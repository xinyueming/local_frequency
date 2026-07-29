"""流程编排与 CLI 入口"""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

from tqdm import tqdm

from local_frequency.config import load_config
from local_frequency.db_builder import build_frequency_db
from local_frequency.file_collector import FileCollector
from local_frequency.file_copier import FileCopier
from local_frequency.ssh_client import SSHClient
from local_frequency.stats import stat_cnv, stat_dnafusion, stat_rnafusion, stat_snv

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Pipeline 运行时配置"""
    ssh: dict = field(default_factory=dict)
    remote: dict = field(default_factory=dict)
    local: dict = field(default_factory=dict)
    stats: dict = field(default_factory=lambda: {"max_workers": 4})
    db_builder: dict = field(default_factory=lambda: {"processes": 2})


class Pipeline:
    """主流程编排器

    执行流程:
    1. SSH 连接 → 收集文件列表
    2. 创建本地目录 → SFTP 拷贝文件
    3. 并行统计 (SNV | DNA fusion | RNA fusion)
    4. 等待统计完成 → 构建频率库
    5. CNV 统计 → 生成 CNV 频率库
    """

    # 需要收集的 SNV/Fusion 文件后缀
    SNV_SUFFIXES: ClassVar[list[str]] = [".filter.xls", ".filter.germline.xls"]
    FUSION_SUFFIXES: ClassVar[list[str]] = [".total.fusion.xls", ".tsv.redup.xls"]
    CNV_SUFFIXES: ClassVar[list[str]] = [".cnv.vcf"]

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.source_dict: dict = {}

    def run(
        self,
        dry_run: bool = False,
        collect_only: bool = False,
        project_ids: list[str] | None = None,
    ) -> bool:
        """执行完整流程

        Args:
            dry_run: 仅预检查，不执行实际操作
            collect_only: 仅收集文件，不统计
            project_ids: 指定项目 ID 过滤

        Returns:
            是否成功
        """
        if dry_run:
            return self._dry_run()

        steps = [
            ("步骤 1/5: SSH 连接与文件收集", self._step_collect),
            ("步骤 2/5: SFTP 拷贝文件", self._step_copy),
            ("步骤 3/5: 并行突变统计", self._step_parallel_stats),
            ("步骤 4/5: 构建频率库", self._step_build_db),
            ("步骤 5/5: CNV 统计与索引", self._step_cnv),
        ]

        if collect_only:
            steps = steps[:2]

        for name, step_fn in tqdm(steps, desc="Pipeline", unit="step"):
            logger.info(">>> %s", name)
            try:
                if not step_fn(project_ids=project_ids):
                    logger.error("%s 失败", name)
                    return False
            except Exception as e:  # noqa: BLE001 — Pipeline 层需捕获所有异常
                logger.error("%s 异常: %s", name, e)
                return False

        logger.info("Pipeline 执行完成")
        return True

    def _dry_run(self) -> bool:
        """预检查模式"""
        logger.info("=== Dry-Run 预检查 ===")

        ssh_cfg = self.config.ssh
        if not ssh_cfg.get("host"):
            logger.error("配置缺少 ssh.host")
            return False

        if not self.config.remote.get("base_path"):
            logger.error("配置缺少 remote.base_path")
            return False

        if not self.config.local.get("base_path"):
            logger.error("配置缺少 local.base_path")
            return False

        logger.info("SSH: %s@%s:%d", ssh_cfg.get("user", "root"), ssh_cfg["host"], ssh_cfg.get("port", 22))
        logger.info("远程目录: %s", self.config.remote["base_path"])
        logger.info("本地目录: %s", self.config.local["base_path"])
        logger.info("最大并行: %d", self.config.stats.get("max_workers", 4))
        logger.info("Dry-Run 检查通过")
        return True

    def _step_collect(self, project_ids: list[str] | None = None) -> bool:
        """SSH 连接并收集文件列表"""
        ssh_cfg = self.config.ssh
        remote_path = self.config.remote["base_path"]

        all_suffixes = self.SNV_SUFFIXES + self.FUSION_SUFFIXES + self.CNV_SUFFIXES
        collector = FileCollector(base_path=remote_path, suffixes=all_suffixes)

        with SSHClient(
            host=ssh_cfg.get("host", ""),
            port=ssh_cfg.get("port", 22),
            user=ssh_cfg.get("user", "root"),
            key_path=ssh_cfg.get("key_path"),
            password=ssh_cfg.get("password"),
        ) as client:
            files = client.walk_remote_dir(remote_path)
            self.source_dict = collector.classify_files(files)

        if not self.source_dict:
            logger.warning("未收集到任何文件")
            return False

        n_projects = len(self.source_dict)
        n_files = sum(len(f) for p in self.source_dict.values() for f in p.values())
        logger.info("收集完成: %d 个项目, %d 个文件", n_projects, n_files)
        return True

    def _step_copy(self, project_ids: list[str] | None = None) -> bool:
        """SFTP 拷贝文件到本地"""
        ssh_cfg = self.config.ssh
        remote_path = self.config.remote["base_path"]
        local_base = self.config.local["base_path"]

        with SSHClient(
            host=ssh_cfg.get("host", ""),
            port=ssh_cfg.get("port", 22),
            user=ssh_cfg.get("user", "root"),
            key_path=ssh_cfg.get("key_path"),
            password=ssh_cfg.get("password"),
        ) as client:
            copier = FileCopier(client.sftp, local_base)
            copied = copier.copy_files(self.source_dict, remote_path, project_ids=project_ids)

        logger.info("拷贝完成: %d 个文件", len(copied))
        return True

    def _step_parallel_stats(self, project_ids: list[str] | None = None) -> bool:
        """并行执行 SNV / DNA fusion / RNA fusion 统计"""
        local_base = self.config.local["base_path"]
        max_workers = self.config.stats.get("max_workers", 4)

        projects = project_ids if project_ids else list(self.source_dict.keys())
        if not projects:
            logger.warning("无项目可统计")
            return True

        stat_tasks = []
        for proj_id in projects:
            proj_dir = str(Path(local_base) / proj_id / "data")
            out_dir = str(Path(local_base) / proj_id / "mutation_frequency_result")
            stat_tasks.append((proj_id, proj_dir, out_dir))

        results = {}
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for proj_id, data_dir, out_dir in stat_tasks:
                # 提交 SNV 统计
                futures[executor.submit(stat_snv, data_dir, out_dir)] = f"{proj_id}_snv"
                # 提交 DNA fusion 统计
                futures[executor.submit(stat_dnafusion, data_dir, out_dir)] = f"{proj_id}_dnafusion"
                # 提交 RNA fusion 统计
                futures[executor.submit(stat_rnafusion, data_dir, out_dir)] = f"{proj_id}_rnafusion"

            for future in tqdm(as_completed(futures), total=len(futures), desc="统计中"):
                task_name = futures[future]
                try:
                    result = future.result()
                    results[task_name] = result
                    logger.info("%s: %d 个样本", task_name, result.get("total", 0))
                except Exception as e:  # noqa: BLE001 — 并行任务需捕获所有异常
                    logger.error("%s 失败: %s", task_name, e)
                    results[task_name] = {"total": 0, "error": str(e)}

        logger.info("并行统计完成: %d 个任务", len(results))
        return True

    def _step_build_db(self, project_ids: list[str] | None = None) -> bool:
        """构建频率库"""
        local_base = self.config.local["base_path"]
        processes = self.config.db_builder.get("processes", 2)

        projects = project_ids if project_ids else list(self.source_dict.keys())
        for proj_id in projects:
            input_dir = str(Path(local_base) / proj_id / "mutation_frequency_result")
            output_dir = str(Path(local_base) / proj_id / "frequency_db")
            stat_path = str(Path(input_dir) / "stat.json")

            build_frequency_db(input_dir, output_dir, stat_path=stat_path, processes=processes)

        return True

    def _step_cnv(self, project_ids: list[str] | None = None) -> bool:
        """CNV 统计与索引"""
        local_base = self.config.local["base_path"]
        projects = project_ids if project_ids else list(self.source_dict.keys())

        for proj_id in projects:
            cnv_dir = str(Path(local_base) / proj_id / "data")
            out_dir = str(Path(local_base) / proj_id / "mutation_frequency_result")
            stat_cnv(cnv_dir, out_dir, sample_prefix=proj_id)

        return True


def create_pipeline(config_path: str) -> Pipeline:
    """从配置文件创建 Pipeline"""
    cfg = load_config(config_path)
    pc = PipelineConfig()

    if "ssh" in cfg:
        pc.ssh = cfg["ssh"]
    if "remote" in cfg:
        pc.remote = cfg["remote"]
    if "local" in cfg:
        pc.local = cfg["local"]
    if "stats" in cfg:
        pc.stats = cfg["stats"]
    if "db_builder" in cfg:
        pc.db_builder = cfg["db_builder"]

    return Pipeline(pc)


def main():
    """CLI 入口"""
    parser = argparse.ArgumentParser(description="Local Frequency - 本地突变频率分析")
    parser.add_argument("--config", required=True, help="配置文件路径")
    parser.add_argument("--dry-run", action="store_true", help="预检查模式")
    parser.add_argument("--collect-only", action="store_true", help="仅收集文件")
    parser.add_argument("--project", nargs="+", help="指定项目 ID 列表")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    pipeline = create_pipeline(args.config)
    success = pipeline.run(
        dry_run=args.dry_run,
        collect_only=args.collect_only,
        project_ids=args.project,
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
