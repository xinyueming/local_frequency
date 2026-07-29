"""端到端集成测试（全链路 mock）"""

import json
from unittest.mock import MagicMock, patch

from local_frequency.pipeline import Pipeline, PipelineConfig


def _create_mock_ssh():
    """创建完整的 mock SSH 客户端"""
    mock = MagicMock()
    mock.__enter__ = MagicMock(return_value=mock)
    mock.__exit__ = MagicMock(return_value=False)
    mock.sftp = MagicMock()
    return mock


def _create_test_files(data_dir):
    """在临时目录创建完整的测试数据"""

    data_dir.mkdir(parents=True, exist_ok=True)

    # SNV filter file
    (data_dir / "sample1.hg38_multianno.filter.xls").write_text(
        "Chr\tPos\tEnd_pos\tRef\tAlt\tGene.refGeneWithVer\tFunc.refGeneWithVer\t"
        "ExonicFunc.refGeneWithVer\tAAChange.refGeneWithVer\tpopulation\t"
        "mutation_frequency\tTAF\n"
        "chr1\t100\t100\tA\tT\tTP53\texonic\tnonsynonymous SNV\t"
        "TP53:NM_000546:exon5:c.500A>T\t1\t1.0\t0.5\n"
    )

    # DNA fusion file
    (data_dir / "sample1.total.fusion.xls").write_text(
        "gene1_chr\tgene1_pos\tgene1\tgene2_chr\tgene2_pos\tgene2\texon\tfreq\n"
        "chr1\t1000\tALK\tchr2\t2000\tROS1\t1-3\t0.3\n"
    )

    # RNA fusion file
    (data_dir / "sample1.tsv.redup.xls").write_text(
        "gene1_chr\tgene1_pos\tgene1\tgene2_chr\tgene2_pos\tgene2\texon\tfreq\n"
        "chr1\t1000\tEML4\tchr2\t2000\tALK\t1-3\t0.2\n"
    )


class TestEndToEndPipeline:
    """端到端 Pipeline 测试"""

    @patch("local_frequency.pipeline.SSHClient")
    @patch("local_frequency.pipeline.stat_cnv")
    @patch("local_frequency.pipeline.build_frequency_db")
    def test_full_pipeline(self, mock_build_db, mock_cnv, MockClient, tmp_path):
        """完整流程：收集 → 拷贝 → 统计 → 建库 → CNV"""
        # 准备 mock SSH
        mock_ssh = _create_mock_ssh()
        mock_ssh.walk_remote_dir.return_value = [
            "/data/a/b/c/d/proj001/sample1.hg38_multianno.filter.xls",
            "/data/a/b/c/d/proj001/sample1.total.fusion.xls",
            "/data/a/b/c/d/proj001/sample1.tsv.redup.xls",
        ]
        MockClient.return_value = mock_ssh

        # 创建本地测试数据
        data_dir = tmp_path / "proj001" / "data"
        _create_test_files(data_dir)

        # stat.json
        stat = {"84": 1}
        result_dir = tmp_path / "proj001" / "mutation_frequency_result"
        result_dir.mkdir(parents=True)
        (result_dir / "stat.json").write_text(json.dumps(stat))

        # 创建 Pipeline
        cfg = PipelineConfig(
            ssh={"host": "test", "port": 22, "user": "test"},
            remote={"base_path": "/data"},
            local={"base_path": str(tmp_path)},
            stats={"max_workers": 1},
            db_builder={"processes": 1},
        )
        p = Pipeline(cfg)
        # 手动注入 source_dict（跳过实际 SSH）
        p.source_dict = {"proj001": {
            ".filter.xls": ["/data/a/b/c/d/proj001/sample1.hg38_multianno.filter.xls"],
            ".total.fusion.xls": ["/data/a/b/c/d/proj001/sample1.total.fusion.xls"],
            ".tsv.redup.xls": ["/data/a/b/c/d/proj001/sample1.tsv.redup.xls"],
        }}

        # 执行统计步骤
        assert p._step_parallel_stats(project_ids=["proj001"]) is True
        assert p._step_build_db(project_ids=["proj001"]) is True
        assert p._step_cnv(project_ids=["proj001"]) is True

        # 验证统计结果文件存在
        assert (result_dir / "mutation_frequency.xls").exists()

    @patch("local_frequency.pipeline.SSHClient")
    def test_collect_and_copy(self, MockClient, tmp_path):
        """仅收集+拷贝流程"""
        mock_ssh = _create_mock_ssh()
        mock_ssh.walk_remote_dir.return_value = [
            "/data/a/b/c/d/proj001/sample.filter.xls",
        ]
        MockClient.return_value = mock_ssh

        cfg = PipelineConfig(
            ssh={"host": "test", "port": 22, "user": "test"},
            remote={"base_path": "/data"},
            local={"base_path": str(tmp_path / "local")},
        )
        p = Pipeline(cfg)
        assert p.run(collect_only=True) is True
        assert "proj001" in p.source_dict
