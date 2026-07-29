"""SFTP 文件拷贝模块"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)


class FileCopier:
    """通过 SFTP 拷贝文件到本地，支持增量拷贝"""

    def __init__(self, sftp, local_base: str):
        """
        Args:
            sftp: paramiko SFTP 客户端
            local_base: 本地存储根目录
        """
        self.sftp = sftp
        self.local_base = os.path.normpath(local_base)

    def copy_files(
        self,
        source_dict: dict[str, dict[str, list[str]]],
        remote_path: str,
        project_ids: list[str] | None = None,
    ) -> dict[str, str]:
        """批量拷贝文件

        Args:
            source_dict: FileCollector.classify_files() 的输出
            remote_path: 远程根目录（用于计算相对路径）
            project_ids: 仅拷贝指定项目，None 表示全部

        Returns:
            {远程路径: 本地路径} 成功拷贝的文件映射
        """
        copied: dict[str, str] = {}
        norm_base = remote_path.rstrip("/")

        projects = project_ids if project_ids else list(source_dict.keys())

        for proj_id in projects:
            if proj_id not in source_dict:
                logger.warning("项目 %s 无数据，跳过", proj_id)
                continue

            proj_dir = source_dict[proj_id]
            for file_list in proj_dir.values():
                for remote_file in file_list:
                    local_path = self._resolve_local_path(
                        remote_file, norm_base, proj_id
                    )
                    if self._skip_copy(local_path):
                        continue
                    if self.copy_single(remote_file, local_path):
                        copied[remote_file] = local_path

        logger.info("拷贝完成: %d 个文件", len(copied))
        return copied

    def copy_single(self, remote_path: str, local_path: str) -> bool:
        """拷贝单个文件

        Args:
            remote_path: 远程文件路径
            local_path: 本地目标路径

        Returns:
            是否成功
        """
        try:
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            self.sftp.get(remote_path, local_path)
            logger.debug("拷贝: %s → %s", remote_path, local_path)
            return True
        except Exception as e:  # noqa: BLE001 — SFTP 拷贝需捕获所有异常
            logger.error("拷贝失败 %s: %s", remote_path, e)
            return False

    def _resolve_local_path(
        self, remote_file: str, remote_base: str, project_id: str
    ) -> str:
        """解析本地目标路径

        目录结构: local_base/项目编号/data/文件名
        """
        filename = os.path.basename(remote_file)
        return os.path.join(self.local_base, project_id, "data", filename)

    def _skip_copy(self, local_path: str) -> bool:
        """检查是否跳过（文件已存在）"""
        if os.path.exists(local_path):
            logger.debug("跳过已存在: %s", local_path)
            return True
        return False
