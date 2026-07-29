"""文件收集与分类模块"""

from __future__ import annotations

import logging
import os
from collections import defaultdict

logger = logging.getLogger(__name__)

# 支持的文件后缀 → 分类键 映射
# 顺序重要：更具体的后缀须排在前面
KNOWN_SUFFIXES = [
    ".filter.germline.xls",
    ".total.fusion.xls",
    ".tsv.redup.xls",
    ".cnv.vcf",
    ".filter.xls",
]


class FileCollector:
    """按文件后缀分类，构建项目编号 → 文件字典"""

    def __init__(
        self,
        base_path: str,
        suffixes: list[str] | None = None,
        project_level: int = 5,
    ):
        """
        Args:
            base_path: 远程数据根目录（用于计算相对路径）
            suffixes: 需要收集的文件后缀列表。None 时使用默认列表。
            project_level: 从 base_path 开始的相对目录层级，该层为项目编号
        """
        self.base_path = base_path.rstrip("/")
        self.suffixes = suffixes if suffixes is not None else KNOWN_SUFFIXES
        self.project_level = project_level

    def classify_files(self, file_paths: list[str]) -> dict[str, dict[str, list[str]]]:
        """将文件路径按项目编号和后缀分类

        Args:
            file_paths: 远程文件绝对路径列表

        Returns:
            {项目编号: {分类键: [文件路径列表]}}
        """
        result: dict[str, dict[str, list[str]]] = defaultdict(
            lambda: defaultdict(list)
        )

        for fpath in file_paths:
            rel_path = os.path.relpath(fpath, self.base_path)
            parts = rel_path.split(os.sep)

            # 按 project_level 提取项目编号
            if len(parts) <= self.project_level:
                continue
            project_id = parts[self.project_level]

            category = self._match_suffix(os.path.basename(fpath))
            if category:
                result[project_id][category].append(fpath)

        logger.info("文件分类完成: %d 个项目", len(result))
        return dict(result)

    def _match_suffix(self, filename: str) -> str | None:
        """匹配文件后缀，返回分类键"""
        for suffix in self.suffixes:
            if filename.endswith(suffix):
                return suffix
        return None
