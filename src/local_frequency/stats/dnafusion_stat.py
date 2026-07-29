"""DNA Fusion 突变频率统计模块"""

from __future__ import annotations

import logging
import os
import re
from collections import Counter

import pandas as pd

logger = logging.getLogger(__name__)


def read_fusion_files(sour_path: str, pattern: str = "fusion") -> list[str]:
    """递归查找 fusion xls 文件

    Args:
        sour_path: 搜索根目录
        pattern: "dna" 匹配 total.fusion.xls，"rna" 匹配 rna.starfusion.tsv.redup.xls

    Returns:
        匹配文件路径列表
    """
    files = []
    if pattern == "dna":
        regex = re.compile(r".*fusion\.xls$")
    elif pattern == "rna":
        regex = re.compile(r".*rna\.starfusion\.tsv\.redup\.xls$")
    else:
        regex = re.compile(r".*redup\.xls$|.*fusion\.xls$")

    for root, _dirs, filenames in os.walk(sour_path):
        for fname in filenames:
            if regex.match(fname):
                files.append(os.path.join(root, fname))

    logger.info("找到 %d 个 %s fusion 文件", len(files), pattern)
    return files


def _parse_xls(filepath: str) -> list[dict]:
    if pd.__version__ >= "1.3.0":
        tb = pd.read_table(filepath, encoding="utf-8", encoding_errors="ignore", index_col=False)
    else:
        tb = pd.read_table(filepath, encoding="utf-8", index_col=False)
    return tb.to_dict(orient="records")


def _make_fusion_key(line: dict) -> str:
    return f'{line["gene1_chr"]}:{line["gene1_pos"]}_{line["gene2_chr"]}:{line["gene2_pos"]}'


def stat_dnafusion(data_dir: str, output_dir: str, sample_prefix: str) -> dict:
    """执行 DNA Fusion 频率统计

    Args:
        data_dir: 输入文件所在目录
        output_dir: 输出目录
        sample_prefix: 项目编号

    Returns:
        {"total": 样本数, "output_file": 输出路径}
    """
    files = read_fusion_files(data_dir, pattern="dna")
    if not files:
        logger.warning("未找到 DNA fusion 文件: %s", data_dir)
        return {"total": 0, "output_file": ""}

    header = [
        "gene1_chr", "gene1_pos", "gene1",
        "gene2_chr", "gene2_pos", "gene2",
        "exon1", "exon2", "FusionChange",
        "population", "mutation_frequency", "TAF",
    ]

    di: dict[str, list] = {}
    li: list[str] = []
    taf: dict[str, list[str]] = {}

    for fi in files:
        content = _parse_xls(fi)
        temp_list: list[str] = []
        for line in content:
            fusion_change = _make_fusion_key(line)
            if fusion_change in temp_list:
                continue
            temp_list.append(fusion_change)

            exon_parts = re.split(r"[-:]+", str(line.get("exon", "")))
            ls1 = [
                str(line.get("gene1_chr", "")),
                str(line.get("gene1_pos", "")),
                str(line.get("gene1", "")),
                str(line.get("gene2_chr", "")),
                str(line.get("gene2_pos", "")),
                str(line.get("gene2", "")),
                exon_parts[0] if len(exon_parts) > 0 else "",
                exon_parts[1] if len(exon_parts) > 1 else "",
            ]

            di[fusion_change] = ls1
            li.append(fusion_change)
            taf.setdefault(fusion_change, []).append(str(line.get("freq", "")))

    counter = Counter(li).most_common()
    total = len(files)

    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"{sample_prefix}_dnafusion.mutation_frequency.xls")

    with open(out_file, "w", encoding="utf-8") as f_out:
        f_out.write("\t".join(header) + "\n")
        for k, v in counter:
            pop = int(v)
            mut_freq = round(pop / total, 4)
            f_out.write("\t".join(di[k]) + f"\t{k}\t{pop}\t{mut_freq}\t{','.join(taf[k])}\n")

    logger.info("DNA Fusion 统计完成: %d 个样本 → %s", total, out_file)
    return {"total": total, "output_file": out_file}
