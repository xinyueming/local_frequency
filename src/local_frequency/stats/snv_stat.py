"""SNV 突变频率统计模块"""

from __future__ import annotations

import logging
import os
import re
from collections import Counter

import pandas as pd

logger = logging.getLogger(__name__)


def read_filter_files(sour_path: str, germline: bool = False) -> list[str]:
    """递归查找 filter.xls 文件

    Args:
        sour_path: 搜索根目录
        germline: True 查找 germline 文件，False 查找体细胞文件

    Returns:
        匹配文件路径列表
    """
    files = []
    pattern = re.compile(r".*hg(19|38)_multianno\.filter.*germline.*") if germline else re.compile(
        r".*hg(19|38)_multianno\.filter(?!.*germline).*"
    )

    for root, _dirs, filenames in os.walk(sour_path):
        for fname in filenames:
            if pattern.match(fname):
                files.append(os.path.join(root, fname))

    logger.info("找到 %d 个 %s 文件", len(files), "germline" if germline else "somatic")
    return files


def _parse_xls(filepath: str) -> list[dict]:
    """读取 xls 文件为 records 列表"""
    if pd.__version__ >= "1.3.0":
        tb = pd.read_table(filepath, encoding="utf-8", encoding_errors="ignore", index_col=False)
    else:
        tb = pd.read_table(filepath, encoding="utf-8", index_col=False)
    return tb.to_dict(orient="records")


def _extract_cdna_info(line: dict) -> str:
    """从记录中提取 cDNA 变异信息"""
    aachange = line.get("AAChange.refGeneWithVer", "")

    if "c." in aachange:
        # 清理 fs* 后缀
        if "fs*" in aachange:
            aachange = aachange.split("*")[0]
        return aachange

    # 回退：使用 new_* 字段构造
    if all(line.get(f"new_{k}") is not None for k in ("gene", "transcript", "exon", "cdna")):
        return f'{line["new_gene"]}:{line["new_transcript"]}:{line["new_exon"]}:{line["new_cdna"]}'

    # 最终回退：用 Chr_Pos_Ref_Alt
    chrom = str(line.get("Chr", "")).strip()
    pos = line.get("Pos", "")
    ref = line.get("Ref", "")
    alt = line.get("Alt", "")
    return f"{chrom}_{pos}_{ref}_{alt}"


def stat_snv(data_dir: str, output_dir: str, sample_prefix: str, germline: bool = False) -> dict:
    """执行 SNV 频率统计

    Args:
        data_dir: 输入文件所在目录（递归搜索）
        output_dir: 输出目录
        sample_prefix: 项目编号（如 84）
        germline: 是否统计胚系变异

    Returns:
        {"total": 样本数, "output_file": 输出路径}
    """
    files = read_filter_files(data_dir, germline)
    if not files:
        logger.warning("未找到匹配的 SNV 文件: %s", data_dir)
        return {"total": 0, "output_file": ""}

    header = [
        "Chr", "Pos", "End_pos", "Ref", "Alt",
        "Gene.refGeneWithVer", "Func.refGeneWithVer",
        "ExonicFunc.refGeneWithVer", "AAChange.refGeneWithVer",
        "population", "mutation_frequency", "TAF",
    ]

    di: dict[str, list] = {}
    li: list[str] = []
    taf: dict[str, list[str]] = {}

    for fi in files:
        content = _parse_xls(fi)
        for line in content:
            # 跳过多基因行
            if "\\x3b" in str(line.get("Gene.refGeneWithVer", "")):
                chrom = str(line.get("Chr", "")).strip()
                pos = line.get("Pos", "")
                ref = line.get("Ref", "")
                alt = line.get("Alt", "")
                key = f"{chrom}_{pos}_{ref}_{alt}"
            else:
                key = _extract_cdna_info(line)

            ls1 = [
                str(line.get("Chr", "")),
                str(line.get("Pos", "")),
                str(line.get("End_pos", "")),
                str(line.get("Ref", "")),
                str(line.get("Alt", "")),
                str(line.get("Gene.refGeneWithVer", "")),
                str(line.get("Func.refGeneWithVer", "")),
                str(line.get("ExonicFunc.refGeneWithVer", "")),
                key,
            ]

            di[key] = ls1
            li.append(key)

            taff = str(line.get("TAF", ""))
            taf.setdefault(key, []).append(taff)

    counter = Counter(li).most_common()
    total = len(files)

    os.makedirs(output_dir, exist_ok=True)
    suffix = "snv_germline" if germline else "snv_somatic"
    out_file = os.path.join(output_dir, f"{sample_prefix}_{suffix}.mutation_frequency.xls")

    with open(out_file, "w", encoding="utf-8") as f_out:
        f_out.write("\t".join(header) + "\n")
        for k, v in counter:
            pop = int(v)
            mut_freq = round(pop / total, 4)
            f_out.write("\t".join(di[k]) + f"\t{k}\t{pop}\t{mut_freq}\t{','.join(taf[k])}\n")

    logger.info("SNV 统计完成: %d 个样本 → %s", total, out_file)
    return {"total": total, "output_file": out_file}
