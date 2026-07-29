"""频率库构建模块 — 将统计结果转换为 annovar 数据库格式"""

from __future__ import annotations

import json
import logging
import os
import shlex
import shutil
import subprocess
from multiprocessing import Pool

import numpy as np

logger = logging.getLogger(__name__)

# 文件名 → stat.json key 映射
FILENAME_TO_KEY = {
    "120.mutation_frequency.xls": "120",
    "120.mutation_frequency_germline.xls": "120_germline",
    "180.mutation_frequency.xls": "180",
    "180.mutation_frequency_germline.xls": "180_germline",
    "680.mutation_frequency.xls": "680",
    "680.mutation_frequency_germline.xls": "680_germline",
    "1100.mutation_frequency.xls": "1100",
    "1100.mutation_frequency_germline.xls": "1100_germline",
    "WES.mutation_frequency.xls": "wes",
    "WES.mutation_frequency_germline.xls": "wes_germline",
    "84.mutation_frequency.xls": "84",
    "84.mutation_frequency_germline.xls": "84_germline",
    "624.mutation_frequency.xls": "624",
    "624.mutation_frequency_germline.xls": "624_germline",
    "cml206.mutation_frequency.xls": "cml206",
    "kywes.mutation_frequency.xls": "kywes",
    "rnaseq.mutation_frequency.xls": "rnaseq",
    "rnapanel.245.mutation_frequency.xls": "rnapanel245",
    "rnapanel.606.mutation_frequency.xls": "rnapanel606",
    "rnaseqfusion.mutation_frequency.xls": "rnaseqfusion",
    "rnapanelfusion.245.mutation_frequency.xls": "rnapanelfusion245",
    "rnapanelfusion.606.mutation_frequency.xls": "rnapanelfusion606",
    "120fusion.mutation_frequency.xls": "120fusion",
    "1100fusion.mutation_frequency.xls": "1100fusion",
    "wesfusion.mutation_frequency.xls": "wesfusion",
    "84fusion.mutation_frequency.xls": "84fusion",
    "624fusion.mutation_frequency.xls": "624fusion",
}


def _parse_taf_values(taf_str: str) -> list[float]:
    """解析 TAF 列中的百分比字符串为浮点数列表"""
    values = []
    for v in taf_str.split(","):
        v = v.strip()
        if not v or v == ".":
            continue
        if v.endswith("%"):
            values.append(float(v.strip("%")) / 100.0)
        else:
            values.append(float(v))
    return values


def xls_to_vcf(xls_path: str, version: str, history_sample_count: int) -> str:
    """将 mutation_frequency.xls 转为 VCF 格式

    Returns:
        VCF 文件路径
    """
    vcf_path = xls_path.replace(".xls", ".vcf")
    with open(vcf_path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]) + "\n")
        with open(xls_path, "r", encoding="utf-8") as fxls:
            for line in fxls:
                line = line.strip()
                if not line or line.startswith("Chr") or line.startswith("gene"):
                    continue

                arr = line.split("\t")
                chrom, pos, _, ref, alt = arr[:5]
                aachange = arr[8]
                local_frequency = float(arr[10])
                local_af_list_str = arr[11] if len(arr) > 11 else ""

                data = np.array(_parse_taf_values(local_af_list_str))
                pct = {
                    "LOCAL_FREQUENCY": local_frequency,
                    "LOCAL_AF_AVG": float(np.mean(data)) if len(data) > 0 else 0,
                    "LOCAL_AF_PCT0": float(np.min(data)) if len(data) > 0 else 0,
                    "LOCAL_AF_PCT25": float(np.percentile(data, 25)) if len(data) > 0 else 0,
                    "LOCAL_AF_PCT50": float(np.percentile(data, 50)) if len(data) > 0 else 0,
                    "LOCAL_AF_PCT75": float(np.percentile(data, 75)) if len(data) > 0 else 0,
                    "LOCAL_AF_PCT100": float(np.max(data)) if len(data) > 0 else 0,
                    "LOCAL_AN": len(data),
                    "LOCAL_AF_VERSION": f"{version}({history_sample_count})",
                }

                info = (
                    f"LOCAL_FREQUENCY={pct['LOCAL_FREQUENCY']:.4f};"
                    f"LOCAL_AF_AVG={pct['LOCAL_AF_AVG']:.4f};"
                    f"LOCAL_AF_PCT0={pct['LOCAL_AF_PCT0']:.4f};"
                    f"LOCAL_AF_PCT25={pct['LOCAL_AF_PCT25']:.4f};"
                    f"LOCAL_AF_PCT50={pct['LOCAL_AF_PCT50']:.4f};"
                    f"LOCAL_AF_PCT75={pct['LOCAL_AF_PCT75']:.4f};"
                    f"LOCAL_AF_PCT100={pct['LOCAL_AF_PCT100']:.4f};"
                    f"LOCAL_AN={pct['LOCAL_AN']};"
                    f"LOCAL_AF_VERSION={pct['LOCAL_AF_VERSION']};"
                    f"AAChange={aachange}"
                )
                fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")

    logger.debug("VCF 生成: %s", vcf_path)
    return vcf_path


def vcf_to_avinput(vcf_path: str) -> str:
    """将 VCF 转为 AVINPUT 格式（纯 Python 实现，替代 perl convert2annovar.pl）"""
    avinput_path = vcf_path.replace(".vcf", ".avinput")
    with open(avinput_path, "w", encoding="utf-8") as fh, open(vcf_path, "r", encoding="utf-8") as fvcf:
        for line in fvcf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            arr = line.split("\t")
            if len(arr) < 8:
                continue
            chrom, pos, _, ref, alt = arr[:5]
            info = arr[7]
            # 输出格式: chrom\tstart\tend\tref\talt\tinfo
            fh.write(f"{chrom}\t{pos}\t{pos}\t{ref}\t{alt}\t{info}\n")
    logger.debug("AVINPUT 生成: %s", avinput_path)
    return avinput_path


def avinput_to_txt(avinput_path: str) -> str:
    """将 AVINPUT 转为 annovar TXT 格式"""
    txt_path = avinput_path.replace(".avinput", ".raw.txt")
    with open(txt_path, "w", encoding="utf-8") as fh:
        fh.write(
            "\t".join([
                "#Chr", "Start", "AAChange", "Ref", "Alt",
                "LOCAL_FREQUENCY", "LOCAL_AF_AVG", "LOCAL_AF_PCT0",
                "LOCAL_AF_PCT25", "LOCAL_AF_PCT50", "LOCAL_AF_PCT75",
                "LOCAL_AF_PCT100", "LOCAL_AN", "LOCAL_AF_VERSION",
            ]) + "\n"
        )
        with open(avinput_path, "r", encoding="utf-8") as fai:
            for line in fai:
                line = line.strip()
                if not line:
                    continue
                arr = line.split("\t")
                chrom, start, _, ref, alt = arr[:5]
                info = arr[5] if len(arr) > 5 else ""

                if info == "INFO":
                    continue

                info_dict = {}
                for i in info.split(";"):
                    if "=" in i:
                        k, v = i.split("=", 1)
                        info_dict[k.strip()] = v.strip()

                fh.write("\t".join([
                    chrom,
                    start,
                    info_dict.get("AAChange", "."),
                    ref,
                    alt,
                    info_dict.get("LOCAL_FREQUENCY", "."),
                    info_dict.get("LOCAL_AF_AVG", "."),
                    info_dict.get("LOCAL_AF_PCT0", "."),
                    info_dict.get("LOCAL_AF_PCT25", "."),
                    info_dict.get("LOCAL_AF_PCT50", "."),
                    info_dict.get("LOCAL_AF_PCT75", "."),
                    info_dict.get("LOCAL_AF_PCT100", "."),
                    info_dict.get("LOCAL_AN", "."),
                    info_dict.get("LOCAL_AF_VERSION", "."),
                ]) + "\n")

    logger.debug("TXT 生成: %s", txt_path)
    return txt_path


def index_txt(txt_path: str) -> str:
    """对 TXT 文件建立索引（排序）

    Returns:
        索引后文件路径
    """
    indexed_path = txt_path.replace(".raw.txt", ".txt")
    try:
        sort_cmd = (
            f"head -n 1 {shlex.quote(txt_path)} > {shlex.quote(indexed_path)}; "
            f"tail -n +2 {shlex.quote(txt_path)} | sort -k1,1 -k2,2n >> {shlex.quote(indexed_path)}"
        )
        subprocess.run(sort_cmd, shell=True, check=True)
        logger.debug("索引完成: %s", indexed_path)
    except subprocess.CalledProcessError as e:
        logger.error("索引失败: %s", e)
        # 回退：直接复制
        shutil.copy2(txt_path, indexed_path)

    return indexed_path


def process_one_xls(xls_path: str, history_sample_count: int) -> dict:
    """处理单个 xls 文件 → VCF → AVINPUT → TXT → 索引

    Returns:
        {"input": xls_path, "output": indexed_path, "version": str}
    """
    version = os.path.basename(os.path.dirname(xls_path))
    logger.info("开始处理: %s (样本数=%d)", xls_path, history_sample_count)

    vcf = xls_to_vcf(xls_path, version, history_sample_count)
    avinput = vcf_to_avinput(vcf)
    raw_txt = avinput_to_txt(avinput)
    indexed = index_txt(raw_txt)

    # 清理中间文件
    for f in (vcf, avinput, raw_txt):
        try:
            os.remove(f)
        except OSError:
            pass

    logger.info("处理完成: %s", indexed)
    return {"input": xls_path, "output": indexed, "version": version}


def _process_wrapper(args: tuple) -> dict:
    """多进程包装器"""
    xls_path, count = args
    return process_one_xls(xls_path, count)


def build_frequency_db(
    input_dir: str,
    output_dir: str,
    stat_path: str | None = None,
    processes: int = 4,
) -> list[str]:
    """构建本地频率数据库

    Args:
        input_dir: mutation_frequency.xls 文件所在目录
        output_dir: 输出目录（项目编号目录）
        stat_path: stat.json 路径，用于获取样本数
        processes: 并行进程数

    Returns:
        生成的文件路径列表
    """
    os.makedirs(output_dir, exist_ok=True)

    # 加载 stat.json
    stat: dict[str, int] = {}
    if stat_path and os.path.exists(stat_path):
        with open(stat_path, encoding="utf-8") as f:
            stat = json.load(f)
        logger.info("加载 stat.json: %d 个条目", len(stat))
    else:
        logger.warning("未找到 stat.json，样本数设为 1")

    # 收集 xls 文件
    xls_files = []
    for fn in os.listdir(input_dir):
        if not fn.endswith(".xls"):
            continue
        key = FILENAME_TO_KEY.get(fn)
        if not key:
            logger.warning("未知文件，跳过: %s", fn)
            continue

        count = stat.get(key, 1)
        xls_files.append((os.path.join(input_dir, fn), count))

    if not xls_files:
        logger.warning("未找到可处理的 xls 文件: %s", input_dir)
        return []

    # 多进程处理
    results = []
    with Pool(processes=processes) as pool:
        results = pool.map(_process_wrapper, xls_files)

    # 移动到输出目录
    output_files = []
    for r in results:
        if not r or not r.get("output"):
            continue
        src = r["output"]
        basename = os.path.basename(src)
        dest = os.path.join(output_dir, basename)
        shutil.copy2(src, dest)

        # 生成 .idx 索引文件
        idx_dest = dest + ".idx"
        try:
            sort_cmd = (
                f"head -n 1 {shlex.quote(dest)} > {shlex.quote(idx_dest)}; "
                f"tail -n +2 {shlex.quote(dest)} | sort -k1,1 -k2,2n >> {shlex.quote(idx_dest)}"
            )
            subprocess.run(sort_cmd, shell=True, check=True)
        except subprocess.CalledProcessError:
            shutil.copy2(dest, idx_dest)

        output_files.append(dest)
        output_files.append(idx_dest)
        logger.info("输出: %s", dest)

    # 清理临时文件
    for r in results:
        if r and r.get("output"):
            try:
                os.remove(r["output"])
            except OSError:
                pass

    logger.info("频率库构建完成: %d 个文件", len(output_files))
    return output_files
