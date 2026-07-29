"""CNV 突变频率统计模块"""

from __future__ import annotations

import gzip
import logging
import os
import re
import shlex
import shutil
import subprocess

logger = logging.getLogger(__name__)


def _process_vcf(input_path: str) -> dict:
    """解析单个 VCF 文件，返回每个区间的 DUP/DEL 出现情况

    Returns:
        {(chrom, start, end): {"DUP": int, "DEL": int, "DUP_CN": [...], "DEL_CN": [...]}}
    """
    records: dict = {}
    basename = os.path.basename(input_path)
    opener = gzip.open if basename.endswith((".gz", ".bgz")) else open

    with opener(input_path, "rt", encoding="ISO-8859-1") as reader:
        for line in reader:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            chrom = cols[0]
            try:
                start = int(cols[1])
            except ValueError:
                continue

            alt = cols[4]
            info = cols[7]

            m = re.search(r"END=(\d+)", info)
            end = int(m.group(1)) if m else start

            # 提取 CN 值（sample 列第三个字段）
            cn = None
            if len(cols) > 9:
                sample_parts = cols[9].split(":")
                if len(sample_parts) >= 3:
                    cn_str = sample_parts[2]
                    try:
                        cn = float(cn_str) if "." in cn_str else int(cn_str)
                    except ValueError:
                        pass

            key = (chrom, start, end)

            if alt == "<DUP>":
                if key not in records:
                    records[key] = {"DUP": 0, "DEL": 0, "DUP_CN": [], "DEL_CN": []}
                if records[key]["DUP"] == 0:
                    records[key]["DUP"] = 1
                    if cn is not None:
                        records[key]["DUP_CN"].append(cn)
            elif alt == "<DEL>":
                if key not in records:
                    records[key] = {"DUP": 0, "DEL": 0, "DUP_CN": [], "DEL_CN": []}
                if records[key]["DEL"] == 0:
                    records[key]["DEL"] = 1
                    if cn is not None:
                        records[key]["DEL_CN"].append(cn)

    return records


def _records_to_df(records: dict, total_files: int) -> list[dict]:
    """将聚合记录转换为行字典列表"""
    import pandas as pd

    rows = []
    for (chrom, start, end), cnt in records.items():
        local_gain_ac = int(cnt.get("DUP", 0))
        local_loss_ac = int(cnt.get("DEL", 0))
        dup_cn_values = cnt.get("DUP_CN", [])
        del_cn_values = cnt.get("DEL_CN", [])

        dup_cn_str = ",".join(str(x) for x in dup_cn_values) if dup_cn_values else "."
        del_cn_str = ",".join(str(x) for x in del_cn_values) if del_cn_values else "."

        local_an = int(total_files)
        if local_an > 0:
            local_gain_af = local_gain_ac / local_an
            local_loss_af = local_loss_ac / local_an
        else:
            local_gain_af = 0.0
            local_loss_af = 0.0

        rows.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "local_gain_ac": local_gain_ac,
            "local_loss_ac": local_loss_ac,
            "local_an": local_an,
            "local_gain_af": local_gain_af,
            "local_loss_af": local_loss_af,
            "dup_cn": dup_cn_str,
            "del_cn": del_cn_str,
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    return df.to_dict(orient="records")


def stat_cnv(data_dir: str, output_dir: str, sample_prefix: str | None = None) -> dict:
    """执行 CNV 频率统计

    Args:
        data_dir: 输入 VCF 文件所在目录
        output_dir: 输出目录
        sample_prefix: 样本前缀，用于生成 local_freq.<prefix>.cnv.txt.gz

    Returns:
        {"total": 样本数, "output_file": 输出路径（.gz 文件）}
    """
    # 查找 VCF 文件
    file_list = []
    if os.path.isdir(data_dir):
        for fn in os.listdir(data_dir):
            if fn.endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
                file_list.append(os.path.join(data_dir, fn))
    else:
        file_list.append(data_dir)

    file_list = [f for f in file_list if os.path.exists(f)]
    if not file_list:
        logger.warning("未找到 CNV VCF 文件: %s", data_dir)
        return {"total": 0, "output_file": ""}

    # 聚合统计
    records: dict = {}
    for f in file_list:
        per_file = _process_vcf(f)
        for key, cnt in per_file.items():
            if key not in records:
                records[key] = {"DUP": 0, "DEL": 0, "DUP_CN": [], "DEL_CN": []}
            records[key]["DUP"] += cnt.get("DUP", 0)
            records[key]["DEL"] += cnt.get("DEL", 0)
            records[key]["DUP_CN"].extend(cnt.get("DUP_CN", []))
            records[key]["DEL_CN"].extend(cnt.get("DEL_CN", []))

    total_files = len(file_list)
    rows = _records_to_df(records, total_files)

    # 写输出
    os.makedirs(output_dir, exist_ok=True)

    if sample_prefix:
        out_txt = os.path.join(output_dir, f"local_freq.{sample_prefix}.cnv.txt")
    else:
        out_txt = os.path.join(output_dir, "local_freq.result.txt")

    import pandas as pd

    df = pd.DataFrame(rows)
    df.to_csv(out_txt, sep="\t", index=False, float_format="%.9f")
    logger.info("CNV 统计完成: %d 个样本 → %s (%d 条记录)", total_files, out_txt, len(rows))

    # bgzip + tabix 索引
    gz_path = _compress_and_index(out_txt)

    return {"total": total_files, "output_file": gz_path}


def _compress_and_index(txt_path: str) -> str:
    """对输出文件进行 bgzip 压缩和 tabix 索引

    Returns:
        .gz 文件路径
    """
    # 优先使用 conda 路径，回退到 PATH
    bgzip_bin = shutil.which("bgzip") or os.path.expanduser("~/miniconda3/bin/bgzip")
    tabix_bin = shutil.which("tabix") or os.path.expanduser("~/miniconda3/bin/tabix")

    sorted_path = txt_path + ".sorted"
    gz_path = txt_path + ".gz" if not txt_path.endswith(".gz") else txt_path

    try:
        sort_cmd = (
            f"head -n 1 {shlex.quote(txt_path)} > {shlex.quote(sorted_path)}; "
            f"tail -n +2 {shlex.quote(txt_path)} | sort -k1,1 -k2,2n >> {shlex.quote(sorted_path)}"
        )
        subprocess.run(sort_cmd, shell=True, check=True)

        bgzip_cmd = f"{shlex.quote(bgzip_bin)} -c {shlex.quote(sorted_path)} > {shlex.quote(gz_path)}"
        subprocess.run(bgzip_cmd, shell=True, check=True)

        subprocess.run([tabix_bin, "-s", "1", "-b", "2", "-e", "3", "-S", "1", gz_path], check=True)

        # 清理临时文件
        try:
            os.remove(txt_path)
        except OSError:
            pass
        try:
            os.remove(sorted_path)
        except OSError:
            pass

        logger.info("生成 bgzip + tabix: %s", gz_path)
    except subprocess.CalledProcessError as e:
        logger.error("压缩/索引失败: %s", e)
    except FileNotFoundError:
        logger.warning("bgzip/tabix 未安装，跳过压缩索引")

    return gz_path
