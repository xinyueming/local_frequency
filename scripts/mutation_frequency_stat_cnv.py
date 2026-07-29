import argparse
import gzip
import pandas as pd
import numpy as np
import re
import os
import sys
import subprocess
import shlex


def process_vcf(input_path: str):
    """Parse a single VCF (compressed or not) and return presence of DUP/DEL per interval for this file.

    For each variant, extracts SV type from INFO field and copy number (CN) from sample column.
    CN is extracted from the second field in the sample column format (e.g., 0/1:3:1 -> 3)
    Returns dict: key -> {'DUP': 0/1, 'DEL': 0/1, 'DUP_CN': [], 'DEL_CN': []}
    """
    records = {}
    basename = os.path.basename(input_path)
    opener = gzip.open if basename.endswith('.gz') or basename.endswith('.bgz') else open
    
    # First, find the header line to locate sample columns
    header_lines = []
    with opener(input_path, 'rt', encoding='ISO-8859-1') as reader:
        for line in reader:
            line = line.strip()
            if line.startswith('#CHROM'):
                header_cols = line.split('\t')
                # Sample columns start from index 9 (0-based)
                # We'll use the first sample column (index 9) for CN extraction
                break
    
    with opener(input_path, 'rt', encoding='ISO-8859-1') as reader:
        for line in reader:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:  # Need at least until sample column
                continue
            
            chrom = cols[0]
            try:
                start = int(cols[1])
            except ValueError:
                continue
            alt = cols[4]
            info = cols[7]
            
            # Extract END position
            m = re.search(r'END=(\d+)', info)
            end = int(m.group(1)) if m else start
            
            # Extract copy number (CN) from sample column
            cn = None
            if len(cols) > 9:  # Has sample column
                sample_field = cols[9]  # First sample column
                # Split by colon and take the second value (index 1)
                sample_parts = sample_field.split(':')
                if len(sample_parts) >= 2:
                    cn_str = sample_parts[2]  # Third field contains CN
                    try:
                        # Convert to appropriate numeric type
                        if '.' in cn_str:
                            cn = float(cn_str)
                        else:
                            cn = int(cn_str)
                    except ValueError:
                        # If conversion fails, skip this CN
                        pass
            
            key = (chrom, start, end)
            
            if alt == '<DUP>':
                if key not in records:
                    records[key] = {'DUP': 0, 'DEL': 0, 'DUP_CN': [], 'DEL_CN': []}
                if records[key]['DUP'] == 0:  # First time seeing this DUP in this file
                    records[key]['DUP'] = 1
                    if cn is not None:
                        records[key]['DUP_CN'].append(cn)
                    
            elif alt == '<DEL>':
                if key not in records:
                    records[key] = {'DUP': 0, 'DEL': 0, 'DUP_CN': [], 'DEL_CN': []}
                if records[key]['DEL'] == 0:  # First time seeing this DEL in this file
                    records[key]['DEL'] = 1
                    if cn is not None:
                        records[key]['DEL_CN'].append(cn)
            else:
                continue

    return records


def records_to_df(records: dict, total_files: int = None):
    """Convert aggregated records to DataFrame.

    If total_files is provided, use it as the denominator (local_an) for AF calculations.
    Otherwise, local_an defaults to local_gain_ac + local_loss_ac.
    Also adds CN columns showing copy number values across samples.
    """
    rows = []
    for (chrom, start, end), cnt in records.items():
        local_gain_ac = int(cnt.get('DUP', 0))
        local_loss_ac = int(cnt.get('DEL', 0))
        
        # Get CN values
        dup_cn_values = cnt.get('DUP_CN', [])
        del_cn_values = cnt.get('DEL_CN', [])
        
        # Create string representations of CN values
        dup_cn_str = ','.join(str(x) for x in dup_cn_values) if dup_cn_values else '.'
        del_cn_str = ','.join(str(x) for x in del_cn_values) if del_cn_values else '.'

        if total_files is not None:
            local_an = int(total_files)
        else:
            local_an = local_gain_ac + local_loss_ac

        if local_an > 0:
            local_gain_af = local_gain_ac / local_an
            local_loss_af = local_loss_ac / local_an
        else:
            local_gain_af = 0.0
            local_loss_af = 0.0

        rows.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'local_gain_ac': local_gain_ac,
            'local_loss_ac': local_loss_ac,
            'local_an': local_an,
            'local_gain_af': local_gain_af,
            'local_loss_af': local_loss_af,
            'dup_cn': dup_cn_str,  # CN values for DUP variants
            'del_cn': del_cn_str,  # CN values for DEL variants
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser(description='Count DUP/DEL from a VCF and output summary (tab-separated)')
    parser.add_argument('-i', '--input', required=True, nargs='+', help='Input VCF file(s) or directory path')
    parser.add_argument('-o', '--output', default=None, help='Output TSV path (default: input_basename.result.txt)')
    parser.add_argument('-S', '--sample-prefix', default=None, help='If provided, use this numeric prefix to name outputs: local_freq.<prefix>.cnv.txt')
    args = parser.parse_args()

    inputs = args.input

    # expand if a single directory is provided
    file_list = []
    for p in inputs:
        if os.path.isdir(p):
            for fn in os.listdir(p):
                if fn.endswith('.vcf') or fn.endswith('.vcf.gz') or fn.endswith('.vcf.bgz'):
                    file_list.append(os.path.join(p, fn))
        else:
            file_list.append(p)

    # validate
    file_list = [f for f in file_list if os.path.exists(f)]
    if not file_list:
        print('No input VCF files found or paths do not exist.', file=sys.stderr)
        sys.exit(1)

    # aggregate across files; per-file presence counts as 1
    records = {}
    for f in file_list:
        per_file = process_vcf(f)
        for key, cnt in per_file.items():
            if key not in records:
                records[key] = {'DUP': 0, 'DEL': 0, 'DUP_CN': [], 'DEL_CN': []}
            
            # Update presence counts
            records[key]['DUP'] += cnt.get('DUP', 0)
            records[key]['DEL'] += cnt.get('DEL', 0)
            
            # Update CN values
            records[key]['DUP_CN'].extend(cnt.get('DUP_CN', []))
            records[key]['DEL_CN'].extend(cnt.get('DEL_CN', []))

    # total number of vcf files = number of samples (assumption)
    total_files = len(file_list)
    df = records_to_df(records, total_files=total_files)

    # decide output naming
    if args.sample_prefix:
        # use the requested sample prefix naming
        prefix = str(args.sample_prefix)
        base_name = f'local_freq.{prefix}.cnv.txt'
        out_unsorted = base_name + '.unsorted'
        # write unsorted first, then we'll sort into final basename
        df.to_csv(out_unsorted, sep='\t', index=False, float_format='%.9f')
        print(f'Processed {len(file_list)} files. Wrote unsorted temporary file to {out_unsorted}')
        out = out_unsorted
    else:
        if args.output:
            out = args.output
        else:
            # default name when multiple inputs: local_freq.result.txt
            if len(file_list) == 1:
                out = os.path.splitext(os.path.basename(file_list[0]))[0] + '.result.txt'
            else:
                out = 'local_freq.result.txt'

        df.to_csv(out, sep='\t', index=False, float_format='%.9f')
        print(f'Processed {len(file_list)} files. Wrote {len(df)} records to {out}')

    # 如果服务器有 bgzip/tabix，可在此处对输出文件进行排序、压缩并建立索引
    try:
        if args.sample_prefix:
            # when sample prefix provided, out is the unsorted temp file; final sorted filename per user requirement:
            prefix = str(args.sample_prefix)
            final_name = f'local_freq.{prefix}.cnv.txt'
            out_sorted = final_name
            source_file = out  # unsorted temp
        else:
            base, ext = os.path.splitext(out)
            out_sorted = f"{base}.sorted{ext}"
            source_file = out

        # 使用 head/tail/sort 来保留表头并按 chrom,start 数值排序
        sort_cmd = (
            f"head -n 1 {shlex.quote(source_file)} > {shlex.quote(out_sorted)}; "
            f"tail -n +2 {shlex.quote(source_file)} | sort -k1,1 -k2,2n >> {shlex.quote(out_sorted)}"
        )
        print('Sorting output with command:', sort_cmd)
        subprocess.run(sort_cmd, shell=True, check=True)

        # 压缩：生成 out_sorted + '.gz'（例如 local_freq.120.cnv.txt.gz）
        gz_path = out_sorted + '.gz' if not out_sorted.endswith('.gz') else out_sorted
        bgzip_cmd = f"bgzip -c {shlex.quote(out_sorted)} > {shlex.quote(gz_path)}"
        print('Compressing with bgzip...')
        subprocess.run(bgzip_cmd, shell=True, check=True)

        # 建立 tabix 索引，chrom=col1 start=col2 end=col3
        # 使用 -f 强制覆盖、-S 1 跳过 1 行 header（因为文件有表头）
        print('Indexing with tabix (skip 1 header line)...')
        subprocess.run(['tabix', '-s', '1', '-b', '2', '-e', '3', '-S', '1', gz_path], check=True)

        # remove temporary unsorted file when we created one for sample-prefix flow
        if args.sample_prefix:
            try:
                os.remove(source_file)
            except Exception:
                pass

        print(f'Created compressed + indexed file: {gz_path} and {gz_path}.tbi')
    except subprocess.CalledProcessError as e:
        print('Error running external command:', e, file=sys.stderr)
    except Exception as e:
        print('Unexpected error during compress/index step:', e, file=sys.stderr)


if __name__ == '__main__':
    main()
