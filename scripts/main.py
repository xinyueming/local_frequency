""" 本地频率数据库注释处理程序 """

import os
import json
import numpy as np

import argparse
import logging
import multiprocessing


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
def init_logger():
    # 配置子进程的日志记录器
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    

class LocalAfProcessor:
    """本地频率数据库注释处理程序"""

    def __init__(self, local_af_file, history_sample_count):
        
        # 本次注释的历史样本数
        self.history_sample_count = history_sample_count
        
        # 需要转换的 xls 文件
        self.local_af_file = local_af_file
        self.input_dir = os.path.dirname(self.local_af_file)
        self.version = os.path.basename(self.input_dir)
        
        self.local_af_file_name = os.path.basename(self.local_af_file)

        # 工作目录
        self.work_dir = os.path.dirname(self.local_af_file)

        # 中间结果目录
        self.tmp_dir = os.path.join(self.work_dir, "tmp")

        # 输出结果目录
        self.output_dir = os.path.join(self.work_dir, "output")

        # 中间文件-vcf文件
        self.vcf_file = os.path.join(
            self.tmp_dir,
            self.local_af_file_name.replace(
                ".xls", ".vcf"
            ),
        )
        # 中间文件-avinput文件
        self.avinput = self.vcf_file.replace("vcf", "avinput")
        # 中间文件-raw.txt文件
        self.rawtxtfile = self.vcf_file.replace("vcf", "raw.txt")
        # 中间文件-txt文件
        tmp_txtfile = self.vcf_file.replace("vcf", "txt")
        # 中间文件-txt文件-文件名
        tmp_txtfile_name = os.path.basename(tmp_txtfile)

        # 输出文件-txt文件
        self.txtfile = os.path.join(self.output_dir, tmp_txtfile_name)

        # 检测文件是否存在
        self._check_file(self.local_af_file)
        self._check_dir(self.vcf_file)
        self._check_dir(self.txtfile)

    def _check_file(self, file):
        """检查文件是否存在"""
        if not os.path.exists(file):
            raise FileNotFoundError(f"{file} not found")

    def _check_dir(self, file):
        """检查输出目录是否存在"""
        dirname = os.path.dirname(file)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    def tr_file_to_vcf(self):
        """把本地频率文件转为vcf文件"""
        with open(self.vcf_file, "w", encoding="utf-8") as fhO:
            fhO.write(
                "\t".join(
                    ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                )
                + "\n"
            )
            with open(self.local_af_file, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    # 移除标题行
                    if line.startswith("Chr") or line.startswith("gene"):
                        continue
                    if line:
                        arr = line.split("\t")
                        chrom, pos, _, ref, alt = arr[:5]
                        aachange = arr[8]
                        # 出现概率
                        local_frequency = float(arr[10])
                        # 本地频率数据（百分比形式）
                        local_af_list_str = arr[11]
                        r = []
                        for value in local_af_list_str.split(","):
                            r.append(float(value.strip("%")) / 100.0)
                        data = np.array(r)
                        # 计算位置值
                        pct_local_af = {
                            "LOCAL_FREQUENCY": local_frequency,
                            "LOCAL_AF_AVG": np.mean(data),
                            "LOCAL_AF_PCT0": np.min(data),
                            "LOCAL_AF_PCT25": np.percentile(data, 25),
                            "LOCAL_AF_PCT50": np.percentile(data, 50),
                            "LOCAL_AF_PCT75": np.percentile(data, 75),
                            "LOCAL_AF_PCT100": np.max(data),
                            "LOCAL_AN": len(data),
                            "LOCAL_AF_VERSION": f'{self.version}({self.history_sample_count})',
                        }
                        new_line = "\t".join(
                            [
                                chrom,
                                pos,
                                ".",
                                ref,
                                alt,
                                ".",
                                ".",
                                f"LOCAL_FREQUENCY={pct_local_af['LOCAL_FREQUENCY']:.4f};\
                                LOCAL_AF_AVG={pct_local_af['LOCAL_AF_AVG']:.4f};\
                                LOCAL_AF_PCT0={pct_local_af['LOCAL_AF_PCT0']:.4f};\
                                LOCAL_AF_PCT25={pct_local_af['LOCAL_AF_PCT25']:.4f};\
                                LOCAL_AF_PCT50={pct_local_af['LOCAL_AF_PCT50']:.4f};\
                                LOCAL_AF_PCT75={pct_local_af['LOCAL_AF_PCT75']:.4f};\
                                LOCAL_AF_PCT100={pct_local_af['LOCAL_AF_PCT100']:.4f};\
                                LOCAL_AN={pct_local_af['LOCAL_AN']};\
                                LOCAL_AF_VERSION={pct_local_af['LOCAL_AF_VERSION']};\
                                AAChange={aachange}"
                                + "\n",
                            ]
                        )
                        fhO.write(new_line)

    def run_cmd(self, cmd):
        """运行命令"""
        logging.info(cmd)
        os.system(cmd)

    def tr_vcf_to_avinput(self):
        """把vcf文件转为avinput文件"""
        cmd = f"perl script/convert2annovar.pl --format vcf4 --includeinfo {self.vcf_file} > {self.avinput}"
        self.run_cmd(cmd)

    def tr_avinput_to_txt(self):
        """把avinput文件转为txt文件"""
        with open(self.rawtxtfile, "w", encoding="utf-8") as fhO:
            fhO.write(
                "\t".join(
                    [
                        "#Chr",
                        "Start",
                        "AAChange",
                        "Ref",
                        "Alt",
                        "LOCAL_FREQUENCY",
                        "LOCAL_AF_AVG",
                        "LOCAL_AF_PCT0",
                        "LOCAL_AF_PCT25",
                        "LOCAL_AF_PCT50",
                        "LOCAL_AF_PCT75",
                        "LOCAL_AF_PCT100",
                        "LOCAL_AN",
                        "LOCAL_AF_VERSION",
                    ]
                )
                + "\n"
            )

            with open(self.vcf_file, encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if line:
                        arr = line.split("\t")
                        chrom, start, end, ref, alt = arr[:5]
                        info = arr[7]

                        info_dict = {}
                        if info == "INFO":
                            continue
                        for i in info.split(";"):
                            if "=" in i:
                                key, value = i.split("=")
                                key = key.strip()
                                value = value.strip()
                                info_dict[key] = value
                        local_frequency = info_dict["LOCAL_FREQUENCY"]
                        local_af_avg = info_dict["LOCAL_AF_AVG"]
                        local_af_pct0 = info_dict["LOCAL_AF_PCT0"]
                        local_af_pct25 = info_dict["LOCAL_AF_PCT25"]
                        local_af_pct50 = info_dict["LOCAL_AF_PCT50"]
                        local_af_pct75 = info_dict["LOCAL_AF_PCT75"]
                        local_af_pct100 = info_dict["LOCAL_AF_PCT100"]
                        local_an = info_dict["LOCAL_AN"]
                        local_af_version = info_dict["LOCAL_AF_VERSION"]
                        aachange = info_dict["AAChange"]
                        
                        new_line = "\t".join(
                            [
                                chrom,
                                start,
                                aachange,
                                ref,
                                alt,
                                local_frequency,
                                local_af_avg,
                                local_af_pct0,
                                local_af_pct25,
                                local_af_pct50,
                                local_af_pct75,
                                local_af_pct100,
                                local_an,
                                local_af_version,
                            ]
                        )

                        fhO.write(new_line + "\n")

    def index_txt(self):
        """索引txt文件"""
        cmd = f"perl script/index_annovar.pl {self.rawtxtfile} --outfile {self.txtfile}"
        self.run_cmd(cmd)

    def process(self):
        """处理"""
        self.tr_file_to_vcf()
        self.tr_vcf_to_avinput()
        self.tr_avinput_to_txt()
        self.index_txt()

def one_process(local_af_file, history_sample_count):
    logger = multiprocessing.get_logger()
    logger.info(f'start local af file: {local_af_file}')
    LocalAfProcessor(local_af_file=local_af_file, history_sample_count=history_sample_count).process()
    logger.info(f'finish local af file: {local_af_file}')

def mul_main():
    # 多进程
    parser = argparse.ArgumentParser(description="转换本地频率文件为 annovar db 的注释格式")
    parser.add_argument("-d", "--xls_dir", required=True, type=str, help="本地频率文件的目录")
    parser.add_argument("-t", "--thread", default=4, type=str, help="本地频率文件的目录")
    
    args = parser.parse_args()
    local_af_xls_dir = args.xls_dir
    thread_size = args.thread

    if not os.path.exists(local_af_xls_dir):
        logging.info(f"本地频率结果目录不存在：{local_af_xls_dir}")
        return

    # 创建进程池，指定最大进程数
    with multiprocessing.Pool(processes=thread_size, initializer=init_logger) as pool:
        for local_af_file_name in os.listdir(local_af_xls_dir):
            if local_af_file_name.endswith(".xls"):
                # 使用进程池执行任务
                local_af_file = os.path.join(local_af_xls_dir, local_af_file_name)
                logging.info(f'add one task: {local_af_file}')
                pool.apply_async(one_process, args=(local_af_file,))

    # 关闭进程池
    pool.close()

    # 等待所有进程完成
    pool.join()   
    
    logging.info('all task has been finished')
    
    
def main():
    parser = argparse.ArgumentParser(description="转换本地频率文件为 annovar db 的注释格式")
    parser.add_argument("-d", "--xls_dir", required=True, type=str, help="本地频率文件的目录")
    
    args = parser.parse_args()
    local_af_xls_dir = args.xls_dir

    if not os.path.exists(local_af_xls_dir):
        logging.info(f"本地频率结果目录不存在：{local_af_xls_dir}")
        return

    stat_file = os.path.join(local_af_xls_dir, "stat.json")
    if not os.path.exists(stat_file):
        logging.info(f"本地频率结果统计文件不存在：{stat_file}")
        return
    
    with open(stat_file, encoding="utf-8") as fh:
        stat = json.load(fh)

    mapping = {
        '120.mutation_frequency.xls': '120',
        '120.mutation_frequency_germline.xls': '120_germline',
        '180.mutation_frequency.xls': '180',
        '180.mutation_frequency_germline.xls': '180_germline',
        '680.mutation_frequency.xls': '680',
        '680.mutation_frequency_germline.xls': '680_germline',
        '1100.mutation_frequency.xls': '1100', 
        '1100.mutation_frequency_germline.xls': '1100_germline',
        'WES.mutation_frequency.xls': 'wes',
        'WES.mutation_frequency_germline.xls': 'wes_germline',
        '84.mutation_frequency.xls': '84',
        '84.mutation_frequency_germline.xls': '84_germline',
        '624.mutation_frequency.xls': '624',
        '624.mutation_frequency_germline.xls': '624_germline',
        'cml206.mutation_frequency.xls': 'cml206',
        'kywes.mutation_frequency.xls': 'kywes',
        'rnaseq.mutation_frequency.xls': 'rnaseq',
        'rnapanel.245.mutation_frequency.xls': 'rnapanel245',
        'rnapanel.606.mutation_frequency.xls': 'rnapanel606',
        'rnaseqfusion.mutation_frequency.xls': 'rnaseqfusion',
        'rnapanelfusion.245.mutation_frequency.xls': 'rnapanelfusion245',
        'rnapanelfusion.606.mutation_frequency.xls': 'rnapanelfusion606',
        '120fusion.mutation_frequency.xls': '120fusion',
        '1100fusion.mutation_frequency.xls': '1100fusion',
        'wesfusion.mutation_frequency.xls': 'wesfusion',
        '84fusion.mutation_frequency.xls': '84fusion',
        '624fusion.mutation_frequency.xls': '624fusion',
    }

    for local_af_file_name in os.listdir(local_af_xls_dir):
        if local_af_file_name.endswith(".xls"):
            
            key = mapping.get(local_af_file_name)
            if not key:
                raise ValueError(f'未能获取到统计样本数的文件名：{local_af_file_name}')
            
            history_sample_count = stat.get(key)
            if not history_sample_count:
                raise ValueError(f'stat.json 文件中，无法相关统计{key}')
            
            local_af_file = os.path.join(local_af_xls_dir, local_af_file_name)
            logging.info(f'start task {local_af_file}')
            one_process(local_af_file, history_sample_count)
    
    logging.info('all task has been finished')

if __name__ == "__main__":
    main()
