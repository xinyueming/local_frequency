import time
import sys
import os
import re
import json
import pandas as pd
from collections import Counter

# 获取当前脚本的绝对路径
script_path = os.path.abspath(__file__)

# 获取脚本所在的目录
script_dir = os.path.dirname(script_path)

def read_R1file(sourpath):
    l = os.listdir(sourpath)
    R1 = [] #用于存放文件路径的list
    for i in l:
        path = os.path.join(sourpath,i)
        if os.path.isdir(path):
            R1.extend(read_R1file(path))
        if os.path.isfile(path):
            pathfile = re.split('[/]+',path)[-1]
            js = re.match('(.*redup.xls)',pathfile)
            js2 = re.match('(.*fusion.xls)',pathfile)
            if js != None or js2 != None:
                R1.append(path)
    return R1

def common_tb(tablefile):
    tb = pd.read_table(tablefile,encoding="utf-8",encoding_errors="ignore",index_col=False)
    tbot = tb.to_dict(orient="records")
    return tbot

def exc_main(s,t):
    li = []
    di = {}
    taf = {}
    header = ['gene1_chr',
              'gene1_pos',
              'gene1',
              'gene2_chr',
              'gene2_pos',
              'gene2',
              'exon1',
              'exon2',
              'FusionChange',
              'population',
              'mutation_frequency',
              'TAF']
    f_out = open(f'{t}/mutation_frequency.xls', 'w')
    f_out.write('\t'.join(header)+'\n')
    files = read_R1file(s)
    for fi in files:
        content = common_tb(fi)
        temp_list = []
        for line in content:
            ls1 = [str(line['gene1_chr']),str(line['gene1_pos']),line['gene1'],str(line['gene2_chr']),str(line['gene2_pos']),line['gene2'],str(re.split('[-:]+',line['exon'])[0]),str(re.split('[-:]+',line['exon'])[1])]
            fusionchange = f'{line["gene1_chr"]}:{line["gene1_pos"]}_{line["gene2_chr"]}:{line["gene2_pos"]}'
            if fusionchange in temp_list:
                continue
            temp_list.append(fusionchange)
            di[fusionchange] = ls1
            li.append(fusionchange)
            taff = ''
            taff = str(line['freq'])
            if fusionchange in taf:
                taf[fusionchange].append(taff)
            else:
                taf[fusionchange] = []
                taf[fusionchange].append(taff)
    counter = Counter(li).most_common()
    total = len(files)
    for k,v in counter:
        pop = int(v)
        mut_freq = round(pop/total,4)
        f_out.write('\t'.join(di[k])+'\t'+k+'\t'+str(pop)+'\t'+str(mut_freq)+'\t'+','.join(taf[k])+'\n')
    f_out.close()
    return total

if __name__=='__main__':
    n = 0
    while n <1:
        outdir = f'{script_dir}/mutation_frequency_result'
        t_rnapanel_fusion = f'{script_dir}/fusion_filter/rnapanel'
        if not os.path.exists(t_rnapanel_fusion):
            os.system(f'mkdir -p {t_rnapanel_fusion}')
        t_rnaseq_fusion = f'{script_dir}/fusion_filter/rnaseq'
        if not os.path.exists(t_rnaseq_fusion):
            os.system(f'mkdir -p {t_rnaseq_fusion}')
        t_120_fusion = f'{script_dir}/fusion_filter/120'
        if not os.path.exists(t_120_fusion):
            os.system(f'mkdir -p {t_120_fusion}')
        t_1100_fusion = f'{script_dir}/fusion_filter/1100'
        if not os.path.exists(t_1100_fusion):
            os.system(f'mkdir -p {t_1100_fusion}')
        t_wes_fusion = f'{script_dir}/fusion_filter/wes'
        if not os.path.exists(t_wes_fusion):
            os.system(f'mkdir -p {t_wes_fusion}')
        t_84_fusion = f'{script_dir}/fusion_filter/84'
        if not os.path.exists(t_84_fusion):
            os.system(f'mkdir -p {t_84_fusion}')
        t_624_fusion = f'{script_dir}/fusion_filter/624'
        if not os.path.exists(t_624_fusion):
            os.system(f'mkdir -p {t_624_fusion}')
        temp_path = f'{script_dir}/fusion_filter_temp'
        now_time = time.strftime("%y%m%d",time.localtime())
        os.system(f'mkdir -p {outdir}/{str(now_time)}')
        if os.path.exists(f'{outdir}/latest'):
            pass
        else:
            os.system(f'mkdir -p {outdir}/latest')
        total_files = read_R1file(temp_path)
        for fl in total_files:
            if 'rnapanel/' in fl:
                os.system(f'mv {fl} {t_rnapanel_fusion}')
            elif 'rnaseq/' in fl:
                os.system(f'mv {fl} {t_rnaseq_fusion}')
            elif '120/' in fl:
                os.system(f'mv {fl} {t_120_fusion}')
            elif '1100/' in fl:
                os.system(f'mv {fl} {t_1100_fusion}')
            elif 'wes/' in fl:
                os.system(f'mv {fl} {t_wes_fusion}')
            elif '84/' in fl:
                os.system(f'mv {fl} {t_84_fusion}')
            elif '624/' in fl:
                os.system(f'mv {fl} {t_624_fusion}')

        total_rnapanel = exc_main(t_rnapanel_fusion, t_rnapanel_fusion)
        total_rnaseq = exc_main(t_rnaseq_fusion, t_rnaseq_fusion)
        total_120 = exc_main(t_120_fusion, t_120_fusion)
        total_1100 = exc_main(t_1100_fusion, t_1100_fusion)
        total_wes = exc_main(t_wes_fusion, t_wes_fusion)
        total_84 = exc_main(t_84_fusion, t_84_fusion)
        total_624 = exc_main(t_624_fusion, t_624_fusion)
        total_json = {
                      'rnapanelfusion': total_rnapanel,
                      'rnaseqfusion': total_rnaseq,
                      '120fusion': total_120,
                      '1100fusion': total_1100,
                      'wesfusion': total_wes,
                      '84fusion': total_84,
                      '624fusion': total_624
                      }
        json.dump(total_json, open("stat.json", "w"), indent=4)
        os.system(f'cp {script_dir}/*.xls {outdir}/latest/')
        os.system(f'cp {script_dir}/stat.json {outdir}/latest/')
        os.system(f'cp {script_dir}/*.xls {outdir}/{str(now_time)}/')
        os.system(f'cp {script_dir}/stat.json {outdir}/{str(now_time)}/')
        #time.sleep(2592000)
        n += 1
