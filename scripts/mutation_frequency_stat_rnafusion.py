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
            js = re.match('(.*rna.starfusion.tsv.redup.xls)',pathfile)
            if js != None:
                R1.append(path)
    return R1

def common_tb(tablefile):
    if pd.__version__ >= "1.3.0":
        tb = pd.read_table(tablefile,encoding="utf-8",encoding_errors="ignore",index_col=False)
    else:
        tb = pd.read_table(tablefile,encoding="utf-8",index_col=False)
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
            ls1 = [str(line['gene1_chr']),str(line['gene1_pos']),line['gene1'],str(line['gene2_chr']),str(line['gene2_pos']),line['gene2'],str(line['exon'].split('-')[0]),str(line['exon'].split('-')[1])]
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
        t_rnapanel_fusion_245 = f'{script_dir}/fusion_filter/rnapanel/245'
        t_rnapanel_fusion_606 = f'{script_dir}/fusion_filter/rnapanel/606'
        t_rnaseq_fusion = f'{script_dir}/fusion_filter/rnaseq'
        temp_path = f'{script_dir}/fusion_filter_temp'
        now_time = time.strftime("%y%m%d",time.localtime())
        os.system(f'mkdir -p {outdir}/{str(now_time)}')
        if os.path.exists(f'{outdir}/latest'):
            pass
        else:
            os.system(f'mkdir -p {outdir}/latest')
        total_files = read_R1file(temp_path)
        for fl in total_files:
            if 'rnapanel/245' in fl:
                os.system(f'mv {fl} {t_rnapanel_fusion_245}')
            elif 'rnapanel/606' in fl:
                os.system(f'mv {fl} {t_rnapanel_fusion_606}')
            elif 'rnaseq/' in fl:
                os.system(f'mv {fl} {t_rnaseq_fusion}')

        total_rnapanel_245 = exc_main(t_rnapanel_fusion_245, t_rnapanel_fusion_245)
        total_rnapanel_606 = exc_main(t_rnapanel_fusion_606, t_rnapanel_fusion_606)
        total_rnaseq = exc_main(t_rnaseq_fusion, t_rnaseq_fusion)
        total_json = {
                      'rnapanelfusion245': total_rnapanel_245,
                      'rnapanelfusion606': total_rnapanel_606,
                      'rnaseqfusion': total_rnaseq,
                      }
        json.dump(total_json, open("stat.json", "w"), indent=4)
        os.system(f'cp {script_dir}/*.xls {outdir}/latest/')
        os.system(f'cp {script_dir}/stat.json {outdir}/latest/')
        os.system(f'cp {script_dir}/*.xls {outdir}/{str(now_time)}/')
        os.system(f'cp {script_dir}/stat.json {outdir}/{str(now_time)}/')
        #time.sleep(2592000)
        n += 1
