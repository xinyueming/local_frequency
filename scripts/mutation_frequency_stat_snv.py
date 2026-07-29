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
            js = re.match('(.*hg38_multianno.filter.*)',pathfile)
            if js != None:
                R1.append(path)
    return R1

def read_info(tablefile):
    code_dict = dict()
    df = pd.read_excel(tablefile)
    for ki in df.iloc:
        if ki.get('样本编号'):
            code_dict[ki['样本编号']] = ki['项目名称']
    return code_dict

def common_tb(tablefile):
    if pd.__version__ >= '1.3.0':
        tb = pd.read_table(tablefile,encoding="utf-8",encoding_errors="ignore",index_col=False)
    else:
        tb = pd.read_table(tablefile,encoding="utf-8",index_col=False)
    tbot = tb.to_dict(orient="records")
    return tbot

def exc_main(s,t,g=''):
    li = []
    di = {}
    taf = {}
    header = ['Chr',
              'Pos',
              'End_pos',
              'Ref',
              'Alt',
              'Gene.refGeneWithVer',
              'Func.refGeneWithVer',
              'ExonicFunc.refGeneWithVer',
              'AAChange.refGeneWithVer',
              'population',
              'mutation_frequency',
              'TAF']
    f_out = open(f'{t}/mutation_frequency.xls', 'w')
    f_out.write('\t'.join(header)+'\n')
    files = read_R1file(s)
    if not g:
        files = [x for x in files if 'germline' not in x]
    else:
        files = [x for x in files if 'germline' in x]
    for fi in files:
        #file_name = re.split('[/]+',fi)[-1].split('.')[0]
        content = common_tb(fi)
        for line in content:
            ls1 = [line['Chr'],str(line['Pos']),str(line['End_pos']),line['Ref'],line['Alt']]
            ls2 = [line['Gene.refGeneWithVer'],line['Func.refGeneWithVer'],line['ExonicFunc.refGeneWithVer']]
            aachange = line['AAChange.refGeneWithVer']
            info = aachange
            try:
                infos = info.split(':')
            except:
                print(info)
                print(fi)
                files.remove(fi)
                break
            c_info = '.'
            try:
                c_info = infos[3]
            except:
                c_info = infos[-1]
            if 'c.' in c_info:
                if '>' in c_info:
                    pass
                else:
                    tem_new_cdna = c_info.lstrip('c.')
                    a_new_cdna = re.findall('[A-Z]+',tem_new_cdna)
                    n_new_cdna = re.findall('\d+',tem_new_cdna)
                    if len(a_new_cdna) == 2:
                        try:
                            if len(n_new_cdna) == 1 and len(a_new_cdna[0]) == 1 and len(a_new_cdna[1]) == 1:
                                c_info = f'c.{n_new_cdna[0]}{a_new_cdna[0]}>{a_new_cdna[1]}'
                        except:
                            pass
                try:
                    infos[3] = c_info
                except:
                    infos[-1] = c_info
            # aaa = infos[-1]
            # if 'p.' in aaa and 'fs' in aaa:
            #     aaa = aaa.replace('X', '*')
            #     if 'fs*' in aaa:
            #         aaas = aaa.split('fs*')
            #         if len(aaas) > 1:
            #             aaas[-1] = str(int(aaas[-1]) + 1)
            #             aaa = 'fs*'.join(aaas)
            # infos[-1] = aaa
            # info = ':'.join(infos)
            if 'fs*' in aachange:
                aachange = aachange.split('*')[0]
            ls1.extend(ls2)
            #he = [lines[0],lines[1],lines[2],lines[3],lines[4],lines[6],lines[7],lines[8]]
            if 'c.' not in aachange:
                aachange = f'{line["new_gene"]}:{line["new_transcript"]}:{line["new_exon"]}:{line["new_cdna"]}'
            if 'c.' not in aachange or '\\x3b' in line['Gene.refGeneWithVer']:
                chrom = line['Chr'].strip()
                pos1 = line['Pos']
                ref = line['Ref']
                alt = line['Alt']
                aachange = f'{chrom}_{pos1}_{ref}_{alt}'
            di[aachange] = ls1
            li.append(aachange)
            taff = ''
            taff = str(line['TAF'])
            if aachange in taf:
                taf[aachange].append(taff)
            else:
                taf[aachange] = []
                taf[aachange].append(taff)
    counter = Counter(li).most_common()
    total = len(files)
    for k,v in counter:
        gene_name = k.split(':')[0]
        pop = int(v)
        mut_freq = round(pop/total,4)
        f_out.write('\t'.join(di[k])+'\t'+k+'\t'+str(pop)+'\t'+str(mut_freq)+'\t'+','.join(taf[k])+'\n')
    f_out.close()
    return total

if __name__=='__main__':
    n = 0
    while n <1:
        outdir = f'{script_dir}/mutation_frequency_result'
        t_84_snv = f'{script_dir}/SNV_filter/84'
        t_624_snv = f'{script_dir}/SNV_filter/624'
        t_245_snv = f'{script_dir}/SNV_filter/245'
        t_606_snv = f'{script_dir}/SNV_filter/606'
        t_cml206_snv = f'{script_dir}/SNV_filter/cml206'
        t_kywes_snv = f'{script_dir}/SNV_filter/kywes'
        t_rnaseq_snv = f'{script_dir}/SNV_filter/RNASeq'
        g_84_snv = f'{script_dir}/SNV_germline_filter/84'
        g_624_snv = f'{script_dir}/SNV_germline_filter/624'
        temp_path = f'{script_dir}/SNV_filter_temp'
        now_time = time.strftime("%y%m%d",time.localtime())
        os.system(f'mkdir -p {outdir}/{str(now_time)}')
        if os.path.exists(f'{outdir}/latest'):
            pass
        else:
            os.system(f'mkdir -p {outdir}/latest')
        total_files = read_R1file(temp_path)
        for fl in total_files:
            if '624/' in fl:
                if 'germline' not in fl:
                    os.system(f'mv {fl} {t_624_snv}')
                else:
                    os.system(f'mv {fl} {g_624_snv}')
            elif '84/' in fl:
                if 'germline' not in fl:
                    os.system(f'mv {fl} {t_84_snv}')
                else:
                    os.system(f'mv {fl} {g_84_snv}')
            elif '245/' in fl:
                if 'germline' not in fl:
                    os.system(f'mv {fl} {t_245_snv}')
            elif '606/' in fl:
                if 'germline' not in fl:
                    os.system(f'mv {fl} {t_606_snv}')
            elif 'cml206/' in fl:
                if 'germline' not in fl:
                    os.system(f'mv {fl} {t_cml206_snv}')
            elif 'kywes/' in fl:
                if 'germline' not in fl:
                    os.system(f'mv {fl} {t_kywes_snv}')
            elif 'RNASeq/' in fl:
                if 'germline' not in fl:
                    os.system(f'mv {fl} {t_rnaseq_snv}')

        total84 = exc_main(t_84_snv, t_84_snv)
        total624 = exc_main(t_624_snv, t_624_snv)
        total245 = exc_main(t_245_snv, t_245_snv)
        total606 = exc_main(t_606_snv, t_606_snv)
        total206 = exc_main(t_cml206_snv, t_cml206_snv)
        total_kywes = exc_main(t_kywes_snv, t_kywes_snv)
        total_rnaseq = exc_main(t_rnaseq_snv, t_rnaseq_snv)
        total84_germline = exc_main(g_84_snv, g_84_snv, 'germline')
        total624_germline = exc_main(g_624_snv, g_624_snv, 'germline')
        total_json = {'84': total84,
                      '624': total624,
                      'rnapanel245': total245,
                      'rnapanel606': total606,
                      'cml206': total206,
                      'kywes': total_kywes,
                      'rnaseq': total_rnaseq,
                      '84_germline': total84_germline,
                      '624_germline': total624_germline,
                      }
        json.dump(total_json, open("stat.json", "w"), indent=4)
        os.system(f'cp {script_dir}/*.xls {outdir}/latest/')
        os.system(f'cp {script_dir}/stat.json {outdir}/latest/')
        os.system(f'cp {script_dir}/*.xls {outdir}/{str(now_time)}/')
        os.system(f'cp {script_dir}/stat.json {outdir}/{str(now_time)}/')
        #time.sleep(2592000)
        n += 1
