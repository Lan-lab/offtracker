import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.description='这算一个小彩蛋'

# 2022.10.21. v3.0: 文件名长度 chromap -> filtered
# 2022.10.26. v4.0: f,r 改成 fw,rv
# 2022.01.11. v4.5: 只取 common chromosomes (chr1-chr22, chrX, chrY, chrM)

# 单文件处理脚本，配合snakemake使用

parser.add_argument("-b", "--bed", type=str, metavar="dir_bed" , required=True, help="dir of bed file")

args = parser.parse_args()

bed_file = pd.read_csv( args.bed, sep='\t', header=None)

common_chr = pd.Series(['chr']*22).str[:] + pd.Series(range(1,23)).astype(str).str[:]
common_chr = pd.concat([common_chr, pd.Series(['chrX','chrY','chrM'])]).to_numpy()

bed_file = bed_file[bed_file[0].isin(common_chr)]

bed_f = bed_file[bed_file[5]=='+']
bed_r = bed_file[bed_file[5]=='-']
bed_f.to_csv(args.bed[:-13] + '.fw.bed',sep='\t',header=False,index=False)
bed_r.to_csv(args.bed[:-13] + '.rv.bed',sep='\t',header=False,index=False)





