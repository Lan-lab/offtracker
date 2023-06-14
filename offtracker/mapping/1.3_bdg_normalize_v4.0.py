import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.description='GG!'

# 2022.10.26. v4.0:	整合了 scale 和翻转 bdg_rv 的功能

# 单文件处理脚本，配合snakemake使用

parser.add_argument("--bdg", type=str, metavar="dir_bdg" , required=True, help="dir of bdg file")
parser.add_argument("--bed", type=str, metavar="dir_bed" , required=True, help="dir of bed file")

args = parser.parse_args()

# total reads
total_reads = sum(1 for _ in open(args.bed))
scale_factor = 1e7/total_reads
if args.bdg[-13:-11] == 'rv':
	scale_factor = -scale_factor

# normalize
bdg_file = pd.read_csv( args.bdg, sep='\t', header=None )
bdg_file[3] = bdg_file[3]*scale_factor
bdg_file.to_csv( args.bdg[:-11]+'.scaled.bdg', sep='\t', header=False, index=False )




