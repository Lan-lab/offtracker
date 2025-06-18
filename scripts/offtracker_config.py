#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2023.08.11. adding a option for not normalizing the bw file
# 2025.05.22. refine the structure
# 2025.06.05. 增加 ignore_chr 选项，默认只取 common chromosomes，用于 1.1_bed2fr.py

import argparse
import os, glob, yaml
import pandas as pd
import shutil, re
import offtracker
import offtracker.X_sequence as xseq
script_dir = os.path.abspath(os.path.dirname(offtracker.__file__))
utility_dir = os.path.join(script_dir, 'utility')
os.chmod( os.path.join(utility_dir, 'bedGraphToBigWig'), 0o755)

###
parser = argparse.ArgumentParser()
parser.description='Mapping fastq files of Tracking-seq.'
parser.add_argument('-f','--folder', type=str, required=True,  help='Directory of the input folder' )
parser.add_argument('-r','--ref'   , type=str, required=True,  help='The fasta file of reference genome')
parser.add_argument('-i','--index' , type=str, required=True,  help='The index file of chromap')
parser.add_argument('-g','--genome', type=str, required=True,  help='File of chromosome sizes, or "hg38", "mm10" ')
parser.add_argument('-o','--outdir', type=str, default='same', help='The output folder')
parser.add_argument('--subfolder'  , type=int, default=0,      help='subfolder level')
parser.add_argument('-t','--thread', type=int, default=4,      help='Number of threads to be used')
parser.add_argument('--blacklist'  , type=str, default='same', help='Blacklist of genome regions in bed format. "none" for no filter')
parser.add_argument('--binsize'    , type=str, default=100,    help='Bin size for calculating bw residue')
parser.add_argument('--normalize'  , type=str, default='True', help='Whether to normalize the BigWig file. "True" or "False"')
parser.add_argument('--ignore_chr' , action='store_true', help='If not set, only chr1-chr22, chrX, chrY, chrM will be analyzed.')


args = parser.parse_args()

if (args.genome == 'hg38') or (args.genome == 'mm10'):
    dir_chrom_sizes = os.path.join(utility_dir, f'{args.genome}.chrom.sizes')
else:
    dir_chrom_sizes = args.genome

if (args.normalize != 'True') & (args.normalize != 'False'):
    raise ValueError('Please provide "True" or "False" for "--normalize"')

if args.blacklist == 'same':
    assert ((args.genome == 'hg38') or (args.genome == 'mm10')), 'Please provide blacklist file, or "--blacklist none" to skip'
    args.blacklist = args.genome
    
if (args.blacklist == 'hg38') or (args.blacklist == 'mm10'):
    blacklist = os.path.join(utility_dir, f'offtracker_blacklist_{args.blacklist}.merged.bed')
else:
    blacklist = args.blacklist

if args.outdir == 'same':
    args.outdir = args.folder
else:
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

if args.ignore_chr:
    args.ignore_chr = '--ignore_chr'
else:
    args.ignore_chr = ''

# 搜索 folder 的 n级子目录下的所有 fastq/fastq.gz/fq/fq.gz 文件
sample_names, files_R1, files_R2 = xseq.detect_fastq(args.folder, n_subfolder=args.subfolder)

assert not isinstance(sample_names, str), 'No fastq file is detected!'

dict_yaml = {
    # fastq 信息
    'files_R1':dict(zip(sample_names,files_R1)),
    'files_R2':dict(zip(sample_names,files_R2)), # 单端 files_R2=[] 结果会自动为 {}
    # 输入输出文件夹
    'input_dir':args.folder,
    'output_dir':args.outdir,
    # 运行参数
    'thread':args.thread,
    'index':args.index,
    'fasta':args.ref,
    'binsize':args.binsize,
    'blacklist':blacklist,
    'genomelen':dir_chrom_sizes,
    'normalize':args.normalize,
    'utility_dir':utility_dir,
    'ignore_chr':args.ignore_chr,
    }

with open( os.path.join(args.outdir,'config.yaml'), 'w') as outfile:
    yaml.dump(dict_yaml, outfile, default_flow_style=False)

snakefile = os.path.join(script_dir, 'snakefile/Snakefile_offtracker.smk')
shutil.copy(snakefile, os.path.join(args.outdir,'Snakefile'))


