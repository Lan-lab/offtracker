#!/usr/bin/env python
# -*- coding: utf-8 -*-

THIS_VERSION = '0.4.2'

# 2025.10.05. 0.4.2. 添加 threads 监测

import argparse
import os, glob, yaml
import pandas as pd
import shutil, re
import offtracker
import offtracker.X_sequence as xseq

script_dir = os.path.abspath(os.path.dirname(offtracker.__file__))
utility_dir = os.path.join(script_dir, 'utility')
file_path = os.path.join(utility_dir, 'bedGraphToBigWig')
file_stat = os.stat(file_path)
file_mode = oct(file_stat.st_mode & 0o777)
if file_mode != '0o755':
    try:
        os.chmod( os.path.join(utility_dir, 'bedGraphToBigWig'), 0o755)
    except:
        print('offtracker may be installed in root but not initialized. Please run "offtracker_init.py" with root permission first.')

###
def main():
    parser = argparse.ArgumentParser()
    parser.description=f'xbulk_qc v{THIS_VERSION}. QC and trim fastq files.'
    parser.add_argument('-f','--folder', type=str, required=True,        help='Directory of the input folder' )
    parser.add_argument('-o','--outdir', type=str, default='same',       help='The output folder')
    parser.add_argument('--subfolder'  , type=int, default=0,            help='subfolder level')
    parser.add_argument('-t','--thread', type=int, default=8,            help='Number of threads to be used')
    parser.add_argument('--include_trimmed', action='store_true',        help='Do not skip trimmed fastq files')
    
    args = parser.parse_args()

    # 自动化的参数调整和报错
    if args.outdir == 'same':
        args.outdir = os.path.join(args.folder,'Trimmed_data')
        if not os.path.exists( args.outdir ):
            os.makedirs( args.outdir )
    else:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

    # 搜索 folder 的 n级子目录下的所有 fastq/fastq.gz/fq/fq.gz 文件
    if args.include_trimmed:
        sample_names, files_R1, files_R2 = xseq.detect_fastq(args.folder, n_subfolder=args.subfolder)
    else:
        sample_names, files_R1, files_R2 = xseq.detect_fastq(args.folder, n_subfolder=args.subfolder, skip_trimmed=True)

    assert not isinstance(sample_names, str), 'No fastq file is detected!'

    ################
    # threads 监测 #
    ################
    import psutil
    n_threads = args.thread
    assert n_threads > 0, f'n_threads should be greater than 0, while {n_threads} is given.'
    cpu_count_total = psutil.cpu_count(logical=True)  # 逻辑 CPU 总数（包括超线程）
    if n_threads > cpu_count_total:
        n_threads = cpu_count_total-1
        print(f'n_threads is reset to {n_threads} due to the total number of threads ({cpu_count_total}).')

    dict_yaml = {
        # fastq 信息
        'files_R1':dict(zip(sample_names,files_R1)),
        'files_R2':dict(zip(sample_names,files_R2)), # 单端 files_R2=[] 结果会自动为 {}
        # 输入输出文件夹
        'input_dir':args.folder,
        'output_dir':args.outdir,
        # 运行参数
        'thread':n_threads,
        'utility_dir':utility_dir
        }


    with open( os.path.join(args.outdir,'config.yaml'), 'w', encoding='utf-8') as outfile:
        yaml.dump(dict_yaml, outfile, default_flow_style=False)

    snakefile = os.path.join(script_dir, 'snakefile/Snakefile_QC.smk')
    shutil.copy(snakefile, os.path.join(args.outdir,'Snakefile'))

    return 'config_qc finished'

if __name__ == '__main__' :
    result = main()
    print(result)

