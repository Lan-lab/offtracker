#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os, glob, yaml
import pandas as pd
import shutil, re
import offtracker
script_dir = os.path.abspath(os.path.dirname(offtracker.__file__))
script_folder= os.path.join(script_dir, 'mapping')
os.chmod( os.path.join(script_folder, 'bedGraphToBigWig'), 0o755)

###
parser = argparse.ArgumentParser()
parser.description='Mapping fastq files of Track-seq.'
parser.add_argument('-f','--folder', type=str, required=True, help='Directory of the input folder' )
parser.add_argument('-r','--ref'   , type=str, required=True, help='The fasta file of reference genome')
parser.add_argument('-i','--index' , type=str, required=True, help='The index file of chromap')
parser.add_argument('-g','--genome', type=str, required=True, help='File of chromosome sizes, or "hg38", "mm10" ')
parser.add_argument('-o','--outdir', type=str, default='same', help='The output folder')
parser.add_argument('--subfolder'  , type=int, default=0,        help='subfolder level')
parser.add_argument('-t','--thread', type=int, default=4,     help='Number of threads to be used')
parser.add_argument('--blacklist'  , type=str, default='same', help='Blacklist of genome regions in bed format.')
parser.add_argument('--binsize'    , type=str, default=10, help='Bin size for calculating bw ratio')

args = parser.parse_args()


if (args.genome == 'hg38') or (args.genome == 'mm10'):
    dir_chrom_sizes = os.path.join(script_folder, f'{args.genome}.chrom.sizes')
else:
    dir_chrom_sizes = args.genome


if args.blacklist == 'same':
    assert ((args.genome == 'hg38') or (args.genome == 'mm10')), 'Please provide blacklist file, or "--blacklist none" to skip'
    args.blacklist = args.genome
    
if (args.blacklist == 'hg38') or (args.blacklist == 'mm10'):
    blacklist = os.path.join(script_folder, f'offtracker_blacklist_{args.blacklist}.merged.bed')
else:
    blacklist = args.blacklist

if args.outdir == 'same':
    args.outdir = args.folder
else:
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

gz_R2 = []
for fastq in ['*2.*fq','*2.*fastq','*2.*fq.gz','*2.*fastq.gz']:
    fq_files = glob.glob( os.path.join(args.folder, args.subfolder*'*/', fastq ) )
    print('{} {} samples detected'.format( len(fq_files), fastq[4:] ) )
    gz_R2.extend( fq_files )

gz_R2.sort()
gz_R2 = pd.Series(gz_R2)
suffix = gz_R2.str.extract('(fastq.*|fq.*)',expand=False)
prefix = gz_R2.str.extract('(.*)(?:.fq|.fastq)',expand=False)

nametype = None
for a_type in ['_trimmed_2', '_2_val_2','_R2_val_2','_R2','_2']:
    len_type = len(a_type)
    if prefix[0][-len_type:] == a_type:
        nametype = a_type
        sample_dir = prefix.str[:-len_type]
        break


if nametype is None:
    # pattern 搜索模式，可能会出 bug
    # find "_R2." or "_2." in prefix[0]
    pattern = re.compile(r'(_R2\.|_2\.)')
    m = pattern.search(prefix[0])
    if m:
        nametype = prefix[0][m.span()[0]:]
        len_type = len(nametype)
        sample_dir = prefix.str[:-len_type]

assert nametype is not None, 'No fastq detected or the file name is invaild!'

sample_name = sample_dir.apply(os.path.basename)

dict_yaml = {
    'suffix':suffix[0],
    'sample':dict(zip(sample_name,sample_dir)),
    'input_dir':args.folder,
    'output_dir':args.outdir,
    'thread':args.thread,
    'index':args.index,
    'fasta':args.ref,
    'binsize':args.binsize,
    'blacklist':blacklist,
    'nametype':nametype,
    'genomelen':dir_chrom_sizes,
    'script_folder':script_folder
    }

with open( os.path.join(args.outdir,'config.yaml'), 'w') as outfile:
    yaml.dump(dict_yaml, outfile, default_flow_style=False)

snakefile = os.path.join(script_dir, 'mapping/Snakefile_Trackseq')
shutil.copy(snakefile, os.path.join(args.outdir,'Snakefile'))


