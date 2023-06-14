#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys,re,time

if sys.version_info < (3,0):
    import platform
    raise Exception(f'python3 is needed, while running {platform.python_version()} now')

import offtracker
from offtracker.X_general import *
script_dir = os.path.abspath(os.path.dirname(offtracker.__file__))
script_folder= os.path.join(script_dir, 'mapping')

import argparse
import pandas as pd
import pybedtools
import multiprocessing as mp
from Bio.Blast.Applications import NcbiblastnCommandline 

def main():
    parser = argparse.ArgumentParser()
    parser.description='Generate candidate regions by sgRNA sequence'
    parser.add_argument('--sgrna' ,       type=str, required=True, help='sgRNA sequence without PAM' )
    parser.add_argument('--pam'   ,       type=str, required=True, help='The protospacer adjacent motif' )
    parser.add_argument('--name'  ,       type=str, required=True, help='custom name of the sgRNA' )
    parser.add_argument('-r','--ref'    , type=str, required=True, help='The fasta file of reference genome')
    parser.add_argument('-b','--blastdb', type=str, required=True, help='blast database')
    parser.add_argument('-o','--outdir' , type=str, required=True, help='The output folder')
    parser.add_argument('-g','--genome' , type=str, default='hg38', help='File of chromosome sizes, or "hg38", "mm10" ')
    parser.add_argument('-t','--thread' , type=int, default=4,     help='Number of threads to be used')
    parser.add_argument('--quick_mode'  , action='store_true',  help='BLAST faster but less candidates.')
    parser.add_argument('--regions'     , type=str, default='auto', nargs='+', help='Regions around candidate sites.' )
   
    args = parser.parse_args()
    
    
    if (args.genome == 'hg38') or (args.genome == 'mm10'):
        dir_chrom_sizes = os.path.join(script_folder, f'{args.genome}.chrom.sizes')
    else:
        dir_chrom_sizes = args.genome
    
    sgRNA_name = args.name
    sgRNA_seq  = args.sgrna
    PAM = args.pam
    n_threads  = args.thread
    dir_output = args.outdir
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    dir_ref_fa = args.ref
    blast_db   = args.blastdb
    quick_mode = args.quick_mode
    if args.regions == 'auto':
        regions = [500, 1000, 2000, 3000]
    else:
        regions = list(map(int, args.regions))
    common_chr = pd.Series(['chr']*23).str[:] + pd.Series(range(23)).astype(str).str[:]
    common_chr = pd.concat([common_chr, pd.Series(['chrX','chrY'])]).to_numpy()
    
    # parameters for alignment
    half_width = 100
    pct_params = 1.0
    frag_len= half_width*2
    location_len = regions[-1]
    dir_df_alignment = os.path.join(dir_output, f'df_alignment_{sgRNA_name}_{location_len}.csv')
    
    
    sgRNA_seq = sgRNA_seq.upper()
    PAM = PAM.upper()
    dir_sgRNA_fasta = os.path.join(dir_output, f'{sgRNA_name}_PAM.fasta')
    dir_sgRNA_blast = os.path.join(dir_output, f'{sgRNA_name}_PAM.blast')
    dir_sgRNA_bed   = os.path.join(dir_output, f'{sgRNA_name}_PAM.bed')
    
    
    possible_sgRNA_PAM = list(product([sgRNA_seq],possible_seq(PAM)))
    possible_sgRNA_PAM = [''.join(combination) for combination in possible_sgRNA_PAM]
    n_seq = len(possible_sgRNA_PAM)
    
    ID = pd.Series(['seq']*n_seq) + pd.Series(range(1,n_seq+1)).astype(str)
    df_sgRNA_PAM = pd.DataFrame({'ID':ID,'sequence':possible_sgRNA_PAM})
    write_fasta(df_sgRNA_PAM, dir_sgRNA_fasta)
    
    
    
    #########
    # BLAST #
    #########
    if os.path.isfile(dir_sgRNA_blast):
        print(f'{dir_sgRNA_blast} exists, skipped.')
    else:
        if quick_mode:
            print('Using quick mode for BLAST')
            blastx_cline = NcbiblastnCommandline(query=dir_sgRNA_fasta, task='blastn-short',out=dir_sgRNA_blast,
                                                db=blast_db, evalue=10000,outfmt=6, num_threads=n_threads,
                                                gapopen=4, gapextend=2, reward=2, word_size=5, dust='no', soft_masking=False)
        else:
            blastx_cline = NcbiblastnCommandline(query=dir_sgRNA_fasta, task='blastn-short',out=dir_sgRNA_blast,
                                                db=blast_db, evalue=100000,outfmt=6, num_threads=n_threads,
                                                gapopen=4, gapextend=2, reward=2, word_size=4, dust='no', soft_masking=False)   
        print(f'BLAST for candidate off-target sites of {sgRNA_name}.')
        blastx_cline()
        print(f'BLAST finished.')
    
    ##############
    # Output bed #
    ##############
    
    blast_regions = pd.read_csv(dir_sgRNA_blast, sep='\t',header=None)
    blast_regions.columns = ['query acc.','chr','% identity','alignment length','mismatches','gap opens','q. start','q. end','st','ed','evalue','bit score']
    blast_regions = blast_regions[blast_regions.evalue<10000]
    
    # reverse strand 
    blast_regions['reverse'] = (blast_regions['st']>blast_regions['ed']).astype(int)
    blast_regions_f = blast_regions[blast_regions.reverse==0].copy()
    blast_regions_r = blast_regions[blast_regions.reverse==1].copy()
    temp = blast_regions_r['st'].copy()
    blast_regions_r['st'] = blast_regions_r['ed']
    blast_regions_r['ed'] = temp
    blast_regions = pd.concat([blast_regions_f, blast_regions_r])
    # sort and add location
    blast_regions = blast_regions.sort_values('evalue').reset_index(drop=True)
    blast_regions['location']=blast_regions['chr'].str[:] + ':' + blast_regions['st'].astype(str).str[:] + '-' + blast_regions['ed'].astype(str).str[:]
    blast_regions = blast_regions.drop_duplicates(subset='location').copy()
    
    # alignment length 筛选
    len_sgRNA=len(sgRNA_seq)
    min_len = len_sgRNA-8
    blast_regions = blast_regions[blast_regions['alignment length']>=min_len].copy().reset_index(drop=True)
    blast_regions = blast_regions.reindex(columns = ['chr', 'st', 'ed' , 'query acc.', '% identity', 'alignment length', 'mismatches',
        'gap opens', 'q. start', 'q. end', 'evalue', 'bit score', 'reverse', 'location'] )
    
    # 输出 bed 用于后续 coverage 计算
    blast_regions_bed = blast_regions[['chr','st','ed']]
    writebed(blast_regions_bed, dir_sgRNA_bed)
    # 对 bed 进行排序但不合并
    a = pybedtools.BedTool(dir_sgRNA_bed)
    a.sort(g=dir_chrom_sizes).saveas( dir_sgRNA_bed )
    print(f'Output {sgRNA_name}_PAM.bed')
    
    
    ############################
    # Output candidate regions #
    ############################
    
    blast_regions_bed = X_readbed(dir_sgRNA_bed)
    blast_regions_bed = blast_regions_bed[blast_regions_bed['chr'].isin(common_chr)]
    blast_regions_bed['midpoint'] = ((blast_regions_bed['st'] + blast_regions_bed['ed'])/2).astype(int)
    blast_regions_bed = blast_regions_bed.drop_duplicates(subset=['chr','midpoint']).copy()
    for a_region in regions:
        candidate_region_left = blast_regions_bed.copy()
        candidate_region_left['ed'] = candidate_region_left['midpoint']
        candidate_region_left['st'] = candidate_region_left['midpoint']-a_region
        candidate_region_left.loc[candidate_region_left['st']<0,'st'] = 0
        # 储存并排序
        left_region =os.path.join(dir_output, f'{sgRNA_name}_candidate_left_{a_region}.bed')
        writebed(candidate_region_left.iloc[:,:3], left_region)
        a = pybedtools.BedTool(left_region)
        a.sort(g=dir_chrom_sizes).saveas( left_region )  
        
        candidate_region_right = blast_regions_bed.copy()
        candidate_region_right['st'] = candidate_region_right['midpoint']
        candidate_region_right['ed'] = candidate_region_right['midpoint']+a_region
        # 储存并排序
        right_region = os.path.join(dir_output, f'{sgRNA_name}_candidate_right_{a_region}.bed')
        writebed(candidate_region_right.iloc[:,:3], right_region)
        a = pybedtools.BedTool(right_region)
        a.sort(g=dir_chrom_sizes).saveas( right_region ) 
    
    # background noise
    for i in range(1,4):
        candidate_region_left = blast_regions_bed.copy()
        candidate_region_left['ed'] = candidate_region_left['midpoint']-5000*i
        candidate_region_left['st'] = candidate_region_left['midpoint']-5000*(i+1)
        candidate_region_left.loc[candidate_region_left['st']<0,'st'] = 0
        candidate_region_left.loc[candidate_region_left['ed']<5000,'ed'] = 5000
        # 储存并排序
        left_region =os.path.join(dir_output, f'{sgRNA_name}_candidate_left_bkg{i}.bed')
        writebed(candidate_region_left.iloc[:,:3], left_region)
        a = pybedtools.BedTool(left_region)
        a.sort(g=dir_chrom_sizes).saveas( left_region )  
        
        candidate_region_right = blast_regions_bed.copy()
        candidate_region_right['st'] = candidate_region_right['midpoint']+5000*i
        candidate_region_right['ed'] = candidate_region_right['midpoint']+5000*(i+1)
        # 储存并排序
        right_region = os.path.join(dir_output, f'{sgRNA_name}_candidate_right_bkg{i}.bed')
        writebed(candidate_region_right.iloc[:,:3], right_region)
        a = pybedtools.BedTool(right_region)
        a.sort(g=dir_chrom_sizes).saveas( right_region )     
    
    print(f'Output candidate regions of {sgRNA_name}.')     
    
    ###################
    # alignment score #
    ###################
    if os.path.isfile(dir_df_alignment):
        print(f'{dir_df_alignment} exists, skipped.')
    else:
        #########
        # 读取 blast bed
        #########
        bed_short = X_readbed(dir_sgRNA_bed)
        bed_short = bed_short[bed_short['chr'].isin(common_chr)].copy()
        bed_short['midpoint'] = ((bed_short['st'] + bed_short['ed'])/2).astype(int)
        bed_short['st'] = bed_short['midpoint'] - half_width 
        bed_short['ed'] = bed_short['midpoint'] + half_width
        bed_short.loc[bed_short['st']<0,'st']=0
        bed_short = bed_short.drop_duplicates()        
        
        #########
        # 根据 bed_f 位点 ed 前后 half_width 取基因组序列
        #########
        
        temp_bed = os.path.join(dir_output, 'temp.bed')
        writebed(bed_short.iloc[:,:3], temp_bed)
        a = pybedtools.BedTool(temp_bed)
        fasta = pybedtools.example_filename(dir_ref_fa)
        a = a.sequence(fi=fasta)
        with open(a.seqfn) as f:
            fasta = {} 
            for line in f:
                line = line.strip() # 去除末尾换行符
                if line[0] == '>':
                    header = line[1:]
                else:
                    sequence = line
                    fasta[header] = fasta.get(header,'') + sequence
        
        # pybedtools 得到位置 chrA:X-Y 时，X数字会往左多1bp
        
        #########
        # local alignment
        #########
        DNA_matrix = {('A','A'): 2, ('A','T'):0.01, ('A','C'):0.01, ('A','G'):0.01, ('A','N'):0.01,
                    ('T','T'): 2, ('T','A'):0.01, ('T','C'):0.01, ('T','G'):0.01, ('T','N'):0.01,
                    ('G','G'): 2, ('G','A'):0.01, ('G','C'):0.01, ('G','T'):0.01, ('G','N'):0.01,
                    ('C','C'): 2, ('C','A'):0.01, ('C','G'):0.01, ('C','T'):0.01, ('C','N'):0.01,
                    ('N','N'): 2, ('N','C'):2, ('N','A'): 2, ('N','G'): 2, ('N','T'): 2}
        mismatch_score = 0.01
        # 添加 PAM
        sgRNA_PAM_fw = sgRNA_seq + PAM
        sgRNA_PAM_rv = reverse_complement(sgRNA_PAM_fw)
        
        list_args_fw=[]
        list_args_rv=[]
        for a_key in fasta.keys():
            seq = re.sub('[^ATCG]','N',fasta[a_key])
            list_args_fw.append( [a_key, sgRNA_PAM_fw, seq, frag_len, DNA_matrix, mismatch_score] )
            list_args_rv.append( [a_key, sgRNA_PAM_rv, seq, frag_len, DNA_matrix, mismatch_score] )
        st = time.time()
        with mp.Pool(n_threads) as p:
            list_align_forward = p.starmap(sgRNA_alignment, list_args_fw)
        ed = time.time()
        print('align_forward:{:.2f}'.format(ed-st))
        st = time.time()
        with mp.Pool(n_threads) as p:
            list_align_reverse = p.starmap(sgRNA_alignment, list_args_rv)
        ed = time.time()
        print('align_reverse:{:.2f}'.format(ed-st))
        #
        df_align_forward = pd.DataFrame(list_align_forward, columns= ['fw_score','fw_pct','fw_target','fw_location','fw_deletion','fw_insertion','fw_mismatch'])
        df_align_reverse = pd.DataFrame(list_align_reverse, columns= ['rv_score','rv_pct','rv_target','rv_location','rv_deletion','rv_insertion','rv_mismatch'])
        df_align_reverse['rv_target'] = df_align_reverse['rv_target'].apply(reverse_complement)
        df_alignment = pd.concat([df_align_forward,df_align_reverse],axis=1)
        df_alignment['location'] = fasta.keys()
        df_alignment['alignment_score'] = df_alignment[['fw_score','rv_score']].max(axis=1)
        df_alignment['fw_score_2'] = df_alignment['fw_score']*(pct_params-df_alignment['fw_pct'].abs())
        df_alignment['rv_score_2'] = df_alignment['rv_score']*(pct_params-df_alignment['rv_pct'].abs())
        df_alignment['best_seq_score'] = df_alignment[['fw_score_2', 'rv_score_2']].max(axis=1)
        df_alignment['best_strand'] = df_alignment[['fw_score_2', 'rv_score_2']].idxmax(axis='columns').replace({'fw_score_2':'+', 'rv_score_2':'-'})
        df_alignment.loc[df_alignment['fw_score_2']==df_alignment['rv_score_2'],'best_strand']='equal_score'
        
        # GG check
        list_best_target = []
        list_best_location = []
        list_delete = []
        list_insert = []
        list_mismat = []
        list_GG = []
        for a_row in df_alignment.iterrows():
            if a_row[1]['best_strand']=='+':
                list_best_target.append(a_row[1]['fw_target'])
                list_best_location.append(a_row[1]['fw_location'])
                list_delete.append(a_row[1]['fw_deletion'])
                list_insert.append(a_row[1]['fw_insertion'])
                list_mismat.append(a_row[1]['fw_mismatch'])
                if a_row[1]['fw_target'][-2:]=='GG':
                    list_GG.append('OK')
                else:
                    list_GG.append('NO')                     
            elif a_row[1]['best_strand']=='-':
                list_best_target.append(a_row[1]['rv_target'])
                list_best_location.append(a_row[1]['rv_location'])
                list_delete.append(a_row[1]['rv_deletion'])
                list_insert.append(a_row[1]['rv_insertion'])
                list_mismat.append(a_row[1]['rv_mismatch'])
                if a_row[1]['rv_target'][-2:]=='GG':
                    list_GG.append('OK')
                else:
                    list_GG.append('NO')  
            else:
                if a_row[1]['fw_target'][-2:]=='GG':
                    list_best_target.append(a_row[1]['fw_target'])
                    list_best_location.append(a_row[1]['fw_location'])
                    list_delete.append(a_row[1]['fw_deletion'])
                    list_insert.append(a_row[1]['fw_insertion'])
                    list_mismat.append(a_row[1]['fw_mismatch'])
                    list_GG.append('OK_same_score')
                # 发现没有 GG 则看 RC
                elif a_row[1]['rv_target'][-2:]=='GG':
                    list_best_target.append(a_row[1]['rv_target'])
                    list_best_location.append(a_row[1]['rv_location'])
                    list_delete.append(a_row[1]['rv_deletion'])
                    list_insert.append(a_row[1]['rv_insertion'])
                    list_mismat.append(a_row[1]['rv_mismatch'])
                    list_GG.append('OK_same_score')
                else:
                    list_best_target.append(a_row[1]['fw_target'])
                    list_best_location.append(a_row[1]['fw_location'])
                    list_delete.append(a_row[1]['fw_deletion'])
                    list_insert.append(a_row[1]['fw_insertion'])
                    list_mismat.append(a_row[1]['fw_mismatch'])                    
                    list_GG.append('NO_same_score')
        # 记入 df_alignment
        df_alignment['deletion'] = list_delete
        df_alignment['insertion'] = list_insert
        df_alignment['mismatch'] = list_mismat
        df_alignment['GG'] = list_GG
        df_alignment['best_target'] = list_best_target
        df_alignment['target_location'] = list_best_location
        
        # 和 df_pivot 一致，左右各 location_len
        bed_short['st'] = bed_short['midpoint'] - location_len 
        bed_short['ed'] = bed_short['midpoint'] + location_len
        bed_short.loc[bed_short['st']<0,'st']=0
        df_alignment.index = igvfmt(bed_short)
        df_alignment.to_csv(dir_df_alignment)
        print(f'Output df_alignment_{sgRNA_name}_{location_len}.csv')
        os.remove(temp_bed)
    
    return 'Done!'
    

if __name__ == '__main__' :
    result = main()
    print(result)



