#!/usr/bin/env python
# -*- coding: utf-8 -*-

import polars as pl
import pandas as pd
import numpy as np
import offtracker
import argparse
import os, glob
import shlex, subprocess
from scipy.stats import norm

def main():
    parser = argparse.ArgumentParser()
    parser.description='New function in 2026. Check and correct potential incorrect target locations.'
    parser.add_argument('-f','--folder'  , type=str, required=True,    nargs='+', help='Directory of the data folder.' )
    parser.add_argument('--name'         , type=str, required=True,    help='custom name of the sgRNA' )
    parser.add_argument('--exp'          , type=str, default='all',    nargs='+', help='A substring mark in the name of experimental samples. The default is to use all samples other than control' )
    parser.add_argument('--control'      , type=str, default='none',   nargs='+', help='A substring mark in the name of control samples. The default is no control. "others" for all samples other than --exp.' )
    parser.add_argument('--fdr'          , type=float, default=0.05,   help='FDR threshold for the final result. Default is 0.05.')
    parser.add_argument('--score'        , type=float, default=1.9,    help='Track score threshold for the final result. Default is 1.9.')
    parser.add_argument('--smooth'       , type=int, default=1,        help='Smooth strength for the signal.')
    parser.add_argument('--window'       , type=int, default=3,        help='Window size for smoothing the signal.')
    parser.add_argument('--binsize'      , type=int, default=100,      help='Window size for smoothing the signal.')
    parser.add_argument('--flank_max'    , type=int, default=100000,   help='Maximun flanking distance from the candidate site.')
    parser.add_argument('--flank_regions', type=int, default=[1000,2000,3000,5000], nargs='+',help='flanking regions for calculating signal.')
    parser.add_argument('--SeqScorePower', type=float, default=4,      help='The seq score power' )
    parser.add_argument('--CtrClip'      , type=float, default=-0.5,   help='The lower clip for control samples' )
    parser.add_argument('-t','--thread'  , type=int, default=4,        help='Number of threads for parallel computing')
    parser.add_argument('-g','--genome'  , type=str, default='hg38',   help='File of chromosome sizes, or "hg38", "mm10" ')
    parser.add_argument('-o','--outdir'  , type=str, default='first',  help='The output folder. Default is the first folder of --folder' )
    parser.add_argument('--outname'      , type=str, default='same',   help='The suffix of output files. Default is the same --exp' )
    # new argument
    parser.add_argument('-r','--ref'     , type=str, required=True,    help='The fasta file of reference genome')
    parser.add_argument('--sgrna' ,        type=str, required=True,    help='One sgRNA sequence without PAM' )
    parser.add_argument('--pam'   ,        type=str, required=True,    help='The protospacer adjacent motif' )
    parser.add_argument('--pam_location',  type=str, default='downstream', help='Upstream or downstream, default is downstream (Cas9)' )
    # not used
    parser.add_argument('--seqfolder'    , type=str, default='none',    help='Actually not used in this script.Only in case you forget to remove this argument.')

    args = parser.parse_args()
    # 2025.08.08. 增加对阳性位点的 target_location 重比对功能，避免 blast 比对后的 realign 在更大范围内的存在不准确的情况
    # 实验性功能，如果 exp 有多个样本的话目前只取第一个 bdg 来分析

    ##########################
    ## parameter initiation ##
    ##########################

    folders = args.folder
    sgRNA_name = args.name + '_loc_correction'
    pattern_exp = args.exp
    pattern_ctr = args.control
    fdr_thresh = args.fdr
    score_thresh = args.score
    binsize = args.binsize
    flank_max = args.flank_max
    flank_regions = args.flank_regions # 如果 analysis 时修改了这个参数没有写 1000 的话会出bug，暂时懒得改了
    smooth_times = args.smooth
    window_size = args.window
    seq_score_power = args.SeqScorePower
    ctr_clip = args.CtrClip


    if args.outname == 'same':
        if isinstance(pattern_exp, list):
            outname = '_'.join(pattern_exp)
        else:
            outname = pattern_exp
    else:
        outname = args.outname

    outdir = args.outdir
    if outdir == 'first':
        outdir = folders[0]
    os.chdir(outdir)
    # out temp folder
    if not os.path.exists( os.path.join(outdir,'temp') ):
        os.makedirs(os.path.join(outdir,'temp'))
    # data temp folder
    for a_folder in folders:
        temp_dir = os.path.join(a_folder, 'temp')
        if not os.path.exists( temp_dir ):
            os.makedirs(temp_dir)

    ##################
    ## glob samples ##
    ##################
    all_sample_names = []
    all_sample_files = []
    for a_folder in folders:    
        bdg_files = pd.Series(glob.glob(os.path.join( a_folder, '*.add.bdg' ))).sort_values().reset_index(drop=True)
        sample_names = bdg_files.apply(os.path.basename).str.extract(r'(.*)\.\d+\.add\.bdg',expand=False)
        all_sample_names.extend( sample_names )
        all_sample_files.extend( bdg_files )
    all_sample_files = pd.Series(all_sample_files)
    all_sample_names = pd.Series(all_sample_names)
    print('all sample names in the folders:')
    print(all_sample_names)
    print('your string pattern for experimental groups: ', pattern_exp)
    ctr_samples = []
    if pattern_ctr == 'none':
        if pattern_exp == 'all':
            exp_samples = list( all_sample_names )
        else:
            exp_samples = []
            for a_mark in pattern_exp:
                exp_samples.extend( list( all_sample_names[all_sample_names.str.contains(a_mark)] ) )
    elif pattern_ctr == 'others':
        if pattern_exp == 'all':
            exp_samples = list( all_sample_names )
        else:
            exp_samples = []
            for a_mark in pattern_exp:
                exp_samples.extend( list( all_sample_names[all_sample_names.str.contains(a_mark)] ) )
            ctr_samples = list( all_sample_names[~all_sample_names.isin(exp_samples)] )
    else:
        for a_mark in pattern_ctr:
            ctr_samples.extend( list( all_sample_names[all_sample_names.str.contains(a_mark)] ) )
        if pattern_exp == 'all':
            exp_samples = list( all_sample_names[~all_sample_names.isin(ctr_samples)] )
        else:
            exp_samples = []
            for a_mark in pattern_exp:
                exp_samples.extend( list( all_sample_names[all_sample_names.str.contains(a_mark)] ) )
    n_exp = len(exp_samples)
    n_ctr = len(ctr_samples)
    print(f'Experimental group has {n_exp} samples:\n{exp_samples}')
    print(f'Control group has {n_ctr} samples:\n{ctr_samples}')

    # mark 错误时
    assert n_exp > 0, 'No experimental sample is found. Please check the name pattern.'
    if (n_ctr==0)&(pattern_ctr != 'none'):
        print('Name pattern for control sample(s) was given, but no file meet the pattern.')
        return 'Program terminated'

    # summarize samples
    bool_exp = all_sample_names.isin(exp_samples)
    bool_ctr = all_sample_names.isin(ctr_samples)
    exp_sample_files = all_sample_files[bool_exp]
    ctr_sample_files = all_sample_files[bool_ctr]
    exp_sample_names = all_sample_names[bool_exp]
    ctr_sample_names = all_sample_names[bool_ctr]
    # selected_sample_files = pd.concat([exp_sample_files,ctr_sample_files])
    # selected_sample_names = pd.concat([exp_sample_names,ctr_sample_names]) # no use



    ####################
    ## run correction ##
    ####################

    # new parameters
    ref_fasta = args.ref
    sgRNA_seq = args.sgrna
    PAM = args.pam
    PAM_loc = args.pam_location
    # read result
    dp_result = pl.read_csv(f'./temp/df_result_{outname}.csv')
    # negative for next section
    bool_fdr_bkg = dp_result['fdr']>fdr_thresh
    bool_score_bkg = dp_result['track_score']<score_thresh
    dp_result_bkg = dp_result.filter(bool_fdr_bkg & bool_score_bkg)
    # positive
    bool_fdr = pl.col('fdr')<=fdr_thresh
    bool_score = pl.col('track_score')>=score_thresh
    dp_result = dp_result.filter(bool_fdr & bool_score)
    # bdg
    dp_bdg = pl.read_csv(exp_sample_files.iloc[0], separator='\t', has_header=False,
                             schema_overrides={'chr':pl.String,'start':pl.Int32,'end':pl.Int32,'residual':pl.Float32})
    # check and realign
    bool_left_neg=(dp_result['exp_L_neg_1000']<-5)&(dp_result['exp_R_neg_1000']==0)
    bool_right_neg=(dp_result['exp_R_neg_1000']<-5)&(dp_result['exp_L_neg_1000']==0)
    list_good_result = []
    list_bad_left = []
    list_bad_right = []
    n_left_for_correct = 0
    n_right_for_correct = 0
    for a_left_bool, a_right_bool, a_row in zip(bool_left_neg, bool_right_neg, dp_result.iter_rows(named=True)):
        if a_left_bool & a_right_bool:
            raise ValueError('abnormal on both left and right')
        if a_left_bool:
            n_left_for_correct += 1
            loc_shift_left = a_row['chr'] + ':' + str(a_row['st']-1000) + '-' + str(a_row['ed']-20)
            region_index = a_row['region_index']
            dp_bdg_chr = dp_bdg.filter(pl.col('chr') == a_row['chr'])
            sr_candidate = offtracker.left_realign(dp_bdg_chr, loc_shift_left, ref_fasta, sgRNA_seq, PAM, PAM_loc, n_iter=0)
            sr_candidate.loc['region_index'] = region_index
            list_bad_left.append(sr_candidate)
        elif a_right_bool:
            n_right_for_correct += 1
            loc_shift_right = a_row['chr'] + ':' + str(a_row['st']+20) + '-' + str(a_row['ed']+1000)
            region_index = a_row['region_index']
            dp_bdg_chr = dp_bdg.filter(pl.col('chr') == a_row['chr'])
            sr_candidate = offtracker.right_realign(dp_bdg_chr, loc_shift_right, ref_fasta, sgRNA_seq, PAM, PAM_loc, n_iter=0)
            sr_candidate.loc['region_index'] = region_index
            list_bad_right.append(sr_candidate)
        else:
            list_good_result.append(a_row)
    dp_result_good = pl.DataFrame(list_good_result)
    df_cand_left = pd.DataFrame(list_bad_left)
    df_cand_right = pd.DataFrame(list_bad_right)
    df_cand_realign = pd.concat([df_cand_left, df_cand_right])
    if len(df_cand_realign) == 0:
        print('No candidate is found for realignment.')
        return 'finished'

    # 情况判断
    n_success_realign = sum(df_cand_realign['realign']=='success')
    n_fail_realign = sum(df_cand_realign['realign']!='success')
    if (n_success_realign == 0) and (n_fail_realign > 0):
        print(f'{n_fail_realign} candidates are found for realignment, but all failed.')
        return 'finished'
    elif (n_success_realign > 0) and (n_fail_realign > 0):
        print(f'{n_success_realign} candidates succeeded, and {n_fail_realign} candidates failed.')
    else:
        print(f'{n_success_realign} candidates succeeded.')

    df_cand_realign = df_cand_realign[df_cand_realign['realign']=='success']
    seqfile = rf'correction_df_candidate_{outname}_realign.csv'
    df_cand_realign.to_csv(seqfile)
    
    # run offtracker_analysis with check_loc mode
    running_log = rf'correction_analysis_{outname}.log'
    # list 转空格分割参数
    if isinstance(pattern_exp, list):
        param_pattern_exp = ' '.join(pattern_exp)
    else:
        param_pattern_exp = pattern_exp
    if isinstance(pattern_ctr, list):
        param_pattern_ctr = ' '.join(pattern_ctr)
    else:
        param_pattern_ctr = pattern_ctr
    if isinstance(flank_regions, list):
        param_flank_regions = ' '.join([str(x) for x in flank_regions])
    else:
        param_flank_regions = flank_regions
    if isinstance(folders, list):
        param_folders = ' '.join([str(x) for x in folders])
    else:
        param_folders = folders

    with open(running_log, "w+") as running_log:
        command   = f'offtracker_analysis.py -t {args.thread} -g {args.genome} --seqfile {seqfile} --name {sgRNA_name} \
        --exp {param_pattern_exp} --control {param_pattern_ctr} --outname {outname}_loc_correction -f {param_folders} -o {outdir} \
        --fdr {fdr_thresh} --window {window_size} --smooth {smooth_times} --SeqScorePower {seq_score_power} \
        --score {score_thresh} --binsize {binsize} --flank_max {flank_max} --flank_regions {param_flank_regions} --CtrClip {ctr_clip} \
        --check_loc'
        command2  = shlex.split('bash -c "{}"'.format(command))
        process_1 = subprocess.Popen(command2, stdout=running_log, stderr=subprocess.STDOUT )
        process_1.wait(timeout=100000)
        retc  = process_1.returncode
        if retc==0:
            print((f'correction_analysis {outname} is done!'))
        else:
            print((f'correction_analysis {outname} is failed!'))


    #######################
    ## recalculate score ##
    #######################

    dp_result_realign = pl.read_csv(f'./temp/df_result_{outname}_loc_correction.csv')

    # 兼容旧版输出列名
    list_col = dp_result_realign.columns[:-5]
    dp_result_new = pl.concat([dp_result_realign[list_col], dp_result_good[list_col], dp_result_bkg[list_col]])

    # 标准化分布, polars 版   
    target_std=0.15
    n_outliers = int(np.ceil(len(dp_result_new)*0.01))
    score_bkg = dp_result_new['raw_score'][n_outliers:-n_outliers]
    mean_score_bkg = score_bkg.mean()
    std_score_bkg = score_bkg.std()
    dp_result_new = dp_result_new.with_columns(
        (pl.col('raw_score').sub(mean_score_bkg)/std_score_bkg).alias('track_score')
    )
    dp_result_new = dp_result_new.with_columns(
        pl.col('track_score').mul(target_std).add(1).alias('track_score')
    )
    dp_result_new = dp_result_new.with_columns(
        pl.col('track_score').clip(lower_bound=0.5).log(base=2).alias('log2_track_score')
    )
    dp_result_new = dp_result_new.sort('track_score', descending=True)

    # pv and fdr
    score_for_fitting = dp_result_new['log2_track_score'][n_outliers:-n_outliers]
    mu, std = norm.fit(score_for_fitting) 
    print('mean_score:{:.3f};std:{:.3f}'.format(mu,std))
    dp_result_new = dp_result_new.with_columns(
        pl.col('log2_track_score').map_elements( lambda x: norm.sf(x,loc=mu,scale=std), return_dtype=pl.Float64 ).clip(lower_bound=1e-320).alias('pv')
    )
    dp_result_new = dp_result_new.with_columns(
        fdr=offtracker.fdr(dp_result_new['pv']).alias('fdr'),
        rank=pl.Series(range(1,len(dp_result_new)+1))
    ) #.with_row_index(name='rank',offset=1)
    dp_result_new.write_csv(f'./temp/df_result_{outname}.csv') # 覆盖原结果

    # ouput Offtracker result
    bool_fdr = pl.col('fdr')<=fdr_thresh
    bool_score = pl.col('track_score')>=score_thresh
    dp_output = dp_result_new.filter(bool_fdr|bool_score)
    if pattern_ctr != 'none':
        dp_output = dp_output[['target_location', 'best_strand','best_target','deletion','insertion','mismatch',
                            'exp_L_length', 'exp_R_length','ctr_L_length','ctr_R_length','L_length','R_length','signal_length',
                            'norm_best_seq_score','track_score', 'log2_track_score','fdr','rank']]
        dp_output.columns = ['target_location', 'strand', 'target', 'deletion', 'insertion', 'mismatch',
                            'exp_L_length', 'exp_R_length','ctr_L_length','ctr_R_length','L_length','R_length','signal_length',
                            'seq_score', 'track_score', 'log2_track_score','FDR', 'rank']
    else:
        dp_output = dp_output[['target_location', 'best_strand','best_target','deletion','insertion','mismatch',
                            'L_length', 'R_length','signal_length',
                            'norm_best_seq_score','track_score', 'log2_track_score','fdr','rank']]
        dp_output.columns = ['target_location', 'strand', 'target', 'deletion', 'insertion', 'mismatch',
                            'L_length', 'R_length','signal_length',
                            'seq_score', 'track_score', 'log2_track_score','FDR', 'rank']
    dp_output.write_csv(f'Offtracker_result_{outname}.csv')

    return 'correction finished'

if __name__ == '__main__' :
    result = main()
    print(result)


