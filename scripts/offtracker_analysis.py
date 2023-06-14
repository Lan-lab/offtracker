#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,glob,sys,re,time,shutil

if sys.version_info < (3,0):
    import platform
    raise Exception(f'python3 is needed, while running {platform.python_version()} now')

import offtracker
from offtracker.X_analysis import *
script_dir = os.path.abspath(os.path.dirname(offtracker.__file__))
script_folder= os.path.join(script_dir, 'mapping')

import argparse
import pandas as pd
import numpy as np
import pybedtools
from scipy.stats import norm

def main():
    parser = argparse.ArgumentParser()
    parser.description='Analyze the ssChIP-seq data.'
    parser.add_argument('-f','--folder' , type=str, required=True,  nargs='+', help='Directory of the data folder.' )
    parser.add_argument('--seqfolder'   , type=str, required=True,  help='Directory of sgRNA information created by seq_cadidates.')
    parser.add_argument('--name'        , type=str, required=True,  help='custom name of the sgRNA' )
    parser.add_argument('--exp'         , type=str, default='all',  nargs='+', help='A substring mark in the name of experimental samples. The default is to use all samples other than control' )
    parser.add_argument('--control'     , type=str, default='none', nargs='+', help='A substring mark in the name of control samples. The default is no control. "others" for all samples other than --exp.' )
    parser.add_argument('--regions'     , type=str, default='auto', nargs='+', help='Regions around candidate sites.' )
    parser.add_argument('--rgO'         , type=str, default='mean', help='Regoins operation. mean/min/max for regions' )
    parser.add_argument('--repO'        , type=str, default='mean', help='Replicate operation for multiple experimental samples (mean/min/max)' )
    parser.add_argument('--repOC'       , type=str, default='mean', help='Replicate operation for multiple control samples (mean/max)' )
    parser.add_argument('-g','--genome' , type=str, default='hg38', help='File of chromosome sizes, or "hg38", "mm10" ')
    parser.add_argument('-o','--outdir' , type=str, default='first',help='The output folder. Default is the first folder of --folder' )
    parser.add_argument('--outname'     , type=str, default='same', help='The suffix of output files. Default is the same --exp' )
    parser.add_argument('--overwrite'   , action='store_true', help='Whether to overwrite existed dataframes.' )
    parser.add_argument('--NoShapeBonus', action='store_true', help='Disable shape bonus for each site.' )
    parser.add_argument('--SignalFormula',type=int, default=4, help='The signal formula' )
    parser.add_argument('--ScoreFormula', type=int, default=3, help='The score formula' )
    parser.add_argument('--ShapeFormula', type=int, default=2, help='The shape formula' )
    parser.add_argument('--ShapeThresh',  type=str, default='auto', help='The threshold for shape formula' )
    parser.add_argument('--ShapeWeight',  type=str, default='auto', help='The weight for shape bonus' )
    parser.add_argument('--SeqScorePower',type=float, default=2, help='The seq score power' )
    parser.add_argument('--ScorePseudocount',type=float, default=0.9, help='The score pseudocount' )
    parser.add_argument('--RescueFactor',type=float, default=1.2, help='The rescue factor' )
    parser.add_argument('--RescueThresh1',type=float, default=1.5, help='The rescue threshold for min raw' )
    parser.add_argument('--RescueThresh2',type=float, default=1.8, help='The rescue threshold for raw score' )
    parser.add_argument('--clean', action='store_true', help='Whether to remove temp files')
    
    print(f'Runing offtracker verision: {offtracker.__version__}')
    # main parameters
    args = parser.parse_args()
    folders = args.folder
    seq_folder = args.seqfolder
    sgRNA_name = args.name
    control_mark = args.control
    exp_mark = args.exp
    ctr_samples = []
    outdir = args.outdir
    if outdir == 'first':
        outdir = folders[0]
    os.chdir(outdir)
    for a_folder in folders:
        temp_dir = os.path.join(a_folder, 'temp')
        if not os.path.exists( temp_dir ):
            os.makedirs(temp_dir)
    
    # parameters for algorithm
    operator = 'p'
    
    dict_signal_formula = {1:signal_formula1, 2:signal_formula2, 3:signal_formula3, 4:signal_formula4}
    dict_score_formula = {1:score_formula1, 2:score_formula2, 3:score_formula3, 4:score_formula4, 5:score_formula5}
    dict_shape_formula = {1:shape_formula1, 2:shape_formula2}

    the_signal_formula = dict_signal_formula[args.SignalFormula]
    if control_mark == 'none':
        print('No control samples, all samples in the folder are regarded as experimental samples.')
        the_score_formula = score_formula5
    else:
        the_score_formula = dict_score_formula[args.ScoreFormula]
    if args.NoShapeBonus:
        the_shape_formula = None
        print('the_shape_formula is None')
    else:
        the_shape_formula = dict_shape_formula[args.ShapeFormula]
    
    if args.ShapeFormula == 1:
        ratio_thresh = 1.5
        exp_weight = 0.3
    elif args.ShapeFormula == 2:
        ratio_thresh = 1
        exp_weight = 0.5
    
    if args.ShapeThresh != 'auto':
        ratio_thresh = float(args.ShapeThresh)
    if args.ShapeWeight != 'auto':
        exp_weight = float(args.ShapeWeight)
    
    seq_score_power = args.SeqScorePower
    if args.regions == 'auto':
        regions = [500, 1000, 2000, 3000]
    else:
        regions = list(map(int, args.regions))
    print('regions:',regions)
    bkgs = ['bkg1','bkg2','bkg3']
    noise_length = 5000
    max_region = max(regions)
    region_op = args.rgO
    repO = args.repO
    repOC = args.repOC
    score_pseudo = args.ScorePseudocount
    rescue_factor = args.RescueFactor
    assert repO in ['mean','min','max'], "--repO only accepts \"mean\",\"min\",or \"max\" "
    assert repOC in ['mean','max'], "--repO only accepts \"mean\" or \"max\" "
    print('repO',repO)
    print('repOC',repOC)
    print('ShapeThresh:',ratio_thresh)
    print('ShapeWeight:',exp_weight)
    print('SignalFormula:',args.SignalFormula)
    print('ScoreFormula:',args.ScoreFormula)
    print('ShapeFormula:',args.ShapeFormula)
    print('SeqScorePower:',args.SeqScorePower)
    print('ScorePseudocount:',args.ScorePseudocount)
    print('RescueFactor:',args.RescueFactor)
    # glob samples, check paired
    all_samples_name = []
    all_samples_pref = []
    for a_folder in folders:    
        forward_bed_files = glob.glob( os.path.join( a_folder, '*.fw.bed' ) )
        forward_bed_files.sort()
        reverse_bed_files = glob.glob( os.path.join( a_folder, '*.rv.bed' ) )
        reverse_bed_files.sort()
        forward_samples = pd.Series(forward_bed_files).apply(os.path.basename).str[:-7]
        reverse_samples = pd.Series(reverse_bed_files).apply(os.path.basename).str[:-7]
        assert (forward_samples == reverse_samples).all()
        all_samples_name.extend( forward_samples )
        all_samples_pref.extend( pd.Series(forward_bed_files).str[:-7] )
    all_samples_pref = pd.Series(all_samples_pref)
    all_samples_name = pd.Series(all_samples_name)

    if control_mark == 'none':
        if exp_mark == 'all':
            exp_samples = list( all_samples_name )
        else:
            exp_samples = []
            for a_mark in exp_mark:
                exp_samples.extend( list( all_samples_name[all_samples_name.str.contains(a_mark)] ) )
    elif control_mark == 'others':
        if exp_mark == 'all':
            exp_samples = list( all_samples_name )
        else:
            exp_samples = []
            for a_mark in exp_mark:
                exp_samples.extend( list( all_samples_name[all_samples_name.str.contains(a_mark)] ) )
            ctr_samples = list( all_samples_name[~all_samples_name.isin(exp_samples)] )
    else:
        for a_mark in control_mark:
            ctr_samples.extend( list( all_samples_name[all_samples_name.str.contains(a_mark)] ) )
        if exp_mark == 'all':
            exp_samples = list( all_samples_name[~all_samples_name.isin(ctr_samples)] )
        else:
            exp_samples = []
            for a_mark in exp_mark:
                exp_samples.extend( list( all_samples_name[all_samples_name.str.contains(a_mark)] ) )
    n_exp = len(exp_samples)
    n_ctr = len(ctr_samples)
    print(f'Experimental group has {n_exp} samples:\n{exp_samples}')
    print(f'Control group has {n_ctr} samples:\n{ctr_samples}')

    # mark 错误时
    assert n_exp > 0, 'No experimental sample is found. Please check the name pattern.'
    if (n_ctr==0)&(control_mark != 'none'):
        print('Name pattern for control sample(s) was given, but no file meet the pattern.')
        return 'Program terminated'
    
    # sequence score
    try:
        df_alignment = pd.read_csv(os.path.join(seq_folder, f'df_alignment_{sgRNA_name}_{max_region}.csv'), index_col=0)
    except FileNotFoundError:
        return 'Please run offtracker_candidates.py first and provide the correct directory with --seqfolder'
    
    # chromosome sizes
    if (args.genome == 'hg38') or (args.genome == 'mm10'):
        dir_chrom_sizes = os.path.join(script_folder,f'{args.genome}.chrom.sizes')
    else:
        dir_chrom_sizes = args.genome
    
    # sgRNA length
    sgRNA_fa = os.path.join(seq_folder, sgRNA_name + '_PAM.fasta')
    if not os.path.isfile(sgRNA_fa):
        raise Exception(f'{sgRNA_name}_PAM.fasta is not in the seqfolder')
    with open(sgRNA_fa,'r') as f:
        temp = f.readlines()
        len_sgRNA_PAM = len(temp[1].strip())
    #
    selected_samples_pref = all_samples_pref[all_samples_name.isin(exp_samples+ctr_samples)]
    selected_samples_name = all_samples_name[all_samples_name.isin(exp_samples+ctr_samples)]
    # read counts on candidate regions
    for a_pref in selected_samples_pref:
        cand_count(a_pref, sgRNA_name, regions, seq_folder, dir_chrom_sizes, overwrite=False)
        bkg_count(a_pref, sgRNA_name, bkgs, seq_folder, dir_chrom_sizes, overwrite=False)
    

    ############
    ## signal ##
    ############
    if args.outname == 'same':
        outname = exp_mark
    else:
        outname = args.outname
    output = f'./temp/df_pivot_{outname}.csv'
    if (os.path.isfile(output))&(not args.overwrite):
        print(f'skip {output}')
    else:
        # 计算每个样本的 df_counts
        dict_df_counts = {}
        for a_pref, a_sample in zip(selected_samples_pref, selected_samples_name):
            dirname = os.path.dirname(a_pref)
            basename = os.path.basename(a_pref)
            temp_dir = os.path.join(dirname, 'temp')
            a_pref = os.path.join(temp_dir, basename)
            dir_df_counts = f'{a_pref}_{sgRNA_name}_df_counts_{the_signal_formula.__name__}_{operator}.csv'
            if (os.path.isfile(dir_df_counts))&(not args.overwrite):
                print(f'df_counts of {a_sample}_{sgRNA_name}_{the_signal_formula.__name__}_{operator} exists, loading.\n')
                df_counts_temp = pd.read_csv(dir_df_counts, index_col=0)
            else:
                print(f'For {a_sample}')
                df_counts_temp, df_noise_temp = load_count(a_sample, regions, bkgs, sgRNA_name = sgRNA_name, signal_formula = the_signal_formula, noise_length=noise_length,
                                                    ratio_thresh=ratio_thresh, exp_weight=exp_weight, shape_formula=the_shape_formula, operator = operator,
                                                    region_op = region_op, dirname = temp_dir, pseudo_count=1)
                print('\n')
                df_counts_temp['noise_1kb'] = df_noise_temp['noise_bp']*1000
                df_counts_temp.to_csv(dir_df_counts)
            dict_df_counts[a_sample] = df_counts_temp
        # 简化 df_counts 再合并
        dict_df_counts_sub = {}
        for a_sample in selected_samples_name:
            df_counts_sub = dict_df_counts[a_sample].copy()
            df_counts_sub = df_counts_sub[['chr','st','ed','midpoint','location','ID_1','ID_2','left_signal','right_signal','signal_FC', 'signal_min','shape_bonus','score']].copy()
            df_counts_sub['sample'] = a_sample
            dict_df_counts_sub[a_sample] = df_counts_sub

        # 合成 pivot 矩阵
        df_pivot = pd.concat(dict_df_counts_sub).pivot(index='location', columns='sample', values='score').fillna(0)
        df_shape = pd.concat(dict_df_counts_sub).pivot(index='location', columns='sample', values='shape_bonus').fillna(0)
        df_signal_FC = pd.concat(dict_df_counts_sub).pivot(index='location', columns='sample', values='signal_FC').fillna(0)
        df_signal_min = pd.concat(dict_df_counts_sub).pivot(index='location', columns='sample', values='signal_min').fillna(0)

        df_pivot['std_all'] = df_pivot.std(axis=1)
        if control_mark != 'none':
            df_pivot['control_mean'] = df_pivot[ctr_samples].mean(axis=1)
            df_pivot['control_max']  = df_pivot[ctr_samples].max(axis=1)
            df_pivot['control_std']  = df_pivot[ctr_samples].std(axis=1)
            df_pivot['control_signal_FC'] = df_signal_FC[ctr_samples].mean(axis=1)
            df_pivot['control_signal_min'] = df_signal_min[ctr_samples].mean(axis=1)
            if repOC == 'mean':
                df_pivot['control_raw'] = df_pivot['control_mean']
                df_pivot['control_shape'] = df_shape[ctr_samples].mean(axis=1)
            elif repOC == 'max':
                df_pivot['control_raw'] = df_pivot['control_max']
                df_pivot['control_shape'] = df_shape[ctr_samples].max(axis=1)
            df_pivot['control_signal'] = df_pivot['control_raw']*df_pivot['control_shape']
        else:
            df_pivot['control_signal'] = 'no_control'

        # group mean/min/max
        if repO == 'mean':
            df_pivot['exp_raw'] = df_pivot[exp_samples].mean(axis=1)
            df_pivot['exp_shape'] = df_shape[exp_samples].mean(axis=1)
        elif repO == 'min':
            df_pivot['exp_raw'] = df_pivot[exp_samples].min(axis=1)
            df_pivot['exp_shape'] = df_shape[exp_samples].min(axis=1)
        elif repO == 'max':
            df_pivot['exp_raw'] = df_pivot[exp_samples].max(axis=1)
            df_pivot['exp_shape'] = df_shape[exp_samples].max(axis=1)
        df_pivot['exp_signal_FC'] = df_signal_FC[exp_samples].mean(axis=1)
        df_pivot['exp_signal_min'] = df_signal_min[exp_samples].mean(axis=1)
        df_pivot['exp_signal'] = df_pivot['exp_raw']*df_pivot['exp_shape']

        # 添加坐标 和 ID
        bed_pivot = bedfmt(df_pivot.index)
        bed_pivot.index = df_pivot.index
        df_pivot = pd.concat([bed_pivot,df_pivot], axis=1)
        df_pivot['midpoint'] = df_pivot['st']+max_region

        point_head = (df_pivot['midpoint']/1000).astype(int)
        df_pivot['ID_1'] = df_pivot['chr'] + ':' + point_head.astype(str)
        point_tail = df_pivot['midpoint'] % 1000
        df_pivot.loc[point_tail<500,'ID_2'] = df_pivot['chr'] + ':' + (point_head-1).astype(str)
        df_pivot.loc[point_tail>=500,'ID_2'] = df_pivot['chr'] + ':' + (point_head+1).astype(str)

        df_alignment_sub = df_alignment.loc[df_pivot.index, ['best_strand','best_target','target_location','deletion','insertion','mismatch','GG','alignment_score','best_seq_score']]
        df_pivot = pd.concat([df_pivot,df_alignment_sub],axis=1)
        df_pivot.to_csv(output)
    ###########
    ## score ##
    ###########
    
    output = f'./temp/df_pivot_score_{outname}.csv'
    if (os.path.isfile(output))&(not args.overwrite):
        print(f'skip {output}')
        df_pivot =  pd.read_csv(output, index_col=0)
    else:
        df_pivot = pd.read_csv(f'./temp/df_pivot_{outname}.csv', index_col=0)
        df_pivot[f'raw_score'] = df_pivot[['exp_signal','control_signal']].apply(lambda x: the_score_formula(x['exp_signal'], x['control_signal'], score_pseudo), axis=1 )
        mean_seq_score = round(df_pivot['best_seq_score'].mean(),3)
        df_pivot['norm_best_seq_score'] = np.power(df_pivot['best_seq_score']/mean_seq_score, seq_score_power)
        df_pivot[f'final_score'] = df_pivot[f'raw_score']*df_pivot['norm_best_seq_score']
        # record the final score before rescue
        df_pivot['final_score_before_rescue'] = df_pivot['final_score'].copy()
        
        # raw rescue
        df_raw_rescue =  df_pivot[df_pivot['raw_score']>1.5].sort_values(by='raw_score', ascending=False)
        df_raw_rescue['min_raw'] = df_raw_rescue[exp_samples].min(axis=1)/df_raw_rescue['control_max'].clip(lower=0.5)
        # 对 min_raw rescue
        index_raw_rescue = df_raw_rescue[df_raw_rescue['min_raw']>args.RescueThresh1].index
        # 直接对 raw 进行 rescue
        index_raw_rescue_2 = df_raw_rescue[df_raw_rescue['raw_score']>args.RescueThresh2].index
        index_raw_rescue = np.union1d(index_raw_rescue,index_raw_rescue_2)
        print('rescue',len(index_raw_rescue),'sites')
        df_pivot['rescue_factor'] = 1
        df_pivot.loc[index_raw_rescue,'rescue_factor'] = rescue_factor
        df_pivot.loc[index_raw_rescue,f'final_score'] = df_pivot.loc[index_raw_rescue,f'final_score']*rescue_factor
        
        # dedup
        df_pivot = df_pivot.sort_values(by=f'final_score', ascending=False)
        list_nondup = dedup_two(df_pivot,'ID_1','ID_2')
        df_pivot = df_pivot[list_nondup]
        
        target_std=0.15
        df_pivot = df_pivot[df_pivot[f'final_score']>0].copy()
        n_bkg_sites = int(len(df_pivot)*0.99)
        score_bkg = df_pivot['final_score'][-n_bkg_sites:]
        mean_score_bkg = score_bkg.mean()
        std_score_bkg = score_bkg.std()
        df_pivot['norm_final_score'] = (df_pivot[f'final_score'] - mean_score_bkg) / std_score_bkg
        df_pivot['norm_final_score'] = df_pivot[f'norm_final_score']*target_std + 1
        df_pivot['norm_final_score'] = df_pivot['norm_final_score'].clip(lower=0)
        df_pivot['log2_norm_final_score'] = np.log2(df_pivot[f'norm_final_score']+1)
        
        # fitting normal distribution
        score_for_fitting = df_pivot['log2_norm_final_score'][-n_bkg_sites:]
        mu, std = norm.fit(score_for_fitting) 
        print('mean_score:{:.3f};std:{:.3f}'.format(mu,std))
        # pv and fdr
        df_pivot['pv'] = df_pivot[f'log2_norm_final_score'].apply( lambda x: norm.sf(x,loc=mu,scale=std) )
        df_pivot['pv'].clip(lower=1e-320,inplace=True)
        df_pivot.to_csv(output)
    
    ############
    ## filter ##
    ############
    
    signal_min_thresh = -0.1
    signal_min_thresh_2 = -2
    search_distance = 40000
    seq_score_thresh = len_sgRNA_PAM*2 - 5*2
    
    df_result = df_pivot.copy()
    ## 
    candidate_dup = df_result[ df_result[f'exp_signal_min']<=signal_min_thresh ].index[:500]
    list_dedup = []
    list_seq_score = []
    for a_loc in candidate_dup:
        temp_chr = df_result.loc[a_loc,'chr']
        temp_seq_score = df_result.loc[a_loc,'best_seq_score']
        # exp_signal_min 太低直接跳过得分判断
        if df_result.loc[a_loc,'exp_signal_min'] <= signal_min_thresh_2 :
            pass
        else:
            # 高于 X 分则跳过
            if temp_seq_score > seq_score_thresh:
                continue
        # 取排其前且同chr者
        temp_df_result = df_result.loc[:a_loc]
        temp_df_result = temp_df_result[temp_df_result['chr'] == temp_chr].iloc[:-1]
        if len(temp_df_result)==0:
            continue
        else:
            # 距离小于 search_distance 者记为信号溢出假阳性
            if (temp_df_result['midpoint']-df_result.loc[a_loc,'midpoint']).abs().min()<search_distance :
                list_dedup.append(a_loc)
                list_seq_score.append(temp_seq_score)
    
    df_result = df_result[~df_result.index.isin(list_dedup)].copy()
    print(f'filter {len(list_dedup)} sites')
    
    df_result['rank'] = range(1,len(df_result)+1)
    df_result['fdr'] = fdr(df_result['pv'])
    df_result.to_csv(f'./temp/df_result_{outname}.csv')
    
    df_output = df_result[df_result['log2_norm_final_score']>=1.3].copy()
    df_output = df_output[['target_location', 'best_strand','best_target','deletion','insertion','mismatch','norm_best_seq_score','norm_final_score', 'log2_norm_final_score','fdr','rank']]
    df_output.columns = ['target_location', 'strand', 'target', 'deletion', 'insertion', 'mismatch', 'seq_score', 'track_score', 'log2(track_score+1)','FDR', 'rank']
    df_output.to_csv(f'Trackseq_result_{outname}.csv', index=False)
    
    if args.clean:
        shutil.rmtree('./temp')

    return 'Done!'
    

if __name__ == '__main__' :
    result = main()
    print(result)


