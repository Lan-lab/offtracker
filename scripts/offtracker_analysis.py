#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,glob,sys,time,shutil

if sys.version_info < (3,0):
    import platform
    raise Exception(f'python3 is needed, while running {platform.python_version()} now')

import offtracker
import offtracker.X_sequence as xseq
script_dir = os.path.abspath(os.path.dirname(offtracker.__file__))
script_folder= os.path.join(script_dir, 'mapping')

import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
from scipy.stats import norm

def main():
    parser = argparse.ArgumentParser()
    parser.description='Analyze the Tracking-seq data.'
    parser.add_argument('-f','--folder'  , type=str, required=True,    nargs='+', help='Directory of the data folder.' )
    parser.add_argument('--seqfolder'    , type=str, required=True,    help='folder containing df_candidate created by offtracker_cadidates.py.')
    parser.add_argument('--name'         , type=str, required=True,    help='custom name of the sgRNA' )
    parser.add_argument('--exp'          , type=str, default='all',    nargs='+', help='A substring mark in the name of experimental samples. The default is to use all samples other than control' )
    parser.add_argument('--control'      , type=str, default='none',   nargs='+', help='A substring mark in the name of control samples. The default is no control. "others" for all samples other than --exp.' )
    parser.add_argument('--fdr'          , type=float, default=0.05,     help='FDR threshold for the final result. Default is 0.05.')
    parser.add_argument('--score'        , type=float, default=1.9,        help='Track score threshold for the final result. Default is 1.9.')
    parser.add_argument('--smooth'       , type=int, default=1,        help='Smooth strength for the signal.')
    parser.add_argument('--window'       , type=int, default=3,        help='Window size for smoothing the signal.')
    parser.add_argument('--binsize'      , type=int, default=100,      help='Window size for smoothing the signal.')
    parser.add_argument('--flank_max'    , type=int, default=100000,   help='Maximun flanking distance from the candidate site.')
    parser.add_argument('--flank_regions', type=int, default=[1000,2000,3000,5000], nargs='+',help='flanking regions for calculating signal.')
    parser.add_argument('--SeqScorePower', type=float, default=4,      help='The seq score power' )
    parser.add_argument('--CtrClip'      , type=float, default=-0.5,     help='The lower clip for control samples' )
    parser.add_argument('-t','--thread'  , type=int, default=4,        help='Number of threads for parallel computing')
    parser.add_argument('-g','--genome'  , type=str, default='hg38',   help='File of chromosome sizes, or "hg38", "mm10" ')
    parser.add_argument('-o','--outdir'  , type=str, default='first',  help='The output folder. Default is the first folder of --folder' )
    parser.add_argument('--outname'      , type=str, default='same',   help='The suffix of output files. Default is the same --exp' )
    parser.add_argument('--signal_only'  , action='store_true', help='A developer option: stop before group analysis. ' )
    parser.add_argument('--overwrite'    , action='store_true', help='Whether to overwrite existed dataframes.' )
    parser.add_argument('--clean'        , action='store_true', help='Whether to remove temp files')


    args = parser.parse_args()

    print(f'Runing offtracker verision: {offtracker.__version__}')
    # main parameters
    folders = args.folder
    sgRNA_name = args.name
    pattern_exp = args.exp
    pattern_ctr = args.control
    fdr_thresh = args.fdr
    score_thresh = args.score
    binsize = args.binsize
    flank_max = args.flank_max
    flank_regions = args.flank_regions
    smooth_times = args.smooth
    window_size = args.window
    seq_score_power = args.SeqScorePower
    n_threads  = args.thread

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
    
    # load df_candidate
    try:
        df_candidate = pd.read_csv(os.path.join(args.seqfolder,f'df_candidate_{sgRNA_name}.csv'), index_col=0)
        df_candidate.index = df_candidate['target_location']
        df_candidate_brief = df_candidate[['chr','st','ed','best_strand','best_target','best_seq_score',
                                 'deletion', 'insertion','mismatch', 'GG', 
                                 'target_location', 'cleavage_site', 'ID_1','ID_2']]
        df_candidate_sub = df_candidate[['chr','cleavage_site']]
    except FileNotFoundError:
        return 'Please run offtracker_candidates.py first and provide the correct directory with --seqfolder'

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
    selected_sample_files = pd.concat([exp_sample_files,ctr_sample_files])
    selected_sample_names = pd.concat([exp_sample_names,ctr_sample_names])


    ##########################
    ## calculate the signal ##
    ##########################

    for a_file, a_name in zip(selected_sample_files, selected_sample_names):
        st = time.time()
        output = os.path.join(outdir, 'temp', a_name + f'.{sgRNA_name}.signal.csv')
        if (os.path.isfile(output))&(not args.overwrite):
            print(output, 'exists, skipped')
            continue
        df_bdg = xseq.read_bed(a_file)
        df_bdg.columns = ['chr','start','end','residual']
        # 将 df_bdg 按照染色体分组
        sample_groups = df_bdg.groupby('chr')
        # 2024.06.03. fix a bug that df_bdg has less chr than df_candidate
        total_chr = df_bdg['chr'].unique()
        df_candidate_sub_temp = df_candidate_sub[df_candidate_sub['chr'].isin(total_chr)]
        # 将 df_candidate_sub 按照染色体分组
        candidate_groups = df_candidate_sub_temp.groupby('chr')

        # 定义一个空的列表，用于存储每个染色体的数据
        chrom_list = []
        # 遍历分组后的数据
        list_index = []
        for chr_name, chr_candidate in candidate_groups:
            # 获取当前染色体对应的 df_sample 数据
            chr_sample = sample_groups.get_group(chr_name)
            # 保留 index
            list_index.extend(list(chr_candidate.index))
            # 将参数数据存储到列表中
            chrom_list.append([chr_sample, chr_candidate, flank_max, smooth_times, window_size, binsize, flank_regions])

        # 多线程运行
        with mp.Pool(n_threads) as p:
            signal_all = p.starmap(offtracker.target_signal_chunk, chrom_list)
        ed = time.time()
        print(f'{ed-st}s for {a_name} with {n_threads} threads')
        df_signal = pd.concat(signal_all)
        df_signal.index = list_index
        df_signal.to_csv(output)

    if args.signal_only:
        return 'signal_only is on, stop here.'
    
    ####################
    ## group analysis ##
    ####################
    if args.outname == 'same':
        if isinstance(pattern_exp, list):
            outname = '_'.join(pattern_exp)
        else:
            outname = pattern_exp
    else:
        outname = args.outname

    output = f'./temp/df_score_{outname}.csv'
    if (os.path.isfile(output))&(not args.overwrite):
        print(f'skip {output}')
        df_score =  pd.read_csv(output, index_col=0)
    else:
        signal_files = pd.Series(glob.glob( os.path.join(outdir, 'temp', f'*{sgRNA_name}.signal.csv') ))
        signal_names = signal_files.apply(os.path.basename).str.extract(rf'(.*)\.{sgRNA_name}\.signal\.csv',expand=False)

        # 读取并合并 samples
        list_df_exp_samples = []
        list_df_ctr_samples = []
        for a_file, a_name in zip(signal_files, signal_names):
            if a_name in exp_samples:
                df_temp = pd.read_csv(a_file, index_col=0)
                df_temp = df_temp.drop(['chr_cleavage'], axis=1)
                list_df_exp_samples.append(df_temp)
            elif a_name in ctr_samples:
                df_temp = pd.read_csv(a_file, index_col=0)
                df_temp = df_temp.drop(['chr_cleavage'], axis=1)
                list_df_ctr_samples.append(df_temp)
            else:
                pass
        
        # 计算每个组内的平均信号
        # 2023.12.07. exp 和 ctr 的信号分开展示
        df_score = df_candidate_brief.copy()
        df_exp = xseq.combine_df(list_df_exp_samples)
        if pattern_ctr != 'none':
            df_ctr = xseq.combine_df(list_df_ctr_samples)
            # 2023.12.10. 给 control 除了 'neg' 特征外的负数范围 clip，防止 exp-ctr 因此出现假阳性
            # 2023.12.31. 将 clip 范围由 -5 改为 -1
            # 2024.01.02. clip 模块移动到 filter and normalize
            # cols_clip = df_ctr.columns[~df_ctr.columns.str.contains('neg_')]
            # df_ctr[cols_clip] = df_ctr[cols_clip].clip(lower=-1)
            # df_exp[cols_clip] = df_exp[cols_clip].clip(lower=-1)
            # df_group_signal = df_exp - df_ctr
            df_exp.columns = 'exp_' + df_exp.columns
            df_ctr.columns = 'ctr_' + df_ctr.columns
            df_score = pd.concat([df_score, df_exp, df_ctr], axis=1)
        else:
            df_score = pd.concat([df_score, df_exp], axis=1)
        # 2024.06.03. 跑样例数据时，只有一个 chr6, 其他都是 nan, 不删除会导致后续计算出错
        df_score = df_score.dropna().copy()
        df_score.to_csv(output)
    
    ##########################
    ## filter and normalize ##
    ##########################
    output = f'./temp/df_result_{outname}.csv'
    if (os.path.isfile(output))&(not args.overwrite):
        print(f'skip {output} as the result exists')
        df_result =  pd.read_csv(output, index_col=0)
    else:
        if pattern_ctr != 'none':
            # 重算 proximal_signal 和 pct_score，因为 clip 了
            cols_exp_L = list('exp_L_' + pd.Series(flank_regions).astype(str))
            cols_exp_R = list('exp_R_' + pd.Series(flank_regions).astype(str))
            cols_ctr_L = list('ctr_L_' + pd.Series(flank_regions).astype(str))
            cols_ctr_R = list('ctr_R_' + pd.Series(flank_regions).astype(str))
            cols_exp_L_pct_score = list('exp_L_pct_score_' + pd.Series(flank_regions).astype(str))
            cols_exp_R_pct_score = list('exp_R_pct_score_' + pd.Series(flank_regions).astype(str))
            cols_ctr_L_pct_score = list('ctr_L_pct_score_' + pd.Series(flank_regions).astype(str))
            cols_ctr_R_pct_score = list('ctr_R_pct_score_' + pd.Series(flank_regions).astype(str))
            df_score['exp_L_mean'] = df_score[cols_exp_L].mean(axis=1)
            df_score['exp_R_mean'] = df_score[cols_exp_R].mean(axis=1)
            df_score['ctr_L_mean'] = df_score[cols_ctr_L].clip(lower=args.CtrClip).mean(axis=1)
            df_score['ctr_R_mean'] = df_score[cols_ctr_R].clip(lower=args.CtrClip).mean(axis=1)
            df_score['exp_L_mean_pct_score'] = df_score[cols_exp_L_pct_score].mean(axis=1)
            df_score['exp_R_mean_pct_score'] = df_score[cols_exp_R_pct_score].mean(axis=1)
            df_score['ctr_L_mean_pct_score'] = df_score[cols_ctr_L_pct_score].clip(lower=args.CtrClip).mean(axis=1)
            df_score['ctr_R_mean_pct_score'] = df_score[cols_ctr_R_pct_score].clip(lower=args.CtrClip).mean(axis=1)
            df_score['L_mean'] = df_score['exp_L_mean'] - df_score['ctr_L_mean']
            df_score['R_mean'] = df_score['exp_R_mean'] - df_score['ctr_R_mean']
            df_score['L_mean_pct_score'] = df_score['exp_L_mean_pct_score'] - df_score['ctr_L_mean_pct_score']
            df_score['R_mean_pct_score'] = df_score['exp_R_mean_pct_score'] - df_score['ctr_R_mean_pct_score']
            df_score['L_length'] = df_score['exp_L_length'] - df_score['ctr_L_length']
            df_score['R_length'] = df_score['exp_R_length'] - df_score['ctr_R_length']
            df_score['signal_length'] = df_score['L_length'] + df_score['R_length']
            df_score['proximal_signal'] = df_score['L_mean'] + df_score['R_mean']
            df_score['pct_score'] = df_score['L_mean_pct_score'] + df_score['R_mean_pct_score']

        # 整理表格
        mean_seq_score = round(df_score['best_seq_score'].mean(),3)
        df_score['norm_best_seq_score'] = np.power(df_score['best_seq_score']/mean_seq_score, seq_score_power)        
        df_score['final_score_1'] = df_score['proximal_signal']*df_score['norm_best_seq_score']
        df_score['final_score_2'] = df_score['pct_score']*df_score['norm_best_seq_score']
        #df_score['final_score_2'] = df_score[f'overall_signal']*df_score['norm_best_seq_score']
        df_score['raw_score'] = df_score['final_score_1'] + df_score['final_score_2']
        df_score = df_score.sort_values('raw_score', ascending=False)    

        # local dedup
        list_nondup = offtracker.dedup_two(df_score,'ID_1','ID_2')
        df_result = df_score[list_nondup].copy()

        # 标准化分布      
        target_std=0.15
        n_outliers = int(np.ceil(len(df_result)*0.01))
        score_bkg = df_result['raw_score'][n_outliers:-n_outliers]
        mean_score_bkg = score_bkg.mean()
        std_score_bkg = score_bkg.std()
        df_result['track_score'] = (df_result['raw_score'] - mean_score_bkg) / std_score_bkg
        df_result['track_score'] = df_result['track_score']*target_std + 1
        df_result = df_result.sort_values(by='track_score', ascending=False)
        df_result['log2_track_score'] = np.log2(df_result['track_score'].clip(lower=0.5))   

        # 单边信号周围有更高分的，去掉
        # v2.1 后 cols_L, cols_R 要手动
        # 2024.01.26. 只看 1kb 了，但这个办法还是无法解决约 100-500 bp 以内有两个相似位点的问题
        if pattern_ctr != 'none':
            cols_L = ['exp_L_1000']
            cols_R = ['exp_R_1000']
        else:
            cols_L = ['L_1000'] # df_score.columns[df_score.columns.str.contains('^L_\d+')]
            cols_R = ['R_1000'] # df_score.columns[df_score.columns.str.contains('^R_\d+')]
        seq_score_thresh = np.power(1.25, seq_score_power)
        search_distance = 100000
        candidate_dup = list(df_result[((df_result[cols_R].max(axis=1)<=0)|(df_result[cols_L].max(axis=1)<=0))&(df_result['log2_track_score']>0.8)].index)
        list_dedup = []
        for a_loc in candidate_dup:
            temp_chr = df_result.loc[a_loc,'chr']
            # 如果序列特别像就不过滤
            temp_seq_score = df_result.loc[a_loc,'norm_best_seq_score']
            if temp_seq_score > seq_score_thresh:
                    continue
            # 取排其前且同chr者
            temp_df_result = df_result.loc[:a_loc]
            temp_df_result = temp_df_result[temp_df_result['chr'] == temp_chr].iloc[:-1]
            if len(temp_df_result)==0:
                continue
            else:
                # 距离小于 search_distance 者记为信号溢出假阳性
                if (temp_df_result['cleavage_site']-df_result.loc[a_loc,'cleavage_site']).abs().min()<search_distance :
                    list_dedup.append(a_loc)
        # 去除重复
        df_result = df_result[~df_result.index.isin(list_dedup)].copy()
        # print(f'filter {len(list_dedup)} sites')

        # fitting normal distribution
        score_for_fitting = df_result['log2_track_score'][n_outliers:-n_outliers]
        mu, std = norm.fit(score_for_fitting) 
        print('mean_score:{:.3f};std:{:.3f}'.format(mu,std))
        # pv and fdr
        df_result['pv'] = df_result['log2_track_score'].apply( lambda x: norm.sf(x,loc=mu,scale=std) )
        df_result['pv'] = df_result['pv'].clip(lower=1e-320)
        df_result['fdr'] = offtracker.fdr(df_result['pv'])
        df_result['rank'] = range(1,len(df_result)+1)
        df_result.to_csv(output)

    output = f'Offtracker_result_{outname}.csv'
    if (os.path.isfile(output))&(not args.overwrite):
        print(f'skip {output} as the result exists')
    else:
        # 2024.06.03. 以防 fdr<=fdr_thresh 滤掉了 track_score>=2 的位点
        bool_fdr = df_result['fdr']<=fdr_thresh
        bool_score = df_result['track_score']>=score_thresh
        # 2025.06.05. BE可能会形成单边信号，导致 track_score 为负数，也保留
        bool_neg_score = df_result['track_score']< -1
        df_output = df_result[bool_fdr|bool_score|bool_neg_score].copy()
        if pattern_ctr != 'none':
            df_output = df_output[['target_location', 'best_strand','best_target','deletion','insertion','mismatch',
                                'exp_L_length', 'exp_R_length','ctr_L_length','ctr_R_length','L_length','R_length','signal_length',
                                'norm_best_seq_score','track_score', 'log2_track_score','fdr','rank']]
            df_output.columns = ['target_location', 'strand', 'target', 'deletion', 'insertion', 'mismatch', 
                                'exp_L_length', 'exp_R_length','ctr_L_length','ctr_R_length','L_length','R_length','signal_length',
                                'seq_score', 'track_score', 'log2_track_score','FDR', 'rank']
        else:
            df_output = df_output[['target_location', 'best_strand','best_target','deletion','insertion','mismatch',
                                'L_length', 'R_length','signal_length',
                                'norm_best_seq_score','track_score', 'log2_track_score','fdr','rank']]
            df_output.columns = ['target_location', 'strand', 'target', 'deletion', 'insertion', 'mismatch', 
                                'L_length', 'R_length','signal_length',
                                'seq_score', 'track_score', 'log2_track_score','FDR', 'rank']
        df_output.to_csv(f'Offtracker_result_{outname}.csv', index=False)

        if args.clean:
            shutil.rmtree('./temp')

    return 'Done!'
    

if __name__ == '__main__' :
    result = main()
    print(result)


