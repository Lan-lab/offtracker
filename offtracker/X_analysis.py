
import pandas as pd
import numpy as np
import os, sys, pybedtools
sys.path.append( os.path.abspath(os.path.dirname(__file__)) )
from X_general import *

def signal_formula1(signal, nonsignal, pseudocount):
    # 背景均值略小于 1
    assert pseudocount>0
    return signal/(nonsignal+pseudocount)

def signal_formula2(signal, nonsignal, pseudocount):
    # 背景均值为 1
    assert pseudocount>0
    return (signal+pseudocount)/(nonsignal+pseudocount)

def signal_formula3(signal, nonsignal, pseudocount):
    # 调整背景均值为 1
    assert pseudocount>0
    out = (signal - nonsignal + 1)
    out.loc[out<0] = 0
    return out

def signal_formula4(signal, nonsignal, pseudocount):
    # 背景均值为 1
    assert pseudocount>0
    return signal

def shape_formula1(normed_left_signal, normed_left_nonsignal, normed_right_signal, normed_right_nonsignal, bkg_noise, ratio_thresh=1.5, exp_weight = 0.3):
    # 防止破坏性的极小值出现
    normed_left_signal = normed_left_signal.clip(lower=1)
    normed_left_nonsignal = normed_left_nonsignal.clip(lower=1)
    normed_right_signal = normed_right_signal.clip(lower=1)
    normed_right_nonsignal = normed_right_nonsignal.clip(lower=1) 
    left_ratio = normed_left_signal/normed_left_nonsignal
    right_ratio = normed_right_signal/normed_right_nonsignal
    #
    good_shape = pd.concat([left_ratio,right_ratio],axis=1).min(axis=1)>ratio_thresh
    left_ratio.loc[~good_shape]=1
    right_ratio.loc[~good_shape]=1
    bonus_coef = left_ratio*right_ratio
    bonus_coef = np.power(bonus_coef,exp_weight)
    return bonus_coef

def shape_formula2(normed_left_signal, normed_left_nonsignal, normed_right_signal, normed_right_nonsignal, bkg_noise, ratio_thresh=1, exp_weight = 0.3):
    left_residual = normed_left_signal - normed_left_nonsignal
    right_residual = normed_right_signal - normed_right_nonsignal
    # 这里先取 max，也可以取 min 或者 mean
    good_left = left_residual>ratio_thresh
    good_right = right_residual>ratio_thresh
    good_shape = good_left&good_right
    left_residual.loc[~good_shape]=0
    right_residual.loc[~good_shape]=0
    bonus_coef = (left_residual+right_residual)/(2*ratio_thresh)
    bonus_coef = np.power(bonus_coef,exp_weight)
    bonus_coef = bonus_coef.clip(lower=1)
    return bonus_coef

def score_formula1(exp, ctrl, pseudocount):
    assert pseudocount>0
    return exp/(max(ctrl,0)+pseudocount)

def score_formula2(exp, ctrl, pseudocount):
    assert pseudocount>0
    return (exp+pseudocount)/(max(ctrl,0)+pseudocount)

def score_formula3(exp, ctrl, pseudocount):
    assert pseudocount>0
    return exp/max(ctrl, pseudocount)

def score_formula4(exp, ctrl, pseudocount):
    assert pseudocount>0
    return max((exp - ctrl), 0)

def score_formula5(exp, ctrl, pseudocount):
    assert pseudocount>0
    return exp

def fdr(p_vals):
    # Benjamini-Hochberg
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr

def dedup_two( df_loc, col_ID_1='ID_1', col_ID_2='ID_2'):
    # 会根据 df_loc 的排序保留第一个 location
    # dedup 结束后，剩下的 ID_1 + ID_2 并集可能会小于 dedup 前的并集
    list_nondup = []
    set_IDs = set()
    df_IDs = df_loc[[col_ID_1,col_ID_2]]
    for a_row in df_IDs.iterrows():
        temp = a_row[1]
        if (temp[col_ID_1] in set_IDs) or (temp[col_ID_2] in set_IDs):
            # 只要有一ID出现过，即便另一ID没出现过，也不更新 set_IDs
            list_nondup.append(False)
        else:
            set_IDs.add(temp[col_ID_1])
            set_IDs.add(temp[col_ID_2])
            list_nondup.append(True)
    return list_nondup

def cand_count(a_pref, the_sgRNA, regions, seq_folder, dir_chrom_sizes, overwrite=False):
    # 
    forward_bed = f'{a_pref}.fw.bed'
    reverse_bed = f'{a_pref}.rv.bed'
    # put into temp_dir
    dirname = os.path.dirname(a_pref)
    basename = os.path.basename(a_pref)
    temp_dir = os.path.join(dirname, 'temp')
    if not os.path.exists( temp_dir ):
        os.makedirs(temp_dir)
    a_pref = os.path.join(temp_dir, basename)
    #
    for a_region in regions:
        work = 1
        if os.path.isfile(f'{a_pref}_candidate_{the_sgRNA}_right_rstrand_{a_region}.count'):
            #print(f'Candidates for {the_sgRNA} within {a_region} of {basename} exists.\n')
            work = 0
        if overwrite:
            print('overwrite mode')
            work = 1
        if work == 1:
            print(f'Working for {basename} within {a_region} bp.\n')
            # 左侧区域
            left_region=os.path.join(seq_folder, f'{the_sgRNA}_candidate_left_{a_region}.bed')
            a = pybedtools.BedTool(left_region)
            b = pybedtools.BedTool(forward_bed)
            c = a.coverage(b,sorted=True,g=dir_chrom_sizes)
            c.saveas(f'{a_pref}_candidate_{the_sgRNA}_left_fstrand_{a_region}.count')
            b = pybedtools.BedTool(reverse_bed)
            c = a.coverage(b,sorted=True,g=dir_chrom_sizes)
            c.saveas(f'{a_pref}_candidate_{the_sgRNA}_left_rstrand_{a_region}.count')
            # 右侧区域
            right_region=os.path.join(seq_folder, f'{the_sgRNA}_candidate_right_{a_region}.bed')
            a = pybedtools.BedTool(right_region)
            b = pybedtools.BedTool(forward_bed)
            c = a.coverage(b,sorted=True,g=dir_chrom_sizes)
            c.saveas(f'{a_pref}_candidate_{the_sgRNA}_right_fstrand_{a_region}.count')
            b = pybedtools.BedTool(reverse_bed)
            c = a.coverage(b,sorted=True,g=dir_chrom_sizes)
            c.saveas(f'{a_pref}_candidate_{the_sgRNA}_right_rstrand_{a_region}.count')

def bkg_count(a_pref, the_sgRNA, bkgs, seq_folder, dir_chrom_sizes, overwrite=False):
    # no need for LF and RR
    forward_bed = f'{a_pref}.fw.bed'
    reverse_bed = f'{a_pref}.rv.bed'
    # put into temp_dir
    dirname = os.path.dirname(a_pref)
    basename = os.path.basename(a_pref)
    temp_dir = os.path.join(dirname, 'temp')
    if not os.path.exists( temp_dir ):
        os.makedirs(temp_dir)
    a_pref = os.path.join(temp_dir, basename)
    #
    work = 1
    for a_region in bkgs:
        if os.path.isfile(f'{a_pref}_candidate_{the_sgRNA}_right_fstrand_{a_region}.count'):
            #print(f'Candidates for {the_sgRNA} within {a_region} of {basename} exists.\n')
            work = 0
        if overwrite:
            print('overwrite mode')
            work = 1
        if work == 1:
            print(f'Working for {basename} within {a_region} bp.\n')
            # 左侧区域
            left_region=os.path.join(seq_folder, f'{the_sgRNA}_candidate_left_{a_region}.bed')
            a = pybedtools.BedTool(left_region)
            b = pybedtools.BedTool(reverse_bed)
            c = a.coverage(b,sorted=True,g=dir_chrom_sizes)
            c.saveas(f'{a_pref}_candidate_{the_sgRNA}_left_rstrand_{a_region}.count')
            # 右侧区域
            right_region=os.path.join(seq_folder, f'{the_sgRNA}_candidate_right_{a_region}.bed')
            a = pybedtools.BedTool(right_region)
            b = pybedtools.BedTool(forward_bed)
            c = a.coverage(b,sorted=True,g=dir_chrom_sizes)
            c.saveas(f'{a_pref}_candidate_{the_sgRNA}_right_fstrand_{a_region}.count')


def load_count(    a_sample, regions, bkgs, sgRNA_name, signal_formula, noise_length, 
    ratio_thresh, exp_weight, shape_formula=None,
    operator='p', region_op='mean', dirname='./',  pseudo_count=1
    ):
    list_noise = []
    for a_bkg in bkgs:
        left_count_rstrand = os.path.join(dirname, f'{a_sample}_candidate_{sgRNA_name}_left_rstrand_{a_bkg}.count')
        left_count_rstrand = readbed(left_count_rstrand)
        right_count_fstrand= os.path.join(dirname, f'{a_sample}_candidate_{sgRNA_name}_right_fstrand_{a_bkg}.count')
        right_count_fstrand = readbed(right_count_fstrand)
        list_noise.append(right_count_fstrand[3])
        list_noise.append(left_count_rstrand[3])
    df_noise = pd.concat(list_noise,axis=1, ignore_index=True)
    df_noise['mean'] = df_noise.mean(axis=1)
    mean_all = df_noise['mean'].mean()
    df_noise['std'] = df_noise.std(axis=1)
    df_noise['outlier'] = df_noise['mean']+2*df_noise['std']
    for i in range(len(bkgs)*2):
        df_noise.loc[df_noise[i]>df_noise['outlier'],i] = np.nan
    df_noise['mean2'] = df_noise.iloc[:,:6].mean(axis=1, skipna = True)
    n_0bkg = sum(df_noise.mean2==0)
    if n_0bkg > 0:
        print(f'{n_0bkg} region(s) with 0 count in background.')
    # 由于有些位置可能出现无法 mapping 而产生大量空白，导致局部噪音过低，因此这里主要是防局部高噪音造成假阳性
    df_noise.loc[df_noise['mean2']<mean_all, 'mean2'] = mean_all
    df_noise['noise_bp'] = df_noise['mean2']/noise_length
    noise_5kb = df_noise['noise_bp'].mean()*5000
    print('Average noise within 5kb on a single strand: {:.2f}'.format(noise_5kb))
    if noise_5kb < 10:
        print('The sequencing depth might be too shallow')
    list_df_counts = []
    for a_region in regions:
        left_count_fstrand= os.path.join(dirname, f'{a_sample}_candidate_{sgRNA_name}_left_fstrand_{a_region}.count')
        left_count_fstrand = readbed(left_count_fstrand)
        left_count_fstrand.columns=[f'chr_{a_region}',f'st_left_{a_region}',f'ed_left_{a_region}',f'counts_left_F_{a_region}',f'cover_left_F_bp_{a_region}',f'length_{a_region}',f'cover_left_F_pct_{a_region}']
        right_count_fstrand=os.path.join(dirname, f'{a_sample}_candidate_{sgRNA_name}_right_fstrand_{a_region}.count')
        right_count_fstrand = readbed(right_count_fstrand)
        right_count_fstrand.columns=[f'chr_{a_region}',f'st_right_{a_region}',f'ed_right_{a_region}',f'counts_right_F_{a_region}',f'cover_right_F_bp_{a_region}',f'length_{a_region}',f'cover_right_F_pct_{a_region}']
        left_count_rstrand= os.path.join(dirname, f'{a_sample}_candidate_{sgRNA_name}_left_rstrand_{a_region}.count')
        left_count_rstrand = readbed(left_count_rstrand)
        left_count_rstrand.columns=[f'chr_{a_region}',f'st_left_{a_region}',f'ed_left_{a_region}',f'counts_left_R_{a_region}',f'cover_left_R_bp_{a_region}',f'length_{a_region}',f'cover_left_R_pct_{a_region}']
        right_count_rstrand=os.path.join(dirname, f'{a_sample}_candidate_{sgRNA_name}_right_rstrand_{a_region}.count')
        right_count_rstrand = readbed(right_count_rstrand)
        right_count_rstrand.columns=[f'chr_{a_region}',f'st_right_{a_region}',f'ed_right_{a_region}',f'counts_right_R_{a_region}',f'cover_right_R_bp_{a_region}',f'length_{a_region}',f'cover_right_R_pct_{a_region}']
        
        df_counts = pd.concat([left_count_fstrand[[f'chr_{a_region}',f'st_left_{a_region}', f'ed_left_{a_region}',
                                                f'counts_left_F_{a_region}',f'cover_left_F_bp_{a_region}',f'cover_left_F_pct_{a_region}']],
                            left_count_rstrand[[f'counts_left_R_{a_region}',f'cover_left_R_bp_{a_region}',f'cover_left_R_pct_{a_region}']],
                            right_count_fstrand[[f'ed_right_{a_region}',f'counts_right_F_{a_region}',f'cover_right_F_bp_{a_region}',f'cover_right_F_pct_{a_region}']],
                            right_count_rstrand[[f'counts_right_R_{a_region}',f'cover_right_R_bp_{a_region}',f'cover_right_R_pct_{a_region}']]],axis=1)
        
        df_counts = df_counts.reindex(columns=[f'chr_{a_region}', f'st_left_{a_region}', f'ed_right_{a_region}', f'ed_left_{a_region}',
                                            f'counts_left_F_{a_region}',f'cover_left_F_bp_{a_region}',f'cover_left_F_pct_{a_region}',
                                            f'counts_left_R_{a_region}',f'cover_left_R_bp_{a_region}',f'cover_left_R_pct_{a_region}',
                                            f'counts_right_F_{a_region}',f'cover_right_F_bp_{a_region}',f'cover_right_F_pct_{a_region}',
                                            f'counts_right_R_{a_region}',f'cover_right_R_bp_{a_region}',f'cover_right_R_pct_{a_region}'])
        df_counts.columns = [f'chr_{a_region}', f'st_{a_region}', f'ed_{a_region}',  f'midpoint_{a_region}',
                            f'counts_left_F_{a_region}',f'cover_left_F_bp_{a_region}',f'cover_left_F_pct_{a_region}',
                            f'counts_left_R_{a_region}',f'cover_left_R_bp_{a_region}',f'cover_left_R_pct_{a_region}',
                            f'counts_right_F_{a_region}',f'cover_right_F_bp_{a_region}',f'cover_right_F_pct_{a_region}',
                            f'counts_right_R_{a_region}',f'cover_right_R_bp_{a_region}',f'cover_right_R_pct_{a_region}']
    
        # signal enrichment = Cs/(Cn+B)
        bkg_noise = df_noise['noise_bp']*a_region
        normed_left_signal = df_counts[f'counts_left_F_{a_region}']/bkg_noise
        normed_left_nonsignal = df_counts[f'counts_left_R_{a_region}']/bkg_noise
        normed_right_signal = df_counts[f'counts_right_R_{a_region}']/bkg_noise
        normed_right_nonsignal = df_counts[f'counts_right_F_{a_region}']/bkg_noise
        df_counts[f'N_LF_{a_region}'] = normed_left_signal
        df_counts[f'N_LR_{a_region}'] = normed_left_nonsignal
        df_counts[f'N_RR_{a_region}'] = normed_right_signal
        df_counts[f'N_RF_{a_region}'] = normed_right_nonsignal
        # 可变公式区
        df_counts[f'left_signal_{a_region}'] = signal_formula(normed_left_signal, normed_left_nonsignal, pseudo_count)
        df_counts[f'right_signal_{a_region}'] = signal_formula(normed_right_signal, normed_right_nonsignal, pseudo_count)
        if shape_formula:
            df_counts[f'shape_bonus_{a_region}'] = shape_formula(normed_left_signal, normed_left_nonsignal,
                                                                 normed_right_signal, normed_right_nonsignal, 
                                                                 bkg_noise, ratio_thresh, exp_weight)
        else:
            df_counts[f'shape_bonus_{a_region}'] = 1
        list_df_counts.append(df_counts)
    df_counts = pd.concat(list_df_counts,axis=1)
    left_signal_cols = df_counts.columns[ df_counts.columns.str.contains('left_signal_') ]
    right_signal_cols = df_counts.columns[ df_counts.columns.str.contains('right_signal_') ]
    
    print('region_op',region_op)
    if region_op == 'mean':
        df_counts['left_signal'] = df_counts[left_signal_cols].mean(axis=1)
        df_counts['right_signal'] = df_counts[right_signal_cols].mean(axis=1)        
    elif region_op == 'max':
        df_counts['left_signal'] = df_counts[left_signal_cols].max(axis=1)
        df_counts['right_signal'] = df_counts[right_signal_cols].max(axis=1)
    elif region_op == 'min':
        df_counts['left_signal'] = df_counts[left_signal_cols].min(axis=1)
        df_counts['right_signal'] = df_counts[right_signal_cols].min(axis=1)    
    else:
        raise Exception('region_op should be "mean", "max", or "min" ')
    
    max_region = max(regions)
    df_counts['chr'] = df_counts[f'chr_{max_region}']
    df_counts['st'] = df_counts[f'st_{max_region}']
    df_counts['ed'] = df_counts[f'ed_{max_region}']
    df_counts['midpoint'] = df_counts[f'midpoint_{max_region}']
    df_counts = df_counts.reindex( columns= ['chr','st','ed','midpoint'] + list(df_counts.columns) )
    df_counts = df_counts.loc[:,~df_counts.columns.duplicated()].copy()
    print('Raw regions:', len(df_counts))
    
    # operator 选择
    if operator == 'p':
        df_counts['score'] = (df_counts['left_signal']+df_counts['right_signal'])/2
    elif operator == 'm':
        df_counts['score'] = np.power(df_counts['left_signal']*df_counts['right_signal'],0.5)
    
    #
    min_region = min(regions)
    
    if shape_formula:
        # 
        shape_bonus_cols = df_counts.columns[ df_counts.columns.str.contains('shape_bonus_') ]
        df_counts['shape_bonus'] = df_counts[shape_bonus_cols].mean(axis=1)

        
    df_counts = df_counts.sort_values(by='score',ascending=False).reset_index(drop=True)
    
    ### 一些其他特征
    
    # 左右最小范围信号绝对强度比值与差值
    max_adjacent_ratio = df_counts[[f'left_signal_{min_region}',f'right_signal_{min_region}']].max(axis=1)
    min_adjacent_ratio = df_counts[[f'left_signal_{min_region}',f'right_signal_{min_region}']].min(axis=1)
    df_counts['signal_FC'] = max_adjacent_ratio/(min_adjacent_ratio+0.001)

    # 左右500bp小信号边若为负数，可能是假的
    df_counts['left_signal_residual'] = df_counts[f'N_LF_{min_region}'] - df_counts[f'N_LR_{min_region}']
    df_counts['right_signal_residual'] = df_counts[f'N_RR_{min_region}'] - df_counts[f'N_RF_{min_region}']
    min_adjacent_signal = df_counts[['left_signal_residual','right_signal_residual']].min(axis=1)
    df_counts['signal_min'] = min_adjacent_signal

    # 具体位置去重
    df_counts['location'] = igvfmt(df_counts)
    df_counts = df_counts.drop_duplicates(subset='location').reset_index(drop=True).copy()
    
    # 第二版 unique_ID 
    point_head = (df_counts['midpoint']/1000).astype(int)
    df_counts['ID_1'] = df_counts['chr'] + ':' + point_head.astype(str)
    point_tail = df_counts['midpoint'] % 1000
    df_counts.loc[point_tail<500,'ID_2'] = df_counts['chr'] + ':' + (point_head-1).astype(str)
    df_counts.loc[point_tail>=500,'ID_2'] = df_counts['chr'] + ':' + (point_head+1).astype(str)
    
    return df_counts, df_noise

