
import pandas as pd
import numpy as np
import os, sys
sys.path.append( os.path.abspath(os.path.dirname(__file__)) )

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

def window_smooth(sr_smooth, window_size=3, times=1):
    window  = np.ones(window_size) / window_size

    bool_index = False
    if isinstance(sr_smooth, pd.Series):
        sr_index = sr_smooth.index 
        bool_index = True
    
    for i in range(times):
        sr_smooth = pd.Series(np.convolve(sr_smooth, window, mode='same'))
    
    if bool_index:
        sr_smooth.index = sr_index
    
    return sr_smooth

# 每 n 个数取平均
def segmental_mean(vector, n, drop='last'):
    # Identify the length and remainder
    length = len(vector)
    rem = length % n
    # If there is a remainder
    if rem != 0:
        if drop=='last':
            main_part = vector[:-rem]  # Part that could be reshaped
            #last_part = vector[-rem:]  # Excessive part in the end
            array = np.array(main_part).reshape(-1, n)
            result = array.mean(axis=1)
        elif drop=='first':
            main_part = vector[rem:]  # Part that could be reshaped
            #first_part = vector[:rem]  # Excessive part in the start
            array = np.array(main_part).reshape(-1, n)
            result = array.mean(axis=1)
    else:
        # If there's no remainder, proceed as usual
        array = np.array(vector).reshape(-1, n)
        result = array.mean(axis=1)

    return result

# v2.1 版本的计算信号长度存在一个问题：bw add 后，不是严格按 binsize 分割的，连续几个区域同数值，会被合并，但是这里都按 binsize 算，导致长度可能偏小
# v2.6 按 flank regions 分别取子集，避免根据 binsize 推测行数的误差，并且加入 trackseq v4 版本的正负bin占比计算
def target_signal(df_bdg_chr, chrom, cleavage_site, flank_max=100000, smooth_times = 1, window_size = 3, 
                  binsize=100, flank_regions=[500,1000,2000,5000], 
                  length_bkg = 20000, length_binsize=1000, length_min_noise=0.2, n_std=1, 
                  end='end',start='start',value='residual', pct_offset=0.0):
    # 输入数据必须是同一条染色体内的
    # 统计 flank regions 的个数
    # n_regions = len(flank_regions)
    ## 根据 binsize 计算每个 flank region 对应长度的 row 个数， 会有偏差
    ## flank_bins = [int(x/binsize) for x in flank_regions] # 取消这个算法
    
    assert flank_max >= max(flank_regions), 'flank_max must be larger than max(flank_regions)'
    assert length_binsize >= binsize, 'length_binsize must be larger than binsize'
    n_merge = int(length_binsize/binsize)
    n_bkg = int(length_bkg/length_binsize)

    # Left
    # 新版增加 list_signal_pct_L 会带 (pos_pct_left, left_pos_sum, left_neg_sum)*n_regions 外加一个 list_pct_score_L
    # list_signal_residual_L 数值和之前类似
    list_signal_pct_L = []
    list_pct_score_L = []
    list_signal_residual_L = []
    df_bdg_chr_L = df_bdg_chr[ (df_bdg_chr[end] >= cleavage_site-flank_max) & (df_bdg_chr[end]<=cleavage_site) ]
    if len(df_bdg_chr_L)<=window_size:
        L_length = 0
        L_overall_signal = 0
        for flank in flank_regions:
            list_signal_pct_L.extend([0,0,0])
            list_pct_score_L.append(0)
            list_signal_residual_L.append(0)
    else:
        ##################
        ## 先算 overall ##
        ##################
        L_length = 0
        # 信号平滑
        if smooth_times > 0:
           signal_residual_L = window_smooth(df_bdg_chr_L[value], window_size=window_size, times=smooth_times)
        else:
           signal_residual_L = df_bdg_chr_L[value]
        # 信号长度
        # 增大 binsize 减少随机波动，增大到 length_binsize 长度
        signal_residual_L_merged = segmental_mean(signal_residual_L, n_merge, drop='first')
        # 防止出现长度不足的情况，然后二次平滑
        if len(signal_residual_L_merged)<=window_size:
            L_length = 0
            L_overall_signal = 0
        else:
            signal_residual_L_merged = window_smooth(signal_residual_L_merged, window_size=3, times=3)
            # 计算背景和阈值
            bkg_L_mean = signal_residual_L_merged[:n_bkg].mean()
            bkg_L_std = max(length_min_noise, signal_residual_L_merged[:n_bkg].std())
            # 平移到均值为0
            # signal_residual_L_merged = signal_residual_L_merged - bkg_L_mean
            # 最后一个小于 threshold 的位点
            signal_start_index = signal_residual_L_merged.index[signal_residual_L_merged<bkg_L_mean+n_std*bkg_L_std].max()
            # 计算信号长度
            L_n_bins = signal_residual_L_merged.index.max()-signal_start_index
            if L_n_bins == 0:
                L_length = 0
            else:
                df_bdg_chr_L_good = df_bdg_chr_L[-n_merge*L_n_bins:]
                L_length = df_bdg_chr_L_good[end].iloc[-1]-df_bdg_chr_L_good[start].iloc[0]
            # 计算 overall 信号强度
            L_overall_signal = signal_residual_L_merged.sum()

        ###################
        ## 再算 proximal ##
        ###################
        # left_region_sum_norm 应该约等于 v2.5 以前的单一数值
        for flank in flank_regions:
            bool_flank = (df_bdg_chr_L[end] >= cleavage_site-flank)
            df_bdg_chr_L_flank = df_bdg_chr_L[ bool_flank ]
            signal_residual_L_flank = signal_residual_L[ bool_flank ]
            if df_bdg_chr_L_flank.empty:
                list_signal_pct_L.extend( [0,0,0] )
                list_pct_score_L.append(0)
                list_signal_residual_L.append(0)
                continue
            # pos and neg
            df_bdg_chr_L_flank_pos = df_bdg_chr_L_flank[df_bdg_chr_L_flank[value] > 0]
            df_bdg_chr_L_flank_neg = df_bdg_chr_L_flank[df_bdg_chr_L_flank[value] <= 0]
            n_pos_left = len(df_bdg_chr_L_flank_pos)
            n_neg_left = len(df_bdg_chr_L_flank_neg)
            # avoid zero
            if n_pos_left == 0:
                pos_pct_left = 0
            else:
                pos_pct_left = n_pos_left/(n_pos_left+n_neg_left)
            # pos/neg value sum
            left_pos_sum = df_bdg_chr_L_flank_pos[value].sum()
            left_neg_sum = df_bdg_chr_L_flank_neg[value].sum()
            list_signal_pct_L.extend( [pos_pct_left,left_pos_sum,left_neg_sum] ) 
            # 平滑 sum
            left_region_sum_norm = 1000*signal_residual_L_flank.sum()/flank
            list_signal_residual_L.append(left_region_sum_norm)            
            # pct_score
            left_pct_score = left_region_sum_norm*max(0,(pos_pct_left-pct_offset))
            list_pct_score_L.append(left_pct_score)

    # Right
    list_signal_pct_R = []
    list_pct_score_R = []
    list_signal_residual_R = []
    df_bdg_chr_R = df_bdg_chr[ (df_bdg_chr[start] <= cleavage_site+flank_max) & (df_bdg_chr[start]>=cleavage_site) ].copy()
    if len(df_bdg_chr_R)<=window_size:
        R_length = 0
        R_overall_signal = 0
        for flank in flank_regions:
            list_signal_pct_R.extend([0,0,0])
            list_pct_score_R.append(0)
            list_signal_residual_R.append(0)
    else:
        ##################
        ## 先算 overall ##
        ##################
        R_length = 0
        # 右侧信号反向
        df_bdg_chr_R[value] = -df_bdg_chr_R[value]
        # 信号平滑
        if smooth_times > 0:
           signal_residual_R = window_smooth(df_bdg_chr_R[value], window_size=window_size, times=smooth_times)
        else:
           signal_residual_R = df_bdg_chr_R[value]
        # 信号长度
        # 增大 binsize 减少随机波动，增大到 length_binsize 长度
        signal_residual_R_merged = segmental_mean(signal_residual_R, n_merge, drop='last')
        # 防止出现长度不足的情况
        if len(signal_residual_R_merged)<=window_size:
            R_length = 0
            R_overall_signal = 0
        else:
            signal_residual_R_merged = window_smooth(signal_residual_R_merged, window_size=3, times=3)
            # 计算背景和阈值
            bkg_R_mean = signal_residual_R_merged[-n_bkg:].mean()
            bkg_R_std = max(length_min_noise, signal_residual_R_merged[-n_bkg:].std())
            # 平移到均值为0
            # signal_residual_R_merged = signal_residual_R_merged - bkg_R_mean
            # 第一个小于 threshold 的位点
            signal_end_index = signal_residual_R_merged.index[signal_residual_R_merged<bkg_R_mean+n_std*bkg_R_std].min()
            # 计算信号长度
            R_n_bins = signal_end_index
            if R_n_bins == 0:
                R_length = 0
            else:
                df_bdg_chr_R_good = df_bdg_chr_R[:n_merge*R_n_bins]
                R_length = df_bdg_chr_R_good[end].iloc[-1]-df_bdg_chr_R_good[start].iloc[0]
            # 计算 overall 信号强度
            R_overall_signal = signal_residual_R_merged.sum()
        ###################
        ## 再算 proximal ##
        ###################
        # 注意，上面的 df_bdg_chr_R[value] 是反向的，因此这边还是 pos 为有信号
        for flank in flank_regions:
            bool_flank = (df_bdg_chr_R[start] <= cleavage_site+flank)
            df_bdg_chr_R_flank = df_bdg_chr_R[ bool_flank ]
            signal_residual_R_flank = signal_residual_R[ bool_flank ]
            if df_bdg_chr_R_flank.empty:
                list_signal_pct_R.extend( [0,0,0] )
                list_pct_score_R.append(0)
                list_signal_residual_R.append(0)
                continue
            # pos and neg
            df_bdg_chr_R_flank_pos = df_bdg_chr_R_flank[df_bdg_chr_R_flank[value] > 0]
            df_bdg_chr_R_flank_neg = df_bdg_chr_R_flank[df_bdg_chr_R_flank[value] <= 0]
            n_pos_right = len(df_bdg_chr_R_flank_pos)
            n_neg_right = len(df_bdg_chr_R_flank_neg)
            # avoid zero
            if n_pos_right == 0:
                pos_pct_right = 0
            else:
                pos_pct_right = n_pos_right/(n_pos_right+n_neg_right)
            # pos/neg value sum
            right_pos_sum = df_bdg_chr_R_flank_pos[value].sum()
            right_neg_sum = df_bdg_chr_R_flank_neg[value].sum()
            list_signal_pct_R.extend( [pos_pct_right,right_pos_sum,right_neg_sum] ) 
            # 平滑 sum
            right_region_sum_norm = 1000*signal_residual_R_flank.sum()/flank
            list_signal_residual_R.append(right_region_sum_norm)            
            # pct_score
            right_pct_score = right_region_sum_norm*max(0,(pos_pct_right-pct_offset))
            list_pct_score_R.append(right_pct_score)


    # calculate proximal_signal
    mean_signal_residual_L = np.mean(list_signal_residual_L)
    mean_signal_residual_R = np.mean(list_signal_residual_R)
    proximal_signal = mean_signal_residual_L+mean_signal_residual_R
    # calculate pct_score
    mean_pct_score_L = np.mean(list_pct_score_L)
    mean_pct_score_R = np.mean(list_pct_score_R)
    pct_score = mean_pct_score_L+mean_pct_score_R
    # calculate length and overall_signal
    signal_length = L_length + R_length
    #pct_signal_length = L_pct_length + R_pct_length
    # 有时候极远处有真编辑位点或者大噪音，会导致 overall_signal 的 bkg 不正确
    if L_overall_signal > 2*(mean_pct_score_L+mean_signal_residual_L):
        L_overall_signal = (mean_pct_score_L+mean_signal_residual_L)/2
    if R_overall_signal > 2*(mean_pct_score_R+mean_signal_residual_R):
        R_overall_signal = (mean_pct_score_R+mean_signal_residual_R)/2
    overall_signal = L_overall_signal + R_overall_signal
    list_return = list_signal_pct_L + list_signal_pct_R + \
                  list_pct_score_L + list_pct_score_R + \
                  list_signal_residual_L + list_signal_residual_R + \
                  [mean_signal_residual_L, mean_signal_residual_R] + \
                  [mean_pct_score_L, mean_pct_score_R] + \
                  [chrom+':'+str(cleavage_site)] + \
                  [L_length, R_length, L_overall_signal, R_overall_signal, signal_length, overall_signal, proximal_signal, pct_score]
                  # [L_pct_length, R_pct_length, pct_signal_length] 暂时不加这里了，额外写一个函数
                  # 2*3*n_regions
                  # 2*n_regions
                  # 2*n_regions
                  # 2+2+1+8

    return list_return

def target_signal_chunk(df_bdg_chr, df_alignment_chr, flank_max=100000, smooth_times = 1, window_size = 3, binsize=100, flank_regions=[500,1000,2000,5000], 
                        length_bkg = 20000, length_binsize=1000, length_min_noise=0.2, n_std=1, pct_offset=0.0):
    # 输入数据必须是同一条染色体内的
    list_target_all = []
    for a_row in df_alignment_chr.iterrows():
        chrom, cleavage_site = a_row[1]
        list_target = target_signal(df_bdg_chr, chrom, cleavage_site, flank_max, smooth_times = smooth_times, window_size = window_size, binsize=binsize, flank_regions=flank_regions, 
                                    length_bkg = length_bkg, length_binsize=length_binsize, length_min_noise=length_min_noise, n_std=n_std, pct_offset=pct_offset)
        list_target_all.append(list_target)
    df_result = pd.DataFrame(list_target_all)
    pct_features_L = [['L_pos_pct_'+x,'L_pos_'+x,'L_neg_'+x] for x in pd.Series(flank_regions).astype(str)]
    pct_features_L = [item for sublist in pct_features_L for item in sublist]
    pct_features_R = [['R_pos_pct_'+x,'R_pos_'+x,'R_neg_'+x] for x in pd.Series(flank_regions).astype(str)]
    pct_features_R = [item for sublist in pct_features_R for item in sublist]
    df_result.columns = pct_features_L + pct_features_R + \
                        list('L_pct_score_' + pd.Series(flank_regions).astype(str)) + list('R_pct_score_' + pd.Series(flank_regions).astype(str)) + \
                        list('L_' + pd.Series(flank_regions).astype(str)) + list('R_' + pd.Series(flank_regions).astype(str)) + \
                        ['L_mean', 'R_mean','L_mean_pct_score','R_mean_pct_score','chr_cleavage',
                         'L_length', 'R_length', 'L_overall_signal', 'R_overall_signal', 'signal_length', 'overall_signal','proximal_signal','pct_score']
    return df_result




# 2024.01.22. 额外写一个 signal length 算法，增加基于 pos_pct 而非 smooth 后的 overall_signal 的 length，叫 singal_length
def signal_length(df_bdg_chr, chrom, cleavage_site, end='end',start='start',value='residual', 
                  flank_max=100000, binsize=100):
    # 输入数据必须是同一条染色体内的
    # Left
    df_bdg_chr_L = df_bdg_chr[ (df_bdg_chr[end] >= cleavage_site-flank_max) & (df_bdg_chr[end]<=cleavage_site) ].copy()
    
    # pos and neg
    df_bdg_chr_L_flank_pos = df_bdg_chr_L_flank[df_bdg_chr_L_flank[value] > 0]
    df_bdg_chr_L_flank_neg = df_bdg_chr_L_flank[df_bdg_chr_L_flank[value] <= 0]
    n_pos_left = len(df_bdg_chr_L_flank_pos)
    n_neg_left = len(df_bdg_chr_L_flank_neg)
    # avoid zero
    if n_pos_left == 0:
        pos_pct_left = 0
    else:
        pos_pct_left = n_pos_left/(n_pos_left+n_neg_left)    
    
    
    df_bdg_chr_R = df_bdg_chr[ (df_bdg_chr[start] <= cleavage_site+flank_max) & (df_bdg_chr[start]>=cleavage_site) ].copy()
    # list_signal_residual_L 数值和之前类似
    list_signal_pct_L = []
    list_pct_score_L = []
    list_signal_residual_L = []
    

    return list_return