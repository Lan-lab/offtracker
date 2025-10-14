
import pandas as pd
import polars as pl
import numpy as np
import os, sys, re
import offtracker.X_sequence as xseq
sys.path.append( os.path.abspath(os.path.dirname(__file__)) )

def fdr(p_vals):
    # Benjamini-Hochberg
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr_value = p_vals * len(p_vals) / ranked_p_values
    fdr_value[fdr_value > 1] = 1
    return fdr_value


def mark_regions_single_chr(dp, min_distance=1000):
    unique_chr = dp['chr'].unique()
    assert len(unique_chr) == 1
    unique_chr = unique_chr[0]

    # Initialize variables for marking regions
    region_id = 1
    current_start = None
    current_end = None
    marked_regions = []

    for row in dp.iter_rows(named=True):
        start, end = row['st'], row['ed']

        if current_start is None:
            # First region
            current_start = start
            current_end = end
            marked_regions.append(f'{unique_chr}_region_{region_id}')
        else:
            if start <= current_end + min_distance:
                # Mark as the same region
                marked_regions.append(f'{unique_chr}_region_{region_id}')
            else:
                # New region
                region_id += 1
                marked_regions.append(f'{unique_chr}_region_{region_id}')
                current_start = start
                current_end = end

        current_end = max(current_end, end)

    return dp.with_columns(region_index=pl.Series(marked_regions))





# def dedup_two( df_loc, col_ID_1='ID_1', col_ID_2='ID_2'):
#     # 会根据 df_loc 的排序保留第一个 location
#     # dedup 结束后，剩下的 ID_1 + ID_2 并集可能会小于 dedup 前的并集
#     list_nondup = []
#     set_IDs = set()
#     df_IDs = df_loc[[col_ID_1,col_ID_2]]
#     for a_row in df_IDs.iterrows():
#         temp = a_row[1]
#         if (temp[col_ID_1] in set_IDs) or (temp[col_ID_2] in set_IDs):
#             # 只要有一ID出现过，即便另一ID没出现过，也不更新 set_IDs
#             list_nondup.append(False)
#         else:
#             set_IDs.add(temp[col_ID_1])
#             set_IDs.add(temp[col_ID_2])
#             list_nondup.append(True)
#     return list_nondup

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
    """_summary_

    Args:
        df_bdg_chr (_type_): .bdg table with the same chromosome
        chrom (_type_): chr name
        cleavage_site (_type_): cleavage site
        flank_max (int, optional): _description_. Defaults to 100000.
        smooth_times (int, optional): _description_. Defaults to 1.
        window_size (int, optional): _description_. Defaults to 3.
        binsize (int, optional): _description_. Defaults to 100.
        flank_regions (list, optional): _description_. Defaults to [500,1000,2000,5000].
        length_bkg (int, optional): _description_. Defaults to 20000.
        length_binsize (int, optional): _description_. Defaults to 1000.
        length_min_noise (float, optional): _description_. Defaults to 0.2.
        n_std (int, optional): _description_. Defaults to 1.
        end (str, optional): _description_. Defaults to 'end'.
        start (str, optional): _description_. Defaults to 'start'.
        value (str, optional): _description_. Defaults to 'residual'.
        pct_offset (float, optional): _description_. Defaults to 0.0.

    Returns:
        _type_: _description_
    """
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
    """

    Args:
        df_bdg_chr (_type_): .bdg table with the same chromosome
        df_alignment_chr (_type_): candidate sites
        flank_max (int, optional): _description_. Defaults to 100000.
        smooth_times (int, optional): _description_. Defaults to 1.
        window_size (int, optional): _description_. Defaults to 3.
        binsize (int, optional): _description_. Defaults to 100.
        flank_regions (list, optional): _description_. Defaults to [500,1000,2000,5000].
        length_bkg (int, optional): _description_. Defaults to 20000.
        length_binsize (int, optional): _description_. Defaults to 1000.
        length_min_noise (float, optional): _description_. Defaults to 0.2.
        n_std (int, optional): _description_. Defaults to 1.
        pct_offset (float, optional): _description_. Defaults to 0.0.

    Returns:
        _type_: _description_
    """
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











############################################################################
# 2025.08.08 新写的的基于 Bio.Align 的 local realign，用于局部修正脱靶位点坐标 #
############################################################################
def create_substitution_matrix(mismatch_score=0.01):
    """
    Create substitution matrix for DNA alignment using Bio.Align format.
    """
    alphabet = 'ATGCN'
    matrix = np.full((len(alphabet), len(alphabet)), mismatch_score)
    
    # Set match scores
    for i in range(len(alphabet)):
        matrix[i][i] = 2.0
    
    # N matches with everything
    n_idx = alphabet.index('N')
    matrix[n_idx, :] = 2.0
    matrix[:, n_idx] = 2.0
    
    return matrix, alphabet

def sgRNA_alignment_new(a_key, sgRNA, seq, substitution_matrix=None, alphabet=None, 
                   mismatch_score=0.01):
    """
    Perform local alignment using Bio.Align instead of deprecated pairwise2.
    """
    from Bio import Align
    if substitution_matrix is None or alphabet is None:
        substitution_matrix, alphabet = create_substitution_matrix(mismatch_score)
    
    # Create aligner
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = Align.substitution_matrices.Array(
        alphabet=alphabet, dims=2, data=substitution_matrix
    )
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2
    aligner.mode = 'local'
    
    try:
        # Perform alignment
        alignments = aligner.align(sgRNA, seq)

        if not alignments:
            # No alignment found, return default values
            return [0, 0, '', f"{a_key.split(':')[0]}:0-0", 0, 0, len(sgRNA)]
        
        # Convert to list for indexing
        alignments = list(alignments)
         
        # Extract alignment information
        coords = alignments[0].coordinates
        start_target = coords[1][0]
        end_target = coords[1][-1]
        
        # Extract target sequence directly from coordinates
        # target = seq[start_target:end_target]
        
        # Get aligned sequences for detailed analysis
        alignment_str = str(alignments[0])
        alignment_lines = alignment_str.split('\n')
        if len(alignment_lines) >= 3:
            aligned_sgrna = [x for x in alignment_lines[0].split(' ') if x != '']
            aligned_genome = [x for x in alignment_lines[2].split(' ') if x != '']
        else:
            raise ValueError("Unexpected alignment format")
        
        assert int(aligned_sgrna[-1]) == len(sgRNA)
        
        # Calculate indels and mismatches
        # deletion = RNA bulge
        # insertion = DNA bulge
        aligned_sgrna_seq = aligned_sgrna[-2]
        aligned_genome_seq = aligned_genome[-2]
        insertion = aligned_sgrna_seq.count('-') if '-' in aligned_sgrna_seq else 0
        deletion = aligned_genome_seq.count('-') if '-' in aligned_genome_seq else 0
        
        # Count mismatches by comparing sequences directly
        # mismatch = 0
        # assert len(aligned_sgrna_seq) == len(aligned_genome_seq)
        # for i in range(len(aligned_sgrna_seq)):
        #     if (aligned_sgrna_seq[i] != aligned_genome_seq[i]) & (aligned_sgrna_seq[i] != 'N') & (aligned_genome_seq[i] != 'N'):
        #         mismatch += 1

        mismatch = round((alignments[0].score % 1)/mismatch_score)
        
        # Calculate target location
        pos_st = int(a_key.split('-')[0].split(':')[1]) + 1
        chr_name = a_key.split(':')[0]
        target_st = pos_st + start_target
        target_ed = pos_st + end_target - 1
        target_location = f"{chr_name}:{target_st}-{target_ed}"
        
        score = alignments[0].score
        
        return [score, aligned_genome_seq, target_location, deletion, insertion, mismatch]
            
    except Exception as e:
        print(f"Alignment error for {a_key}: {e}")
        return [0, 0, '', f"{a_key.split(':')[0]}:0-0", 0, 0, len(sgRNA)]

def local_realign(sgRNA_seq, fasta, PAM='NGG', PAM_loc='downstream'):
    # 添加 PAM
    if PAM_loc == 'downstream':
        sgRNA_PAM_fw = sgRNA_seq + PAM
    else:
        sgRNA_PAM_fw = PAM + sgRNA_seq
    sgRNA_PAM_rv = xseq.reverse_complement(sgRNA_PAM_fw)
    list_args_fw=[]
    list_args_rv=[]
    for a_key, a_seq in fasta.items():
        # 2025.04.25 修正大小写问题
        a_seq = re.sub('[^ATCG]','N',a_seq.upper())
        list_args_fw.append( [a_key, sgRNA_PAM_fw, a_seq])
        list_args_rv.append( [a_key, sgRNA_PAM_rv, a_seq])
    list_align_forward = [sgRNA_alignment_new(*args) for args in list_args_fw]
    list_align_reverse = [sgRNA_alignment_new(*args) for args in list_args_rv]
    #
    df_align_forward = pd.DataFrame(list_align_forward, columns= ['fw_score', 'fw_target','fw_location','fw_deletion','fw_insertion','fw_mismatch'])
    df_align_reverse = pd.DataFrame(list_align_reverse, columns= ['rv_score', 'rv_target','rv_location','rv_deletion','rv_insertion','rv_mismatch'])
    df_align_reverse['rv_target'] = df_align_reverse['rv_target'].apply(xseq.reverse_complement)
    df_candidate = pd.concat([df_align_forward,df_align_reverse],axis=1)
    df_candidate['location'] = fasta.keys()
    df_candidate['alignment_score'] = df_candidate[['fw_score','rv_score']].max(axis=1)
    df_candidate['best_seq_score'] = df_candidate[['fw_score', 'rv_score']].max(axis=1)
    df_candidate['best_strand'] = df_candidate[['fw_score', 'rv_score']].idxmax(axis='columns').replace({'fw_score':'+', 'rv_score':'-'})
    df_candidate.loc[df_candidate['fw_score']==df_candidate['rv_score'],'best_strand']='equal_score'
            
    # GG check
    # 2023.12.05 增加 cleavage_site 推测
    list_best_target = []
    list_best_location = []
    list_cleavage_site = []
    list_delete = []
    list_insert = []
    list_mismat = []
    list_GG = []
    for a_row in df_candidate.iterrows():
        if a_row[1]['best_strand']=='+':
            list_best_target.append(a_row[1]['fw_target'])
            list_best_location.append(a_row[1]['fw_location'])
            list_cleavage_site.append(int(a_row[1]['fw_location'].split('-')[1]) - 6)
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
            list_cleavage_site.append(int(a_row[1]['rv_location'].split('-')[0].split(':')[1]) + 5)
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
                list_cleavage_site.append(int(a_row[1]['fw_location'].split('-')[1]) - 6)
                list_delete.append(a_row[1]['fw_deletion'])
                list_insert.append(a_row[1]['fw_insertion'])
                list_mismat.append(a_row[1]['fw_mismatch'])
                list_GG.append('OK_same_score')
            # 发现没有 GG 则看 RC
            elif a_row[1]['rv_target'][-2:]=='GG':
                list_best_target.append(a_row[1]['rv_target'])
                list_best_location.append(a_row[1]['rv_location'])
                list_cleavage_site.append(int(a_row[1]['rv_location'].split('-')[0].split(':')[1]) + 5)
                list_delete.append(a_row[1]['rv_deletion'])
                list_insert.append(a_row[1]['rv_insertion'])
                list_mismat.append(a_row[1]['rv_mismatch'])
                list_GG.append('OK_same_score')
            else:
                list_best_target.append(a_row[1]['fw_target'])
                list_best_location.append(a_row[1]['fw_location'])
                list_cleavage_site.append(int(a_row[1]['fw_location'].split('-')[1]) - 6)
                list_delete.append(a_row[1]['fw_deletion'])
                list_insert.append(a_row[1]['fw_insertion'])
                list_mismat.append(a_row[1]['fw_mismatch'])                    
                list_GG.append('NO_same_score')
    # 记入 df_candidate
    df_candidate['deletion'] = list_delete
    df_candidate['insertion'] = list_insert
    df_candidate['mismatch'] = list_mismat
    df_candidate['GG'] = list_GG
    df_candidate['best_target'] = list_best_target
    df_candidate['target_location'] = list_best_location
    df_candidate['cleavage_site'] = list_cleavage_site
    df_candidate = pd.concat([xseq.bedfmt(df_candidate['target_location']), df_candidate], axis=1)

    return df_candidate

def left_realign(dp_bdg_chr, loc_shift_left, ref_fasta, sgRNA_seq, PAM, PAM_loc, n_iter):
    # print(loc_shift_left)
    fasta = xseq.get_seq(loc_shift_left, ref_fasta)
    df_candidate = local_realign(sgRNA_seq, fasta, PAM, PAM_loc)
    sr_candidate = df_candidate.iloc[0].copy()
    chrom = sr_candidate['chr']
    cleavage_site = sr_candidate['cleavage_site']
    flank_regions = [500]
    signals = target_signal(dp_bdg_chr.to_pandas(), chrom, cleavage_site, flank_regions=flank_regions)
    L_neg_1000 = signals[2]
    R_neg_1000 = signals[5]
    # 如果右侧范围变负数了，说明过头了
    if R_neg_1000 < 0:
        sr_candidate.loc['realign'] = 'fail'
        return sr_candidate

    # 计算左移后的 L_neg_1000，如果还是负数则迭代，最多迭代 10 次
    if L_neg_1000 < 0:
        st = sr_candidate['st']
        ed = sr_candidate['ed']
        loc_shift_left = f'{chrom}:{int(st)-1000}-{int(ed)-20}'
        n_iter += 1
        if n_iter < 10:
            return left_realign(dp_bdg_chr, loc_shift_left, ref_fasta, sgRNA_seq, PAM, PAM_loc, n_iter)
        else:
            sr_candidate.loc['realign'] = 'fail'
            return sr_candidate
    else:
        sr_candidate.loc['realign'] = 'success'
        return sr_candidate

def right_realign(dp_bdg_chr, loc_shift_right, ref_fasta, sgRNA_seq, PAM, PAM_loc, n_iter):
    # print(loc_shift_right)
    fasta = xseq.get_seq(loc_shift_right, ref_fasta)
    df_candidate = local_realign(sgRNA_seq, fasta, PAM, PAM_loc)
    sr_candidate = df_candidate.iloc[0].copy()
    chrom = sr_candidate['chr']
    cleavage_site = sr_candidate['cleavage_site']
    flank_regions = [500]
    signals = target_signal(dp_bdg_chr.to_pandas(), chrom, cleavage_site, flank_regions=flank_regions)
    L_neg_1000 = signals[2]
    R_neg_1000 = signals[5]
    # 如果左侧范围变负数了，说明过头了
    if L_neg_1000 < 0:
        sr_candidate.loc['realign'] = 'fail'
        return sr_candidate

    # 计算右移后的 R_neg_1000，如果还是负数则迭代，最多迭代 10 次
    if R_neg_1000 < 0:
        st = sr_candidate['st']
        ed = sr_candidate['ed']
        loc_shift_right = f'{chrom}:{int(st)+20}-{int(ed)+1000}'
        n_iter += 1
        if n_iter < 10:
            return right_realign(dp_bdg_chr, loc_shift_right, ref_fasta, sgRNA_seq, PAM, PAM_loc, n_iter)
        else:
            sr_candidate.loc['realign'] = 'fail'
            return sr_candidate
    else:
        sr_candidate.loc['realign'] = 'success'
        return sr_candidate


