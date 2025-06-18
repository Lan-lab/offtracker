
import math
import pandas as pd
from itertools import product
import numpy as np
import os, glob

ambiguous_nt = {'A': ['A'],
                'T': ['T'],
                'C': ['C'],
                'G': ['G'],
                'U': ['T'],
                'R': ['A', 'G'],
                'W': ['A', 'T'],
                'M': ['A', 'C'],
                'S': ['C', 'G'],
                'Y': ['C', 'T'],
                'K': ['G', 'T'],
                'V': ['A', 'C', 'G'],
                'H': ['A', 'C', 'T'],
                'D': ['A', 'G', 'T'],
                'B': ['C', 'G', 'T'],
                'N': ['A', 'C', 'G', 'T']}

def is_seq_valid(sequence, extra=True, ambiguous_nt=ambiguous_nt):
    if extra:
        valid_nucleotides = list(ambiguous_nt.keys())
    else:
        valid_nucleotides = ['A', 'T', 'C', 'G']
    for nucleotide in sequence:
        if nucleotide not in valid_nucleotides:
            return nucleotide
    return True

def possible_seq(sequence):
    valid_check = is_seq_valid(sequence)
    if valid_check == True:
        possible_nucleotides = []
        for x in sequence:
            possible_nucleotides.append(ambiguous_nt[x])
        possible_combinations = list(product(*possible_nucleotides))
        sequences = [''.join(combination) for combination in possible_combinations]
    else:
        raise KeyError(f'Unvalid character \'{valid_check}\' in sequence')
    return sequences

# 包含 degenerate base pairs
def get_base_score(base1, base2, exact_score=2, partial_match=2, mismatch_score=0.01):
    base1 = ambiguous_nt[base1]
    base2 = ambiguous_nt[base2]
    if base1 == base2:
        return exact_score
    if list(np.union1d(base1,base2)) == base1 or list(np.union1d(base1,base2)) == base2:
        # 其中一个是子集，注意顺序不一致会导致不等，所以必须排好序
        return partial_match
    return mismatch_score


def complement(seq):
    dict_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-':'-',
                  'M': 'K', 'R': 'Y', 'W': 'W', 'S': 'S', 'Y': 'R', 'K':'M',
                  'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V'} 
    bases = list(seq) 
    letters = [dict_complement[base] for base in bases] 
    return ''.join(letters)

def reverse(seq):
    return seq[::-1]

def reverse_complement(seq):
    return complement(seq[::-1])


def write_fasta(df, output, name_col = 'ID', sequence_col='sequence', line_len = 60): 
    with open(output,'w') as f: 
        for i in range(len(df)):
            f.write( '>' + df[name_col].iloc[i] + '\n')
            n_line = math.ceil( len(df[sequence_col].iloc[i])/line_len)
            for j in range(n_line):
                f.write( df[sequence_col].iloc[i][ j*line_len : (j+1)*line_len ]   + '\n')
    return 'fasta is written.'

def write_bed(df, bed_dir):
    return df.to_csv(bed_dir, sep='\t', header=None, index=False)

def read_bed(bed_dir):
    return pd.read_csv(bed_dir,sep='\t',header=None)

def X_readbed(bed_dir):
    bed = pd.read_csv(bed_dir,sep='\t',header=None)
    bed.columns = ['chr','st','ed'] + list(bed.columns[3:])
    return bed

def igvfmt(bed):
    sr_igv = bed.iloc[:,0].str[:] + ':' + bed.iloc[:,1].astype(str).str[:] + '-' + bed.iloc[:,2].astype(str).str[:]
    return sr_igv

def bedfmt(igv):
    igv = pd.Series(igv)
    igv = igv.str.extract('(.*):(.*)-(.*)')
    igv.columns = ['chr','st','ed']
    igv['st'] = igv['st'].astype(int)
    igv['ed'] = igv['ed'].astype(int)
    return igv

def add_ID(df, chr_col=0, midpoint='cleavage_site'):#, midpoint='midpoint'):
	chr_col_name = df.columns[chr_col]
	print(f'chromosome col = {chr_col_name}')
	point_head = (df[midpoint]/1000).astype(int)
	df['ID_1'] = df[chr_col_name] + ':' + point_head.astype(str)
	point_tail = df[midpoint] % 1000
	df.loc[point_tail<500,'ID_2'] = df[chr_col_name] + ':' + (point_head-1).astype(str)
	df.loc[point_tail>=500,'ID_2'] = df[chr_col_name] + ':' + (point_head+1).astype(str)
	return df



def detect_fastq(folder, n_subfolder, NGS_type='paired-end'):
    """
    搜索 folder 的 n级子目录下的所有 fastq/fastq.gz/fq/fq.gz 文件
    paired-end 模式 : 识别 2.fq/2.fastq 为 paired-end 的 R2 文件，并验证对应 R1 文件
    single-end 模式 : 所有 fastq/fastq.gz/fq/fq.gz 文件都视为 single-end 文件
    
    不建议 2. 和 fq/fastq 之间有其他字符，如 2.trimmed.fq.gz，因为中间字符不确定，使用通配符容易误判文件名其他的2.
    样本名不要带点，建议用_分割特征，同特征内分割不要用_可以用-，如 sample_day-hour_type_batch_rep_1.fq.gz

    Input
    ----------
    folder : 根目录
    n_subfolder : n级子目录

    Parameter
    ----------
    NGS_type : 'paired-end' or 'single-end'
    
    Output
    ----------
    sample_names : 识别的样品名
    files_R1 : R1文件的完整路径
    files_R2 : R2文件的完整路径

    """
    # import os, sys, glob
    # import pandas as pd
    if NGS_type == 'paired-end':
        print('paired-end mode')
        files_R2 = []
        # 支持四种文件扩展名
        # 个人习惯包含绝对路径
        for fastq in ['*2.fq','*2.fastq','*2.fq.gz','*2.fastq.gz']:
            fq_files = glob.glob( os.path.join(folder, n_subfolder*'*/', fastq ) )
            print(f'{len(fq_files)} {fastq[2:]} samples detected')
            files_R2.extend( fq_files )
        #
        if len(files_R2) > 0:
            files_R2 = pd.Series(files_R2).sort_values().reset_index(drop=True)
            # 拆分文件名
            suffix = files_R2.str.extract(r'(\.fastq.*|\.fq.*)',expand=False)
            prefix = files_R2.str.extract(r'(.*)(?:.fq|.fastq)',expand=False)
            # 将 prefix 进一步拆分为 sample_dir （真样品名） 和 nametype （某种统一后缀），支持五种样本名后缀
            nametype = []
            sample_dir = []
            for a_prefix in prefix:
                for a_type in ['_trimmed_2', '_2_val_2','_R2_val_2','_R2','_2']:
                    len_type = len(a_type)
                    if a_prefix[-len_type:] == a_type:
                        nametype.append(a_type)
                        sample_dir.append(a_prefix[:-len_type])
                        break
            assert len(nametype) == len(files_R2), 'The file name pattern is invaild!'
            nametype = pd.Series(nametype)
            sample_dir = pd.Series(sample_dir)
            # 根据 R2 文件，检查 R1 文件是否存在
            files_R1 = sample_dir + nametype.str.replace('2','1') + suffix
            for i in range(len(files_R1)):
                assert os.path.exists(files_R1[i]), f'{files_R1[i]} not found!'
            sample_names = sample_dir.apply(os.path.basename)
        else:
            print('No paired-end samples detected!')
            sample_names = 'no sample'
            files_R1 = []

    elif NGS_type == 'single-end':
        print('single-end mode')
        files_R1 = []
        files_R2 = [] # 占位
        # 支持四种文件扩展名
        # 个人习惯包含绝对路径
        for fastq in ['*.fq','*.fastq','*.fq.gz','*.fastq.gz']:
            fq_files = glob.glob( os.path.join(folder, n_subfolder*'*/', fastq ) )
            print(f'{len(fq_files)} {fastq[1:]} samples detected')
            files_R1.extend( fq_files )
        files_R1 = pd.Series(files_R1).sort_values()
        #
        if len(files_R1) > 0:
            # 拆分文件名
            suffix = files_R1.str.extract(r'(\.fastq.*|\.fq.*)',expand=False)
            prefix = files_R1.str.extract(r'(.*)(?:.fq|.fastq)',expand=False)
            # 单端模式下，所有前缀都视为样品名
            sample_names = prefix.apply(os.path.basename)
        else:
            print('No single-end samples detected!')
            sample_names = 'no sample'
            files_R1 = []

    return sample_names, files_R1, files_R2


def sgRNA_alignment(a_key, sgRNA, seq, frag_len, DNA_matrix=None, mismatch_score = 0.01, return_align=False):
    from Bio import pairwise2
    import numpy as np
    if DNA_matrix is None:
        DNA_matrix = {('A','A'): 2, ('A','T'):0.01, ('A','C'):0.01, ('A','G'):0.01, ('A','N'):2,
                    ('T','T'): 2, ('T','A'):0.01, ('T','C'):0.01, ('T','G'):0.01, ('T','N'):2,
                    ('G','G'): 2, ('G','A'):0.01, ('G','C'):0.01, ('G','T'):0.01, ('G','N'):2,
                    ('C','C'): 2, ('C','A'):0.01, ('C','G'):0.01, ('C','T'):0.01, ('C','N'):2,
                    ('N','N'): 2, ('N','C'):2, ('N','A'): 2, ('N','G'): 2, ('N','T'): 2}        
    # a_key 是 pybedtools 得到的位置 chrA:X-Y 而 X 数字会往左多1bp
    alignments = pairwise2.align.localds( sgRNA, seq, DNA_matrix, -2, -2, penalize_extend_when_opening=False)
    # 有时会存在得分相同的不同 alignment 方式，选取最接近中点的
    position_pct = [alignments[x].start/(frag_len-len(sgRNA)-4) - 0.5 for x in range(len(alignments))]
    align_score  = [alignments[x].score for x in range(len(alignments))]
    mid_aligment = np.argmin(np.abs(position_pct))
    position_pct = position_pct[mid_aligment]
    best_alignment = alignments[mid_aligment]
    target = best_alignment.seqB[best_alignment.start:best_alignment.end]
    deletion = target.count('-')
    insertion = best_alignment.end - best_alignment.start - len(sgRNA)
    if insertion<0:
        #当比对到边缘有悬空时发生
        insertion=0
    # 用小数点记录，默认 DNA_matrix 下，如果 mismatch 大于 100 会出错
    mismatch = round((best_alignment.score % 1)/mismatch_score)
    # 推算 target_location
    pos_st = int(a_key.split('-')[0].split(':')[1]) + 1 # 减去多的1bp
    pos_ed = int(a_key.split('-')[1])
    chr_name = a_key.split(':')[0]
    target_st = pos_st + best_alignment.start
    target_ed = pos_st + best_alignment.end - 1 - deletion # 2023.12.05 修正 deletion 错位
    target_location = chr_name + ':' + str(target_st) + '-' + str(target_ed)
    if return_align:
        return [best_alignment.score, position_pct, target, target_location, deletion, insertion, mismatch, best_alignment.seqB]
    else:
        return [best_alignment.score, position_pct, target, target_location, deletion, insertion, mismatch]


def combine_df(list_df, op = 'mean'):
    # df 行列、结构必须一模一样，非数字部分也一模一样，只有数字不同
    df_nondigit = list_df[0].select_dtypes(exclude=[float, int])
    if op=='mean':
        df_combined = pd.concat(list_df).groupby(level=0).mean(numeric_only=True)
    elif op=='max':
        df_combined = pd.concat(list_df).groupby(level=0).max(numeric_only=True)
    elif op=='min':
        df_combined = pd.concat(list_df).groupby(level=0).min(numeric_only=True)
    else:
        print('op must be mean, max or min')
    #
    df_combined = pd.concat([df_nondigit, df_combined], axis=1)
    return df_combined
