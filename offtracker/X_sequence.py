
import math
import pandas as pd
from itertools import product
import numpy as np

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
                'N': ['A', 'T', 'C', 'G']}

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

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-':'-',
                  'M': 'K', 'R': 'Y', 'W': 'W', 'S': 'S', 'Y': 'R', 'K':'M',
                  'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V'} 
    bases = list(seq) 
    letters = [complement[base] for base in bases] 
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

def sgRNA_alignment(a_key, sgRNA, seq, frag_len, DNA_matrix=None, mismatch_score = 0.01, return_align=False):
    from Bio import pairwise2
    import numpy as np
    if DNA_matrix is None:
        DNA_matrix = {('A','A'): 2, ('A','T'):0.01, ('A','C'):0.01, ('A','G'):0.01, ('A','N'):0.01,
                    ('T','T'): 2, ('T','A'):0.01, ('T','C'):0.01, ('T','G'):0.01, ('T','N'):0.01,
                    ('G','G'): 2, ('G','A'):0.01, ('G','C'):0.01, ('G','T'):0.01, ('G','N'):0.01,
                    ('C','C'): 2, ('C','A'):0.01, ('C','G'):0.01, ('C','T'):0.01, ('C','N'):0.01,
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
