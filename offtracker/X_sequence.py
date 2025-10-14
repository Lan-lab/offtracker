
import math
import pandas as pd
from itertools import product
import numpy as np
import os, glob, sys
import polars as pl

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



def detect_fastq(folder, n_subfolder, NGS_type='paired-end', skip_trimmed=False):
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

        if skip_trimmed:
            files_R2 = [f for f in files_R2 if '_trimmed_2.fq.gz' not in f]
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

def sgRNA_alignment_new(a_key, sgRNA, seq, substitution_matrix=None, alphabet=None, 
                   mismatch_score=0.01):
    """
    Perform local alignment using Bio.Align instead of deprecated pairwise2.
    """
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


def get_seq(location, ref_fasta, return_df=False) -> dict:
    """
    根据 genome location 取序列
    location 如果是 list, str 则形式为 "chr1:123456-123458" 或 ["chr1:123456-123458", "chr2:123456-123458"]
    location 如果是 pd.DataFrame, pl.DataFrame 则默认前三列为 bed 格式
    ref_fasta 是参考基因组 fasta 文件的路径
    默认返回字典, key 是位置, value 是序列
    如果 return_df 为 True, 则返回 pl.DataFrame, 第一列为位置, 第二列为序列

    pybedtools 返回的序列实际上没有包括坐标 start 的那个碱基，这一点和 twoBitToFa 一样
    但是 IGV/UCSC 等序列是包括 start 的碱基，blast 的结果也是包括 start 的
    所以在后续分析时要注意这一点
    """
    if sys.platform[:3]=='win':
        # windows 似乎装不了 pybedtools
        raise ValueError('windows 似乎装不了 pybedtools')
    else:
        import pybedtools
    #########
    # 根据 genome location 取序列
    #########
    if isinstance(location,(list,str)):
        bed_loc = bedfmt(location)
    elif isinstance(location,pd.DataFrame):
        bed_loc = location.iloc[:,:3]
    elif isinstance(location,pl.DataFrame):
        bed_loc = location[:,:3]
    else:
        raise ValueError('location must be a list, str or pd.DataFrame')
    
    fasta = pybedtools.example_filename(ref_fasta)
    temp_bed = './temp_amp_loc.bed'
    write_bed(bed_loc, temp_bed)
    a = pybedtools.BedTool(temp_bed)
    a = a.sequence(fi=fasta)
    with open(a.seqfn, encoding='utf-8') as f:
        dict_seq = {} # 定义一个空的字典
        for line in f:
            line = line.strip() # 去除末尾换行符
            if line[0] == '>':
                header = line[1:]
            else:
                sequence = line
                dict_seq[header] = dict_seq.get(header,'') + sequence
    
    # remove temp_amp_loc.bed
    os.remove(temp_bed)
    if return_df:
        return pl.DataFrame(list(dict_seq.items()), orient='row', schema={'location':pl.String,'sequence':pl.String})
    else:
        return dict_seq

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
