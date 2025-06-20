# Offtracker

Offtracker is an end to end pipeline of Tracking-seq data analysis for detecting off-target sites of any genome editing tools that generate double-strand breaks (DSBs) or single-strand breaks (SSBs).

## System requirements

* Linux/Unix 
* Python >= 3.6

## Dependency

```bash
# We recommend creating a new environment using mamba/conda to avoid compatibility problems
# If you don't use mamba, just replace the code with conda 
# Windows systems may not be compatible with pybedtools.
mamba create -n offtracker -c bioconda blast snakemake pybedtools deeptools chromap
```


## Installation 


```bash
# Activate the environment
conda activate offtracker

# Direct installation with pip
pip install offtracker

# (Alternative) Download the offtracker from github
git clone https://github.com/Lan-lab/offtracker.git 
cd offtracker
pip install .
```


## Before analyzing samples

**Important: Do not use hard-masked genome.fa**, in which repeats are masked by capital Ns and reads should have been mapped to these region (e.g. MHC region) will be lost. Besides, the genome.fa **should not** contain alternate loci like chr2_KI270776v1_alt and chr6_GL000256v2_alt, which may cause multi-mappings and the reads may be discarded.

For example, https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz is soft-masked genome with alternate loci. https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.masked.gz is hard-masked genome. **Do not** use these two as reference genome.

http://cistrome.org/~galib/MAESTRO/references/scATAC/Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz is the genome used for the example data.

```bash
# The following command can be used to check whether alternate loci of chr6 are present in the reference genome.
grep "^>chr6" genome.fa
```

```bash
# Build chromap index (only need once for each genome)
chromap -i -r /Your_Path_To_Reference/hg38_genome.fa \
-o /Your_Path_To_Reference/hg38_genome.chromap.index

# Build blast index (only need once for each genome)
makeblastdb -input_type fasta -title hg38 -dbtype nucl -parse_seqids \
-in /Your_Path_To_Reference/hg38_genome.fa \
-out /Your_Path_To_Reference/hg38_genome.blastdb \
-logfile /Your_Path_To_Reference/hg38_genome.blastdb.log

# Generate candidate regions by sgRNA sequence (need once for each genome and sgRNA)
# --name: a user-defined name of the sgRNA, which will be used in the following analysis.
offtracker_candidates.py -t 8 -g hg38 \
-r /Your_Path_To_Reference/hg38_genome.fa \
-b /Your_Path_To_Reference/hg38_genome.blastdb \
--name 'VEGFA2' --sgrna 'GACCCCCTCCACCCCGCCTC' --pam 'NGG' \
-o /Your_Path_To_Candidates_Folder

```


## Quality control and adapter trimming

```bash
# Generate snakemake config file for quality control and adapter trimming.
offtracker_qc.py -t 4 \
-f /Your_Path_To_Input_Folder \
--subfolder 0

cd /Your_Path_To_Input_Folder/Trimmed_data
snakemake -np # dry run to check whether everything is alright
nohup snakemake --cores 16 1>${outdir}/sm_qc.log 2>&1 &

"""
Set “--subfolder 0” if the file structure is like: 
| - Input_Folder
  | - sample1_R1.fastq.gz
  | - sample1_R2.fastq.gz
  | - sample2_R1.fastq.gz
  | - sample2_R2.fastq.gz
Set “--subfolder 1” if the file structure is like:
| - Input_Folder
  | - Sample1_Folder
    | - sample1_R1.fastq.gz
    | - sample1_R2.fastq.gz
  | - Sample2_Folder
    | - sample2_R1.fastq.gz
    | - sample2_R2.fastq.gz

The script “offtracker_qc.py” will create a “Trimmed_data” folder under /Your_Path_To_Input_Folder. 
If “-o /Your_Path_To_Output” is set, the output will be redirected to /Your_Path_To_Output.
"""
```

## Strand-specific mapping of Tracking-seq data 

```bash

# Generate snakemake config file for mapping
# Results will be generated in /Your_Path_To_Output, if -o is not set, the output will be in the same folder as the fastq files
offtracker_config.py -t 8 -g hg38 --blacklist hg38 \
-r /Your_Path_To_Reference/hg38_genome.fa \
-i /Your_Path_To_Reference/hg38_genome.chromap.index \
-f /Your_Path_To_Trimmed_Data \
-o /Your_Path_To_Output \ 
--subfolder 0 

# Warning: Do not contain "fastq" or "fq" in the folder name, otherwise the program may treat the folder as a fastq file
# This problem may be fixed in the future

# Run the snakemake program
cd /Your_Path_To_Fastq
snakemake -np # dry run
nohup snakemake --cores 16 1>sm_mapping.log 2>sm_mapping.err &

## about cores
# --cores of snakemake must be larger than -t of offtracker_config.py
# parallel number = cores/t

## about output
# This part will generate "*.fw.scaled.bw" and ".rv.scaled.bw" for IGV visualization
# "*.fw.bed" and "*.rv.bed" are used in the next part.
```


## Analyzing the genome-wide off-target sites

```bash
# In this part, multiple samples in the same condition can be analyzed in a single run by pattern recognition of sample names

offtracker_analysis.py -g hg38 --name "VEGFA2" \
--exp 'Cas9_VEGFA2' \
--control 'WT' \
--outname 'Cas9_VEGFA_293' \
-f /Your_Path_To_Output \
--seqfolder /Your_Path_To_Candidates

# --name: the same gRNA name you set when running offtracker_candidates.py
# --exp/--control: add one or multiple patterns of file name in regular expressions
# If multiple samples meet the pattern, their signals will be averaged. Thus, only samples with the same condition should be included in a single analysis.

# This step will generate Offtracker_result_{outname}.csv
# Default FDR is 0.05, which can be changed by --fdr. This will empirically make the threshold of Track score around 2.
# Sites with Track score >=2, which is a empirical threshold, are output regardless of FDR.
# Intermediate files are saved in ./temp folder, which can be deleted.
# Keeping the intermediate files can make the analysis faster if involving previously analyzed samples (e.g. using the same control samples for different analyses)
```

## Off-target sequences visualization

```bash
# After get the Offtracker_result_{outname}.csv, you can visualize the off-target sites with their genomic sequence with the following command:

offtracker_plot.py --result Your_Offtracker_Result_CSV \
--sgrna 'GACCCCCTCCACCCCGCCTC' --pam 'NGG'

# The default output is a pdf file with Offtracker_result_{outname}.pdf
# Assigning a specific output file with another suffix can change the format. e.g., "--output Offtracker_plot.png" will generate a png file.
# The orange dash line indicates the empirical threshold of Track score = 2
# Empirically, the off-target sites with Track score < 2 are less likely to be real off-target sites.
```


## Note1, when not using hg38 or mm10

The default setting only includes chr1-chr22, chrX, chrY, and chrM. (only suitable for human and mouse) \
If you are using reference genomes without "chr" at the beginning, or want to analyze all chromosomes or other species, you can set "--ignore_chr" when running offtracker_config.py to skip chromosome filter.

Currently, this software is only ready-to-use for mm10 and hg38. For any other genome, e.g., hg19, please add a genome size file named "hg19.chrom.sizes" to .\offtracker\utility. Besides, add "--blacklist none" or "--blacklist Your_Blacklist" (e.g., ENCODE blacklist) when running offtracker_config.py, because we only include blacklists for mm10 and hg38.

## Note2

The FDRs in the Tracking-seq result do not reflect the real off-target probability.
It is strongly recommended to observe the "fw.scaled.bw" and "rv.scaled.bw" using genome browser like IGV to visually inspect each target location from the Tracking-seq result.



# Example Data

Here are example data that contains reads of chr6 from HEK293T cells edited with Cas9 + sgRNA VEGFA_site_2 (VEGFA2) and reads of chr6 from wild type HEK293T cells:

https://figshare.com/articles/dataset/WT_HEK239T_chr6/25956034

It takes about 5-10 minutes to run the mapping (offtracker_config.py & snakemake) of example data with -t 8 and --cores 16 (2 parallel tasks)

## Signal visualization

After mapping, there will be 4 .bw files in the output folder:
```bash
Cas9_VEGFA2_chr6.fw.scaled.bw

Cas9_VEGFA2_chr6.rv.scaled.bw

WT_chr6.fw.scaled.bw

WT_chr6.rv.scaled.bw
```
These files can be visualized in genome browser like IGV:

![signal](https://github.com/Lan-lab/offtracker/blob/main/example_output/signals_example.png?raw=true)

The signal (coverage) for each sample is normalized to 1e7/total_reads. As only reads mapping to chr6 were extracted in the example data, the signal range is much higher than that of the whole genome samples. 

## Whole genome off-target analysis

For analyzing the signals (offtracker_analysis.py), it takes about 3-5 minutes and outputs a file named "Offtracker_result_{outname}.csv"

After that, you can visualize the off-target sites with their genomic sequence (offtracker_plot.py) and get an image like this:

![offtarget](https://github.com/Lan-lab/offtracker/blob/main/example_output/sequences_example.png?raw=true)


After finishing the pipeline, if “chr6:31400832-31400854” and “chr6:31495044-31495066” are missing in the plot, it is most likely due to either:

•	Using a hard-masked reference genome (where repeats are replaced with 'N's)

•	The presence of alternate loci (e.g., chr6_GL000256v2_alt) in the genome.

These two off-target sites locate in the region of MHC class I chain-related protein A and B (MICA and MICB), which is polymorphic (resulting in alternate loci in contigs like “chr6_GL000256v2_alt”) and contains interspersed repeats (resulting in sequences masked by capital 'N's in a hard-masked genome). Please try again with unmasked or soft-masked genome without alternate loci.



# Citation

If you use Tracking-seq or OFF-TRACKER in your research, please cite the following paper:

Zhu, M., Xu, R., Yuan, J., Wang, J. et al. Tracking-seq reveals the heterogeneity of off-target effects in CRISPR–Cas9-mediated genome editing. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02307-y

The signal visualization of .bw file here was generated by the Integrative Genomics Viewer (IGV) software. The signal visualization in the Tracking-seq article above was generated by either IGV or pyGenomeTracks:

Robinson, J., Thorvaldsdóttir, H., Winckler, W. et al. Integrative genomics viewer. Nat Biotechnol 29, 24–26 (2011). https://doi.org/10.1038/nbt.1754

Lopez-Delisle L, Rabbani L, Wolff J, Bhardwaj V, Backofen R, Grüning B, Ramírez F, Manke T. pyGenomeTracks: reproducible plots for multivariate genomic data sets. Bioinformatics. 2020 Aug 3:btaa692. doi: 10.1093/bioinformatics/btaa692.

