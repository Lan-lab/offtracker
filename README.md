# OFF-TRACKER

OFF-TRACKER is an end to end pipeline of Tracking-seq data analysis for detecting off-target sites of any genome editing tools that generate double-strand breaks (DSBs) or single-strand breaks (SSBs).

## System requirements

* Linux/Unix 
* Python >= 3.6

## Dependency

```bash
# We recommend creating a new enviroment using mamba/conda to avoid compatibility problems
# If you don't use mamba, just replace the code with conda 
mamba create -n offtracker -c bioconda blast snakemake pybedtools
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

```bash
# Build blast index (only need once for each genome)
makeblastdb -input_type fasta -title hg38 -dbtype nucl -parse_seqids \
-in /Your_Path_To_Reference/hg38_genome.fa \
-out /Your_Path_To_Reference/hg38_genome.blastdb \
-logfile /Your_Path_To_Reference/hg38_genome.blastdb.log

# Build chromap index (only need once for each genome)
chromap -i -r /Your_Path_To_Reference/hg38_genome.fa \
-o /Your_Path_To_Reference/hg38_genome.chromap.index

# Generate candidate regions by sgRNA sequence (need once for each genome and sgRNA)
# --name: the name of the sgRNA, which will be used in the following analysis
offtracker_candidates.py -t 8 -g hg38 \
-r /Your_Path_To_Reference/hg38_genome.fa \
-b /Your_Path_To_Reference/hg38_genome.blastdb \
--name 'VEGFA2' --sgrna 'GACCCCCTCCACCCCGCCTC' --pam 'NGG' \
-o /Your_Path_To_Candidates

```

## Strand-specific mapping of Tracking-seq data 

```bash
# Generate snakemake config file 
# --subfolder: If different samples are in seperate folders, set this to 1
# Results will be generated in /Your_Path_To_Output, if -o is not set, the output will be in the same folder as the fastq files
offtracker_config.py -t 8 -g hg38 --blacklist hg38 \
-r /Your_Path_To_Reference/hg38_genome.fa \
-i /Your_Path_To_Reference/hg38_genome.chromap.index \
-f /Your_Path_To_Fastq \
-o /Your_Path_To_Output \ 
--subfolder 0 

# Warning: Do not contain "fastq" or "fq" in the folder name, otherwise the program will treat the folder as a fastq file
# This problem will be fixed in the future version

# Run the snakemake program
cd /Your_Path_To_Fastq
snakemake -np # dry run
nohup snakemake --cores 16 1>snakemake.log 2>snakemake.err &

## about cores
# --cores of snakemake must be larger than -t of offtracker_config.py
# parallel number = cores/t

## about output
# This part will generate "*.fw.scaled.bw" and ".rv.scaled.bw" for IGV visualization
# "*.fw.bed" and "*.rv.bed" are used in the next part.
```


## Analyzing the genome-wide off-target sites

```bash
# In this part, multiple samples in the same condition can be analyzed in a single run by pattern recogonization of sample names

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
# Change the suffix of the output file to change the format (e.g.: .png)
# The orange dash line indicates the empirical threshold of Track score = 2
# Empirically, the off-target sites with Track score < 2 are less likely to be real off-target sites.
```


## Note1

The default setting only includes chr1-chr22, chrX, chrY, and chrM. Please make sure the reference genome contains "chr" at the beginning. 

Currently, this software is only ready-to-use for mm10 and hg38. For any other genome, e.g., hg19, please add genome size file named "hg19.chrom.sizes" to .\offtracker\mapping and instal manually. Besides, add "--blacklist none" or "--blacklist Your_Blacklist" (e.g., ENCODE blacklist) when running offtracker_config.py, because we only provide blacklists for mm10 and hg38.

If you have a requirement for species other than human/mouse, please post an issue.

## Note2

The FDRs in the Tracking-seq result do not reflect the real off-target probability.
It is strongly recommended to observe the "fw.scaled.bw" and "rv.scaled.bw" using genome browser like IGV to visually inspect each target location from the Tracking-seq result.



# Example Data

Here are example data that contains reads of chr6 from HEK293T cells edited with Cas9 + sgRNA VEGFA2 and wild type cells:

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

The signal (coverage) for each sample is normalized to 1e7/total_reads. As only reads mapping to chr6 were extracted in the example data, the signal range is higher than that of the whole genome. 

## Whole genome off-target analysis

For analyzing the signals (offtracker_analysis.py), it takes about 3-5 minutes and outputs a file named "Offtracker_result_{outname}.csv"

After that, you can visualize the off-target sites with their genomic sequence (offtracker_plot.py) and get an image like this:

![offtarget](https://github.com/Lan-lab/offtracker/blob/main/example_output/sequences_example.png?raw=true)

# Citation

If you use Tracking-seq or OFF-TRACKER in your research, please cite the following paper:

Zhu, M., Xu, R., Yuan, J., Wang, J. et al. Tracking-seq reveals the heterogeneity of off-target effects in CRISPR–Cas9-mediated genome editing. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02307-y

The signal visualization of .bw file here was generated by the Integrative Genomics Viewer (IGV) software. The signal visualization in the Tracking-seq article above was generated by either IGV or pyGenomeTracks:

Robinson, J., Thorvaldsdóttir, H., Winckler, W. et al. Integrative genomics viewer. Nat Biotechnol 29, 24–26 (2011). https://doi.org/10.1038/nbt.1754

Lopez-Delisle L, Rabbani L, Wolff J, Bhardwaj V, Backofen R, Grüning B, Ramírez F, Manke T. pyGenomeTracks: reproducible plots for multivariate genomic data sets. Bioinformatics. 2020 Aug 3:btaa692. doi: 10.1093/bioinformatics/btaa692.

