OFF-TRACKER
=======================

OFF-TRACKER is an end to end pipeline of Track-seq data analysis for detecting off-target sites of any genome editing tools that generate double-strand breaks (DSBs) or single-strand breaks (SSBs).

System requirements
-----
* Linux/Unix 
* Python >= 3.6

Dependency
-----

```bash
# We recommend creating a new enviroment using mamba/conda to avoid compatibility problems
# If you don't use mamba, just replace the code with conda 
mamba create -n offtracker -c bioconda blast snakemake pybedtools
```


Installation 
-----

```bash
# activate the environment
conda activate offtracker

# Direct installation with pip
pip install offtracker

# (Alternative) Download the offtracker from github
git clone https://github.com/Lan-lab/offtracker.git 
cd offtracker
pip install .
```


Before analyzing samples
-----

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
offtracker_candidates.py -t 8 -g hg38 \
-r /Your_Path_To_Reference/hg38_genome.fa \
-b /Your_Path_To_Reference/hg38_genome.blastdb \
--name 'HEK4' --sgrna 'GGCACTGCGGCTGGAGGTGG' --pam 'NGG' \
-o /Your_Path_To_Candidates

```

Strand-specific mapping of Track-seq data 
-----

```bash
# Generate snakemake config file 
offtracker_config.py -t 8 -g hg38 --blacklist hg38 \
-r /Your_Path_To_Reference/hg38_genome.fa \
-i /Your_Path_To_Reference/hg38_genome.chromap.index \
-f /Your_Path_To_Fastq \
-o /Your_Path_To_Output \ 
--subfolder 0 

# --subfolder: If different samples are in seperate folders, set this to 1
# -o: Default is outputting to /Your_Path_To_Fastq

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


Analyzing the off-target sites
-----

```bash
# In this part, multiple samples in the same condition can be analyzed in a single run by pattern recogonization of sample names

offtracker_analysis.py -g hg38 --name "HEK4" \
--exp 'Cas9_HEK4.*293' \
--control 'control' \
--outname 'Cas9_HEK4_293' \
-f /Your_Path_To_Output \
--seqfolder /Your_Path_To_Candidates

# --name: the same as that in offtracker_candidates.py
# --exp/--control: add one or multiple patterns of file name in regex


# This step will generate Trackseq_result_{outname}.csv
# Intermediate files are saved in ./temp folder, which can be deleted 
# Keeping the intermediate files can make the analysis faster if involving previously analyzed samples (e.g. using the same control samples for different analyses)
```


Note1
--------------
The default setting only includes chr1-chr22, chrX, chrY, and chrM.

Please make sure the reference genome contains "chr" at the beginning. 

If you have requirement for other chromosomes or species other than human/mouse, please post an issue.

Note2
--------------
Currently, this software is only ready-to-use for mm10 and hg38. 

For any other genome, say hg19, please add genome size file named "hg19.chrom.sizes" to .\offtracker\mapping before install.

Besides, add "--blacklist none" or "--blacklist Your_Blacklist" when running offtracker_config.py

Note3
--------------
Instead of setting a fixed threhold for track score or FDR, it is strongly recommended to observe the "fw.scaled.bw" and "rv.scaled.bw" using IGV to check top target locations from the Track-seq result.

