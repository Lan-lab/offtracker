# 2023.08.11. adding a option for not normalizing the bw file
# 2024.01.23. add --fixedStep to bigwigCompare for not merging neighbouring bins with equal values.
# 2025.05.22. refine the structure

configfile: "config.yaml"

# # fastq 信息
_files_R1 = config['files_R1'] # dict型, key 为 sample
_files_R2 = config['files_R2'] # dict型, key 为 sample
# # 运行参数
_output_dir = config["output_dir"]
_thread = config['thread']
_BinSize = str(config["binsize"])
_normalize = config["normalize"]


import os

if _normalize == "True":
    rule all:
        input:
            expand( os.path.join(_output_dir,"{sample}.fw.bed"), sample=_files_R1 ),
            expand( os.path.join(_output_dir,"{sample}.rv.bed"), sample=_files_R1 ),
            expand( os.path.join(_output_dir,"{sample}.fw.scaled.bw"), sample=_files_R1 ),
            expand( os.path.join(_output_dir,"{sample}.rv.scaled.bw"), sample=_files_R1 ),
            expand( os.path.join(_output_dir,"{sample}." + _BinSize + ".add.bdg"),sample=_files_R1 ),
elif _normalize == "False":
    rule all:
        input:
            expand( os.path.join(_output_dir,"{sample}.fw.bed"), sample=_files_R1 ),
            expand( os.path.join(_output_dir,"{sample}.rv.bed"), sample=_files_R1 ),
            expand( os.path.join(_output_dir,"{sample}.fw.raw.bw"), sample=_files_R1 ),
            expand( os.path.join(_output_dir,"{sample}.rv.raw.bw"), sample=_files_R1 ),
else:
    raise ValueError('Please provide "True" or "False" for "--normalize" when running offtracker_config.py')


rule chromap:
    input:
        R1=lambda w: _files_R1[w.sample],
        R2=lambda w: _files_R2[w.sample]
    threads:
        _thread
    params:
        index=config["index"],
        fasta=config["fasta"]
    output:
        temp(os.path.join(_output_dir,"{sample}.chromapx.bed"))
    shell:
        """
        chromap -l 3000 --low-mem --BED --remove-pcr-duplicates \
        --min-read-length 10 --allocate-multi-mappings \
        -x {params.index} -r {params.fasta} -t {threads} -1 {input.R1} -2 {input.R2} -o {output}
        """

if config["blacklist"] != 'none':
    rule remove_blacklist:
        input:
            os.path.join(_output_dir,"{sample}.chromapx.bed")
        threads:
            _thread
        params:
            blacklist=config["blacklist"]
        output:
            temp(os.path.join(_output_dir,"{sample}.filtered.bed"))
        shell:
            "bedtools intersect -a {input} -b {params.blacklist} -v > {output}"

    rule bed2fr:
        input:
            os.path.join(_output_dir,"{sample}.filtered.bed")
        threads:
            _thread
        params:
            dir_script=config["utility_dir"],
            ignore_chr=config["ignore_chr"],
        output:
            fw=os.path.join(_output_dir,"{sample}.fw.bed"),
            rv=os.path.join(_output_dir,"{sample}.rv.bed")
        shell:
            "python {params.dir_script}/1.1_bed2fr.py -b {input} {params.ignore_chr}"
else:
    rule bed2fr:
        input:
            os.path.join(_output_dir,"{sample}.chromapx.bed")
        threads:
            _thread
        params:
            dir_script=config["utility_dir"],
            ignore_chr=config["ignore_chr"],
        output:
            fw=os.path.join(_output_dir,"{sample}.fw.bed"),
            rv=os.path.join(_output_dir,"{sample}.rv.bed")
        shell:
            "python {params.dir_script}/1.1_bed2fr.py -b {input} {params.ignore_chr}"    

rule bed2bdg_fw:
    input:
        os.path.join(_output_dir,"{sample}.fw.bed")
    threads:
        _thread
    params:
        gl=config["genomelen"]
    output:
        temp(os.path.join(_output_dir,"{sample}.fw.bdg"))
    shell:
        "bedtools genomecov -bg -i {input} -g {params.gl} > {output}"

rule bed2bdg_rv:
    input:
        os.path.join(_output_dir,"{sample}.rv.bed")
    threads:
        _thread
    params:
        gl=config["genomelen"]
    output:
        temp(os.path.join(_output_dir,"{sample}.rv.bdg"))
    shell:
        "bedtools genomecov -bg -i {input} -g {params.gl} > {output}"

rule bdg_sort_fw:
    input:
        fw=os.path.join(_output_dir,"{sample}.fw.bdg")
    threads:
        _thread
    output:
        temp(os.path.join(_output_dir,"{sample}.fw.sorted.bdg"))
    shell:
        "bedtools sort -i {input.fw} > {output}"

rule bdg_sort_rv:
    input:
        rv=os.path.join(_output_dir,"{sample}.rv.bdg")
    threads:
        _thread
    output:
        temp(os.path.join(_output_dir,"{sample}.rv.sorted.bdg"))
    shell:
        "bedtools sort -i {input.rv} > {output}"

if _normalize == "True":
    rule bdg_normalize_fw:
        input:
            bdg=os.path.join(_output_dir,"{sample}.fw.sorted.bdg"),
            bed=os.path.join(_output_dir,"{sample}.fw.bed")
        threads:
            _thread
        params:
            dir_script=config["utility_dir"]
        output:
            temp(os.path.join(_output_dir,"{sample}.fw.scaled.bdg"))
        shell:
            "python {params.dir_script}/1.3_bdg_normalize_v4.0.py --bdg {input.bdg} --bed {input.bed}"
    
    rule bdg_normalize_rv:
        input:
            bdg=os.path.join(_output_dir,"{sample}.rv.sorted.bdg"),
            bed=os.path.join(_output_dir,"{sample}.rv.bed")
        threads:
            _thread
        params:
            dir_script=config["utility_dir"]
        output:
            temp(os.path.join(_output_dir,"{sample}.rv.scaled.bdg"))
        shell:
            "python {params.dir_script}/1.3_bdg_normalize_v4.0.py --bdg {input.bdg} --bed {input.bed}"
    
    rule bdg2bw_fw:
        input:
            os.path.join(_output_dir,"{sample}.fw.scaled.bdg")
        threads:
            _thread
        params:
            gl=config["genomelen"],
            dir_script=config["utility_dir"]
        output:
            os.path.join(_output_dir,"{sample}.fw.scaled.bw")
        shell:
            "{params.dir_script}/bedGraphToBigWig {input} {params.gl} {output}"
    
    rule bdg2bw_rv:
        input:
            os.path.join(_output_dir,"{sample}.rv.scaled.bdg")
        threads:
            _thread
        params:
            gl=config["genomelen"],
            dir_script=config["utility_dir"]
        output:
            os.path.join(_output_dir,"{sample}.rv.scaled.bw")
        shell:
            "{params.dir_script}/bedGraphToBigWig {input} {params.gl} {output}"
    
    rule bwAdd:
        input:
            fw=os.path.join(_output_dir,"{sample}.fw.scaled.bw"),
            rv=os.path.join(_output_dir,"{sample}.rv.scaled.bw")
        threads:
            _thread
        output:
            os.path.join(_output_dir,"{sample}." + _BinSize + ".add.bdg")
        shell:
            """
            bigwigCompare --binSize {_BinSize} -p {threads} --verbose -o {output} \
            --outFileFormat bedgraph --fixedStep \
            --bigwig1 {input.fw} \
            --bigwig2 {input.rv} \
            --operation add 
            """
else:
    rule bdg_reverse_rv:
        input:
            os.path.join(_output_dir,"{sample}.rv.sorted.bdg")
        threads:
            _thread
        output:
            temp(os.path.join(_output_dir,"{sample}.rv.sorted_r.bdg"))
        shell:
            "awk -F '\t' -v OFS='\t' '{{$4=-$4; print}}' {input} > {output}"
    
    rule bdg2bw_fw:
        input:
            os.path.join(_output_dir,"{sample}.fw.sorted.bdg")
        threads:
            _thread
        params:
            gl=config["genomelen"],
            dir_script=config["utility_dir"]
        output:
            os.path.join(_output_dir,"{sample}.fw.raw.bw")
        shell:
            "{params.dir_script}/bedGraphToBigWig {input} {params.gl} {output}"
    
    rule bdg2bw_rv:
        input:
            os.path.join(_output_dir,"{sample}.rv.sorted_r.bdg")
        threads:
            _thread
        params:
            gl=config["genomelen"],
            dir_script=config["utility_dir"]
        output:
            os.path.join(_output_dir,"{sample}.rv.raw.bw")
        shell:
            "{params.dir_script}/bedGraphToBigWig {input} {params.gl} {output}"




