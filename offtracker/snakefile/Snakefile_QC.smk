# 更新记录:
# 2022.05.04. v1.0:    初步运行, fastp + multiqc
# 2024.01.17. v2.0:    翻新结构，匹配 X_NGS 框架

# 参数列表
configfile: "config.yaml"

### config['files_R1'], config['files_R2'] 为 dict型

# # fastq 信息
_files_R1 = config['files_R1'] # dict型, key 为 sample
_files_R2 = config['files_R2'] # dict型, key 为 sample
# # 输入输出文件夹
# config['input_dir']
_output_dir = config["output_dir"]
# # 运行参数
_thread = config['thread']
# config['utility_dir']

import os

############################
# conditional output_files #
############################
output_HT = expand( os.path.join(_output_dir,"{sample}_fastp.html"), sample=_files_R1)
output_JS = expand( os.path.join(_output_dir,"{sample}_fastp.json"), sample=_files_R1)
output_MQC = os.path.join(_output_dir,"MultiQC_Report_Raw.html")
output_R1 = expand( os.path.join(_output_dir,"{sample}_trimmed_1.fq.gz"), sample=_files_R1) # dict 会自动迭代 keys
output_R2 = expand( os.path.join(_output_dir,"{sample}_trimmed_2.fq.gz"), sample=_files_R1)

output_files = output_HT + output_JS + [output_MQC] + output_R1 + output_R2

rule all:
    input:
        output_files

#######################
## fastp and multiQC ##
#######################
rule QCtrim:
    input:
        R1=lambda w: _files_R1[w.sample],
        R2=lambda w: _files_R2[w.sample]
    threads:
        _thread
    output:
        R1=os.path.join(_output_dir,"{sample}_trimmed_1.fq.gz"),
        R2=os.path.join(_output_dir,"{sample}_trimmed_2.fq.gz"),
        HT=os.path.join(_output_dir,"{sample}_fastp.html"),
        JS=os.path.join(_output_dir,"{sample}_fastp.json")
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} \
        -h {wildcards.sample}_fastp.html -j {wildcards.sample}_fastp.json \
        --length_required 10 --thread {threads} --detect_adapter_for_pe --disable_quality_filtering
        """

rule multiqc:
    input:
        expand( os.path.join(_output_dir,"{sample}_fastp.html"), sample=_files_R1 )
    threads:
        _thread
    output:
        os.path.join(_output_dir,"MultiQC_Report_Raw.html")
    shell:
        "multiqc {_output_dir} -n MultiQC_Report_Raw --outdir {_output_dir}"
