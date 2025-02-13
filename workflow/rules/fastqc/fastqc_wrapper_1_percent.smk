# Snakefile

import re
import os


#def get_fq_names_CGF_ID_ds_1_perc(wildcards):
#    fastq_R1 = []
#    fastq_R2 = []
#    for filename in os.listdir(config['results_folder'] + "/" + wildcards.sample + "/split_fastq/bam_to_fastq_" + wildcards.bin_start + "-" + wildcards.bin_end):
#        if re.search(r'R1', filename) and filename.endswith('.fastq.gz') and 'trimmed' not in filename:
#            fastq_R1.append(os.path.join(config['results_folder'] + "/" + wildcards.sample, filename))
#        elif re.search(r'R2', filename) and filename.endswith('.fastq.gz') and 'trimmed' not in filename:
#            fastq_R2.append(os.path.join(config['results_folder'] + "/" + wildcards.sample, filename))
#    
#    return fastq_R1, fastq_R2


rule fastqc_R1_ds_1_perc:
    input:
        fastq=config['results_folder'] + "/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R1.fastq.gz"
    output:
        html=config['results_folder'] + "/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/fastqc/{sample}_R1.html",
        zip=config['results_folder'] + "/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/fastqc/{sample}_R1_fastqc.zip"
    params:
        extra="--extract"
    threads:
        config['fastqc_threads']
    resources:
        mem_mb=4096
    wrapper:
        "v5.5.0/bio/fastqc"


rule fastqc_R2_ds_1_perc:
    input:
        fastq=config['results_folder'] + "/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R2.fastq.gz"
    output:
        html=config['results_folder'] + "/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/fastqc/{sample}_R2.html",
        zip=config['results_folder'] + "/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/fastqc/{sample}_R2_fastqc.zip"
    params:
        extra="--extract"
    threads:
        config['fastqc_threads']
    resources:
        mem_mb=4096
    wrapper:
        "v5.5.0/bio/fastqc"
