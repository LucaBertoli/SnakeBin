# Snakefile

import re
import os


def get_fq_names_CGF_ID(wildcards):
    fastq_R1 = []
    fastq_R2 = []
    for filename in os.listdir(config['results_folder'] + "/" + wildcards.sample):
        if re.search(r'R1', filename) and filename.endswith('.fastq.gz') and 'trimmed' not in filename:
            fastq_R1.append(os.path.join(config['results_folder'] + "/" + wildcards.sample, filename))
        elif re.search(r'R2', filename) and filename.endswith('.fastq.gz') and 'trimmed' not in filename:
            fastq_R2.append(os.path.join(config['results_folder'] + "/" + wildcards.sample, filename))
    
    return fastq_R1, fastq_R2

rule fastqc_R1:
    input:
        fastq=lambda wildcards: get_fq_names_CGF_ID(wildcards)[0]
    output:
        html=config['results_folder'] + "/{sample}/fastqc/{sample}_R1.html",
        zip=config['results_folder'] + "/{sample}/fastqc/{sample}_R1_fastqc.zip"
    params:
        extra="--extract"
    log:
        log=config['results_folder'] + "/{sample}/logs/fastqc/{sample}_R1.log"
    threads:
        config['fastqc_threads']
    resources:
        mem_mb=4096
    wrapper:
        "v5.5.0/bio/fastqc"


rule fastqc_R2:
    input:
        fastq=lambda wildcards: get_fq_names_CGF_ID(wildcards)[1]
    output:
        html=config['results_folder'] + "/{sample}/fastqc/{sample}_R2.html",
        zip=config['results_folder'] + "/{sample}/fastqc/{sample}_R2_fastqc.zip"
    params:
        extra="--extract"
    log:
        log=config['results_folder'] + "/{sample}/logs/fastqc/{sample}_R2.log"
    threads:
        config['fastqc_threads']
    resources:
        mem_mb=4096
    wrapper:
        "v5.5.0/bio/fastqc"
