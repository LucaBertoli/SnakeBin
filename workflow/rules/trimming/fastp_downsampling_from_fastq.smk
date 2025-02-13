# Snakefile
import glob as glob
import os

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    return sorted(glob.glob(config['results_folder'] +'/'+ wildcards.sample + '/*R1*fastq.gz'))

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    return sorted(glob.glob(config['results_folder'] +'/'+ wildcards.sample + '/*R2*fastq.gz'))

rule trimming:
    input:
        fastq_R1=get_fq1, #fastq_grezzi
        fastq_R2=get_fq2 #fastq_grezzi
    output:
        trimmed_R1=temp(config['results_folder']+"/{sample}/{sample}.trimmed.R1.fastq.gz"), #trimmed1.fastq.gz
        trimmed_R2=temp(config['results_folder']+"/{sample}/{sample}.trimmed.R2.fastq.gz") #trimmed2.fastq.gz
    threads: 
        config['fastp_threads']
    params:
        threads=config['fastp_threads']
    log:
       config['results_folder']+"/{sample}/logs/fastp/{sample}.log" 
    shell:
        """
        {config[path_fastp]} \\
        -i {input.fastq_R1[0]} \\
        -I {input.fastq_R2[0]} \\
        -o {output.trimmed_R1} \\
        -O {output.trimmed_R2} \\
        --failed_out {config[results_folder]}/{wildcards.sample}/{wildcards.sample}_failed.fastq.gz \\
        -h {config[results_folder]}/{wildcards.sample}/{wildcards.sample}_adapters_removal_report.html \\
        -j {config[results_folder]}/{wildcards.sample}/{wildcards.sample}_fastp.json \\
        -w {params.threads}
        """

