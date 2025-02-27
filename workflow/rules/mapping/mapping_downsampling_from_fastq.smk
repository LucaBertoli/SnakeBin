# Snakefile
import glob as glob
from datetime import datetime
import os

def get_fq1_CGF_ID(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    fq1_files = sorted(glob.glob(config['results_folder'] +'/'+ wildcards.sample + '/*R1*fastq.gz'))
    if not fq1_files:
        raise ValueError(f"No fastq files found for sample {wildcards.sample}")
    return fq1_files[0].split('/')[-1:][0].split('_')[0]

def get_time(wildcards):
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y-%H:%M:%S")
    return dt_string

rule bwa_mem2_mem:
    input:
        trimmed1=config['results_folder']+"/{sample}/{sample}.trimmed.R1.fastq.gz", #trimmed R1
        trimmed2=config['results_folder']+"/{sample}/{sample}.trimmed.R2.fastq.gz" #trimmed R2     
    output:
        start_sorted_bam=temp(config['results_folder']+"/{sample}/{sample}_start_sorted.bam"), #bwa-mapped samples
    threads: 
        max(config['bwa_threads'], config['samtools_threads'])
    params:
        bwa_threads=config['bwa_threads'],
        samtools_threads=config['samtools_threads'],
        CGF_ID_R1=get_fq1_CGF_ID,
    shell:
        """
        {config[path_bwa]} mem \\
        -R '@RG\\tID:{params.CGF_ID_R1}\\tPU:lane\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tCN:CGF-ddlab\\tPL:ILLUMINA' \\
        -t {params.bwa_threads} \\
        {config[reference_fasta]} \\
        {input.trimmed1} \\
        {input.trimmed2} \\
        | {config[path_samtools]} sort \\
        --threads {params.samtools_threads} \\
        -o {output.start_sorted_bam}
        """

rule samtools_index:
    input:
        start_sorted_bam=config['results_folder']+"/{sample}/{sample}_start_sorted.bam"
    output:
        start_sorted_bam_bai=temp(config['results_folder']+"/{sample}/{sample}_start_sorted.bam.bai"),
    threads: 
        config['samtools_threads']
    params:
        samtools_threads=config['samtools_threads']
    shell: 
        """
        {config[path_samtools]} index \\
        --threads {params.samtools_threads} \\
        -o {output.start_sorted_bam_bai} \\
        {input.start_sorted_bam}
        """
