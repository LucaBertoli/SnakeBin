# Snakefile

import re
import os


rule fastqc_R1_ds_1_perc_merged:
    input:
        fastq=config['results_folder'] + "/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/merged.R1.fastq.gz"
    output:
        html=config['results_folder'] + "/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/fastqc/merged_R1.html",
        zip=config['results_folder'] + "/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/fastqc/merged_R1_fastqc.zip"
    params:
        extra="--extract"
    threads:
        config['fastqc_threads']
    resources:
        mem_mb=4096
    wrapper:
        "v5.5.0/bio/fastqc"


rule fastqc_R2_ds_1_perc_merged:
    input:
        fastq=config['results_folder'] + "/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/merged.R2.fastq.gz"
    output:
        html=config['results_folder'] + "/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/fastqc/merged_R2.html",
        zip=config['results_folder'] + "/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/fastqc/merged_R2_fastqc.zip"
    params:
        extra="--extract"
    threads:
        config['fastqc_threads']
    resources:
        mem_mb=4096
    wrapper:
        "v5.5.0/bio/fastqc"
