# Snakefile
import glob as glob
import os

def get_fq(samples):
    # Code that returns a list of fastq files for read 1 and read 2 based on *wildcards.sample*
    l_fastq_R1 = []
    l_fastq_R2 = []
    for sample in samples:
        l_fastq_R1.append([f for f in glob.glob(config['results_folder'] + '/' + sample + '/*R1*fastq.gz') if "trimmed" not in f])
        l_fastq_R2.append([f for f in glob.glob(config['results_folder'] + '/' + sample + '/*R2*fastq.gz') if "trimmed" not in f])
    return l_fastq_R1, l_fastq_R2

l_fastq_R1, l_fastq_R2 = get_fq(samples)

rule merge_raw_fastq:
    input:
        fastq_R1=l_fastq_R1, #fastq_grezzi
        fastq_R2=l_fastq_R2 #fastq_grezzi
    output:
        merged_sample_R1=temp(config['results_folder']+"/merged_raw_fastq/merged.R1.fastq.gz"),
        merged_sample_R2=temp(config['results_folder']+"/merged_raw_fastq/merged.R2.fastq.gz")
    shell:
        """
        cat {input.fastq_R1} > {output.merged_sample_R1}
        cat {input.fastq_R2} > {output.merged_sample_R2}
        """

