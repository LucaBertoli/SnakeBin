# Snakefile


rule merge_reads_from_fastq:
    input:
        downsampled_fastq_R1=expand(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{{bin_start}}-{{bin_end}}/{sample}_R1.fastq.gz",sample=samples),
        downsampled_fastq_R2=expand(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{{bin_start}}-{{bin_end}}/{sample}_R2.fastq.gz",sample=samples),
    output:
        merged_fastq_R1=temp(config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/merged.R1.fastq.gz"),
        merged_fastq_R2=temp(config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/merged.R2.fastq.gz"),
    shell:
        """
        cat {input.downsampled_fastq_R1} > {output.merged_fastq_R1}
        cat {input.downsampled_fastq_R2} > {output.merged_fastq_R2}
        """
