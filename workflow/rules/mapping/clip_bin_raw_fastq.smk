# Snakefile
import os

rule clip_bam_bin_raw_fastq:
    input:
        recalibrated_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam",
    output:
        clipped_bam=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam"),
    shell:
        """
        {config[path_bamutil]} clipOverlap \\
        --in {input.recalibrated_bam} \\
        --out {output.clipped_bam}
        """

rule samtools_index_clipped_bin_raw_fastq:
    input:
        clipped_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
    output:
        clipped_bam_index=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai"),
    params:
        samtools_threads=config['samtools_threads']
    shell:
        """
        {config[path_samtools]} index \\
        --threads {params.samtools_threads} \\
        -o {output.clipped_bam_index} \\
        {input.clipped_bam}
        """