# Snakefile

rule baserecalibrator_downsampled_fastq_per_bin:
    input:
        bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.bam",
    output:
        recal_table=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/recal_data.table",
    log:
        config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/recal.log",
    shell: 
        """
        gatk BaseRecalibrator \\
        -R {config[reference_fasta]} \\
        -I {input.bam} \\
        --known-sites {config[mills_vcf]} \\
        --known-sites {config[dbsnp_vcf]} \\
        -O {output.recal_table}
        """

rule applyBQSR_downsampled_fastq_per_bin:
    input:
        recal_table=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/recal_data.table",
        bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.bam",
    output:
        recal_bam=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam")
    log:
        config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/recal.log",
    shell:
        """
        gatk ApplyBQSR \\
        --static-quantized-quals 10 \\
        --static-quantized-quals 20 \\
        --static-quantized-quals 30 \\
        -R {config[reference_fasta]} \\
        -I {input.bam} \\
        -O {output.recal_bam} \\
        --bqsr-recal-file {input.recal_table}
        """

rule samtools_index_recal_bam_downsampled_fastq_per_bin:
    input:
        recal_sorted_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam"
    output:
        recal_sorted_bam_bai=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam.bai")
    params:
        samtools_threads=config['samtools_threads']
    shell:
        """
        {config[path_samtools]} index \\
        --threads {params.samtools_threads} \\
        -o {output.recal_sorted_bam_bai} \\
        {input.recal_sorted_bam}
        """

