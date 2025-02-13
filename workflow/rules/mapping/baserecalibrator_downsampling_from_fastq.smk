# Snakefile

rule baserecalibrator_downsampled_fastq:
    input:
        mark_dup_bam=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.bam",
    output:
        recal_table=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/recal_data.table",
    shell: 
        """
        gatk BaseRecalibrator \\
        -R {config[reference_fasta]} \\
        -I {input.mark_dup_bam} \\
        --known-sites {config[mills_vcf]} \\
        --known-sites {config[dbsnp_vcf]} \\
        -O {output.recal_table}
        """

rule applyBQSR_downsampled_fastq:
    input:
        recal_table=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/recal_data.table",
        mark_dup_bam=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.bam",
    output:
        recal_bam=temp(config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam")
    shell:
        """
        gatk ApplyBQSR \\
        --static-quantized-quals 10 \\
        --static-quantized-quals 20 \\
        --static-quantized-quals 30 \\
        -R {config[reference_fasta]} \\
        -I {input.mark_dup_bam} \\
        -O {output.recal_bam} \\
        --bqsr-recal-file {input.recal_table}
        """

rule samtools_index_recal_downsampled_fastq:
    input:
        recal_bam=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam",
    output:
        recal_bam_bai=temp(config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam.bai"),
    params:
        samtools_threads=config['samtools_threads']
    shell:
        """
        {config[path_samtools]} index \\
        --threads {params.samtools_threads} \\
        -o {output.recal_bam_bai} \\
        {input.recal_bam}
        """

