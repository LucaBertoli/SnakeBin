#Snakemake

rule flagstat_1_position_based:
    input:
        start_bam = config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/start_sorted.bam"
    output:
        flagstat = config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/flagstat.txt"
    threads: 
        config['samtools_threads']
    params:
        samtools_threads = config['samtools_threads']
    shell:
        """
        samtools flagstat {input.start_bam} --threads {params.samtools_threads} > {output.flagstat}
        """

rule flagstat_1_position_based_dedup:
    input:
        mark_dup_bam = config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/alignment.dedup.bam"
    output:
        flagstat_dedup = config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/flagstat_dedup.txt"
    threads: 
        config['samtools_threads']
    params:
        samtools_threads = config['samtools_threads']
    shell:
        """
        samtools flagstat {input.mark_dup_bam} --threads {params.samtools_threads} > {output.flagstat_dedup}
        """