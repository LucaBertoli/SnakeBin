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