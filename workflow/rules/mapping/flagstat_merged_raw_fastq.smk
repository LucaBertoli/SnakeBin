# Snakefile

rule samtools_flagstat_merged_fastq:
    input:
        start_sorted_bam=config['results_folder']+"/merged_raw_fastq/start_sorted.bam",
        start_sorted_bam_bai=config['results_folder']+"/merged_raw_fastq/start_sorted.bam.bai",
    output:
        flagstat=config['results_folder']+"/merged_raw_fastq/flagstat.txt",
    threads: 
        config['samtools_threads']
    params:
        samtools_threads=config['samtools_threads']
    shell:
        """
        {config[path_samtools]} flagstat \\
        --threads {params.samtools_threads} \\
        {input.start_sorted_bam} > {output.flagstat}
        """ 
