# Snakefile

rule samtools_flagstat_1_percent:
    input:
        start_sorted_bam=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/start_sorted.bam",
        start_sorted_bam_bai=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/start_sorted.bam.bai",
    output:
        flagstat=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/flagstat.txt",
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