# Snakefile

rule samtools_flagstat_position_based_per_sample:
    input:
        start_sorted_bam=config['results_folder']+"/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/start_sorted.bam",
        start_sorted_bam_bai=config['results_folder']+"/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/start_sorted.bam.bai",
    output:
        flagstat=config['results_folder']+"/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/flagstat.txt",
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