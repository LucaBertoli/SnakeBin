# Snakefile

rule samtools_flagstat_recall_downsampling_from_fastq:
    input:
        start_sorted_bam=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        start_sorted_bam_bai=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai",
    output:
        flagstat=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/flagstat_recall.txt",
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

