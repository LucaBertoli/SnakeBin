# Snakefile

rule markduplicates_bam_from_merged_fastq_raw_fastq:
    input:
        start_sorted_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/start_sorted.bam",
    output:
        mark_dup_bam=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.bam"),
        metric_dup=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/duplicates.txt",
    shell:
        """
        java -Xmx8G -Djava.io.tmpdir=/home/temp -jar {config[path_picard]} MarkDuplicates \\
        I={input.start_sorted_bam} \\
        O={output.mark_dup_bam} \\
        M={output.metric_dup} \\
        REMOVE_DUPLICATES=true \\
        VALIDATION_STRINGENCY=SILENT \\
        CREATE_INDEX=true \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
        """
