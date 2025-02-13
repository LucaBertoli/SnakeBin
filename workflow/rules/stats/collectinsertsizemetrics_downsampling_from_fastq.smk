# Snakefile

rule collect_insert_metrics_downsampling_from_fastq:
    input:
        clipped_bam=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
    output:
        insert_hist=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.hist.pdf",
        insert_hist_out=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.output"
    shell:
        """
        java -jar {config[path_picard]} CollectInsertSizeMetrics \\
        I={input.clipped_bam} \\
        H={output.insert_hist} \\
        O={output.insert_hist_out} \\
        AS=true \\
        VALIDATION_STRINGENCY=SILENT
        """
