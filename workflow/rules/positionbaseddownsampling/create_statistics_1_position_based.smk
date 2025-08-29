# Snakefile

rule stats_downsampled_1_position_based_dup:
    input:
        duplicates=config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/duplicates.txt",
    output:
        duplicates_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/stats/duplicates_clip.tsv",
    shell:
        """
        python3 workflow/scripts/metricsParser.py {input.duplicates} > {output.duplicates_tsv}
        """

rule stats_downsampled_1_position_based_flagstat:
    input:
        flagstat=config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/flagstat.txt",
    output:
        flagstat_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/stats/flagstat.tsv",
    shell:
        """
        python3 workflow/scripts/metricsParser.py {input.flagstat} > {output.flagstat_tsv}
        """

rule aggregate_stats_1_position_based:
    input:
        duplicates_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/stats/duplicates_clip.tsv",
        flagstat=config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/stats/flagstat.tsv",
    output:
        sample_stats=config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/stats.tsv",
    shell:
        """
        echo -e "{wildcards.bin_start}-{wildcards.bin_end}\\t$(paste -d'\\t' \\
            {input.flagstat} \\
            {input.duplicates_tsv})" \\
            > {output.sample_stats}
        """

rule merge_stats_1_position_based:
    input:
        sample_stats=expand(config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/stats.tsv", zip , bin_start=l_bin_start, bin_end=l_bin_end),
    output:
        aggregate_stats=config['results_folder']+"/merged_raw_fastq/statistics_ds_1_pos.tsv"
    shell:
        """
        cat {input.sample_stats} >> {output.aggregate_stats}
        """

