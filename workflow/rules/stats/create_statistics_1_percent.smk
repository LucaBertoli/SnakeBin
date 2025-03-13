# Snakefile

rule stats_downsampled_ds_1_perc:
    input:
        duplicates=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/duplicates.txt",
        flagstat=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/flagstat.txt",
    output:
        duplicates_tsv=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/stats/duplicates_clip.tsv",
        flagstat_tsv=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/stats/flagstat_tsv",
    shell:
        """
        python3 workflow/scripts/metricsParser.py {input.duplicates} > {output.duplicates_tsv}
        python3 workflow/scripts/metricParser.py {input.flagstat} > {output.flagstat_tsv}
        """


rule aggregate_stats_ds_1_perc:
    input:
        duplicates_tsv=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/stats/duplicates_clip.tsv",
        flagstat_tsv=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/stats/flagstat.tsv",
    output:
        sample_stats=config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/stats.tsv",
    shell:
        """
        echo -e "{wildcards.bin_start}-{wildcards.bin_end}\\t$(paste -d'\\t' \\
            {input.flagstat_tsv} \\
            {input.duplicates_tsv})" \\
            > {output.sample_stats}
        """

rule merge_stats_ds_1_perc:
    input:
        sample_stats=expand(config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/stats.tsv", zip , bin_start=l_bin_start, bin_end=l_bin_end),
    output:
        aggregate_stats=config['results_folder']+"/merged_downsampled_1_percent_sequencial/statistics_ds.tsv"
    shell:
        """
        cat {input.sample_stats} >> {output.aggregate_stats}
        """

