# Snakefile

rule per_sample_stats_bin_raw_fastq:
    input:
        insert_size=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.output",
        flagstat_recall=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/flagstat_recall.txt",
        target_stats=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10-stats.tsv",
        on_near_off_target=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_HsMetrics.txt",
        unif_cov_target=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/unif_of_coverage.txt",

    output:
        insert_size_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/insert_des_size_clip.tsv",
        flagstat_recall_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/map_dedup_clip.tsv",
        target_stats_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/Statistics_clip.tsv",
        on_near_off_target_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/HsMetrics_clip.tsv",
        unif_cov_target_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/unif_of_coverage.tsv",
    shell:
        """
        python3 workflow/scripts/metricsParser.py {input.insert_size} > {output.insert_size_tsv}
        python3 workflow/scripts/metricsParser.py {input.flagstat_recall} > {output.flagstat_recall_tsv}
        python3 workflow/scripts/metricsParser.py {input.target_stats} > {output.target_stats_tsv}
        python3 workflow/scripts/metricsParser.py {input.on_near_off_target} > {output.on_near_off_target_tsv}
        grep "Uniformity of coverage" {input.unif_cov_target} | sed 's/^.*: //' > {output.unif_cov_target_tsv}
        """


rule per_sample_aggregate_stats_bin_raw_fastq:
    input:
        insert_size_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/insert_des_size_clip.tsv",
        flagstat_recall_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/map_dedup_clip.tsv",
        target_stats_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/Statistics_clip.tsv",
        on_near_off_target_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/HsMetrics_clip.tsv",
        unif_cov_target_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats/unif_of_coverage.tsv",
    output:
        sample_stats=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/stats.tsv",
    shell:
        """
        echo -e "{wildcards.bin_start}-{wildcards.bin_end}\\t{wildcards.coverage_downsampling_levels}\\t$(paste -d'\\t' \\
            {input.insert_size_tsv} \\
            {input.flagstat_recall_tsv} \\
            {input.target_stats_tsv} \\
            {input.on_near_off_target_tsv} \\
            {input.unif_cov_target_tsv})" \\
            > {output.sample_stats}
        """

rule aggregate_stats_bin_raw_fastq:
    input:
        sample_stats=expand(expand(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{levels}/insert_{bin_start}-{bin_end}/stats.tsv", zip , bin_start=l_bin_start, bin_end=l_bin_end, allow_missing=True), levels=coverage_downsampling_levels),
    output:
        aggregate_stats=config['results_folder']+"/merged_raw_fastq/statistics_ds.tsv"
    shell:
        """
        cat {input.sample_stats} >> {output.aggregate_stats}
        """

