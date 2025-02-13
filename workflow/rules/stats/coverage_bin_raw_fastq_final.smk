# Snakefile


rule exome_coverage_bin_raw_fastq_final:
    input:
        clipped_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        target_bed=config["target_regions"]["design"]["bed"]
    output:
        coverage_hist_design=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-capture.hist.coverage.gz",
    params:
        levels=lambda wildcards: config["coverage_downsampling_levels"][wildcards.coverage_downsampling_levels]
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -abam {input.clipped_bam} \\
        -b {input.target_bed} \\
        | gzip > {output.coverage_hist_design}
        """

rule coverage_pass_region_bin_raw_fastq_final:
    input:
        callable_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/callable.bed",
        target_bed=config["target_regions"]["design"]["bed"]
    output:
        coverage_callable_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable.bed"
    params:
        levels=lambda wildcards: config["coverage_downsampling_levels"][wildcards.coverage_downsampling_levels]
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -a {input.callable_bed} \\
        -b {input.target_bed} \\
        > {output.coverage_callable_bed}
        """


rule coverage_pass_dp10_region_bin_raw_fastq_final:
    input:
        callable_dp10_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/callable_DP10.bed",
        target_bed=config["target_regions"]["design"]["bed"]
    output:
        coverage_callable_dp10_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10.bed"
    params:
        levels=lambda wildcards: config["coverage_downsampling_levels"][wildcards.coverage_downsampling_levels]
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -a {input.callable_dp10_bed} \\
        -b {input.target_bed} \\
        > {output.coverage_callable_dp10_bed}
        """


rule gene_coverage_bin_raw_fastq_final:
    input:
        coverage_hist=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-capture.hist.coverage.gz",
        coverage_callable_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable.bed",
        coverage_callable_dp10_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10.bed",
    output:
        tsv_design=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10-stats.tsv"
    params:
        levels=lambda wildcards: config["coverage_downsampling_levels"][wildcards.coverage_downsampling_levels]
    shell:
        """
        python workflow/scripts/geneCoverage.py \\
        {input.coverage_hist} \\
        {input.coverage_callable_bed} \\
        {input.coverage_callable_dp10_bed}
        """


