# Snakefile


rule exome_coverage:
    input:
        clipped_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        target_bed=config["target_regions"]["design"]["bed"]
    output:
        coverage_hist_design=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-capture.hist.coverage.gz",
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -abam {input.clipped_bam} \\
        -b {input.target_bed} \\
        | gzip > {output.coverage_hist_design}
        """

rule coverage_pass_region:
    input:
        callable_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/callable.bed",
        target_bed=config["target_regions"]["design"]["bed"]
    output:
        coverage_callable_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable.bed"
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -a {input.callable_bed} \\
        -b {input.target_bed} \\
        > {output.coverage_callable_bed}
        """


rule coverage_pass_dp10_region:
    input:
        callable_dp10_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/callable_DP10.bed",
        target_bed=config["target_regions"]["design"]["bed"]
    output:
        coverage_callable_dp10_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10.bed"
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -a {input.callable_dp10_bed} \\
        -b {input.target_bed} \\
        > {output.coverage_callable_dp10_bed}
        """


rule gene_coverage:
    input:
        coverage_hist=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-capture.hist.coverage.gz",
        coverage_callable_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable.bed",
        coverage_callable_dp10_bed=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10.bed",
    output:
        tsv_design=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10-stats.tsv"
    shell:
        """
        python workflow/scripts/geneCoverage.py \\
        {input.coverage_hist} \\
        {input.coverage_callable_bed} \\
        {input.coverage_callable_dp10_bed}
        """


