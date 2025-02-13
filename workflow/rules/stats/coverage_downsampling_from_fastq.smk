# Snakefile


rule exome_coverage_downsampling_from_fastq:
    input:
        clipped_bam=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
    output:
        coverage_hist_design=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-capture.hist.coverage.gz",
    params:
        target_bed=lambda wildcards: config["target_regions"]["design"]["bed"]
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -abam {input.clipped_bam} \\
        -b {params.target_bed} \\
        | gzip > {output.coverage_hist_design}
        """

rule coverage_pass_region_downsampling_from_fastq:
    input:
        callable_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable.bed"
    output:
        coverage_callable_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable.bed"
    params:
        target_bed=lambda wildcards: config["target_regions"]["design"]["bed"]
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -a {input.callable_bed} \\
        -b {params.target_bed} \\
        > {output.coverage_callable_bed}
        """


rule coverage_pass_dp10_region_downsampling_from_fastq:
    input:
        callable_dp10_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_DP10.bed"
    output:
        coverage_callable_dp10_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10.bed"
    params:
        target_bed=lambda wildcards: config["target_regions"]["design"]["bed"]
    shell:
        """
        {config[path_bedtools]} coverage \\
        -hist \\
        -a {input.callable_dp10_bed} \\
        -b {params.target_bed} \\
        > {output.coverage_callable_dp10_bed}
        """


rule gene_coverage_downsampling_from_fastq:
    input:
        coverage_hist=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-capture.hist.coverage.gz",
        coverage_callable_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable.bed",
        coverage_callable_dp10_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10.bed",
    output:
        tsv_design=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10-stats.tsv"
    shell:
        """
        python workflow/scripts/geneCoverage.py \\
        {input.coverage_hist} \\
        {input.coverage_callable_bed} \\
        {input.coverage_callable_dp10_bed}
        """


