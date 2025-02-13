# Snakefile
import csv
import os

def near_distance(wildcards):
    with open(config['results_folder'] + "/merged_raw_fastq/split_bam_merged/" + wildcards.coverage_downsampling_levels + "/insert_" + wildcards.bin_start + "-" + wildcards.bin_end + "/alignment.dedup.recal.clip.output", newline='') as csvfile:
        tsv_reader = csv.reader(csvfile, delimiter='\t')
        is_mean_insert_size_row = False
        for row in tsv_reader:
            if is_mean_insert_size_row:
                mean_insert_size = int(float(row[5]))
                return mean_insert_size
            if "MEAN_INSERT_SIZE" in row:
                is_mean_insert_size_row = True
    return None

rule hsmetrics_bin_raw_fastq_final:
    input:
        clipped_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        insert_hist_out=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.output"
    output:
        output_metrics=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_HsMetrics.txt",
        output_target=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_PER_TARGET_COVERAGE.txt",
        output_base=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_PER_BASE_COVERAGE.txt",
    params:
        near_distance=near_distance,
        bait_intervals=lambda wildcards: config["target_regions"]["design"]["intervals"],
        target_intervals=lambda wildcards: config["target_regions"]["design"]["intervals"],
        levels=lambda wildcards: config["coverage_downsampling_levels"][wildcards.coverage_downsampling_levels]
    shell:
        """
        java -Xmx80G -jar {config[path_picard]} CollectHsMetrics \\
        I={input.clipped_bam} \\
        R={config[reference_fasta]} \\
        O={output.output_metrics} \\
        BAIT_INTERVALS={params.bait_intervals} \\
        TARGET_INTERVALS={params.target_intervals} \\
        PER_TARGET_COVERAGE={output.output_target} \\
        PER_BASE_COVERAGE={output.output_base} \\
        VALIDATION_STRINGENCY=SILENT \\
        NEAR_DISTANCE={params.near_distance} 
        """