# Snakefile
import csv
import os

def near_distance(wildcards):
    with open(config['results_folder'] + "/merged_downsampled_bam/" + wildcards.million_fragment_downsampling + "/insert_" + wildcards.bin_start + "-" + wildcards.bin_end + "/alignment.dedup.recal.clip.output", newline='') as csvfile:
        tsv_reader = csv.reader(csvfile, delimiter='\t')
        is_mean_insert_size_row = False
        for row in tsv_reader:
            if is_mean_insert_size_row:
                mean_insert_size = int(float(row[5]))
                return mean_insert_size
            if "MEAN_INSERT_SIZE" in row:
                is_mean_insert_size_row = True
    return None

rule hsmetrics_downsampling_from_bam:
    input:
        clipped_bam=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        insert_hist_out=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.output"
    output:
        output_metrics=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_HsMetrics.txt",
        output_target=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_PER_TARGET_COVERAGE.txt",
        output_base=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_PER_BASE_COVERAGE.txt",
    params:
        near_distance=near_distance,
        bait_intervals=lambda wildcards: config["target_regions"]["design"]["intervals"],
        target_intervals=lambda wildcards: config["target_regions"]["design"]["intervals"]
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

# Definizione delle regole con utilizzo dei wildcards
#rule hsmetrics:
#    input:
#        clipped_bam=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip.bam",
#        insert_hist_out=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip.output"
#    output:
#        output_metrics=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['design'][:-4])+"_HsMetrics.txt",
#        output_target=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['design'][:-4])+"_PER_TARGET_COVERAGE.txt",
#        output_base=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['design'][:-4])+"_PER_BASE_COVERAGE.txt",
#    params:
#        near_distance=near_distance,
#    shell:
#        """
#        java -Xmx80G -jar {config[path_picard]} CollectHsMetrics \\
#        I={input.clipped_bam} \\
#        R={config[reference_fasta]} \\
#        O={output.output_metrics} \\
#        BAIT_INTERVALS= {config[design_intervals]} \\
#        TARGET_INTERVALS= {config[design_intervals]} \\
#        PER_TARGET_COVERAGE= {output.output_target} \\
#        PER_BASE_COVERAGE= {output.output_base} \\
#        VALIDATION_STRINGENCY=SILENT \\
#        NEAR_DISTANCE={params.near_distance}
#        """
#
#rule hsmetrics_refseq:
#    input:
#        clipped_bam=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip.bam",
#        insert_hist_out=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip.output"
#    output:
#        hsmetrix=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['refseq'][:-4])+"_HsMetrics.txt",
#        per_target=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['refseq'][:-4])+"_PER_TARGET_COVERAGE.txt",
#        per_base=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['refseq'][:-4])+"_PER_BASE_COVERAGE.txt",
#    params:
#        near_distance=near_distance,
#    shell:
#        """
#        java -Xmx80G -jar {config[path_picard]} CollectHsMetrics \\
#        I={input.clipped_bam} \\
#        R={config[reference_fasta]} \\
#        O={output.hsmetrix} \\
#        BAIT_INTERVALS= {config[refseq_intervals]} \\
#        TARGET_INTERVALS= {config[refseq_intervals]} \\
#        PER_TARGET_COVERAGE= {output.per_target} \\
#        PER_BASE_COVERAGE= {output.per_base} \\
#        VALIDATION_STRINGENCY=SILENT \\
#        NEAR_DISTANCE={params.near_distance}
#        """
#
#rule hsmetrics_refseq_slop20:
#    input:
#        clipped_bam=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip.bam",
#        insert_hist_out=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip.output"
#    output:
#        hsmetrix=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['refseq_slop20'][:-4])+"_HsMetrics.txt",
#        per_target=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['refseq_slop20'][:-4])+"_PER_TARGET_COVERAGE.txt",
#        per_base=config['results_folder']+"/{sample}/{sample}_alignment.dedup.recal.clip_"+os.path.basename(config['refseq_slop20'][:-4])+"_PER_BASE_COVERAGE.txt",
#    params:
#        near_distance=near_distance,
#    shell:
#        """
#        java -Xmx80G -jar {config[path_picard]} CollectHsMetrics \\
#        I={input.clipped_bam} \\
#        R={config[reference_fasta]} \\
#        O={output.hsmetrix} \\
#        BAIT_INTERVALS= {config[refseq_intervals_slop20]} \\
#        TARGET_INTERVALS= {config[refseq_intervals_slop20]} \\
#        PER_TARGET_COVERAGE= {output.per_target} \\
#        PER_BASE_COVERAGE= {output.per_base} \\
#        VALIDATION_STRINGENCY=SILENT \\
#        NEAR_DISTANCE={params.near_distance}
#        """
