#Snakemake


rule downsampling_1_position_based:
    input:
        split_bam = config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/start_sorted.bam",
        split_bai = config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/start_sorted.bam.bai"
    output:
        split_ds_bam = temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/start_sorted.bam"),
        split_ds_bai = temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/ds_1_percent_position_based/insert_{bin_start}-{bin_end}/start_sorted.bai"),
    params:
        fraction=0.1
    shell:
        """
        gatk PositionBasedDownsampleSam -I {input.split_bam} -O {output.split_ds_bam} --FRACTION {params.fraction}
        """