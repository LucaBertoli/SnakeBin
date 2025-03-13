#Snakemake


rule downsampling_1_position_based_per_sample:
    input:
        split_bam = config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/start_sorted.bam",
        split_bai = config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/start_sorted.bam.bai"
    output:
        split_ds_bam = temp(config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/start_sorted.bam"),
        split_ds_bai = temp(config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/start_sorted.bam.bai"),
    params:
        fraction=0.01
    shell:
        """
        gatk PositionBasedDownsampleSam -I {input.split_bam} -O {output.split_ds_bam} --FRACTION {params.fraction}
        samtools index {output.split_ds_bam}
        """

rule sort_split_bam_samples_1_position_based_per_sample:
    input:
        split_bam=config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/start_sorted.bam",
        split_bai=config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/start_sorted.bam.bai"
    output:
        sorted_split_bam=temp(config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/sorted_start_sorted.bam"),
    shell:
        """
        {config[path_samtools]} sort -T $(dirname {input.split_bam}) -n {input.split_bam} -o {output.sorted_split_bam}
        """


rule sort_split_bam_samples_to_fastq_1_position_based_per_sample:
    input:
        sorted_split_bam=config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/sorted_start_sorted.bam",
    output:
        bam_to_fastq_R1=temp(config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/{sample}_R1.fastq.gz"),
        bam_to_fastq_R2=temp(config['results_folder']+"/{sample}/split_bam/ds_1_percent_position_based_per_sample/insert_{bin_start}-{bin_end}/{sample}_R2.fastq.gz")
    shell:
        """
        {config[path_samtools]} bam2fq -1 {output.bam_to_fastq_R1} -2 {output.bam_to_fastq_R2} {input.sorted_split_bam}
        """