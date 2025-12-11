# Snakefile

rule split_bam_samples:
    input:
        start_sorted_bam=config['results_folder']+"/{sample}/{sample}_start_sorted.bam",
        start_sorted_bam_bai=config['results_folder']+"/{sample}/{sample}_start_sorted.bam.bai",
    output:
        split_bam=temp(config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/start_sorted.bam"),
        split_bai=temp(config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/start_sorted.bam.bai")
    shell:
        """
        {config[path_samtools]} view -h {input.start_sorted_bam} | \\
        awk 'substr($0,1,1)=="@" || ($9>={wildcards.bin_start} && $9<={wildcards.bin_end}) || ($9<=-{wildcards.bin_start} && $9>=-{wildcards.bin_end})' | \\
        {config[path_samtools]} view -b > {output.split_bam}
        {config[path_samtools]} index {output.split_bam}
        """

rule sort_split_bam_samples:
    input:
        split_bam=config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/start_sorted.bam",
        split_bai=config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/start_sorted.bam.bai"
    output:
        sorted_split_bam=temp(config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/sorted_start_sorted.bam"),
    shell:
        """
        {config[path_samtools]} sort -T $(dirname {input.split_bam}) -n {input.split_bam} -o {output.sorted_split_bam}
        """

rule sort_split_bam_samples_to_fastq:
    input:
        sorted_split_bam=config['results_folder']+"/{sample}/split_bam/insert_{bin_start}-{bin_end}/sorted_start_sorted.bam",
    output:
        bam_to_fastq_R1=temp(config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R1.fastq.gz"),
        bam_to_fastq_R2=temp(config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R2.fastq.gz")
    shell:
        """
        {config[path_samtools]} bam2fq -1 {output.bam_to_fastq_R1} -2 {output.bam_to_fastq_R2} {input.sorted_split_bam}
        """