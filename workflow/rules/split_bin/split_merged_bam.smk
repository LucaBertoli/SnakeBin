# Snakefile

rule split_bam:
    input:
        start_sorted_bam=config['results_folder']+"/merged_raw_fastq/start_sorted.bam",
        start_sorted_bam_bai=config['results_folder']+"/merged_raw_fastq/start_sorted.bam.bai",
    output:
        split_bam=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/start_sorted.bam"),
        split_bai=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/start_sorted.bam.bai")
    shell:
        """
        {config[path_samtools]} view -h {input.start_sorted_bam} | \\
        awk 'substr($0,1,1)=="@" || ($9>={wildcards.bin_start} && $9<={wildcards.bin_end}) || ($9<=-{wildcards.bin_start} && $9>=-{wildcards.bin_end})' | \\
        {config[path_samtools]} view -b > {output.split_bam}
        {config[path_samtools]} index {output.split_bam}
        """

