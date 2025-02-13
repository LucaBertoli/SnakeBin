# Snakefile

rule callable_loci_downsampling_from_fastq:
    input:
        clipped_bam=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        clipped_bam_index=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai"
    output:
        callable_table=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_table.txt",
        callable_status=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_status.bed",
    shell:
        """
        {config[path_java1.8]} -jar {config[path_gatk3.8]} \\
        -T CallableLoci \\
        -R {config[reference_fasta]} \\
        -I {input.clipped_bam} \\
        -summary {output.callable_table} \\
        -o {output.callable_status}
        """

rule callable_bed_downsampling_from_fastq:
    input:
        callable_status=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_status.bed"
    output:
        callable_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable.bed"
    shell:
        """
        awk -F "\\t" '$4=="CALLABLE" {{print $1"\\t"$2"\\t"$3}}' {input.callable_status} > {output.callable_bed}
        """



rule callable_loci_DP10_downsampling_from_fastq:
    input:
        clipped_bam=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        clipped_bam_index=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai"
    output:
        callable_table=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_table_DP10.txt",
        callable_status=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_status_DP10.bed",
    shell:
        """
        {config[path_java1.8]} -jar {config[path_gatk3.8]} \\
        -T CallableLoci \\
        -R {config[reference_fasta]} \\
        -minDepth 10 \\
        -I {input.clipped_bam} \\
        -summary {output.callable_table} \\
        -o {output.callable_status}
        """

rule callable_bed_DP10_downsampling_from_fastq:
    input:
        callable_status=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_status_DP10.bed"
    output:
        callable_bed=config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/callable_DP10.bed"
    shell:
        """
        awk -F "\\t" '$4=="CALLABLE" {{print $1"\\t"$2"\\t"$3}}' {input.callable_status} > {output.callable_bed}
        """
