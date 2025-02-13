# Snakefile

rule bwa_mem2_mem_merged_fastq:
    input:
        trimmed1=config['results_folder']+"/merged_raw_fastq/trimmed.R1.fastq.gz", #trimmed R1
        trimmed2=config['results_folder']+"/merged_raw_fastq/trimmed.R2.fastq.gz" #trimmed R2     
    output:
        start_sorted_bam=temp(config['results_folder']+"/merged_raw_fastq/start_sorted.bam"), #bwa-mapped samples
    threads: 
        max(config['bwa_threads'], config['samtools_threads'])
    params:
        bwa_threads=config['bwa_threads'],
        samtools_threads=config['samtools_threads'],
    shell:
        """
        {config[path_bwa]} mem \\
        -R '@RG\\tID:merged\\tPU:lane\\tLB:merged\\tSM:merged\\tCN:CGF-ddlab\\tPL:ILLUMINA' \\
        -t {params.bwa_threads} \\
        {config[reference_fasta]} \\
        {input.trimmed1} \\
        {input.trimmed2} \\
        | {config[path_samtools]} sort \\
        --threads {params.samtools_threads} \\
        -o {output.start_sorted_bam}
        """

rule samtools_index_merged_fastq:
    input:
        start_sorted_bam=config['results_folder']+"/merged_raw_fastq/start_sorted.bam"
    output:
        start_sorted_bam_bai=temp(config['results_folder']+"/merged_raw_fastq/start_sorted.bam.bai"),
    threads: 
        config['samtools_threads']
    params:
        samtools_threads=config['samtools_threads']
    shell: 
        """
        {config[path_samtools]} index \\
        --threads {params.samtools_threads} \\
        -o {output.start_sorted_bam_bai} \\
        {input.start_sorted_bam}
        """
