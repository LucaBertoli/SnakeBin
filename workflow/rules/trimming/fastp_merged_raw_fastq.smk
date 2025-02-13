# Snakefile

rule trimming_merged_raw_fastq:
    input:
        fastq_R1=config['results_folder']+"/merged_raw_fastq/merged.R1.fastq.gz",
        fastq_R2=config['results_folder']+"/merged_raw_fastq/merged.R2.fastq.gz"
    output:
        trimmed_R1=temp(config['results_folder']+"/merged_raw_fastq/trimmed.R1.fastq.gz"),
        trimmed_R2=temp(config['results_folder']+"/merged_raw_fastq/trimmed.R2.fastq.gz")
    threads: 
        config['fastp_threads']
    params:
        threads=config['fastp_threads']
    shell:
        """
        {config[path_fastp]} \\
        -i {input.fastq_R1} \\
        -I {input.fastq_R2} \\
        -o {output.trimmed_R1} \\
        -O {output.trimmed_R2} \\
        --failed_out {config[results_folder]}/merged_raw_fastq/trimmed_failed.fastq.gz \\
        -h {config[results_folder]}/merged_raw_fastq/trimmed_adapters_removal_report.html \\
        -j {config[results_folder]}/merged_raw_fastq/trimmed_fastp.json \\
        -w {params.threads}
        """

