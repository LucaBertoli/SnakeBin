# Snakefile
#in the fragment-based downsampling workflor, downsamples each insert size bin of each sample to a specific fragment number, ensuring that the final downsampled dataset contains a consistent number of fragments across all samples

rule extract_reads_from_fastq:
    input:
        bam_to_fastq_R1=config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R1.fastq.gz",
        bam_to_fastq_R2=config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R2.fastq.gz",
    output:
        downsampled_fastq_R1=temp(config['results_folder']+"/{sample}/downsampled_fastq/{million_fragment_downsampling}/fastq_{bin_start}-{bin_end}/{sample}_R1.fastq.gz"),
        downsampled_fastq_R2=temp(config['results_folder']+"/{sample}/downsampled_fastq/{million_fragment_downsampling}/fastq_{bin_start}-{bin_end}/{sample}_R2.fastq.gz"),
    threads: 
        config['bgzip_threads']
    params:
        num_frag = lambda wildcards: (config["million_fragment_downsampling"].get(str(wildcards.million_fragment_downsampling), 0) // len(samples))*4,
        num_samples=len(samples),
        bgzip_threads=config["bgzip_threads"]
    shell:
        """
        echo "Extracting {params.num_frag} fragments from {input.bam_to_fastq_R1} and {input.bam_to_fastq_R2} to {output.downsampled_fastq_R1} and {output.downsampled_fastq_R2}"
        set -euxo pipefail
        mkdir -p $(dirname {output.downsampled_fastq_R1})
        zcat {input.bam_to_fastq_R1} | awk 'NR <= {params.num_frag}' | bgzip -@ {params.bgzip_threads} > {output.downsampled_fastq_R1}
        zcat {input.bam_to_fastq_R2} | awk 'NR <= {params.num_frag}' | bgzip -@ {params.bgzip_threads} > {output.downsampled_fastq_R2}
        """

