# Snakefile

# Fragments per samples for downsampling
million_fragment_downsampling = config['million_fragment_downsampling']

def num_farg(million_fragment_downsampling, samples):
    # Aggiunto controllo per evitare divisione per zero
    if len(samples) == 0:
        raise ValueError("La lista dei campioni è vuota, impossibile calcolare i frammenti per campione.")
    return int(million_fragment_downsampling // len(samples))

num_fragments_per_sample = num_farg(million_fragment_downsampling, samples)

rule extract_reads_from_fastq:
    input:
        bam_to_fastq_R1=config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R1.fastq.gz",
        bam_to_fastq_R2=config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R2.fastq.gz",
    output:
        downsampled_fastq_R1=temp(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{bin_start}-{bin_end}/{sample}_R1.fastq.gz"),
        downsampled_fastq_R2=temp(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{bin_start}-{bin_end}/{sample}_R2.fastq.gz"),
    threads: 
        config['bgzip_threads']
    params:
        num_frag=num_fragments_per_sample * 4,
        bgzip_threads=config["bgzip_threads"]
    shell:
        """
        set -euxo pipefail
        mkdir -p $(dirname {output.downsampled_fastq_R1})
        zcat {input.bam_to_fastq_R1} | awk 'NR <= {params.num_frag}' | bgzip -@ {params.bgzip_threads} > {output.downsampled_fastq_R1}
        zcat {input.bam_to_fastq_R2} | awk 'NR <= {params.num_frag}' | bgzip -@ {params.bgzip_threads} > {output.downsampled_fastq_R2}
        """

