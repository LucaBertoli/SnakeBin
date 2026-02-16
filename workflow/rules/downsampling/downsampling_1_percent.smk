# Snakefile
#in the fragment-based downsampling workflow, downsamples each insert size bin of each sample to 1% of the original number of fragments, then merges the downsampled fastq files across samples for each bin.

rule extract_reads_from_fastq_ds_1_perc:
    input:
        bam_to_fastq_R1=config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R1.fastq.gz",
        bam_to_fastq_R2=config['results_folder']+"/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/{sample}_R2.fastq.gz",
        html_R1=config['results_folder'] + "/{sample}/split_fastq/bam_to_fastq_{bin_start}-{bin_end}/fastqc/{sample}_R1.html",
    output:
        downsampled_fastq_R1=temp(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{bin_start}-{bin_end}/{sample}_1_percent_sequencial_R1.fastq.gz"),
        downsampled_fastq_R2=temp(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{bin_start}-{bin_end}/{sample}_1_percent_sequencial_R2.fastq.gz"),
    threads: 
        config['bgzip_threads']
    params:
        bgzip_threads=config["bgzip_threads"]
    shell:
        """
        set -euxo pipefail
        mkdir -p $(dirname {output.downsampled_fastq_R1})
        max_lines=$(python3 workflow/scripts/metricsParser.py {input.html_R1} | awk -F "\\t" '{{printf "%.0f\\n", ($1/100)}}' | awk -F "\\t" '{{print ($1*4)}}')
        zcat {input.bam_to_fastq_R1} | awk -v max_lines=$max_lines 'NR <= max_lines' | bgzip -@ {params.bgzip_threads} > {output.downsampled_fastq_R1}
        zcat {input.bam_to_fastq_R2} | awk -v max_lines=$max_lines 'NR <= max_lines' | bgzip -@ {params.bgzip_threads} > {output.downsampled_fastq_R2}
        """

rule merge_reads_from_fastq_ds_1_perc:
    input:
        downsampled_1_percent_fastq_R1=expand(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{{bin_start}}-{{bin_end}}/{sample}_1_percent_sequencial_R1.fastq.gz",sample=samples),
        downsampled_1_percent_fastq_R2=expand(config['results_folder']+"/{sample}/downsampled_fastq/fastq_{{bin_start}}-{{bin_end}}/{sample}_1_percent_sequencial_R2.fastq.gz",sample=samples),
    output:
        merged_1_percent_fastq_R1=temp(config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/merged.R1.fastq.gz"),
        merged_1_percent_fastq_R2=temp(config['results_folder']+"/merged_downsampled_1_percent_sequencial/insert_{bin_start}-{bin_end}/merged.R2.fastq.gz"),
    shell:
        """
        cat {input.downsampled_1_percent_fastq_R1} > {output.merged_1_percent_fastq_R1}
        cat {input.downsampled_1_percent_fastq_R2} > {output.merged_1_percent_fastq_R2}
        """
   

