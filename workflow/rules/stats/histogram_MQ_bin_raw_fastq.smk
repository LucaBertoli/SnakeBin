#Snakefile

rule MQ_histogram:
    input:
        sorted_split_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        sorted_split_bam_bai=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai"
    output:
        MQ_histogram_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_histogram.tsv",
        mapq_text=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/mapq_values.txt")
    params:
        targets=lambda wildcards: config["target_regions"]["design"]["bed"]
    shell:
        """
        echo "Starting to create the histogram"
        samtools view {input.sorted_split_bam} | awk '{{print $5}}' > {output.mapq_text}
        sort -n {output.mapq_text} | uniq -c | awk '{{print $2 "\\t" $1}}' > {output.MQ_histogram_tsv}
        echo "Histogram created"
        """

rule MQ_histogram_target:
    input:
        sorted_split_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        sorted_split_bam_bai=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai"
    output:
        MQ_histogram_target_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_histogram_target.tsv",
        filtered_target=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/filtered_target.bam"),
        mapq_text_target=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/mapq_values_target.txt")
    params:
        targets=lambda wildcards: config["target_regions"]["design"]["bed"]
    shell:
        """
        echo "Starting to create the histogram"
        bedtools intersect -abam {input.sorted_split_bam} -b {params.targets} > {output.filtered_target}
        samtools view {output.filtered_target} | awk '{{print $5}}' > {output.mapq_text_target}
        sort -n {output.mapq_text_target} | uniq -c | awk '{{print $2 "\\t" $1}}' > {output.MQ_histogram_target_tsv}
        echo "Histogram created"
        """

rule MQ_histogram_LM:
    input:
        sorted_split_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        sorted_split_bam_bai=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai"
    output:
        MQ_histogram_LM_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_histogram_LM.tsv",
        filtered_LM=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/filtered_LM.bam"),
        mapq_text_LM=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/mapq_values_LM.txt")
    params:
        targets=lambda wildcards: config["low_mappability_regions"]["bed"]
    shell:
        """
        echo "Starting to create the histogram"
        bedtools intersect -abam {input.sorted_split_bam} -b {params.targets} > {output.filtered_LM}
        samtools view {output.filtered_LM} | awk '{{print $5}}' > {output.mapq_text_LM}
        sort -n {output.mapq_text_LM} | uniq -c | awk '{{print $2 "\\t" $1}}' > {output.MQ_histogram_LM_tsv}
        echo "Histogram created"
        """

rule MQ_stats:
    input:
        MQ_histogram_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_histogram.tsv"
    output:
        MQ_stats=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_stats.tsv"
    run:
        import pandas as pd
        import sys
        from statistics import median, mode
        df = pd.read_csv(input.MQ_histogram_tsv, sep='\t', header=None, names=['MQ', 'count'])
        expanded = df.loc[df.index.repeat(df['count'])]['MQ'].reset_index(drop=True)
        media = round(expanded.mean(), 2)
        median = round(median(expanded), 2)
        mode = round(mode(expanded), 2)
        with open(output.MQ_stats, 'w') as f:
            f.write("Statistica\tValore\n")
            f.write(f"Media\t{media}\n")
            f.write(f"Mediana\t{median}\n")
            f.write(f"Moda\t{mode}\n")

rule MQ_stats_target:
    input:
        MQ_histogram_target_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_histogram_target.tsv"
    output:
        MQ_stats_target=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_target_stats.tsv"
    run:
        import pandas as pd
        import sys
        from statistics import median, mode
        df = pd.read_csv(input.MQ_histogram_target_tsv, sep='\t', header=None, names=['MQ', 'count'])
        expanded = df.loc[df.index.repeat(df['count'])]['MQ'].reset_index(drop=True)
        media_target = round(expanded.mean(), 2)
        median_target = round(median(expanded), 2)
        mode_target = round(mode(expanded), 2)
        with open(output.MQ_stats_target, 'w') as f:
            f.write("Statistica\tValore\n")
            f.write(f"Media\t{media_target}\n")
            f.write(f"Mediana\t{median_target}\n")
            f.write(f"Moda\t{mode_target}\n")

rule MQ_stats_LM:
    input:
        MQ_histogram_LM_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_histogram_LM.tsv"
    output:
        MQ_stats_LM=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/MQ_LM_stats.tsv"
    run:
        import pandas as pd
        import sys
        from statistics import median, mode
        df = pd.read_csv(input.MQ_histogram_LM_tsv, sep='\t', header=None, names=['MQ', 'count'])
        expanded = df.loc[df.index.repeat(df['count'])]['MQ'].reset_index(drop=True)
        media_LM = round(expanded.mean(), 2)
        median_LM = round(median(expanded), 2)
        mode_LM = round(mode(expanded), 2)
        with open(output.MQ_stats_LM, 'w') as f:
            f.write("Statistica\tValore\n")
            f.write(f"Media\t{media_LM}\n")
            f.write(f"Mediana\t{median_LM}\n")
            f.write(f"Moda\t{mode_LM}\n")