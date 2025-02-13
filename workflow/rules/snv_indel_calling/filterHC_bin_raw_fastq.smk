# Snakefile
import pandas as pd


rule genotypeGVCFs:
    input:
        gcvf_file = config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.g.vcf.gz",
    output:
        raw_HC_vcf = temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw.vcf.gz"),
        raw_HC_vcf_tbi = temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw.vcf.gz.tbi"),
        raw_HC_split_vcf = temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.split.raw.vcf.gz"),
        raw_HC_split_vcf_tbi = temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.split.raw.vcf.gz.tbi")
    shell: 
        """
        gatk GenotypeGVCFs \\
        -O {output.raw_HC_vcf} \\
        -R {config[reference_fasta]} \\
        -V {input.gcvf_file} \\
        && \\
        gatk IndexFeatureFile \\
        -I {output.raw_HC_vcf} \\
        && \\
        bcftools norm \\
        {output.raw_HC_vcf} \\
        -f {config[reference_fasta]} \\
        -m - both | bgzip > {output.raw_HC_split_vcf} \\
        && \\
        gatk IndexFeatureFile \\
        -I {output.raw_HC_split_vcf}
        """

rule divide_SNP_INDEL:
    input:
        raw_HC_vcf=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.split.raw.vcf.gz",
        raw_HC_vcf_tbi=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.split.raw.vcf.gz.tbi",
    output:
        raw_snps_vcf=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_snps.vcf.gz"),
        raw_snps_vcf_tbi=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_snps.vcf.gz.tbi"),
        raw_indels_vcf=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_indels.vcf.gz"),
        raw_indels_vcf_tbi=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_indels.vcf.gz.tbi"),
    shell:
        """
        gatk SelectVariants \\
        --select-type-to-include SNP \\
        --output {output.raw_snps_vcf} \\
        -V {input.raw_HC_vcf} \\
        && \\
        gatk SelectVariants \\
        --select-type-to-exclude SNP \\
        --output {output.raw_indels_vcf} \\
        -V {input.raw_HC_vcf}
        """

rule applyGATKHardFilters:
    input:
        raw_snps_vcf=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_snps.vcf.gz",
        raw_snps_vcf_tbi=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_snps.vcf.gz.tbi",
        raw_indels_vcf=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_indels.vcf.gz",
        raw_indels_vcf_tbi=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_indels.vcf.gz.tbi",
    output:
        raw_filtered_snps_vcf=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_filtered_snps.vcf.gz"),
        raw_filtered_snps_vcf_tbi=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_filtered_snps.vcf.gz.tbi"),
        raw_filtered_indels_vcf=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_filtered_indels.vcf.gz"),
        raw_filtered_indels_vcf_tbi=temp(config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.raw_filtered_indels.vcf.gz.tbi"),
        variants_filtered_vcf=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.variants.filtered.vcf.gz",
        variants_filtered_vcf_tbi=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.variants.filtered.vcf.gz.tbi",
        variants_selected_vcf=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.variants.selected.vcf.gz",
        variants_selected_vcf_tbi=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.variants.selected.vcf.gz.tbi",
    log:
        config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.GATKfilter.log"
    shell: 
        """
        gatk VariantFiltration \\
        -O {output.raw_filtered_snps_vcf} \\
        -R {config[reference_fasta]} \\
        -V {input.raw_snps_vcf} \\
        --filter-expression "QD < 2.0" \\
        --filter-name "Broad_QD_filter" \\
        --filter-expression "MQ < 40.0" \\
        --filter-name "Broad_MQ_filter" \\
        --filter-expression "FS > 60.0" \\
        --filter-name "Broad_FS_filter" \\
        --filter-expression "SOR > 3.0" \\
        --filter-name "Broad_SOR_filter" \\
        --filter-expression "MQRankSum < -12.5" \\
        --filter-name "Broad_MQRankSum_filter" \\
        --filter-expression "ReadPosRankSum < -8.0" \\
        --filter-name "Broad_ReadPosRankSum_filter" \\
        && \\
        gatk VariantFiltration \\
        -O {output.raw_filtered_indels_vcf} \\
        -R {config[reference_fasta]} \\
        -V {input.raw_indels_vcf} \\
        --filter-expression "QD < 2.0" \\
        --filter-name "Broad_QD_filter" \\
        --filter-expression "FS > 200.0" \\
        --filter-name "Broad_FS_filter" \\
        -filter-expression "ReadPosRankSum < -20.0" \\
        --filter-name "Broad_ReadPosRankSum_filter" \\
        && \\
        gatk MergeVcfs \\
        -I {output.raw_filtered_snps_vcf} \\
        -I {output.raw_filtered_indels_vcf} \\
        -O {output.variants_filtered_vcf} \\
        && \\
        gatk SelectVariants \\
        -R {config[reference_fasta]} \\
        --variant {output.variants_filtered_vcf} \\
        --exclude-filtered \\
        -O {output.variants_selected_vcf}
        """
