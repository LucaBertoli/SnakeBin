# Snakefile

rule callHC_downsampling_from_fastq:
    input:
        clipped_bam = config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        clipped_bam_bai = config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai",
    output:
        raw_HC_gvcf = config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.g.vcf.gz",
    log:
        config['results_folder']+"/merged_downsampled_bam/insert_{bin_start}-{bin_end}/insert_{bin_start}-{bin_end}.log"
    shell:
        """
        gatk HaplotypeCaller \\
        -R {config[reference_fasta]} \\
        -I {input.clipped_bam} \\
        -ERC GVCF \\
        --output {output.raw_HC_gvcf} \\
        --standard-min-confidence-threshold-for-calling 30.0 \\
        --dont-use-soft-clipped-bases true \\
        --sample-ploidy 2
        """
