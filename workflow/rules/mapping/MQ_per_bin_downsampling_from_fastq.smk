#Snakefile

rule MQ_extraction_ds_millions:
    input: 
        recalibrated_bam=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam",
        recalibrated_bai=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam.bai",
    output:
        /merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/
