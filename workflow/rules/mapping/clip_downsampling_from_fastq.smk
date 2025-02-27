# Snakefile

rule clip_bam_downsampling_fastq:
    input:
        recalibrated_bam=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam",
        recalibrated_bai=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.bam.bai",
    output:
        clipped_bam=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
    shell:
        """
        {config[path_bamutil]} clipOverlap \\
        --in {input.recalibrated_bam} \\
        --out {output.clipped_bam}
        """

rule samtools_index_downsampling_fastq:
    input:
        clipped_bam=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
    output:
        clipped_bam_index=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai",
    params:
         samtools_threads=config['samtools_threads']
    shell:
        """
        {config[path_samtools]} index \\
        --threads {params.samtools_threads} \\
        -o {output.clipped_bam_index} \\
        {input.clipped_bam}
        """

