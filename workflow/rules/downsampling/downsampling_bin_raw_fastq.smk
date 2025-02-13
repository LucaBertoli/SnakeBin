# Snakefile
import gzip
import os
import subprocess


rule downsampling_raw_bam_per_bin:
    input:
        split_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        split_bai=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai",
        split_bam_cov_tsv=config['results_folder']+"/merged_raw_fastq/split_bam_merged/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10-stats.tsv",
    output:
        downsampled_bam=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam",
        downsampled_bai=config['results_folder']+"/merged_raw_fastq/split_bam_merged/{coverage_downsampling_levels}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip.bam.bai"
    threads: 
        config['bwa_threads']
    params:
        levels=lambda wildcards: config["coverage_downsampling_levels"][wildcards.coverage_downsampling_levels]
    run:
        try:
            with open(input.split_bam_cov_tsv, 'r') as tsv_file:
                real_coverage = float(tsv_file.readline().strip().split('\t')[3])
        except (IndexError, ValueError, FileNotFoundError) as e:
            print(f"Errore nella lettura del TSV: {e}")
            return
        target_cov = int(params.levels)
        print(target_cov)
        five_perc = target_cov * 0.05 # 5% di target_cov, equivale a 3X su 60X
        print(f"real: {real_coverage}")
        print(f"target: {target_cov}")
        if real_coverage == 0 or real_coverage < target_cov - five_perc or target_cov + five_perc >= real_coverage >= target_cov - five_perc:
            # Crea link simbolici se la copertura reale è sufficientemente vicina
            os.symlink(os.path.abspath(input.split_bam), output.downsampled_bam)
            os.symlink(os.path.abspath(input.split_bai), output.downsampled_bai)
            return
        ratio = round(target_cov / real_coverage, 3)
        print(f"ratio: {ratio}")
        # Esegui downsampling
        subprocess.run([
            config["sambamba"], 
            'view', 
            '-h', 
            '-t', str(config["bwa_threads"]), 
            '-s', str(ratio), 
            '-f', 'bam', 
            input.split_bam,  
            '-o', output.downsampled_bam  
        ], check=True)

