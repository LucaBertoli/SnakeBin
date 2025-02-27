# Snakefile


rule unif_coverage_downsampling_from_fastq:
    input:
        tsv=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip-callable_DP10-stats.tsv",
        per_base_cov=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/alignment.dedup.recal.clip_PER_BASE_COVERAGE.txt"
    output:
        unif_cov=config['results_folder']+"/merged_downsampled_bam/{million_fragment_downsampling}/insert_{bin_start}-{bin_end}/unif_of_coverage.txt"
    params:
        targets=lambda wildcards: config["target_regions"]["design"]["bed"]
    shell:
        """
        tsv={input.tsv}; #coverage file from callable loci analysis
        target={input.per_base_cov}; #PER_BASE_COVERAGE.txt from collectHSMetrics
        kit={params.targets};
        t=$(head -1 $tsv | cut -f4);
        t2=$(echo | awk '{{print "'"$t"'"*0.2}}');
        threshold=$(printf "%.0f\\n" $t2);
        echo "Threshold: $threshold" > {output};

        total=`cat $target | wc -l | cut -d ' ' -f1`;
        echo "Total #of bases: $total" >> {output} ;

        bases=`cat $target | awk -v t=$threshold '$4>t' | wc -l`;
        echo "Total bases: $bases" >> {output};

        pct=`printf "%.2f\n" $(echo | awk "{{print ($bases/$total)*100}}")`;
        echo "Uniformity of coverage (Pct > 0.2*mean): $pct" >> {output};
        """

