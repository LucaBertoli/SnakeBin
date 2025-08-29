Pipeline in snakemake generata per comparare le performance di frammenti in base all'insert size.
Gli insert size vengono comparati raggruppando campioni WES e dividendo i relativi frammenti in bin di inserto da 1 a 1000 cons step di 50.

Al momento la pipeline esegue i seguenti moduli:
	-unione dei campioni, divisione in bin e calcolo delle statistiche sia a copertura grezza che sottocampionata in base alla copertura (con sambamba, soglie impostabili dal config.yaml)
	-divisione dei campioni in bin, downsampling a parità di frammenti presi sequenzialmente dai fastq (prendendo un numero uniforme di frammenti per ogni campione, soglie impostabili dal config.yaml), unione dei bin sottocampionati, calcolo statistiche.
	-divisione dei campioni in bin, downsampling all'1% di frammenti sequenziali partendo dai fastq mergiati (super fastq creati unendo i fastq dei singoli campioni, prendendo un numero uniforme di frammenti per ogni campione), unione dei bin sottocampionati, calcolo statistiche.

conda activate /home/morlandi/miniconda3/envs/snakemake
snakemake --cores 1 --dry-run ##test
nohup snakemake --cores 10 --use-conda & ##run
snakemake --delete-temp-output #cleanup temp files
