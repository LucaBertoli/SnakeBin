conda activate /home/morlandi/anaconda3/envs/snakemake
snakemake --cores 1 --dry-run ##test
nohup snakemake --cores 10 --use-conda & ##run
