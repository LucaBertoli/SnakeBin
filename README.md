# SnakeBin
<img width="2000" height="758" alt="immagine" src="https://github.com/user-attachments/assets/a98b03c8-09f3-42b3-b343-6bfae4947a01" />

Repository developed as part of my PhD thesis. 
SnakeBin is a Snakemake pipeline for insert-size-resolved Whole Exome Sequencing (WES) analysis of NGS data, stratifying reads into normalized Insert Size Bins (ISBs) to quantify fragment length effects on coverage efficiency, mappability, genotypability, and variant calling.

Starting from raw fastq files, it aligns them against the reference genome and extracts per-read insert sizes (TLEN), bins fragments (e.g., 50bp intervals), and computes ISB-stratified metrics including depth/bredth/uniformity of coverage, on/near/off-target fractions, mapping quality histograms, duplication rates, callable bases (as defined by GATK CallableLoci).

- Coverage-Based Downsampling 
  - Merges sample's fastq, trims and aligns reads, generates ISBs, subsamples reads per ISB to fixed target depth across bins (e.g., 50x average coverage, eliminating depth confounding effects), postprocesses BAMs, computes relevant metrics and optional variant calling.

- Fragment-Based Downsampling 
  - Aligns samples, generates sample-specific ISBs, converts BAMs to FASTQs, sequentially subsamples ISBs by selecting the first n° fragments (ensuring an equal contribution of each sample), merges sample-specific ISBs to cohort-specific ISBs, realigns the ISBs, postprocesses BAMs, computes relevant metrics and optional variant calling.
​

---

## 🚀 Features
- Insert-size resolved WES Analysis
- Fragment-based normalization
- Coverage-based normalization
- Variant calling
- Scalability
---

## 📁 Repository Structure

SnakeBin \
├── config # config YAML file \
├── workflow \
│ ├── SnakeFile # workflow \
│ ├── rules \
│ │ ├── downsampling # snakefiles with the downsampling rules \
│ │ ├── fastqc # snakefile with fastqc rules \
│ │ ├── mapping #snakefiles with mapping, base recalibration, clipping, duplicate removal and flagstat rules \
│ │ ├── merging # snakefiles with the fastq merging rules \
│ │ ├── positionbaseddownsampling # snakefiles with the picard position based 1% fragment downsampling with the relative duplicate removal and flagstat rules \
│ │ ├── positionbaseddownsampling_per_sample # like above, but per-sample \
│ │ ├── snv_calling # snakefiles with the vairant calling and filtering rules \
│ │ ├── split_bin # snakefiles with the insert size bin BAM generation rules \
│ │ ├── stats # snakefiles with the metrics computation rules \
│ │ ├── trimming # snakefiles with the trimming rules \
│ ├── scripts # additional python scripts for the metrics collection \
└── README.md


---

## 🛠️ Installation and Usage

```bash
#Clone locally the repository: 
git clone https://github.com/Lab-Delledonne-bioinfo/SnakeBin
cd SnakeBin

▶️ Usage
#If not installed globally, activate the snakemake conda environment:
conda activate snakemake

#Example of pipeline testing in dry-run mode:
snakemake --cores 1 --dry-run

#Example of pipeline running in nohup:
nohup snakemake --cores 10 > Example_Snakemake.log &

#Temporary data cleanup after pipeline computation:
snakemake --delete-temp-output #cleanup temp files
```
⚙️ Configuration

#Before running the pipeline edit following file. In particular, edit the sample/result folder path, insert size binning scheme, downsampling levels, target regions files for metrics computation, tools executable files.
```
config.yaml
```

📦 Dependencies
```bash
#Orchestration
Snakemake 
#main language
Python 3, java 1.8, java (latest)
#Libraries (python)
os, sys, collections, gzip, glob, csv, pandas, datetime, re, subprocess, shutil
#Tools
fastqc, fastp, BWA-MEM2, samtools, picard, bamutil, bedtools, GATK v3.8, GATK (latest), bcftools, RTG tools, sambamba
```

👤 Author \
Email: lucabertoli10@yahoo.it








