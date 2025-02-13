import os
import sys

def help():
    print("""
SYNTAX:
 - script.py <filename>

The program retrieves metrics from files that match:

 - html == fastqc (number of frags and %GC)
 - *.output == CollectInsertSizeMetrics (average insert size)
 - dupl*.txt == MarkDuplicated (% duplicates)
 - flagstat == samtools flagstats (# mapped fragments)
 - *stats.tsv == first line of gene stats table
 - HsMetrics.txt == picard HsMetrics
   (ON_BAIT_BASES, NEAR_BAIT_BASES, OFF_BAIT_BASES,
    FOLD_ENRICHMENT, FOLD_80_BASE_PENALTY)
""", file=sys.stderr)

def get_fastqc(html):
    sample = os.path.basename(html).split('/')[0]
    reads = 0
    length = 0
    gc = 0
    with os.popen(f'lynx --dump {html}') as f:
        for line in f:
            if "Total Sequences" in line:
                reads = int(line.split()[-1])
            elif "%GC" in line:
                gc = int(line.split()[-1])
    print(f"{reads}\t{gc}")

def get_coverage(file):
    sample = os.path.basename(os.path.dirname(file))
    with open(file) as f:
        line = f.readline().strip()
        if line.startswith("all"):
            a = line.split('\t')
            a[0] = sample
            print("\t".join(a[2:]))

def get_flagstats(file):
    sample = os.path.basename(os.path.dirname(file))
    supplementary = 0
    mate_mapped = 0
    singletons = 0
    with open(file) as f:
        for line in f:
            if "with itself and mate mapped" in line:
                mate_mapped = int(int(line.split()[0]) / 2)
            elif "singletons" in line:
                singletons = int(line.split()[0])
    print(f"{supplementary + mate_mapped + singletons}")

def get_duplicates(file):
    sample = os.path.basename(os.path.dirname(file))
    with open(file) as f:
        for line in f:
            fields = line.split('\t')
            if len(fields) > 8 and fields[8].startswith("0."):
                print(f"{float(fields[8]) * 100:.2f}")
                break

def insert_size(file):
    sample = os.path.basename(os.path.dirname(file))
    index = 0
    with open(file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if "MEDIAN_INSERT_SIZE" in fields:
                index = fields.index("MEAN_INSERT_SIZE")
                dim = next(f).strip().split('\t')
                print(f"{float(dim[index]):.2f}\t{float(dim[index-5]):.2f}\t{float(dim[index-4]):.2f}")
                break

def get_hs_metrics(file):
    sample = os.path.basename(os.path.dirname(file))
    with open(file) as f:
        for line in f:
            if "ON_BAIT_BASES" in line:
                headers = line.strip().split('\t')
                pos = {h: i for i, h in enumerate(headers)}
                values = next(f).strip().split('\t')
                on_bait = float(values[pos['ON_BAIT_BASES']])
                near_bait = float(values[pos['NEAR_BAIT_BASES']])
                off_bait = float(values[pos['OFF_BAIT_BASES']])
                total = on_bait + near_bait + off_bait
                print(f"{on_bait/total*100:.2f}\t{near_bait/total*100:.2f}\t{off_bait/total*100:.2f}\t"
                      f"{float(values[pos['AT_DROPOUT']]):.2f}\t{float(values[pos['GC_DROPOUT']]):.2f}\t"
                      f"{values[pos['FOLD_80_BASE_PENALTY']]}")
                break

def main():
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
        file = sys.argv[1]
        if file.endswith('html'):
            get_fastqc(file)
        elif file.endswith('duplicates.txt'):
            get_duplicates(file)
        elif file.endswith('flagstat.txt') or file.endswith('flagstat_recall.txt'):
            get_flagstats(file)
        elif file.endswith('output'):
            insert_size(file)
        elif file.endswith('HsMetrics.txt'):
            get_hs_metrics(file)
        elif file.endswith('stats.tsv'):
            get_coverage(file)
        else:
            help()
    else:
        help()

if __name__ == "__main__":
    main()

