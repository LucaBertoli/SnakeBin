import gzip
import sys
from collections import defaultdict
import os

HELP = """SYNTAX:

 - python3 script.py <coverageBed.bed.gz> <callable1.bed> ... [callable_n.bed]

NOTE:
The destination file will be named: callable_n-stats.tsv
"""

## deifnizione funzione di help
def print_help_and_exit():
    print(HELP, file=sys.stderr)
    sys.exit()

#definizione funzione per aggiunta di una regione
def add_region(ss, chr_name, gene_name, region, length, coverage, thresholds, index, update_length):
    #se il nome del cromosoma non è presente in ss, crealo
    if chr_name not in ss:
        ss[chr_name] = {}
    #se il nome del gene non è presente in ss[chr_name], crealo
    if gene_name not in ss[chr_name]:
        ss[chr_name][gene_name] = {}
    #se la regione non è presente in ss, inizializzala inserendo 0 in tutti i campi (7 campi + un campo per ogni file callable in input)
    if region not in ss[chr_name][gene_name]:
        ss[chr_name][gene_name][region] = [0] * (7 + len(callable))

    #aggiorna i primi due elementi della lista sommando la lunghezza e la copertura della regione
    if update_length:
        ss[chr_name][gene_name][region][0] += length
        ss[chr_name][gene_name][region][1] += coverage * length

    #aggiorna il numero di basi con almeno un certo coverage (1,5,10,20,30) o basi genotipizzate, aggiungendo al contatore la lunghezza della regione se il suo coverage è maggiore o uguale alla soglia. Per aggiornare la genotipizzabilità, la soglia è 1 e la copertura sarà 1 se è callable o 0 se non è callable
    for i, threshold in enumerate(thresholds):
        if coverage >= threshold:
            ss[chr_name][gene_name][region][index + i] += length

#definizione per la stampa della regione nel file in output
def print_region(ss, chr_name, gene_name, coord):
    #genera il nome delle regioni, che gene o gene::posizione
    key = gene_name
    if gene_name != coord:
        coord_to_add= str(int(coord) - 10**10)
        key += f"::{coord_to_add}"

    region_data = ss[chr_name][gene_name][coord]
    #calcola la percentuale della regione coperta ad una serta soglia
    row = [
        chr_name, key, region_data[0],
        f"{region_data[1] / region_data[0]:.2f}",
        *(f"{100 * region_data[i] / region_data[0]:.2f}" for i in range(2, 7))
    ]

    #calcola la percentuale della regione genotipizzata per ogni callable
    for i in range(len(callable)):
        row.append(f"{100 * region_data[7 + i] / region_data[0]:.2f}")

    #mette tutti i dati della regione in una riga separata da tab
    return "\t".join(map(str, row))

#inizio main script
if len(sys.argv) < 3:
    print_help_and_exit()

#definizione nome output file
dest = sys.argv[-1].replace(".bed", "-stats.tsv")

#se il file è gia presente e non è vuoto, non fare nulla
if os.path.isfile(dest) and os.path.getsize(dest) > 0:
    print("Nothing to do", file=sys.stderr)
    sys.exit()
else:
    print(f"Producing stats file: {dest}", file=sys.stderr)

#definizione variabili
depth = ""
callable = []
ss = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

#lettura dei file di input, il primo file sarà il file di copertura, il secondo e terzo saranno callable
for f in sys.argv[1:]:
    if os.path.isfile(f):
        if not depth:
            depth = f
        else:
            callable.append(f)

#apertura del file di copertura (bedtools coverage del BAM sul design) e lettura delle righe
with gzip.open(depth, 'rt') as file:
    for line in file:
        if line.startswith("all") or line.startswith("genome"):
            continue

        fields = line.strip().split("\t")
        #esempio fields = ['chr1', '1693320', '1693482', 'NADK', '84', '1', '162', '0.0061728']
        add_region(ss, fields[0], fields[3], str(int(fields[1]) + 10**10), int(fields[5]), int(fields[4]), [1, 5, 10, 20, 30], 2, 1)
        add_region(ss, fields[0], fields[3], fields[3], int(fields[5]), int(fields[4]), [1, 5, 10, 20, 30], 2, 1)
        add_region(ss, 'all', 'all', 'all', int(fields[5]), int(fields[4]), [1, 5, 10, 20, 30], 2, 1)

#aggiorna le genotipizzabilità delle regioni usando i le coperture dei callable (generati da bedtools coverage di callable.bed e callable_DP10.bed sul design) senza aggiornare le lunghezze e coverage delle regioni
for i, c_file in enumerate(callable):
    with open(c_file, 'r') as file:
        for line in file:
            if line.startswith("all") or line.startswith("genome"):
                continue

            fields = line.strip().split("\t")
            if int(fields[-4]) != 1:
                continue

            add_region(ss, fields[0], fields[3], str(int(fields[1]) + 10**10), int(fields[5]), int(fields[4]), [1], 7 + i, 0)
            add_region(ss, fields[0], fields[3], fields[3], int(fields[5]), int(fields[4]), [1], 7 + i, 0)
            add_region(ss, 'all', 'all', 'all', int(fields[5]), int(fields[4]), [1], 7 + i, 0)

#apertura e scrittura del file tsv in output
with open(dest, 'w') as out:
    out.write(print_region(ss, 'all', 'all', 'all') + "\n")

    for chr_name in sorted(ss):
        #stampa i geni
        for gene_name in sorted(ss[chr_name]):
            if gene_name == 'all':
                continue
            out.write(print_region(ss, chr_name, gene_name, gene_name) + "\n")

        #stampa le regioni
        for gene_name in sorted(ss[chr_name]):
            for pos in sorted(ss[chr_name][gene_name]):
                if pos == gene_name:
                    continue
                out.write(print_region(ss, chr_name, gene_name, pos) + "\n")


## ss è un dizionario di un dizionario di una lista, a differenza di un dizionario normale se una key non esiste essa viene creata automaticamente. La struttura della variabile è:
## ss = {
##     'chr1': {
##         'geneA': {
##             'region1': [length, coverage, ...],
##             'region2': [length, coverage, ...]
##         },
##         'geneB': {
##             'region1': [length, coverage, ...]
##         }
##     },
##     'chr2': {
##         'geneC': {
##             'region1': [length, coverage, ...]
##         }
##     }
## }