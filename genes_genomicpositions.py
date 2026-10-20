#!/usr/bin/env python3

from collections import defaultdict

############################################################
# Input
############################################################

GFF_FILE = "genes_annotations_612.gff3"

############################################################
# Helper function
############################################################

def parse_attributes(attr):

    d = {}

    for item in attr.strip().split(";"):

        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v

    return d


############################################################
# Data structures
############################################################

gene_transcripts = defaultdict(list)

transcripts = {}

############################################################
# First pass
# Collect genes and transcripts
############################################################

print("Reading GFF...")

with open(GFF_FILE) as f:

    for line in f:

        if line.startswith("#"):
            continue

        cols = line.rstrip().split("\t")

        if len(cols) != 9:
            continue

        feature = cols[2]

        attr = parse_attributes(cols[8])

        if feature == "mRNA":

            tr = attr.get("ID")
            gene = attr.get("Parent")

            if tr is None or gene is None:
                continue

            gene_transcripts[gene].append(tr)

            transcripts[tr] = {
                "blocks": [],
                "cds": [],
                "utr5": [],
                "utr3": []
            }

############################################################
# Second pass
# Collect CDS and UTRs
############################################################

with open(GFF_FILE) as f:

    for line in f:

        if line.startswith("#"):
            continue

        cols = line.rstrip().split("\t")

        if len(cols) != 9:
            continue

        feature = cols[2]

        start = int(cols[3])
        end = int(cols[4])

        attr = parse_attributes(cols[8])

        parent = attr.get("Parent")

        if parent is None:
            continue

        for tr in parent.split(","):

            if tr not in transcripts:
                continue

            if feature == "CDS":

                transcripts[tr]["cds"].append((start, end))
                transcripts[tr]["blocks"].append((start, end))

            elif feature in (
                    "five_prime_UTR",
                    "5UTR",
                    "five_prime_utr"):

                transcripts[tr]["utr5"].append((start, end))
                transcripts[tr]["blocks"].append((start, end))

            elif feature in (
                    "three_prime_UTR",
                    "3UTR",
                    "three_prime_utr"):

                transcripts[tr]["utr3"].append((start, end))
                transcripts[tr]["blocks"].append((start, end))

############################################################
# Analyse
############################################################

cds_bp = 0
utr5_bp = 0
utr3_bp = 0
intron_bp = 0

genes_used = 0

for gene, tr_list in gene_transcripts.items():

    ########################################################
    # Longest transcript
    ########################################################

    best = None
    best_len = -1

    for tr in tr_list:

        length = 0

        for s, e in transcripts[tr]["blocks"]:
            length += e - s + 1

        if length > best_len:
            best_len = length
            best = tr

    if best is None:
        continue

    genes_used += 1

    ########################################################
    # Sum CDS / UTR lengths
    ########################################################

    for s, e in transcripts[best]["cds"]:
        cds_bp += e - s + 1

    for s, e in transcripts[best]["utr5"]:
        utr5_bp += e - s + 1

    for s, e in transcripts[best]["utr3"]:
        utr3_bp += e - s + 1

    ########################################################
    # Reconstruct exons
    ########################################################

    blocks = sorted(transcripts[best]["blocks"])

    if len(blocks) == 0:
        continue

    merged = [list(blocks[0])]

    for s, e in blocks[1:]:

        if s <= merged[-1][1] + 1:

            merged[-1][1] = max(merged[-1][1], e)

        else:

            merged.append([s, e])

    ########################################################
    # Introns = gaps between merged exons
    ########################################################

    if len(merged) > 1:

        for i in range(len(merged)-1):

            intron_start = merged[i][1] + 1
            intron_end = merged[i+1][0] - 1

            if intron_end >= intron_start:

                intron_bp += intron_end - intron_start + 1

############################################################
# Output
############################################################

total = cds_bp + intron_bp + utr5_bp + utr3_bp

print("\n-----------------------------------------------------")
print("Genomic annotation (longest transcript per gene)")
print("-----------------------------------------------------")
print(f"Genes analysed : {genes_used}")
print("")

print(f"CDS      : {cds_bp/1e6:8.2f} Mb   ({100*cds_bp/total:6.2f}%)")
print(f"Introns  : {intron_bp/1e6:8.2f} Mb   ({100*intron_bp/total:6.2f}%)")
print(f"5' UTR   : {utr5_bp/1e6:8.2f} Mb   ({100*utr5_bp/total:6.2f}%)")
print(f"3' UTR   : {utr3_bp/1e6:8.2f} Mb   ({100*utr3_bp/total:6.2f}%)")

print("-----------------------------------------------------")
print(f"Total    : {total/1e6:8.2f} Mb")
print("-----------------------------------------------------")