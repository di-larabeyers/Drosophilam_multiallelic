#!/usr/bin/env python3

############################################################
# annotationfinder.py
#
# Annotates multiallelic loci using FlyBase GFF3.
#
# Gene assignment follows the same logic as genefinder.py:
#   BED -> bedtools intersect -> GFF feature "gene"
#   -> Name attribute
#
# Additional annotations:
#   Gene
#   Gene_ID
#   Transcript
#   Feature
#   Region
#   Protein_coding
#   Strand
#
# Important:
#   - Multiple genes at the same genomic position are preserved.
#   - Multiple transcripts are preserved.
#   - Positions without a GFF3 overlap are retained as NA.
#
# Input:
#   hcfilter.bed
#   linvilla_allc_hcf_converted.tsv
#   genes_annotations_612.gff3
#
# Output:
#   annotated_multiallelic_sites.tsv
############################################################


import subprocess
import pandas as pd
from collections import defaultdict


############################################################
# Files
############################################################

BED_FILE = "hcfilter.bed"
INPUT_FILE = "linvilla_allc_hcf_converted.tsv"
GFF_FILE = "genes_annotations_612.gff3"

TEMP_FILE = "temp_intersect2.tsv"
OUTPUT_FILE = "annotated_multiallelic_sites2.tsv"


############################################################
# 1. Run bedtools
#
# -wao:
#   Write the original BED entry for every overlap.
#   Importantly, also keep BED entries that have NO
#   overlap with the GFF3.
#
# This allows positions without annotation to be
# retained and later assigned "NA".
############################################################

cmd = (
    f'bedtools intersect '
    f'-wao '
    f'-a "{BED_FILE}" '
    f'-b "{GFF_FILE}" '
    f'> "{TEMP_FILE}"'
)

print("Running bedtools...")
print(cmd)

subprocess.run(cmd, shell=True, check=True)


############################################################
# 2. Read intersect output
############################################################

print("Reading intersect output...")

intersect = pd.read_csv(
    TEMP_FILE,
    sep="\t",
    header=None,
    dtype=str
)

print("Intersect rows:", len(intersect))


############################################################
# 3. Parse GFF
#
# This is used for:
#   - transcript relationships
#   - features
#   - protein-coding status
#
# Gene assignment itself is taken directly from the
# GFF "gene" feature in the bedtools intersect output.
############################################################

print("Parsing GFF...")

gene_info = {}
transcript_to_gene = {}


with open(GFF_FILE) as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")

        if len(fields) != 9:
            continue

        chrom = fields[0]
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes = fields[8]

        ####################################################
        # Attribute dictionary
        ####################################################

        attr = {}

        for item in attributes.split(";"):

            if "=" in item:

                k, v = item.split("=", 1)
                attr[k] = v


        ####################################################
        # Gene
        ####################################################

        if feature == "gene":

            gene_id = attr.get("ID")

            if gene_id is None:
                continue

            gene_info[gene_id] = {

                "gene_name": attr.get("Name", "NA"),

                "chrom": chrom,

                "start": start,

                "end": end,

                "strand": strand,

                "transcripts": set(),

                "features": set(),

                "protein_coding": False

            }


        ####################################################
        # mRNA
        ####################################################

        elif feature == "mRNA":

            transcript = attr.get("ID")
            parent = attr.get("Parent")

            if transcript is None:
                continue

            transcript_to_gene[transcript] = parent

            if parent in gene_info:

                gene_info[parent]["transcripts"].add(transcript)


        ####################################################
        # Transcript-dependent features
        ####################################################

        elif feature in {
            "CDS",
            "exon",
            "five_prime_UTR",
            "three_prime_UTR",
            "start_codon",
            "stop_codon"
        }:

            parent = attr.get("Parent")

            if parent is None:
                continue

            parents = parent.split(",")

            for transcript in parents:

                if transcript not in transcript_to_gene:
                    continue

                gene = transcript_to_gene[transcript]

                if gene not in gene_info:
                    continue

                gene_info[gene]["features"].add(feature)

                if feature == "CDS":

                    gene_info[gene]["protein_coding"] = True


print("--------------------------------------")
print("Genes in GFF:", len(gene_info))
print("Transcripts:", len(transcript_to_gene))

protein = sum(
    g["protein_coding"]
    for g in gene_info.values()
)

print("Protein coding genes:", protein)
print("Non-coding genes:", len(gene_info) - protein)
print("--------------------------------------")


############################################################
# 4. Collect annotations directly from BED/intersect
#
# IMPORTANT:
#
# The annotation is attached to the actual BED interval:
#
#   key = (chromosome, start, end)
#
# Multiple GFF features belonging to the same BED position
# are collected together.
#
# Multiple genes are stored in a SET and therefore are NOT
# overwritten.
############################################################

print("Collecting annotations...")


annotations = defaultdict(lambda: {

    "Gene": set(),

    "Gene_ID": set(),

    "Transcript": set(),

    "Features": set(),

    "Feature": "NA",

    "Region": "NA",

    "Protein_coding": set(),

    "Strand": set()

})


############################################################
# BED has 3 columns:
#
#   0 = chromosome
#   1 = start
#   2 = end
#
# GFF starts at column 3.
############################################################

BED_COLS = 3


############################################################
# Process every bedtools intersection
############################################################

for _, row in intersect.iterrows():

    ########################################################
    # BED coordinates
    ########################################################

    chrom = str(row[0])
    start = int(row[1])
    end = int(row[2])

    key = (chrom, start, end)

    ann = annotations[key]


    ########################################################
    # Check whether a GFF feature was found
    #
    # With -wao, a position without overlap receives "."
    # in the GFF columns.
    ########################################################

    feature = str(row[BED_COLS + 2])

    if feature == ".":

        # No GFF3 annotation for this BED position.
        # Keep the default NA values.
        continue


    ########################################################
    # GFF information
    ########################################################

    strand = str(row[BED_COLS + 6])
    attributes = str(row[BED_COLS + 8])


    ########################################################
    # Store feature
    ########################################################

    ann["Features"].add(feature)


    ########################################################
    # Gene assignment
    #
    # Only a GFF "gene" feature assigns a gene.
    #
    # IMPORTANT:
    # Multiple genes are stored using .add()
    # instead of replacing the previous gene.
    ########################################################

    if feature == "gene":

        attr = {}

        for item in attributes.split(";"):

            if "=" in item:

                k, v = item.split("=", 1)
                attr[k] = v


        gene_name = attr.get(
            "Name",
            attr.get("ID", "NA")
        )

        gene_id = attr.get(
            "ID",
            "NA"
        )


        # Store ALL genes
        ann["Gene"].add(gene_name)

        # Store ALL gene IDs
        ann["Gene_ID"].add(gene_id)

        # Store strand
        if strand != ".":
            ann["Strand"].add(strand)


    ########################################################
    # Transcript
    #
    # Multiple transcripts are preserved.
    ########################################################

    elif feature == "mRNA":

        attr = {}

        for item in attributes.split(";"):

            if "=" in item:

                k, v = item.split("=", 1)
                attr[k] = v


        transcript = attr.get("ID")

        if transcript is not None:

            ann["Transcript"].add(transcript)

        if strand != ".":
            ann["Strand"].add(strand)


############################################################
# 5. Determine feature and region
#
# This does NOT influence gene assignment.
#
# Feature priority:
#
#   CDS
#   exon
#   five_prime_UTR
#   three_prime_UTR
#   gene
#
############################################################

print("Assigning regions...")


priority = [

    "CDS",

    "exon",

    "five_prime_UTR",

    "three_prime_UTR",

    "gene"

]


region_name = {

    "CDS": "CDS",

    "exon": "Exon",

    "five_prime_UTR": "5UTR",

    "three_prime_UTR": "3UTR",

    "gene": "Intron"

}


for ann in annotations.values():

    feature = "NA"


    ########################################################
    # Select highest-priority feature
    ########################################################

    for f in priority:

        if f in ann["Features"]:

            feature = f

            break


    ann["Feature"] = feature


    ########################################################
    # Assign region
    ########################################################

    if feature == "NA":

        ann["Region"] = "NA"

    else:

        ann["Region"] = region_name[feature]


    ########################################################
    # Protein-coding status
    #
    # If at least one assigned gene is protein-coding,
    # the position is classified as protein-coding.
    #
    # If genes are assigned but none are protein-coding,
    # the position is classified as non-protein-coding.
    #
    # If no gene is assigned:
    #   NA
    ########################################################

    if len(ann["Gene_ID"]) == 0:

        ann["Protein_coding"] = {"NA"}

    else:

        protein_status = set()

        for gene_id in ann["Gene_ID"]:

            if gene_id in gene_info:

                if gene_info[gene_id]["protein_coding"]:

                    protein_status.add("Yes")

                else:

                    protein_status.add("No")


        if "Yes" in protein_status:

            ann["Protein_coding"] = {"Yes"}

        elif "No" in protein_status:

            ann["Protein_coding"] = {"No"}

        else:

            ann["Protein_coding"] = {"NA"}


print("Finished assigning regions.")


############################################################
# 6. Read original multiallelic table
############################################################

print("Reading allele numbers...")


alleles = pd.read_csv(

    INPUT_FILE,

    sep="\t",

    usecols=[
        "chr",
        "pos",
        "n_alleles"
    ]

)


############################################################
# 7. Keep one row per genomic position
#
# The input contains 15 rows per position.
############################################################

alleles = alleles.iloc[::15].reset_index(drop=True)


alleles = alleles.sort_values(

    ["chr", "pos"]

)


print("Positions:", len(alleles))


############################################################
# 8. Read BED coordinates
############################################################

print("Reading BED positions...")


bed_positions = pd.read_csv(

    BED_FILE,

    sep="\t",

    header=None,

    names=[
        "chr",
        "start",
        "end"
    ],

    dtype={
        "chr": str
    }

)


print("BED positions:", len(bed_positions))


############################################################
# 9. Build mapping from TSV position to BED interval
#
# IMPORTANT:
#
# This assumes that the TSV "pos" corresponds to the
# BED "start".
#
# No automatic pos-1 or pos+1 correction is applied.
############################################################

bed_map = {}


for _, row in bed_positions.iterrows():

    key = (

        str(row["chr"]),

        int(row["start"])

    )


    bed_map[key] = (

        str(row["chr"]),

        int(row["start"]),

        int(row["end"])

    )


############################################################
# 10. Create output
############################################################

output = []


for _, row in alleles.iterrows():

    chr_ = str(row["chr"])

    pos = int(row["pos"])


    ########################################################
    # Find corresponding BED interval
    ########################################################

    bed_key = (chr_, pos)


    if bed_key in bed_map:

        annotation_key = bed_map[bed_key]


        if annotation_key in annotations:

            ann = annotations[annotation_key]


        else:

            ann = {

                "Gene": set(),

                "Gene_ID": set(),

                "Transcript": set(),

                "Feature": "NA",

                "Region": "NA",

                "Protein_coding": {"NA"},

                "Strand": set()

            }


    else:

        ####################################################
        # No corresponding BED position
        ####################################################

        ann = {

            "Gene": set(),

            "Gene_ID": set(),

            "Transcript": set(),

            "Feature": "NA",

            "Region": "NA",

            "Protein_coding": {"NA"},

            "Strand": set()

        }


    ########################################################
    # Convert sets to strings
    #
    # Multiple values are separated by ";"
    ########################################################

    genes = (

        ";".join(sorted(ann["Gene"]))

        if len(ann["Gene"]) > 0

        else "NA"

    )


    gene_ids = (

        ";".join(sorted(ann["Gene_ID"]))

        if len(ann["Gene_ID"]) > 0

        else "NA"

    )


    transcripts = (

        ";".join(sorted(ann["Transcript"]))

        if len(ann["Transcript"]) > 0

        else "NA"

    )


    strands = (

        ";".join(sorted(ann["Strand"]))

        if len(ann["Strand"]) > 0

        else "NA"

    )


    protein_coding = (

        ";".join(sorted(ann["Protein_coding"]))

        if len(ann["Protein_coding"]) > 0

        else "NA"

    )


    ########################################################
    # Add output row
    ########################################################

    output.append({

        "chr": row["chr"],

        "pos": row["pos"],

        "n_alleles": row["n_alleles"],

        "Gene": genes,

        "Gene_ID": gene_ids,

        "Transcript": transcripts,

        "Feature": ann["Feature"],

        "Region": ann["Region"],

        "Protein_coding": protein_coding,

        "Strand": strands

    })


############################################################
# 11. Save output
############################################################

output = pd.DataFrame(output)


output = output.sort_values(

    ["chr", "pos"]

).reset_index(drop=True)


output.to_csv(

    OUTPUT_FILE,

    sep="\t",

    index=False

)


############################################################
# 12. Summary
############################################################

print("Finished.")

print("--------------------------------------")

print("Output file:", OUTPUT_FILE)

print("Sites:", len(output))


print(
    "Sites with gene:",
    (output["Gene"] != "NA").sum()
)


print(
    "Sites without gene:",
    (output["Gene"] == "NA").sum()
)


############################################################
# Count positions with multiple genes
############################################################

multiple_gene_sites = (

    output["Gene"]
    .str.contains(";", regex=False)
    .sum()

)


print(
    "Sites with multiple genes:",
    multiple_gene_sites
)


############################################################
# Count unique genes
#
# Because several genes can now occur in one cell,
# split them before counting.
############################################################

unique_genes = set()


for genes in output.loc[
    output["Gene"] != "NA",
    "Gene"
]:

    for gene in genes.split(";"):

        unique_genes.add(gene)


print(
    "Unique genes:",
    len(unique_genes)
)


############################################################
# Count sites without GFF annotation
############################################################

print(
    "Sites without GFF annotation:",
    (output["Region"] == "NA").sum()
)

print("--------------------------------------")