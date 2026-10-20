import pandas as pd
import numpy as np

####################
# Filter parameters (to adapt to hf filter)
min_reads = 1
max_n_alleles = 10

####################
# Input Data (file After variant calling)
files = [
    "filtered_sites_allc.tsv"
]

results = []

####################
for file in files:

    print(f"Bearbeite {file}...")

    # TSV ohne Header einlesen
    df = pd.read_csv(file, sep="\t", header=None, dtype=str)

    # Spaltennamen
    n_cols = df.shape[1]

    names = (
        ["chr", "pos", "ref", "alternative allele"]
        + [f"tp{i}" for i in range(1, n_cols - 3)]
    )

    df.columns = names

    count_columns = df.columns[4:]

    for _, row in df.iterrows():

        ####################
        # Allees per time
        allele_lists = []

        for col in count_columns:

            counts = [
                int(x) for x in row[col].split(",")
                if x != ""
            ]

            allele_lists.append(counts)

        ####################
        # Number of Alleles
        n_alleles = len(allele_lists[0])

        # Filter: maximale Allelanzahl
        if n_alleles > max_n_alleles:
            continue

        ####################
        # HARD FILTER

        passes_filter = all(
            all(elem >= min_reads for elem in counts)
            for counts in allele_lists
        )

        if not passes_filter:
            continue

        ####################
        # Calculate coverage
        coverages = [
            sum(counts)
            for counts in allele_lists
        ]

        ####################
        # save results
        results.append({
            "chr": row["chr"],
            "pos": int(row["pos"]),
            "ref": row["ref"],
            "alternative allele": row["alternative allele"],
            "n_alleles": n_alleles,
            "coverage_total": np.sum(coverages),
            "coverage_mean_over_time_points": np.mean(coverages),
            "coverage_max_time_point": np.max(coverages)
        })

####################
# new table
out = pd.DataFrame(results)

# Sortieren
out = out.sort_values(["chr", "pos"])

####################
# save
output_file = (
    "filtered_sites_coverage_multiallelic_sites_"
    "minreads1_maxalleles10.tsv"
)

out.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nFertig!")
print(f"Anzahl Sites: {len(out)}")
print(f"Gespeichert als: {output_file}")