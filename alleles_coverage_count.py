import pandas as pd
import numpy as np

####################
# Input files (file after variant calling)
files = [
    "filtered_sites_allc_mincov_611_skipcov_3663.tsv"]

results = []
####################
# Create new file with:

for file in files:

    print(f"Bearbeite {file}...")

# Read TSV without header
    df = pd.read_csv(file, sep="\t", header=None, dtype=str)

    ### Column names
    n_cols = df.shape[1]
    names = ["chr", "pos", "ref", "alternative allele"] + [f"tp{i}" for i in range(1, n_cols - 3)]
    df.columns = names

    count_columns = df.columns[4:]

    for _, row in df.iterrows():

        ### Anzahl der Allele (Referenz + Alternative)
        if row["alternative allele"] == ".":
            n_alleles = 1
        else:
            n_alleles = len(row["alternative allele"].split(",")) + 1

        coverages = []

        for col in count_columns:
            counts = [int(x) for x in row[col].split(",") if x != ""]
            coverages.append(sum(counts))

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

out = pd.DataFrame(results)

## Sort by chromosome and position
out = out.sort_values(["chr", "pos"])

####################
# save
out.to_csv(
    "filtered_sites_coverage_multiallelic_sites_covfiltertsv",
    sep="\t",
    index=False
)

print(f"\nFertig!")
print(f"Anzahl Sites: {len(out)}")
print("Gespeichert als: filtered_sites_coverage_multiallelic_sites_covfilter.tsv")