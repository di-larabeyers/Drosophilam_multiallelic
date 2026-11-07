#!/usr/bin/env python3

import pandas as pd

############################################################
# Files
############################################################

INPUT = "filtered_sites_coverage_multiallelic_sites_covfilter.tsv"
OUTPUT = "coverage_summary_100kb_windows.tsv"

WINDOW_SIZE = 100000

############################################################
# Read input
############################################################

print("Reading input...")

df = pd.read_csv(INPUT, sep="\t")

print(f"{len(df)} multiallelic sites loaded.\n")

############################################################
# Assign genomic windows
############################################################

df["window_start"] = (df["pos"] // WINDOW_SIZE) * WINDOW_SIZE
df["window_end"] = df["window_start"] + WINDOW_SIZE - 1

############################################################
# Summarize windows
############################################################
summary = (
    df
    .groupby(["chr", "window_start", "window_end"], as_index=False)
    .agg(
        n_sites=("pos", "count"),

        # Summe der Coverage über alle Zeitpunkte
        tot_cov_per_window=("coverage_total", "sum"),

        # Durchschnitt der bereits über die Zeitpunkte gemittelten Coverage
        # aller multiallelischen Sites im jeweiligen 100-kb-Fenster
        mean_cov_per_time_window=("coverage_mean_over_time_points", "mean"),

        # Maximale Coverage eines Locus in einem einzelnen Zeitpunkt
        max_cov_per_site=("coverage_max_time_point", "max"),

        # Mittlere Allelanzahl
        mean_n_alleles=("n_alleles", "mean")
    )
)
############################################################
# Site density
############################################################

summary["sites_per_kb"] = (
    summary["n_sites"] /
    (WINDOW_SIZE / 1000)
)

############################################################
# Sort output
############################################################

summary = summary.sort_values(
    ["chr", "window_start"]
).reset_index(drop=True)

############################################################
# Save
############################################################

summary.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

############################################################
# Summary
############################################################

print("-------------------------------------------")
print(f"Window size          : {WINDOW_SIZE:,} bp")
print(f"Number of windows    : {len(summary)}")
print("-------------------------------------------")

print(summary.head())

print("\nOutput written to:")
print(OUTPUT)