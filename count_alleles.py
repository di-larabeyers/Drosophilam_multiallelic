
import pandas as pd

#import tsv-file (change path here)
# that is the tsv file you originally got from the data_conversion script and also used as input for the F-Matrix script; hcf means hard and coverage filter
inp = pd.read_csv("linvilla_allc_hcf_converted.tsv", sep="\t", usecols=["n_alleles"])

# Häufigkeit der Werte zählen
counts = inp["n_alleles"].value_counts()

# Umwandeln in einen neuen DataFrame
result_df = counts.reset_index()

result_df.columns = ['n_alleles', 'frequency']

# Teilen der Werte in der Spalte "Häufigkeit" durch 15
result_df['frequency'] = result_df['frequency'] / 15

result_df.to_csv("n_alleles.tsv", sep="\t", index=False)