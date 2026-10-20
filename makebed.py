import pandas as pd

inp = pd.read_csv("linvilla_1_10_allelefreq.tsv", sep="\t")

# Behalte nur jede fünfzehnte Zeile (pro genomischer Position 15 Timepoints, somit jede 15. Zeile eine neue SNP-Position)
inpf = inp.iloc[::15]

# Nur die ersten beiden Spalten behalten (Chr, Pos)
df = inpf.iloc[:, 1:3]  # Auswahl der ersten zwei Spalten mit iloc

# Neue dritte Spalte hinzufügen: Wert der zweiten Spalte (Index 1) + 1
df['C'] = df.iloc[:, 1] + 1

# Als BED-Datei speichern (tabulatorgetrennt)
df.to_csv('hardfilter.bed', sep='\t', index=False, header=False)  # Kein Header für BED-Format