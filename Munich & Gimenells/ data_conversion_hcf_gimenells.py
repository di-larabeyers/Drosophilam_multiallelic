import pandas as pd
import numpy as np


#change the name of your population here:
pop="gimenells"

#estimation of generations between the timepoints based on sample dates and Drosophila generations for summer (10) and winter (2); needs to be changed when looking at another data set (other sampling dates)
gens = [0, 10, 12, 22, 24, 34, 35, 36, 44, 46, 48, 60, 62]

print("importing and merging tsv-files")

#import tsv-file (change path here)
zl = pd.read_csv("filtered_sites_cov_resistence_genes.tsv", sep="\t", header=None)


# DataFrames untereinander hängen
inp = pd.concat([zl], ignore_index=True)



#count the number of columns in your input
column_count = inp.shape[1]

print("converting allele frequencies to lists")
#convert allele frequencies at all timepoints in lists
for column in inp.columns[4:column_count]:
    inp[column] = inp[column].apply(lambda x: [int(i) for i in x.split(",")])


print("generating lists")

def generate_lists(row, column_count):

    #list to save all generated lists for this row
    rowlists = []

    for column_index in range(4, column_count):  # start at column 4
        current_column = row[column_index]     # current column
        n_alleles = len(row[4])                 # number of alleles in current row
        
        # make dynamic list based on current row
        new_list = [
        pop,                   # name of population (should be the same for the whole dataframe)
        row[0],                # name of chromosome
        int(row[1]),         # position on the chromosome
        n_alleles,           # number of alleles (should be the same for the whole row)
        column_index - 3      # index for the number of alleles
        ] + [elem / sum(current_column) if sum(current_column) != 0 else 0 for elem in current_column] # allele frequencies
        rowlists.append(new_list) # save all generated lists for this row in rowlists
    
    return rowlists


# generate lists for each row
all_generated_lists = inp.apply(lambda row: generate_lists(row, column_count), axis=1)

# merge all the lists into a new dataframe
print("merging lists to new dataframe")
flattened_data = [item for sublist in all_generated_lists for item in sublist]
new_df = pd.DataFrame(flattened_data)
new_df[len(new_df.columns)] = new_df.iloc[:, 5:].sum(axis=1)

# add column for the generation
new_df[len(new_df.columns)] = np.tile(gens, len(new_df) // len(gens) + 1)[:len(new_df)]


# make a header for the new dataframe
header = ["pop", "chr", "pos", "n_alleles",	"time_point"]

# starting wih column 5 names depend on number of columns in the whole dataframe
header += [f"af_{i-4}" for i in range(5, len(new_df.columns) - 2)]

header += ["sum", "gen"]

# give the new dataframe the generated header
new_df.columns = header

# name the file based on the population name
filename = pop +"_all_resistencegenes_hcf_converted.tsv"

# generate new tsv file from the new converted dataframe
print("generating new tsv-file")
new_df.to_csv(filename, sep="\t", index=False)
