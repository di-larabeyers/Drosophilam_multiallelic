# pandas package is needed to handle the dataframes
# numpy is needed for adding the generations column
import pandas as pd
import numpy as np

###################################################################################################################
#your settings:

#change the name of your population here:
pop="linvilla"

#estimation of generations between the timepoints for Linvilla based on sample dates and Drosophila generations for summer (10) and winter (2); needs to be changed when looking at another data set (other sampling dates)
gens = [0, 10, 12, 22, 24, 34, 35, 36, 44, 46, 48, 60, 70, 72, 82]

# set a hard filter (change here):
min_reads = 1
max_n_alleles = 10

###################################################################################################################

print("importing tsv-file") # to see status while the program is running

#import tsv-file (change path here)
inp= pd.read_csv("filtered_sites_allc.tsv", sep="\t" , header=None)

#count the number of columns in your input
column_count = inp.shape[1]


print("converting allele frequencies to lists") # to see status while the program is running

#convert allele frequencies at all timepoints in lists (needed since they come comma seperated in a tab seperated file)
for column in inp.columns[4:column_count]:
    inp[column] = inp[column].apply(lambda x: [int(i) for i in x.split(",")])


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

def apply_with_condition(row, column_count):
    if all(all(elem >= min_reads for elem in lst) for lst in row[4:]) and len(row[4]) <= max_n_alleles:   # hard filter
        return generate_lists(row, column_count)  # if condition is true: list is generated and returned
    # function returns None value if condition is false

print("generating lists for each row") # to see status while the program is running

# generate lists for each row
all_generated_lists = inp.apply(lambda row: apply_with_condition(row, column_count), axis=1)

# delete all None values that might have been generated due to the if for the hard filter
filtered_lists = all_generated_lists.dropna().tolist()  # Konvertiert das Ergebnis in eine Liste ohne None-Werte

# merge all the lists into a new dataframe
print("merging lists to new dataframe") # to see status while the program is running

flattened_data = [item for sublist in filtered_lists for item in sublist] # flatten the lists in lists
new_df = pd.DataFrame(flattened_data)

###################################################################################################################

# add column with the sum of all allele frequencies in one row of the new dataframe (quality control, should be 1.0)
new_df[len(new_df.columns)] = new_df.iloc[:, 5:].sum(axis=1)


# add column for the generation
new_df[len(new_df.columns)] = np.tile(gens, len(new_df) // len(gens) + 1)[:len(new_df)]

n_pos = len(new_df)/15

print(f"number of positions: {n_pos}")



###################################################################################################################
# make a header for the new dataframe

# first 5 headers are always the same
header = ["pop", "chr", "pos", "n_alleles","time_point"]

# starting wih column 5 names depend on number of columns in the whole dataframe
header += [f"af_{i-4}" for i in range(5, len(new_df.columns) - 2)]

# last 2 headers are always the same
header += ["sum", "gen"]

# give the new dataframe the generated header
new_df.columns = header




# name the file based on the population and used filters
filename = pop + "_" + str(min_reads) + "_" + str(max_n_alleles) + "_allelefreq.tsv"

# generate new tsv file from the new converted dataframe
print("generating new tsv-file") # to see status while the program is running
new_df.to_csv(filename, sep="\t", index=False)