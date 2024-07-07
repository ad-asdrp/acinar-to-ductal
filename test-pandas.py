#imports
import scanpy as sc
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

#reading the CSV file
df = pd.read_csv('hopefullydataset.csv.gz', index_col=0)

#creating 2 separate dataframes
df_2wks = df.iloc[:,:-1116] 
df_4wks = df.iloc[:,1471:] 

#converting dataframs to AnnData objects
df_2wks = df_2wks.T
adata_2wks = sc.AnnData(X=df_2wks.values)

adata_2wks.var_names = df_2wks.columns
adata_2wks.obs_names = df_2wks.index

print(adata_2wks)

df_4wks = df_4wks.T
adata_4wks = sc.AnnData(X=df_4wks.values)

adata_4wks.var_names = df_4wks.columns
adata_4wks.obs_names = df_4wks.index
print(adata_4wks)

def cell_proportion(gene_of_interest, data):
    if gene_of_interest in data.var_names:
        # Get the expression values for the gene of interest
        gene_expression = data[:, gene_of_interest].X

        # Calculate the proportion of cells expressing the gene
        total_cells = data.shape[0]
        expressing_cells = np.sum(gene_expression > 0)  # Count cells where expression > 0
        proportion_expressing = expressing_cells / total_cells
        return proportion_expressing
    else:
        return -1

#proportion of cells expressing gene
gene_of_interest = 'Kras'

two_wks = cell_proportion(gene_of_interest, adata_2wks)
four_wks = cell_proportion(gene_of_interest, adata_4wks)

print('2 weeks proportion of ' + gene_of_interest + ' is ' +  str(two_wks) + '.')
print('4 weeks proportion of ' + gene_of_interest + ' is ' +  str(four_wks) + '.')

proportion_2wks = {}
for i in adata_2wks.var_names:
    proportion_2wks[i] = cell_proportion(i, adata_2wks)

proportion_4wks = {}
for i in adata_4wks.var_names:
    proportion_4wks[i] = cell_proportion(i, adata_4wks)

proportion_differences = {}

for i in adata_2wks.var_names:
    proportion_differences[i] = cell_proportion(i, adata_2wks) - cell_proportion(i, adata_4wks)

sorted_proportion_differences = sorted(proportion_differences.items(),key=lambda x:x[1], reverse=True)
print(sorted_proportion_differences)

sorted_proportion_differences = sorted(proportion_differences.items(),key=lambda x:x[1])
print(sorted_proportion_differences)

gene_of_interest = 'Kras'

difference = cell_proportion(gene_of_interest, adata_2wks) - cell_proportion(gene_of_interest, adata_4wks)
print(difference)

