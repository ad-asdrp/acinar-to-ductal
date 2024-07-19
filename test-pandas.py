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

#make copy of initial AnnData Object
import copy
adata_2wks_final = copy.deepcopy(adata_2wks)
adata_4wks_final = copy.deepcopy(adata_4wks)

#Filter cells - 2wks
sc.pp.filter_cells(adata_2wks, min_genes=20)

#Filter genes
sc.pp.filter_genes(adata_2wks, min_cells=20)

#Calculate Quality Control metrics
adata_2wks.var['mt'] = adata_2wks.var_names.str.startswith('MT-')  # assuming 'MT-' prefix for mitochondrial genes
sc.pp.calculate_qc_metrics(adata_2wks, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#Filter cells based on QC metrics
adata_2wks = adata_2wks[adata_2wks.obs.n_genes_by_counts < 5000, :]
adata_2wks = adata_2wks[adata_2wks.obs.pct_counts_mt < 5, :]

#Normalize the data
sc.pp.normalize_total(adata_2wks, target_sum=1e4)

#Filter cells - 4wks
sc.pp.filter_cells(adata_4wks, min_genes=20)

#Filter genes
sc.pp.filter_genes(adata_4wks, min_cells=20)

#Calculate Quality Control metrics
adata_4wks.var['mt'] = adata_4wks.var_names.str.startswith('MT-')  # assuming 'MT-' prefix for mitochondrial genes
sc.pp.calculate_qc_metrics(adata_4wks, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#Filter cells based on QC metrics
adata_4wks = adata_4wks[adata_4wks.obs.n_genes_by_counts < 5000, :]
adata_4wks = adata_4wks[adata_4wks.obs.pct_counts_mt < 5, :]

#Normalize the data
sc.pp.normalize_total(adata_4wks, target_sum=1e4)

#Perform PCA
sc.tl.pca(adata_2wks, svd_solver='arpack')
sc.tl.pca(adata_4wks, svd_solver='arpack')

#Compute the neighborhood graph
sc.pp.neighbors(adata_2wks, n_neighbors=10, n_pcs=40)
sc.pp.neighbors(adata_4wks, n_neighbors=10, n_pcs=40)

#Compute UMAP for visualization
sc.tl.umap(adata_2wks)
sc.tl.umap(adata_4wks)

#Cluster the cells
sc.tl.leiden(adata_2wks, resolution=0.1)  # You can adjust the resolution parameter
sc.tl.leiden(adata_4wks, resolution=0.1)  # You can adjust the resolution parameter

#Visualize the clusters
sc.pl.umap(adata_2wks, color=['leiden'])
sc.pl.umap(adata_4wks, color=['leiden'])

sc.tl.dendrogram(adata_2wks, 'leiden')
sc.tl.rank_genes_groups(adata_2wks, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata_2wks, n_genes=20, sharey=False)
sc.pl.rank_genes_groups_heatmap(adata_2wks, n_genes=20, groupby='leiden', show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata_2wks, n_genes=10, groupby='leiden')
sc.pl.rank_genes_groups_matrixplot(adata_2wks, n_genes=10, groupby='leiden')
ranked_genes = adata_2wks.uns['rank_genes_groups']
groups = ranked_genes['names'].dtype.names
for group in groups:
    print(f"\nCluster {group}")
    print(ranked_genes['names'][group][:20])

# Define known acinar/ductal cell markers
acinar_markers = ['Amy2a', 'Ela1']
ductal_markers = ['Pyy','Dclk1', 'Foxq1', 'Krt19', 'S100a10', 'Perp', 'Wfdc2', 'Ctsh', 'Krt7', 'Mmp7', 'Ambp', 'Mmp7', 'Anxa4', 'Slc4a4', 'Krt19', 'Ceacam1']

# Plot the expression of known acinar cell markers to identify acinar cell clusters
sc.pl.umap(adata_2wks, color=ductal_markers)
sc.pl.umap(adata_4wks, color=ductal_markers)

sc.pl.umap(adata_2wks, color=acinar_markers)
sc.pl.umap(adata_4wks, color=acinar_markers)

def calculate_mean_with_pseudocount(data, gene_of_interest, pseudocount):
    # Check if gene_of_interest exists in data.var_names
    if gene_of_interest not in data.var_names:
        raise ValueError(f"Gene {gene_of_interest} not found in data")
    
    # Find column index of gene_of_interest
    idx = np.where(data.var_names == gene_of_interest)[0][0]
    
    # Calculate mean expression with pseudocount
    mean_expr = np.mean(data.X[:, idx] + pseudocount)
    
    return mean_expr

def log_fold_change(data1, data2, gene_of_interest, pseudocount=1e-6):
    # Calculate mean expression with pseudocount
    mean1 = calculate_mean_with_pseudocount(data1, gene_of_interest, pseudocount)
    mean2 = calculate_mean_with_pseudocount(data2, gene_of_interest, pseudocount)
    
    # Calculate fold change
    fold_change = mean2 / mean1
    
    # Calculate log2 fold change
    loggy = np.log2(fold_change)
    
    return loggy

log_fold_change(adata_2wks_final, adata_4wks_final, 'KrasG12D')

#rank LFC
lfc_value = {}

for i in adata_2wks.var_names:
    if i in adata_4wks.var_names:
        lfc_value[i] = log_fold_change(adata_2wks, adata_4wks, i)

sorted_lfc_value = sorted(lfc_value.items(),key=lambda x:x[1], reverse=False)
sorted_lfc_value_2 = sorted(lfc_value.items(),key=lambda x:x[1], reverse=True)

print(sorted_lfc_value)
print("")
print(sorted_lfc_value_2)
