
# Imports, data loading, and data processing like the template

# Determine cells to keep and filter
keep = (adata.obs['pct_counts_mt'] > 2) & (adata.obs['pct_counts_mt'] < 18) & (adata.obs['total_counts'] > 10000) & (adata.obs['total_counts'] < 80000)
print("Removed cells: %d"%(adata.n_obs - sum(keep))) #7353

adata = adata[keep, :]
adata.shape
#(7797, 21566)

# Normalize and log like template, save raw data prior to filtering for highly variable genes

# Find highly variable genes and subset adata object
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
adata = adata[:, adata.var.highly_variable]

# Regress out, scale, PCA like template

# Compute the neighborhood graph. 15 is default so you can also write this as: sc.pp.neighbors(adata)
sc.pp.neighbors(adata, n_neighbors=15)

# Visualize cells via UMAP - use defualt parameters. Just to see what it looks like
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CLDN5', 'CD3D', 'VEGFA', 'GAPDH', 'TOP2A'], use_raw = False, save ='.pdf')
## these marker genes were given - nice to visualize genes that might show up in data (GBM)
## the genes replaced the tutorial (pbmc) genes

##### Change number of clusters (i.e. resolution paramter) - test 0.1, 0.2, 0.3, 0.4
# Clustering
sc.tl.leiden(adata, resolution = 0.3) # Increasing resolution paramter increases the number of clusters
adata.obs['leiden'].describe()['unique'] # This will tell you how many clusters you have with the given resolution
sc.pl.umap(adata, color=['leiden', 'CLDN5', 'CD3D', 'TOP2A'], use_raw = False, save='_leiden_0.3.pdf') # change file name based on resolution used

# Find marker genes - this step takes some time sometimes. Run it by itself and let if finish out before next step.
sc.tl.rank_genes_groups(adata, 'leiden', corr_method='benjamini-hochberg', method='wilcoxon')

# Show the 10 top ranked genes per cluster
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(10) ## this is top 10 genes, you can change 10 to any number of genes you want to look at
    # You can also replace "pvals" with "logfoldchanges" to look at a different statistic. Might help.
tmp = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}) # all genes in your adata object
tmp.to_csv("marker_genes_leiden_0.1.csv") # dumps a csv file into your working directory so you can see all of the marker genes and their p-values for each cluster
# Use this csv file to compare genes to the list of marker genes given
# Tentatively have an idea of what each cluster belongs to but DO NOT LABEL yet.

########################## STOP HERE ###############################################

### Add all your output (UMAP and marker genes) from resolution = 0.1 to PowerPoint so you have it organized together.
### Research marker genes for each cluster - note their function, if it is a marker gene for a specific cell type, etc. This will be helpful for cluster identification.
### Start to map each cluster with a cell type (does not have to be specific). Can be as simple as "cluster 0 has a lot of immune cell markers"
### Keep notes!!

### Once you have everything for the first resolution tested in a PowerPoint, go back to line 31 and change your resolution.
### Run lines 31-47 again. Make sure to change your pdf file names with the correct resolution
### Save all of your output to the PowerPoint and do the same thing as before.

##### After you complete all of your resolutions, look at all your plots and notes together to choose the optimal number of clusters for this data.
##### Note your resolution value that gave you that number of clusters and RERUN your code again at line 31 now with your optimal resolution.

### Proceed with the rest of the template

### By picking an optimal resolution, you should know: 1) your number of clusters and your cell identities (it is okay if you have an unknown cell identity)


# Identify cell types - the length of this list should equal your number of clusters
# The order of this list matches the order of numbers. If you think cluster 0 is T-cells, T-cells should be the first name in this list.
new_cluster_names = [
    'A', 'B',
    'C', 'D',
    'E', ....]
adata.rename_categories('leiden', new_cluster_names)

# Plot with labels
sc.pl.umap(adata, color='leiden', use_raw = False, save='_new_idents.pdf')


# Play with the other plots - given with template
