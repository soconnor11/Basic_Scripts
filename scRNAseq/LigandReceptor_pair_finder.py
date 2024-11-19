########################################################################
## Mono- and Triculture scRNA-seq analysis                           ##
##  ______     ______     __  __                                      ##
## /\  __ \   /\  ___\   /\ \/\ \                                     ##
## \ \  __ \  \ \___  \  \ \ \_\ \                                    ##
##  \ \_\ \_\  \/\_____\  \ \_____\                                   ##
##   \/_/\/_/   \/_____/   \/_____/                                   ##
## @Developed by: Plaisier Lab                                        ##
##   (https://plaisierlab.engineering.asu.edu/)                       ##
##   Arizona State University                                         ##
##   242 ISTB1, 550 E Orange St                                       ##
##   Tempe, AZ  85281                                                 ##
## @github: https://github.com/plaisier-lab/gbmTriculture             ##
## @dockerhub:  https://hub.docker.com/r/cplaisier/scrna_seq_velocity ##
## @Author:  Samantha O'Connor and Chris Plaisier                     ##
## @License:  GNU GPLv3                                               ##
##                                                                    ##
## If this program is used in your analysis please                    ##
## mention who built it. Thanks. :-)                                  ##
########################################################################

## Uses the cplaisier/scrna_seq_velocity Docker image from Docker Hub:  https://hub.docker.com/r/cplaisier/scrna_seq_velocity
# docker run -it -v '<local file path>:/files' cplaisier/scrna_seq_velocity
# Be sure to replace the <local file path> with the path to your scRNA-seq data

#---------------------------------- Updated 05/10/22 -----------------------------------#
# Identify upregulated receptors in triculture GSCs relative to monoculture GSCs
# Filter receptors based off presence of secreted ligands in other triculture cell populations (HUVECS, ASTROCYTES)
# Filter receptors to those associated with at least one upregulated pathway in triculture GSCs

#------------------------
# Setup section
#-----------------------

# Imports
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom
import mygene
mg = mygene.MyGeneInfo()
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Set up file structure
tag = 'mono'
newpath = 'results'
newpath1 = newpath+'/'+tag
if not os.path.exists(newpath):
    os.makedirs(newpath)

if not os.path.exists(newpath1):
    os.makedirs(newpath1)

#------------------------
# Load monoculture data
#-----------------------

print('Loading monoculture scRNA-seq...')
# Build scanpy data object
Mono_adata = sc.read_10x_mtx(
        './MN1_monoculture/filtered_feature_bc_matrix/',
        var_names = 'gene_ids',
        cache=True)
Mono_adata.var_names_make_unique()
Mono_adata
Mono_adata.shape #(2573, 36601)

# Add copy number data from CONICSmap analysis to adata object
Mono_PP_CONICS = pd.read_csv("MN1_posterior_probabilities_CONICSmap.csv", index_col = 0)
Mono_adata.obs = pd.concat([Mono_adata.obs,Mono_PP_CONICS], axis=1)
Mono_adata.obs['4q_new']= 1-Mono_adata.obs['4q'].abs()
Mono_adata.obs

# Initial default filtering of cells and genes
sc.pp.filter_cells(Mono_adata, min_genes=200)
sc.pp.filter_genes(Mono_adata, min_cells=3)

# Generate quality control data
mito_genes = Mono_adata.var['gene_symbols'].str.startswith('MT-')
Mono_adata.obs['percent_mito'] = np.sum(Mono_adata[:, mito_genes].X, axis=1).A1 / np.sum(Mono_adata.X, axis=1).A1
Mono_adata.obs['n_counts'] = Mono_adata.X.sum(axis=1).A1
Mono_adata.obs['n_genes'] = [i.count_nonzero() for i in Mono_adata.X]

# Secondary data dependent filtering
Mono_adata = Mono_adata[Mono_adata.obs.n_genes > 2000, :]
Mono_adata = Mono_adata[Mono_adata.obs.percent_mito < 0.18, :]
Mono_adata = Mono_adata[Mono_adata.obs.n_counts > 5000, :]
print(Mono_adata.shape) #(2109, 21466)

# Downstream processing
sc.pp.normalize_total(Mono_adata, target_sum=1e4)
sc.pp.log1p(Mono_adata)
Mono_adata.raw = Mono_adata
sc.pp.highly_variable_genes(Mono_adata, n_top_genes=4000)
Mono_adata2 = Mono_adata[:, Mono_adata.var.highly_variable] # save as new variable so can go back if needed
sc.pp.regress_out(Mono_adata2, ['n_counts', 'percent_mito'])
sc.pp.scale(Mono_adata2, max_value=10)
sc.tl.pca(Mono_adata2, svd_solver='arpack')


sc.pp.neighbors(Mono_adata2) #, n_neighbors=10, n_pcs=40)
sc.tl.umap(Mono_adata2)
# Cluster
sc.tl.leiden(Mono_adata2, resolution = 0.13)

# Label clusters with marker genes
marker_genes = ['ENSG00000261371', 'ENSG00000101463', 'ENSG00000147889', 'ENSG00000160307']
new_cluster_names = ['HUVECS1', 'HUVECS2' , 'GSCs', 'ASTROCYTES']
Mono_adata2.rename_categories('leiden', new_cluster_names)

# Combine HUVEC clusters together
old_to_new = dict(
    HUVECS1='HUVECS',
    HUVECS2='HUVECS',
    GSCs='GSCS',
    ASTROCYTES='ASTROCYTES'
)
Mono_adata2.obs['new_clusters'] = (
    Mono_adata2.obs['leiden']
    .map(old_to_new)
    .astype('category')
)

# Differential expression between clusters
sc.tl.rank_genes_groups(Mono_adata2, 'new_clusters', corr_method='benjamini-hochberg', method='wilcoxon')
# Save out marker gene names
pd.DataFrame(Mono_adata2.uns['rank_genes_groups']['names']).to_csv(newpath1+'/'+tag+'_marker_genes_names.csv')

# Get a table with scores and results
result = Mono_adata2.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']}).to_csv(newpath1+'/'+tag+'_marker_genes.csv')


# Plot UMAP of merged clusters
os.chdir(newpath1)

sc.pl.umap(Mono_adata2, color='new_clusters', legend_loc='right margin', title='', frameon=True, save='_'+ tag+'_clusters_labeled.pdf')

# Put all figures together on same pdf
genes1 = ['PECAM1','S100B','CDKN2A']
color1 = 'binary'
convGenes = [list(Mono_adata2.var.index[Mono_adata2.var['gene_symbols']==i])[0] for i in genes1]
sc.pl.umap(Mono_adata2, color=convGenes, color_map = color1, legend_loc='right margin', frameon=True, ncols=1, vmin=0, vmax=1, save='_'+tag+'_gene_expression_binary.pdf')
chroms = ['17q','4q_new']
sc.pl.umap(Mono_adata2, color=chroms, color_map = color1,legend_loc='right margin', frameon=True, ncols=1, vmin=0, vmax=1, save='_'+tag+'_chrom_expression_binary.pdf')

"""
# Generate separate pdfs for marker genes
for gene in genes1:
    convGene = Mono_adata2.var.index[Mono_adata2.var['gene_symbols']==gene][0]
    sc.pl.umap(Mono_adata2, color=convGene, color_map = color1, vmin=0, vmax=1, save='_'+tag+'_'+gene+'.pdf')

# Generate separate pdfs for copy numbers
chroms = ['4q_new', '17q']
for chrom in chroms:
    sc.pl.umap(Mono_adata2, color=chrom, color_map = color1, vmin=0, vmax=1, save = '_'+tag+'_'+chrom+'.pdf')
"""

#------------------------
# Load triculture data
#-----------------------

# Set up file structure
os.chdir('../../')
tag = 'tri'
newpath2 = newpath+'/'+tag
if not os.path.exists(newpath2):
    os.makedirs(newpath2)

print('Loading triculture scRNA-seq...')

Tri_adata = sc.read_10x_mtx(
        './MN2_triculture/filtered_feature_bc_matrix/',
        var_names = 'gene_ids',
        cache=True)
Tri_adata.var_names_make_unique()
Tri_adata
Tri_adata.shape #(2645, 36601)

# Add copy number data from CONICSmap analysis
Tri_PP_CONICS = pd.read_csv("MN2_posterior_probabilities_CONICSmap.csv", index_col = 0)
Tri_adata.obs = pd.concat([Tri_adata.obs,Tri_PP_CONICS], axis=1)
Tri_adata.obs['4q_new']= 1-Tri_adata.obs['4q'].abs()
Tri_adata.obs

# Initial default filtering of cells and genes
sc.pp.filter_cells(Tri_adata, min_genes=200)
sc.pp.filter_genes(Tri_adata, min_cells=3)

# Generate quality control metrics
mito_genes = Tri_adata.var['gene_symbols'].str.startswith('MT-')
Tri_adata.obs['percent_mito'] = np.sum(Tri_adata[:, mito_genes].X, axis=1).A1 / np.sum(Tri_adata.X, axis=1).A1
Tri_adata.obs['n_counts'] = Tri_adata.X.sum(axis=1).A1
Tri_adata.obs['n_genes'] = [i.count_nonzero() for i in Tri_adata.X]

# Data dependent filtering
Tri_adata = Tri_adata[Tri_adata.obs.n_genes > 2000, :]
Tri_adata = Tri_adata[Tri_adata.obs.percent_mito < 0.18, :]
Tri_adata = Tri_adata[Tri_adata.obs.n_counts > 5000, :]
print(Tri_adata.shape) #(1947, 21307)

# Downstream processing
sc.pp.normalize_total(Tri_adata, target_sum=1e6)
sc.pp.log1p(Tri_adata)
Tri_adata.raw = Tri_adata
sc.pp.highly_variable_genes(Tri_adata, n_top_genes=4000)
Tri_adata2 = Tri_adata[:, Tri_adata.var.highly_variable]
sc.pp.regress_out(Tri_adata2, ['n_counts', 'percent_mito'])
sc.pp.scale(Tri_adata2, max_value=10)
sc.tl.pca(Tri_adata2, svd_solver='arpack')


sc.pp.neighbors(Tri_adata2) #, n_neighbors=10, n_pcs=40)
sc.tl.umap(Tri_adata2)

# Cluster
sc.tl.leiden(Tri_adata2, resolution = 0.15)

# Label clusters with marker genes
marker_genes = ['ENSG00000261371', 'ENSG00000148773', 'ENSG00000078018', 'ENSG00000147889', 'ENSG00000171848']
new_cluster_names = ['HUVECS1', 'HUVECS2', 'ASTROCYTES', 'GSCs', 'UNKNOWN']
Tri_adata2.rename_categories('leiden', new_cluster_names)

# Combine HUVECS clusters together
old_to_new = dict(
    HUVECS1='HUVECS',
    HUVECS2='HUVECS',
    ASTROCYTES='ASTROCYTES',
    GSCs='GSCS',
    UNKNOWN='UNKNOWN'
)
Tri_adata2.obs['new_clusters'] = (
    Tri_adata2.obs['leiden']
    .map(old_to_new)
    .astype('category')
)


# Differential expression between clusters
sc.tl.rank_genes_groups(Tri_adata2, 'new_clusters', corr_method='benjamini-hochberg', method='wilcoxon')
# Save out marker gene names
pd.DataFrame(Tri_adata2.uns['rank_genes_groups']['names']).to_csv(newpath2+'/'+tag+'_marker_genes_names.csv')

# Get a table with names and results
result = Tri_adata2.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']}).to_csv(newpath2+'/'+tag+'_marker_genes.csv')


# plot UMAP of merged clusters
os.chdir(newpath2)
sc.pl.umap(Tri_adata2, color='new_clusters', legend_loc='right margin', title='', frameon=True, save='_'+tag+'_clusters_labeled.pdf')


# Put all figures together on same pdf
genes1 = ['PECAM1','S100B','CDKN2A']
color1 = 'binary'
convGenes = [list(Tri_adata2.var.index[Tri_adata2.var['gene_symbols']==i])[0] for i in genes1]
sc.pl.umap(Tri_adata2, color=convGenes, color_map = color1, legend_loc='right margin', frameon=True, ncols=1, vmin=2, vmax=8, save='_'+tag+'_gene_expression_binary.pdf')
chroms = ['17q','4q_new']
sc.pl.umap(Tri_adata2, color=chroms, color_map = color1,legend_loc='right margin', frameon=True, ncols=1, vmin=0, vmax=1, save='_'+tag+'_chrom_expression_binary.pdf')



#-----------------------------------------------------------------------------------
# Integrate monoculture and triculture data to look at differential expression
#-----------------------------------------------------------------------------------

# Set up file structure
os.chdir('../../')
tag = 'integrated'
newpath3 = newpath+'/'+tag
if not os.path.exists(newpath3):
    os.makedirs(newpath3)


print('Integrate scRNA-seq datasets...')

# Begin integration
var_names = Mono_adata.var_names.intersection(Tri_adata.var_names)
len(var_names) # 20457 common genes
# Subset data to common genes
Mono_adata3 = Mono_adata[:, var_names]
Tri_adata3 = Tri_adata[:, var_names]
Mono_adata3.obs = Mono_adata2.obs
Tri_adata3.obs = Tri_adata2.obs

adata_concat = Mono_adata3.concatenate(Tri_adata3, batch_categories=['Mono', 'Tri'])
adata_concat.obs.new_clusters = adata_concat.obs.new_clusters.astype('category')
adata_concat.shape #(4056, 20457)
# 4056 cells; 20457 genes

# Using BBKNN
#sc.tl.pca(adata_concat)
#sc.external.pp.bbknn(adata_concat, batch_key='batch')
#sc.tl.umap(adata_concat)
#sc.pl.umap(adata_concat, color=['batch', 'new_clusters'], save='_MN1_MN2_concat_BBKNN.pdf')

## Create new observation
tmp1 = adata_concat.obs['batch']
tmp2 = adata_concat.obs['new_clusters']
tmp3 = []
for i in range(len(tmp1)):
    tmp3.append(tmp1.iloc[i]+'_'+tmp2.iloc[i])

adata_concat.obs['batch_new_clusters'] = tmp3
adata_concat.obs['batch_new_clusters'].value_counts()
tmp = pd.DataFrame(adata_concat.obs['batch_new_clusters'].value_counts())
tmp.rename(columns = {'batch_new_clusters':'cell_counts'}).to_csv(newpath3+'/'+tag+'_cell_counts.csv')

# Remove triculture unknown cluster
adata_concat = adata_concat[adata_concat.obs['batch_new_clusters']!='Tri_UNKNOWN']
adata_concat.shape #(4024, 20457)

# Normalize/prep adata_concat
sc.pp.normalize_total(adata_concat, target_sum=1e6)
sc.pp.log1p(adata_concat)
sc.tl.pca(adata_concat, svd_solver='arpack')
sc.pp.neighbors(adata_concat)
sc.tl.umap(adata_concat)

"""
# Duplicate adata_concat to calculate expression across ALL cells in culture condition
new_adata_concat = adata_concat.copy()
tmp4 = adata_concat.obs['batch']
tmp5 = []
for i in range(len(tmp4)):
    tmp5.append(tmp4.iloc[i]+'_all')

new_adata_concat.obs['batch_new_clusters'] = tmp5
new_new_adata_concat = adata_concat.concatenate(new_adata_concat)

# Plot receptors
os.chdir(newpath3)
receptors = ['PDGFRA', 'LGR6', 'FPR1', 'FGFR4', 'LRP8', 'F3']
dict = {'PDGFRA ligands': ['PDGFA', 'PDGFB', 'PDGFC', 'PDGFD'], 'LGR6 ligand': ['RSPO3'], 'FPR1 ligand': ['SAA1'], 'FGFR4 ligand': ['FGF5'], 'LRP8 ligands': ['APOE', 'LRPAP1', 'RELN'], 'F3 ligand': ['IL6', 'TFPI']}

sc.pl.matrixplot(adata_concat, receptors, gene_symbols = 'gene_symbols', groupby = 'batch_new_clusters', dendrogram=True, cmap='Blues', standard_scale='var', swap_axes = True, save = '_'+tag+'_receptors.pdf')
sc.pl.matrixplot(adata_concat, dict, gene_symbols = 'gene_symbols', groupby = 'batch_new_clusters', dendrogram=True, cmap='Purples', standard_scale='var', save = '_'+tag+'_receptors_and_ligands_axes_swap.pdf')


sc.pl.dotplot(new_new_adata_concat, receptors, gene_symbols = 'gene_symbols', groupby='batch_new_clusters', save = '_'+tag+'_test.pdf')
sc.pl.matrixplot(new_new_adata_concat, receptors, gene_symbols = 'gene_symbols', groupby = 'batch_new_clusters', dendrogram=True, cmap='Reds', standard_scale='var', save = '_'+tag+'_test.pdf')
"""


#------------------------ Gene conversion (ensembl & gene symbols)--------------------------------#

gene_symbols_conv = pd.DataFrame(adata_concat.var['gene_symbols']).reset_index()
geneconversion = mg.querymany(list(gene_symbols_conv['index']), scopes='ensembl.gene')
df = pd.DataFrame(geneconversion)
df1 = df[['query', 'symbol', 'entrezgene']].rename(columns={'query':'ensembl'})
df2 = df.set_index('query')['entrezgene']

# Integrate entrez gene IDs into adata_concat
adata_concat.var['entrez_gene'] = df2


#----------------------------------------------------
# Find genes DE in triculture data (compared to mono)
#-----------------------------------------------------

# Find upregulated genes in triculture GSCs compared to monoculture GSCs
sc.tl.rank_genes_groups(adata_concat, groupby='batch_new_clusters', groups = ['Tri_GSCS'], reference = 'Mono_GSCS', use_raw=False, n_genes = False)

# Get a table with names and results
result = adata_concat.uns['rank_genes_groups']
groups = result['names'].dtype.names
tmp = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']})


# Format differentially expressed GSC genes
Tri_DE_GSC_all = tmp.rename(columns ={'Tri_GSCS_n':'genes', 'Tri_GSCS_l':'fold_change', 'Tri_GSCS_p':'pvals_adj'})
info1 = pd.merge(Tri_DE_GSC_all, df1, how='left',right_on='ensembl',left_on='genes')[['genes', 'symbol','entrezgene','fold_change','pvals_adj']]
info1.to_csv('Differentially_expressed_tri_GSC_genes.csv')
Tri_genes_with_cutoffs = info1[(info1['fold_change']>=0.3) & (info1['pvals_adj']<=0.05)]
len(Tri_genes_with_cutoffs) #525
# 5/10/22: 1553
Tri_genes_with_cutoffs.to_csv('Upregulated_tri_GSC_genes.csv')


# Identify receptor presence in triculture GSCs
adata_concat_Tri_GSCs = adata_concat[adata_concat.obs['batch_new_clusters'] =='Tri_GSCS']
len(adata_concat_Tri_GSCs) #217
adata_concat_Tri_GSCs_2 = adata_concat_Tri_GSCs[:,Tri_genes_with_cutoffs['genes']]
GSC_DEGs_present = adata_concat_Tri_GSCs_2.var_names[(np.sum(adata_concat_Tri_GSCs_2.X.todense() > 0, axis=0)/(len(adata_concat_Tri_GSCs_2.obs))>=0.10).tolist()[0]].tolist()
len(GSC_DEGs_present)
# 491 present
# 5/10/22: 1481

# Differentially expressed genes based on FC >= 0.3 or 1.0 (used later on in analysis - downstream pathways)
Tri_GSC_DE_genes = info1[info1['fold_change']>=0.3]
len(Tri_GSC_DE_genes) # 4281
# 5/10/22: 7837

Tri_GSC_DE_genes_FC1 = info1[info1['fold_change']>=1.0]
len(Tri_GSC_DE_genes_FC1) # 2087
# 5/10/22/: 3192

Tri_genes_with_cutoffs_strict = Tri_genes_with_cutoffs[Tri_genes_with_cutoffs['fold_change']>=1.0]


###############################################
## Identify putative ligand/receptor pairs   ##
## from Ramilowski_et_al_2015                ##
## PMID = 26198319                           ##
## Hypothesis:  Receptor in GSC and secreted ##
## ligand from HUVEC or Astrocyte            ##
###############################################

os.chdir('../')
print('Identify putative ligand-receptor pairs...')

## Load up receptor and ligand Pairs
# Load receptor-ligand pairs from Ramilowski
LR_pairs = pd.read_excel("Ramilowski_et_al_2015.xlsx", sheet_name = 'All.Pairs')

# Literature supported interactions
LR_pairs_LS = LR_pairs[LR_pairs['Pair.Evidence']=='literature supported']
len(LR_pairs_LS) #1894
receptors_LS = pd.DataFrame(LR_pairs_LS['Receptor.ApprovedSymbol'])
receptors1_LS = list(set(list(receptors_LS['Receptor.ApprovedSymbol'])))
receptors_df_LS = pd.DataFrame(data = receptors1_LS, columns = ['receptors'])
len(receptors_df_LS) #589 unique receptors

# Ligand info
ligands_LS = pd.DataFrame(LR_pairs_LS['Ligand.ApprovedSymbol'])
len(ligands_LS['Ligand.ApprovedSymbol'].unique()) #642


# Subset pairs to pairs wtih secreted ligands ( move secreted ligand filter up 01/04/22 )
secreted_ligands = pd.read_csv('protein_class_Secreted.csv') #2943
secr_ligands = np.intersect1d(list(set(list(ligands_LS['Ligand.ApprovedSymbol']))), list(secreted_ligands['Gene']))
len(secr_ligands) #519
np.setdiff1d(ligands_LS['Ligand.ApprovedSymbol'].unique(), secr_ligands)

secr_ligand_pairs_Ram = LR_pairs_LS[LR_pairs_LS['Ligand.ApprovedSymbol'].isin(secr_ligands)]
len(secr_ligand_pairs_Ram) # 1504 pairs

receptors_LS = pd.DataFrame(secr_ligand_pairs_Ram['Receptor.ApprovedSymbol'])
receptors1_LS = list(set(list(receptors_LS['Receptor.ApprovedSymbol'])))
receptors_df_LS = pd.DataFrame(data = receptors1_LS, columns = ['receptors'])
len(receptors_df_LS) # 496 unique receptors

# Ligand info
ligands_LS = pd.DataFrame(secr_ligand_pairs_Ram['Ligand.ApprovedSymbol'])
len(ligands_LS['Ligand.ApprovedSymbol'].unique()) #519

# Literature supporting receptors
Tri_GSC_DE_Receps_LS = pd.merge(Tri_genes_with_cutoffs, receptors_df_LS, how = 'left', right_on='receptors', left_on = 'symbol').sort_values(by=['fold_change'], ascending=False)
Tri_GSC_DE_Receps_LS_2 = Tri_GSC_DE_Receps_LS[Tri_GSC_DE_Receps_LS['receptors'].notna()].drop(columns = ['receptors'])
len(Tri_GSC_DE_Receps_LS_2) #30
Tri_GSC_DE_Receps_LS_2.to_csv(newpath3+'/'+tag+'_upregulated_literature_supp_receptors_triculture_GSCs.csv')



## Subset receptors into only receptors present in >10% triculture upregulated GSCs
# Subset adata_concat so only looking at Triculture GSCs
adata_concat_Tri = adata_concat[adata_concat.obs['batch'] =='Tri']
adata_concat_Tri_GSCS = adata_concat_Tri[adata_concat_Tri.obs['new_clusters']=='GSCS']

## Identify presence of receptors in triculture GSCs
adata_concat_Tri_receptors = adata_concat_Tri_GSCS[:,Tri_GSC_DE_Receps_LS_2['genes']]
GSC_receptors_present = adata_concat_Tri_receptors.var_names[(np.sum(adata_concat_Tri_receptors.X.todense() > 0, axis=0)/(len(adata_concat_Tri_receptors.obs))>=0.10).tolist()[0]].tolist()
len(GSC_receptors_present)
# 29 receptors > 0.1
tmp = pd.merge(df1, pd.DataFrame(GSC_receptors_present, columns = ['present_receptors']), how='left', right_on='present_receptors', left_on = 'ensembl')
GSC_receptors_present_sym = tmp[tmp['present_receptors'].notna()].drop(columns = ['present_receptors'])
np.setdiff1d(list(Tri_GSC_DE_Receps_LS_2['symbol']), list(GSC_receptors_present_sym['symbol']))
# filtered out GLP2R
Tri_GSC_DE_Receps_LS_2 = Tri_GSC_DE_Receps_LS_2[Tri_GSC_DE_Receps_LS_2['symbol']!='GLP2R']

### Look at receptors corresponding ligands
## See if ligands are present in other clusters
# But some receptors have many different ligands!!
y = secr_ligand_pairs_Ram[secr_ligand_pairs_Ram['Receptor.ApprovedSymbol'].isin(Tri_GSC_DE_Receps_LS_2['symbol'])]
len(y) #78
another_y = secr_ligand_pairs_Ram[secr_ligand_pairs_Ram['Pair.Evidence']=='literature supported'][secr_ligand_pairs_Ram[secr_ligand_pairs_Ram['Pair.Evidence']=='literature supported']['Receptor.ApprovedSymbol'].isin(Tri_GSC_DE_Receps_LS_2['symbol'])]
len(another_y) #78

ligands1 = list(set(list(y['Ligand.ApprovedSymbol'])))
len(ligands1) #71
y1 = y[['Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol']]
y1['pairs'] = y1.values.tolist()

## Narrow down to secreted ligands
#secreted_ligands = pd.read_csv('protein_class_Secreted.csv')
#common_secr_ligands = np.intersect1d(list(set(list(y1['Ligand.ApprovedSymbol']))), list(secreted_ligands['Gene']))
#len(common_secr_ligands) #91
#y2 = y1[y1['Ligand.ApprovedSymbol'].isin(common_secr_ligands)]
# len(list(set(list(y2['Receptor.ApprovedSymbol']))))
# 103 secreted ligands that match up to 28 unique receptors
#np.setdiff1d(np.array(Tri_GSC_DE_Receps_LS_2['symbol']), np.array(list(set(list(y2['Receptor.ApprovedSymbol'])))))
# array(['DCHS1', 'IL18R1', 'PLXNC1', 'TNFRSF10D'], dtype=object) do not come up as secreted

## Narrow down to ligands actually present in adata_concat
secrLigandAC = np.intersect1d(list(set(list(y1['Ligand.ApprovedSymbol']))),list(adata_concat.var['gene_symbols']))
print("secreted ligands in adata_concat dataset:", secrLigandAC)
y3 = y1[y1['Ligand.ApprovedSymbol'].isin(secrLigandAC)]


#y3.to_csv('Tri_Unique_Receptor_secreted_Ligand_Pairs.csv')
#np.setdiff1d(np.array(Tri_GSC_DE_Receps_LS_2['symbol']), np.array(list(set(list(y3['Receptor.ApprovedSymbol'])))))
# ['DCHS1', 'IL18R1', 'IL1R1', 'PLA2R1', 'PLXNC1','TNFRSF10D'] are

newreceps = np.intersect1d(np.array(Tri_GSC_DE_Receps_LS_2['symbol']), np.array(list(set(list(y3['Receptor.ApprovedSymbol'])))))
len(newreceps) #27
## 6/8/21
# 'CD46', 'CELSR1', 'CMKLR1', 'DPP4', 'F3', 'FGFR4', 'FPR1', 'FXYD6', 'GPC1', 'IL20RA', 'IL6R', 'IL27RA', 'ITGA10', 'ITGA3', 'KDR', 'LDLR', 'LGR6', 'LRP8', 'MCAM', 'NOTCH3', 'PDGFRA', 'PLSCR4', 'PLXNA2','PTPRM', 'RAMP1', 'SCARB1', 'VIPR1'

# Plot all receptors with dotplot
sc.pl.dotplot(adata_concat, list(newreceps), gene_symbols = 'gene_symbols', groupby='batch_new_clusters', save = '_060821_dotplot_26_receps_new.pdf')

Tri_GSC_DE_Receps_LS_3 = Tri_GSC_DE_Receps_LS_2[Tri_GSC_DE_Receps_LS_2['symbol'].isin(newreceps)]
#['IL20RA', 'PDGFRA', 'LGR6', 'FPR1', 'LRP8', 'FGFR4', 'CMKLR1', 'F3', 'IL27RA', 'PLSCR4', 'KDR', 'NOTCH3', 'PLXNA2', 'ITGA10', 'FXYD6', 'RAMP1', 'GPC1', 'CELSR1', 'DPP4', 'PTPRM', 'VIPR1', 'SCARB1', 'LDLR', 'MCAM', 'CD46', 'ITGA3']
#info4 = info4.sort_values(by=['fold_change'], ascending=False)
#info4_w_entrez = info3_w_entrez[info3_w_entrez['gene_symbols'].isin(newreceps)]
Tri_GSC_DE_Receps_LS_4 = pd.merge(Tri_GSC_DE_Receps_LS_3, pd.DataFrame(y3['Receptor.ApprovedSymbol'].value_counts()).reset_index().rename(columns={"Receptor.ApprovedSymbol":"total_num_of_ligands"}), how = 'left', right_on = 'index', left_on = 'symbol').drop(columns = ['index'])
#info6 = info5[['gene_symbols', 'fold_change', 'total_num_of_ligands']]
print(Tri_GSC_DE_Receps_LS_4)
print("Receptors with secreted ligands in adata_concat:", list(Tri_GSC_DE_Receps_LS_4['symbol']))

### Investigate ligand expression in other cell types
secretedligand_list = list(set(list(y3['Ligand.ApprovedSymbol'])))
tmp1 = pd.merge(df1, pd.DataFrame(secretedligand_list, columns = ['secreted_ligand']), how='left', right_on='secreted_ligand', left_on = 'symbol')
secretedligands_to_investigate = tmp1[tmp1['secreted_ligand'].notna()].drop(columns = ['secreted_ligand'])
# 49 secreted ligands
secretedligands_to_investigate_ensembl_list = secretedligands_to_investigate['ensembl']

## Subset cell types into unique variables
adata_concat_Tri_ASTRO_HUVECS = adata_concat_Tri[adata_concat_Tri.obs['new_clusters']!='GSCS']
adata_concat_Tri_ASTRO = adata_concat_Tri_ASTRO_HUVECS[adata_concat_Tri_ASTRO_HUVECS.obs['new_clusters']!='HUVECS']
adata_concat_Tri_HUVECS = adata_concat_Tri_ASTRO_HUVECS[adata_concat_Tri_ASTRO_HUVECS.obs['new_clusters']!='ASTROCYTES']


###### Presence (>10% of cells) in cell type filter  - presence (not mean evaluated)
adata_concat_Tri_ASTRO_HUVECS_ligand = adata_concat_Tri_ASTRO_HUVECS[:,secretedligands_to_investigate_ensembl_list]
HUV_ASTRO_ligand = adata_concat_Tri_ASTRO_HUVECS_ligand.var_names[(np.sum(adata_concat_Tri_ASTRO_HUVECS_ligand.X.todense() > 0, axis=0)/(len(adata_concat_Tri_ASTRO_HUVECS_ligand.obs))>=0.1).tolist()[0]].tolist()
len(HUV_ASTRO_ligand) #39 ligands

adata_concat_Tri_ASTRO_ligand = adata_concat_Tri_ASTRO[:,secretedligands_to_investigate_ensembl_list]
ASTRO_ligand = adata_concat_Tri_ASTRO_ligand.var_names[(np.sum(adata_concat_Tri_ASTRO_ligand.X.todense() > 0, axis=0)/(len(adata_concat_Tri_ASTRO_ligand.obs))>=0.1).tolist()[0]].tolist()
len(ASTRO_ligand) # 36 ligands

adata_concat_Tri_HUVECS_ligand = adata_concat_Tri_HUVECS[:,secretedligands_to_investigate_ensembl_list]
HUVEC_ligand = adata_concat_Tri_HUVECS_ligand.var_names[(np.sum(adata_concat_Tri_HUVECS_ligand.X.todense() > 0, axis=0)/(len(adata_concat_Tri_HUVECS_ligand.obs))>=0.1).tolist()[0]].tolist()
len(HUVEC_ligand)#35 ligands

top_secreted_ligands = list(set(HUV_ASTRO_ligand + ASTRO_ligand + HUVEC_ligand))
len(top_secreted_ligands)# 41 ligands
tmp2 = pd.merge(df1, pd.DataFrame(top_secreted_ligands, columns = ['secreted_ligand']), how='left', right_on='secreted_ligand', left_on = 'ensembl')
top_secreted_ligands_2 = tmp2[tmp2['secreted_ligand'].notna()].drop(columns = ['secreted_ligand'])

top_pairs_2 = pd.merge(y3, top_secreted_ligands_2['symbol'], how = 'left', right_on = 'symbol', left_on = 'Ligand.ApprovedSymbol').dropna()
len(top_pairs_2) # 51
top_pairs_2.sort_values(by=['Receptor.ApprovedSymbol']).to_csv('121421_top_LR_pairs.csv')

Tri_GSC_DE_Receps_LS_5 = pd.merge(Tri_GSC_DE_Receps_LS_4, pd.DataFrame(top_pairs_2['Receptor.ApprovedSymbol'].value_counts()).reset_index().rename(columns={"Receptor.ApprovedSymbol":"num_of_ligands_present"}), how = 'left', right_on = 'index', left_on = 'symbol')
Tri_GSC_DE_Receps_LS_6 = Tri_GSC_DE_Receps_LS_5.dropna().drop(columns =['index'])
print(Tri_GSC_DE_Receps_LS_6)
print("Receptors with secreted ligands PRESENT in astrocytes, HUVECS, or both:", list(Tri_GSC_DE_Receps_LS_6['symbol']))
Tri_GSC_DE_Receps_LS_6.to_csv('Tri_Receptor_Ligand_pairs_secreted_and_present_with_info.csv')
len(Tri_GSC_DE_Receps_LS_6) #17 receptors that pass filter
#['PDGFRA', 'LGR6', 'FPR1', 'LRP8', 'FGFR4', 'F3', 'PLSCR4', 'KDR', 'NOTCH3', 'PLXNA2', 'GPC1', 'CELSR1','IL6R', 'VIPR1', 'SCARB1', 'LDLR', 'ITGA3']

sc.pl.dotplot(adata_concat, list(Tri_GSC_DE_Receps_LS_6['symbol']), gene_symbols = 'gene_symbols', groupby='batch_new_clusters', save = '_121421_17_receps_new.pdf')

for recep in list(set(list(top_pairs_2['Receptor.ApprovedSymbol']))):
    print(recep, list(top_pairs_2[top_pairs_2['Receptor.ApprovedSymbol']==recep]['symbol']))


############################################
## Idenitfy activated downstream pathways ##
## from receptor in GSCs.                 ##
############################################
"""
print('Integrate downstream pathways...')
################ Look at downstream pathways of upregulated receptors in Tri GSC cluster #########################
# Database info - put in dictionary
# gene sets as entrez IDs
database = {}
for path in [['c2.cp.pid.v7.2.entrez.gmt', 'pid'], ['c2.cp.wikipathways.v7.2.entrez.gmt', 'wiki'], ['c2.cp.biocarta.v7.2.entrez.gmt', 'biocarta'], ['c2.cp.kegg.v7.2.entrez.gmt', 'kegg']]:
    database[path[1]] = {}
    with open('pathways/'+path[0], 'r') as inFile:
        while 1:
            inLine = inFile.readline()
            if not inLine:
                break
            splitUp = inLine.strip().split('\t')
            Name = splitUp.pop(0)
            url = splitUp.pop(0)
            database[path[1]][Name] = splitUp

database_df = pd.DataFrame.from_dict(database, orient='index').melt()[~pd.isnull(pd.DataFrame.from_dict(database, orient='index').melt().value)].rename(columns={'variable':'pathway', 'value':'genes'})
database_df['database'] = [i.split('_')[0] for i in list(database_df['pathway'])]
database_df = database_df.set_index('database')
database_df_2 = database_df.explode('genes')
# merge with other gene IDs
database_df_3 = pd.merge(database_df_2, df1, how='left', right_on='entrezgene', left_on='genes').drop(columns=['entrez'])

# Receptors that have mapped pathways
tmp = pd.merge(database_df_2, Tri_GSC_DE_Receps_LS_6['entrezgene'], how = 'left', right_on = 'entrezgene', left_on = 'genes')
recep_pathways = tmp[tmp['entrezgene'].notna()].drop(columns = ['entrezgene'])
mapped_receptor_info = Tri_GSC_DE_Receps_LS_6[Tri_GSC_DE_Receps_LS_6['entrezgene'].isin(list(recep_pathways['genes'].unique()))]
mapped_receptor_info = mapped_receptor_info.rename(columns ={'genes':'ensembl'})
mapped_receptor_info_2 = pd.merge(mapped_receptor_info, recep_pathways, how = 'left', right_on = 'genes', left_on = 'entrezgene').drop(columns=['genes'])

print("Candidate receptors that have mapped pathway:", list(mapped_receptor_info['symbol']))
print("Mapped pathways associated with candidate receptors:", list(mapped_receptor_info_2['pathway'].unique()))
print("Candidate receptors that do not have mapped pathways:", list(np.setdiff1d(list(Tri_GSC_DE_Receps_LS_6['symbol']), list(mapped_receptor_info['symbol']))))


# Hypergeometric
pathways = list(mapped_receptor_info_2['pathway'].unique())
len(pathways)
# 128
path = []
M = len(df1['entrezgene'].dropna())
# total number of genes (with entrez IDs): 16481
N = len(Tri_GSC_DE_genes_FC1['entrezgene'].dropna())
# total number of differentially expressed genes in Triculture GSC cluster (FC>=1): 1359
geneNames = []
genesinpathway = []
genesinbigM = []
geneNames_DE = []
genesDE = []
hypergeometric = []
for pathway in pathways:
    path.append(pathway)
    path_subset = database_df_3[database_df_3['pathway']==pathway]
    tmp = list(path_subset['genes'])
    tmp_symbols = list(path_subset['symbol'])
    geneNames.append(",".join([str(i) for i in tmp_symbols]))
    n = len(tmp)
    genesinpathway.append(n)
    tmp2 = np.intersect1d(tmp,list(df1['entrezgene'].dropna()))
    n_new = len(tmp2)
    genesinbigM.append(n_new)
    tmp3 = list(np.intersect1d(tmp,Tri_GSC_DE_genes_FC1['entrezgene'].dropna()))
    if len(tmp3) == 0:
        geneNames_DE.append("NA")
    else:
        tmp3_symbols = list(path_subset[path_subset['genes'].isin(tmp3)]['symbol'])
        geneNames_DE.append(",".join([str(i) for i in tmp3_symbols]))
    k = len(tmp3)
    genesDE.append(k)
    tmp4 = hypergeom.sf(k,M,n_new,N,loc=0)
    hypergeometric.append(tmp4)

hyperinfo = pd.DataFrame({'path':path, 'geneNames': geneNames, 'num_genesinpathway':genesinpathway, 'num_genesinbigM':genesinbigM, 'geneNamesDE': geneNames_DE, 'num_genesDE':genesDE, 'hypergeometric':hypergeometric})
hyperinfo_2 = pd.merge(mapped_receptor_info_2, hyperinfo, how = 'left', left_on = 'pathway', right_on = 'path').rename(columns={'fold_change':"log_fold_change"})
hyperinfo_3 = hyperinfo_2[['ensembl', 'symbol', 'entrezgene', 'log_fold_change', 'p_val', 'pval_adj', 'total_num_of_ligands', 'num_of_ligands_present', 'pathway', 'geneNames','num_genesinpathway', 'num_genesinbigM', 'geneNamesDE', 'num_genesDE', 'hypergeometric']]
hyperinfo_4 = hyperinfo_3.sort_values(by =['hypergeometric'])
enriched_pathway_info = hyperinfo_4[hyperinfo_4['hypergeometric']<=0.05]
enriched_pathway_info.to_csv('061721_enriched_pathway_info.csv')

enriched_pathway_info_2 = enriched_pathway_info[['symbol', 'pathway']].sort_values(by='symbol').groupby(['symbol']).agg({'pathway': lambda x: x.tolist()})
len(enriched_pathway_info_2) #14
print(enriched_pathway_info_2)
print("Receptors with enriched pathways:", list(enriched_pathway_info['symbol'].unique()))

enriched_pathway_info_3 = enriched_pathway_info.sort_values(by='log_fold_change', ascending=False)

# stacked violin plot
genes1 = list(enriched_pathway_info_3['symbol'].unique())
#convGenes = [list(adata_concat.var.index[adata_concat.var['gene_symbols']==i])[0] for i in genes1]
for celltype in ['HUVECS', 'ASTROCYTES', 'GSCS']:
    sc.pl.stacked_violin(adata_concat[adata_concat.obs.new_clusters == celltype], genes1, gene_symbols = 'gene_symbols', groupby='batch', swap_axes=True, figsize = (2,15), save = '_gene_distribution_'+celltype+'.pdf')
"""


# ----------------------------------pathway gene sets as symbols ---------------------------------------#
print('Integrate downstream pathways...')
################ Look at downstream pathways of upregulated receptors in Tri GSC cluster #########################
# Database info - put in dictionary
database = {}
for path in [['c2.cp.pid.v7.4.symbols.gmt', 'pid'], ['c2.cp.wikipathways.v7.4.symbols.gmt', 'wiki'], ['c2.cp.biocarta.v7.4.symbols.gmt', 'biocarta'], ['c2.cp.kegg.v7.4.symbols.gmt', 'kegg']]:
    database[path[1]] = {}
    with open('pathways_V2/'+path[0], 'r') as inFile:
        while 1:
            inLine = inFile.readline()
            if not inLine:
                break
            splitUp = inLine.strip().split('\t')
            Name = splitUp.pop(0)
            url = splitUp.pop(0)
            database[path[1]][Name] = splitUp

database_df = pd.DataFrame.from_dict(database, orient='index').melt()[~pd.isnull(pd.DataFrame.from_dict(database, orient='index').melt().value)].rename(columns={'variable':'pathway', 'value':'genes'})
database_df['database'] = [i.split('_')[0] for i in list(database_df['pathway'])]
database_df = database_df.set_index('database')
database_df_2 = database_df.explode('genes')
# merge with other gene IDs
database_df_3 = pd.merge(database_df_2, df1, how='left', right_on='symbol', left_on='genes')

# Receptors that have mapped pathways
tmp = pd.merge(database_df_2, Tri_GSC_DE_Receps_LS_6['symbol'], how = 'left', right_on = 'symbol', left_on = 'genes')
recep_pathways = tmp[tmp['symbol'].notna()]
mapped_receptor_info = Tri_GSC_DE_Receps_LS_6[Tri_GSC_DE_Receps_LS_6['symbol'].isin(list(recep_pathways['genes'].unique()))]
# 15
mapped_receptor_info_2 = pd.merge(mapped_receptor_info, recep_pathways, how = 'left', right_on = 'genes', left_on = 'symbol')

print("Candidate receptors that have mapped pathway:", list(mapped_receptor_info['symbol']))
print("Mapped pathways associated with candidate receptors:", list(mapped_receptor_info_2['pathway'].unique()))
print("Candidate receptors that do not have mapped pathways:", list(np.setdiff1d(list(Tri_GSC_DE_Receps_LS_6['symbol']), list(mapped_receptor_info['symbol']))))


# Hypergeometric
pathways = list(mapped_receptor_info_2['pathway'].unique())
len(pathways)
# 128
M = len(df1['symbol'].dropna())
# total number of genes (with entrez IDs): 16481
# total number of genes (with symbols) : 17280
#N = len(Tri_GSC_DE_genes_FC1['symbol'].dropna())
# total number of differentially expressed genes in Triculture GSC cluster (FC>=1): 1359
# symbol: 1506
N = len(Tri_genes_with_cutoffs['symbol']) #525 DE genes
#N = len(Tri_genes_with_cutoffs_strict['symbol']) #DE genes with higher fold change cutoff

path = []
geneNames = []
genesinpathway = []
genesinbigM = []
geneNames_DE = []
genesDE = []
hypergeometric = []
for pathway in pathways:
    path.append(pathway)
    path_subset = database_df_3[database_df_3['pathway']==pathway]
    tmp = list(path_subset['genes'])
    #tmp_symbols = list(path_subset['symbol'])
    geneNames.append(",".join([str(i) for i in tmp]))
    n = len(tmp)
    genesinpathway.append(n)
    tmp2 = np.intersect1d(tmp,list(df1['symbol'].dropna()))
    n_new = len(tmp2)
    genesinbigM.append(n_new)
    #tmp3 = list(np.intersect1d(tmp,Tri_GSC_DE_genes_FC1['symbol'].dropna()))
    tmp3 = list(np.intersect1d(tmp,Tri_genes_with_cutoffs['symbol'].dropna()))
    #tmp3 = list(np.intersect1d(tmp,Tri_genes_with_cutoffs_strict['symbol'].dropna()))
    if len(tmp3) == 0:
        geneNames_DE.append("NA")
    else:
        #tmp3_symbols = list(path_subset[path_subset['genes'].isin(tmp3)]['symbol'])
        geneNames_DE.append(",".join([str(i) for i in tmp3]))
    k = len(tmp3)
    genesDE.append(k)
    tmp4 = hypergeom.sf(k,M,n_new,N,loc=0)
    hypergeometric.append(tmp4)

hyperinfo = pd.DataFrame({'path':path, 'geneNames': geneNames, 'num_genesinpathway':genesinpathway, 'num_genesinbigM':genesinbigM, 'geneNamesDE': geneNames_DE, 'num_genesDE':genesDE, 'hypergeometric':hypergeometric})

hyperinfo_2 = pd.merge(mapped_receptor_info_2, hyperinfo, how = 'left', left_on = 'pathway', right_on = 'path').rename(columns={'fold_change':"log2_fold_change"})
hyperinfo_3 = hyperinfo_2.sort_values(by =['hypergeometric'])
enriched_pathway_info = hyperinfo_3[hyperinfo_3['hypergeometric']<=0.05]
len(enriched_pathway_info['symbol_x'].unique())

enriched_pathway_info.sort_values(by=['log2_fold_change'], ascending=False).to_csv('121621_enriched_pathway_info_DEGs_match.csv') # 16


enriched_pathway_info.sort_values(by=['log2_fold_change'], ascending=False).to_csv('121521_enriched_pathway_info_DEGs_match_strict.csv') #15

enriched_pathway_info_2 = enriched_pathway_info[['symbol_x', 'pathway']].sort_values(by='symbol_x').groupby(['symbol_x']).agg({'pathway': lambda x: x.tolist()})
len(enriched_pathway_info_2) #14
print(enriched_pathway_info_2)
print("Receptors with enriched pathways:", list(enriched_pathway_info['symbol_x'].unique()))

enriched_pathway_info_3 = enriched_pathway_info.sort_values(by='log2_fold_change', ascending=False)

# stacked violin plot
genes1 = list(enriched_pathway_info_3['symbol_x'].unique())
#convGenes = [list(adata_concat.var.index[adata_concat.var['gene_symbols']==i])[0] for i in genes1]
for celltype in ['HUVECS', 'ASTROCYTES', 'GSCS']:
    sc.pl.stacked_violin(adata_concat[adata_concat.obs.new_clusters == celltype], genes1, gene_symbols = 'gene_symbols', groupby='batch', swap_axes=True, figsize = (2,15), save = '_121421_gene_distribution_'+celltype+'.pdf')
