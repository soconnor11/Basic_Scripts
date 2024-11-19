

#------------------------------------------------------
# RNA velocity analysis
#---------------------------------------------------

docker run -it -v '/home/soconnor/U5_hNSC_Neural_G0:/files' cplaisier/scrna_seq_velocity_6_7_2021
# Install for data proessing
pip3 install -U umap-learn==0.3.10
#pip3 install tqdm
#pip3 install ipywidgets

#------------------------------------------------------
# Setup
#---------------------------------------------------

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

#------------------------------------------------------
# Load files into scanpy

tag = 'hNSC'
resdir = 'scVelo'
# Files need: velocyto loom file, loom file from Seurat completed analysis
# Load data from Seurat

adata_Seurat = scv.read_loom('data/cellcycle_int_integrated.loom')
adata_Seurat.shape #(6171, 19216)
adata_Seurat = adata_Seurat[adata_Seurat.obs['orig_ident']=='WT']

ccAF1_umap = pd.read_csv(resdir+"/"+tag+'_UMAPCoordinates.csv', header = 0, index_col = 0)
#ccAF1_umap = ccAF1_umap.loc[[True if not i.find('-1_1')==-1 else False for i in ccAF1_umap.index]]
#ccAF1_umap.index = [i.rstrip('-1_1') for i in list(ccAF1_umap.index)]
adata_Seurat.obsm['umap_cell_embeddings'] = np.array(ccAF1_umap.loc[adata_Seurat.obs_names])

# Load U5 data from velocyto
adata_vc = scv.read_loom('data/U5_hNSC/WT/U5_velocyto.loom')
adata_vc.shape


#------------------------------------------------------
# Merge data
#---------------------------------------------------

adata = scv.utils.merge(adata_Seurat, adata_vc)
adata.shape #(3015, 4567)

#adata.obs["seurat_new_clusters"]=adata.obs["seurat_clusters"]-1
#adata.obs.seurat_new_clusters = adata.obs.seurat_new_clusters.astype('category')
#adata.obs.seurat_new_clusters = adata.obs.seurat_new_clusters.rename_categories({'0':'Late G1', '1':'M', '2':'M/EarlyG1', '3':'G1', '4':'S/G2', '5':'G2', '6':'G1/other', '7':'G0', '8':'S'})
#adata.obs.seurat_new_clusters= adata.obs.seurat_new_clusters.cat.reorder_categories(['S', 'S/G2', 'G2', 'M', 'M/EarlyG1', 'Late G1', 'G0', 'G1/other', 'G1'], ordered=True)
#adata.rename_categories('seurat_new_clusters', new_cluster_names)

# Colorize by UMAP colors
#ident_colors = pd.read_csv(resdir+"/"+tag+"_seurat_umap_colors.csv")
#adata.uns['ClusterName_colors']=list(ident_colors["ucols"])
​
#------------------------------------------------------
# Preprocess the data
#---------------------------------------------------

#scv.pp.filter_genes(adata, min_shared_counts=10)
#scv.pp.normalize_per_cell(adata)
#scv.pp.filter_genes_dispersion(adata, n_top_genes=3000)
#scv.pp.log1p(adata)
#scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=3000)
#scv.pp.moments(adata, n_pcs=19, n_neighbors=10, use_rep='umap_cell_embeddings')
#scv.pp.pca(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=3000)
scv.pp.moments(adata) # can play with how many cells used to calculate
​
#------------------------------------------------------
# Estimate RNA velocity
#---------------------------------------------------

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.latent_time(adata)
adata.obs['latent_time'].to_csv('data/U5_hNSC/WT/latent_time.csv')
