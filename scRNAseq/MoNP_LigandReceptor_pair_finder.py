########################################################################
## MoNP ligand receptor pair finder across all cell types             ##
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

# docker run -it -v '/home/soconnor/old_home/MoNP/:/files' cplaisier/scrna_seq_velocity
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
resdir = 'analysis'
newpath = 'analysis_LR'
resdir2 = resdir+'/'+newpath
if not os.path.exists(resdir2):
    os.makedirs(resdir2)


# Load mouse MNN ligand receptor database
lr_network = pd.read_csv(resdir +'/mnn_lr_network.csv', index_col=0)

# Load merged data
data = sc.read_loom(resdir+'/MoNP_MoNP_VP_merged_022024.loom')

# Load in cluster markers
markers = {}
for cluster1 in range(0,10):
    tmp = pd.read_csv(resdir+'/cluster_'+str(cluster1)+'_marker_genes_MoNP_VP_vs_MoNP.csv', index_col =0)
    tmp.dropna(inplace=True)
    tmp = tmp[tmp['p_val_adj']<=0.05]
    if cluster1 == 0:
        cluster1 = 'VSMC/Fibro'
    if cluster1 == 1:
        cluster1 = 'M2'
    if cluster1 == 2:
        cluster1 = 'M1'
    if cluster1 == 3:
        cluster1 = 'Mo'
    if cluster1 == 4:
        cluster1 = 'EC'
    if cluster1 == 5:
        cluster1 = 'DC'
    if cluster1 == 6:
        cluster1 = 'M3'
    if cluster1 == 7:
        cluster1 = 'T cell'
    if cluster1 == 8:
        cluster1 = 'M4'
    if cluster1 == 9:
        cluster1 = 'Div'
    markers[cluster1] = tmp

# Separate by source
con = data[data.obs['source'] == 'MoNP']
vp = data[data.obs['source'] == 'MoNP_VP']

# Find LR infor in control data (MoNP)
# Initialize dataframe
pairs = pd.DataFrame(columns=['Ligand', 'Receptor', 'Sender', 'Receiver', 'Ligand_Counts', 'Ligand_Cells_Exp', 'Ligand_Avg_Exp', 'Ligand_Diff_Exp', 'Receptor_Counts', 'Receptor_Cells_Exp', 'Receptor_Avg_Exp', 'Receptor_Diff_Exp'])
# Set celltypes
cell_types = list(con.obs['integrated_cluster_named'].unique())
for pair_ind in range(len(lr_network)):
    lig1 = lr_network.iloc[pair_ind]['from']
    rec1 = lr_network.iloc[pair_ind]['to']
    if lig1 in data.var_names and rec1 in data.var_names:
        for celltype1 in cell_types:
            # Subset data to particular cell type
            receiver1 = celltype1
            subset = con[con.obs['integrated_cluster_named']==receiver1]
            # Find receptor details
            counts_rec = subset[:,rec1].layers['counts'].todense().sum()
            perc_cells_rec = round(np.count_nonzero(subset[:,rec1].layers['counts'].todense())/len(subset.obs), 4)
            avg_exp_rec = round(subset[:,rec1].X.todense().mean(),4)
            if len(markers[receiver1][markers[receiver1]['symbol'] == rec1]) != 0:
                diff_exp_rec = round(float(markers[receiver1][markers[receiver1]['symbol'] == rec1]['avg_log2FC']),4)
            else:
                diff_exp_rec = np.nan
            # Find matched ligand details in every cluster
            for celltype2 in cell_types:
                sender1 = celltype2
                subset2 = con[con.obs['integrated_cluster_named']==sender1]
                counts_lig = subset2[:,lig1].layers['counts'].todense().sum()
                perc_cells_lig = round(np.count_nonzero(subset2[:,lig1].layers['counts'].todense())/len(subset2.obs), 4)
                avg_exp_lig = round(subset2[:,lig1].X.todense().mean(),4)
                if len(markers[sender1][markers[sender1]['symbol'] == lig1]) != 0:
                    diff_exp_lig = round(float(markers[sender1][markers[sender1]['symbol'] == lig1]['avg_log2FC']),4)
                else:
                    diff_exp_lig = np.nan
                pairs = pairs.append({'Ligand': lig1, 'Receptor': rec1, 'Sender': sender1, 'Receiver': receiver1, 'Ligand_Counts': counts_lig, 'Ligand_Cells_Exp': perc_cells_lig, 'Ligand_Avg_Exp': avg_exp_lig, 'Ligand_Diff_Exp': diff_exp_lig, 'Receptor_Counts': counts_rec, 'Receptor_Cells_Exp': perc_cells_rec, 'Receptor_Avg_Exp': avg_exp_rec, 'Receptor_Diff_Exp': diff_exp_rec}, ignore_index=True)


pairs.to_csv(resdir+'/MoNP_LR_counts_and_expression_per_cluster_030124.csv')

df = pairs.dropna()
df2 = df[df['Receptor_Diff_Exp'] > 0]
df3 = df2[df2['Ligand_Diff_Exp'] > 0]
df4 = df3[df3['Receptor_Avg_Exp'] > 0]
df5 = df4[df4['Ligand_Avg_Exp'] > 0]
df6 = df5[df5['Ligand_Cells_Exp']>0.1]
df7 = df6[df6['Receptor_Cells_Exp']>0.1]

df7.to_csv(resdir+'/MoNP_ligand_and_receptor_positive_diff_exp.csv')


df = pairs[pairs['Receptor_Diff_Exp'].notna()]
df2 = df[df['Receptor_Diff_Exp'] > 0]

df4 = df3[df3['Receptor_Avg_Exp'] > 0]
df5 = df4[df4['Ligand_Avg_Exp'] > 0]
df6 = df5[df5['Ligand_Cells_Exp']>0.1]
df7 = df6[df6['Receptor_Cells_Exp']>0.1]


df8 = df[df['Receptor_Diff_Exp'] < 0]
df9 = df8[df8['Ligand_Diff_Exp'] < 0]
df10 = df9[df9['Receptor_Avg_Exp'] > 0]
df11 = df10[df10['Ligand_Avg_Exp'] > 0]
df12 = df11[df11['Ligand_Cells_Exp']>0.1]
df13 = df12[df12['Receptor_Cells_Exp']>0.1]

df13.to_csv(resdir+'/MoNP_ligand_and_receptor_negative_diff_exp.csv')




# Find LR info in VP data (MoNP_VP)
# Initialize dataframe
pairs2 = pd.DataFrame(columns=['Ligand', 'Receptor', 'Sender', 'Receiver', 'Ligand_Counts', 'Ligand_Cells_Exp', 'Ligand_Avg_Exp', 'Ligand_Diff_Exp', 'Receptor_Counts', 'Receptor_Cells_Exp', 'Receptor_Avg_Exp', 'Receptor_Diff_Exp'])
# Set celltypes
cell_types = list(vp.obs['integrated_cluster_named'].unique())
for pair_ind in range(len(lr_network)):
    lig1 = lr_network.iloc[pair_ind]['from']
    rec1 = lr_network.iloc[pair_ind]['to']
    if lig1 in data.var_names and rec1 in data.var_names:
        for celltype1 in cell_types:
            # Subset data to particular cell type
            receiver1 = celltype1
            subset = vp[vp.obs['integrated_cluster_named']==receiver1]
            # Find receptor details
            counts_rec = subset[:,rec1].layers['counts'].todense().sum()
            perc_cells_rec = round(np.count_nonzero(subset[:,rec1].layers['counts'].todense())/len(subset.obs), 4)
            avg_exp_rec = round(subset[:,rec1].X.todense().mean(),4)
            if len(markers[receiver1][markers[receiver1]['symbol'] == rec1]) != 0:
                diff_exp_rec = round(float(markers[receiver1][markers[receiver1]['symbol'] == rec1]['avg_log2FC']),4)
            else:
                diff_exp_rec = np.nan
            # Find matched ligand details in every cluster
            for celltype2 in cell_types:
                sender1 = celltype2
                subset2 = vp[vp.obs['integrated_cluster_named']==sender1]
                counts_lig = subset2[:,lig1].layers['counts'].todense().sum()
                perc_cells_lig = round(np.count_nonzero(subset2[:,lig1].layers['counts'].todense())/len(subset2.obs), 4)
                avg_exp_lig = round(subset2[:,lig1].X.todense().mean(),4)
                if len(markers[sender1][markers[sender1]['symbol'] == lig1]) != 0:
                    diff_exp_lig = round(float(markers[sender1][markers[sender1]['symbol'] == lig1]['avg_log2FC']),4)
                else:
                    diff_exp_lig = np.nan
                pairs2 = pairs2.append({'Ligand': lig1, 'Receptor': rec1, 'Sender': sender1, 'Receiver': receiver1, 'Ligand_Counts': counts_lig, 'Ligand_Cells_Exp': perc_cells_lig, 'Ligand_Avg_Exp': avg_exp_lig, 'Ligand_Diff_Exp': diff_exp_lig, 'Receptor_Counts': counts_rec, 'Receptor_Cells_Exp': perc_cells_rec, 'Receptor_Avg_Exp': avg_exp_rec, 'Receptor_Diff_Exp': diff_exp_rec}, ignore_index=True)

pairs2.to_csv(resdir+'/MoNP_VP_LR_counts_and_expression_per_cluster_030124.csv')


df = pairs2.dropna()
df2 = df[df['Receptor_Diff_Exp'] > 0]
df3 = df2[df2['Ligand_Diff_Exp'] > 0]
df4 = df3[df3['Receptor_Avg_Exp'] > 0]
df5 = df4[df4['Ligand_Avg_Exp'] > 0]
df6 = df5[df5['Ligand_Cells_Exp']>0.1]
df7 = df6[df6['Receptor_Cells_Exp']>0.1]

df7.to_csv(resdir+'/MoNP_VP_ligand_and_receptor_positive_diff_exp.csv')

df8 = df[df['Receptor_Diff_Exp'] < 0]
df9 = df8[df8['Ligand_Diff_Exp'] < 0]
df10 = df9[df9['Receptor_Avg_Exp'] > 0]
df11 = df10[df10['Ligand_Avg_Exp'] > 0]
df12 = df11[df11['Ligand_Cells_Exp']>0.1]
df13 = df12[df12['Receptor_Cells_Exp']>0.1]

df13.to_csv(resdir+'/MoNP_VP_ligand_and_receptor_negative_diff_exp.csv')




























# Old code
#---------------------------
# Load in integrated data
#---------------------------

#data1 = sc.read_loom(resdir + '/MoNP_MoNP_VP_integrated.loom')
# Load mouse MNN ligand receptor database
lr_network = pd.read_csv(resdir +'/mnn_lr_network.csv', index_col=0)
#pairs = pd.DataFrame(columns = range(11), index = range(11))
pairs = pd.DataFrame()
for i in range(11):
    print(i)
    # Load in differentially expressed genes for MoNP_VP vs. MoNP at each cluster
    mgenes = pd.read_csv(resdir +'/cluster_'+str(i)+'_marker_genes_MoNP_VP_vs_MoNP.csv', index_col=0)
    # Find upregulated receptors
    mgenes = mgenes[(mgenes['avg_log2FC']>=0.25) & (mgenes['p_val_adj']<=0.05)].dropna()
    enr_receptors = list(set(list(lr_network['to'])) & set(list(mgenes['symbol'])))
    tmp = lr_network[lr_network['to'].isin(enr_receptors)]
    ligands_to_check = list(tmp['from'])
    for j in range(11):
        if j!=i:
            print(j)
            other_mgenes = pd.read_csv(resdir +'/cluster_'+str(j)+'_marker_genes_MoNP_VP_vs_MoNP.csv', index_col=0)
            # Find upregulated ligands too
            other_mgenes = other_mgenes[(other_mgenes['avg_log2FC']>=0.25) & (other_mgenes['p_val_adj']<=0.05)].dropna()
            # Find just present ligands
            #other_mgenes = other_mgenes[(other_mgenes['avg_log2FC']>=0)].dropna()
            tmp2 = other_mgenes[other_mgenes['symbol'].isin(ligands_to_check)]
            tmp3 = tmp[tmp['from'].isin(list(tmp2['symbol']))]
            tmp3 = tmp3.rename(columns = {'from': 'ligand', 'to':'receptor'})
            tmp3['from'] = j
            tmp3['to'] = i
            pairs = pairs.append(tmp3, ignore_index=True)

pairs.to_csv(resdir2 + '/lr_pairs_both_upregulated_all_clusters_MoNP_VP_vs_MoNP.csv')

### presence filter - use expression data from integrated - > average expression (RNA)

# DOWNREGULATED (negative log2fc differentially expressed VP vs. no VP)


############################################
## Idenitfy activated downstream pathways ##
## from receptor in GSCs.                 ##
############################################

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
