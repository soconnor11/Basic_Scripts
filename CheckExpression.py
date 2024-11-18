
# docker run -it -v '/home/soconnor/scGlioma_1:/files' cplaisier/scrna_seq_velocity_6_7_2021
#------------------------
# Import packages
#-----------------------

import pandas as pd
import numpy as np
import scanpy as sc
import os

#------------------------
# Setup section
#-----------------------

dirs = ['GSE84465', 'GSE131928', 'GSE70630', 'GSE89567']
datasets = ['GSE84465_all.loom', 'GSE131928_10X.loom', 'GSE70630.loom', 'GSE89567.loom']
sets = ['Darmanis', 'Neftel', 'Tirosh', 'Venteicher']

#receptors = ['FPR1', 'FPR2', 'SCARB1', 'LGR6']
receptors = ['PDGFRA', 'LGR6', 'FPR1', 'FGFR4', 'LRP8', 'F3']
#ligands = ['SAA1', 'RSPO3']
#both = receptors + ligands
#dict = {'Receptors' : receptors, 'Ligands' : ligands}


# Select dataset to investigate
n = 3
dir = dirs[n]
dataset = datasets[n]
set = sets[n]
data = sc.read_loom(dirs[n]+'/'+datasets[n])

# Find average expression of receptors of interest
avgExp_all = {}
perCell_all = {}
if n == 0:
    print(sets[n])
    #patient_info = ['BT_S1', 'BT_S2', 'BT_S4', 'BT_S6']
    data = data[data.obs['cell_type'] == 'Neoplastic']
    P1 = data[data.obs['patient'] == 'BT_S1'] #268
    P2 = data[data.obs['patient'] == 'BT_S2'] #531
    P3 = data[data.obs['patient'] == 'BT_S4'] #163
    P4 = data[data.obs['patient'] == 'BT_S6'] #129
    tags = ['P1', 'P2', 'P3', 'P4']
    samples = [P1, P2, P3, P4]
if n == 1:
    print(sets[n])
    P1 = data[data.obs['patient'] == 'MGH143'] #2182
    P2 = data[data.obs['patient'] == 'MGH124'] #1332
    P3 = data[data.obs['patient'] == 'MGH102'] #1237
    P4 = data[data.obs['patient'] == 'MGH125'] #901
    P5 = data[data.obs['patient'] == 'MGH115'] #661
    P6 = data[data.obs['patient'] == 'MGH126'] #152
    tags = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']
    samples = [P1, P2, P3, P4, P5, P6]
if n == 2:
    print(sets[n])
    P1 = data[data.obs['patient'] == 'MGH36'] #694
    P2 = data[data.obs['patient'] == 'MGH53'] #726
    P3 = data[data.obs['patient'] == 'MGH54'] #1174
    P4 = data[data.obs['patient'] == 'MGH60'] #430
    P5 = data[data.obs['patient'] == 'MGH93'] #439
    P6 = data[data.obs['patient'] == 'MGH97'] #584
    tags = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']
    samples = [P1, P2, P3, P4, P5, P6]
if n == 3:
    print(sets[n])
    P1 = data[data.obs['patient'] == 'MGH42'] #474
    P2 = data[data.obs['patient'] == 'MGH43'] #251
    P3 = data[data.obs['patient'] == 'MGH44'] #567
    P4 = data[data.obs['patient'] == 'MGH45'] #423
    P5 = data[data.obs['patient'] == 'MGH56'] #839
    P6 = data[data.obs['patient'] == 'MGH57'] #343
    P7 = data[data.obs['patient'] == 'MGH103'] #113
    tags = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
    samples = [P1, P2, P3, P4, P5, P6, P7]
    # Find expression for genes of interest
for gene1 in receptors:
    print(gene1)
    avgExp = {}
    perCell = {}
    genExp_all = {}
    for celltype in list(data.obs['patient'].unique()):
        # normalized expression data of genes of interest
        tmp1 = data[data.obs['patient'] == celltype][:,gene1].layers['norm_data'].todense()
        avgExp[celltype] = np.mean(tmp1)
        perCell[celltype] = int(tmp1.astype(bool).sum(axis=0))/len(tmp1)
        # normalized expression data of all genes
        tmp = pd.DataFrame(data[data.obs['patient'] == celltype].layers['norm_data'].todense())
        # remove full columns of 0
        tmp = tmp.loc[:, (tmp != 0).any(axis=0)]
        genExp = {}
        for metric in ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']:
            genExp[metric] = tmp.mean(axis=0).describe().loc[metric]
        genExp_all[celltype] = genExp
    avgExp_all[gene1] = avgExp
    perCell_all[gene1] = perCell
df = pd.DataFrame(avgExp_all)
print(df)
pd.DataFrame(avgExp_all).to_csv('for_GBMtriculture/'+set+'_receptor_expression.csv')
#pd.DataFrame(perCell_all).to_csv('for_GBMtriculture/'+set+'_receptor_percentage.csv')
pd.DataFrame(genExp_all).to_csv('for_GBMtriculture/'+set+'_general_expression_stats.csv')
for metric in ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']:
    # Combine dataset metric information with avg gene expression for each receptor of interest
    df[metric] = pd.DataFrame(genExp_all).loc[metric]
    df.to_csv('for_GBMtriculture/'+set+'_combined.csv')
passedFilt = {}
# Set filter cutoff
filter = '25%'
for patient in list(data.obs['patient'].unique()):
    #print(patient)
    tmp = df.loc[patient]
    d1 = pd.DataFrame(pd.DataFrame(tmp)[0:6][patient]>=tmp.loc[filter])
    passed = list(d1[d1].dropna().index)
    passedFilt[patient] = passed
pd.DataFrame.from_dict(passedFilt, orient='index')
pd.DataFrame.from_dict(passedFilt, orient='index').to_csv('for_GBMtriculture/'+set+'_receptors_greater_than_expression.csv')
# Generate dotplots of receptors / ligands of interest
os.chdir('for_GBMtriculture')
#for i in range(len(samples)):
    #sc.pl.dotplot(samples[i], receptors, groupby='patient', save=set+'_'+tags[i]+'.pdf')
sc.pl.dotplot(data, receptors, groupby='patient', save=set+'.pdf')
os.chdir("../")



#------------------------
# Find expression of genes of interest
#-----------------------

# mean expression
avgExp_all = {}
for gene1 in receptors:
    print(gene1)
    avgExp = {}
    genExp = {}
    for celltype in list(datas[dirs[n]].obs['patient'].unique()):
        #print(celltype)
        avgExp[celltype] = np.mean(datas[dirs[n]][datas[dirs[n]].obs['patient'] == celltype][:,gene1].X.todense())
        #avgExp[celltype] = pd.DataFrame(datas[dirs[n]][datas[dirs[n]].obs['patient'] == celltype][:,gene1].X.todense()).describe().loc['mean']
        genExp[celltype] = np.mean(datas[dirs[n]][datas[dirs[n]].obs['patient'] == celltype].X.todense())
        genExp[celltype] = pd.DataFrame(datas[dirs[n]][datas[dirs[n]].obs['patient'] == celltype].X.todense()).describe().loc['mean']
    avgExp_all[gene1] = avgExp

pd.DataFrame(avgExp_all).to_csv('for_GBMtriculture/'+set+'_receptor_expression.csv')













#for n in [0,1,2,3]:
#for n in [0,1,2,3,4,5]:
for n in [0,1,2,3,4,5,6]:
    sc.pl.dotplot(samples[n], dict, groupby='patient', save=tags[n]+'_'+set+'.pdf')

# mean expression
avgExp_all = {}
for gene1 in receptors:
    print(gene1)
    avgExp = {}
    for celltype in list(data.obs['patient'].unique()):
        #print(celltype)
        avgExp[celltype] = np.mean(data[data.obs['patient'] == celltype][:,gene1].X.todense())
    avgExp_all[gene1] = avgExp

pd.DataFrame(avgExp_all).to_csv('for_GBMtriculture/'+set+'_receptor_expression.csv')






# quartile / tertile score - legend
# red = 0, >0 is yellow; > median, green? only captures expression level - not how many cells have expression







receptors = ['FPR1', 'FPR2', 'LGR6', 'LRP8', 'PDGFRA', 'F3', 'FGFR4']
ligands = ['SAA1', 'RSPO3', 'PDGFA', 'PDGFB', 'PDGFC', 'PDGFD', 'IL6', 'APOE', 'LRPAP1']

receptors = ['TLR2', 'TLR4', 'SCARB1', 'AGER', 'P2RX7', 'CD36', 'SELS', 'GRM7', 'ADRA2A', 'MTNR1A']
receptors = ['TLR2', 'TLR4', 'SCARB1', 'AGER', 'P2RX7', 'CD36', 'GRM7', 'ADRA2A', 'MTNR1A']

receptors = ['FPR1', 'FPR2']
ligands = ['SAA1']
both = receptors + ligands
dict = {'Receptors' : receptors, 'Ligands' : ligands}

# Other SAA1 receptors / ligand interactions (literature)

new_receptors = ['SAA1', 'FPR1', 'FPR2', 'TLR2', 'TLR4', 'SCARB1', 'AGER', 'P2RX7', 'CD36', 'SELS']
new_receptors = ['TLR2', 'TLR4', 'SCARB1', 'AGER', 'P2RX7', 'CD36', 'SELS']

new_ligands = ['CXCL8'] #synergize with TLRs to enhance neutrophil recruitment through activation of FPR2
new_receptors = ['ENSG00000137462', 'ENSG00000136869', 'ENSG00000073060', 'ENSG00000204305', 'ENSG00000089041', 'ENSG00000135218', 'ENSG00000131871']
new_ligands = ['ENSG00000169429']
new_both = new_receptors + new_ligands
new_dict = {'New receptors' : new_receptors, 'New ligands' : new_ligands}

# SAA1 more receptors from CellTalk
new_receptors_2 = ['GRM7', 'ADRA2A', 'MTNR1A']

# Darmanis
sc.pl.stacked_violin(data, receptors, groupby='cell_type', swap_axes=True,  save = '.pdf')
sc.pl.dotplot(S6, dict, groupby='cell_type', save='_S6_.pdf')
sc.pl.dotplot(data, dict, groupby='patient', save='_V2.pdf')
sc.pl.dotplot(data, new_receptors, groupby='cell_type', save='_new.pdf')

data[:,both]
data[:,'FPR1']
data[data[:,'FPR1'].obs['cell_type'] == 'Neoplastic']

avgExp_all = {}
for gene1 in new_receptors:
    print(gene1)
    avgExp = {}
    for celltype in list(data.obs['cell_type'].unique()):
        #print(celltype)
        avgExp[celltype] = np.mean(data[data.obs['cell_type'] == celltype][:,gene1].X.todense())
    avgExp_all[gene1] = avgExp

pd.DataFrame(avgExp_all).to_csv('avg_lig_recep_exprs_celltypes_new_2.csv')


avgExp_all = {}
for gene1 in receptors:
    print(gene1)
    avgExp = {}
    for celltype in list(data.obs['patient'].unique()):
        #print(celltype)
        avgExp[celltype] = np.mean(data[data.obs['patient'] == celltype][:,gene1].X.todense())
    avgExp_all[gene1] = avgExp

pd.DataFrame(avgExp_all)

# Triculture
sc.pl.dotplot(adata_concat, new_dict, groupby='batch_new_clusters', save='_SAA1_new_receptors.pdf')


avgExp_all = {}
for gene1 in new_both:
    print(gene1)
    avgExp = {}
    for celltype in list(adata_concat.obs['batch_new_clusters'].unique()):
        #print(celltype)
        avgExp[celltype] = np.mean(adata_concat[adata_concat.obs['batch_new_clusters'] == celltype][:,gene1].X.todense())
    avgExp_all[gene1] = avgExp

pd.DataFrame(avgExp_all).to_csv('avg_lig_recep_exprs_celltypes.csv')
