
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

receptors = ['FPR1', 'FPR2', 'SCARB1', 'TLR4', 'LGR6']
ligands = ['SAA1', 'RSPO3']
both = receptors + ligands
dict = {'Receptors' : receptors, 'Ligands' : ligands}

receptors = ['SCARB1']

# Select dataset to investigate
n = 3
dir = dirs[n]
dataset = datasets[n]
set = sets[n]
data = sc.read_loom(dirs[n]+'/'+datasets[n])

# Find average expression of receptors of interest
if n == 0:
    print(sets[n])
    #patient_info = ['BT_S1', 'BT_S2', 'BT_S4', 'BT_S6']
    data = data[data.obs['cell_type'] == 'Neoplastic']
    P1 = data[data.obs['patient'] == 'BT_S1']
    P2 = data[data.obs['patient'] == 'BT_S2']
    P3 = data[data.obs['patient'] == 'BT_S4']
    P4 = data[data.obs['patient'] == 'BT_S6']
    tags = ['P1', 'P2', 'P3', 'P4']
    samples = [P1, P2, P3, P4]
if n == 1:
    print(sets[n])
    P1 = data[data.obs['patient'] == 'MGH143']
    P2 = data[data.obs['patient'] == 'MGH124']
    P3 = data[data.obs['patient'] == 'MGH102']
    P4 = data[data.obs['patient'] == 'MGH125']
    P5 = data[data.obs['patient'] == 'MGH115']
    P6 = data[data.obs['patient'] == 'MGH126']
    tags = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']
    samples = [P1, P2, P3, P4, P5, P6]
if n == 2:
    print(sets[n])
    P1 = data[data.obs['patient'] == 'MGH36']
    P2 = data[data.obs['patient'] == 'MGH53']
    P3 = data[data.obs['patient'] == 'MGH54']
    P4 = data[data.obs['patient'] == 'MGH60']
    P5 = data[data.obs['patient'] == 'MGH93']
    P6 = data[data.obs['patient'] == 'MGH97']
    tags = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']
    samples = [P1, P2, P3, P4, P5, P6]
if n == 3:
    print(sets[n])
    P1 = data[data.obs['patient'] == 'MGH42']
    P2 = data[data.obs['patient'] == 'MGH43']
    P3 = data[data.obs['patient'] == 'MGH44']
    P4 = data[data.obs['patient'] == 'MGH45']
    P5 = data[data.obs['patient'] == 'MGH56']
    P6 = data[data.obs['patient'] == 'MGH57']
    P7 = data[data.obs['patient'] == 'MGH103']
    tags = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
    samples = [P1, P2, P3, P4, P5, P6, P7]
    # Find expression for genes of interest
avgExp_all = {}
perCell_all = {}
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
#filter = '25%'
filter = '50%'
for patient in list(data.obs['patient'].unique()):
    #print(patient)
    tmp = df.loc[patient]
    #d1 = pd.DataFrame(pd.DataFrame(tmp)[0:5][patient]>=tmp.loc[filter])
    d1 = pd.DataFrame(pd.DataFrame(tmp)[0:1][patient]>=tmp.loc[filter])
    passed = list(d1[d1].dropna().index)
    passedFilt[patient] = passed
print(pd.DataFrame.from_dict(passedFilt, orient='index'))
pd.DataFrame.from_dict(passedFilt, orient='index').to_csv('for_GBMtriculture/'+set+'_receptors_greater_than_expression.csv')
# Generate dotplots of receptors / ligands of interest
os.chdir('for_GBMtriculture')
for i in range(len(samples)):
    sc.pl.dotplot(samples[i], dict, groupby='patient', save=set+'_'+tags[i]+'.pdf')
os.chdir("../")
