##########################################################
## ccAF:  classifyPrimaryCells_gliomas.py               ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

# General
from os.path import exists
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from copy import deepcopy
from multiprocessing import Pool
import pickle

# Single Cell Packages
import scanpy as sc

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import classification_report

# Custom classes for classification
import classifiersV3 as cl

def confusionMatrix(a1, b1):
    v1 = set(a1)
    v2 = set(b1)
    dfTmp = pd.DataFrame([pd.Series(list(a1)),pd.Series(list(b1))]).T
    df1 = pd.DataFrame(index=v1, columns = v2)
    for i in v1:
        for j in v2:
            df1.loc[i,j] = dfTmp.loc[(dfTmp[0]==i) & (dfTmp[1]==j)].shape[0]
    return df1

# Load up classifier as ccAF1 -> Genes are Ensmbl gene IDs
with open('results/ACTINN/ccAF_2672_smaller.pkl','rb') as pklFile:
    ccAF1 = pickle.load(pklFile)

# HGNC -> downlaoded from HGNC website (https://www.genenames.org/download/custom/)
hgncEnsmbl = pd.read_csv('Whitfield/geneConversions/hgnc_geneSymbols_ensmbl.txt', index_col=1, header=0, sep='\t')
hgncEnsmbl = hgncEnsmbl.loc[~hgncEnsmbl['Ensembl ID(supplied by Ensembl)'].isnull()]

ensmblHgnc = pd.Series(hgncEnsmbl.index)
ensmblHgnc.index = list(hgncEnsmbl['Ensembl ID(supplied by Ensembl)'])

hgncPrevEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Previous symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Previous symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncPrevEnsmbl[j] = ensmbl

hgncAliasEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Alias symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Alias symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncAliasEnsmbl[j] = ensmbl

# Load up loom file
set1_scanpy = sc.read_loom('NSC/cellcycle_int_integrated_V3_6_18_2020.loom')
set1_scanpy = set1_scanpy[:,[True if i in hgncEnsmbl.index or i in hgncPrevEnsmbl or i in hgncAliasEnsmbl else False for i in set1_scanpy.var_names]]
# Convert to ensmbl
convertMe_hgnc = pd.DataFrame(np.nan, index=set1_scanpy.var_names, columns=['ensembl'])
numConverted_hgnc = 0
missed = []
for j in set1_scanpy.var_names:
    if j in hgncEnsmbl.index:
        convertMe_hgnc.loc[j,'ensembl'] = hgncEnsmbl.loc[j,'Ensembl ID(supplied by Ensembl)']
    elif j in hgncPrevEnsmbl:
        convertMe_hgnc.loc[j,'ensembl'] = hgncPrevEnsmbl[j]
    elif j in hgncAliasEnsmbl:
        convertMe_hgnc.loc[j,'ensembl'] = hgncAliasEnsmbl[j]
    else:
        missed.append(j)

set1_scanpy.var_names = pd.Index(convertMe_hgnc['ensembl'])

#############################
### ACTINN Classification ###
#############################
    # Cross validation within dataset
print('\nACTINN '+set1+'...')
set1PredLbls = ccAF1.predict_labels(set1_scanpy)
#set1PredProbs = ccAF1.predict_probs(set1_scanpy)
#set1PredProbs.columns = set1_scanpy.obs_names
set1_scanpy.obs['ccAF_new'] = set1PredLbls
# Dataframe of true labels, predictions, probabilities for all iterations
DF = pd.concat([set1_scanpy.obs, pd.DataFrame({'Predictions':set1PredLbls},index=set1_scanpy.obs_names)], axis=1)
DF.to_csv('results/ccAF_results_NSC.csv')
# Make Table 1 data
t1 = pd.DataFrame(pd.Series(set1PredLbls).value_counts()).T
t1.index = ['Counts']
t2 = confusionMatrix(set1_scanpy.obs['clusts_named'], set1PredLbls)
table1 = pd.concat([t1, t2[list(t1.columns)]])
table1.to_csv('results/table1_NSC.csv')
