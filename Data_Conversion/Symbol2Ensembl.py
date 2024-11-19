import numpy as np
import pandas as pd
import scanpy as sc
import os

# Load up conversion file
symEnsmbl = pd.read_csv('geneConversion/U5/filtered_feature_bc_matrix/features.tsv.gz', header = None, index_col=1, sep='\t')
tmp1 = pd.Series(symEnsmbl.index)
tmp1.loc[symEnsmbl.index.duplicated()] = [i+'.1' for i in symEnsmbl.loc[symEnsmbl.index.duplicated()].index]
symEnsmbl.index = pd.Index(tmp1)

# HGNC -> downlaoded from HGNC website (https://www.genenames.org/download/custom/)
hgncEnsmbl = pd.read_csv('geneConversion/Whitfield/hgnc_geneSymbols_ensmbl.txt', index_col=1, header=0, sep='\t')
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

# Load data
datasets = {}
for set1 in sets:
    # Load in normalized U5-hNSC data
    print('\nLoading U5-hNSC scRNA-seq data...')
    datasets[set1] = sc.read_loom(resdir+'/'+set1+'_normalized_gene_symbols.loom')
    datasets[set1].obs['dataset'] = set1
    # Convert genes to Ensembl
    print('\nConverting gene symbols to Ensembl...')
    common = [i for i in datasets[set1].var_names if i in symEnsmbl[0]]
    datasets[set1]= datasets[set1][:,common]
    datasets[set1].var_names = pd.Index(symEnsmbl.loc[datasets[set1].var_names,0], name='Ensembl')
