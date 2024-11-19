
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom

import mygene
mg = mygene.MyGeneInfo()

# Database info - put in dictionary
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
database_df_3 = pd.merge(database_df_2, df1, how='left', right_on='entrezgene', left_on='genes').drop(columns=['entrezgene'])

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
print(enriched_pathway_info_2)
print("Receptors with enriched pathways:", list(enriched_pathway_info['symbol'].unique()))
