
# Find enriched gene ontology terms for specific single cell clusters

docker run -it -v '/home/soconnor:/files' soconnor/scrna_seq_velocity_10062021
#------------------------
# Set up section / Load packages
#-----------------------

import numpy as np
import pandas as pd
import gseapy as gp
import os

#------------------------------------------------------
# Set up working directory / Organize files to load
#---------------------------------------------------

resdir = 'ccAF_newData/redo_analysis'
res1 = 0.8
for res1 in np.arange(0.1, 0.9, 0.1):
    res1 = str(res1)
    newpath = resdir+'/'+res1+'/enrichr'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    # Read in marker genes for each cluster
    mgenes = pd.read_csv(resdir+'/'+res1+'/'+res1+'_scTransform_Markers_together.csv')
    mgenes['cluster'] = mgenes['cluster'].astype(str)
    # Run enrichr
    enr_results = {}
    for clust1 in mgenes['cluster'].unique():
        enr = gp.enrichr(gene_list=mgenes[mgenes['cluster'] == clust1]['gene'], gene_sets=['GO_Biological_Process_2021'], organism='Human', description='test_name', outdir='test/enrichr_kegg', cutoff=0.05, no_plot=True)
        if 'Adjusted P-value' in enr.results.columns:
            enr_results[clust1] = enr.results.loc[enr.results['Adjusted P-value']<=0.05,:]
            enr_results[clust1].to_csv(newpath+'/'+clust1+'_enrichr_results.csv')
