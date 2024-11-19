
#------------------------------------------------------
# ccAF application
#-----------------------------------------------------

# for KAT5 Project
# Samantha O'Connor
# Updated 03/17/22
# ccAF description and build information: O'Connor et al., 2021
# Example: GSC827_inVitro_combined


#------------------------------------------------------
# Setup
#----------------------------------------------------

# File needed before begin: pre-normalized scRNA-seq data with genes as ensembl IDs as LOOM file.

#docker pull cplaisier/scrna_seq_velocity # if haven't already pulled
# Run the docker container using the following command (replace with home directory / parent folder of folder where scRNA-seq loom file is):
#docker run -it -v '/home/soconnor/old_home/GSC_bam/:/files' cplaisier/scrna_seq_velocity
#docker run -it -v '/media/omics1/ccNN/usftp21.novogene.com/01.RawData/:/files' cplaisier/scrna_seq_velocity
#docker run -it -v '/home/soconnor/old_home/GBM_tumors/:/files' cplaisier/scrna_seq_velocity
docker run -it -v '/home/soconnor/old_home/GSE159929_RAW/:/files' cplaisier/scrna_seq_velocity
cd /files
python3

# Imports
import pandas as pd
import scanpy as sc
import ccAF
import mygene
mg = mygene.MyGeneInfo()

#------------------------------------------------------
# Read in file and apply ccAF
#---------------------------------------------------

tags = ['Bladder', 'Blood', 'Common_bile_duct', 'Esophagus', 'Heart', 'Liver', 'Lymph_node', 'Marrow', 'Muscle', 'Rectum', 'Skin', 'Small_intestine', 'Spleen', 'Stomach', 'Trachea']
tag = tags[0]
dir1 = tag+'/analysis_output'

# Load in data
data = sc.read_loom(dir1+'/'+tag+'_data.loom')
#data = sc.read_loom(dir1+'/'+tag+'_data_2.loom') # BT324
#data = sc.read_loom(dir1+'/'+tag+'_data_3.loom') # BT333
data
data.var_names_make_unique()
data.obs_names_make_unique()

# Conver gene symbols to ensembl IDs if necessary
tmp = mg.querymany(data.var_names, scopes='symbol', fields='ensembl.gene', species = 'human', as_dataframe=True)
tmp2 = pd.DataFrame(tmp[tmp['notfound']!=True]['ensembl.gene']).dropna()
tmp3 = tmp2.reset_index().drop_duplicates(subset=['query']).set_index('query')
data_subset = data[:,data.var_names.isin(tmp3.index)]
data_subset.var['symbol'] = data_subset.var.index
data_subset.var['ensembl'] = tmp3['ensembl.gene']
data_subset.var.index = data_subset.var['ensembl']


# Apply ccAF
data_subset.obs['ccAF'] = ccAF.ccAF.predict_labels(data_subset)
data_subset.obs['ccAF'].value_counts()
pd.DataFrame(data_subset.obs['ccAF']).to_csv(dir1+'/'+tag+'_ccAF_calls.csv')
