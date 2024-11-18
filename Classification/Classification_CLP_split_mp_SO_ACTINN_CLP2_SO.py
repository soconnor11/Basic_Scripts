##########################################################
## OncoMerge:  classifiers.py                           ##
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

# ssh soconnor@plaisierlab1.sbhse.dhcp.asu.edu
# docker run -it -v '/home/soconnor/ccAF:/files' cplaisier/scrna_seq_velocity
# cd /files
# ls
# python3

# General
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
import pickle

# Single Cell Packages
import scanpy as sc

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.preprocessing import scale, robust_scale, quantile_transform
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim

# Custom classes for classification
import classifiersV3 as cl

# Classifier loading and writing
#import torch

# Plotting
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Function to translate predictions
# exp1, whit1 = compute_whitfield_experimental(whitfield_pred, experimental)
def compute_whitfield_experimental(whitfield_pred, conversion = {'M/Early G1':'M/Early G1','G2':'G2/M','S':'S','S/G2':'S/G2','Late G1':'Late G1','G1':'G1', 'M':'G2/M', 'G0':'Neural G0', 'G1/other':'G1/other'}):
    #exp1 = list(experimental.dropna())
    #whit1 = list(whitfield_pred[experimental.dropna().index])
    whit2 = [conversion[i] for i in whitfield_pred]
    return whit2

#compute_whitfield_experimental(whitfield.obs['randomForest'], meta['Cycle'])
def compute_whitfield_translation(whitfield_pred, conversion = {'M/Early G1':'M/G1','G2/M':'G2/M','S':'S','S/G2':'G2','Late G1':'S','G1':'G1/S', 'Unknown':'Unknown', 'Neural G0':'G1/S', 'G1/other':'Unknown'}):
    # Subset to only those time points with experimental cell cycle state data
    #whit1 = pd.Series(list(whitfield_pred[meta['Cycle'].dropna().index]))
    #
    ## Convert predictions to ['M','S']
    whit2 = [conversion[i] for i in whitfield_pred]
    #
    # Compute and return accuracy
    return whit2

# Load up classifier as ccAF1 -> Genes are Ensmbl gene IDs
with open('results/ACTINN/ccAF_2672_smaller.pkl','rb') as pklFile:
    ccAF1 = pickle.load(pklFile)

g2e = pd.read_csv('Whitfield/geneConversions/ensembl_entrez.csv', index_col=1, header=0)
g2e = g2e.loc[g2e.index.dropna()]
g2e = g2e.loc[~g2e.index.duplicated()]

for whit1 in [['whitfield_dataPlusScores_6_16_2020_', '.csv', 'no_manip'], ['whitfield_dataPlusScores_pow10_6_16_2020_', '.csv', 'pow10'], ['whitfield_dataPlusScores_pow2_6_22_2020_', '.csv', 'pow2'], ['whitfield_dataPlusScores_FC_6_29_2020_', '.csv', 'FC'], ['whitfield_dataPlusScores_6_30_2020_', '_1334.csv', '1334'], ['whitfield_dataPlusScores_pow10_6_30_2020_', '_1334.csv', 'pow10_1334'], ['whitfield_dataPlusScores_pow2_6_30_2020_', '_1334.csv', 'pow2_1334'], ['whitfield_dataPlusScores_FC_6_30_2020_', '_1334.csv', 'FC_1334']]:
    for mod1 in ['none', 'scale', 'robust', 'quantile']:
        # Load up whitfield dataset & incorporate ideal vector data from Figure 1 to meta
        experiments = {'TT1':[0,12], 'TT2':[12,38], 'TT3':[38,86], 'TN':[86,105], 'SHAKE':[105,114]}
        whitfield = {}
        print(whit1, mod1)
        for exp1 in experiments:
            whitfield[exp1] = sc.read_csv('Whitfield/'+whit1[0]+exp1+whit1[1], first_column_names=True).T
            var_names = [str(g2e.loc[float(i),'Gene stable ID']) for i in whitfield[exp1].var_names if float(i) in g2e.index]
            whitfield[exp1] = whitfield[exp1][:,[True if float(i) in g2e.index else False for i in whitfield[exp1].var_names]]
            whitfield[exp1].var_names = pd.Index(var_names)
            if mod1 == 'scale':
                tmp = whitfield[exp1].X
                whitfield[exp1].X = scale(tmp, axis = 1)
            if mod1 == 'robust':
                tmp = whitfield[exp1].X
                whitfield[exp1].X = robust_scale(tmp, axis = 1)
            if mod1 == 'quantile':
                tmp = whitfield[exp1].X
                whitfield[exp1].X = quantile_transform(tmp, axis = 1)
            whitfield[exp1].var_names = [i.rstrip('.0') for i in whitfield[exp1].var_names]
            whitfield[exp1].var_names_make_unique()

        whitfield_mg = pd.read_csv('Whitfield/markergenes_ForPlotting.csv', header=0, index_col=0, skiprows=[1,2,3])
        whitfield_phase = whitfield_mg['Phase']
        whitfield_mg_1 = whitfield_mg.drop('Phase', axis=1)
        whitfield_vectors = whitfield_mg_1.iloc[20:25]
        maxWhitfield_vector = [i.replace(' phase Vector','').replace(' vector','').replace(' Vector','') for i in whitfield_vectors.idxmax()]
        maxWhitfield_vector = [k.replace('_', '/') for k in maxWhitfield_vector]
        maxWhitfield_vector_2 = whitfield_mg.iloc[0:20].groupby('Phase').median().idxmax()
        maxWhitfield_vector_2 = [k.replace('_', '/') for k in maxWhitfield_vector_2]

        # Experimental and exprs3
        # Load up experimental meta data
        meta = pd.read_csv('Whitfield/metaInformation.csv', header=0, index_col=0)
        experimental = meta['Cycle']
        maxWhitfield_vector_3 = [maxWhitfield_vector_2[i] if not isinstance(meta['Cycle'].iloc[i], str) else {'M':'G2/M','S':'S'}[str(meta['Cycle'].iloc[i])] for i in range(len(list(maxWhitfield_vector_2)))]

        # Reconstitute output
        testPredLbls = {'NN':[]}

        #############
        ### NN CV ###
        #############

        # Classification with ground truth dataset
        for exp1 in ['TT1','TT2','TT3','TN','SHAKE']:
            testPredLbls['NN'] += list(ccAF1.predict_labels(whitfield[exp1]))


        # Set true labels
        truelab = maxWhitfield_vector_3
        # truelab = whitfield.obs["experimental"]

        # Compute scores for each classifier method
        pred = compute_whitfield_translation(testPredLbls['NN'])
        pred_exp = compute_whitfield_experimental(testPredLbls['NN'])

        # Dataframe of true labels, predictions, probabilities for all iterations
        DF = pd.DataFrame({'Experimental': meta['Cycle'],'maxWhitfield_mean':maxWhitfield_vector,'maxWhitfield_median': maxWhitfield_vector_2, 'maxWhitfield_median_exp': maxWhitfield_vector_3, 'Translated_Predictions_experimental':pred_exp, 'Translated_Predictions':pred, 'Predictions':testPredLbls['NN']})
        DF.to_csv('method/NN/whitfield/ccAF_CV_results_Whitfield_'+whit1[2]+'_pad_'+mod1+'.csv')

        # Get classification report for each iteration
        performanceResults = classification_report(truelab, pred, output_dict=True, zero_division=0)

        # Convert into a dataframe
        performDF = pd.DataFrame(performanceResults).T
        states1 = ['G1/S','G2', 'G2/M', 'M/G1', 'S']
        performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
        performDF['Classifier'] = 'NN'
        performDF.to_csv('method/NN/whitfield/CV_classification_report_Whitfield_'+whit1[2]+'_pad_'+mod1+'.csv')
