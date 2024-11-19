##########################################################
## ccAF:  plottingClassifiers.py                        ##
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

# Plotting
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

## Load datasets and concatenate
SVMrej_CR = pd.read_csv('results/SVMrej/CV_classification_report_2672.csv')
NN_CR = pd.read_csv('results/ACTINN/CV_classification_report_2672.csv')
RF_CR = pd.read_csv('results/RFpy/CV_classification_report_2672.csv')
KNN_CR = pd.read_csv('results/KNN/CV_classification_report_2672.csv')

Classifiers_CR = pd.concat([SVMrej_CR, RF_CR, KNN_CR, NN_CR], axis = 0)
Classifiers_CR.rename(columns={'Unnamed: 0':'Cell Cycle State'}, inplace=True)

# Write out means
CR = Classifiers_CR.groupby(['Cell Cycle State','Classifier']).mean()
CR['Subset'] = 2672
CR.to_csv('results/cell_cycle_states_CV_stats.csv')

## Plot f1-score boxplot
# hue = cell cycle state
sns.set(style="darkgrid")
fig, ax = plt.subplots(figsize=(20,10))
sns.boxplot(x="Classifier", y="f1-score", hue = "Cell Cycle State", data=Classifiers_CR)
plt.savefig('results/F1_scores_for_each_state_per_classifier.pdf')
plt.clf()
# hue = classifier
sns.set(style="darkgrid")
fig, ax = plt.subplots(figsize=(20,10))
sns.boxplot(hue="Classifier", y="f1-score", x = "Cell Cycle State", data=Classifiers_CR)
plt.savefig('results/F1_scores_for_each_state_per_classifier_V2.pdf')
plt.clf()
