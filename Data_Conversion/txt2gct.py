

docker run -it -v '/home/soconnor/old_home/GBM_tumors/:/files' cplaisier/scrna_seq_velocity_6_7_2021

import pandas as pd
import os

# Load txt file and convert to gct file
tags = ["BT322", "BT333", "BT368"]
for tag in tags:
    resdir = tag+'/analysis_output'
    file = tag+".txt"
    txt = pd.read_csv(resdir+'/'+file, sep="\t", header=0, index_col=0)
    txt.insert(0,"Description",["Blah"]*txt.shape[0])
    gct = []
    gct.append("#1.2")
    gct.append(str(txt.shape[0])+"\t"+str(txt.shape[1]))
    gct.append("Name")
    #gct.append(txt.to_string())
    os.chdir(resdir)
    with open(tag+".gct", "w") as gctFile:
        gctFile.write("\n".join(gct))
        txt.to_csv(gctFile,sep="\t")
    os.chdir("../..")






"""
resdir = "ccAF_improve/gsea"
# Load txt file
#file = "U5_hNSC_expression_data.txt"
file = "integrated_expression_data.txt"
txt = pd.read_csv(resdir+"/"+file, sep="\t", header=0, index_col=0)
txt.insert(0,"Description",["Blah"]*txt.shape[0])

gct = []
gct.append("#1.3")
gct.append(str(txt.shape[0])+"\t"+str(txt.shape[1]-1)+"\t1\t0")
gct.append("Name")
with open(resdir+"/integrated_expression_data.gct", "w") as gctFile:
    gctFile.write("\n".join(gct))
    txt.to_csv(gctFile,sep="\t")
"""
