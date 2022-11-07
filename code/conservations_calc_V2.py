# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 10:36:47 2022

@author: janve
"""
from joblib import Parallel, delayed


def isGlycolationSite(seq,index):
    #print(seq[index-3:index+3])
    try:
        if list(seq)[index] == 'N' and list(seq)[index+1] != 'P' and (list(seq)[index+2] == 'S' or list(seq)[index+2] == 'T'):
            return True
        else:
            return False
    except:
        pass        
    

import pandas as pd
dbPTM = pd.read_csv('../data/dpPTMextended.tsv', sep ="\t")

multiz = pd.read_csv("../data/compiledMultizDatabase.csv", names=["id", "seq"])
multiz["organism"] = [item.split(".")[1].split("_")[1] for item in multiz["id"]]
multiz["id"] = [item.split(".")[0] for item in multiz["id"]]
dbPTM['conservation'] = 0

def calculateConservationForProtein(rowIndex, glycoIndex, subMultiz):
    conservationCounter = 0
    # print(subMultiz.loc[subMultiz["organism"]=="hg38"]["seq"].item()[glycoIndex-10:glycoIndex+10])
    if isGlycolationSite(subMultiz.loc[subMultiz["organism"]=="hg38"]["seq"].item(), glycoIndex-1) == True:
        conservationCounter += 1
        subMultiz = subMultiz.drop(subMultiz[subMultiz["organism"]=="hg38"].index, axis=0)
        # print(list(humanSequence)[index-10:index+10])
        for subMultizIndex, subMultizRow in subMultiz.iterrows():
            if isGlycolationSite(subMultizRow["seq"], glycoIndex-1) == True:
                conservationCounter += 1
    print(conservationCounter)
    return rowIndex, conservationCounter


results = Parallel(n_jobs=60, verbose=10)(
        delayed(calculateConservationForProtein)(rowIndex, row["index"], multiz[multiz["id"] == row["ensId"]].reset_index(drop=True))
        for rowIndex, row in dbPTM.iterrows()
)

for result in results:
    try:
        dbPTM.loc[result[0], "conservation"] = result[1]
    except:
        print("failed")
dbPTM.to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)


