"""
Created on Mon Nov 21 15:10:10 2022

@author: janve
"""

import pandas as pd
import re
from joblib import Parallel, delayed
from seqfold import fold

def dotBracket(seq):
    structs = fold(seq)
    desc = ["."] * len(seq)
    for s in structs:
        if len(s.ij) == 1:
            i, j = s.ij[0]
            desc[i] = "("
            desc[j] = ")"
    print("".join(desc))
    return("".join(desc))

def calculatederivative(y,window):
    derivative =  []
    if window%2 == 0:
        window = window/2
    else:
        window = (window-1)/2
    window = int(window)
    for i in range(0,window):
        derivative.append(abs(y[0]-y[i+window]))
    for i in range(0+window,len(y)-window):
        derivative.append(abs(y[i-window]-y[i+window]))
    for i in range(len(y)-window,len(y)):
        derivative.append(abs(y[i-window]-y[len(y)-1]))
    return derivative

def foldedness(seq, window):
    print(seq)
    result = dotBracket(seq) 
    print(result + '  1')
    split_result = re.split('',result)
    foldedness = [0]*len(split_result)
    for i in range(0,len(split_result)):
        if split_result[i] == '(':
            try:  
                foldedness[i] = foldedness[i-1]+1
            except: 
                foldedness[i] = 1
        if split_result[i] == '.':
            try:  
                foldedness[i] = foldedness[i-1]+0
            except: 
                foldedness[i] = 0
        if split_result[i] == ')':
                foldedness[i] = foldedness[i-1]-1
    deriv =  calculatederivative(foldedness, window)
    return deriv


dbPTM = pd.read_csv('dpPTMextended.tsv',sep = "\t")
window = 9

results = Parallel(n_jobs= None, verbose=5)(delayed(foldedness)(seq,window) for seq in dbPTM["cdna"])
for result in results:
    try:
        dbPTM.loc[result[0], "foldedness"] = result[1]
    except:
        print("* Failed for entry: "+result[0])

dbPTM.to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)

    