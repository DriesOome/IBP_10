"""
Created on Mon Nov 21 15:10:10 2022

@author: janve
"""

import pandas as pd
import re
from joblib import Parallel, delayed
from seqfold import fold
from math import floor

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

def dotBracketlocal(seq,splitwindow):
    seqlist= []
    n=0
    output = ""
    #print(len(seq)/splitwindow)
    for i in range(floor(len(seq)/splitwindow)): 
        seqlist.append(seq[n:n+splitwindow])
        #print(seq[n:n+splitwindow])
        n+=splitwindow
    for seqpart in seqlist:
        structs = fold(seqpart)
        desc = ["."] * len(seqpart)
        for s in structs:
            if len(s.ij) == 1:
                i, j = s.ij[0]
                desc[i] = "("
                desc[j] = ")"
        #print("".join(desc))
        output = output + "".join(desc)
    #print(output)
    return output

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

def foldedness(seq, window,splitwindow):
    #print(seq)
    result = dotBracketlocal(seq,splitwindow) 
    # print(result)
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
    # print(deriv)
    return deriv

def foldednessForProtein(rowIndex, row):
    try:
        seq = row["cdna"]
        if row["ensId"] in memoData:
            fold = memoData["ensId"]
        else:
            fold = foldedness(seq, window, splitwindow)
        return rowIndex, fold
    except Exception as e:
        print(e)

def extractFoldedness():
    dbPTM = pd.read_csv('../data/dpPTMextended.tsv', sep="\t")
    dbPTM["foldedness"] = 0
    results = Parallel(n_jobs=-5, verbose=5)(
        delayed(foldednessForProtein)(rowIndex, row)
        for rowIndex, row in dbPTM.iterrows())
    for result in results:
        try:
            dbPTM.loc[result[0], "foldedness"] = str(result[1])
        except Exception as e:
            print(e)
            print("* Failed for entry: ")

    dbPTM.to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)

window = 9
splitwindow = 50
memoData = {}

# extractFoldedness()

    