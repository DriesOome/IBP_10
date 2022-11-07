# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 10:36:47 2022

@author: janve
"""

def isGlycolationSite(seq,index):
    try: 
        if list(seq)[index+1] == 'N' and list(seq)[index+2] != 'P' and (list(seq)[index+3] == 'S' or list(seq)[index+3] == 'T'):
            return True
        else:
            return False
    except:
        pass        
    

import pandas as pd
data1 = pd.read_csv('dpPTMextended.tsv',sep = "\t")

data2 = pd.read_csv("compiledMultizDatabase.csv", names=["id", "seq"])
data2["organism"] = [item.split(".")[1].split("_")[1] for item in data2["id"]]
data2["id"] = [item.split(".")[0] for item in data2["id"]]
data1['conservation'] = 0

for i in range(0,len(data1["ensId"])):
    for j in range(0, len(data2['id'])):
        if data1['ensId'][i] == data2['id'][j] and data2['organism'][j] == 'hg38':
            print('next prot')
            index = data1['index'][i]
            seq = data2['seq'][j]
            print(list(seq)[index+1])
            if isGlycolationSite(seq, index) == True:
                data1['conservation'][i] += 1
        if data1['ensId'][i] == data2['id'][j]:
            seq = data2['seq'][j]
            if isGlycolationSite(seq, index):
                data1['conservation'][i] += 1
