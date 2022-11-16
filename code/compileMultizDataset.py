import os.path

import pandas as pd
import csv
import numpy as np
from readFasta import read_fasta

def compileMultizDataset():
    if not os.path.exists("../data/compiledMultizDatabase.csv"):
        print(">Compiling multiz dataset")
        csvFile = open("../data/multizDatabaseFull.csv", "w+")
        writer = csv.DictWriter(csvFile, fieldnames=["id","seq"])
        data = pd.DataFrame(columns = ["seq"])
        with open('../data/startingDatasets/knownCanonical.multiz100way.exonAA.fa') as fp:
            progress = 0
            for name, seq in read_fasta(fp):
                progress += 1
                idsplit = (name.split(' ')[0]).split('_')
                id = "_".join(idsplit[0:2])
                if not id in data.index:
                    data.loc[id] = [[-1]*(int(idsplit[3])-1)]
                    # print(data.loc[id])
                if not -1 in data.loc[id]["seq"]:
                    writer.writerow({"id": id[1:], "seq": ''.join(data.loc[id]["seq"])})
                    data.drop(id, axis=0,inplace=True)
                else:
                    try:
                        data.loc[id]["seq"][int(idsplit[2])-1] = seq
                    except:
                        print("index out of range")
                        data.loc[id]["seq"].append(seq)
                progress += 1
                print(str(progress) + " --- " + str(len(data)))
    else:
        print(">Multiz dataset already compiled")


