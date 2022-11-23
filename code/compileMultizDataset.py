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
        data = {}
        with open('../data/startingDatasets/knownCanonical.multiz100way.exonAA.fa') as fp:
            progress = 0
            fastaDict = read_fasta(fp)
            for name in fastaDict:
                progress += 1
                idsplit = (name.split(' ')[0]).split('_')
                id = "_".join(idsplit[0:2])
                subSequenceIndex = int(idsplit[2])
                totalSubSequences = int(idsplit[3])
                if not id in data:
                    data[id] = [-1]*(totalSubSequences)
                    print(">New entry made for: "+id)
                try:
                    data[id][subSequenceIndex-1] = fastaDict[name]
                except Exception as e:
                    print(e)
                if not -1 in data[id]:
                    writer.writerow({"id": id[1:], "seq": ''.join(data[id])})
                    print(">Completed for: "+id)
                    del data[id]
                progress += 1
                print(str(progress/len(fastaDict)/2))
    else:
        print(">Multiz dataset already compiled")


