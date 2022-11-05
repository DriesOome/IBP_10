import sys

import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from readFasta import read_fasta

def extractCDNA():
    dbPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
    dbPTM["cdna"] = np.nan
    with open('../data/startingDatasets/Homo_sapiens.cdna.all.fa') as fp:
        progress = 0
        fastaDict = read_fasta(fp)
        fastaSize = len(fastaDict)
        print(">Extracting cDNA")
        for name in fastaDict:
            id = name.split(" ")[0].split(".")[0][1:]
            dbPTM.loc[dbPTM['ensId'] == id, 'cdna'] = fastaDict[name]
            progress += 1
            sys.stdout.write('\r'+str("progress: ")+str(progress/fastaSize*100)+"%")
            sys.stdout.flush()
    #print(dbPTM["cdna"])
    dbPTM.reset_index(drop=True).to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)
    print(">DONE")

