import pandas as pd
import numpy as np
from joblib import Parallel, delayed

def calculateConservationForProtein(ensId: str, conservationData):
    try:
        #print("calculating conservation for: "+ ensId)
        conservationData["seq"] = conservationData["seq"].astype("string")
        countVector = [0]*len(str(conservationData.loc[0]["seq"]))
        for rowIndex in range(1, len(conservationData["seq"])):
            #print(str(rowIndex) + " - comparing to: " + conservationData.loc[rowIndex]["seq"])
            for i in range(1, len(str(conservationData.loc[0]["seq"]))):
                if str(conservationData.loc[0]["seq"])[i] == conservationData.loc[rowIndex]["seq"][i]:
                    countVector[i] += 1
        countVector = [i/100 for i in countVector]
        #print("countvector: " + str(countVector))
        return ensId, countVector
    except Exception as e:
        print("Failed for protein: "+ensId)
        #print(e)

def calculateConservation():
    print(">Calculating conservation - this might take a while... (20-40 min est.)")
    data = pd.read_csv("../data/compiledMultizDatabase.csv", names=["id", "seq"])
    data["organism"] = [item.split(".")[1].split("_")[1] for item in data["id"]]
    data["id"] = [item.split(".")[0] for item in data["id"]]

    dbPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
    dbPTM["conservationArray"] = np.nan

    results = Parallel(n_jobs=-5, verbose=10)( # all cpu's except for 1 will be used
        delayed(calculateConservationForProtein)(ensId, data.loc[data["id"] == ensId].reset_index(drop=True))
        for ensId in dbPTM["ensId"].unique()
    )
    for result in results:
        try:
            dbPTM.loc[dbPTM["ensId"] == result[0], "conservationArray"] = str(result[1])
        except:
            print("")
    dbPTM.to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)
    print(">DONE")
