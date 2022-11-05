import pandas as pd
import numpy as np
# reading in multiz dataset
def extractProteinSequence():
    print(">Extracting protein sequence from multiz")
    multizData = pd.read_csv("../data/compiledMultizDatabase.csv", names=["id", "seq"])
    multizData["organism"] = [item.split(".")[1].split("_")[1] for item in multizData["id"]]
    multizData["id"] = [item.split(".")[0] for item in multizData["id"]]
    multizData = multizData[multizData["organism"]=="hg38"]

    dbPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
    dbPTM["proteinSequence"] = np.nan
    for ensId in dbPTM["ensId"]:
        try:
            multizEntry = multizData[multizData["id"] == ensId]
            sequence = str(multizEntry.iloc[0]["seq"])
            dbPTM.loc[dbPTM["ensId"]==ensId, "proteinSequence"] = sequence
        except Exception as e:
            print("* Failed for protein: "+ensId)
    dbPTM.drop(dbPTM[dbPTM["proteinSequence"] == "Nan"].index).to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)
    print(">DONE")
