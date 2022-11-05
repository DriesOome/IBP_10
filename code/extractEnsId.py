import pandas as pd

def extractEnsId():
    print(">Linking uniprotIds to ensIds")
    fromToData = pd.read_csv("../data/startingDatasets/ensembleToUniprot.tsv", sep="\t", names=["ensId", "uniProtId", "protName"])
    dpPTM = pd.read_csv("../data/startingDatasets/dbPTMdataset.tsv", sep="\t")
    # select all ids for which we have a glyco site
    fromToData = fromToData[fromToData["uniProtId"].isin(dpPTM["uniProtId"])].reset_index(drop=True).drop(labels=["protName"], axis=1)
    dpPTM = pd.merge(dpPTM, fromToData, on="uniProtId", how="left").dropna().reset_index(drop=True)
    dpPTM[["proteinName", "uniProtId", "ensId", "index", "sequence"]].to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)
    print(">DONE")