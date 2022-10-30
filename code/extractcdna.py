import pandas as pd
import numpy as np

from readFasta import read_fasta


fromToData = pd.read_csv("../data/ensToUniMultiz.tsv", sep="\t", names=["ensId", "uniProtId", "protName"])
glycoSitesData = pd.read_csv("../data/humanGlycoSitesScored.tsv", sep="\t")
# select all ids for which we have a glyco site
fromToData = fromToData[fromToData["uniProtId"].isin(glycoSitesData["uniProtId"])].reset_index(drop=True).drop(labels=["protName"], axis=1)
glycoSitesData = pd.merge(glycoSitesData, fromToData, on="uniProtId", how="left").dropna().reset_index(drop=True)

glycoSitesData["cdna"] = np.nan
with open('../data/Homo_sapiens.cdna.all.fa') as fp:
    progress = 0
    for name, seq in read_fasta(fp):
        id = name.split(" ")[0].split(".")[0][1:]
        glycoSitesData.loc[glycoSitesData['ensId'] == id, 'cdna'] = seq
        progress += 1
        print(progress)
print(glycoSitesData["cdna"])
glycoSitesData.reset_index(drop=True).to_csv("../data/glycoSitesCdna.tsv", sep="\t")


