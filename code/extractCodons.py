import pandas as pd
import numpy as np
def clipCDNA(seq: str):
    maxPos = len(seq)-1
    pos = 0
    while maxPos-pos > 3:
        if seq[pos] == "A" and seq[pos+1] == "T" and seq[pos+2] == "G":
            return seq[pos:]
        else:
            pos += 1

    print("no atg found")
    return -1

def splitCodons(seq: str):
    maxPos = len(seq)-1
    pos = 0
    codons = []
    while maxPos-pos > 3:
        codons.append(seq[pos:pos+3])
        pos += 3
    return codons

codonFrequencies = {
    "TTT": 17.6, "TCT": 15.2, "TAT": 12.2, "TGT": 10.6,
    "TTC": 20.3, "TCC": 17.7,  "TAC": 15.3,  "TGC": 12.6,
    "TTA":  7.7,  "TCA": 12.2,  "TAA":  1.0,  "TGA":  1.6,
    "TTG": 12.9, "TCG":  4.4,  "TAG":  0.8, "TGG": 13.2,

    "CTT": 13.2, "CCT": 17.5,  "CAT": 10.9, "CGT":  4.5,
    "CTC": 19.6,  "CCC": 19.8, "CAC": 15.1,  "CGC": 10.4,
    "CTA":  7.2,  "CCA": 16.9,  "CAA": 12.3,  "CGA":  6.2,
    "CTG": 39.6,  "CCG":  6.9, "CAG": 34.2,  "CGG": 11.4,

    "ATT": 16.0,  "ACT": 13.1,  "AAT": 17.0,  "AGT": 12.1,
    "ATC": 20.8,  "ACC": 18.9,  "AAC": 19.1, "AGC": 19.5,
    "ATA":  7.5, "ACA": 15.1,  "AAA": 24.4,  "AGA": 12.2,
    "ATG": 22.0,  "ACG":  6.1,  "AAG": 31.9, "AGG": 12.0,

    "GTT": 11.0,  "GCT": 18.4,  "GAT": 21.8,  "GGT": 10.8,
    "GTC": 14.5,  "GCC": 27.7, "GAC": 25.1, "GGC": 22.2,
    "GTA":  7.1,  "GCA": 15.8,  "GAA": 29.0,  "GGA": 16.5,
    "GTG": 28.1,  "GCG":  7.4,  "GAG": 39.6,  "GGG": 16.5
}
fractionFrequencies = {
    "TTT": 0.45, "TCT": 0.18, "TAT": 0.43, "TGT": 0.45,
    "TTC": 0.55, "TCC": 0.22, "TAC": 0.57, "TGC": 0.55,
    "TTA": 0.07, "TCA": 0.15, "TAA": 0.28, "TGA": 0.52,
    "TTG": 0.13, "TCG": 0.06, "TAG": 0.20, "TGG": 1.00,

    "CTT": 0.13, "CCT": 0.28, "CAT": 0.41, "CGT": 0.08,
    "CTC": 0.20, "CCC": 0.33, "CAC": 0.59, "CGC": 0.19,
    "CTA": 0.07, "CCA": 0.27, "CAA": 0.25, "CGA": 0.11,
    "CTG": 0.41, "CCG": 0.11, "CAG": 0.75, "CGG": 0.21,

    "ATT": 0.36, "ACT": 0.24, "AAT": 0.46, "AGT": 0.15,
    "ATC": 0.48, "ACC": 0.36, "AAC": 0.54, "AGC": 0.24,
    "ATA": 0.16, "ACA": 0.28, "AAA": 0.42, "AGA": 0.20,
    "ATG": 1,    "ACG": 0.12, "AAG": 0.58, "AGG": 0.2,

    "GTT": 0.18, "GCT": 0.26, "GAT": 0.46, "GGT": 0.16,
    "GTC": 0.24, "GCC": 0.40, "GAC": 0.54, "GGC": 0.34,
    "GTA": 0.11, "GCA": 0.23, "GAA": 0.42, "GGA": 0.25,
    "GTG": 0.47, "GCG": 0.11, "GAG": 0.58, "GGG": 0.25
}

def assignTranslationSpeedScores(codons: [str]):
    frequencies = []
    for codon in codons:
        frequencies.append(fractionFrequencies.get(codon))
    return frequencies


glycoSitesData = pd.read_csv("../data/glycoSitesCdna.csv", sep=",").drop(labels=["Unnamed: 0.1", "Unnamed: 0"], axis=1)
glycoSitesData = glycoSitesData.dropna().reset_index(drop=True)

glycoSitesData["cdna"] = [clipCDNA(cdna) for cdna in glycoSitesData["cdna"]]
glycoSitesData["codons"] = np.nan
glycoSitesData["codons"] = [splitCodons(cdna) for cdna in glycoSitesData["cdna"]]
glycoSitesData["frequencies"] = np.nan
glycoSitesData["frequencies"] = [assignTranslationSpeedScores(codons) for codons in glycoSitesData["codons"]]
glycoSitesData.to_csv("../data/codonsExtract.tsv", sep="\t")