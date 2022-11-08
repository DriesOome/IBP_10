import sys

import pandas as pd
import numpy as np

def checkIfCDNAisProtein(cdnaSequence: str, protSequence):
    translatedProt = translateCDNA(splitCodons(cdnaSequence[0:30]))
    return translatedProt == protSequence[0:9]

def clipCDNA(entry):
    cdnaSeq = entry["cdna"]
    protSeq = entry["proteinSequence"]
    maxPos = len(cdnaSeq) - 1
    pos = 0
    while maxPos-pos > 3:
        if cdnaSeq[pos] == "A" and cdnaSeq[pos + 1] == "T" and cdnaSeq[pos + 2] == "G":
            if checkIfCDNAisProtein(cdnaSeq[pos:], protSeq):
                return cdnaSeq[pos:(pos+len(protSeq)*3)]
            else:
                pos+=1
        else:
            pos += 1

    print("no atg found")
    return -1

geneticCode = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAA": "stop", "TAG": "stop", "TGA": "stop"
}

def translateCDNA(codons):
    seq = []
    for codon in codons:
        seq.append(geneticCode.get(codon))
    return ''.join(seq)

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
# https://www.cs.tau.ac.il/~tamirtul/MTDR/mu_vals.html
decodingTime = {
"AAA":0.14788,
"AAC":0.14034,
"AAG":0.12458,
"AAT":0.15498,
"ACA":0.15344,
"ACC":0.12808,
"ACG":0.13043,
"ACT":0.15349,
"AGA":0.14815,
"AGC":0.13401,
"AGG":0.14382,
"AGT":0.15579,
"ATA":0.16686,
"ATC":0.12393,
"ATG":0.14733,
"ATT":0.13501,
"CAA":0.17789,
"CAC":0.14571,
"CAG":0.13826,
"CAT":0.16813,
"CCA":0.15097,
"CCC":0.11139,
"CCG":0.10966,
"CCT":0.14516,
"CGA":0.15941,
"CGC":0.08957,
"CGG":0.10514,
"CGT":0.12995,
"CTA":0.15401,
"CTC":0.13435,
"CTG":0.13518,
"CTT":0.16164,
"GAA":0.16382,
"GAC":0.14829,
"GAG":0.13094,
"GAT":0.16581,
"GCA":0.15615,
"GCC":0.1045,
"GCG":0.08279,
"GCT":0.13372,
"GGA":0.15957,
"GGC":0.09991,
"GGG":0.11195,
"GGT":0.12311,
"GTA":0.14702,
"GTC":0.13617,
"GTG":0.12615,
"GTT":0.14843,
"TAC":0.15038,
"TAT":0.15317,
"TCA":0.16693,
"TCC":0.13484,
"TCG":0.1222,
"TCT":0.15756,
"TGC":0.13681,
"TGG":0.15278,
"TGT":0.16268,
"TTA":0.17189,
"TTC":0.14625,
"TTG":0.16597,
"TTT":0.15944,
"TAA": 0,
"TGA": 0,
"TAG": 0
}

#Function that calculates codon frequency within a protein
def calcFrequencyInProt(codons):
    #input = list of codons for a protein
    dict = {}
    #total number of codons
    total = len(codons)
    
    for codon in codons:
    # Add 1 to the value corresponding to the 'codon' key (and set to 1 if it doesn't exist yet)
        dict[codon] = dict.get(codon, 0) + 1
    #dict.items() returns a list of occurrences, so we can use that alongside 'total' to compute relative frequency
    relativeFreq = [(codon, count / total * 100) for codon, count in dict.items()]
    
    #printing result in %
    for codon, pct in sorted(relativeFreq, key=lambda x: x[0]):
        print(codon, f"{pct:.2f}%")
#    print ("Count of all codons in the current protein is : \n" + str(dict)) 
    
    return dict

def assignTranslationSpeedScores(codons: [str]):
    frequencies = []
    dict = calcFrequencyInProt(codons)
    for codon in codons:
        frequencies.append(dict[codon])
    return frequencies

def extractCodons():
    print(">Extracting codons")
    dpPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
    dpPTM = dpPTM.dropna().reset_index(drop=True)
    dpPTM["cdna"] = [clipCDNA(row) for index, row in dpPTM.iterrows()]
    dpPTM = dpPTM.drop(dpPTM[dpPTM["cdna"] == -1.0].index)
    dpPTM["codons"] = np.nan
    dpPTM["codons"] = [splitCodons(cdna) for cdna in dpPTM["cdna"]]

    print(">Assigning translation speed scores")
    dpPTM["frequencies"] = np.nan
    dpPTM["frequencies"] = [assignTranslationSpeedScores(codons) for codons in dpPTM["codons"]]
    dpPTM.to_csv("../data/dpPTMextended.tsv", sep="\t", index=False)
    print(">DONE")

extractCodons()