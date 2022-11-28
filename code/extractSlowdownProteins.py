from calculateSpeedVectors import calculateSpeedVectors, calculateAverageSpeed, calculateMedianSpeed
from extractCodons import splitCodons, translateCDNA
import pandas as pd
import seqlogo
import numpy as np
from PIL import EpsImagePlugin
EpsImagePlugin.gs_windows_binary = "C:/Program Files (x86)/gs32bit/gs10.00.0/bin/gswin32c.exe"
dpPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
windowSize = 5
offstetStart = 200
speedVectors = calculateSpeedVectors(dpPTM, windowSize, offstetStart)
averageSpeed = calculateMedianSpeed(speedVectors.values(), offstetStart)
targetPos = offstetStart+70
averageAtTargetPos = averageSpeed[targetPos]
print(">Threshold value: "+str(averageAtTargetPos))
slowdownVectors = [k.split("/")[0] for (k, v) in speedVectors.items() if v[targetPos]>averageAtTargetPos*1.15]
try:
    with open("../data/slowdownProtein.txt", 'w') as file:
        for protId in slowdownVectors:
            file.write(protId+"\n")
except IOError:
    print("I/O error")
def cDNAPM(dpPTM, slowdownVectors):
    PM = np.zeros(shape=(9, 4))
    totalCount = 0
    for protId in slowdownVectors:
        try:
            group = dpPTM.loc[dpPTM["ensId"] == protId]
            for rowIndex, row in group.iterrows():
                ptmIndex = row["index"] - 1
                codons = splitCodons(row["cdna"])
                if translateCDNA([codons[ptmIndex]]) == "N":
                    totalCount += 1
                    ptmSurrounding = ''.join(codons[ptmIndex - 1:ptmIndex + 2])
                    for m in range(0, 9):
                        if "A" == ptmSurrounding[m]:
                            PM[m][0] += 1
                        elif "C" == ptmSurrounding[m]:
                            PM[m][1] += 1
                        elif "G" == ptmSurrounding[m]:
                            PM[m][2] += 1
                        elif "T" == ptmSurrounding[m]:
                            PM[m][3] += 1
        except Exception as e:
            print(e)
    PM = PM / totalCount
    return PM

def proteinPM(dpPTM, slowdownVectors):
    PM = np.zeros(shape=(9, 20))
    totalCount = 0
    for protId in slowdownVectors:
        try:
            group = dpPTM.loc[dpPTM["ensId"] == protId]
            for rowIndex, row in group.iterrows():
                ptmIndex = row["index"] - 1
                seq = splitCodons(row["proteinSequence"])
                if seq[ptmIndex] == "N":
                    totalCount += 1
                    ptmSurrounding = seq[ptmIndex - 1:ptmIndex + 2]
                    # yes yes very poorly made excuse me
                    for m in range(0, 9):
                        if "A" == ptmSurrounding[m]:
                            PM[m][0] += 1
                        elif "C" == ptmSurrounding[m]:
                            PM[m][1] += 1
                        elif "D" == ptmSurrounding[m]:
                            PM[m][2] += 1
                        elif "E" == ptmSurrounding[m]:
                            PM[m][3] += 1
                        elif "F" == ptmSurrounding[m]:
                            PM[m][4] += 1
                        elif "G" == ptmSurrounding[m]:
                            PM[m][5] += 1
                        elif "H" == ptmSurrounding[m]:
                            PM[m][6] += 1
                        elif "I" == ptmSurrounding[m]:
                            PM[m][7] += 1
                        elif "K" == ptmSurrounding[m]:
                            PM[m][8] += 1
                        elif "L" == ptmSurrounding[m]:
                            PM[m][9] += 1
                        elif "M" == ptmSurrounding[m]:
                            PM[m][10] += 1
                        elif "N" == ptmSurrounding[m]:
                            PM[m][11] += 1
                        elif "P" == ptmSurrounding[m]:
                            PM[m][12] += 1
                        elif "Q" == ptmSurrounding[m]:
                            PM[m][13] += 1
                        elif "R" == ptmSurrounding[m]:
                            PM[m][14] += 1
                        elif "S" == ptmSurrounding[m]:
                            PM[m][15] += 1
                        elif "T" == ptmSurrounding[m]:
                            PM[m][16] += 1
                        elif "V" == ptmSurrounding[m]:
                            PM[m][17] += 1
                        elif "W" == ptmSurrounding[m]:
                            PM[m][18] += 1
                        elif "Y" == ptmSurrounding[m]:
                            PM[m][19] += 1
        except Exception as e:
            print(e)
    PM = PM / totalCount
    return PM
PM = proteinPM(dpPTM, slowdownVectors)
print(PM)
PM = seqlogo.Ppm(PM, alphabet_type="AA")
seqlogo.seqlogo(PM, filename="../figures/seqlogo.pdf",format="pdf")



