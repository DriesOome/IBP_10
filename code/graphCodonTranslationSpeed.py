import math
from math import floor, ceil

import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt

# parameters
windowSize = 5
offsetStart = 110
sequonEfficiencyLower = -1
sequonEfficiencyUpper = 5
conservationUpper = 20
conservationLower = -1

def calculateSpeedVectors():
    averagedSpeed = [0] * offsetStart * 2
    speedVectors = []
    countPerPosition = [0] * offsetStart * 2
    failedCount = 0
    for index, row in dpPTM.iterrows():
        speedVector = [0] * offsetStart * 2
        try:
            index = row["index"]
            if row["proteinSequence"][index-1] != "N":
                raise Exception("Index does not point to an 'N' residue")
            frequencies = [float(f) for f in row["frequencies"][1:-1].split(", ")]
            for offset in range(-offsetStart, offsetStart):
                posInGene = index+offset
                if floor(windowSize/2) < posInGene < len(frequencies)-ceil(windowSize/2):
                    speedVector[offset+offsetStart] = sum(frequencies[posInGene-floor(windowSize/2):posInGene+ceil(windowSize/2)])/windowSize

                    countPerPosition[offset+offsetStart] += 1
            averagedSpeed = np.add(averagedSpeed, speedVector)
            speedVectors.append(speedVector[windowSize:len(averagedSpeed)-windowSize])
        except Exception as e:
            failedCount += 1
            # print(e)
            # print("failed for protein")
    print(">Failed for "+str(failedCount)+"/"+str(len(dpPTM["ensId"]))+" entries")
    averagedSpeed = np.divide(averagedSpeed, countPerPosition)
    averagedSpeed = averagedSpeed[windowSize:len(averagedSpeed)-windowSize]
    # for i in range(0, len(averagedSpeed)):
    #     averagedSpeed[i] = averagedSpeed[i]/countPerPosition[i]
    return averagedSpeed, speedVectors
# read dataframe
dpPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
dpPTM = dpPTM[dpPTM["sequonEfficiency"]<sequonEfficiencyUpper]
dpPTM = dpPTM[dpPTM["sequonEfficiency"]>sequonEfficiencyLower]
dpPTM = dpPTM[dpPTM["conservation"]<conservationUpper]
dpPTM = dpPTM[dpPTM["conservation"]>conservationLower]
# dpPTM = dpPTM.head(100)
# building speed vector
averagedSpeed, speedVectors = calculateSpeedVectors()
# saving output
data = pd.DataFrame(columns=["x", "speed"])
data["x"] = list(range(-offsetStart+windowSize, offsetStart-windowSize))
data["averageSpeed"] = averagedSpeed
data.to_csv("../data/outputSpeed.txt", sep=",")

# plot
plt.figure(figsize=(10,5))
for speedVector in speedVectors:
    plt.plot(list(range(-floor(len(speedVector)/2), ceil(len(speedVector)/2))),speedVector,linewidth=0.005)
plt.axhline(y = np.mean([a for a in averagedSpeed[0:floor(len(averagedSpeed)/3)]]), color = 'b', linestyle = '-')
plt.plot(list(range(-floor(len(averagedSpeed)/2), ceil(len(averagedSpeed)/2))), averagedSpeed, color = 'r', linewidth=1)
plt.xlabel("Positions relative to sequon")
plt.ylabel("DecodingTime")
title = "Sequon position vs Decoding time\nConservation interval "+str(conservationLower)+"-"+str(conservationUpper)
plt.title(title)
plt.savefig("../figures/latest.png", format="png", dpi=150)
plt.show()