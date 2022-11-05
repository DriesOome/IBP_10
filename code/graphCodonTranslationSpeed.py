import math
from math import floor, ceil

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

# parameters
windowSize = 11
offsetStart = 200
sequonEfficiencyLower = -2
sequonEfficiencyUpper = 14

def calculateSpeedVectors():
    averagedSpeed = [0] * offsetStart * 2
    speedVectors = {}
    countPerPosition = [0] * offsetStart * 2
    for index, row in dpPTM.iterrows():
        speedVector = [0] * offsetStart * 2
        try:
            index = row["index"]
            frequencies = [float(f) for f in row["frequencies"][1:-1].split(", ")]
            for offset in range(-offsetStart, offsetStart):
                posInGene = index+offset
                if floor(windowSize/2) < posInGene < len(frequencies)-ceil(windowSize/2):
                    speedVector[offset+offsetStart] = sum(frequencies[posInGene-floor(windowSize/2):posInGene+ceil(windowSize/2)])/windowSize
                    countPerPosition[offset+offsetStart] += 1
                    #averagedSpeed[offset+offsetStart] += sum(frequencies[posInGene-floor(windowSize/2):posInGene+ceil(windowSize/2)])/windowSize
            averagedSpeed = np.add(averagedSpeed, speedVector)
            speedVectors[row["ensId"] + "/" + str(index)] = speedVector[windowSize:len(averagedSpeed)-windowSize]
        except Exception as e:
            print(e)
            print("failed for protein")
    averagedSpeed = np.divide(averagedSpeed, countPerPosition)
    averagedSpeed = averagedSpeed[windowSize:len(averagedSpeed)-windowSize]
    # for i in range(0, len(averagedSpeed)):
    #     averagedSpeed[i] = averagedSpeed[i]/countPerPosition[i]
    return averagedSpeed, speedVectors
# read dataframe
dpPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
dpPTM = dpPTM[dpPTM["sequonEfficiency"]<sequonEfficiencyUpper]
dpPTM = dpPTM[dpPTM["sequonEfficiency"]>sequonEfficiencyLower]
# building speed vector
averagedSpeed, speedVectors = calculateSpeedVectors()

# saving output
data = pd.DataFrame(columns=["x", "speed"])
data["x"] = list(range(-offsetStart+windowSize, offsetStart-windowSize))
data["averageSpeed"] = averagedSpeed
data.to_csv("../data/outputSpeed.txt", sep=",")

# plot
fig, ax = plt.subplots()
for speedVector in speedVectors.values():
    speedVector = [speed for speed in speedVector if speed != 0]
    ax.plot(list(range(-floor(len(speedVector)/2), ceil(len(speedVector)/2))), speedVector, linewidth=0.001)
#plt.axhline(y = np.mean([a for a in averagedSpeed[0:floor(len(averagedSpeed)/4)]]), color = 'r', linestyle = '-')
ax.plot(list(range(-floor(len(averagedSpeed)/2), ceil(len(averagedSpeed)/2))), averagedSpeed, linewidth=1)
plt.show()