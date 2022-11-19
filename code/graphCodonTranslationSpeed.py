import math
from math import floor, ceil

import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt
import statistics
from calculateSpeedVectors import calculateAverageSpeed, calculateSpeedVectors, calculateMedianSpeed
# parameters
windowSize = 5
offsetStart = 110
sequonEfficiencyLower = -1
sequonEfficiencyUpper = 14
conservationUpper = 101
conservationLower = -1




# read dataframe
dpPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
dpPTM = dpPTM[dpPTM["sequonEfficiency"]<sequonEfficiencyUpper]
dpPTM = dpPTM[dpPTM["sequonEfficiency"]>sequonEfficiencyLower]
dpPTM = dpPTM[dpPTM["conservation"]<conservationUpper]
dpPTM = dpPTM[dpPTM["conservation"]>conservationLower]
# dpPTM = dpPTM.head(100)
# building speed vector
speedVectors = calculateSpeedVectors(dpPTM, windowSize, offsetStart).values()
averageSpeed = calculateAverageSpeed(speedVectors, offsetStart)
medianSpeed = calculateMedianSpeed(speedVectors, offsetStart)
# saving output
data = pd.DataFrame(columns=["x", "speed"])
data["x"] = list(range(-offsetStart+windowSize, offsetStart-windowSize))
data["averageSpeed"] = averageSpeed[windowSize:len(averageSpeed) - windowSize]
data.to_csv("../data/outputSpeed.txt", sep=",")

# plot
plt.figure(figsize=(10,5))
for speedVector in speedVectors:
    speedVector = [speed for speed in speedVector if speed !=0]
    plt.plot(list(range(-floor(len(speedVector)/2), ceil(len(speedVector)/2))),speedVector,linewidth=0.005)
plt.axhline(y = np.mean([a for a in averageSpeed[0:floor(len(averageSpeed)/3)]]), color = 'b', linestyle = '-')
plt.plot(list(range(-floor(len(averageSpeed)/2), ceil(len(averageSpeed)/2))), averageSpeed, color = 'r', linewidth=1)
plt.plot(list(range(-floor(len(medianSpeed)/2), ceil(len(medianSpeed)/2))), medianSpeed, color = 'g', linewidth=1)
plt.xlabel("Positions relative to sequon")
plt.ylabel("DecodingTime")
title = "Sequon position vs Decoding time\nConservation interval "+str(conservationLower)+"-"+str(conservationUpper)
plt.title(title)
plt.savefig("../figures/latest.png", format="png", dpi=150)
plt.show()