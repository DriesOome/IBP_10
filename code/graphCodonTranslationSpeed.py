import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

glycoSitesData = pd.read_csv("../data/codonsExtract.tsv", sep="\t")
offsetStart = 200
averagedSpeed = [0]*offsetStart*2
countPerPosition = [0]*offsetStart*2
sequonEfficiencyLower = 5
sequonEfficiencyUpper = 8
for index, row in glycoSitesData.iterrows():
    print(row)
    index = row["index"]
    frequencies = row["frequencies"][1:-2].split(", ")
    if sequonEfficiencyLower < row["sequonEfficiency"] < sequonEfficiencyUpper:
        try:
            for offset in range(-offsetStart, offsetStart):
                posInGene = index+offset
                if 0 < posInGene < len(frequencies):
                    averagedSpeed[offset+offsetStart] += float(frequencies[posInGene])
                    countPerPosition[offset+offsetStart] += 1
        except:
            print("failed for row: " + str(row["ensId"]))
for i in range(0, len(averagedSpeed)):
    averagedSpeed[i] = averagedSpeed[i]/countPerPosition[i]
print(averagedSpeed)


# plot
fig, ax = plt.subplots()

ax.plot(np.linspace(-offsetStart, offsetStart, len(averagedSpeed)), averagedSpeed, linewidth=2.0)

plt.show()