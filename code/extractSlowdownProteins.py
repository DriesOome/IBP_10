from calculateSpeedVectors import calculateSpeedVectors, calculateAverageSpeed, calculateMedianSpeed
import pandas as pd
dpPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
windowSize = 5
offstetStart = 200
speedVectors = calculateSpeedVectors(dpPTM, windowSize, offstetStart)
averageSpeed = calculateMedianSpeed(speedVectors.values(), offstetStart)
targetPos = offstetStart+70
averageAtTargetPos = averageSpeed[targetPos]
print(">Threshold value: "+str(averageAtTargetPos))
slowdownVectors = [k.split("/")[0] for (k, v) in speedVectors.items() if v[targetPos]>averageAtTargetPos*1.1]
print(slowdownVectors)
try:
    with open("../data/slowdownProtein.txt", 'w') as file:
        for protId in slowdownVectors:
            file.write(protId+"\n")
except IOError:
    print("I/O error")