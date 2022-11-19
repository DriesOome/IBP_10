import pandas as pd
import numpy as np
from math import floor, ceil
import statistics

def calculateSpeedVectors(dpPTM, windowSize, offsetStart):
    speedVectors = {}
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
            speedVectors[row["ensId"]+"/"+str(index)] = speedVector
        except Exception as e:
            failedCount += 1
            # print(e)
            # print("failed for protein")
    print(">Failed for "+str(failedCount)+"/"+str(len(dpPTM["ensId"]))+" entries")
    return speedVectors

def calculateAverageSpeed(speedVectors, offsetStart):
    averageSpeed = [0] * (offsetStart) * 2
    for pos in range(0, len(averageSpeed)):
        totalSpeedForPos = 0
        count = 0
        for speedVector in speedVectors:
            if speedVector[pos] != 0:
                count += 1
                totalSpeedForPos += speedVector[pos]
        averageSpeed[pos] = totalSpeedForPos/count
    return averageSpeed

def calculateMedianSpeed(speedVectors, offsetStart):
    medianSpeed = [0] * (offsetStart) * 2
    for pos in range(0, len(medianSpeed)):
        positionVector = []
        for speedVector in speedVectors:
            if speedVector[pos] != 0:
                positionVector.append(speedVector[pos])
        medianSpeed[pos] = statistics.median(positionVector)
    return medianSpeed