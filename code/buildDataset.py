from compileMultizDataset import compileMultizDataset
from extractEnsId import extractEnsId
from extractProteinSequence import extractProteinSequence
from sequonEfficiencyCalc import  extractSequonEfficiency
from extractcdna import extractCDNA
from extractCodons import extractCodons
from calculateConservation import calculateConservation
# CODE full analysis
extractEnsId()
compileMultizDataset()
extractProteinSequence()
extractSequonEfficiency()
extractCDNA()
extractCodons()
calculateConservation()
