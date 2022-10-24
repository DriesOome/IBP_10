import pandas

# POU => G, Y, N, S, C, Q, T
POU = ["G","Y","N","S","C","Q","T"]
# non positive => NOT H, R, K
POS = ["H","R","K"]
# function for sequon efficiency
def calcSequonEfficiency(sequence: str, glycoIndex: int):
    print(sequence)
    score: float = -0.5070
    # check if index is correct
    if sequence[glycoIndex] != "N":
        print("given index does not point to N residue")
        return -1
    # calculate score
    # if nonpositive at i-2 add 1.3292
    if sequence[glycoIndex-2] not in POS:
        score += 1.3292
    # if F, Y, W, C, M at i-2 add 2.1615
    if sequence[glycoIndex-2] in ["F", "Y", "W", "C", "M"]:
        score += 2.1615
    # if POU at i-1 add 3.6061
    if sequence[glycoIndex - 1] in POU:
        score += 3.6061
    # if C, M, S, T, A, V, I, L, G at i-1 add +3.0451
    if sequence[glycoIndex - 1] in ["C", "M", "S", "T", "C", "Q", "T"]:
        score += 3.0451
    # if POU at i+1 add 3.9373
    if sequence[glycoIndex+1] in POU:
        score += 3.9373

    return score


# read in data
data = pandas.read_csv('../data/GlycosylationSites',sep='\t')
# fetch human only data
data = data.loc[data["proteinName"].str.contains("HUMAN")].dropna()
print(data.head())
# calculate sequon efficiency test
print("TEST 1")
print(calcSequonEfficiency("ARGLTNYSKIL", 5) == 11.4107)
print("TEST 2")
print(calcSequonEfficiency("ARGADNGTKIL", 5) == 4.7595)
# calculate for the whole dataset
data["sequonEfficiency"] = [calcSequonEfficiency(sequence, 10) for sequence in data["sequence"]]
# write to file
data = data.drop(data[data["sequonEfficiency"] == -1.0].index)
data = data.reset_index()
data.to_csv("../data/humanGlycoSitesScored.tsv", sep="\t")



