
import pandas as pd
import re
from Sequon_efficiency import calcSequonEfficiency

data = open('knownCanonical.multiz100way.exonAA.fa')
lines = data.readlines()
current_protein = 'start'
output = pd.DataFrame(columns=['protein','glycosite','amount','efficiency'])
select = 'init'
glycosite = []
amount = []
index = 0
efficiency = []
endCheckpos1 = False
endCheckpos2 = False
amountpostpos1 = [0]
amountpostpos2 = [0]
amountprepos1 = [0]
amountprepos2 = [0]

for line in lines:
    splitline = re.split(' |_', line)
    
    if splitline[0].startswith('>'):
        if splitline[0] != current_protein:
            if select == 'check':
                for i in range(0,len(glycosite)):
                    if amount[i] == 1 and endCheckpos1 == True and i == 0:
                        for x in range(1,len(amountpostpos1)):
                            if (amountpostpos1[x] and amountprepos1[x]) == 1:
                                amount[i] += 1
                        print([current_protein,glycosite[i]+index,amount[i],efficiency[i]])
                    if amount[i] == 1 and endCheckpos2 == True and i == (0 or 1):
                        for x in range(1,len(amountpostpos1)):
                            if (amountpostpos2[x] and amountprepos2[x]) == 1:
                                amount[i] += 1
                        print([current_protein,glycosite[i]+index,amount[i],efficiency[i]])
                    #print(current_protein,glycosite[i],amount[i])
                    output.loc[len(output)]=[current_protein,glycosite[i]+index,amount[i],efficiency[i]]
                    #output.append([current_protein,glycosite[i],amount[i]])
                    #pd.concat(output,pd.Series([current_protein,glycosite[i],amount[i]]))
            current_protein = splitline[0]
            
            index_after = splitline[4]
            index = 0
            select = 'first'
        elif splitline[0] == current_protein and splitline[1] == 'hg38':
            if select == 'check':
                for i in range(0,len(glycosite)):
                    if amount[i] == 1 and endCheckpos1 == True and i == 0:
                        for x in range(1,len(amountpostpos1)):
                            if (amountpostpos1[x] and amountprepos1[x]) == 1:
                                amount[i] += 1
                        print([current_protein,glycosite[i]+index,amount[i],efficiency[i]])
                    if amount[i] == 1 and endCheckpos2 == True and i == (0 or 1):
                        for x in range(1,len(amountpostpos1)):
                            if (amountpostpos2[x] and amountprepos2[x]) == 1:
                                amount[i] += 1
                        print([current_protein,glycosite[i]+index,amount[i],efficiency[i]])
                    #print(current_protein,glycosite[i],amount[i])
                    #print(current_protein,glycosite[i],amount[i])
                    output.loc[len(output)]=[current_protein,glycosite[i]+index,amount[i],efficiency[i]]
            index = index + int(index_after)
            index_after = splitline[4]
            select = 'first'
    elif select == 'first':
        glycosite = []
        if endCheckpos1 == True:
            if list(line)[1] == 'S' or list(line)[1] == 'T':
                glycosite.append(-2)
        if endCheckpos2 == True:
            if list(line)[1] != 'P':
                if list(line)[2] == 'S' or list(line)[2] == 'T':
                    glycosite.append(-1)
        for i in range(0,len(list(line))-2):
            if list(line)[i] == 'N':
                if list(line)[i+1] != 'P':
                    if list(line)[i+2] == 'T' or list(line)[i+2] == 'S':
                        glycosite.append(i)
                        efficiency.append(calcSequonEfficiency(line, i))
                        select = 'check'
        if len(glycosite) == 0:
            select = 'skip'
        amount = [1]*len(glycosite)
        if list(line)[len(list(line))-2] == 'N' and list(line)[len(list(line))-1] != 'P':
            endCheckpos1 = True
            amountprepos1 = [1]
            amountpostpos1 = [1]
            select ='check'
        elif list(line)[len(list(line))-1] == 'N':
            endCheckpos2 = True
            amountpos2 = [1]
            amountpostpos2 = [1]
            select ='check'
        else:
            endCheckpos1 = False
            endCheckpos2 = False
    elif select == 'check':
        if endCheckpos1 == True :
            if list(line)[len(list(line))-2] == 'N' and list(line)[len(list(line))-1] != 'P':
                amountprepos1.append(1)
            else:
                amountprepos1.append(0)
        if endCheckpos2 == True : 
            if list(line)[len(list(line))-1] == 'N':
                amountprepos2.append(1)
            else:
                amountprepos2.append(0)
        else: 
            for i in range(0, len(glycosite)):
                try: 
                    if glycosite[i] > 0: 
            
                        if list(line)[glycosite[i]] == 'N':
                            if list(line)[glycosite[i]+1] != 'P':
                                if list(line)[glycosite[i]+2] == 'S' or list(line)[glycosite[i]+2] == 'T':
                                    amount[i] = amount[i] + 1
                    elif glycosite[i] == -1:
                        if list(line)[1] != 'P':
                            if list(line)[2] == 'S' or list(line)[2] == 'T':
                                amountpostpos1.append(1)
                        else: 
                            amountpostpos1.append(0)
                    elif glycosite[i] == -2:
                        if list(line)[2] == 'S' or list(line)[2] == 'T':
                            amountpostpos2.append(1)
                        else:
                            amountpostpos2.append(0)
                except:
                    pass
                    
    
    output.to_csv('efficiency&abundanceOfSequonsMultiz',sep="\t")
            
            