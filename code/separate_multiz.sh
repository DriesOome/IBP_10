#!/bin/bash

#Separating entries (delimiter = 2 new lines):
awk -v RS='\n\n\n' '{print > ("enst-" NR ".txt")}' knownCanonical.multiz100way.exonAA.fa

#Renaming files:
for i in enst-*;
do
 ensemblid=$(grep -o ENST[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9] "$i" | sort -u);
 new_filename="${ensemblid}.txt";
 if [ ! -f "$new_filename" ]; then
    mv "$i" "$new_filename";
 fi;
done
