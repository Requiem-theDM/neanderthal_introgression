#!/bin/bash
echo -e "variantChrom\tvariantPosition\tID\tREF\tvariantAlt\tSEGMENT\tALLELE\tSCORE\tNMATCH\tDMATCH\tPOPGROUP"> ../data/data_clean/Sprime_results/combinedSprimes.tsv
for file in ../data/data_clean/Sprime_results/*.ND_match
do
    popgroup="${file#*/}"
    popgroup="${popgroup#*/}"
    popgroup="${popgroup#*/}"
    popgroup="${popgroup#*/}"
    popgroup="${popgroup%%.*}"
    python3 appendPopGroup.py ${file} ${popgroup} >> ../data/data_clean/Sprime_results/combinedSprimes.tsv
done