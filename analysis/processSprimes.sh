#!/bin/bash
echo -e "CHROM\tPOS\tID\tREF\ALT\tSEGMENT\tALLELE\tSCORE\tNMATCH\tDMATCH\tPOSGROUP"> ../data/data_clean/Sprime_results/combinedSprimes.tsv
for file in ../data/data_clean/Sprime_results/*.ND_match
do
    popgroup="${file#*/}"
    popgroup="${popgroup#*/}"
    popgroup="${popgroup#*/}"
    popgroup="${popgroup#*/}"
    popgroup="${popgroup%%.*}"
    python3 appendPopGroup.py ${file} ${popgroup} >> ../data/data_clean/Sprime_results/combinedSprimes.tsv
done