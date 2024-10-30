#!/usr/bin/env python3

# Establish Environment
import sys
import argparse as ap
from datetime import datetime

import pandas as pd

# Log Alignments to both Console and FIle
class Logging:
    def __init__(self, filename):
        self.out_file = open(filename, "a")
        self.old_stdout = sys.stdout
        sys.stdout = self
    def write(self, text): 
        self.old_stdout.write(text)
        self.out_file.write(text)
    def __enter__(self): 
        return self
    def __exit__(self, type, value, traceback):
        sys.stdout = self.old_stdout

def float_range(min_val, max_val):
    def _float_range(value):
        try:
            value = float(value)
        except ValueError:
            raise ap.ArgumentTypeError(f"Invalid float value: {value}")

        if value < min_val or value > max_val:
            raise ap.ArgumentTypeError(f"Value must be in range [{min_val}, {max_val}]")
        return value
    return _float_range

# Process Arguments from Command Line to set desired PIP threshold, match type, and association file to search.
def processArgs():
    parser = ap.ArgumentParser(description='Finds matching variant calls between SPrime Introgression data and eQTL finemapping.')
    parser.add_argument('pip',default=0.5,help='Sets the desired PIP threshold.',type=float_range(0.0, 1.0))
    parser.add_argument('matchtype',default="either",help='Sets the desired match type.',choices=['d', 'dmatch', 'denisovan','n', 'nmatch', 'neanderthal','b', 'bmatch', 'both','e', 'ematch', 'either'])
    parser.add_argument('matchin',default="sig",help='Sets the desired eQTL associations file to search.',choices=['sig','all','any'])
    parser.add_argument('-ni','--nonIntrogressed',action='store_true',help='Returns non-introgressed eQTL finemaps.')
    args = parser.parse_args()

    pip = args.pip
    matchtype = args.matchtype.lower()
    matchin = args.matchin.lower()
    global nonIntrogressed
    nonIntrogressed = args.nonIntrogressed

    if matchtype in ["d", "dmatch", "denisovan"]:
        matchtype = "dmatch"
    elif matchtype in ["n", "nmatch", "neanderthal"]:
        matchtype = "nmatch"
    elif matchtype in ["b", "bmatch", "both"]:
        matchtype = "both"
    elif matchtype in ["e", "ematch", "either"]:
        matchtype = "either"

    if matchin in ["all", "any"]:
        matchin = "all"

    with Logging("log.txt"):
        if matchin == "sig":
            if matchtype == "dmatch":
                print(f"## Performing Alignment with PIP = {pip}, looking for Denisovan matches in SNPs with significant associations.")
            elif matchtype == "nmatch":
                print(f"## Performing Alignment with PIP = {pip}, looking for Neaderthal matches in SNPs with significant associations.")
            elif matchtype == "both":
                print(f"## Performing Alignment with PIP = {pip}, looking for both Denisovan and Neaderthal matches in SNPs with significant associations.")
            elif matchtype == "either":
                print(f"## Performing Alignment with PIP = {pip}, looking for either Denisovan or Neaderthal matches in SNPs with significant associations.")
        else:
            if matchtype == "dmatch":
                print(f"## Performing Alignment with PIP = {pip}, looking for Denisovan matches in SNPs with any association.")
            elif matchtype == "nmatch":
                print(f"## Performing Alignment with PIP = {pip}, looking for Neaderthal matches in SNPs with any association.")
            elif matchtype == "both":
                print(f"## Performing Alignment with PIP = {pip}, looking for both Denisovan and Neaderthal matches in SNPs with any association.")
            elif matchtype == "either":
                print(f"## Performing Alignment with PIP = {pip}, looking for either Denisovan or Neaderthal matches in SNPs with any association.")
    
    return pip, matchtype, matchin

def importeQTLs(pip, matchin):
    if matchin == "all":
        with Logging("log.txt"):
            print(f"## Importing eQTL for SNPs with any association.")
        eqtls = pd.read_csv("../data/data_clean/eQTL_finemapping/eQTL_finemapping.allAssociations.MAGE.v1.0.txt",sep="\t")

    if matchin == "sig":
        with Logging("log.txt"):
            print(f"## Importing eQTL for SNPs with significant associations.")
        eqtls = pd.read_csv("../data/data_clean/eQTL_finemapping/eQTL_finemapping.significantAssociations.MAGE.v1.0.txt",sep="\t")

    eqtls = eqtls.drop(columns=["variantRef","variant_kgpID","geneSymbol","variantCredibleSet"])
    eqtls = eqtls[eqtls.variantChrom != "chrX"]
    eqtls["variantChrom"] = eqtls["variantChrom"].str[3:]
    eqtls["variantChrom"] = eqtls["variantChrom"].astype("int64")
    eqtls = eqtls[eqtls.variantPIP>pip]
    baseschecked = eqtls.shape[0]
    return eqtls, baseschecked

def importSprimes(matchtype):
    with Logging("log.txt"):
            print(f"## Importing SPrime data for all Population Groups.")
    
    introgressions = pd.read_csv("../data/data_clean/Sprime_results/combinedSprimes.tsv",sep="\t")
    introgressions = introgressions.drop(columns=["REF","SEGMENT","ALLELE","SCORE"])
    if matchtype == "nmatch":
        introgressions = introgressions[introgressions.NMATCH=="match"]
    elif matchtype == "dmatch":
        introgressions = introgressions[introgressions.DMATCH=="match"]
    elif matchtype == "both":
        introgressions = introgressions[((introgressions["NMATCH"]=="match") & (introgressions["DMATCH"]=="match"))]
    elif matchtype == "either":
        introgressions = introgressions[((introgressions["NMATCH"]=="match") | (introgressions["DMATCH"]=="match"))]

    popgroups = [y for x, y in introgressions.groupby("POPGROUP")]
    
    totalintrogressions = 0
    for popgroup in range(len(popgroups)):
        totalintrogressions += popgroups[popgroup].shape[0]
    averageintrogressioncount = totalintrogressions / len(popgroups)
    return introgressions, averageintrogressioncount

# This actually performs the alignment, prints the alignment to a file keyed by start time and search parametes, then prints the Alignment Density and Introgression Enrichment in Selected Genes to both terminal and the log file for the run.
def alignSNPs(start_time, pip, matchin, matchtype, baseschecked, averageintrogressioncount, eqtls, introgressions):
   
    if nonIntrogressed == True:
        with Logging("log.txt"):
            print(f"## Looking for mismatches between eQTLs and Introgressions.")

        alignments = pd.merge(eqtls, introgressions,how="outer",on=['variantChrom','variantPosition','variantAlt','variant_rsID'],indicator=True).query('_merge == "left_only"').drop(columns=["NMATCH","DMATCH","POPGROUP","_merge"])
        alignments.to_csv(f"../data/data_clean/Alignments/{start_time}_alignment_{pip}_{matchin}_{matchtype}_mismatch.tsv",sep="\t")
        alignmentcount = alignments.shape[0]
        with Logging("log.txt"):
            print(f"## {alignmentcount} Non-introgressed SNPs were found out of {baseschecked} SNPs.")

    else:
        with Logging("log.txt"):
            print(f"## Looking for matches between eQTLs and Introgressions.")

        alignments = pd.merge(eqtls, introgressions,how='inner',on=['variantChrom','variantPosition','variantAlt','variant_rsID'])
        alignments.to_csv(f"../data/data_clean/Alignments/{start_time}_alignment_{pip}_{matchin}_{matchtype}.tsv",sep="\t")
        alignmentcount = alignments.shape[0]
    
        with Logging("log.txt"):
            print(f"## {alignmentcount} Alignments were found out of {averageintrogressioncount} introgressions across {baseschecked} SNPs.")
            print(f"## Alignment Density: {alignmentcount/baseschecked}")
            print(f"## Introgression Enrichment in Searched Genes: {alignmentcount/baseschecked/(averageintrogressioncount/3000000000)}")

def main():
    # Checks start time for keying of output file and runtime calculation.
    start_time = datetime.now()
    with Logging("log.txt"):
        print(f"New Alignment Session")
        print(f"Start Time: {start_time}")

    pip, matchtype, matchin = processArgs()

    eqtls, baseschecked = importeQTLs(pip, matchin)

    introgressions, averageintrogressioncount = importSprimes(matchtype)

    alignSNPs(start_time, pip, matchin, matchtype, baseschecked, averageintrogressioncount, eqtls, introgressions)

    # Logs stop time and total runtime to both console and the log file, then appends a trailing newline to visually separate runtime logs in the log file.
    with Logging("log.txt"):
        print(f"Stop Time: {datetime.now()}")
        print(f"Total Runtime: {datetime.now() - start_time}")
        print(f"")

# Only calls the main function if this script is being called as expected, otherwise the script just does nothing.
if __name__ == "__main__":
    main()