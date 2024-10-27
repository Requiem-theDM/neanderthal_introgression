#!/usr/bin/env python3

# Establish Environment
import sys
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

# Process Arguments from Command Line to set desired PIP threshold, match type, and association file to search.
def processArgs():
    try: sys.argv[1].isdigit()
    except:
        pip = 0.5
        with Logging("log.txt"):
            print(f"## Using Default 0.5 PIP.")
    else: pip = float(sys.argv[1])

    try: sys.argv[2]
    except:
        matchtype = "nmatch"
        with Logging("log.txt"):
            print(f"## Using Default Neanderthal Match.")
    else:
        if sys.argv[2].lower() in ["d", "dmatch", "denisovan"]:
            matchtype = "dmatch"
        elif sys.argv[2].lower() in ["n", "nmatch", "neanderthal"]:
            matchtype = "nmatch"
        elif sys.argv[2].lower() in ["b", "bmatch", "both"]:
            matchtype = "both"
        elif sys.argv[2].lower() in ["e", "ematch", "either"]:
            matchtype = "either"
        else:
            matchtype = "nmatch"

    try: sys.argv[3]
    except:
        matchin = "sig"
        with Logging("log.txt"):
            print(f"## Using Default Search in Significant Associations.")
    else:
        if sys.argv[3].lower() in ["all", "any"]:
            matchin = "all"
        else:
            matchin = "sig"

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
    with Logging("log.txt"):
        if matchin == "sig":
            print(f"## Importing eQTL for SNPs with significant associations.")
        else:
            print(f"## Importing eQTL for SNPs with any association.")

    if matchin == "all":
        eqtls = pd.read_csv("../data/data_clean/eQTL_finemapping/eQTL_finemapping.allAssociations.MAGE.v1.0.txt",sep="\t")
        eqtls = eqtls.drop(columns=["variantRef","variant_kgpID","variant_rsID","geneSymbol","variantCredibleSet"])
        eqtls = eqtls[eqtls.variantChrom != "chrX"]
        eqtls["variantChrom"] = eqtls["variantChrom"].str[3:]
        eqtls["variantChrom"] = eqtls["variantChrom"].astype("int64")
        eqtls = eqtls[eqtls.variantPIP>pip]
        baseschecked = eqtls.shape[0]

    if matchin == "sig":
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
    with Logging("log.txt"):
            print(f"## Looking for matches between eQTLs and Introgressions.")

    alignments = pd.merge(eqtls, introgressions)
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