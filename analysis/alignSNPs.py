#!/usr/bin/env python3

import sys
from datetime import datetime

class Sprime:
    def __init__(self,chrom,pos,snp,nmatch,dmatch,popgroup):
        self.chrom = chrom
        self.pos = pos
        self.snp = snp
        self.nmatch = nmatch
        self.dmatch = dmatch
        self.popgroup = popgroup

class eQTL:
    def __init__(self,chrom,pos,snp,geneID,pip):
        self.chrom = chrom
        self.pos = pos
        self.snp = snp
        self.geneID = geneID
        self.pip = pip

class SNP:
    def __init__(self,chrom,pos,snp,geneID,pip,nmatch,dmatch,popgroup):
        self.chrom = chrom
        self.pos = pos
        self.snp = snp
        self.geneID = geneID
        self.pip = pip
        self.nmatch = nmatch
        self.dmatch = dmatch
        self.popgroup = popgroup

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

def processArgs():
    try: sys.argv[1].isdigit()
    except: pip = 0.5
    else: pip = float(sys.argv[1])

    try: sys.argv[2]
    except: matchtype = "nmatch"
    else:
        if sys.argv[2].lower() == "d" or sys.argv[2].lower() == "dmatch" or sys.argv[2].lower() == "denisovan":
            matchtype = "dmatch"
        elif sys.argv[2].lower() == "n" or sys.argv[2].lower() == "nmatch" or sys.argv[2].lower() == "neanderthal":
            matchtype = "nmatch"
        elif sys.argv[2].lower() == "b" or sys.argv[2].lower() == "bmatch" or sys.argv[2].lower() == "both":
            matchtype = "both"
        elif sys.argv[2].lower() == "e" or sys.argv[2].lower() == "ematch" or sys.argv[2].lower() == "either":
            matchtype = "either"
        else:
            matchtype = "nmatch"

    try: sys.argv[3]
    except: matchin = "sig"
    else:
        if sys.argv[3].lower() == "all" or sys.argv[3].lower() == "any":
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

def importeQTL(pip, matchin, eqtls):
    with Logging("log.txt"):
        if matchin == "sig":
            print(f"## Importing eQTL for SNPs with significant associations.")
        else:
            print(f"## Importing eQTL for SNPs with any association.")

    if matchin == "sig": 
        with open(f"../data/data_clean/eQTL_finemapping/eQTL_finemapping.significantAssociations.MAGE.v1.0.txt") as file:
            file.readline()
            for line in file:
                fields = line.split()
                if fields[0].lstrip("chr") != "X" and (fields[3] == "A" or fields[3] == "T" or fields[3] == "C" or fields[3] == "G"):
                    if float(fields[8]) >= pip:
                        eqtls.append(eQTL(fields[0].lstrip("chr"), fields[1], fields[3], fields[6], float(fields[8])))

    if matchin == "all":
        with open(f"../data/data_clean/eQTL_finemapping/eQTL_finemapping.allAssociations.MAGE.v1.0.txt") as file:
            file.readline()
            for line in file:
                fields = line.split()
                if float(fields[8]) >= pip:
                    eqtls.append(eQTL(fields[0].lstrip("chr"), fields[1], fields[3], fields[6], float(fields[8])))

def importSprime(matchtype, introgressions):
    with Logging("log.txt"):
            print(f"## Importing SPrime data for all Population Groups.")

    with open(f"../data/data_clean/Sprime_results/combinedSprimes.tsv") as file:
        file.readline()
        for line in file:
            fields = line.split()
            if matchtype == "both" and fields[8] == "match" and fields[9] == "match":
                introgressions.append(Sprime(fields[0], fields[1], fields[4], fields[8], fields[9], fields[10]))
            
            elif matchtype == "either" and (fields[8] == "match" or fields[9] == "match"):
                introgressions.append(Sprime(fields[0], fields[1], fields[4], fields[8], fields[9], fields[10]))

            elif matchtype == "dmatch" and fields[9] == "match":
                introgressions.append(Sprime(fields[0], fields[1], fields[4], fields[8], fields[9], fields[10]))

            elif matchtype == "nmatch" and fields[8] == "match":
                introgressions.append(Sprime(fields[0], fields[1], fields[4], fields[8], fields[9], fields[10]))

def importDebugger(eqtls, introgressions):
    eqchromosomes = []
    eqsnp = []
    for eqtl in eqtls:
        eqchromosomes.append(f"{eqtl.chrom}")
        eqsnp.append(f"{eqtl.snp}")
    intrchromosomes = []
    intrsnp = []
    for introgression in introgressions:
        intrchromosomes.append(f"{introgression.chrom}")
        intrsnp.append(f"{introgression.snp}")

    with Logging("log.txt"):
        print(f"## The set of chromosomes in the eQTLs is: {set(eqchromosomes)}")
        print(f"## The set of chromosomes in the introgressions is: {set(intrchromosomes)}")
        print(f"## The set of SNPs in the eQTLs is: {set(eqsnp)}")
        print(f"## The set of SNPs in the introgressions is: {set(intrsnp)}")

def splitbyChromosome(eqtls, introgressions):
    with Logging("log.txt"):
        print(f"## Splitting eQTLs by chromosome.")

    eqtls_chrom1 = []
    eqtls_chrom2 = []
    eqtls_chrom3 = []
    eqtls_chrom4 = []
    eqtls_chrom5 = []
    eqtls_chrom6 = []
    eqtls_chrom7 = []
    eqtls_chrom8 = []
    eqtls_chrom9 = []
    eqtls_chrom10 = []
    eqtls_chrom11 = []
    eqtls_chrom12 = []
    eqtls_chrom13 = []
    eqtls_chrom14 = []
    eqtls_chrom15 = []
    eqtls_chrom16 = []
    eqtls_chrom17 = []
    eqtls_chrom18 = []
    eqtls_chrom19 = []
    eqtls_chrom20 = []
    eqtls_chrom21 = []
    eqtls_chrom22 = []

    for eqtl in eqtls:
        if eqtl.chrom == "1": eqtls_chrom1.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "2": eqtls_chrom2.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "3": eqtls_chrom3.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "4": eqtls_chrom4.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "5": eqtls_chrom5.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "6": eqtls_chrom6.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "7": eqtls_chrom7.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "8": eqtls_chrom8.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "9": eqtls_chrom9.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "10": eqtls_chrom10.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "11": eqtls_chrom11.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "12": eqtls_chrom12.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "13": eqtls_chrom13.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "14": eqtls_chrom14.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "15": eqtls_chrom15.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "16": eqtls_chrom16.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "17": eqtls_chrom17.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "18": eqtls_chrom18.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "19": eqtls_chrom19.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "20": eqtls_chrom20.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "21": eqtls_chrom21.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))
        elif eqtl.chrom == "22": eqtls_chrom22.append(eQTL(eqtl.chrom,eqtl.pos,eqtl.snp,eqtl.geneID,eqtl.pip))

    eqtls = [eqtls_chrom1,eqtls_chrom2,eqtls_chrom3,eqtls_chrom4,eqtls_chrom5,eqtls_chrom6,eqtls_chrom7,eqtls_chrom8,eqtls_chrom9,eqtls_chrom10,eqtls_chrom11,eqtls_chrom12,eqtls_chrom13,eqtls_chrom14,eqtls_chrom15,eqtls_chrom16,eqtls_chrom17,eqtls_chrom18,eqtls_chrom19,eqtls_chrom20,eqtls_chrom21,eqtls_chrom22]
    with Logging("log.txt"):
        print(f"## eQTLs were found for {len(eqtls)} chromosomes.")

    with Logging("log.txt"):
        print(f"## Splitting introgressions by chromosome.")

    introgressions_chrom1 = []
    introgressions_chrom2 = []
    introgressions_chrom3 = []
    introgressions_chrom4 = []
    introgressions_chrom5 = []
    introgressions_chrom6 = []
    introgressions_chrom7 = []
    introgressions_chrom8 = []
    introgressions_chrom9 = []
    introgressions_chrom10 = []
    introgressions_chrom11 = []
    introgressions_chrom12 = []
    introgressions_chrom13 = []
    introgressions_chrom14 = []
    introgressions_chrom15 = []
    introgressions_chrom16 = []
    introgressions_chrom17 = []
    introgressions_chrom18 = []
    introgressions_chrom19 = []
    introgressions_chrom20 = []
    introgressions_chrom21 = []
    introgressions_chrom22 = []

    for introgression in introgressions:
        if introgression.chrom == "1": introgressions_chrom1.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "2": introgressions_chrom2.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "3": introgressions_chrom3.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "4": introgressions_chrom4.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "5": introgressions_chrom5.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "6": introgressions_chrom6.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "7": introgressions_chrom7.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "8": introgressions_chrom8.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "9": introgressions_chrom9.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "10": introgressions_chrom10.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "11": introgressions_chrom11.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "12": introgressions_chrom12.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "13": introgressions_chrom13.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "14": introgressions_chrom14.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "15": introgressions_chrom15.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "16": introgressions_chrom16.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "17": introgressions_chrom17.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "18": introgressions_chrom18.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "19": introgressions_chrom19.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "20": introgressions_chrom20.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "21": introgressions_chrom21.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        elif introgression.chrom == "22": introgressions_chrom22.append(Sprime(introgression.chrom,introgression.pos,introgression.snp,introgression.nmatch,introgression.dmatch,introgression.popgroup))
        
    introgressions = [introgressions_chrom1,introgressions_chrom2,introgressions_chrom3,introgressions_chrom4,introgressions_chrom5,introgressions_chrom6,introgressions_chrom7,introgressions_chrom8,introgressions_chrom9,introgressions_chrom10,introgressions_chrom11,introgressions_chrom12,introgressions_chrom13,introgressions_chrom14,introgressions_chrom15,introgressions_chrom16,introgressions_chrom17,introgressions_chrom18,introgressions_chrom19,introgressions_chrom20,introgressions_chrom21,introgressions_chrom22]
    with Logging("log.txt"):
        print(f"## Introgressions were found in {len(introgressions)} chromosomes.")
    
    return eqtls, introgressions

def alignSNPs(eqtls, introgressions, alignedSNPs):
    with Logging("log.txt"):
            print(f"## Attempting to align Introgressions with eQTLs.")

    alignments = 0
    baseschecked = 0
    introgressioncount = 0
    for chromosome in range(len(eqtls)):
        with Logging("log.txt"):
            print(f"## Aligning Chromosome {chromosome + 1} out of {len(eqtls)}.")
        for eqtl in eqtls[chromosome]:
            for introgression in introgressions[chromosome]:
                if eqtl.pos == introgression.pos and eqtl.snp == introgression.snp:
                    alignedSNPs.append(SNP(eqtl.chrom, eqtl.pos, eqtl.snp, eqtl.geneID, eqtl.pip, introgression.nmatch, introgression.dmatch, introgression.popgroup))
                    alignments +=1
                introgressioncount += len(introgressions[chromosome])

        baseschecked += len(eqtls[chromosome])
    
    with Logging("log.txt"):
        print(f"## {alignments} Alignments were found out of {introgressioncount} introgressions across {baseschecked} SNPs.")
        print(f"## Alignment Density: {alignments/baseschecked}")
        print(f"## Introgression Enrichment in Searched Genes: {alignments/baseschecked/(introgressioncount/3000000000)}")
            
def outputAlignments(start_time, pip, matchtype, matchin, alignedSNPs):
    header = "Chomosome\tPosition\tSNP\tgeneID\tPIP\tNeanderthal_Match\tDenisovan_Match\tPopulation_Group"
    with open(f"../data/data_clean/Alignments/{start_time}_alignment_{pip}_{matchin}_{matchtype}.tsv","w") as file:
        file.write(f"{header}\n")
        for alignedSNP in alignedSNPs:
            file.write(f"{alignedSNP.chrom}\t{alignedSNP.pos}\t{alignedSNP.snp}\t{alignedSNP.geneID}\t{alignedSNP.pip}\t{alignedSNP.nmatch}\t{alignedSNP.dmatch}\t{alignedSNP.popgroup}\n")

def main():
    start_time = datetime.now()
    with Logging("log.txt"):
        print(f"New Alignment Session")
        print(f"Start Time: {start_time}")

    pip, matchtype, matchin = processArgs()

    eqtls = []
    importeQTL(pip, matchin, eqtls)

    introgressions = []
    importSprime(matchtype, introgressions)
    
    # importDebugger(eqtls, introgressions)

    eqtls, introgressions = splitbyChromosome(eqtls, introgressions)

    alignedSNPs = []
    alignSNPs(eqtls, introgressions, alignedSNPs)

    outputAlignments(start_time, pip, matchtype, matchin, alignedSNPs)

    stop_time = datetime.now()
    runtime = stop_time - start_time
    with Logging("log.txt"):
        print(f"Stop Time: {stop_time}")
        print(f"Total Runtime: {runtime}")
        print(f"")

if __name__ == "__main__":
    main()