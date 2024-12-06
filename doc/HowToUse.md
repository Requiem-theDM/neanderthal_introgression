## Steps to Reproduce Findings

1. Clone this git repository. 
2. Download all files from the following links to the data/data_raw directory:
- [eQTL Files](https://www.dropbox.com/scl/fo/3filfi60qhxgokmjb6weq/AGLqwdegrXF5NUYqAl7EI6Q/QTL_results/eQTL_results/eQTL_finemapping_results?dl=0&rlkey=ihnv8hpbj1ilb96cxsm8fc42p&subfolder_nav_tracking=1)
- [SPrime Data](https://data.mendeley.com/datasets/y7hyt83vxr/1)

3. You should have two .zip folders in your data_raw directory.

4. In terminal, navigate to analysis and run the following command:
- `./processArchives.sh`

5. In terminal, run the following command:
- `./processSprimes.sh`

6. In terminal, run the following commands to generate datasets for downstream analysis:
- `./alignSNPs.py 0 e all`
- `./alignSNPs.py 0 e all -ni`
- alignSNPs takes the following parameters: desired PIP search value, desired match type, and desired search file, **in that order.**
    - PIP search value defaults to 0.5, and accepts any number
    - Match type defaults to Neanderthals and accepts any of the following in a case-insensitive manner:
        - Neanderthal match type: n, nmatch, neanderthal
        - Denisovan match type: d, dmatch, denisovan
        - Both match types: b, bmatch, both
        - Either match type: e, ematch, either
    - Search file defaults to eQTL_finemapping.significantAssociations.MAGE.v1.0.txt and accepts either of the following in a case-insensitive manner:
        - SNPs with Significant Associations: sig
        - SNPs with Any Association: all, any
    - `-ni` can be optionally specified to perform an inverted search and look for SNPs that were not found in the specified introgression group.
- alignSNPs prints logs of each run to log.txt, stored in the analysis folder.
- alignSNPs outputs matched SNPs stamped with search start date, start time, and parameters in a tab separated values format in the data/data_clean/Alignments directory.

## Using R script

1. Once alignments have been generated, open R files (Number_of_SNPs_per_POPGROUP_Neanderthals_and_Denisovan(SIGNIFICANTVARIANTS).R,
Number_of_SNPs_per_POPGROUP_Neanderthals_and_Denisovan.R, UniqueSNPsAssociatedwithEachPopulationGroup(AllVariants).R,
UniqueSNPsAssociatedwithEachPopulationGroup(SignificantVariants).R) within RStudio
2. Ensure that the file used for file_path matches the file you generated using ./alignSNPs.py 0 e all within your data_clean directory
   - This file should be the ONLY the introgressed data (no non-introgressed data needs to be loaded.)
3. Run R code line by line to generate plots

## Identifying SNPs of Interest for Further Analysis

1. Open IDSNPsForLitSearch.R in R Studio.
    - Ensure file path matches file generated and working directory is accurately set. 
    - The code currently is set to identify SNPs matching both Denisovans and Neanderthals with significant associations. 
    - Adjust the PIP threshold if desired (currently set at 0.4). 
3. Run all lines to identify SNPs of interest from the dataset.
