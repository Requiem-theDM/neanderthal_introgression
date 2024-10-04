# Reviving the Past: A Replication of Neanderthal DNA's Role in Human Gene Expression
## Josh Beale, Caroline Moore, Fijare Plous

## Project Description
Our project aims to replicate Figure 4 from McCoy et al., 2017, titled “Impacts of Neanderthal-introgressed sequences on the landscape of human gene expression.” This figure visualizes purifying selection on the gradual movement of genes from one species into the gene pool of another, via the hybridization of humans and Neanderthals. The analysis within the figure looks at allele frequency classification between Neanderthal genes that have moved into the human population and those that have not  Ultimately McCoy et. al. found no significant differences in allele-specific expression between introgressed and non-introgressed variants.

## Figure Example
- https://www.ncbi.nlm.nih.gov/pmc/articles/instance/6219754/bin/nihms-848864-f0004.jpg 

## Datasets
- https://github.com/mccoy-lab/MAGE/tree/main
- https://data.mendeley.com/datasets/y7hyt83vxr/1

## Software Requirements
- Python3
- R version 4.4.1
- GNU bash, version 3.2.57(1)-release

## Proposed Workflow
1. Align credible sets from finemapping MAGE eQTL based on varying PIP thresholds with SNPs from Sprime dataset based on if the Sprime matches either neanderthals or denisovans.
2. Generate figures for neanderthal matches
- Derived Allele Frequency vs. Number of SNPs
- Significant differences in the estimated proportions of introgressed (pN) and non-introgressed (pH) variants showing significant ASE (at 10% FDR), stratified by derived allele frequency
3. Generate figures for Denisovan matches
- Derived Allele Frequency vs. Number of SNPs
- Significant differences in the estimated proportions of introgressed (pN) and non-introgressed (pH) variants showing significant ASE (at 10% FDR), stratified by derived allele frequency


```
Repository Structure
.
├── README.md
├── analysis            <- all things data analysis
├── comm                <- internal communication such as meeting notes
├── data
│   ├── data_clean      <- clean version of the data
│   └── data_raw        <- raw data (don't touch) (Or links to obtaining it)
├── diss
│   └── presentations   <- end of semester presentation
├── doc                 <- documentation
└── misc                <- miscellaneous files that don't fit elsewhere  
```