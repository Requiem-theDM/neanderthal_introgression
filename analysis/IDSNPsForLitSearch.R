# Load required libraries
library(tidyverse)

# Load the TSV file
setwd("~/qb24_project/data/data_clean/Alignments/")
file_path <- "2024-10-30 12:24:57.975162_alignment_0.0_sig_both.tsv"
df <- read_tsv(file_path)

#Select just variantPIP and variant_rsID 
df_subset <- select(df, variant_rsID, variantPIP)


#Remove duplicates
df_no_dup <- distinct(df_subset, variant_rsID, variantPIP)
print(df_no_dup$variant_rsID)

#The SNP variants present in the dataset have rsIDs rs7479832, rs3210908,
#rs6421977, rs6421976, and rs6598024  


#Select for variants that have a PIP score greater than 0.4
df_filtered <- filter(df_no_dup, variantPIP > 0.4)
print(df_filtered$variant_rsID)

#Only rs7479832 has a PIP value greater than 0.4 (0.4468349)







