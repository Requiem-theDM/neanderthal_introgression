# Load required libraries
library(tidyverse)
library(ggthemes)
library(cowplot)
# Load the TSV file
setwd("~/qb24_project/data/data_clean/Alignments/")
file_path <- "2024-10-30 12:27:25.625848_alignment_0.0_all_either.tsv"
introgressions = read_tsv(file_path)
non_introgressed = read_tsv("2024-10-30 12:25:38.798494_alignment_0.0_sig_either_mismatch.tsv")
# Define derived allele frequency bins (1% bins)
bins <- seq(0, 1, by = 0.01)  # Create bins from 0 to 1 with a width of 0.01 (1% bins)
labels <- paste0(head(bins, -1), "-", tail(bins, -1))  # Generate labels for the bins
# Filter data for introgressions unique to denisovans (DMATCH = "match" & NMATCH != "match"),
# introgressions unique to neanderthals (NMATCH = "match" & DMATCH != "match"),
# and introgressions from both (DMATCH = "match" and NMATCH = "match")
denisovans <- introgressions %>%
  filter(NMATCH != "match" & DMATCH == "match") %>%
  mutate(Category = "Denisovans") %>% subset(select = -c(NMATCH,DMATCH))
neanderthals <- introgressions %>%
  filter(NMATCH == "match" & DMATCH != "match") %>%
  mutate(Category = "Neanderthals") %>% subset(select = -c(NMATCH,DMATCH))
both <- introgressions %>%
  filter(NMATCH == "match" & DMATCH == "match") %>%
  mutate(Category = "Both") %>% subset(select = -c(NMATCH,DMATCH))
non_introgressed = non_introgressed %>% mutate(POPGROUP = "N/A") %>% mutate(Category = "Neither")
# Combine data
combined_df = bind_rows(denisovans, neanderthals, both, non_introgressed)

inset = ggplot(filter(combined_df,Category=='Both',variantPIP>0.1)) +
  geom_histogram(aes(x=variantPIP),position='dodge',binwidth=0.01) +
  scale_y_continuous() +
  labs(x="Likelihood of Being Causative",
       y="Number of SNPs") +
  theme_base(base_size = 12) +
  scale_fill_manual("black") 

all = ggplot(combined_df) +
  geom_histogram(aes(x=variantPIP,fill = Category),position='dodge',binwidth=0.01) +
  labs(x="Likelihood of Being Causative",
       y="Number of SNPs") +
  theme_base(base_size = 12) +
  scale_fill_colorblind()
  

ggdraw(all + theme_half_open(12)) +
  draw_plot(inset, .25, .25, .5, .5)

           
       