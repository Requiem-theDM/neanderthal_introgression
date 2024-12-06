# Load necessary libraries
library(tidyverse)
library(viridis)  # Ensure viridis is available for colorblind-safe palettes

# Load the data
file_path <- "2024-10-30 12_27_25.625848_alignment_0.0_all_either.tsv"
data <- read.delim(file_path)

# Count unique occurrences of variant_rsID for each POPGROUP
variant_pop_unique <- data %>%
  group_by(POPGROUP) %>%
  summarise(unique_SNPs = n_distinct(variant_rsID), .groups = "drop")

# Inspect the summary
print(variant_pop_unique)

# Plot the unique SNPs associated with each population
ggplot(variant_pop_unique, aes(x = reorder(POPGROUP, -unique_SNPs), y = unique_SNPs, fill = POPGROUP)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Population Group", 
    y = "Number of SNPs",
    title = "Introgressed SNPs of Different Population Groups"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)  # Center title and make it bigger
  ) +
  scale_fill_viridis_d(option = "plasma")  # Using 'plasma' for a visually appealing, perceptible palette

