# Load necessary libraries
library(tidyverse)
library(viridis)  # Ensure viridis is available for colorblind-safe palettes

# Load the data
file_path <- "2024-10-30 12_27_25.625848_alignment_0.0_all_either.tsv"
data <- read.delim(file_path)

# Filter data for variantPIP > 0.1
filtered_data <- data %>%
  filter(variantPIP > 0.1)

# Split filtered data into Neanderthal (NMATCH == "match") and Denisovan (DMATCH == "match")
neanderthal_data <- filtered_data %>%
  filter(NMATCH == "match")

denisovan_data <- filtered_data %>%
  filter(DMATCH == "match")

# Check data after filtering
print(neanderthal_data)
print(denisovan_data)

# Function to count unique occurrences and plot
create_plot <- function(dataset, title, filename) {
  # Check if the dataset is valid
  if (nrow(dataset) == 0) {
    stop("Dataset is empty, no data to plot!")
  }
  
  variant_pop_unique <- dataset %>%
    group_by(POPGROUP) %>%
    summarise(unique_variants = n_distinct(variant_rsID), .groups = "drop")
  
  # Plot
  plot <- ggplot(variant_pop_unique, aes(x = reorder(POPGROUP, -unique_variants), y = unique_variants, fill = POPGROUP)) +
    geom_bar(stat = "identity") +
    labs(
      x = "Population Group", 
      y = "Number of SNPs",
      title = title
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16)  # Center and increase the title size
    ) +
    scale_fill_viridis_d(option = "cividis")  # Using colorblind-safe palette
  
  # Save plot
  ggsave(filename = filename, plot = plot, width = 8, height = 5, dpi = 300)
  
  return(plot)
}

# Create and save plots for Neanderthal and Denisovan data (filtered for matches and variantPIP > 0.1)
neanderthal_plot <- create_plot(neanderthal_data, 
                                "Notably Introgressed Neanderthal SNPs per Population Group", 
                                "neanderthal_plot_variantPIP.png")

denisovan_plot <- create_plot(denisovan_data, 
                              "Notably Introgressed Denisovan SNPs per Population Group", 
                              "denisovan_plot_variantPIP.png")

# Display the plots
print(neanderthal_plot)
print(denisovan_plot)
