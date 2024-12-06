# Load necessary libraries
library(tidyverse)
library(viridis)  # Ensure viridis is available for colorblind-safe palettes

# Load the data
file_path <- "2024-10-30 12_27_25.625848_alignment_0.0_all_either.tsv"
data <- read.delim(file_path)

# Split data into Neanderthal (NMATCH == "match") and Denisovan (DMATCH == "match") based on NMATCH and DMATCH
neanderthal_data <- data %>%
  filter(NMATCH == "match")

denisovan_data <- data %>%
  filter(DMATCH == "match")

# Function to count unique occurrences and plot
create_plot <- function(dataset, title, filename) {
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
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  # Center and make title larger
    ) +
    scale_fill_viridis_d(option = "plasma")  # Using colorblind-safe palette
  
  # Save plot
  ggsave(filename = filename, plot = plot, width = 8, height = 5, dpi = 300)
  
  return(plot)
}

# Create and save plots for Neanderthal and Denisovan data (filtered for matches)
neanderthal_plot <- create_plot(neanderthal_data, 
                                "Introgressed Neanderthal SNPs per Population Group", 
                                "neanderthal_plot_all.png")

denisovan_plot <- create_plot(denisovan_data, 
                              "Introgressed Denisovan SNPs per Population Group", 
                              "denisovan_plot_all.png")

# Display the plots
print(neanderthal_plot)
print(denisovan_plot)
