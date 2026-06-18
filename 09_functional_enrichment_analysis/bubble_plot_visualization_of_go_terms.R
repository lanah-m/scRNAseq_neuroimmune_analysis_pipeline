# ============================================================
# GO term bubble plot visualization (NF1)
#
# Purpose:
# This script generates GO term bubble plots for NF1 using
# enriched GO terms derived from intersection-based analysis.
# ============================================================


# 0. Setup ####
suppressMessages({
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
})

# Load input data ####
file_path <- "path/to/NF1_GO_terms_by_intersection_final.xlsx"

df <- read.xlsx(file_path)


# 1. User-defined function: plot_go_bubble() ####
#
# Purpose:
# Generates GO term bubble plots from enrichment results.
# The plot visualizes GO term enrichment across condition
# intersections using bubble size (gene count) and color
# (-log10 adjusted p-value).
#
# Inputs:
# - df: data frame containing GO enrichment results with columns:
#       intersection, term_name, avg_intersection_size,
#       avg_adjusted_p_value
# - width: numeric, width of output PDF
# - output_file: full file path for saving the plot
#
# Output:
# - Saves a PDF bubble plot to the specified output path
# - No object is returned (side-effect function)
plot_go_bubble <- function(df, width, output_file) {
  
  pdf(output_file, width = width, height = 8)
  
  ggplot(df, aes(
    x = intersection,
    y = term_name,
    size = avg_intersection_size,
    color = -log10(avg_adjusted_p_value)
  )) +
    geom_point() +
    labs(
      x = "Intersection",
      y = "GO Term",
      size = "Gene count",
      color = "-log10(adj p-value)"
    ) +
    scale_size_continuous(limits = c(0, NA)) +
    scale_color_gradient(low = "darkblue", high = "lightblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  dev.off()
}


# 2. Data preprocessing ####
# Extract celltype from filename
file_name <- basename(file_path)
celltype <- sub("_GO_terms.*", "", file_name)
df$celltype <- celltype

# Make variables numeric so it shows as a continous varibale 
df$avg_intersection_size <- as.numeric(df$avg_intersection_size)
df$avg_adjusted_p_value <- as.numeric(df$avg_adjusted_p_value)

# Check available GO terms
unique(df$term_name)


# 3. Define GO term ordering ####
go_terms_order <- c(
  # Initiation (CCI Acute, SNI Acute)
  "regulation of intracellular transport",
  "establishment of mitochondrion localization",
  
  # Resolving (CCI Mid)
  "neuron maturation",
  "regulation of vesicle fusion",
  "synaptic vesicle membrane organization",
  
  # Resolved (CCI Chronic)
  "endoplasmic reticulum membrane",
  "mitochondrial membrane",
  "positive regulation of endocytosis",
  
  # Chronification (SNI Mid, SNI Chronic)
  "oxidative phosphorylation",
  "ATP synthesis coupled electron transport",
  "respiratory electron transport chain"
)

# Check for mismatches between expected and observed terms
setdiff(unique(df$term_name), go_terms_order)

# Convert to factor and order as listed above
df$term_name <- factor(df$term_name, levels = rev(go_terms_order))


# 4. Define 'intersection' ordering ####
intersection_order <- c(
  "CCI Acute,SNI Acute",
  "CCI Mid",
  "CCI Chronic",
  "SNI Chronic,SNI Mid"
)

# Convert to factor and order as listed above
df$intersection <- factor(df$intersection, levels = intersection_order)


# 5. Generate plots ####
output_dir <- "path/to/folder"

plot_go_bubble(
  df,
  width = 6,
  output_file = file.path(output_dir, "go_term_bubble_plot_nf1.pdf")
)
