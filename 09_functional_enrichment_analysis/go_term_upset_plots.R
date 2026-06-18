# ============================================================
# GO term overlap visualization across conditions
#
# Purpose:
# This script reads GO-term enrichment results from multiple
# Excel sheets, maps them to experimental conditions, and
# generates an UpSet plot to visualize overlaps in GO terms
# across conditions.
#
# Inputs:
# - Excel file containing GO-term results in multiple sheets
# - Each sheet corresponds to one experimental condition group
#
# Outputs:
# - UpSet plot PDF showing GO-term intersections across conditions
#
# Additional Notes:
# - Sheet names are mapped to standardized condition IDs
# - GO terms are de-duplicated per condition
# ============================================================


# 0. Preprocessing and setup ####

# Load required libraries
library(readxl)
library(dplyr)
library(tibble)
library(UpSetR)

# Set input file path
file_path <- "path/to/go_terms_file.xlsx"

# Define desired condition order for plotting
desired_order <- c(
  "sni_chronic", "sni_mid", "sni_acute",
  "cci_chronic", "cci_mid", "cci_acute"
)

# Map Excel sheet names to condition identifiers
sheet_to_condition <- c(
  "SNI Chronic" = "sni_chronic",
  "SNI Mid"     = "sni_mid",
  "SNI Acute"   = "sni_acute",
  "CCI Chronic" = "cci_chronic",
  "CCI Mid"     = "cci_mid",
  "CCI Acute"   = "cci_acute"
)

# Set output directory for results
out_dir <- "path/to/output_folder"


# 1. Read Excel sheet structure ####

# Retrieve all sheet names from input file
sheets <- excel_sheets(file_path)

# Initialize list to store GO terms per condition
go_list <- list()


# 2. Extract GO terms per condition ####

for (sheet in sheets) {
  
  # Log progress
  message("Reading sheet: ", sheet)
  
  # Read current sheet
  df <- read_excel(file_path, sheet = sheet)
  
  # Map sheet name to condition ID
  condition <- sheet_to_condition[[sheet]]
  
  # Extract unique GO terms for condition
  go_list[[condition]] <- unique(df$term_name)
}


# 3. Format GO-term list ####

# Reorder list to match biological/experimental structure
go_list <- go_list[desired_order]


# 4. Generate UpSet plot ####

# Define output PDF path
pdf_path <- file.path(out_dir, "upsetplot_GO_terms.pdf")

# Open PDF device
pdf(pdf_path)

# Generate UpSet plot
print(
  UpSetR::upset(
    fromList(go_list),
    order.by   = "freq",
    keep.order = TRUE,
    sets       = desired_order,
    nsets      = length(go_list)
  )
)

# Close PDF device
dev.off()