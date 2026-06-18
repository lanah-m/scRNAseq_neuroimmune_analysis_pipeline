# ============================================================
# Reading SCCAF clusters into R 
#
# Purpose: This script loads an existing CCA-integrated Seurat
# object, imports SCCAF-refined cluster assignments, attaches
# them to the object metadata, and saves an updated Seurat 
# object for downstream analysis and annotation.
# ============================================================


# 0. Preprocessing and set up ####
# Load required packages quietly
suppressMessages({
  library(Seurat)
  library(ggplot2)
})

# Set seed for reproducibility 
set.seed(1234)

# Load CCA-integrated Seurat object
adata <- readRDS( "path/to/allcells_cca_integrated.rds")

# Load SCCAF cluster assignments
adata_sccaf <- read.csv(
  "path/to/allcells_r0.01-0.03_cln.csv",
  header = FALSE, # File has no header
  row.names = 1   # row names correspond to cell barcodes
)


# 1. Attach SCCAF metadata ####
# Add SCCAF cluster labels to Seurat object metadata
adata[["SCCAF"]] <- adata_sccaf

# Set SCCAF clusters as active identity class
Idents(adata) <- "SCCAF"


# 2. Double check the SCCAF metadata added as expected ####
# Visualize SCCAF clusters on UMAP
DimPlot(
  adata,
  reduction = "umap",
  label = TRUE,
  raster = FALSE
)

DimPlot(
  adata,
  reduction = "umap",
  group.by = c("SCCAF", "sampleID"),
  label = TRUE,
  raster = FALSE
)


# 3. Save processed object ####
# Save updated Seurat object with SCCAF annotations
saveRDS(adata, file = "path/to/allcells_cca_integrated_cln.rds")
