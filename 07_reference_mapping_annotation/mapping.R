# ============================================================
# Single-cell RNA-seq mapping to Renthal DRG atlas
#
# Purpose:
# This script maps a preprocessed single-cell RNA-seq dataset (query) 
# that has yet to be annotated onto a published DRG reference atlas.
# There were two main goals:
# 1) Transfer atlas-based cell type annotations to query data
# 2) Project query cells into atlas UMAP space for visualization
#
# References:
# 1) https://doi.org/10.1126/sciadv.adj9173
# 2) https://doi.org/10.1016/j.neuron.2020.07.026
# 3) Seurat mapping: https://satijalab.org/seurat/articles/integration_mapping
# 4) Anchoring example:
#    https://github.com/Renthal-Lab/harmonized_atlas/blob/main/New_Human_Data_Anchoring/Anchoring.R
# 5) Seurat multimodal mapping:
#    https://satijalab.org/seurat/articles/multimodal_reference_mapping
#
# Additional Notes:
# - SCT assay was intentionally NOT used due to complications
#   and inconsistency with reference workflow (Ref 4, Ref 5 issue discussions).
# - RNA assay with NormalizeData() is used for consistency
#   with reference preprocessing.
# ============================================================


# 0. Preprocessing and setup ####
suppressMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(pheatmap)
  library(grid)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(sceasy)
  library(reticulate)
  reticulate::use_condaenv("wscrna_py38", required = TRUE)
  loompy <- reticulate::import("loompy")
  library(scCustomize)
})

set.seed(100)

setwd("path/to/folder")


# 1. Load query dataset (neurons) ####
obj <- readRDS("path/to/neurons_clustered.rds")

# Use RNA assay (not SCT) for mapping consistency
DefaultAssay(obj) <- "RNA"

# Prepare RNA assay structure
# Join layers 
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# split by batch
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch) # important for normalization by batch

# Confirm split by batch
obj[["RNA"]]

# Normalize query using same strategy as reference (RNA-based normalization)
obj <- NormalizeData(obj)


# 2. Load and preprocess reference atlas ####
atlas <- readRDS("path/to/renthal_science_advances_neurons.Rds")

# Subset to include only the data from mouse 
atlas <- subset(atlas, subset = Species %in% c("Mouse"))

# Convert to Seurat v5-compatible structure
atlas[["RNA"]] <- as(object = atlas[["RNA"]], Class = "Assay5")
atlas <- UpdateSeuratObject(atlas)

# Check if atlas structure as expected
DimPlot(atlas, group.by = "Atlas_annotation")

# Use RNA assay (not SCT) for mapping consistency
DefaultAssay(atlas) <- "RNA"

# Remove injured Atf3+ cluster because the goal is to annotate by neuron 
# identity not injury state 
atlas <- subset(atlas, subset = Atlas_annotation != "Atf3")

# Re-check atlas structure
DimPlot(atlas, group.by = "Atlas_annotation")

# Reconstruct PCA space for reference (required because atlas lacks PCA)
VariableFeatures(atlas@assays$RNA) <- VariableFeatures(atlas@assays$integrated)
atlas <- ScaleData(atlas)
atlas <- RunPCA(atlas, reduction.name = "pca.rna")


# 3. Anchor-based mapping (PCA-based integration) ####
# Find anchors between reference atlas and query
transfer.anchors <- FindTransferAnchors(
  reference = atlas,
  query = obj,
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:30,
  reference.reduction = "pca.rna"
)

# Transfer atlas annotations to query
predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = atlas$Atlas_annotation,
  dims = 1:30
)

# Add predictions to query metadata
obj <- AddMetaData(obj, metadata = predictions)


# 4. Build reference UMAP model for projection ####
DefaultAssay(atlas) <- "integrated"

# Save original atlas UMAP
atlas[["science_advances_pca_refumap"]] <- atlas[["umap"]]

# Recompute UMAP WITH model return for projection
atlas <- RunUMAP(
  atlas,
  dims = 1:30,
  reduction = "pca.rna",
  return.model = TRUE
)


# 5. Map query into atlas UMAP space ####
obj <- MapQuery(
  anchorset = transfer.anchors,
  reference = atlas,
  query = obj,
  refdata = list(celltype = "Atlas_annotation"),
  reference.reduction = "pca.rna",
  reduction.model = "umap"
)


# 6. Visualization of mapping results ####
# Start PDF file
pdf_file <- file.path(
  "path/to/folder",
  "renthal_science_advances_neuron_mapping_using_pca.pdf"
)

pdf(pdf_file, width = 10, height = 10)

# Query mapped onto original CCA space
DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = "predicted.id",
  label = TRUE,
  repel = TRUE
) + ggtitle("Query mapped using PCA-based atlas")

# Query projected into atlas UMAP
DimPlot(
  obj,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  repel = TRUE
) + ggtitle("Query projected into atlas UMAP space")

# Atlas in recomputed UMAP space
DimPlot(
  atlas,
  reduction = "umap",
  group.by = "Atlas_annotation",
  label = TRUE,
  repel = TRUE
) + ggtitle("Atlas (recomputed UMAP)")

# Atlas original UMAP
DimPlot(
  atlas,
  reduction = "science_advances_pca_refumap",
  group.by = "Atlas_annotation",
  label = TRUE,
  repel = TRUE
) + ggtitle("Atlas original UMAP")

# Close PDF
dev.off()


# 7. Save processed objects ####
output_dir <- "path/to/folder"

# Save rds
saveRDS(
  atlas,
  file = file.path(output_dir, "science_advances_pca_mapped.rds")
)

# Remove integrated assay before export (reduces compatibility issues)
DefaultAssay(atlas) <- "RNA"
atlas[["integrated"]] <- NULL

# Save h5ad 
as.anndata(
  x = atlas,
  file_path = output_dir,
  file_name = "science_advances_pca_mapped.h5ad",
  assay = "RNA",
  main_layer = "data",
  other_layers = c("counts"),
  transfer_norm.data = TRUE,
  transfer_dimreduc = TRUE,
  verbose = TRUE
)
