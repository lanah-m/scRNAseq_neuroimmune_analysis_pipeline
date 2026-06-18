# ============================================================
# Data integration for the removal of batch effects
#
# Purpose: This script performs batch integration of a preprocessed scRNA-seq 
# dataset using multiple methods (RPCA, Harmony, and CCA) to compare integration performance.
# ============================================================


# 0. Preprocessing and set up ####
# Load necessary packages quietly
suppressMessages({
  library(Seurat)    
  library(ggplot2)     
  library(sctransform)
})

# Set seed for reproducibility 
set.seed(1234)

# Read the pre-integration clustered rds
adata <- readRDS('path/to/clustered_preintegration.rds')

# Set dims to use as was chosen in the previous script  
dims_to_use <- 23
 

# 1. Prepare Seurat v5 layers ####
# Ensure RNA assay is explicitly a v5 Assay 
adata[["RNA"]] <- as(adata[["RNA"]], Class = "Assay5")

# Join all sample layers prior to splitting
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])

# Now Split by batch so integration is performed *between batches*
adata[["RNA"]] <- split(adata[["RNA"]], f = adata$batch)

# Confirm RNA assay is batch-split
adata[["RNA"]]


# 2.a RPCA Integration ####
adata <- IntegrateLayers(
  object = adata,
  method = RPCAIntegration,
  normalization.method = "SCT",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)

adata <- FindNeighbors(
  adata,
  reduction = "integrated.rpca",
  dims = 1:dims_to_use
)

adata <- FindClusters(
  adata,
  resolution = final_res,
  cluster.name = "rpca_clusters"
)

adata <- RunUMAP(
  adata,
  reduction = "integrated.rpca",
  dims = 1:dims_to_use,
  reduction.name  = "umap.rpca"
)

# 2.b RPCA plots ####
pdf("path/to/rpca_umap.pdf", height = 10, width = 30)

DimPlot(
  adata,
  reduction = "umap.rpca",
  group.by = c("sampleID", "batch", "sample_batch", "rpca_clusters"),
  combine = FALSE,
  cols = my_colors,
  raster = FALSE
)

DimPlot(
  adata,
  reduction = "umap.rpca",
  group.by = "sample_batch",
  split.by = "sampleID",
  combine = FALSE,
  cols = my_colors,
  raster = FALSE
)

dev.off()


# 3.a Harmony Integration ####
adata <- IntegrateLayers(
  object = adata,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  new.reduction = "integrated.harmony",
  verbose = FALSE
)

adata <- FindNeighbors(
  adata,
  reduction = "integrated.harmony",
  dims = 1:dims_to_use
)

adata <- FindClusters(
  adata,
  resolution = final_res,
  cluster.name = "harmony_clusters"
)

adata <- RunUMAP(
  adata,
  reduction = "integrated.harmony",
  dims = 1:dims_to_use,
  reduction.name  = "umap.harmony"
)

# 3.b Harmony plots ####
pdf("path/to/harmony_umap.pdf", height = 10, width = 30)

DimPlot(
  adata,
  reduction = "umap.harmony",
  group.by = c("sampleID", "batch", "sample_batch", "harmony_clusters"),
  combine = FALSE,
  cols = my_colors,
  raster = FALSE
)

DimPlot(
  adata,
  reduction = "umap.harmony",
  group.by = "sample_batch",
  split.by = "sampleID",
  combine = FALSE,
  cols = my_colors,
  raster = FALSE
)

dev.off()


# 4.a CCA Integration ####
adata <- IntegrateLayers(
  object = adata,
  method = CCAIntegration,
  normalization.method = "SCT",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

adata <- FindNeighbors(
  adata,
  reduction = "integrated.cca",
  dims = 1:dims_to_use
)

adata <- FindClusters(
  adata,
  resolution = final_res,
  cluster.name = "cca_clusters"
)

adata <- RunUMAP(
  adata,
  reduction = "integrated.cca",
  dims = 1:dims_to_use,
  reduction.name = "umap.cca"
)

# 4.b CCA plots ####
pdf("path/to/cca_umap.pdf", height = 10, width = 30)

DimPlot(
  adata,
  reduction = "umap.cca",
  group.by = c("sampleID", "batch", "sample_batch", "cca_clusters"),
  combine = FALSE,
  cols = my_colors,
  raster = FALSE
)

DimPlot(
  adata,
  reduction = "umap.cca",
  group.by = "sample_batch",
  split.by = "sampleID",
  combine = FALSE,
  cols = my_colors,
  raster = FALSE
)

dev.off()


# 5. Save processed object ####
saveRDS(adata, "path/to/integrated.rds")