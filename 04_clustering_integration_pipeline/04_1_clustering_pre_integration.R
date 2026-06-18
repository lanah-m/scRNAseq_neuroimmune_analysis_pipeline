# ============================================================
# Clustering and UMAP Visualization
#
# Purpose: This script performs normalization, dimensionality reduction, 
# clustering across multiple resolutions, and UMAP visualization of a QC-filtered 
# scRNA-seq Seurat object. 
# ============================================================


# 0. Preprocessing and set up ####
# Load necessary packages quietly
suppressMessages({
  library(Seurat)    
  library(clustree)    
  library(ggplot2)     
  library(patchwork) 
  library(future)
  library(sctransform)
})

# Set seed for reproducibility 
set.seed(1234)

# Read the qc-filtered rds
adata <- readRDS('path/to/qc-filtered.rds')


# 1. Normalization, feature selection, and dimensionality reduction ####
# Increase memory usage limits so you are not hit with an error
options(future.globals.maxSize = 18700 * 1024^2)  # increase if computer crashes

# Split by batch. In doing do, normalization and variable feature identification 
# is performed for each batch independently such that a consensus set of 
# variable features is identified
adata[["RNA"]] <- split(adata[["RNA"]], f = adata$batch)

# Check that its split by batch 
adata[["RNA"]]

# Normalize, scale, and find variable features using SCTranform() 
# SCTransform() is better than Scale(), FindVariableFeatures(), and Normalize() 
# for removing confounds by technical factors including sequencing depth
adata <- SCTransform(
  adata, 
  vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mito'), 
  vst.flavor = 'v2', 
  return.only.var.genes = FALSE # set return.only.var.genes to false to scale all genes
  ) 

# PCA
adata <- RunPCA(adata, verbose = FALSE)

# Elbow plot to choose number of PCs
ElbowPlot(adata, ndims = 50) 

# Choose the number of dims/PCs
dims_to_use <- 23 # This is your chosen number of PCs/dims based on the ElbowPlot

# Find neighbors using PCs 
adata <- FindNeighbors(adata, dims = 1:dims_to_use, verbose = FALSE)


# 2. Cluster across multiple resolutions ####
# Set resolutions you want to iterate through 
resolutions <- seq(0.0, 0.10, by = 0.02)

# Loop through each resolution to find clusters and create UMAP 
for (res in resolutions) {
  adata <- FindClusters(
    adata,
    resolution = res,
    cluster.name = paste0("merged_res_", res), 
    verbose = FALSE
  )
  
  adata <- RunUMAP(
    adata,
    dims = 1:dims_to_use,
    n.neighbors = 30,
    min.dist = 0.3,
    reduction.name  = "umap.merged",
    verbose = FALSE
  )
}


# 3. Cluster tree visualization ####
# Assess cluster stability across resolutions
clustree(adata, prefix = "merged_res_")

# Chose resolution after inspecting clustree output
final_res <- 0.04 


# 4. Final clustering + UMAP based on chosen resolution ####
# Rerun clustering 
adata <- FindClusters(
  adata,
  resolution = final_res,
  cluster.name = paste0("merged_res_", final_res), 
  verbose = FALSE
)

# Rerun making the UMAP
adata <- RunUMAP(
  adata,
  dims = 1:dims_to_use,
  n.neighbors = 30,
  min.dist = 0.3,
  reduction.name = "umap.merged",
  verbose = FALSE
)


# 5. Ensure clusters do not show obvious QC-related artifacts ####
# Start PDF 
pdf("path/to/file.pdf")

# Plot the QC metrics
FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "percent.mito")

FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(
  adata,
  features = "percent.mito",
  group.by = "seurat_clusters",
  pt.size  = 0.1
) +
  ggtitle("Percent Mitochondrial Content by Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Close the PDF
dev.off() # Clusters did not show strong QC-driven separation, so analysis proceeds


# 6. Plot UMAPs ####
# Start PDF
pdf('path/to/file.pdf', height = 15, width = 45)

# Plot the UMAPs colored by metadata and clustering
DimPlot(
  adata,
  reduction = "umap.merged",
  group.by  = c(
    "sampleID",
    "batch",
    "sample_batch",
    paste0("merged_res_", final_res)
  ),
  combine = FALSE,
  cols    = my_colors,
  raster  = FALSE
)

DimPlot(
  adata,
  reduction = "umap.merged",
  group.by  = "sample_batch",
  split.by  = "sampleID",
  combine   = FALSE,
  cols      = my_colors,
  raster    = FALSE
)

DimPlot(
  adata,
  reduction = "umap.merged",
  label     = TRUE,
  group.by  = c(
    paste0("merged_res_", final_res),
    "sampleID",
    "batch",
    "sample_batch"
  ),
  cols   = my_colors,
  raster = FALSE
)

DimPlot(
  adata,
  reduction = "umap.merged",
  label     = TRUE,
  group.by  = "sampleID",
  split.by  = "sample_batch",
  cols      = my_colors,
  raster    = FALSE
)

# Close PDF
dev.off() # In examining the UMAPs, I saw batch effects and as such proceeded 
# with integration


# 7. Save processed object ####
saveRDS(adata, "path/to/clustered_preintegration.rds")
