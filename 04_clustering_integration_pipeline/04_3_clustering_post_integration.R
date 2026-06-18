# ============================================================
# Clustering after integration (CCA-selected)
#
# Purpose: This script performs downstream analysis of a CCA-integrated scRNA-seq 
# dataset, including clustering, resolution selection, and export to 
# RDS and AnnData formats for downstream use.
# ============================================================


# 0. Preprocessing and set up ####
# Load necessary packages quietly
suppressMessages({
  library(Seurat)    
  library(ggplot2)     
  library(sctransform)
  library(sceasy)
  library(reticulate) 
  reticulate::use_condaenv("wscrna_py38", required = TRUE) 
  loompy <- reticulate::import('loompy') 
  library(scCustomize)
})

# Set seed for reproducibility 
set.seed(1234)

# Read the pre-integration clustered rds
adata <- readRDS('path/to/integrated.rds')

# Set dims to use as was chosen in a previous script  
dims_to_use <- 23

# Copy the content of 'umap.cca' into 'umap'and integrated.cca into pca
adata@reductions$pca <- adata@reductions$`integrated.cca`
adata@reductions$umap <- adata@reductions$`umap.cca`


# 1. Keep only CCA-related reductions to avoid ambiguity downstream ####
# CCA was chosen as it showed the best batch correction
keep_reductions <- c(
  "integrated.cca",
  "umap.cca",
  "umap",
  "pca"
)

for (r in names(adata@reductions)) {
  if (!(r %in% keep_reductions)) {
    adata@reductions[[r]] <- NULL
  }
}


# 2. Explore multiple clustering resolutions ####
resolutions <- seq(0.00, 0.10, by = 0.01)

for (res in resolutions) {
  adata <- FindClusters(
    adata,
    resolution = res,
    verbose = FALSE,
    cluster.name = paste0("cca_res_", res)
  )
  
  adata <- RunUMAP(
    adata,
    reduction = "integrated.cca",
    dims = 1:dims_to_use,
    n.neighbors = 30,
    min.dist = 0.3,
    reduction.name = "umap.cca",
    verbose = FALSE
  )
}


# 3. Assess cluster stability ####
# Visualize how clusters change across resolutions
clustree(adata, prefix = "cca_res_")

# Chosen resolution based on its stability from the clustree() inspection
final_res <- 0.02  


# 4. 'Final' clustering (this will later be optimized with SCCAF) ####
adata <- FindClusters(
  adata,
  resolution = final_res,
  verbose = FALSE,
  cluster.name = paste0("cca_res_", final_res)
)

adata <- RunUMAP(
  adata,
  reduction = "integrated.cca",
  dims = 1:dims_to_use,
  n.neighbors = 30,
  min.dist = 0.3,
  reduction.name = "umap.cca",
  verbose = FALSE
)
# NOTE: After running UMAP, the following warning was observed: "Key ‘umapcca_’ 
# taken, using ‘maujg_’ instead" indicating Seurat auto-renamed the embedding columns.


# 5. Sanity check plots ####
# This plot will show the CCA UMAP 
DimPlot(
  adata,
  reduction = "umap.cca",
  group.by  = c(
    "sampleID",
    "batch",
    "sample_batch",
    paste0("cca_res_", final_res)
  ),
  combine = FALSE
)

# This plot should also show the CCA UMAP. If not, there is a problem
DimPlot(
  adata,
  group.by = "sample_batch",
  combine  = FALSE
)


# 6. Save processed object ####
# Switch back to parallel for speed
plan(multisession, workers = 3)

# Save as RDS
saveRDS(adata, file = 'path/to/allcells_cca_integrated.rds')

# Save the object directly from seurat to anndata to avoid gene drop out
as.anndata(
  x = adata,
  file_path = "path/to/folder",
  file_name = "allcells_cca_integrated.h5ad",
  assay = "SCT",                       # Make sure you're using the correct assay
  main_layer = "data",                 # Use normalized data as main layer
  other_layers = c("counts"),          # Also include raw counts
  transfer_norm.data = TRUE,           # Explicitly transfer norm.data layer
  transfer_dimreduc = TRUE,            # Include UMAP, PCA, etc.
  verbose = TRUE                       # Show progress
)

# Save the SCTransform() identified highly variable features to a file
hvg_genes <- adata@assays[["SCT"]]@var.features

write.csv(
  hvg_genes, 
  file = "path/to/sct_hvg.csv", 
  row.names = FALSE
)
