# ============================================================
# Quality Control (QC) Filtering
#
# Purpose: Perform scRNA-seq QC including doublet and statistcal threshold- 
# based filtering to generate a Seurat object with only high-quality cells.
# ============================================================


# 0. Preprocessing and set up ####
suppressPackageStartupMessages({
  library(Seurat)                 
  library(SingleCellExperiment)   
  library(scDblFinder)          
  library(scater)                
  library(BiocParallel)          
  library(Matrix)                
  library(dplyr)                
  library(ggplot2)          
  library(future)
})

# Set parallel processing for speed
plan(multisession, workers = 3)

# Set seed for reproducibility 
set.seed(100)

# read in the object
adata <- readRDS('path/to/merged.rds')


# 1. Define QC thresholds BEFORE doublet removal ####
# We define QC outliers before removing any doublets via scDblFinder because they 
# are based on statistical thresholds that will change once cells are removed 

# Extract metadata from Seurat object
qc <- adata@meta.data  

# Detect outliers using 5 MAD as a threshold 
qc$umi_outlier <- isOutlier(qc$nCount_RNA, nmads = 5, type = "both")
qc$feature_outlier <- isOutlier(qc$nFeature_RNA, nmads = 5, type = "both")
qc$mito_outlier <- isOutlier(qc$percent.mito, nmads = 5, type = "both")

# Write QC flags back to Seurat object
adata@meta.data <- qc


# 2. Doublet detection (scDblFinder) ####
# We now switch to doublet removal via scDblFinder rather than removing the 
# flagged outliers first becuase the default expected doublet rate is calculated 
# on the basis of the cells given, and if you excluded a lot of cells as low 
# quality, scDblFinder might think that the doublet rate should be lower than it is.

# Convert Seurat to SingleCellExperiment object as that is what scDblFinder requires
sce <- as.SingleCellExperiment(adata)

# Ensure reproducibility 
bp <- SnowParam(workers = 3, RNGseed = 1234)

# Run doublet finder 
# Since both the cluster-based or not cluster based doublet identificaton 
# approach of scDblFinder perform very similarly overall in benchmarks, 
# I just went with the default (not cluster based).
sce <- scDblFinder(
  sce,
  samples = "batch",
  BPPARAM = bp 
  # You can also use dbr.per1k =  to set mutliplet rate
) 

# Inspect estimated doublet assignments
# If values don't makes sense, go back and adjust dbr.per1k paramater of scDblFinder()
table(sce$scDblFinder.class) 
# Singlet 68335, doublet 14535


# 3. Transfer scDblFinder results back to Seurat ####
# Extract metadata from SingleCellExperiment object
meta_scdblfinder <- as.data.frame(sce@colData@listData)

# Keep only relevant scDblFinder outputs
meta_scdblfinder <- dplyr::select(meta_scdblfinder, starts_with("scDblFinder"))

# Check to make sure the first few rows look as expected 
head(meta_scdblfinder)

# Ensure rownames match cell barcodes
rownames(meta_scdblfinder) <- sce@colData@rownames

# Check again
head(meta_scdblfinder)

# Add metadata back into Seurat object
adata <- AddMetaData(object = adata, metadata = meta_scdblfinder)


# 4. Visualize doublet classification ####
# Visualize QC metrics stratified by scDblFinder classification
# to confirm that doublets correspond to expected high RNA/feature signals.

# Start PDF
pdf("path/to/file.pdf", width = 11, height = 6)

# Create the violin plots
VlnPlot(
  adata,
  group.by = "sampleID",
  split.by = "scDblFinder.class",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0
) +
  theme(legend.position = "right")

# Close the PDF
dev.off()


# 5. Remove doublets ####
adata <- subset(adata, subset = scDblFinder.class == "singlet")

# Check how many cells remain
adata@assays
# 27653 features, 68335 cells 


# 6. Visualizing remaining outliers to make sure it looks as expected ####
# Now that doublets have been identified and removed with scDblFinder
# we can switch back to visualizing and removing the outliers we identified 
# previously with isOutlier()

# Start the PDF
pdf("path/to/file.pdf", width = 11, height = 11)

# Mitochondrial % vs UMI counts plot
ggplot(
  adata@meta.data,
  aes(
    x = nCount_RNA,
    y = percent.mito,
    color = umi_outlier | feature_outlier | mito_outlier
  )
) +
  geom_point(alpha = 0.6, size = 1) +
  facet_wrap(~sampleID, scales = "free") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(
    title = "% mito vs nCount_RNA",
    x = "nCount_RNA (UMIs)",
    y = "% Mitochondrial",
    color = "Outlier"
  ) +
  theme_minimal()


# Genes vs UMIs plot
ggplot(
  adata@meta.data,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    color = umi_outlier | feature_outlier | mito_outlier
  )
) +
  geom_point(alpha = 0.6, size = 1) +
  facet_wrap(~sampleID, scales = "free") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(
    title = "nCount_RNA vs nFeature_RNA",
    x = "nCount_RNA (UMIs)",
    y = "nFeature_RNA (genes)",
    color = "Outlier"
  ) +
  theme_minimal()

# Close the PDF
dev.off()


# 7. Remove the remaining outliers ####
# Identify high-quality cells by removing any flagged as outliers 
cells_to_keep <- rownames(qc)[
  !(qc$umi_outlier | qc$feature_outlier | qc$mito_outlier)
]

# Remove all low quality cells
adata <- subset(adata, cells = cells_to_keep)


# 8. Inspect and save final qc-filtered object ####
# List the number of genes and features remaining after qc for future reference
adata@assays # 27653 features, 59037 cells 

# Save the object
saveRDS(adata, 'path/to/qc-filtered.rds')
