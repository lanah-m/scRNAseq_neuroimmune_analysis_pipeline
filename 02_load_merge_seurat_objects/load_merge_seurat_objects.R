# ============================================================
# Multi-batch single-cell RNA sequnecing (scRNA seq) preprocessing and Seurat object construction
#
# Purpose: Loads CellBender-filtered 10X HDF5 files from libraries across 3 batches, converts them to Seurat objects, and merges them into a single dataset. Outputs a merged Seurat object with batch/sample metadata and QC metrics for downstream analysis.
#
# Additional Notes: 3 experimental scRNA seq batches were done. After assessment of the sequencing depth from all of the samples from all of the batches, we noticed that there were issues with sequencing depth with some of the samples. As such, we send those samples with low sequencing depth to be sequenced further. The sequencing information was combined from both runs of these samples at the cellranger stage in the 10X cloud analysis platform. After that was done, cellbender was run and then the samples were analyzed as shown below. 
# ============================================================


# 0. Preprocessing and set up ####
# Load necessary packages quietly
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(rhdf5)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(scater)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(future)
})

# Set parallel processing for speed
plan(multisession, workers = 3)

# Set seed for reproducibility 
set.seed(100)


# 1. User-defined function: load_seurat_from_h5() ####
#
# Purpose:
# Reads 10X Genomics HDF5 data, converts it into a sparse
# matrix, constructs a Seurat object, and attaches metadata.
#
# Inputs:
# - file: path to 10X HDF5 file
# - sample_name: project/sample identifier
# - sampleID: biological condition label
# - batch: sequencing batch label
#
# Output:
# - Seurat object with sampleID and batch metadata
load_seurat_from_h5 <- function(file, sample_name, sampleID, batch) {
  
  # Open HDF5 file
  h5_data <- H5File$new(file, mode = "r")
  
  # Convert 10X HDF5 → sparse matrix
  # This is a workaround because of issues with Seurat's read10XH5() function
  feature_matrix <- Matrix::sparseMatrix(
    i = h5_data[["matrix/indices"]][] ,
    p = h5_data[["matrix/indptr"]][] ,
    x = h5_data[["matrix/data"]][] ,
    dimnames = list(
      h5_data[["matrix/features/name"]][] ,
      h5_data[["matrix/barcodes"]][] 
    ),
    dims = h5_data[["matrix/shape"]][] ,
    index1 = FALSE
  )
  
  # Close file connection 
  h5_data$close_all()
  
  # Ensure gene names are unique (required for Seurat)
  rownames(feature_matrix) <- make.unique(rownames(feature_matrix))
  
  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = feature_matrix,
    min.cells = 3,       # Keep only genes expressed in >= 3 cells
    min.features = 200,  # Keep only cells with >= 200 detected genes
    project = sample_name
  )
  
  # Add metadata
  obj$sampleID <- sampleID   # Biological identity
  obj$batch <- batch         # Sequencing batch
  
  # Set active identity class
  Idents(obj) <- "sampleID"
  
  return(obj)
}


# 2. Define sample metadata ####
samples <- data.frame(
  # sample_name: object name + project label
  sample_name = c(
    "naive_b1",
    "cci_mid_b3",
    "cci_chronic_b3",
    "sni_mid_b3",
    "sni_chronic_b3",
    "cci_acute_b6",
    "sni_acute_b6",
    "cci_chronic_b6",
    "naive_b6"
  ),
  
  # sampleID: biological identity (used for grouping cells later)
  sampleID = c(
    "naive",
    "cci_mid",
    "cci_chronic",
    "sni_mid",
    "sni_chronic",
    "cci_acute",
    "sni_acute",
    "cci_chronic",
    "naive"
  ),
  
  # batch: sequencing batch (important for batch correction later)
  batch = c(
    "batch_1",
    "batch_3",
    "batch_3",
    "batch_3",
    "batch_3",
    "batch_6",
    "batch_6",
    "batch_6",
    "batch_6"
  ),
  
  # file: path to filtered 10X HDF5 output file from cellbender
  file = c(
    "path/to/Naive_B1/Naive_B1_filtered.h5",
    "path/to/CCI_Mid_B3/CCI_Mid_B3_filtered.h5",
    "path/to/CCI_Chronic_B3/CCI_Chronic_B3_filtered.h5",
    "path/to/sni_mid/sni_mid_filtered.h5",
    "path/to/sni_chronic/sni_chronic_filtered.h5",
    "path/to/cci_acute/cci_acute_filtered.h5",
    "path/to/SNI_Acute_B6/SNI_Acute_B6_filtered.h5",
    "path/to/CCI_Chronic_B6/CCI_Chronic_B6_filtered.h5",
    "path/to/naive/naive_filtered.h5"
  ),
  
  stringsAsFactors = FALSE
)


# 3. Read in all samples ####
# Initialize an empty list to store your Seurat objects in
seurat_list <- list()

# Loop through all samples 
for (i in seq_len(nrow(samples))) {
  
  # Write a message to indicate which sample is being processed 
  cat("Loading:", samples$sample_name[i], "\n")
  
  # Use the load_seurat_from_h5() function defined earlier on each sample
  seurat_list[[samples$sample_name[i]]] <- load_seurat_from_h5(
    file = samples$file[i],
    sample_name = samples$sample_name[i],
    sampleID = samples$sampleID[i],
    batch = samples$batch[i]
  )
}


# 4. Merge all Seurat objects as one ####
adata <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1]
)

# Join all the RNA layers which currently split by sample together  
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])


# 5. Add combined metadata ####
# Creates a unique identifier per sample-batch combination 
# Useful for observing batch effects later on
adata$sample_batch <- factor(
  paste(adata$sampleID, adata$batch, sep = "_")
)

# Compute mitochondrial gene percentage
adata[["percent.mito"]] <- PercentageFeatureSet(
  adata,
  pattern = "^mt-"
)


# 6. Inspect and save final merged object ####
# List the number of genes and features for future reference
adata@assays # 27653 features for 82870 cells
 
# Save the object 
saveRDS(adata, "path/to/merged.rds")
