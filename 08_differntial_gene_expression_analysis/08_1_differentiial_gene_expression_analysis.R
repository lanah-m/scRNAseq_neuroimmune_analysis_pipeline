# ============================================================
# Differential expression across conditions (cell-type resolved)
#
# Purpose:
# This script performs differential expression (DE) analysis
# comparing experimental conditions against a "naive" baseline
# within each cluster. More specifically, clusters are split by 
# condition such that and within each cluster, each condition is 
# compared to its naive state. Combined DE table across all cell 
# types and conditions were generated. 
#
# Additional Notes:
# - RNA assay is used (NOT SCT) for consistency with prior steps
# - Positive logFC = upregulated in condition vs naive
# - Negative logFC = downregulated in condition vs naive
# ============================================================


# 0. Preprocessing and setup ####
suppressMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(openxlsx)
})

set.seed(1234)

# Load input data
adata <- readRDS(
  "path/to/allcells_cca_integrated_cln.rds"
)


# 1. RNA assay normalization (consistent with previous workflows) ####
DefaultAssay(adata) <- "RNA"

# Join and split RNA layers for batch-aware normalization
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])
adata[["RNA"]] <- split(adata[["RNA"]], f = adata$batch)

# Check that it is split by batch
adata[["RNA"]]

# Normalize using RNA-based approach (NOT SCT)
adata <- NormalizeData(adata)

# Re-join layers after normalization
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])


# 2. Define identity structure for DE analysis ####
# Check what you clusters are at the moment
unique(adata@meta.data[["SCCAF"]])

# Create combined identity: cell type × condition
adata$celltype.condition <- paste(adata$SCCAF, adata$sampleID, sep = "_")
Idents(adata) <- "celltype.condition"


# 3. Differential expression analysis (condition vs naive) ####
cell_types <- unique(adata$SCCAF)
conditions <- unique(adata$sampleID)

de_results <- list()

for (ct in cell_types) {
  
  # Skip if naive baseline does not exist for this cell type
  if (!paste(ct, "naive", sep = "_") %in% Idents(adata)) next
  
  for (cond in conditions) {
    
    # Skip naive vs naive comparison
    if (cond == "naive") next
    
    ident1 <- paste(ct, cond, sep = "_")     # non-naive condition
    ident2 <- paste(ct, "naive", sep = "_")  # naive
    
    # Run DE analysis
    de <- FindMarkers(
      adata,
      ident.1 = ident1,
      ident.2 = ident2,
      assay = "RNA",
      verbose = FALSE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    # Keep only significant genes
    de <- de[de$p_val_adj <= 0.05, ]
    
    # Skip empty results
    if (nrow(de) == 0) next
    
    # Add annotation columns
    de$cell_type <- ct
    de$condition <- cond
    de$gene <- rownames(de)
    
    # Store result
    de_results[[paste(ident1, "vs", ident2, sep = "_")]] <- de
  }
}

# Combine results into single table
all_de <- bind_rows(de_results)


# 4. Export results ####
write.xlsx(
  all_de,
  file = "path/to/condition_vs_naive_DEG_by_celltype.xlsx",
  rowNames = FALSE
)
