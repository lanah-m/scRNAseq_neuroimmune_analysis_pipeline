# ============================================================
# Annoation of Clusters via Cluster DEGs and known cell markers
#
# Purpose: This script identifies conserved cell-type marker
# genes across batches using FindConservedMarkers() on the
# RNA assay. This approach is used as an alternative to
# PrepSCTFindMarkers() due to errors arising from multiple
# SCT models with unequal library sizes.
#
# # Additional Notes: 
# - Seurat issue: https://github.com/satijalab/seurat/issues/9130
# - Seurat integration vignette:
#   https://satijalab.org/seurat/articles/integration_introduction.html
# ============================================================


# 0. Preprocessing and setup ####
# Load required libraries quietly
suppressMessages({
  library(Seurat)     
  library(openxlsx)  
  library(future)   
  library(metap)    
  library(multtest)    
})

# Set parallelization strategy (sequential for stability)
plan(future::sequential)

# Increase future global size to avoid memory errors
options(future.globals.maxSize = 18700 * 1024^2)

# Set output directory
output_dir <- "path/to/folder"

# Load input data (CCA-integrated, SCCAF-annotated Seurat object)
adata <- readRDS("path/to/allcells_cca_integrated_cln.rds")

# Set default assay to RNA for DE analysis
DefaultAssay(adata) <- "RNA"


# 1. User-defined function: plot_feature_pdf() ####
#
# Purpose:
# Generate FeaturePlot visualizations for predefined marker
# sets and export them as a multi-page PDF.
#
# Inputs:
# - adata: Seurat object
# - marker_list: Named list of genes
# - pdf_path: Output PDF file path
#
# Output:
# - PDF file containing FeaturePlot panels
plot_feature_pdf <- function(adata, marker_list, pdf_path) {
  pdf(pdf_path)
  for (genes in marker_list) {
    FeaturePlot(adata, features = genes, combine = FALSE)
  }
  dev.off()
  cat("Plots saved at:", pdf_path, "\n")
}


# 2. RNA assay preparation for conserved DE ####
# Join all RNA layers prior to splitting
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])

# Split RNA assay by batch so each batch is normalized independently
adata[["RNA"]] <- split(adata[["RNA"]], f = adata$batch)

# Confirm RNA assay is split by batch
adata[["RNA"]]

# Normalize RNA data (creates 'data' layer used for DE)
adata <- NormalizeData(adata)

# Re-join RNA layers after normalization
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])


# 3. Identify conserved markers per cell type ####
# Retrieve unique SCCAF cell identities
idents <- unique(adata@meta.data[["SCCAF"]])

# Initialize list to store markers during loop
allcells.markers.list <- list()

# For each SCCAF cluster, identify conserved markers
for (ident in idents) {
  message("Processing cell type: ", ident)
  
  markers <- FindConservedMarkers(
    adata,
    ident.1        = ident,
    grouping.var   = "batch",
    assay          = "RNA",
    slot           = "data",
    only.pos       = TRUE,
    min.pct        = 0.85,
    logfc.threshold = 1,
    verbose        = FALSE
  )
  
  # Add metadata columns
  markers$celltype <- ident
  markers$gene     <- rownames(markers)
  
  # Compute min/max adjusted p-values across batches
  adj_cols <- grep("p_val_adj$", colnames(markers), value = TRUE)
  markers$min_adj_pval <- apply(markers[, adj_cols], 1, min, na.rm = TRUE)
  markers$max_adj_pval <- apply(markers[, adj_cols], 1, max, na.rm = TRUE)
  
  # Compute min/max log2 fold-change across batches
  logfc_cols <- grep("avg_log2FC$", colnames(markers), value = TRUE)
  markers$min_log2FC <- apply(markers[, logfc_cols], 1, min, na.rm = TRUE)
  markers$max_log2FC <- apply(markers[, logfc_cols], 1, max, na.rm = TRUE)
  
  # Compute min/max pct.1 across batches
  pct1_cols <- grep("pct\\.1$", colnames(markers), value = TRUE)
  markers$min_pct1 <- apply(markers[, pct1_cols], 1, min, na.rm = TRUE)
  markers$max_pct1 <- apply(markers[, pct1_cols], 1, max, na.rm = TRUE)
  
  # Compute min/max pct.2 across batches
  pct2_cols <- grep("pct\\.2$", colnames(markers), value = TRUE)
  markers$min_pct2 <- apply(markers[, pct2_cols], 1, min, na.rm = TRUE)
  markers$max_pct2 <- apply(markers[, pct2_cols], 1, max, na.rm = TRUE)
  
  # Store results
  allcells.markers.list[[ident]] <- markers
}


# 4. Post-processing and filtering of DEGs ####
# Keep only summary and identifier columns
allcells.markers.list <- lapply(allcells.markers.list, function(df) {
  df[, c(
    "celltype",
    "gene",
    "min_log2FC",
    "max_log2FC",
    "min_adj_pval",
    "max_adj_pval",
    "min_pct2",
    "max_pct2",
    "min_pct1",
    "max_pct1"
  )]
})

# Combine all cell types into one data frame
adata.markers <- do.call(rbind, allcells.markers.list)

# Check it looks as expected
head(adata.markers)

# Keep only the genes significant in all batches
adata_sig <- adata.markers[adata.markers$max_adj_pval < 0.05, ]

# Remove mitochondrial and ribosomal genes
adata_sig <- adata_sig[!grepl("^mt-", adata_sig$gene, ignore.case = TRUE), ]
adata_sig <- adata_sig[!grepl("^Rp[sl]", adata_sig$gene, ignore.case = TRUE), ]

# Confirm final structure
head(adata_sig)

# Export conserved DEGs to Excel
write.xlsx(
  adata_sig,
  file = file.path(
    output_dir,
    "allcells_conserved_deg_log2fc1_minpct0.85.xlsx"
  )
)


# 5. Plot the expression of known markers ####  
# Define marker sets
marker_sets <- list(
  Sex = c("Xist", "Eif2s3y"),
  Endothelial = c("Cldn5", "Plvap", "Emcn", "Flt1", "Egfl7", "Pecam1", "Cd31", "Ly6c1", "Esam"),
  Mesenchymal = c("Cd34", "Apod", "Pdgfra"),
  Pericytes = c("Notch3", "Pdgfrb", "Acta2", "Kcnj8"),
  Fibroblasts = c("Dcn", "Col1a1", "Ptgds", "Apod", "Fn1", "Col3a"),
  Satellite_Glia = c("Fabp7", "Ednrb", "Cdh19", "Plp1"),
  RBC_Progenitors = c("Hbb-bt", "Hba-a2", "Hba-a1").
  Neurons = c("Rbfox3", "Isl1", "Tubb3", "Gal", "Tac", "Prph"),
  Schwann = c("Pllp", "Prx", "Mag", "Pmp2", "Ncamp"),
  Myelinating_Schwann = c("Erbb3", "S100b", "Mbp", "Pmp22", "Mpz"),
  Nonmyelinating_Schwann = c("Ngfr", "P75", "Cdh2", "L1cam", "Ednrb", "Emp1", "Sema3e"),
  Immune = c("Csf1r", "Itgax", "Ighm", "Trbc2", "S100a9"),
  Proliferating = c("Mki67", "Top2a", "Cks1b", "H2afx", "Cks2")
)

# Generate feature plot PDF
plot_feature_pdf(
  adata,
  marker_sets,
  "path/to/CellMarkersFeaturePlot.pdf"
)


# 6. Review the excel file and plots generated to annotate the SCCAF clusters ####


# 7. Biological annotation of clusters ####
# Based on the excel file and plots, annotate the clusters 
adata$SCCAF <- recode(
  adata$SCCAF,
  '0'= 'Neurons',
  '1'= 'Sat Glia & nSchwann',
  '2'= 'Fibroblasts',
  '3'= 'Mono, Mac, & DC',
  '4'= 'Endothelial',
  '5'= 'mSchwann',
  '6'= 'Pericytes',
  '7'= 'T Cells',
  '8'= 'B Cells',
  '9'= 'Neutrophils',
  '10'= 'RBCs'
)


# 9. Save ####
# Save processed RDS object 
saveRDS(
  adata,
  file.path(output_dir, "adata_cln.RDS")
)

# Save processed H5ad object 
reticulate::use_condaenv("wscrna_py38", required = TRUE)
loompy <- reticulate::import("loompy")

as.anndata(
  x = adata,
  file_path = output_dir,
  file_name = "adata_cln.h5ad",
  assay = "SCT",
  main_layer = "data",
  other_layers = c("counts"),
  transfer_norm.data = TRUE,
  transfer_dimreduc = TRUE,
  verbose = TRUE
)
