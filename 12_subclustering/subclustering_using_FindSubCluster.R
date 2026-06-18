# ============================================================
# B cell subclustering and conserved marker analysis

# Purpose: This script performs downstream subclustering of B cells from an
# integrated Seurat object, identifies conserved markers across batches, and 
# annotates biological B cell states. 
# ============================================================


# 0. Preprocessing and setup ####
# Load packages
suppressMessages({
  library(Seurat)
  library(openxlsx)
  library(dplyr)
  library(sceasy)
  library(reticulate) 
  reticulate::use_condaenv("wscrna_py38", required = TRUE) 
  loompy <- reticulate::import('loompy') 
  library(scCustomize)
})

# Define output directory
output_dir <- "path/to/folder"

# Set seed for reproducibility
set.seed(1234)

# Load data 
immune <- readRDS("path/to/file.rds")


# 1. Subclustering B cells ####
bcell <- FindSubCluster(
  immune,
  cluster = "B Cells",
  graph.name = "SCT_snn",
  subcluster.name = "B_cell_subclustered",
  resolution = 0.3
)

# Visualize subclusters
pdf_file <- file.path(output_dir, "bcell_subclustering.pdf")
pdf(pdf_file, height = 7, width = 7)

DimPlot(
  bcell,
  group.by = "B_cell_subclustered",
  label = TRUE
)

dev.off()

# Set active identity and assay 
Idents(bcell) <- "B_cell_subclustered"
DefaultAssay(bcell) <- "RNA"


# 2. Batch-aware normalization ####
# This is nessary for the conserved marker detection that is performed later

# Join layers before splitting by batch  
bcell[["RNA"]] <- JoinLayers(bcell[["RNA"]])

# Split by batch 
bcell[["RNA"]] <- split(bcell[["RNA"]], f = bcell$batch)

# Normalize 
bcell <- NormalizeData(bcell)

# Join layers again 
bcell[["RNA"]] <- JoinLayers(bcell[["RNA"]])


# 3. Conserved marker detection across batches ####
# This is to to support in the annoation of the identified subclusters

# Define identities (subclusters) to iterate through in a loop 
idents <- c(
  "B Cells_0",
  "B Cells_1",
  "B Cells_2"
)

# Make an empty list to store results in 
bcell.markers.list <- list()

# DEG analysis for each identity (subcluster) defined above 
for (ident in idents) {
  message("Processing: ", ident)
  
  markers <- FindConservedMarkers(
    bcell,
    ident.1 = ident,
    grouping.var = "batch",
    assay = "RNA",
    slot = "data",
    only.pos = TRUE,
    min.pct = 0.5,
    logfc.threshold = 0.5
  )
  
  # Add identity column and gene names
  markers$celltype <- ident
  markers$gene <- rownames(markers)
  
  # Compute min/max adjusted p-values across batches
  adj_cols <- grep("p_val_adj$", colnames(markers), value = TRUE)
  markers$min_adj_pval <- apply(markers[, adj_cols], 1, min, na.rm = TRUE)
  markers$max_adj_pval <- apply(markers[, adj_cols], 1, max, na.rm = TRUE)
  
  # Compute min/max log2FC across batches
  logfc_cols <- grep("avg_log2FC$", colnames(markers), value = TRUE)
  markers$min_log2FC <- apply(markers[, logfc_cols], 1, min, na.rm = TRUE)
  markers$max_log2FC <- apply(markers[, logfc_cols], 1, max, na.rm = TRUE)
  
  # Compute min/max pct.1 across batches
  pct1_cols <- grep("pct\\.1$", colnames(markers), value = TRUE) # note double escape for dot
  markers$min_pct1 <- apply(markers[, pct1_cols], 1, min, na.rm = TRUE)
  markers$max_pct1 <- apply(markers[, pct1_cols], 1, max, na.rm = TRUE)
  
  # Compute min/max pct.2 across batches
  pct2_cols <- grep("pct\\.2$", colnames(markers), value = TRUE)
  markers$min_pct2 <- apply(markers[, pct2_cols], 1, min, na.rm = TRUE)
  markers$max_pct2 <- apply(markers[, pct2_cols], 1, max, na.rm = TRUE)
  
  # Store in list 
  bcell.markers.list[[ident]] <- markers
}

# 4. Combine and clean DEG results 
# Keep only desired columns for each data frame inside the list
bcell.markers.list <- lapply(bcell.markers.list, function(df) {
  df[, c(
    "celltype", "gene",
    "min_log2FC", "max_log2FC",
    "min_adj_pval", "max_adj_pval",
    "min_pct2", "max_pct2",
    "min_pct1", "max_pct1"
  )]
})

# Combine the data from the selects columns chosen into one data frame
adata.markers <- do.call(rbind, bcell.markers.list)

# Filter to include only statistically significant conserved markers 
adata_sig <- adata.markers[adata.markers$max_adj_pval < 0.05, ]

# Remove mitochondrial and ribosomal genes
adata_sig <- adata_sig[!grepl("^mt-", adata_sig$gene, ignore.case = TRUE), ]
adata_sig <- adata_sig[!grepl("^Rp[sl]", adata_sig$gene, ignore.case = TRUE), ]


# 5. Export DEG results ####
write.xlsx(
  adata_sig,
  file = file.path(output_dir, "bcell_conserved_deg_log2fc0.5_minpct0.5.xlsx")
)


# 6. Feature plots showing expression of known marker genes to support annoation ####
pdf_file <- file.path(output_dir, "bcellsubclusteringfeatureplots.pdf")
pdf(pdf_file)

FeaturePlot(
  immune,
  features = c(
    # B cell markers
    "Cd79a", "Cd79b", "Cd19",
    
    # B cell developmental stages
    "Rag1", "Rag2", "Dntt", "Tnfrsf13c", "Pax5", "Tcf3",
    "Ebf1", "Pou2f2", "Igkc", "Iglc1", "Iglc2",
    "Vpreb3", "Igll1", "Ighm", "Ighd",
    
    # Plasmablasts / plasma cells
    "Sdc1", "Xbp1", "Jchain",
    
    # Isotypes
    "Ighm", "Ighd", "Ighg1", "Ighg2b", "Ighg3", "Igha", "Ighe",
    
    # Functional / activation
    "Cd40", "Cd80", "Cd86"
  ),
  combine = FALSE
)

dev.off()


# 7. Review the excel file and plots generated to annotate the B cell clusters ####


# 8. Biological annotation of clusters ####
# Based on the excel file and plots, annotate the clusters 
bcell$B_cell_subclustered <- recode(
  bcell$B_cell_subclustered,
  "B Cells_0" = "Transitional and Immature B Cells",
  "B Cells_1" = "Mature B and Plasma Cells",
  "B Cells_2" = "pre-B Cell"
)


# 9. Save ####
# Save cluster assignments 
write.csv(
  bcell$B_cell_subclustered,
  file.path(output_dir, "b_cell_cluster_assignments.csv")
)

# Save the subsetted B Cell Object 
# Subset to only the B cell clusters
bcell <- subset(
  bcell,
  idents = c(
    "Transitional and Immature B Cells",
    "Mature B and Plasma Cells",
    "pre-B Cell"
  )
)

# Plot to make sure that you subsetted properly 
DimPlot(
  bcell,
  group.by = "B_cell_subclustered",
  label = TRUE
)

# Save processed RDS object 
saveRDS(
  bcell,
  file.path(output_dir, "bcell_clustered.RDS")
)

# Save processed H5ad object 
reticulate::use_condaenv("wscrna_py38", required = TRUE)
loompy <- reticulate::import("loompy")

as.anndata(
  x = bcell,
  file_path = output_dir,
  file_name = "bcells.h5ad",
  assay = "SCT",
  main_layer = "data",
  other_layers = c("counts"),
  transfer_norm.data = TRUE,
  transfer_dimreduc = TRUE,
  verbose = TRUE
)