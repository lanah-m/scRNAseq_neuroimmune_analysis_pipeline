# ============================================================
# CellChat analysis across experimental conditions
#
# Purpose:
# This script performs CellChat-based cell–cell communication
# analysis across multiple experimental conditions derived
# from a Seurat object. It subsets data by condition, refines
# cell-type annotations using external metadata, computes
# ligand–receptor interactions, and generates multiple
# visualization outputs and interaction summaries.
#
# Inputs:
# - Seurat object containing integrated single-cell RNA data
# - External CSV annotation files for B cells, T cells,
#   and neutrophils
# - Defined condition groups for comparative analysis
#
# Outputs:
# - CellChat objects per condition
# - Heatmaps, chord plots, and bubble plots
# - Interaction count matrices
# - Saved RDS and CSV outputs per condition
#
# Additional notes:
# - CellChat version must be >= 2
# - "hyperalgesic priming" is not yet biologically validated and 
#   "resolved" might be more accurate term
# ============================================================

# 0. Preprocessing And Setup ####
#Load packages
library(dplyr)
library(reshape2)
library(Seurat)
library(clustree)
library(hdf5r)
library(SeuratDisk)
library(openxlsx)
library(Nebulosa)
library(CellChat)
library(patchwork)
library(future)
library(presto)
library(tidyr)

# Further set up
options(stringsAsFactors = FALSE)
setwd("path/to/cellchat_directory")
set.seed(100)

# Increase future global size to avoid memory errors
options(future.globals.maxSize = 15700 * 1024^2) 


# 1. User-defined function: save dynamic bubble
#
# Purpose:
# Generate CellChat bubble plots with dynamic sizing based on
# number of source and target cell types.
#
# Inputs:
# - CellChat object
# - sources: vector of source cell types
# - targets: vector of target cell types
# - filename: output PDF path
#
# Output:
# - PDF file containing dynamically sized bubble plot
save_dynamic_bubble <- function(cellchat, sources, targets, filename) {
  
  # Generate Bubble Plot Object
  p <- netVisual_bubble(
    cellchat,
    sources.use = sources,
    targets.use = targets,
    remove.isolate = TRUE,
    return.data = FALSE
  )
  
  # Compute Axis Complexity For Dynamic Scaling
  x_count <- length(unique(ggplot_build(p)$layout$panel_params[[1]]$x$breaks))
  y_count <- length(unique(ggplot_build(p)$layout$panel_params[[1]]$y$breaks))
  
  # Set Dynamic PDF Dimensions
  pdf(
    file = filename,
    width  = max(6, x_count * 0.25),
    height = max(6, y_count * 0.25)
  )
  
  print(p)
  dev.off()
}


# 2. Load Data And Initialize Seurat Object ####
adata <- readRDS("path/to/allcells_cln.rds")

# Set the defult assay to be RNA rather than SCT 
DefaultAssay(adata) <- "RNA"


# 3. RNA Preprocessing For CellChat Compatibility ####
# Preprocessing of my dataset (using the RNA assay; this is what cellchat uses 
# in their vignette)
# Join all the layers together so we can split them as we like in the next step 
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])

# Split by batch
adata[["RNA"]] <- split(adata[["RNA"]], f = adata$batch)

# Normalize - normalization is performed for each batch independently
adata <- NormalizeData(adata)

# Now join them again 
adata[["RNA"]] <- JoinLayers(adata[["RNA"]])


# 4. Set Identities ####
# Set 'sampleID' as the active identity
Idents(adata) <- "sampleID"

# Check that your clusters are as you expect: 
table(adata@meta.data[["merged_metadata_neuron_predicted_id"]])


# 5. Load And Refine B Cell Annotations ####
# Read B-cell metadata
b <- read.csv("path/to/bcell_assignments.csv", header = FALSE, row.names = 1)

# Create a vector of B cell lables to keep 
bcell_labels <- c(
  "Mature B And Plasma Cells",
  "Transitional And Immature B Cells",
  "Pre-B Cell"
)

# Delete all the rows that are not one of the B cell lables defined above 
b <- b[b$V2 %in% bcell_labels, , drop = FALSE]

# Write a new csv with only the b cells that have not been deleted
write.csv(b, "path/to/bcell_subclusters.csv")

# Read that csv in 
b <- read.csv("path/to/bcell_subclusters.csv", header = TRUE, row.names = 1)

# Add the labels as a column in the seurat object
adata[["b"]] <- b

# Override the labels for only the T cells that have not been deleted
adata$merged_metadata_neuron_predicted_id[
  !is.na(adata$b)
] <- adata$b[!is.na(adata$b)]


# 6. Load And Refine T Cell Annotations ####
# Read T-cell metadata
tcell <- read.csv("path/to/tcell_assignments.csv", header = FALSE, row.names = 1)

# Create a vector of T cell lables to keep 
tcell_labels <- c(
  "Cd4 Abt", "Cd8 Abt", "Ifn Responsive Abt",
  "Ilc2", "Klra7+ Gdt", "Rorc+ Gdt", "Nk & Ilc1"
)

# Delete all the rows that are not one of the B cell lables defined above 
tcell <- tcell[tcell$V2 %in% tcell_labels, , drop = FALSE]

# Write a new csv with only the t cells that have not been deleted
write.csv(tcell, "path/to/tcell_subclusters.csv")

# Read that csv in 
tcell <- read.csv("path/to/tcell_subclusters.csv", header = TRUE, row.names = 1)

# Add the labels as a column in the seurat object
adata[["tcell"]] <- tcell

# Override the labels for only the T cells that have not been deleted
adata$merged_metadata_neuron_predicted_id[
  !is.na(adata$tcell)
] <- adata$tcell[!is.na(adata$tcell)]


# 7. Load And Refine Neutrophil Annotations ####
# Read neutrophil cell metadata
neutrophil <- read.csv("path/to/neutrophil_assignments.csv",
                       header = FALSE, row.names = 1)

# Create a vector of neutrophil cell lables to keep 
neutrophil_labels <- c("G0", "G1", "G2", "G3", "G4", "G5a", "G5b", "G5c")

# Delete all the rows that are not one of the neutrophil cell lables defined above 
neutrophil <- neutrophil[neutrophil$V2 %in% neutrophil_labels, , drop = FALSE]

# Write a new csv with only the neutrophil cells that have not been deleted
write.csv(neutrophil, "path/to/neutrophil_subclusters.csv")

# Read that csv in 
neutrophil <- read.csv("path/to/neutrophil_subclusters.csv",
                       header = TRUE, row.names = 1)

# Add the labels as a column in the seurat object
adata[["neutrophil"]] <- neutrophil

# Override the labels for only the neutrophil cells that have not been deleted
adata$merged_metadata_neuron_predicted_id[
  !is.na(adata$neutrophil)
] <- neutrophil[!is.na(adata$neutrophil)]


# 8. Final Filtering And Visualization ####
# Set the active identity so the next line works as inteded
Idents(adata) <- "merged_metadata_neuron_predicted_id"

# Remove any remaining unclassified b and t cell  
adata <- subset(
  adata,
  idents = setdiff(levels(adata),
                   c("T & Nk Cells", "B Cells", "Neutro"))
)

# Check that it looks correct: 
table(adata@meta.data[["merged_metadata_neuron_predicted_id"]]) 

# Plot the clusters and make sure it looks right 
DimPlot(
  adata,
  reduction = "umap",
  label = TRUE,
  group.by = "merged_metadata_neuron_predicted_id",
  raster = FALSE
)


# 9. Rename Cluster Groups ####
# First, make sure it's a character
adata$merged_metadata_neuron_predicted_id <- as.character(
  adata$merged_metadata_neuron_predicted_id
)

# Now, rename to combine clusters as you like 
adata$merged_metadata_neuron_predicted_id[
  adata$merged_metadata_neuron_predicted_id %in% c("Pep1", "Pep2")
] <- "Pep"

adata$merged_metadata_neuron_predicted_id[
  adata$merged_metadata_neuron_predicted_id %in% c("G0", "G1", "G2")
] <- "G0-2"

adata$merged_metadata_neuron_predicted_id[
  adata$merged_metadata_neuron_predicted_id %in% c("G5a", "G5b", "G5c")
] <- "G5"

adata$merged_metadata_neuron_predicted_id <- factor(
  adata$merged_metadata_neuron_predicted_id
)


# 10. Define Experimental Conditions ####
conditions_list <- list(
  Naive = c("naive"),
  SNI = c("sni_acute", "sni_mid", "sni_chronic"),
  CCI = c("cci_acute", "cci_mid", "cci_chronic"),
  Pain_Initiation = c("sni_acute", "cci_acute"),
  Pain_Resolution = c("cci_mid"),
  Pain_Chronification = c("sni_mid", "sni_chronic"),
  Hyperalgesic_Priming = c("cci_chronic")
)

all_idents <- levels(adata$merged_metadata_neuron_predicted_id)


# 11. CellChat Analysis Loop ####
for (condition_name in names(conditions_list)) {
  
  message("Processing: ", condition_name)
  
  # Subset Data By Condition
  adata_subset <- subset(
    adata,
    subset = sampleID %in% conditions_list[[condition_name]]
  )
  
  # Convert assay to v5 if needed
  adata_subset[["RNA"]] <- as(
    object = adata_subset[["RNA"]],
    Class = "Assay5"
  )
  
  UpdateSeuratObject(adata_subset)
  
  # Set identities explicitly and fix factor levels
  Idents(adata_subset) <- adata_subset$merged_metadata_neuron_predicted_id
  Idents(adata_subset) <- factor(Idents(adata_subset), levels = all_idents)
  
  # Create CellChat Object
  cellchat <- createCellChat(
    object = adata_subset,
    group.by = "ident",
    assay = "RNA"
  )
  
  # Load CellChat Database
  CellChatDB <- CellChatDB.mouse
  cellchat@DB <- subsetDB(CellChatDB) # use all of DB except for 
  # "Non-protein Signaling" for cell-cell communication analysis
  
  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Inference
  cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat)
  
  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Define Sources And Targets
  valid_clusters <- levels(cellchat@idents)
  
  sources <- c(
    "Migratory DCs", "pDCs", "cDC1", "cDC2",
    "Prolif DCs", "Mast", "Mono", "Prolif Mono", "MDM",
    "Cd163+ Mac", "Ifn Mac", "MG-like Mac", "Ccr2+ Mac",
    "G0-2", "G3", "G4", "G5",
    "Mature B And Plasma Cells", "Transitional And Immature B Cells",
    "Pre-B Cell", "CD4 AbT", "CD8 AbT", "Ifn Responsive AbT",
    "ILC2", "Klra7+ γδT", "Rorc+ γδT", "NK & ILC1"
  )
  
  targets <- c("NF1", "NF2", "NF3", "PEP", "SST", "NP", "cLTMR1")
  
  # Define the clusters to show in your plots
  sources_to_show <- intersect(sources, valid_clusters)
  targets_to_show <- intersect(targets, valid_clusters)
  
  # Ensure all rows/columns exist in netP
  all_groups <- unique(c(sources, targets))
  missing <- setdiff(all_groups, rownames(cellchat@netP[[1]]))
  
  # Adjust plot Sizing
  x_count <- length(targets_to_show)
  y_count <- length(sources_to_show)
  
  pdf_width  <- max(6, x_count * 0.3)
  pdf_height <- max(6, y_count * 0.3)
  
  # Heatmap (Weight)
  # The sys.sleep() is nessary for the code to work so don't remove it
  heatmap_plot <- netVisual_heatmap(
    cellchat, 
    row.show = sources_to_show, 
    col.show = targets_to_show, 
    measure = "weight", 
    slot.name = "netP")
  
  pdf(
    file = paste0(condition_name, "_cellchat_signalling", "_strength_heatmaps.pdf"), 
    width = pdf_width, 
    height = pdf_height)
  
  Sys.sleep(0.5)
  print(heatmap_plot)  # Explicitly print the plot
  Sys.sleep(0.5)
  dev.off()   #Turn the PDF off 
  Sys.sleep(0.5)
  
  # Heatmap (Count)
  heatmap_plot <- netVisual_heatmap(
    cellchat, 
    row.show = sources_to_show, 
    col.show = targets_to_show, 
    measure = "count", 
    slot.name = "netP")
  
  pdf(
    file = paste0(condition_name, "_cellchat_signalling", "_count_heatmaps.pdf"), 
    width = pdf_width, 
    height = pdf_height)
  
  Sys.sleep(0.5)
  print(heatmap_plot)  # Explicitly print the plot
  Sys.sleep(0.5)
  dev.off()   #Turn the PDF off 
  Sys.sleep(0.5)
  
  # Chord Plot
  pdf(
    file = paste0(condition_name, "_cellchat_signalling", "_cord_plot.pdf"), 
    width = 9, 
    height = 6) # This line goes first to avoid an error. Don't move it
  
  Sys.sleep(0.5)
  
  chord_plot <- netVisual_chord_gene(
    cellchat, 
    sources.use = sources_to_show, 
    targets.use = targets_to_show, 
    slot.name = "netP", 
    reduce = 0.01)
  
  Sys.sleep(0.5)
  print(chord_plot)  # Explicitly print the plot
  Sys.sleep(0.5)
  dev.off()   #Turn the PDF off 
  Sys.sleep(0.5)
  
  # Bubble Plot
  save_dynamic_bubble(
    cellchat,
    sources = sources_to_show,
    targets = targets_to_show,
    filename = paste0(condition_name, "_cellchat_all_immune_bubbleplot.pdf")
  )
  
  # Subset CellChat Object
  cellchats <- subsetCellChat(
    cellchat,
    idents.use = c(sources_to_show, targets_to_show)
  )
  
  # Extract the inferred cellular com network as a data frame 
  df.nets <- subsetCommunication(cellchats)
  
  # Count number of interactions between source-target pairs
  interaction_count_matrix <- df.nets %>%
    filter(source %in% sources_to_show,
           target %in% targets_to_show) %>%
    group_by(source, target) %>%
    summarise(num_interactions = n(), .groups = "drop") %>%
    dcast(source ~ target, value.var = "num_interactions", fill = 0)
  
  # Save number of interactions to excel sheet 
  write.xlsx(interaction_count_matrix,
             file = paste0(condition_name, "_interaction_counts.xlsx"),
             rowNames = FALSE)
  
  # Save remaining outputs
  saveRDS(cellchat,
          file = paste0(condition_name, "_cellchat_allsignalling.rds"))
  
  saveRDS(cellchats,
          file = paste0(condition_name, "_cellchat_immune_to_neuron_signalling.rds"))
  
  write.csv(df.net,
            file = paste0(condition_name, "_cellchat_allsignalling.csv"),
            row.names = FALSE)
  
  write.csv(df.nets,
            file = paste0(condition_name, "_cellchat_immune_to_neuron_signalling.csv"),
            row.names = FALSE)
  
  message("Finished processing: ", condition_name)
}