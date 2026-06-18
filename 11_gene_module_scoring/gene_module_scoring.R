# ============================================================
# Gene set feature plotting in macrophages
#
# Purpose:
# This script generates gene set enrichment feature plots
# from a Seurat object of macrophages using AddModuleScore().
# It visualizes predefined functional gene programs on UMAP.
#
# Inputs:
# - Seurat object with UMAP reduction
# - Predefined macrophage gene sets
#
# Outputs:
# - Combined UMAP feature plot PDF
# - Individual gene set feature plots
#
# Notes:
# - Gene sets are adapted from Sanin et al. (2022)
# - Used for functional annotation of macrophage states
# ============================================================


# 0. Setup ####
suppressMessages({
  library(dplyr)
  library(Seurat)
  library(clustree)
  library(hdf5r)
  library(SeuratDisk)
  library(openxlsx)
  library(Nebulosa)
  library(sctransform)
  library(SeuratObject)
  library(future)
  library(glmGamPoi)
  library(zellkonverter)
  library(Matrix)
  library(rhdf5)
})

# Set working directory
setwd('path/to/gene_set_feature_plots')

# Load macrophage Seurat object
macrophage <- readRDS('path/to/mmd_cln.rds')

# Replace UMAP reduction for consistency
macrophage@reductions[["umap"]] <- macrophage@reductions[["umap.cca"]]


# 1. Define gene sets ####
# These gene sets represent macrophage functional programs
ident.scores <- list(
  "Complement & Phagocytosis" = c('C1qc','F13a1','C1qa','C4b','Cfh','C5ar1',
                                  'Snx2','Tgfbr2','Dab2','Folr2','Cltc','Wwp1',
                                  'Cd209d','Mrc1','Cd209f','Cd209g','Cd36',
                                  'Ctsb','Lgmn','Cltc','Cd63'),
  
  "ECM & Actin Regulation" = c('Cd44','Sdc1','Fn1',
                               'Pfn1','Actg1','Tmsb4x'),
  
  "Antigen Presentation" = c('H2-Ab1','H2-Aa','H2-Eb1',
                             'H2-Oa','H2-DMb2','H2-Ob','H2-DMb1'),
  
  "Innate Immune Response" = c('Tnfaip8l2','Cyba','Rsad2','Anxa1','Ifitm3',
                               'Fcgr1','Fgr','Oasl2','Clec4e','Clec4d',
                               'Pglyrp1','Oas3','Isg20','Samhd1',
                               'Hmgb2','Rnase6','Slpi','Msrb1','Gbp2'),
  
  "Phagosome" = c('Ctss','Cyba','Msr1',
                  'Fcgr1','Coro1a','Thbs1','Ncf4','Fcgr3'),
  
  "ROS Producers" = c('Nox4','Dao','Qsox1','lox','nos1','Nos2','nos3','Ncf2',
                      'Mtnd1','Sdha','Sdhb','Sdhc','Sdhd','Cyc1','Ogdh','Gpd1',
                      'Ddo','Pipox','Hao1','Hao2','Paox','Uox','Duox1','Duoxa1',
                      'Cox1','Cox-2','Xdh','Nrf2','Syk','Acox1','Acox2','Acox3',
                      'Tafa4','Romo1','Thbs1','Ulbp1','mt-Co2','ERo1α and β',
                      'PDI','Nox1','Cybb','Cyba','Nox3','Atpsckmt','Etfa',
                      'Sod1','Sod2'),
  
  "ROS Scavengers" = c('Prdx1','Prdx2','Prdx3','Prdx4','Prdx5','Prdx6',
                       'Txn1','Txn2','Txnrd1','Txnrd2','Txnrd3','Txndc5',
                       'Gpx1','Gpx2','Gpx3','Gpx4','Gpx5','Gpx6','Cat',
                       'Glrx','Glrx2','Glrx3','Glrx5','Nqo2','Scara3',
                       'Gsr','Sod1','Sod2'),
  
  "Senescence & SASP" = c('Cdkn2a','Cdkn1a','Cdkn2d','Casp8','Il1b','Glb1',
                          'Serpine1','Axl','Ccl2','Ccl3','Ccl4','Ccl5',
                          'Csf1','Itgax','Mmp9','Mmp12','Spp1','Tgfb1',
                          'Vim','Gdf15','Tnfrsf1a'),
  
  "ECM Organization" = c('Col1a1','Nid1','Dpt','B4galt1',
                         'Lum','Col3a1','Ccdc80','Serpinh1','Ddr2'),
  
  "Cycling" = c('Racgap1','Cks1b','Stmn1','Ran','Smc4','Top2a',
                'Cks2','Ube2c','Cenpw','Smc2','Cenpa'),
  
  "Microglial Signature" = c('Sall1','Tmem119','Htra1','P2ry12','Olfml3')
)


# 2. Module scoring ####
# Score gene set activity per cell
macrophage <- AddModuleScore(
  object = macrophage,
  features = ident.scores,
  name = "Enrichment"
)


# 3. Feature plotting ####
# Generate UMAP feature plots for enrichment scores
f1c <- FeaturePlot(
  object = macrophage,
  features = paste0('Enrichment', 1:length(ident.scores)),
  reduction = "umap",
  min.cutoff = 'q15',
  max.cutoff = 'q85',
  pt.size = 0.01,
  combine = FALSE
)

# Clean plot appearance for publication
f1c <- lapply(f1c, function(x) {
  x + theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(face = 'plain', size = 20)
    #+ NoLegend() # remove the # and run again to generate a plot w/o a legend
  )
})

# Add gene set titles
for (i in seq_along(f1c)) {
  f1c[[i]] <- f1c[[i]] +
    labs(title = names(ident.scores)[[i]])
}


# 4. Export plots ####
# Combine all feature plots into one figure
combined_plot <- wrap_plots(f1c, ncol = length(f1c))

# Save combined plot
ggsave(
  "path/to/gene_set_plots.pdf",
  combined_plot,
  width = 17.5,
  height = 15
)

# Save each feature plot separately
for (i in seq_along(f1c)) {
  file_name <- paste0(
    "path/to/gene_set_feature_plots/",
    names(ident.scores)[[i]], "_feature_plot.pdf"
  )
  
  ggsave(file_name, plot = f1c[[i]], width = 5.5, height = 5)
}