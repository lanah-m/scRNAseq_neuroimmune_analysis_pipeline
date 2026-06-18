# scRNAseq-neuroimmune-pipeline

## Overview

This repository contains the complete computational workflow developed for my MSc thesis investigating **Neuroimmune Dynamics in the Dorsal Root Ganglia During Peripheral Nerve Injury–Induced Neuropathic Pain** using single-cell RNA sequencing (scRNA-seq). Raw data and thesis results are not included. Instead, the repository documents the code from processing raw sequencing output to downstream biological analyses.

---

### Environment & reproducibility

Single-cell analysis is highly sensitive to software versions. Differences in package versions can lead to changes in clustering, dimensionality reduction, or statistical outputs. To ensure reproducibility, both Python and R environments are included in the renv.lock (for r) and requirements.txt (for python) files. 

---

### 00_preprocessing_alignment

scRNA-seq data are initially generated as raw sequencing reads composed of nucleotide sequences (A, T, C, G). These reads must be aligned to a reference genome in order to quantify gene expression. In this project, a **custom reference genome** was constructed to account for a knock-in reporter gene (tdTomato), ensuring that transgene expression was correctly captured rather than discarded or misaligned. cellranger_custom_genome_creation.md describes how to create this custom reference genome and cellranger_alignment.md describes how to use a reference genome to align the sequencing data. 

---

### 01_ambient_rna_and_empty_droplet_removal

Droplet-based scRNA-seq platforms often capture ambient RNA originating from lysed cells or the surrounding solution. This contamination can lead to false-positive expression. CellBender is used to remove ambient RNA while distinguishing true cells from empty droplets. cellbender_code.bat is the code to input into terminal or command prompt and the README.md describes in mode detail what was done.

---

### 02_load_seurat_objects

After the initial quality control with cellbender, filtered count matrices are loaded into R as Seurat objects and then all samples are merged together as a single object. load_merge_seurat_objects.R is the code to do that. 

---

### 03_qc_filtering

Quality control is performed to remove technical artifacts such as doublets and dying cells. qc_filtering.R is the code to do that. 

---

### 04_clustering_and_integration_pipeline

Following data cleaning with quality control, clustering using Seurat is performed. Initial clustering is conducted (04_1_clustering_pre_integration.R) without batch correction to assess whether integration is necessary. When batch effects are detected, CCA-based integration is applied (04_2_integration.R). Post-integration clustering (04_3_clustering_post_integration.R) allows evaluation of whether batch correction improves or distorts biologically meaningful structure. 

---

### 05_sccaf_ml_optimization

Following the clustering done using Seurat, SCCAF (Single-Cell Clustering Assessment Framework) is used. SCCAF uses machine learning to improve clustering. 05_1_sccaf_clustering_optimization.ipynb is the implementation of SCCAF and 05_2_reading_sccaf_results_into_r.R takes the SCCAF clustering results and reads it into R. 

---

### 06_de_novo_clustering_annotation

After establishing the clusters of cells, you now how to annotate them. You can do this using two parallel methods 1) use the genes that are unregulated in the cluster relative to the other clusters and 2) using known marker genes of certain populations (like Iba1 for macrophage for example). The annotation.R script generates an excel sheet with the unregulated genes in each cluster as well as FeaturePlots showing the expression of markers genes. 

---

### 07_reference_mapping_annotation

In some cases, steps 04-06 are circumvented by taking the quality controlled scRNA sequencing data and mapping it to another already annotated scRNA seq dataset. Here, mapping refers to the projection of query cells onto a shared low-dimensional embedding using a reference scRNA seq/snRNA seq atlas such that the annotations of the atlas
cells are transferred to the query cells based on shared gene expression patterns.

---

### 08_differential_gene_expression

One important comparison that is often of interest is that between sample types. For example, what makes the acute sample different than the chronic sample? To get at this, a differential gene expression analyses is performed to identify transcriptional changes associated with a specific samples (08_1_differential_gene_expression_analysis.R). Results can then be visualized using cluster-by-condition heatmaps (08_2_visualizing_differential_gene_expression.ipynb).

---

### 09_functional_enrichment

While the genes obtained from the differential gene expression can help us understand what's changing from sample to sample on a gene level, it can be difficult to understand whats happening on a larger scale. While differential gene expression tells us which genes change between samples, it doesn’t immediately show what those changes mean biologically. Functional enrichment analysis helps by linking genes to known biological terms. Each gene in the genome is already annotated in databases like Gene Ontology (GO) with terms describing what it does, where it works in the cell, or what processes it is involved in (for example: “inflammatory response” or “axon regeneration”). Enrichment analysis takes your list of unregulated genes and checks which of these GO terms appear more often than expected by chance. In other words, it looks for biological labels that are shared across many of your changed genes, and highlights them. This turns a long gene list into clearer functional themes. To do this, you first go and input your gene list into https://biit.cs.ut.ee/gprofiler/gost with the correct specific and your databases of interest check off. You then download the CSV with all the genes and stats associated with it. To visualize these results, that csv is the input for both the bubble_plot_visualization_of_go_terms.R (to make Bubble plots) and go_term_upset_plots.R (to make UpSet plots) scripts.

---

### 10_intercellular_communication

CellChat (cellchat.R) is used to infer potential ligand–receptor interactions between cell populations. 

---

### 11_gene_module_scoring

gene_module_scoring.R calculates gene module scores by evaluating the average expression of predefined gene sets (modules) within each cell. These modules represent biological pathways or functional programs (e.g., stress response, inflammation, senescence), allowing comparison of pathway activity across conditions or cell types at the single-cell level. The average expression of the gene module is plotted on a UMAP similar to FeaturePlots generated in 06_de_novo_clustering_annotation. 

---

### 12_subclustering

Following an initial clustering of scRNA seq data, you make want to subcluster a cluster (or a few clusters) further. For example, when you have neurons, immune cells, pericytes, endothelial cells, etc. but you are specifically interested only in neurons, you make choose to take all the neuron clusters and then subcluster those. Generally you do this step either after 06_de_novo_clustering_annotation or 07_reference_mapping_annotation. You can do this in two ways. The first way is by using the FindSubCluster() function in Seurat and the code for that is in subclustering_using_FindSubCluster.R. The second way is to use the subset() function in Seurat and then redo 04_clustering_and_integration_pipeline to 06_de_novo_clustering_annotation following all steps following (and including) SCTransform() in the 04_clustering_and_integration_pipeline. Alternatively, you can use the subset() function in Seurat and then do 06_de_novo_clustering_annotation.

---
