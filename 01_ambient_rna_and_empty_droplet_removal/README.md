# CellBender Background Removal (scRNA-seq)

This pipeline applies CellBender remove-background to 10x Genomics scRNA-seq datasets to remove ambient RNA contamination prior to downstream analysis.

## Tool
- CellBender v0.3.2

## Input
- raw_feature_bc_matrix.h5 files from Cell Ranger

## Samples
- Naive_B1
- CCI_Chronic_B3
- CCI_Mid_B3
- CCI_Chronic_B6
- SNI_Acute_B6

## Method

This step was originally executed using a Windows batch (.bat) environment. The command shown below is the equivalent execution command reformatted in Unix-style (Bash) syntax for reproducibility and compatibility with Linux-based bioinformatics environments commonly used in CellBender workflows.

```bash
cellbender remove-background \
  --input "/path/to/sample/raw_feature_bc_matrix.h5" \
  --output "/path/to/sample/sample"