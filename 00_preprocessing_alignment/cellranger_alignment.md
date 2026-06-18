# 10x Genomics Cloud Analysis: Aligning Reads to a Custom Reference Genome

This document describes how to use the 10x Genomics Cloud Analysis platform to align sequencing reads to a custom reference genome generated using `cellranger mkref`.

---

## Overview

The 10x Genomics Cloud Analysis platform enables processing of single-cell sequencing data without requiring local high-performance computing resources. In this workflow, FASTQ files are uploaded and aligned to a custom reference genome that includes an additional protein/transgene sequence (e.g., tdTomato) integrated into a standard reference.

---

## Prerequisites

Before starting, ensure you have:

* FASTQ files generated from sequencing
* A custom reference genome generated using `cellranger mkref`
* A 10x Genomics Cloud Analysis account
* Access to a compatible web browser

---

## Step 1: Access 10x Genomics Cloud Analysis

Navigate to the [10x Genomics Cloud Analysis platform](https://www.10xgenomics.com/products/cloud-analysis) and log in using your 10x Genomics account credentials.

---

## Step 2: Upload FASTQ Files

1. Create a new analysis project
2. Upload FASTQ files generated from your sequencing run. Using the 10x Genomics Cloud command line interface (CLI) is the fastest way to upload data but it will require time to set up on your computer.
3. Confirm that sample metadata is correctly assigned

---

## Step 3: Upload Custom Reference Genome

1. Select **“Custom Reference”** during pipeline setup
2. Upload the entire `mkref` output directory 

---

## Step 4: Configure the Analysis Pipeline

Select the Fastq files you uploaded into the cloud and select the appropriate libray type. Press create new analysis and select the apprioptate pipeline depending on your experiment.

---

## Step 5: Run Alignment

Start the analysis. The cloud platform will align reads to the custom reference genome. Processing time depends on dataset size but will likley take a couple of hours per library.

---

## Step 6: Download Results

Once complete, download all files including:

* Gene expression matrix (filtered and raw)
* Alignment BAM files (optional)
* QC metrics and summary reports

---

## Reference

Detailed steps can also be found in the documention present in the [10x Genomics Cloud Analysis Platform](https://www.10xgenomics.com/products/cloud-analysis)
