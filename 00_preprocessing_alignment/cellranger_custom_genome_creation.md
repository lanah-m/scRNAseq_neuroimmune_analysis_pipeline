# Creation of a Custom Reference Genome for Read Alignment

This document details the procedure used to generate a custom reference genome by augmenting an existing reference genome with an additional protein/transgene sequence and its corresponding gene annotation. In this implementation, the approach was applied to incorporate the tdTomato transgene, a red fluorescent protein commonly used as a reporter for lineage tracing studies. However, the same workflow can be used to add any protein or transgene of interest and for generalization purposes, placeholders were used in the workflow descibed below to facilitate adaptation to other proteins or transgenes of interest. The resulting reference genome was used for downstream read alignment and analysis.

## Background and Rationale

To enable the alignment of sequencing reads to a protein/transgene not present in the standard reference genome, a custom reference genome was constructed. This was achieved by appending a custom FASTA sequence and a corresponding GTF annotation to an existing reference genome prior to reference indexing with Cell Ranger.

## Input Files and Requirements

### Base Reference Genome

The following files from an existing reference genome were used:

- genome.fa — reference genome FASTA file  
- genes.gtf — reference gene annotation file  

Source of reference genome files: 10x Genomics mouse reference “2024-A”, based on the GRCm39 genome assembly with Gencode vM33 gene annotations, obtained from the [10x Genomics reference downloads](https://www.10xgenomics.com/support/software/cell-ranger/downloads).

### Custom Protein / Transgene

The sequence of the protein/transgene of interest was provided in the following formats:

- <protein_of_interest>.fa – FASTA file containing the protien of intrest seqeunce. In the case of tdTomato, the coding sequence was obtained from the [SnapGene plasmid database](https://www.snapgene.com/plasmids/fluorescent_protein_genes_and_plasmids/tdTomato) and the transgene was assigned the custom contig name “TDTOMATO”. 
- <protein_of_interest>.gtf – custom annotation file defining the tdTomato transgene locus, generated for incorporation into the custom reference genome. Naming conventions used: Gene-level annotation used gene_id "tdTomato_g", transcript_id "tdTomato_t", and gene_name "tdTomato_n".

## Directory Structure

A new directory was created to contain all files associated with the custom reference genome:

custom_ref_<genome_name>/
├── genome.fa
├── genes.gtf
├── <protein_of_interest>.fa
└── <protein_of_interest>.gtf

All subsequent steps were performed from within this directory.

## Construction of the Custom Reference Genome

The FASTA sequence corresponding to the protein or transgene of interest was appended to the existing reference genome FASTA file using the following command:

cat <protein_of_interest>.fa >> genome.fa

The corresponding gene annotation was appended to the existing GTF file using the following command:

cat <protein_of_interest>.gtf >> genes.gtf

### Alternative Manual Method

If the above commands were unavailable or unsuccessful, the same result was achieved manually by copying the full contents of protein_of_interest.fa and pasting them at the end of genome.fa, followed by copying the full contents of protein_of_interest.gtf and pasting them at the end of genes.gtf. Care was taken to ensure that file formatting remained intact and that no existing entries were modified.

## Reference Genome Indexing with Cell Ranger

Following construction of the combined FASTA and GTF files, the custom reference genome was indexed using cellranger mkref. This step was performed on a Linux workstation located in the University of Alberta ACE Core facility (Hobman computer).

The following command was executed from the working directory containing the modified reference files:

cellranger mkref \
  --genome=<custom_reference_name> \
  --fasta=genome.fa \
  --genes=genes.gtf

## Output and Downstream Use

The cellranger mkref command generated a reference genome directory containing indexed FASTA files and processed gene annotations compatible with downstream Cell Ranger workflows.

This custom reference genome was subsequently used for the alignment of reads from a single cell RNA sequencing expeirment utilizing trangenic mice expressing tdTomato. 

## Notes and Considerations

Chromosome and gene naming conventions in the custom FASTA and GTF files were matched to those of the base reference genome where applicable. The GTF file formatting was verified to ensure compatibility with Cell Ranger. This procedure appends new entries to the reference genome and does not alter existing annotations.

## Software and Versions

Operating system / compute environment: Linux-based high-performance workstation (Hobman computer, ACE Core facility).