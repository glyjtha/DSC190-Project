# DSC 190 Statistical Genomics Project
## Molecular Mechanisms Based on Breast Cancer Subtypes in TNBC
By: Qingyu Wang, Juntong Ye

## Introduction
This project explores the genetic underpinnings of breast tissue traits by integrating GWAS summary statistics, eQTL data from GTEx, and RNA-seq data from GEO. The analysis focuses on identifying key genes and variants that influence breast tissue traits and their expression patterns.

## Objectives
- Identify SNPs acting as eQTLs in TNBC subtypes.
- Uncover signaling pathways critical to tumor biology and therapy.
- Integrate GWAS, eQTL, and RNA-seq data to advance personalized medicine.
  
## Data Sources
- **GWAS Data**: GCST90454344.tsv containing SNP IDs, effect sizes, p-values, and alleles from GWAS.
- **eQTL Data**: From the [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/qtl), featuring gene IDs, variant IDs, TSS distances, p-values, and effect sizes.
- **RNA-Seq Data**: From the [GEO Dataset GSE180878](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE180878), containing single-cell RNA sequencing data:
  - `GSE180878_Li_Brugge_10XscRNAseq_Metadata_human.csv`
  - `GSE180878_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_human.csv`

## How to Use the Code
1. **Clone the Repository**:
   ```bash
   https://github.com/glyjtha/DSC190-Project.git
