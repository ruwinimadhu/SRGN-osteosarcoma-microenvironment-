# SRGN-osteosarcoma-microenvironment-
SRGN in osteosarcoma — analysis code

Analysis code for the manuscript:

"Serglycin is a conserved osteoclast-enriched marker of the osteosarcoma microenvironment with prognostic relevance across dogs and humans"

Ruwini Madhushika Herath. Preprint on bioRxiv; submitted to npj Precision Oncology.

This repository contains the R scripts used to produce the single-cell, correlation, survival, and transcription-factor analyses reported in the paper, across canine and human osteosarcoma datasets.

Repository contents
Script	What it does	Datasets
SRGN_complete_analysis.R	Main pipeline in five layers: (1) canine scRNA-seq annotation and SRGN expression by cell type, (2) bulk RNA-seq SRGN correlation, (3) canine survival (KM), (4) SRGN transcription-factor regulon activity (DoRothEA/VIPER), (5) human scRNA-seq annotation and SRGN expression, plus the cross-species summary.	GSE252470, GSE238110, COTC022 clinical, GSE162454
srgn_survival_GSE39055.R	Human recurrence-free survival (Fig. 3): median-split KM, optimal-cutpoint KM, and Cox proportional hazards. Self-contained — downloads the data via getGEO().	GSE39055
srgn_target_os_validation.R	Independent prognostic validation: SRGN as a continuous, pre-specified predictor in Cox regression (overall survival). Self-contained — pulls data programmatically.	TARGET-OS

The two survival/validation scripts are standalone and download their own data, so either can be run on its own without the main pipeline.

Data availability

All datasets are public.

Accession	Description	Access
GSE252470	Canine osteosarcoma scRNA-seq (10x)	NCBI GEO
GSE238110	Canine (DOG2) bulk RNA-seq	NCBI GEO
COTC022 (SOC Data Set)	Canine clinical / survival data	COTC022 supplementary data
GSE162454	Human osteosarcoma scRNA-seq (10x)	NCBI GEO
GSE39055	Human osteosarcoma, recurrence-free survival	NCBI GEO (auto-downloaded)
TARGET-OS	Human osteosarcoma expression + clinical	UCSC Xena / GDC (auto-downloaded)

srgn_survival_GSE39055.R and srgn_target_os_validation.R fetch their data automatically. For SRGN_complete_analysis.R, download the GEO files above and place them in your working directory before running (see below).

Requirements
R (developed on 4.5.1)
Bioconductor

R packages used across the scripts:

r
# CRAN
install.packages(c(
  "Seurat", "ggplot2", "ggpubr", "dplyr", "tidyr", "readxl", "readr",
  "ggrepel", "patchwork", "Matrix", "survival", "survminer", "UCSCXenaTools"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "dorothea", "viper", "TCGAbiolinks"))
Running the scripts

Standalone (nothing to download manually):

r
source("srgn_survival_GSE39055.R")     # human RFS — Fig. 3
source("srgn_target_os_validation.R")  # TARGET-OS validation

Main pipeline:

Download the GEO/COTC input files listed under Data availability into one folder.
Open SRGN_complete_analysis.R and set the working directory near the top (setwd(...)) to that folder.
Run the script. Figures and result files (.png, .rds, .csv) are written to the working directory.
Code availability statement

For the manuscript, you can use:

All analysis code is available at https://github.com/<your-username>/<repo>. Scripts reproduce the single-cell, correlation, survival, and transcription-factor analyses across the canine and human osteosarcoma datasets listed above. All input data are publicly available under the accessions provided.

License

MIT (see individual script headers).

Contact

Ruwini Madhushika Herath
