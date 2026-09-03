# liu-embryo-lineage-repro

Reproducibility and re-analysis of single-cell RNA sequencing data from:
> **Liu et al. (2025)** – *"Genome-coverage single-cell histone modifications for embryo lineage tracing"*, *Nature*.  
> DOI: [10.1038/s41586-025-08656-1](https://doi.org/10.1038/s41586-025-08656-1)

---

## 📌 Project Overview
This repository is part of a bachelor's thesis project focused on computational reproducibility and standardizing single-cell RNA-seq workflows according to **FAIR** (Findable, Accessible, Interoperable, Reusable) data principles.

The pipeline reproduces the single-cell transcriptomic processing of the early embryonic reference dataset (`GSE45719`), evaluating data ingestion, legacy matrix sanitization, quality control standards, and UMAP validation to trace early embryonic lineage progression.

---

## 📂 Repository Structure

```text
├── .github/workflows/
│   └── docker-build.yml        # Automated CI/CD workflow to build and push container to GHCR
├── figures/                    # Output plots (QC violin plots, UMAPs, PCA projections)
├── scRNAseq/                   # Core R scripts for matrix sanitization and Seurat standardization
├── Dockerfile                  # Container specification (R environment, Seurat, dependencies)
├── .gitignore                  # Excludes heavy binary objects (*.rds, raw data tables)
└── README.md                   # Project documentation
```
---

## Computational Environment & Container
To ensure exact reproducibility across platforms (such as the DOra high-performance cluster and local environments), all dependencies are containerized.

Base Environment: R (containerized runtime ensuring legacy format compatibility)

Core Libraries: Seurat, SeuratObject, Matrix, tidyverse (dplyr, ggplot2)

Container Registry: GitHub Container Registry (ghcr.io)

## Running via Docker / JupyterHub
The container image is built automatically and published publicly:

```text
docker pull ghcr.io/ilariag04/liu-embryo-lineage-repro:latest
```
---

## Reproducibility Notes & Pipeline Workflow
Data Ingestion & Sanitization: Ingestion of legacy single-cell expression tables dating back to 2014 (GSE45719_scRNA.csv), applying custom sanitization wrappers to handle missing values, trailing characters, and formatting incompatibilities.

Standardization & Normalization: Execution of pure mathematical calculations (such as TPM conversion) and construction of custom S4 assay objects (CreateAssayObject) to bypass Seurat v5 architectural constraints without altering biological counts.

Quality Control & Filtering: Strict quality filtering based on library complexity and cellular feature distributions to eliminate technical artifacts.

Dimensionality Reduction & Validation: Principal component analysis (PCA) and non-linear dimensional reduction (UMAP) performed to capture transcriptomic transitions and isolate cellular clusters across early developmental stages.
