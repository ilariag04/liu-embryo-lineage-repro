# liu-embryo-lineage-repro

Reproducibility of Liu et al. (2025) - Nature "Genome-coverage single-cell histone modifications for embryo lineage tracing". 

Part of the thesis project: **Single cell RNA-seq, A gold standard approach for a FAIR analysis**.

## Reproducibility Notes and Challenges
During the computational replication of the single-cell RNA sequencing (scRNA-seq) analysis, a significant reproducibility barrier was identified. The original Materials and Methods section stated the use of Seurat and R, but **failed to declare the specific minor releases or provide a dependency lockfile**. 

Because single-cell analysis packages change significantly across versions, running the original code "as is" on modern environments causes syntax errors (e.g., changes in how `Idents` are handled).

To guarantee true FAIR reproducibility:
1. All core mathematical algorithms and statistical parameters from the original paper (e.g., number of PCs, scale factor = 10,000, n_neighbors = 30) have been strictly preserved to ensure comparable biological results.
2. The R code has been refactored exclusively to update deprecated Seurat syntax.
3. The entire reproducible workflow has been encapsulated in the provided `Dockerfile` to create a stable, isolated environment deployable on JupyterHub.
