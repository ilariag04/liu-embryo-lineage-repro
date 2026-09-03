setwd("/home/jovyan/TACIT/fig2")
library(Seurat)
library(dplyr)

message("1. Caricamento file...")
raw_df <- read.csv("/home/jovyan/TACIT/scRNA-seq_2014/GSE45719_scRNA.csv", stringsAsFactors = FALSE)

# Sanificazione e isolamento dei nomi
genes <- make.unique(trimws(as.character(raw_df[, 1])))
cells <- make.unique(trimws(as.character(colnames(raw_df)[-1])))

message("2. Conversione numerica e TPM...")
mat <- suppressWarnings(data.matrix(raw_df[, -1]))
mat[is.na(mat)] <- 0

col_sums <- colSums(mat)
col_sums[col_sums == 0] <- 1
mat_tpm <- sweep(mat, 2, col_sums, FUN = "/") * 1e6

rownames(mat_tpm) <- genes
colnames(mat_tpm) <- cells

message("3. Creazione oggetto Seurat...")
rna_assay <- CreateAssayObject(counts = mat_tpm)
ref_obj <- CreateSeuratObject(counts = rna_assay, project = "GSE45719")

message("4. Generazione grafici QC e Pipeline...")
pdf("all_reference.pdf")
# Usiamo solo nFeature e nCount per evitare problemi di percent.mt vuoti
print(VlnPlot(ref_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 1))
print(FeatureScatter(ref_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

ref_obj <- NormalizeData(ref_obj)
ref_obj <- FindVariableFeatures(ref_obj, nfeatures = 3000)
# Regrediamo solo su nFeature_RNA per massima stabilità
ref_obj <- ScaleData(ref_obj, vars.to.regress = c("nFeature_RNA"))
ref_obj <- RunPCA(ref_obj, npcs = 20)
ref_obj <- RunUMAP(ref_obj, dims = 1:20, reduction = "pca", n.neighbors = 5L)
ref_obj <- FindNeighbors(ref_obj, dims = 1:20, reduction = "pca")
ref_obj <- FindClusters(ref_obj, resolution = 0.5)

pdf("ref_stage.pdf")
print(DimPlot(ref_obj, reduction = "umap", label = TRUE, pt.size = 2))
print(DimPlot(ref_obj, reduction = "pca", label = TRUE, pt.size = 2))
dev.off()

saveRDS(ref_obj, file = "./reference.rds")
message("ANALISI COMPLETATA CON SUCCESSO! File reference.rds generato.")
