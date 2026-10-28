library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(stringr)
library(readxl)
library(tibble)
library(BuenColors)

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2.rds")

sgRNA_qc <- read.table(
  file = "../output/PS2_tRNA_gRNA_calls.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
cells_qc <- sgRNA_qc$barcode
so_p2 <- so_p2[,cells_qc]
dim(so_p2) # 18099 103372 (prev. 103853)

sparse_sgRNA <- sparseMatrix(i = as.factor(sgRNA_qc$max_id), j = as.factor(sgRNA_qc$barcode), x = sgRNA_qc$total_umi)
rownames(sparse_sgRNA) <- levels(as.factor(sgRNA_qc$max_id))
colnames(sparse_sgRNA) <- levels(as.factor(sgRNA_qc$barcode))
cells <- intersect(colnames(sparse_sgRNA),colnames(so_p2))
sparse_sgRNA <- sparse_sgRNA[,cells]
so_p2[["sgRNA"]] <- CreateAssayObject(sparse_sgRNA[,cells])

identity_vec <- sgRNA_qc$max_id; names(identity_vec) <- sgRNA_qc$barcode
so_p2@meta.data$identity <- unname(identity_vec[colnames(so_p2)])
so_p2@meta.data$identity_NTC <- case_when(
  grepl("targeting", so_p2@meta.data$identity) ~ "NTC",
  TRUE ~ so_p2@meta.data$identity
)
head(so_p2@meta.data)

DefaultAssay(so_p2) <- "sgRNA"
so_p2 <- NormalizeData(so_p2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p2 <- FindClusters(so_p2, resolution = 0.2)

DefaultAssay(so_p2) <- "RNA"
so_p2 <- NormalizeData(so_p2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p2 <- FindClusters(so_p2, resolution = 0.2)

saveRDS(so_p2, file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
