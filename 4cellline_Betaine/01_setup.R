library(dplyr)
library(Seurat)
library(stringr)
library(readxl)
library(tibble)
library(BuenColors)

filter <- dplyr::filter
select <- dplyr::select

source("../../global_functions/tRNA_import_massive.R")
source("../../global_functions/tRNA_valid_massive.R")

so_4clb <- createSeurat(c("4CLB"), c("1"))
trna_raw_data <- import_tRNA(so_4clb, c("4CLB_massive"), c("1"))
so_4clb <- append_tRNA(so_4clb, trna_raw_data)
dim(so_4clb) # 18129 22682

# filter tRNA probe pair data
df_valid <- tRNA_valid(so_4clb) # Missing anticodons: UAA and UAG. Expected: UAA and UAG
so_4clb[["tRNA_pp_valid"]] <- CreateAssayObject(df_valid[,colnames(so_4clb)])

# normalize data
DefaultAssay(so_4clb) <- "tRNA_pp_valid"
so_4clb <- NormalizeData(so_4clb) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_4clb <- FindClusters(so_4clb, resolution = 0.2)

DefaultAssay(so_4clb) <- "tRNA_ac"
so_4clb <- NormalizeData(so_4clb) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_4clb <- FindClusters(so_4clb, resolution = 0.2)

DefaultAssay(so_4clb) <- "tRNA_aa"
so_4clb <- NormalizeData(so_4clb) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_4clb <- FindClusters(so_4clb, resolution = 0.2)

DefaultAssay(so_4clb) <- "RNA"
so_4clb <- NormalizeData(so_4clb) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_4clb <- FindClusters(so_4clb, resolution = 0.2)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so_4clb <- CellCycleScoring(so_4clb, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(so_4clb, group.by = c("seurat_clusters", "Phase"), label = TRUE)

saveRDS(so_4clb, file = "../../../trna_large_data/output/seurat_object_4cellline_Betaine_deepSeq_261019.rds")
