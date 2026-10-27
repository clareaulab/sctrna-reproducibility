library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(stringr)
library(readxl)
library(tibble)
library(BuenColors)

source("../../global_functions/tRNA_import.R")
source("../../global_functions/tRNA_valid.R")

so_p2 <- createSeurat(c("tRNAPS_FlexV2_Control", "tRNAPS_FlexV2_AAdepr1", "tRNAPS_FlexV2_AAdepr2", "tRNAPS_FlexV2_Torin1"),
                      c("Control", "AAdepr1", "AAdepr2", "Torin1"))
trna_raw_data <- import_tRNA(so_p2, c("tRNAPS_FlexV2_Control", "tRNAPS_FlexV2_AAdepr1", "tRNAPS_FlexV2_AAdepr2", "tRNAPS_FlexV2_Torin1"),
                             c("Control", "AAdepr1", "AAdepr2", "Torin1"))
so_p2 <- append_tRNA(so_p2, trna_raw_data)
dim(so_p2) # 18132 103853

# add metadata for sample
cell_names <- colnames(so_p2)
sample <- ifelse(grepl("-Control$", cell_names), "Control", 
               ifelse(grepl("-AAdepr1$", cell_names), "AAdepr1",
                      ifelse(grepl("-AAdepr2$", cell_names), "AAdepr2",
                             ifelse(grepl("-Torin1$", cell_names), "Torin1", NA))))
so_p2$sample <- sample

# filter tRNA probe pair data
df_valid <- tRNA_valid(so_p2)
so_p2[["tRNA_pp_valid"]] <- CreateAssayObject(df_valid[,colnames(so_p2)])

# normalize data
DefaultAssay(so_p2) <- "tRNA_pp_valid"
so_p2 <- NormalizeData(so_p2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p2 <- FindClusters(so_p2, resolution = 0.2)

DefaultAssay(so_p2) <- "tRNA_ac"
so_p2 <- NormalizeData(so_p2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_p2 <- FindClusters(so_p2, resolution = 0.2)

DefaultAssay(so_p2) <- "tRNA_aa"
so_p2 <- NormalizeData(so_p2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_p2 <- FindClusters(so_p2, resolution = 0.2)

DefaultAssay(so_p2) <- "RNA"
so_p2 <- NormalizeData(so_p2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p2 <- FindClusters(so_p2, resolution = 0.2)

table(so_p2$sample)
# AAdepr1 AAdepr2 Control  Torin1 
# 25322   25733   26448   26350 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so_p2 <- CellCycleScoring(so_p2, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(so_p2, group.by = c("seurat_clusters", "sample", "Phase"), label = TRUE)

saveRDS(so_p2, file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2.rds")
