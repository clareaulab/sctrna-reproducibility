library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(stringr)
library(readxl)
library(tibble)
select <- dplyr::select
filter <- dplyr::filter

source("../../global_functions/tRNA_import_massive.R")
source("../../global_functions/tRNA_valid_massive.R")

so_p <- createSeurat(c("tRNAPS_1","tRNAPS_2","tRNAPS_3","tRNAPS_4"), c("1","2","3","4"))
trna_raw_data <- import_tRNA(so_p, c("tRNAPS_1","tRNAPS_2","tRNAPS_3","tRNAPS_4"), c("1","2","3","4"))
so_p <- append_tRNA(so_p, trna_raw_data)
dim(so_p) #18099 89102

# add metadata for sample
cell_names <- colnames(so_p)
sample <- ifelse(grepl("-1$", cell_names), "1", 
               ifelse(grepl("-2$", cell_names), "2",
                      ifelse(grepl("-3$", cell_names), "3",
                             ifelse(grepl("-4$", cell_names), "4", NA))))
so_p$sample <- sample

# filter tRNA probe pair data
df_valid <- tRNA_valid(so_p)
so_p[["tRNA_pp_valid"]] <- CreateAssayObject(df_valid[,colnames(so_p)])

# normalize data
DefaultAssay(so_p) <- "tRNA_pp_valid"
so_p <- NormalizeData(so_p) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p <- FindClusters(so_p, resolution = 0.2)

DefaultAssay(so_p) <- "tRNA_ac"
so_p <- NormalizeData(so_p) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_p <- FindClusters(so_p, resolution = 0.2)

DefaultAssay(so_p) <- "tRNA_aa"
so_p <- NormalizeData(so_p) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_p <- FindClusters(so_p, resolution = 0.2)

DefaultAssay(so_p) <- "RNA"
so_p <- NormalizeData(so_p) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p <- FindClusters(so_p, resolution = 0.2)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so_p <- CellCycleScoring(so_p, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(so_p, group.by = c("seurat_clusters", "sample", "Phase"), label = TRUE)

saveRDS(so_p, file = "../../../trna_large_data/output/seurat_object_tRNAPS_260901.rds")
