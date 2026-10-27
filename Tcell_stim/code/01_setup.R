library(dplyr)
library(Matrix)
library(Seurat)
library(stringr)
library(readxl)
library(BuenColors)

source("../../global_functions/tRNA_import_massive.R")
source("../../global_functions/tRNA_valid_massive.R")

so_TcE <- createSeurat(c("exhTcell_D2", "D8_acute", "D8_chronic"), c("D2", "D8_acute", "D8_chronic"))
trna_raw_data <- import_tRNA(so_TcE, c("exhTcell_D2", "D8_acute", "D8_chronic"), c("D2", "D8_acute", "D8_chronic"))
so_TcE <- append_tRNA(so_TcE, trna_raw_data)
dim(so_TcE) # 18099 54549

# add metadata for sample
cell_names <- colnames(so_TcE)
cond <- ifelse(grepl("-D2$", cell_names), "D2", 
                  ifelse(grepl("-D8_acute$", cell_names), "D8_acute", 
                         ifelse(grepl("-D8_chronic$", cell_names), "D8_chronic", NA)))
so_TcE$condition <- cond

# filter tRNA probe pair data
df_valid <- tRNA_valid(so_TcE) # Missing anticodons: UAA, UAG. Expected: UAA and UAG
so_TcE[["tRNA_pp_valid"]] <- CreateAssayObject(df_valid[,colnames(so_TcE)])

# normalize data
DefaultAssay(so_TcE) <- "tRNA_pp_valid"
so_TcE <- NormalizeData(so_TcE) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_TcE <- FindClusters(so_TcE, resolution = 0.2)

DefaultAssay(so_TcE) <- "tRNA_ac"
so_TcE <- NormalizeData(so_TcE) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_TcE <- FindClusters(so_TcE, resolution = 0.2)

DefaultAssay(so_TcE) <- "tRNA_aa"
so_TcE <- NormalizeData(so_TcE) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_TcE <- FindClusters(so_TcE, resolution = 0.2)

DefaultAssay(so_TcE) <- "RNA"
so_TcE <- NormalizeData(so_TcE) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_TcE <- FindClusters(so_TcE, resolution = 0.2)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so_TcE <- CellCycleScoring(so_TcE, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(so_TcE, group.by = c("seurat_clusters", "condition", "Phase"), label = TRUE)

saveRDS(so_TcE, file = "../../../trna_large_data/output/seurat_object_exhTcell.rds")
