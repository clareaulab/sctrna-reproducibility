library(dplyr)
library(Seurat)

source("../../global_functions/tRNA_import_massive.R")
source("../../global_functions/tRNA_valid_massive.R")

so_4cl <- createSeurat(c("4CL"), c("1"))
trna_raw_data <- import_tRNA(so_4cl, c("4CL"), c("1"))
so_4cl <- append_tRNA(so_4cl, trna_raw_data)
dim(so_4cl) # 18099 18821

# filter tRNA probe pair data
df_valid <- tRNA_valid(so_4cl) # Missing anticodons: GAU, UAA and UAG. Expected: UAA and UAG
so_4cl[["tRNA_pp_valid"]] <- CreateAssayObject(df_valid[,colnames(so_4cl)])

# normalize data
DefaultAssay(so_4cl) <- "tRNA_pp_valid"
so_4cl <- NormalizeData(so_4cl) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_4cl <- FindClusters(so_4cl, resolution = 0.2)

DefaultAssay(so_4cl) <- "tRNA_ac"
so_4cl <- NormalizeData(so_4cl) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_4cl <- FindClusters(so_4cl, resolution = 0.2)

DefaultAssay(so_4cl) <- "tRNA_aa"
so_4cl <- NormalizeData(so_4cl) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_4cl <- FindClusters(so_4cl, resolution = 0.2)

DefaultAssay(so_4cl) <- "RNA"
so_4cl <- NormalizeData(so_4cl) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_4cl <- FindClusters(so_4cl, resolution = 0.2)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so_4cl <- CellCycleScoring(so_4cl, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(so_4cl, group.by = c("seurat_clusters", "Phase"), label = TRUE)

saveRDS(so_4cl, file = "../../../trna_large_data/output/seurat_object_4cellline.rds")
