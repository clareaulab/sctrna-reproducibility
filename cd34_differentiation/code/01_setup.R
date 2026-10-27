library(dplyr)
library(Seurat)

source("../../global_functions/tRNA_import_massive.R")
source("../../global_functions/tRNA_valid_massive.R")

so_cd34 <- createSeurat(c("CD34_D0", "CD34_D7_Ery", "CD34_D7_Mono", "CD34_D7_Mega"), c("D0", "D7_Ery", "D7_Mono", "D7_Mega"))
dim(so_cd34) #18099 76245
trna_raw_data <- import_tRNA(so_cd34, c("CD34_D0", "CD34_D7_Ery", "CD34_D7_Mono", "CD34_D7_Mega"), c("D0", "D7_Ery", "D7_Mono", "D7_Mega"))
so_cd34 <- append_tRNA(so_cd34, trna_raw_data)
dim(so_cd34) #18099 75984

# add metadata for sample
cell_names <- colnames(so_cd34)
cond <- ifelse(grepl("-D0$", cell_names), "D0", 
                  ifelse(grepl("-D7_Ery$", cell_names), "D7_Ery",
                         ifelse(grepl("-D7_Mono$", cell_names), "D7_Mono",
                                ifelse(grepl("-D7_Mega$", cell_names), "D7_Mega", NA))))
so_cd34$condition <- cond

# filter tRNA probe pair data
df_valid <- tRNA_valid(so_cd34) # Missing anticodons: UAA, UAG. Expected: UAA and UAG
so_cd34[["tRNA_pp_valid"]] <- CreateAssayObject(df_valid[,colnames(so_cd34)])
dim(so_cd34) #18099 60679

# normalize data
DefaultAssay(so_cd34) <- "tRNA_pp_valid"
so_cd34 <- NormalizeData(so_cd34) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_cd34 <- FindClusters(so_cd34, resolution = 0.2)

DefaultAssay(so_cd34) <- "tRNA_ac"
so_cd34 <- NormalizeData(so_cd34) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_cd34 <- FindClusters(so_cd34, resolution = 0.2)

DefaultAssay(so_cd34) <- "tRNA_aa"
so_cd34 <- NormalizeData(so_cd34) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:20)
so_cd34 <- FindClusters(so_cd34, resolution = 0.2)

DefaultAssay(so_cd34) <- "RNA"
so_cd34 <- NormalizeData(so_cd34) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_cd34 <- FindClusters(so_cd34, resolution = 0.2)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so_cd34 <- CellCycleScoring(so_cd34, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(so_cd34, group.by = c("seurat_clusters", "condition", "Phase"), label = TRUE)

#saveRDS(so_cd34, file = "../../../trna_large_data/output/seurat_object_CD34diff_260907.rds")
