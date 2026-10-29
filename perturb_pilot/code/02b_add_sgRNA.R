library(dplyr)
library(Seurat)

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_Perturb_pilot.rds")

sgRNA_qc <- read.table(
  #file = "../data/PerturbPilot_tRNA_gRNA_calls.tsv",
  file = "../output/PerturbPilot_tRNA_gRNA_calls_261027.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
cells_qc <- sgRNA_qc$barcode
so_p <- so_p[,cells_qc]
dim(so_p) # 18099 10899 # 18099 20081

sparse_sgRNA <- sparseMatrix(i = as.factor(sgRNA_qc$max_id), j = as.factor(sgRNA_qc$barcode), x = sgRNA_qc$total_umi)
rownames(sparse_sgRNA) <- levels(as.factor(sgRNA_qc$max_id))
colnames(sparse_sgRNA) <- levels(as.factor(sgRNA_qc$barcode))
cells <- intersect(colnames(sparse_sgRNA),colnames(so_p))
sparse_sgRNA <- sparse_sgRNA[,cells]
so_p[["sgRNA"]] <- CreateAssayObject(sparse_sgRNA[,cells])

identity_vec <- sgRNA_qc$max_id; names(identity_vec) <- sgRNA_qc$barcode
so_p@meta.data$identity <- unname(identity_vec[colnames(so_p)])
so_p@meta.data$identity_NTC <- case_when(
  grepl("targeting", so_p@meta.data$identity) ~ "NTC",
  TRUE ~ so_p@meta.data$identity
)
head(so_p@meta.data)

DefaultAssay(so_p) <- "sgRNA"
so_p <- NormalizeData(so_p) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p <- FindClusters(so_p, resolution = 0.2)

DefaultAssay(so_p) <- "tRNA_pp_valid"
so_p <- NormalizeData(so_p) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p <- FindClusters(so_p, resolution = 0.2)

DefaultAssay(so_p) <- "RNA"
so_p <- NormalizeData(so_p) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_p <- FindClusters(so_p, resolution = 0.2)

#saveRDS(so_p, file = "../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA.rds")
saveRDS(so_p, file = "../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA_261027.rds")
