library(dplyr)
library(Matrix)
library(Seurat)

so_TcE <- readRDS(file = "../../../trna_large_data/output/seurat_object_exhTcell.rds")

# import protein
import_protein <- function(sample_id, barcode_idx = ""){
  mtx <- fread(paste0("../data/protein/",sample_id,"_Protein_full.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/protein/",sample_id,"_Protein_full.barcodes.txt"), header = FALSE)[[1]], "-", barcode_idx)
  colnames(matx) <- sub("\\..*", "", paste0(fread(paste0("../data/protein/",sample_id,"_Protein_full.genes.txt"), header = FALSE)[[1]]))
  return(t(matx))
}

#c("D2" ,"D8acute", "D8chronic"), c("D2", "D8_acute", "D8_chronic")
tcell_mtx <- cbind(import_protein("D2", "D2"),
                   import_protein("D8acute", "D8_acute"),
                   import_protein("D8chronic", "D8_chronic"))

# filter very low counts
tcell_mtx <- tcell_mtx[,colSums(tcell_mtx) >= 100]

# filter out very high / isotype control guys
tcell_mtx_meta <- data.frame(
  total = colSums(tcell_mtx),
  control_tags = colSums(tcell_mtx[grepl("sotypeCtrl", rownames(tcell_mtx)),]) + 1
)

ggplot(tcell_mtx_meta, aes(total, control_tags)) +
  geom_point() + scale_x_log10() + scale_y_log10()

remove_tcell_protein <- tcell_mtx_meta$total > 25000 | tcell_mtx_meta$control_tags > 75
sum(remove_tcell_protein) #518

tcell_mtx <- tcell_mtx[,!remove_tcell_protein]

saveRDS(tcell_mtx, file = "../output/tcell_mtx.rds")

rn_tcell_mtx2 <- make.unique(rownames(tcell_mtx))
tcell_mtx2 <- tcell_mtx
rownames(tcell_mtx2) <- rn_tcell_mtx2

# add protein data to seurat object
colnames(so_TcE)
cells_t2 <- intersect(colnames(tcell_mtx2), colnames(so_TcE)) #54103
dim(so_TcE) #18099 54549
so_TcE <- so_TcE[,cells_t2]
dim(so_TcE) #18099 54103
cells_t3 <- intersect(colnames(so_TcE), colnames(tcell_mtx2)) #54103
dim(tcell_mtx2) #344 66136
tcell_mtx2 <- tcell_mtx2[,cells_t3]
dim(tcell_mtx2) #344 54103
so_TcE[["Protein"]] <- CreateAssayObject(tcell_mtx2)

DefaultAssay(so_TcE) <- "Protein"
so_TcE <- NormalizeData(so_TcE) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_TcE <- FindClusters(so_TcE, resolution = 0.2)

DefaultAssay(so_TcE) <- "RNA"
so_TcE <- NormalizeData(so_TcE) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_TcE <- FindClusters(so_TcE, resolution = 0.2)

saveRDS(so_TcE, file="../../../trna_large_data/output/seurat_object_exhTcell_protein.rds")
