library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(stringr)
library(readxl)
library(tibble)
library(BuenColors)

so_4clb <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline_Betaine_deepSeq_261019.rds")

DimPlot(so_4clb, group.by = "seurat_clusters")

# look for mRNA markers
FeaturePlot(so_4clb, c("CXCR4", "CD3E", "TRBC2", # Jurkat
                       "S100A4", "CD44", "KRT8", # K562
                       "NCAM1", "TTC26", "SPI1", # HeLa
                       "GATA1", "LYZ", "BCR" # HEK293T
                       ))

Jurkat <- c(3,7)
K562 <- c(1,4,5)
HeLa <-  c(0)
HEK293T <- c(2,6)

so_4clb$cell_line <- case_when(
  so_4clb$seurat_clusters %in% Jurkat ~ "Jurkat",
  so_4clb$seurat_clusters %in% HeLa ~ "HeLa",
  so_4clb$seurat_clusters %in% HEK293T ~ "HEK293T",
  so_4clb$seurat_clusters %in% K562 ~ "K562",
  TRUE ~ "help"
)
DimPlot(so_4clb, group.by = c("seurat_clusters", "cell_line", "Phase"), label = TRUE)

col1 = c("HeLa"="#C2CC8E", "HEK293T"="darkblue", "Jurkat"="maroon", "K562"="cyan3")
DimPlot(so_4clb, reduction = "umap", group.by = "cell_line", label = FALSE) + 
  theme_void() +labs(title = NULL) + NoLegend() + scale_color_manual(values = col1)

saveRDS(so_4clb, file="../../../trna_large_data/output/seurat_object_4cellline_Betaine_deepSeq_annotated_261019.rds")
