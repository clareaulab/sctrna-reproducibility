library(Seurat)

so_4cl <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline.rds")

DimPlot(so_4cl, group.by = "seurat_clusters")

# look for mRNA markers
FeaturePlot(so_4cl, c("CXCR4", "CD3E", "TRBC2", # Jurkat
                       "S100A4", "CD44", "KRT8", # K562
                       "NCAM1", "TTC26", "SPI1", # HeLa
                       "GATA1", "LYZ", "BCR" # HEK293T
                       ))

Jurkat <- c(1,6)
HeLa <-  c(0,5)
HEK293T <- c(3,7,8)
K562 <- c(2,4)

so_4cl$cell_line <- case_when(
  so_4cl$seurat_clusters %in% Jurkat ~ "Jurkat",
  so_4cl$seurat_clusters %in% HeLa ~ "HeLa",
  so_4cl$seurat_clusters %in% HEK293T ~ "HEK293T",
  so_4cl$seurat_clusters %in% K562 ~ "K562",
  TRUE ~ "help"
)

DimPlot(so_4cl, group.by = c("seurat_clusters", "cell_line", "Phase"), label = TRUE)

col1 = c("HeLa"="darkseagreen", "HEK293T"="darkblue", "Jurkat"="maroon", "K562"="cyan3")

DimPlot(so_4cl, reduction = "umap", group.by = "cell_line", label = FALSE) + 
  theme_void() + labs(title = NULL) + NoLegend() + scale_color_manual(values = col1)

saveRDS(so_4cl, file="../../../trna_large_data/output/seurat_object_4cellline_annotated.rds")
