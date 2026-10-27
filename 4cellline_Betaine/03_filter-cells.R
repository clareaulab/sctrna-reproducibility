library(dplyr)
library(Seurat)
library(BuenColors)

so_4clb <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline_Betaine_deepSeq_annotated_261019.rds")

col1 = c("HeLa"="darkseagreen", "HEK293T"="darkblue", "Jurkat"="violetred4", "K562"="cyan3")

# filter for cells that are left after clean-up for umap in tRNA space
cells_keep <- readLines("../data/4cellline_filtered_barcodes_keep.txt")
dim(so_4clb) # 18129 22682
so_4clb_filt <- so_4clb[,cells_keep]
dim(so_4clb_filt) # 18129 18528

# tRNA counts for manuscript text
colMeans(so_4clb_filt@meta.data[, c("nCount_RNA", "nCount_tRNA_pp_valid")])
# nCount_RNA nCount_tRNA_pp_valid 
# 2520.354            11321.873
so_4clb_filt@meta.data %>%
  summarise(across(c(nCount_RNA, nCount_tRNA_pp_valid),
                   list(mean = mean, median = median, sd = sd)))

DefaultAssay(so_4clb_filt) <- "RNA"
so_4clb_filt <- NormalizeData(so_4clb_filt) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_4clb_filt <- FindClusters(so_4clb_filt, resolution = 0.2)

p <- DimPlot(so_4clb_filt, reduction = "umap", group.by = "cell_line", label = FALSE, pt.size = 0.005) +
  scale_color_manual(values = col1) + theme_void() +  labs(title = NULL) + NoLegend()
cowplot::ggsave2(p, file = "../output/umap_cellline_figure1.png", width = 5.8, height = 5.2, dpi = 400)

saveRDS(so_4clb, file="../../../trna_large_data/output/seurat_object_4cellline_Betaine_deepSeq_annotated_261019_filt.rds")
