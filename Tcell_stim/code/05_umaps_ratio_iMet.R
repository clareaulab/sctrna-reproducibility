library(dplyr)
library(Seurat)

so_TcE <- readRDS(file = "../../../trna_large_data/output/seurat_object_exhTcell_MM_filtered_0422.rds")

# per-cell tRNA/mRNA ratio
so_TcE$tRNA_mRNA_ratio <- log1p(so_TcE$nCount_tRNA_pp_valid)/log1p(so_TcE$nCount_RNA)

FeaturePlot(so_TcE, features = "tRNA_mRNA_ratio", min.cutoff = "q02", max.cutoff = "q98",
            reduction = "wnn.3.umap") +
  scale_colour_viridis_c(option = "viridis")
#cowplot::ggsave2(p, file = "../output/umap_translational-capacity_viridis.png", width = 6.5, height = 5.2, dpi = 400)
#cowplot::ggsave2(p, file = "../output/umap_translational-capacity_viridis.pdf", width = 6.5, height = 5.2, dpi = 400)

# plot relative umap for iMet
imet_counts <- GetAssayData(so_TcE, assay = "tRNA_aa", layer = "counts")["iMet", ]
so_TcE$reliMet <- imet_counts / so_TcE$nCount_tRNA_pp_valid
p <- FeaturePlot(so_TcE, features = "reliMet", max.cutoff = "q95", reduction = "wnn.3.umap",
                 cols = c("lightgrey", "darkblue"), pt.size = 0.01)
cowplot::ggsave2(p, file = "../output/umap_reliMet2.png", width = 6.5, height = 5.8)

# quantify iMet fold change
res <- FindMarkers(so_TcE, assay = "tRNA_aa", group.by = "condition",
                   ident.1 = "D8_acute", ident.2 = "D2",
                   features = "iMet", logfc.threshold = 0, min.pct = 0)
#p_val avg_log2FC pct.1 pct.2 p_val_adj
#iMet     0  0.7134122 0.983 0.705         0
format(res$p_val, digits = 3)
res <- FindMarkers(so_TcE, assay = "tRNA_aa", group.by = "condition",
                   ident.1 = "D8_acute", ident.2 = "D8_chronic",
                   features = "iMet", logfc.threshold = 0, min.pct = 0)
#p_val avg_log2FC pct.1 pct.2 p_val_adj
#iMet     0   1.124556 0.983  0.94         0
format(res$p_val, digits = 3)

df <- LayerData(so_TcE, assay = "tRNA_aa", layer = "data")["iMet", ]
df <- so_TcE@assays$tRNA_aa$counts["iMet",]
grp <- so_TcE$condition
wilcox_log10p(df[grp == "D8_acute"], df[grp == "D2"]) # -1259.084
wilcox_log10p(df[grp == "D8_acute"], df[grp == "D8_chronic"]) # -3400.926
