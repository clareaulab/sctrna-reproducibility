library(Seurat)

so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_FDL_260907.rds")

# transcriptomic markers
DefaultAssay(so_cd34) <- "RNA"
p <- FeaturePlot(so_cd34, features = c("CD34"), cols = c("lightgray", "violetred4"),
            min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_CD34.png", width = 8, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("GATA1"), cols = c("lightgray", "violetred4"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_GATA1.png", width = 8, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("CD14"), cols = c("lightgray", "violetred4"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_CD14.png", width = 8, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("PF4"), cols = c("lightgray", "violetred4"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_PF4.png", width = 8, height = 7.2, dpi = 400)


# tRNA markers
DefaultAssay(so_cd34) <- "tRNA_pp_valid"
p <- FeaturePlot(so_cd34, features = c("PP94-Leu-UUG"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp94_Leu_UUG.png", width = 8, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("PP163-Val-GUC"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp163_Val_GUC.png", width = 8, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("PP56-Cys-UGC"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp56_Cys_UGC.png", width = 8, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("PP82-Ile-AUC"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp82_Ile_AUC.png", width = 8, height = 7.2, dpi = 400)

# additional tRNA markers for supplement
DefaultAssay(so_cd34) <- "tRNA_pp_valid"
p <- FeaturePlot(so_cd34, features = c("PP44-Asp-GAC"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp44_Asp_GAC.png", width = 8.4, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("PP45-Asp-GAC"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp45_Asp_GAC.png", width = 8.4, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("PP120-Met-AUG"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp120_Met_AUG.png", width = 8.4, height = 7.2, dpi = 400)

p <- FeaturePlot(so_cd34, features = c("PP172-iMet-AUG"), cols = c("lightgray", "darkblue"),
                 min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
cowplot::ggsave2(p, file = "../output/umap_pp172_iMet_AUG.png", width = 8.4, height = 7.2, dpi = 400)
