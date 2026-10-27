library(dplyr)
library(Seurat)

so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_FDL_260907.rds")

DefaultAssay(so_cd34) <- "RNA"
trna_transcription <- c("POLR3A","POLR3B","POLR3C","POLR3D", "POLR3E","POLR3F",
                        "POLR3G","POLR3H","POLR3K", "BDP1", "GTF3C1","GTF3C2",
                        "GTF3C3","GTF3C4","GTF3C5","GTF3C6", "GTF3A")

# Ucell module score
so_cd34 <- AddModuleScore_UCell(so_cd34, features = list(trna_transcription = trna_transcription))
p <- FeaturePlot(so_cd34, "trna_transcription_UCell", cols = c("lightgrey", "violetred4"),
            min.cutoff = "q05", max.cutoff = "q95", reduction =  'fdl')
cowplot::ggsave2(p, file = "../output/umap_module_tRNAtranscription_Ucell.png", width = 8.4, height = 7.2, dpi = 400)
