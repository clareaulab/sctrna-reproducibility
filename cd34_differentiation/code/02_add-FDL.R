library(dplyr)
library(Seurat)

so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_260907.rds")
fdl <- read.csv("../data/FDL_coordinates.csv")

colnames(so_cd34)
tail(colnames(so_cd34))

fdl$barcode <- paste0(
  substr(fdl$X, 1, 16), 
  "-", 
  sub(".*-", "", fdl$X)
)
fdl$barcode <- gsub("ay", "", fdl$barcode)

rownames(fdl) <- fdl$barcode
fdl$barcode <- NULL
fdl$X <- NULL
so_cd34 <- AddMetaData(so_cd34, metadata = fdl)

so_cd34 <- subset(so_cd34, subset = !is.na(FA1) & !is.na(FA2)) # removing 3440 cells
dim(so_cd34) #18099 72544

ggplot(so_cd34@meta.data, aes(x = FA1, y = FA2, color = seurat_clusters)) +
  geom_point(size = 0.5, alpha = 0.8) +
  theme_minimal() +
  labs(x = "FDL1", y = "FDL2", color = "Cluster")

fdl_coords <- as.matrix(so_cd34@meta.data[, c("FA1", "FA2")])
rownames(fdl_coords) <- colnames(so_cd34)

so_cd34[["fdl"]] <- CreateDimReducObject(
  embeddings = fdl_coords,
  key = "FDL_",
  assay = DefaultAssay(so_cd34)
)

saveRDS(so_cd34, file = "../../../trna_large_data/output/seurat_object_CD34_diff_FDL_260907.rds")


# export FDL layout colored by sample
col1 <- c("D0" = "darkblue", "D7_Ery" = "firebrick", "D7_Mono" = "lightblue",
          "D7_Mega" = "mediumpurple")
p <- DimPlot(so_cd34, reduction = "fdl", group.by = "condition", cols = col1)
cowplot::ggsave2(p, file = "../output/umap_FDL_conditions.png", width = 9, height = 8, dpi = 400)

# umap colored by tRNA counts
p <- FeaturePlot(so_cd34, "nCount_tRNA_pp_valid", cols = c("lightgray", "darkblue"), reduction = "fdl",
            min.cutoff = "q02", max.cutoff = "q98")
cowplot::ggsave2(p, file = "../output/umap_counts_tRNApp.png", width = 8.4, height = 7.2, dpi = 400)
