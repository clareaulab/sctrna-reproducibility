library(dplyr)
library(Seurat)
filter <- dplyr::filter
select <- dplyr::select
"%ni%" <- Negate("%in%")

so_TcE <- readRDS(file = "../../../trna_large_data/output/seurat_object_exhTcell_protein.rds")
so_TcE <- JoinLayers(so_TcE)
so_TcE@reductions

# mRNA
DefaultAssay(so_TcE) <- "RNA"
so_TcE <- NormalizeData(so_TcE)
so_TcE <- FindVariableFeatures(so_TcE)
so_TcE <- ScaleData(so_TcE)
so_TcE <- RunPCA(so_TcE, reduction.name = "pca_mRNA")

# tRNA
DefaultAssay(so_TcE) <- "tRNA_pp_valid"
so_TcE <- NormalizeData(so_TcE)
so_TcE <- FindVariableFeatures(so_TcE)
so_TcE <- ScaleData(so_TcE)
so_TcE <- RunPCA(so_TcE, reduction.name = "pca_tRNA")

# protein
DefaultAssay(so_TcE) <- "Protein"
so_TcE <- NormalizeData(so_TcE, normalization.method = "CLR", margin = 2)
so_TcE <- ScaleData(so_TcE)
so_TcE <- RunPCA(so_TcE, reduction.name = "pca_protein")

so_TcE@reductions
Reductions(so_TcE)

# Now run multimodal neighbors and embedding
so_TcE <- FindMultiModalNeighbors(object = so_TcE,
                                    reduction.list = list("pca_mRNA", "pca_tRNA", "pca_protein"),
                                    dims.list = list(1:30, 1:15, 1:15))
so_TcE <- RunUMAP(so_TcE, nn.name = "weighted.nn", reduction.name = "wnn.3.umap", reduction.key = "Uw3_" )
so_TcE <- FindClusters(so_TcE, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)

DefaultAssay(so_TcE) <- "Protein"
FeaturePlot(so_TcE, features = c("nCount_RNA", "nCount_tRNA_pp_valid", "nCount_Protein"), min.cutoff = "q0", max.cutoff = "q90",
            reduction =  'wnn.3.umap',   pt.size = 0.01, ncol=3)
FeaturePlot(so_TcE, features = c("CD19", "CD20", "MS4A1"), min.cutoff = "q0", max.cutoff = "q90",
            reduction =  'wnn.3.umap',   pt.size = 0.01, ncol=3)
FeaturePlot(so_TcE, features = c("CD4", "CD8A", "IL7R"), min.cutoff = "q0", max.cutoff = "q90",
            reduction =  'wnn.3.umap',   pt.size = 0.01, ncol=3)

# filter
so_TcE_filt <- subset(so_TcE, subset = seurat_clusters %ni% c(22, 27, 28))
DimPlot(so_TcE_filt, group.by = c("condition", "seurat_clusters"), label = T)

so_TcE_filt <- FindMultiModalNeighbors(object = so_TcE_filt,
                                  reduction.list = list("pca_mRNA", "pca_tRNA", "pca_protein"),
                                  dims.list = list(1:30, 1:15, 1:15))
so_TcE_filt <- RunUMAP(so_TcE_filt, nn.name = "weighted.nn", reduction.name = "wnn.3.umap", reduction.key = "Uw3_" )
so_TcE_filt <- FindClusters(so_TcE_filt, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)

so_TcE_filt$nCount_norm <- so_TcE_filt$nCount_Protein/(so_TcE_filt$nCount_Protein + so_TcE_filt$nCount_RNA +
                                                          so_TcE_filt$nCount_tRNA_pp_valid)
FeaturePlot(so_TcE_filt, features = c("nCount_RNA", "nCount_tRNA_pp_valid", "nCount_Protein", "nCount_norm"), min.cutoff = "q0", max.cutoff = "q90",
            reduction =  'wnn.3.umap', pt.size = 0.01, ncol=3)

#saveRDS(so_TcE_filt, "../../../trna_large_data/output/seurat_object_exhTcell_MM_filtered_0422.rds")
