library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)

filter <- dplyr::filter
select <- dplyr::select

sample_colors <- c(
  Con = "#fdde8f",
  Tor = "#4FBBD6",
  AAd = "#D65147"
)

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
so_p2 <- subset(so_p2, subset = target != "TRMT5")
table(so_p2$identity_NTC)
so_p2$sample2 <- substr(so_p2$sample, 1, 3)
table(so_p2$sample2)

# embedding all three conditions in one umap
pseudobulk_by_target_sample <- function(so, assay = "RNA") {
  mat <- as.matrix(so@assays[[assay]]$counts)
  targets <- so@meta.data$target
  samples <- so@meta.data$sample2
  group <- paste(targets, samples, sep = "_")
  pb <- rowsum(t(mat), group = group)
  pb <- t(pb)
  return(pb)
}


## embedding all three conditions, mRNA
so_p2 <- JoinLayers(so_p2)
pb_rna <- pseudobulk_by_target_sample(so_p2, "RNA")

guide_sample_meta <- data.frame(
  guide_sample = colnames(pb_rna)
)
guide_sample_meta$target <- sub("_(.*)", "", guide_sample_meta$guide_sample)
guide_sample_meta$sample2 <- sub(".*_", "", guide_sample_meta$guide_sample)

# get CPM / cp10k
pb_rna_norm <- t(t(pb_rna) / colSums(pb_rna) * 1e6)
pb_rna_norm <- log1p(pb_rna_norm)
rna_mat <- t(pb_rna_norm)

gene_var <- apply(rna_mat, 2, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:500]
rna_mat_var <- rna_mat[, top_genes]
rna_pca <- prcomp(rna_mat_var, scale. = TRUE)

#rna_pca <- prcomp(rna_mat, scale. = TRUE)
rna_umap <- umap(rna_pca$x[,1:30])
rna_umap_df <- data.frame(
  UMAP1 = rna_umap[,1],
  UMAP2 = rna_umap[,2],
  guide_sample_meta
)

p <- ggplot(rna_umap_df, aes(UMAP1, UMAP2)) +
  geom_point(data = rna_umap_df %>% dplyr::filter(sample2 %in% c("AAd", "Tor")),
             aes(fill = sample2), shape = 21, size = 1.8, color = "darkgrey") +
  geom_point(data = rna_umap_df %>% dplyr::filter(sample2 == "Con"),
             aes(fill = sample2), shape = 21, size = 1.8, color = "darkgrey", stroke = 0.4) +
  scale_fill_manual(values = sample_colors) + pretty_plot() + L_border()
ggsave(filename = "../output/umap_guide_embedding_mRNA.pdf", plot = p, width = 5.5, height = 5, units = "in")


## embedding all three conditions, tRNA
pb_rna <- pseudobulk_by_target_sample(so_p2, "tRNA_pp_valid")

guide_sample_meta <- data.frame(guide_sample = colnames(pb_rna))
guide_sample_meta$target <- sub("_(.*)", "", guide_sample_meta$guide_sample)
guide_sample_meta$sample2 <- sub(".*_", "", guide_sample_meta$guide_sample)

# get CPM
pb_rna_norm <- t(t(pb_rna) / colSums(pb_rna) * 1e6)
pb_rna_norm <- log1p(pb_rna_norm)

rna_mat <- t(pb_rna_norm)
gene_var <- apply(rna_mat, 2, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:30]
rna_mat_var <- rna_mat[, top_genes]
rna_pca <- prcomp(rna_mat_var, scale. = TRUE)

rna_umap <- umap(rna_pca$x[,1:15])
rna_umap_df <- data.frame(
  UMAP1 = rna_umap[,1],
  UMAP2 = rna_umap[,2],
  guide_sample_meta
)
ggplot(rna_umap_df, aes(UMAP1, UMAP2, color = sample2)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = target), vjust = -0.6, size=2) + pretty_plot() + L_border()

p <- ggplot(rna_umap_df, aes(UMAP1, UMAP2)) +
  geom_point(data = rna_umap_df %>% dplyr::filter(sample2 %in% c("AAd", "Tor")),
             aes(fill = sample2), shape = 21, size = 1.8, color = "darkgrey") +
  geom_point(data = rna_umap_df %>% dplyr::filter(sample2 == "Con"),
             aes(fill = sample2), shape = 21, size = 1.8, color = "darkgrey", stroke = 0.4) +
  scale_fill_manual(values = sample_colors) + pretty_plot() + L_border()
ggsave(filename = "../output/umap_guide_embedding_tRNA.pdf", plot = p, width = 5.5, height = 5, units = "in")

p <- ggplot(rna_umap_df, aes(UMAP1, UMAP2)) +
  geom_point(data = rna_umap_df %>% dplyr::filter(sample2 %in% c("AAd", "Tor")),
             aes(fill = sample2), shape = 21, size = 3, color = "darkgrey") +
  geom_point(data = rna_umap_df %>% dplyr::filter(sample2 == "Con"),
             aes(fill = sample2), shape = 21, size = 3, color = "darkgrey", stroke = 0.4) +
  scale_fill_manual(values = sample_colors) + pretty_plot() + L_border() +
  geom_text(data = rna_umap_df %>% dplyr::filter(sample2 == "Con"),
            aes(label = target), vjust = -0.6, size=4) + pretty_plot() + L_border()
ggsave(filename = "../output/umap_guide_embedding_tRNA_label.pdf", plot = p, width = 15.5, height = 15, units = "in")
