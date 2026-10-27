library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
library(pheatmap)

set.seed(1234)
filter <- dplyr::filter
select <- dplyr::select

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_sgRNA_260901.rds")
so_p$target <- str_remove(so_p$identity_NTC, "_[AB]$")
table(so_p$target)
dim(so_p) #18099 89088
so_p <- subset(so_p, subset = target != "TRMT5")
dim(so_p) #18099 88351

## this plot for mRNA
guides_trnaps <- read_excel("../data/opool_Twist_Perturb-tRNA_manual_250707 copy.xlsx")
guides_chrom <- guides_trnaps[193:260,1]
guides_chrom_vec <- as.character(guides_chrom[[1]])
guides_chrom_vec <- unique(gsub("-[AB]","",guides_chrom_vec))
guides_chrom_vec_ntc <- c(guides_chrom_vec, "NTC")
all_guides <- unique(so_p$target)
guides_get <- setdiff(all_guides, guides_chrom_vec_ntc)

so_p <- NormalizeData(so_p)
so_p <- FindVariableFeatures(so_p, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(so_p)
groups <- so_p[["target"]][,1]
#groups_use <- groups[groups %in% guides_get]
counts <- so_p[["RNA"]]$counts[hvgs,]
split_cells <- split(colnames(so_p), groups)
pb_rna <- lapply(split_cells, function(cells) {
  if (length(cells) == 0) return(NULL)
  Matrix::rowSums(counts[, cells, drop = FALSE]) / length(cells)
})
# convert each pseudobulk vector to a one-column matrix
pb_df_rna <- as.data.frame(do.call(rbind, pb_rna))
pb_df_rna <- t(pb_df_rna)
pb_df_rna2 <- data.frame(pb_df_rna)
pb_df_rna2$cell_id <- seq_len(nrow(pb_df_rna2))
pb_df_collapsed <- pb_df_rna2 %>%
  pivot_longer(cols = -cell_id, names_to = "guide", values_to = "value") %>%
  dplyr::mutate(prefix = sub("_.*", "", guide)) %>%
  group_by(cell_id, prefix) %>%
  dplyr::summarize(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = prefix, values_from = mean_value)
pb_df_collapsed <- pb_df_collapsed %>% 
  select(-cell_id)
pb_z_rna <- t(scale(t(pb_df_collapsed)))
sim_rna <- cor(pb_z_rna, use="pairwise.complete.obs")  # guides x guides
dist_mat_rna <- as.dist(1 - sim_rna)
sim_rna <- as.matrix(sim_rna)
sim_rna[sim_rna > 1] <- 1
sim_rna[sim_rna < -1] <- -1
my_colors <- jdb_palette("solar_flare", type = "continuous")
my_breaks <- seq(-1.05, 1.05, length.out = length(my_colors))
pheatmap(sim_rna, 
         color = my_colors, 
         breaks = my_breaks)
pheatmap::pheatmap(sim_rna, 
                   clustering_distance_rows = dist_mat_trna,
                   clustering_distance_cols = dist_mat_trna,
                   color = my_colors, 
                   breaks = my_breaks)
sim_rna_sub <- sim_rna[guides_get, guides_get, drop = FALSE]
p <- pheatmap::pheatmap(sim_rna_sub, color = my_colors, breaks = my_breaks)
cowplot::ggsave2(p, file = "../output/heatmap_VariableFetaures_mRNA_260902.pdf", width = 20, height = 20)

# save row/column order for tRNA heatmap
row_order <- p$tree_row$order
col_order <- p$tree_col$order
row_names_ordered <- rownames(sim_rna_sub)[row_order]
col_names_ordered <- colnames(sim_rna_sub)[col_order]


## for tRNA
so_p <- NormalizeData(so_p, assay = "tRNA_pp_valid")
so_p <- FindVariableFeatures(so_p, selection.method = "vst", nfeatures = 150,
                             assay = "tRNA_pp_valid")
hvgs <- VariableFeatures(so_p, assay = "tRNA_pp_valid")
groups <- so_p[["target"]][,1]
#groups_use <- groups[groups %in% guides_get]
counts <- so_p[["tRNA_pp_valid"]]$counts[hvgs,]
split_cells <- split(colnames(so_p), groups)
pb_trna <- lapply(split_cells, function(cells) {
  if (length(cells) == 0) return(NULL)
  Matrix::rowSums(counts[, cells, drop = FALSE]) / length(cells)
})
# convert each pseudobulk vector to a one-column matrix
pb_df_trna <- as.data.frame(do.call(rbind, pb_trna))
pb_df_trna <- t(pb_df_trna)
pb_df_trna2 <- data.frame(pb_df_trna)
pb_df_trna2$cell_id <- seq_len(nrow(pb_df_trna2))
pb_df_collapsed <- pb_df_trna2 %>%
  pivot_longer(cols = -cell_id, names_to = "guide", values_to = "value") %>%
  dplyr::mutate(prefix = sub("_.*", "", guide)) %>%
  group_by(cell_id, prefix) %>%
  dplyr::summarize(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = prefix, values_from = mean_value)
pb_df_collapsed <- pb_df_collapsed %>% 
  select(-cell_id)
pb_z_trna <- t(scale(t(pb_df_collapsed)))
sim_trna <- cor(pb_z_trna, use="pairwise.complete.obs")  # guides x guides
dist_mat_trna <- as.dist(1 - sim_trna)
sim_trna <- as.matrix(sim_trna)
sim_trna[sim_trna > 1] <- 1
sim_trna[sim_trna < -1] <- -1
my_colors <- jdb_palette("solar_flare", type = "continuous")
my_breaks <- seq(-1.05, 1.05, length.out = length(my_colors))
#diag(sim_trna) <- NA
pheatmap(sim_trna, 
         color = my_colors, 
         breaks = my_breaks)
pheatmap::pheatmap(sim_trna, 
                   clustering_distance_rows = dist_mat_trna,
                   clustering_distance_cols = dist_mat_trna,
                   color = my_colors, 
                   breaks = my_breaks)
sim_trna_sub <- sim_trna[guides_get, guides_get, drop = FALSE]
p <- pheatmap::pheatmap(sim_trna_sub, color = my_colors, breaks = my_breaks)
#cowplot::ggsave2(p, file = "../output/heatmap_VariableFetaures_tRNA_200features_freeClustering_260902.pdf", width = 20, height = 20)

sim_trna_ordered <- sim_trna[row_names_ordered, col_names_ordered, drop = FALSE]
p <- pheatmap::pheatmap(
  sim_trna_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = my_colors,
  breaks = my_breaks
)
cowplot::ggsave2(p, file = "../output/heatmap_VariableFetaures_tRNA_150features_260902.pdf", width = 20, height = 20)
