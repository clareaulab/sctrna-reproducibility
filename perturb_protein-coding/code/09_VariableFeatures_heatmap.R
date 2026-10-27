library(dplyr)
library(Seurat)
library(ComplexHeatmap)

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_sgRNA_260901.rds")

# plot log2 FCs
guides <- c("GTF3C1_A", "GTF3C1_B", "GTF3C6_A", "GTF3C6_B", "ALKBH8_A",  "ALKBH8_B",
  "POLR3C_A", "POLR3C_B", "TRMT5_A", "TRMT5_B", "MTOR_A", "MTOR_B",
  "ELP3_A", "ELP4_A", "ELP3_B", "ELP4_B")
all_diff_out <- read.csv("../output/all_diff_out_tRNApp_260901.csv")
all_want <- all_diff_out %>% filter(gRNA %in% guides)
best_guides <- all_want %>%
  filter(p_val_adj < 10^-10) %>% filter(pct.2 > 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>% group_by(gRNA) %>% 
  top_n(n=2, wt = abs(avg_log2FC)) %>% pull(probe) %>% unique()
fc_mat <- all_want %>%
  dplyr::filter(probe %in% best_guides) %>%
  dplyr::select(probe, gRNA, avg_log2FC) %>%
  pivot_wider(names_from = gRNA,
    values_from = avg_log2FC,
    values_fill = 0)

fc_mat <- as.data.frame(fc_mat)
rownames(fc_mat) <- fc_mat$probe
fc_mat$probe <- NULL
fc_mat <- as.matrix(fc_mat)
fc_mat[fc_mat > 2] <- 2
fc_mat[fc_mat < -2] <- -2

my_colors <- jdb_palette("solar_flare", type = "continuous")
my_breaks <- seq(-2.05, 2.05, length.out = length(my_colors))
p <- pheatmap(fc_mat, color = my_colors, breaks = my_breaks)
cowplot::ggsave2(p, file = "../output/heatmap_variable-features_mtor2.pdf", width = 5.4, height = 4)
