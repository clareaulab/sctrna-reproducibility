library(dplyr)
library(Seurat)
library(BuenColors)

filter <- dplyr::filter
select <- dplyr::select

so_TcE <- readRDS("../../../trna_large_data/output/seurat_object_exhTcell_MM_filtered_0422.rds")

log2FC_cutoff <- 1
pval_cutoff <- 0.05

markers1 <- c("CD279", "CD155", "CD73", "IRF4")
markers2 <- c("CD3E", "CD3E.1", "CD56", "CD80", "CD86")
selenoprotein_markers <- c("TXNRD1", "GPX4")
p <- ggplot(diff, aes(x = avg_log2FC, y = -1 * log10(p_val_adj))) +
  geom_point(data = filter(diff, !(gene %in% c(markers1, markers2, selenoprotein_markers))),
             color = "lightgrey", alpha = 1.2) +
  geom_point(data = filter(diff, gene %in% markers1),
             color = "violetred4", size = 1.5) +
  geom_point(data = filter(diff, gene %in% markers2),
             color = "#5a88bf", size = 1.5) +
  geom_point(data = filter(diff, gene %in% selenoprotein_markers),
             color = "#009da6", size = 1.5) +
  geom_text_repel(
    data = filter(diff, gene %in% c(markers1, markers2)), 
    aes(label = gene), color = "black", size = 3,
    max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3,
    segment.color = "grey50", segment.size = 0.4
  ) + theme_minimal() + pretty_plot() + L_border() +
  labs(title = "Differential proteins in D8 chronic vs. acute",
       x = "log2 Fold Change", y = "-log10 Adjusted P-value")
cowplot::ggsave2(p, file = "../output/volcano_protein_D8c_D8a_label_reduced4.pdf", width = 3.5, height = 4)
