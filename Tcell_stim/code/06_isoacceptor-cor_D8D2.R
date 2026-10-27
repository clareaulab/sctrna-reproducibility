library(dplyr)
library(Seurat)
library(ggplot2)
library(BuenColors)

filter <- dplyr::filter
select <- dplyr::select

so_TcE <- readRDS(file = "../../../trna_large_data/output/seurat_object_exhTcell_MM_filtered_0422.rds")
so_TcE <- JoinLayers(so_TcE, assay="RNA")

run_fm <- function(ident1, assay) {
  FindMarkers(so_TcE, assay = assay,
              ident.1 = ident1, ident.2 = "D2",
              logfc.threshold = 0, min.pct = 0) %>%
    tibble::rownames_to_column("feature") %>%
    mutate(pass = p_val_adj <= 0.05 & pmax(pct.1, pct.2) >= 0.1)
}

trna_a <- run_fm("D8_acute",   "tRNA_aa")
trna_c <- run_fm("D8_chronic", "tRNA_aa")

trna_fc <- inner_join(
  trna_a %>% select(feature, lfc_acute   = avg_log2FC, pass_acute   = pass),
  trna_c %>% select(feature, lfc_chronic = avg_log2FC, pass_chronic = pass),
  by = "feature") %>%
  filter(pass_acute | pass_chronic)

rho <- cor(trna_fc$lfc_chronic, trna_fc$lfc_acute, method = "spearman")
cor(trna_fc$lfc_chronic, trna_fc$lfc_acute, method = "spearman") # 0.2072276
cor(trna_fc$lfc_chronic, trna_fc$lfc_acute, method = "pearson") # 0.2534636

highlight <- c("SeC", "iMet")

plot_df <- trna_fc %>%
  mutate(group = ifelse(grepl(paste(highlight, collapse = "|"), feature),
                        "highlight", "other"))

p <- ggplot(plot_df, aes(x = lfc_chronic, y = lfc_acute)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = group), size = 3) +
  geom_text_repel(aes(label = feature, color = group), size = 3.2,
                  max.overlaps = Inf, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, size = 3.5,
           label = paste0("rho = ", round(rho, 2))) +
  scale_color_manual(values = c(highlight = "violetred4", other = "grey50")) +
  labs(x = expression(log[2]~"FC D8 chronic / D2"),
       y = expression(log[2]~"FC D8 acute / D2")) +
  xlim(-4.5, 4.1) + ylim(-4.5, 4.1) +
  pretty_plot() + L_border() +
  theme(legend.position = "none")
cowplot::ggsave2(p, file = "../output/cor_tRNAaa2.pdf", width = 4.2, height = 4, dpi = 400)

# quantify SeC fold change
res <- FindMarkers(so_TcE, assay = "tRNA_aa", group.by = "condition",
                   ident.1 = "D8_acute", ident.2 = "D2",
                   features = "SeC", logfc.threshold = 0, min.pct = 0)
#    p_val avg_log2FC pct.1 pct.2 p_val_adj
#SeC     0   2.161247 0.947 0.226         0
format(res$p_val, digits = 3)
res <- FindMarkers(so_TcE, assay = "tRNA_aa", group.by = "condition",
                   ident.1 = "D8_chronic", ident.2 = "D2",
                   features = "SeC", logfc.threshold = 0, min.pct = 0)
#    p_val avg_log2FC pct.1 pct.2 p_val_adj
#SeC     0   3.605195 0.969 0.226         0
format(res$p_val, digits = 3)

# relative umap for SeC
sec_counts <- GetAssayData(so_TcE, assay = "tRNA_aa", layer = "counts")["SeC", ]
so_TcE$relSeC <- sec_counts / so_TcE$nCount_tRNA_pp_valid
p <- FeaturePlot(so_TcE, features = "relSeC", max.cutoff = "q95", reduction = "wnn.3.umap",
                 cols = c("lightgrey", "darkblue"), pt.size = 0.01)
cowplot::ggsave2(p, file = "../output/umap_relSeC2.png", width = 6.5, height = 5.8)
