library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
select <- dplyr::select
filter <- dplyr::filter

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
so_p2$sample2 <- substr(so_p2$sample, 1, 3)
so_p2_ctrl <- subset(so_p2, subset = sample2 == "Con")

# define ISR and non-ISR control guides
guides_isr <- c("MTOR", "RPL14", "POLR1B", "MARS", "EIF2B4", "EIF2B3", "HARS", "AARS", 
                "TARS", "POLR3D", "QARS", "HSPA9", "ATF5", "BDP1")

Idents(so_p2_ctrl) <- so_p2_ctrl$target
DefaultAssay(so_p2_ctrl) <- "RNA"
deg_isr <- FindMarkers(
  so_p2_ctrl,
  ident.1 = guides_isr,
  ident.2 = "NTC",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
deg_isr$gene <- rownames(deg_isr)
log2FC_cutoff <- 1
pval_cutoff <- 0.05
deg_isr <- deg_isr %>%
  mutate(
    significance = case_when(
      abs(avg_log2FC) >= log2FC_cutoff & p_val_adj < pval_cutoff ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )
ggplot(deg_isr, aes(x = avg_log2FC, y = -1*log10(p_val_adj), color = significance)) +
  geom_point(size=3.5) +
  scale_color_manual(values = c("lightgrey", "darkblue")) +
  theme_minimal() +
  labs(title = "Differential mRNAs") +
  #scale_y_continuous(expand = c(0,0))+
  geom_text_repel(
    data = filter(deg_isr, significance == "Significant"), 
    aes(label = gene),
    size = 5.5,
    max.overlaps = 15,
    box.padding = 0.3,
    point.padding = 0.2
  ) + pretty_plot() + L_border()

# label ISR genes
deg_isr <- deg_isr %>%
  dplyr::mutate(
    highlight = dplyr::case_when(
      gene %in% c("NUPR1", "G0S2", "ATF3", "ATF5", "DDIT3", "DDIT4") ~ "highlight",
      significance == "Significant" ~ "significant",
      TRUE ~ "other"
    )
  )

ggplot(deg_isr, aes(x = avg_log2FC, y = -log10(p_val_adj), color = highlight)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c(
    "highlight" = "violetred4",
    "significant" = "#939bce",
    "other" = "lightgrey"
  )) +
  theme_minimal() +
  labs(title = "Differential mRNAs") +
  geom_text_repel(
    data = dplyr::filter(deg_isr, gene %in% c("NUPR1", "G0S2", "ATF3", "ATF5", "DDIT3", "DDIT4")),
    aes(label = gene),
    size = 5.5,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) + pretty_plot() + L_border()

genes_to_label <- c("NUPR1", "G0S2", "ATF3", "ATF5", "DDIT3", "DDIT4")
p <- ggplot(deg_isr, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = dplyr::filter(deg_isr, !(gene %in% genes_to_label)),
    aes(color = highlight), size = 3.5) +
  geom_point(
    data = dplyr::filter(deg_isr, gene %in% genes_to_label), color = "violetred4", size = 3.8) + 
  scale_color_manual(values = c("significant" = "#939bce", "other" = "lightgrey")) +
  geom_text_repel(
    data = dplyr::filter(deg_isr, gene %in% genes_to_label),
    aes(label = gene),
    size = 5.5,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/volcano_control_ISR-vs-NTC.pdf", width = 4, height = 4.5, dpi = 400)
