library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)

filter <- dplyr::filter
select <- dplyr::select

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")

# compare AAdepl and Torin1 treatment
so_p2_AT <- subset(so_p2, subset = sample %in% c("AAdepr1", "AAdepr2", "Torin1"))
so_p2_AT$sample2 <- substr(so_p2_AT$sample, 1, 3)
unique(so_p2_AT$sample2)

so_p2_AT_ntc <- subset(so_p2_AT, subset = identity_NTC == "NTC")

# look at mRNA
Idents(so_p2_AT_ntc) <- "sample2"
so_p2_AT_ntc <- JoinLayers(so_p2_AT_ntc)
diff_atm <- FindMarkers(so_p2_AT_ntc, ident.1 = "Tor", ident.2 = "AAd", assay = "RNA")
diff_atm$gene <- rownames(diff_atm)
log2FC_cutoff <- 1
pval_cutoff <- 0.05
diff_atm <- diff_atm %>%
  mutate(
    significance = case_when(
      abs(avg_log2FC) >= log2FC_cutoff & p_val_adj < pval_cutoff ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    SLC_gene = str_starts(gene, "SLC"),
    MT_gene = str_starts(gene, "MT-"),
    ARS_gene = str_detect(gene, "ARS(1|2)$")
  )
diff_atm <- diff_atm %>%
  mutate(
    gene_group = case_when(
      str_starts(gene, "SLC") ~ "SLC",
      str_detect(gene, "ARS1") ~ "cyto-ARS",
      str_detect(gene, "ARS2") ~ "mito-ARS",
      str_starts(gene, "MT-") ~ "MT",
      TRUE ~ "Other"
    )
  )
diff_atm <- diff_atm %>%
  mutate(gene_group = factor(gene_group, levels = c("Other", "SLC", "MT", "cyto-ARS", "mito-ARS"))) %>%
  arrange(gene_group)
p <- ggplot(diff_atm, aes(x = avg_log2FC, y = -1*log10(p_val_adj), color = gene_group)) +
  geom_point(size=2.4) +
  scale_color_manual(values = c(
    "Other" = "lightgrey",
    "SLC"   = "darkblue",
    "cyto-ARS"   = "skyblue", # darkgreen
    "mito-ARS"   = "#de9aaf", # de9aad #e1bdc0
    "MT"    = "violetred4"
  )) +
  theme_minimal() +
  labs(title = "DEGs in Torin-1 treatment vs. AA depletion") +
  geom_text_repel(
    data = filter(diff_atm, gene_group %in% c("SLC")), aes(label = gene),
    size = 3, max.overlaps = 15, box.padding = 0.3, point.padding = 0.4) + 
  pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/volcano_Torin1_v_AAdepl_geneGroups_6.pdf", width = 4.4, height = 3.8)
