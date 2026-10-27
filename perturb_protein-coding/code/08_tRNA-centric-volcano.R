library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
filter <- dplyr::filter
select <- dplyr::select

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_sgRNA_260901.rds")
dim(so_p)

so_p$target <- str_remove(so_p$identity_NTC, "_[AB]$")
Idents(so_p) <- so_p$target
table(Idents(so_p))
guides <- levels(Idents(so_p))
guides_noNTC <- guides[guides != "NTC"]

# Glu-GAA
gaa_de <- map_dfr(guides_noNTC, function(g) {
  res <- FindMarkers(so_p, ident.1 = g, ident.2 = "NTC",
    assay = "tRNA_ac", features = "GAA",
    test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
  res %>% mutate(guide = g)
})
volcano_df <- gaa_de %>%
  mutate(neg_log10_p = -log10(p_val_adj + 1e-300),
    significant = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
top_guides <- volcano_df %>% arrange(p_val_adj) %>%
  filter(abs(avg_log2FC) > 0.1) %>% slice_head(n = 21)
top_guides2 <- c("TRMT5", "ALKBH8", "ELP3", "ELP4", "ELP5", "BRF2", #"CTU1", "CTU2",
                 "GTF3C1", "POLR1D", "POLR3C", "GTF3C6", "MTOR", "RPTOR", "POLR3D")
top_guides_df <- volcano_df %>% filter(guide %in% top_guides2)
p <- ggplot(volcano_df, aes(avg_log2FC, neg_log10_p)) +
  geom_point(color = "grey70", size=1.4) +
  geom_point(data = top_guides_df, color = "violetred4", size = 1.5) +
  geom_text_repel(data = top_guides_df, aes(label = guide), size = 1.8) +
  theme_classic() + pretty_plot(fontsize=5) + L_border()
cowplot::ggsave2(p, file = "../output/volcano_Glu-GAA_guides_noNTC_260904.pdf", width = 2, height = 2)


# Arg AGA
guides_noNTC <- guides[guides != "NTC"]
Idents(so_p) <- "target"
aga_de <- map_dfr(guides_noNTC, function(g) {
  res <- FindMarkers(so_p, ident.1 = g, ident.2 = "NTC",
    assay = "tRNA_ac", features = "AGA",
    test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
  res %>% mutate(guide = g)
})
volcano_df <- aga_de %>%
  mutate(neg_log10_p = -log10(p_val_adj + 1e-300),
    significant = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
top_guides <- volcano_df %>% arrange(p_val_adj) %>%
  filter(abs(avg_log2FC) > 0.1) %>% slice_head(n = 21)
top_guides2 <- c("TRMT5", "ALKBH8", "ELP3", "ELP4", "ELP5", "BRF2", #"CTU1", "CTU2",
                 "GTF3C1", "POLR1D", "POLR3C", "GTF3C6", "MTOR", "RPTOR", "POLR3D")
top_guides_df <- volcano_df %>% filter(guide %in% top_guides2)
p <- ggplot(volcano_df, aes(avg_log2FC, neg_log10_p)) +
  geom_point(color = "grey70", size=1.4) +
  geom_point(data = top_guides_df, color = "violetred4", size = 1.5) + 
  geom_text_repel(data = top_guides_df, aes(label = guide), size = 1.8) +
  theme_classic() + pretty_plot(fontsize=5) + L_border()
cowplot::ggsave2(p, file = "../output/volcano_Arg-AGA_guides_260904.pdf", width = 2, height = 2)
