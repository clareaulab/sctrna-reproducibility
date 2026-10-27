library(dplyr)
library(Seurat)
library(ggplot2)
library(BuenColors)

filter <- dplyr::filter
select <- dplyr::select

so_TcE <- readRDS(file = "../../../trna_large_data/output/seurat_object_exhTcell_MM_filtered_0422.rds")
so_TcE <- JoinLayers(so_TcE, assay="RNA")

so_TcE$tRNA_mRNA_ratio <- log1p(so_TcE$nCount_tRNA_pp_valid) /
  log1p(so_TcE$nCount_RNA)

prot_feat <- paste0(Key(so_TcE[["Protein"]]), "MYC")
rna_feat  <- paste0(Key(so_TcE[["RNA"]]), "MYC")

df <- FetchData(so_TcE,
                vars = c(prot_feat, rna_feat, "tRNA_mRNA_ratio", "condition",
                         "nCount_RNA", "nCount_tRNA_pp_valid", "nCount_Protein"),
                layer = "counts")
colnames(df)[1:2] <- c("myc_prot", "myc_rna")

df <- df %>%
  mutate(myc_prot_rna = log2((myc_prot + 1) / (myc_rna + 1)))

# Spearman correlations: pooled and within each condition
cor_all <- cor(df$tRNA_mRNA_ratio, df$myc_prot_rna, method = "spearman")
cor_cond <- df %>%
  group_by(condition) %>%
  summarise(rho = cor(tRNA_mRNA_ratio, myc_prot_rna, method = "spearman"),
            n = n(), .groups = "drop")
cor_all; cor_cond

cond_cols <- c("D2" = "lightblue", "D8_acute" = "#5987C1", "D8_chronic" = "violetred4")

p <- ggplot(df, aes(x = tRNA_mRNA_ratio, y = myc_prot_rna)) +
  geom_point(aes(color = condition), size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
  scale_color_manual(values = cond_cols) +
  scale_x_log10() +
  facet_wrap(~ condition, scales = "free_x") +
  geom_text(data = cor_cond,
            aes(label = paste0("rho = ", round(rho, 2))),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3,
            inherit.aes = FALSE) +
  labs(x = "tRNA / mRNA ratio",
       y = expression(log[2]~"(MYC protein / MYC mRNA)")) +
  pretty_plot() + L_border() +
  theme(legend.position = "none")
#cowplot::ggsave2(p, file = "../output/cor_MTratio_MYC.pdf", width = 5, height = 2.2, dpi = 400)
cowplot::ggsave2(p, file = "../output/cor_MTratio_MYC.png", width = 5, height = 2.2, dpi = 400)


# NFKB2
prot_feat <- paste0(Key(so_TcE[["Protein"]]), "NFKB2")
rna_feat  <- paste0(Key(so_TcE[["RNA"]]), "NFKB2")

df <- FetchData(so_TcE,
                vars = c(prot_feat, rna_feat, "tRNA_mRNA_ratio", "condition",
                         "nCount_RNA", "nCount_tRNA_pp_valid", "nCount_Protein"),
                layer = "counts")
colnames(df)[1:2] <- c("myc_prot", "myc_rna")

df <- df %>%
  mutate(myc_prot_rna = log2((myc_prot + 1) / (myc_rna + 1)))

# Spearman correlations: pooled and within each condition
cor_all <- cor(df$tRNA_mRNA_ratio, df$myc_prot_rna, method = "spearman")
cor_cond <- df %>%
  group_by(condition) %>%
  summarise(rho = cor(tRNA_mRNA_ratio, myc_prot_rna, method = "spearman"),
            n = n(), .groups = "drop")
cor_all; cor_cond

cond_cols <- c("D2" = "lightblue", "D8_acute" = "#5987C1", "D8_chronic" = "violetred4")

p <- ggplot(df, aes(x = tRNA_mRNA_ratio, y = myc_prot_rna)) +
  geom_point(aes(color = condition), size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
  scale_color_manual(values = cond_cols) +
  scale_x_log10() +
  facet_wrap(~ condition, scales = "free_x") +
  geom_text(data = cor_cond,
            aes(label = paste0("rho = ", round(rho, 2))),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3,
            inherit.aes = FALSE) +
  labs(x = "tRNA / mRNA ratio",
       y = expression(log[2]~"(NFKB2 protein / NFKB2 mRNA)")) +
  pretty_plot() + L_border() +
  theme(legend.position = "none")
#cowplot::ggsave2(p, file = "../output/cor_MTratio_NFKB2.pdf", width = 5, height = 2.2, dpi = 400)
cowplot::ggsave2(p, file = "../output/cor_MTratio_NFKB2.png", width = 5, height = 2.2, dpi = 400)
