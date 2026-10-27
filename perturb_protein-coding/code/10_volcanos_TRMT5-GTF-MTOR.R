library(dplyr)
library(Seurat)
library(BuenColors)
library(ggplot2)

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_sgRNA_260901.rds")
so_p$target <- str_remove(so_p$identity_NTC, "_[AB]$")

pp <- read.csv("../../global_functions/small_data/tRNA_probes_g3_design_human_pairs_250105.csv")
pp <- pp %>%
  mutate(G37 = str_sub(seq1, -1, -1) == "C")

# plot TRMT5 guide A
fc_a <- FindMarkers(so_p, ident.1 = "TRMT5_B", ident.2 = "NTC", 
                    assay = "tRNA_pp_valid", group.by = "identity_NTC",
                    logfc.threshold = 0)
fc_a$probe <- rownames(fc_a)
fc_a <- fc_a %>%
  mutate(probePairID = str_extract(probe, "^PP\\d+"))
fc_a <- fc_a %>%
  left_join(pp %>% select(probePairID, G37), by = "probePairID")
fc_sig <- fc_a %>%
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0 & pct.2 > 0)
p <- ggplot(fc_sig, aes(x = avg_log2FC, y = -1*log(p_val_adj), color = G37)) +
  geom_point() +
  scale_color_manual(values = c("lightgrey", "violetred4")) +
  theme_minimal() + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/volcano_G37_guideA.pdf", width = 3, height = 2.7)

# look at guide B individually
fc_c <- FindMarkers(so_p, ident.1 = "TRMT5_B", ident.2 = "NTC", 
                    assay = "tRNA_pp_valid", group.by = "identity_NTC",
                    logfc.threshold = 0)
fc_c$probe <- rownames(fc_c)
fc_c <- fc_c %>%
  mutate(probePairID = str_extract(probe, "^PP\\d+"))
fc_c <- fc_c %>%
  left_join(pp %>% select(probePairID, G37), by = "probePairID")
fc_sig <- fc_c %>%
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0 & pct.2 > 0)
p <- ggplot(fc_sig, aes(x = avg_log2FC, y = -1*log(p_val_adj), color = G37)) +
  geom_point() +
  scale_color_manual(values = c("lightgrey", "violetred4")) +
  theme_minimal() + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/volcano_G37_guideB_2.pdf", width = 3, height = 2.7)

# function to get colvano for individual guide
get_volcano <- function(guide_AB){
  fc_a <- FindMarkers(so_p, ident.1 = guide_AB, ident.2 = "NTC", 
                      assay = "tRNA_pp_valid", group.by = "identity_NTC",
                      logfc.threshold = 0)
  fc_a$probe <- rownames(fc_a)
  fc_a <- fc_a %>%
    mutate(probePairID = str_extract(probe, "^PP\\d+"))
  fc_a <- fc_a %>%
    left_join(pp %>% select(probePairID, G37), by = "probePairID")
  fc_sig <- fc_a %>%
    filter(p_val_adj < 0.05) %>%
    filter(pct.1 > 0 & pct.2 > 0)
  p <- ggplot(fc_sig, aes(x = avg_log2FC, y = -1 * log(p_val_adj))) +
    geom_point(data = filter(fc_sig, G37 == FALSE), color = "lightgrey", size=1.4) +
    geom_point(data = filter(fc_sig, G37 == TRUE), color = "violetred4", size=1.4) +
    theme_minimal() + pretty_plot() + L_border()
  p
}

# choose additional example
unique(so_p$target)

p <- get_volcano("GTF3C1_A")
#cowplot::ggsave2(p, file = "../output/volcano_GTF3C1_guideA.pdf", width = 2.2, height = 2.7)

p <- get_volcano("GTF3C1_B")
#cowplot::ggsave2(p, file = "../output/volcano_GTF3C1_guideB.pdf", width = 2.2, height = 2.7)

p <- get_volcano("MTOR_A")
#cowplot::ggsave2(p, file = "../output/volcano_MTOR_guideA.pdf", width = 2.2, height = 2.7)

p <- get_volcano("MTOR_B")
#cowplot::ggsave2(p, file = "../output/volcano_MTOR_guideB.pdf", width = 2.2, height = 2.7)
