library(dplyr)
library(Seurat)
library(BuenColors)
library(ggplot2)

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_260901.rds")
so_p[["RNA"]] <- JoinLayers(so_p[["RNA"]])
so_p$target <- str_remove(so_p$identity_NTC, "_[AB]$")

# dotplot for tRNA counts for selection of guides, including mTor and tRNA modification
target_mtor_ntc <- c("MTOR", "RPTOR", "LAMTOR1", "LAMTOR3", "LAMTOR4",
                     "MIOS", "FLCN", "SEH1L", "NTC")
target_selection <- c(target_mtor_ntc, "SIRT7",
                      "CTU1", "CTU2", "ELP3", "ELP4", "ELP5", 
                      "ALKBH8", "URM1", "TRNT1", "TRMT112", "DKC1", "ADAT3",
                      "POLR1A", "POLR1B", "POLR1C", "POLR1D", "POLR3C", "POLR3D", "GTF3C1", "GTF3C6", "BDP1"
)
#so_p$target <- gsub("_[AB]$", "", so_p$identity_NTC)
df_target <- so_p@meta.data %>%
  filter(target %in% target_selection) %>% 
  group_by(target) %>%
  summarize(
    median_tRNA = median(nCount_tRNA_pp_valid),
    mean_tRNA   = mean(nCount_tRNA_pp_valid),
    n_cells     = n(),
    .groups = "drop"
  )
ntc_value_mean <- df_target %>%
  filter(target == "NTC") %>%
  pull(mean_tRNA)
p <- ggplot(df_target, aes(x = reorder(target, mean_tRNA), y = mean_tRNA)) +
  geom_point(size = 2) +
  geom_hline(yintercept = ntc_value_mean, 
             linetype = "dashed", 
             linewidth = 0.6,
             colour = "violetred4") +
  coord_flip() + theme_bw() +
  labs(y = "Mean nCount_tRNA_pp_valid", x = "Target") + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/dotplot_tRNAcount_guides-selection2_collapsed_mean.pdf", width = 3.3, height = 5.5)
