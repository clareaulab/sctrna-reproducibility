library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
library(pheatmap)

select <- dplyr::select
filter <- dplyr::filter

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
so_p2$sample2 <- substr(so_p2$sample, 1, 3)
so_p2$target_sample <- paste0(so_p2$identity_NTC, "_", so_p2$sample2)
so_p2$target_sample2 <- paste0(so_p2$target, "_", so_p2$sample2)
table(so_p2$target_sample2)

# plot log2 FCs
possible_grnas <- unique(so_p2$target_sample2)[!grepl("^NTC_", unique(so_p2$target_sample2))]
lapply(possible_grnas, function(one_trna_grna) {
  condition <- sub(".*_(Con|Tor|AAd)$", "\\1", one_trna_grna)
  ntc <- paste0("NTC_", condition)
  
  fmout <- FindMarkers(so_p2, ident.1 = one_trna_grna, ident.2 = ntc,
                       assay = "tRNA_pp_valid", group.by = "target_sample2", logfc.threshold = 0, min.pct = 0)
  
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout$NTC <- ntc
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out
write.csv(all_diff_out, "../output/all_diff_out_tRNApp_target.csv")
#all_diff_out <- read_csv("../output/all_diff_out_tRNApp_target.csv")

unique(so_p2$target_sample)
guides_want <- c("CTU1_Con", "CTU1_AAd", "CTU1_Tor", "CTU2_Con", "CTU2_AAd", "CTU2_Tor",
                 "ELP3_Con", "ELP3_AAd", "ELP3_Tor", "ELP4_Con", "ELP4_AAd", "ELP4_Tor",
                 "ELP5_Con", "ELP5_AAd", "ELP5_Tor", "URM1_Con", "URM1_AAd", "URM1_Tor",
                 "ALKBH8_Con", "ALKBH8_AAd", "ALKBH8_Tor", "TRMT112_Con", "TRMT112_AAd", "TRMT112_Tor")
all_diff_out2 <- all_diff_out %>% filter(gRNA %in% guides_want)

all_want <- all_diff_out2 %>%
  filter(p_val < 0.05) %>%
  dplyr::arrange(abs(p_val_adj)) %>% 
  filter(pct.2 > 0.1) %>% filter(pct.1 > 0.1) %>%
  pull(probe) %>% unique() %>% head(20)

fc_mat <- all_diff_out2 %>%
  dplyr::filter(probe %in% all_want) %>%
  dplyr::select(probe, gRNA, avg_log2FC) %>%
  pivot_wider(
    names_from = gRNA,
    values_from = avg_log2FC,
    values_fill = 0
  )

fc_mat <- as.data.frame(fc_mat)
rownames(fc_mat) <- fc_mat$probe
fc_mat$probe <- NULL
fc_mat <- as.matrix(fc_mat)
fc_mat[fc_mat < -1] <- -1
my_colors <- jdb_palette("solar_flare", type = "continuous")
my_breaks <- seq(-0.5, 0.5, length.out = length(my_colors))
pheatmap(fc_mat, 
         color = my_colors, 
         breaks = my_breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

gene      <- sub("_.*", "", colnames(fc_mat))   # ALKBH8, CTU1, ...
treatment <- sub(".*_", "", colnames(fc_mat))   # AAd, Con, Tor

ord <- order(gene, factor(treatment, levels = c("Con", "AAd", "Tor")))
fc_mat_sorted <- fc_mat[, ord]

ann <- data.frame(gene = gene[ord], treatment = treatment[ord],
                  row.names = colnames(fc_mat_sorted))
ann_colors <- list(treatment = c(Con = "#fdde8f", AAd = "#C8553D", Tor = "#5AA9C8"))
p <- pheatmap(fc_mat_sorted,
         color = my_colors,
         breaks = my_breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = ann,
         annotation_colors = ann_colors,
         gaps_col = which(diff(as.integer(factor(gene[ord]))) != 0))
#cowplot::ggsave2(p, file = "../output/heatmap_variable-features_modification_split.pdf", width = 5.8, height = 4.4)
