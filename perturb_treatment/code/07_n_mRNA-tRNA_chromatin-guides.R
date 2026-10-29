library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
library(pheatmap)
filter <- dplyr::filter
select <- dplyr::select

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
so_p2 <- subset(so_p2, subset = target != "TRMT5")
table(so_p2$identity_NTC)
so_p2_ctrl <- subset(so_p2, subset = sample == "Control")
so_p2_AAd <- subset(so_p2, subset = sample %in% c("AAdepr1", "AAdepr2"))
so_p2_tor <- subset(so_p2, subset = sample == "Torin1")

possible_grnas <- unique(so_p2_ctrl$identity_NTC)[!grepl("NTC", unique(so_p2_ctrl$identity_NTC))]
possible_grnas <- setdiff(possible_grnas, c("TRMT5_A", "TRMT5_B"))

lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p2_ctrl, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "tRNA_pp_valid",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_ctrl
#write.csv(all_diff_out_ctrl, "../output/all_diff_out_tRNApp_261008_ctrl.csv")
all_diff_out_ctrl <- read.csv("../output/all_diff_out_tRNApp_261008_ctrl.csv")

grna_t <- all_diff_out_ctrl %>% dplyr::group_by(gRNA) %>%
  dplyr::summarize(n_hits = sum(p_val_adj < 0.05 #& abs(avg_log2FC) > 0.5
                                ), count = n()) %>% arrange(desc(n_hits))

so_p2_ctrl <- JoinLayers(so_p2_ctrl, assay = "RNA")
lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p2_ctrl, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "RNA",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_mRNA_ctrl
#write.csv(all_diff_out_mRNA_ctrl, "../output/all_diff_out_mRNA_261008_ctrl.csv")
all_diff_out_mRNA_ctrl <- read.csv("../output/all_diff_out_mRNA_261008_ctrl.csv")

grna_m <- all_diff_out_mRNA_ctrl %>% group_by(gRNA) %>%
  #dplyr::summarize(n_hits = sum(p_val_adj < 0.01), count = n()) %>% arrange(desc(n_hits))
  dplyr::summarize(n_hits = sum(p_val_adj < 0.05 #& abs(avg_log2FC) > 1
                                ), count = n()) %>% arrange(desc(n_hits))

n_mRNA_df <- all_diff_out_mRNA_ctrl %>% dplyr::group_by(probe) %>% dplyr::summarize(is_diffential = sum(p_val_adj < 0.05))
often_differential_mRNA <- n_mRNA_df %>% #dplyr::filter(is_diffential >= 5) %>% 
  dplyr::pull(probe)

grna_merge <- merge(grna_m, grna_t, by = "gRNA")
grna_merge <- grna_merge %>% 
  dplyr::rename("mRNA"=n_hits.x)
grna_merge <- grna_merge %>% 
  dplyr::rename("tRNA"=n_hits.y)
grna_merge <- grna_merge %>% 
  dplyr::select(gRNA,mRNA,tRNA)


# AAd cells
lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p2_AAd, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "tRNA_pp_valid",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_aad
#write.csv(all_diff_out_aad, "../output/all_diff_out_tRNApp_261008_aad.csv")
all_diff_out_aad <- read.csv("../output/all_diff_out_tRNApp_261008_aad.csv")

grna_t <- all_diff_out_aad %>% dplyr::group_by(gRNA) %>%
  dplyr::summarize(n_hits = sum(p_val_adj < 0.05 #& abs(avg_log2FC) > 0.5
  ), count = n()) %>% arrange(desc(n_hits))

so_p2_AAd <- JoinLayers(so_p2_AAd, assay = "RNA")
lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p2_AAd, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "RNA",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_mRNA_aad
write.csv(all_diff_out_mRNA_aad, "../output/all_diff_out_mRNA_261008_aad.csv")
all_diff_out_mRNA_aad <- read.csv("../output/all_diff_out_mRNA_261008_aad.csv")

grna_m <- all_diff_out_mRNA_aad %>% group_by(gRNA) %>%
  #dplyr::summarize(n_hits = sum(p_val_adj < 0.01), count = n()) %>% arrange(desc(n_hits))
  dplyr::summarize(n_hits = sum(p_val_adj < 0.05 #& abs(avg_log2FC) > 1
  ), count = n()) %>% arrange(desc(n_hits))

n_mRNA_df <- all_diff_out_mRNA_aad %>% dplyr::group_by(probe) %>% dplyr::summarize(is_diffential = sum(p_val_adj < 0.05))
often_differential_mRNA <- n_mRNA_df %>% #dplyr::filter(is_diffential >= 5) %>% 
  dplyr::pull(probe)

grna_merge_aad <- merge(grna_m, grna_t, by = "gRNA")
grna_merge_aad <- grna_merge_aad %>% 
  dplyr::rename("mRNA"=n_hits.x)
grna_merge_aad <- grna_merge_aad %>% 
  dplyr::rename("tRNA"=n_hits.y)
grna_merge_aad <- grna_merge_aad %>% 
  dplyr::select(gRNA,mRNA,tRNA)


# Torin1 cells
lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p2_tor, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "tRNA_pp_valid",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_tor
write.csv(all_diff_out_tor, "../output/all_diff_out_tRNApp_261008_tor.csv")
all_diff_out_tor <- read.csv("../output/all_diff_out_tRNApp_261008_tor.csv")

grna_t <- all_diff_out_tor %>% dplyr::group_by(gRNA) %>%
  dplyr::summarize(n_hits = sum(p_val_adj < 0.05 #& abs(avg_log2FC) > 0.5
  ), count = n()) %>% arrange(desc(n_hits))

so_p2_tor <- JoinLayers(so_p2_tor, assay = "RNA")
lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p2_tor, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "RNA",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_mRNA_tor
#write.csv(all_diff_out_mRNA_tor, "../output/all_diff_out_mRNA_261008_aad.csv")
all_diff_out_mRNA_tor <- read.csv("../output/all_diff_out_mRNA_261008_aad.csv")

grna_m <- all_diff_out_mRNA_tor %>% group_by(gRNA) %>%
  #dplyr::summarize(n_hits = sum(p_val_adj < 0.01), count = n()) %>% arrange(desc(n_hits))
  dplyr::summarize(n_hits = sum(p_val_adj < 0.05 #& abs(avg_log2FC) > 1
  ), count = n()) %>% arrange(desc(n_hits))

n_mRNA_df <- all_diff_out_mRNA_tor %>% dplyr::group_by(probe) %>% dplyr::summarize(is_diffential = sum(p_val_adj < 0.05))
often_differential_mRNA <- n_mRNA_df %>% #dplyr::filter(is_diffential >= 5) %>% 
  dplyr::pull(probe)

grna_merge_tor <- merge(grna_m, grna_t, by = "gRNA")
grna_merge_tor <- grna_merge_tor %>% 
  dplyr::rename("mRNA"=n_hits.x)
grna_merge_tor <- grna_merge_tor %>% 
  dplyr::rename("tRNA"=n_hits.y)
grna_merge_tor <- grna_merge_tor %>% 
  dplyr::select(gRNA,mRNA,tRNA)


# combine all threee
grna_merge_tor$condition <- "Torin1"
grna_merge_aad$condition <- "AAd"
grna_merge$condition <- "Ctrl"

grna_all <- rbind(grna_merge, grna_merge_tor, grna_merge_aad)

ggplot(grna_all, aes(x = mRNA, y = tRNA, color = condition)) +
  labs(x = "# of differentially expressed genes - mRNA", y = "# of differentially expressed genes - tRNA") +
  geom_point() + pretty_plot() +
  #geom_text_repel(data = grna_merge_tor, #data = filter(all_diff_out, significance == "Significant"),
  #  aes(label = gRNA), size = 3,
  #  max.overlaps = 5, box.padding = 0.2, point.padding = 0.1) + 
  pretty_plot() + L_border()


guides_trnaps <- read_excel("../data/opool_Twist_Perturb-tRNA_manual_250707 copy.xlsx")
guides_chrom <- guides_trnaps[193:260,1]
guides_chrom_vec <- as.character(guides_chrom[[1]])
guides_chrom_vec <- gsub("-","_",guides_chrom_vec)
guides_chrom_vec2 <- setdiff(guides_chrom_vec, "NTC")
guides_chrom_vec3 <- gsub("_[A/B]", "", guides_chrom_vec2)
guides_chrom_vec4 <- setdiff(guides_chrom_vec3, "BRF2")
guides_chrom_vec5 <- setdiff(guides_chrom_vec2, c("BRF2_A", "BRF2_B"))

ggplot(grna_all %>% filter(gRNA %in% guides_chrom_vec5), aes(x = mRNA, y = tRNA, color = condition)) +
  labs(x = "# of differentially expressed genes - mRNA", y = "# of differentially expressed genes - tRNA") +
  geom_point() + pretty_plot() +
  #geom_text_repel(data = grna_merge_tor, #data = filter(all_diff_out, significance == "Significant"),
  #  aes(label = gRNA), size = 3,
  #  max.overlaps = 5, box.padding = 0.2, point.padding = 0.1) + 
  pretty_plot() + L_border() #+ xlim(0,700) + ylim(0,120)


# box plot
ggplot(grna_all %>% filter(gRNA %in% guides_chrom_vec5), aes(x = mRNA, y = tRNA, color = condition)) +
  labs(x = "# of differentially expressed genes - mRNA", y = "# of differentially expressed genes - tRNA") +
  geom_point() + pretty_plot() +
  pretty_plot() + L_border()

grna_all %>%
  mutate(in_vec5 = ifelse(gRNA %in% guides_chrom_vec5, "chromatin guide", "all other")) %>%
  ggplot(aes(x = condition, y = tRNA, fill = in_vec5)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
               outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.8, alpha = 0.6) +
  labs(x = "Condition",
       y = "# of differentially expressed genes - tRNA",
       fill = NULL) +
  pretty_plot() + L_border()


# tRNA
grna_all %>%
  mutate(in_vec5 = ifelse(gRNA %in% guides_chrom_vec5, "chromatin guide", "all other"),
         condition = factor(condition, levels = c("Ctrl", "AAd", "Torin1"))) %>%
  ggplot(aes(x = condition, y = tRNA, fill = in_vec5)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
               outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.8, alpha = 0.6) +
  stat_compare_means(aes(group = in_vec5,
                       label = sprintf("p = %.2f", after_stat(p))),
                   method = "wilcox.test",
                   label.y.npc = 0.98, size = 3)
  scale_fill_manual(values = c("all other" = "grey", "chromatin guide" = "violetred4")) +
  labs(x = "Condition",
       y = "# of differentially expressed genes - tRNA",
       fill = NULL) +
  pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/boxplot_n_chrom_tRNA.pdf", width = 5, height = 3.5)

  
# get exact p-values
df <- grna_all %>%
    mutate(in_vec5 = ifelse(gRNA %in% guides_chrom_vec5, "chromatin guide", "all other"),
           condition = factor(condition, levels = c("Ctrl", "AAd", "Torin1")))
compare_means(tRNA ~ in_vec5, data = df, group.by = "condition",
                method = "wilcox.test", p.adjust.method = "BH")

# mRNA
p <- grna_all %>%
  mutate(in_vec5 = ifelse(gRNA %in% guides_chrom_vec5, "chromatin guide", "all other"),
         condition = factor(condition, levels = c("Ctrl", "AAd", "Torin1"))) %>%
  ggplot(aes(x = condition, y = mRNA, fill = in_vec5)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
               outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.8, alpha = 0.6) +
  stat_compare_means(aes(group = in_vec5),
                     method = "wilcox.test", label = "p.format",
                     label.y.npc = 0.98, size = 3) +
  scale_fill_manual(values = c("all other" = "grey", "chromatin guide" = "violetred4")) +
  labs(x = "Condition",
       y = "# of differentially expressed genes - mRNA",
       fill = NULL) +
  pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/boxplot_n_chrom_mRNA.pdf", width = 5, height = 3.5)

# get exact p-values
df <- grna_all %>%
  mutate(in_vec5 = ifelse(gRNA %in% guides_chrom_vec5, "chromatin guide", "all other"),
         condition = factor(condition, levels = c("Ctrl", "AAd", "Torin1")))
compare_means(mRNA ~ in_vec5, data = df, group.by = "condition",
              method = "wilcox.test", p.adjust.method = "BH")
