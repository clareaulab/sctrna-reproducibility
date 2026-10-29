library(dplyr)
library(Seurat)

so_p <- readRDS("../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA.rds")
unique(so_p$identity_NTC)
so_p$target <- str_remove(so_p$identity_NTC, "_guide[123]$")
trna_targets <- unique(so_p$target)[unique(so_p$target) != "NTC"]
grna_targets <- unique(so_p$identity_NTC)[unique(so_p$identity_NTC) != "NTC"]

so_p <- NormalizeData(so_p, assay = "tRNA_pp")

# look at crude proportions
odf <- data.frame(so_p@meta.data, data.matrix((t(so_p@assays$tRNA_aa@counts)/colSums(pilot@assays$tRNA_aa@counts))*100))

# SeC
melt_df <- odf[c("target", rownames(so_p@assays$tRNA_aa@counts))] %>%
  reshape2::melt(id.vars = c("target"))
melt_df$on_target <- str_split_fixed(melt_df[[1]], "-", 3)[,2] == melt_df[,"variable"]
melt_df$target_label <- str_split_fixed(melt_df$target, "-", 3)[, 2]
#melt_df$target_label <- factor(melt_df$target_label, levels = unique(melt_df$target_label))
p <- ggplot(melt_df %>% filter(variable == "SeC"), aes(x = target_label, y = value, color = on_target)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_bw() + coord_cartesian(ylim=c(0,0.8)) + pretty_plot() +
  scale_color_manual(values = c("TRUE" = "violetred4", "FALSE" = "#74ccd4")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#cowplot::ggsave2(p, file = "../output/isodecoder-comp_SeC3.pdf", width = 4.2, height = 1.7, dpi = 400)


# reference for AA and AC mapping
mapping <- readxl::read_xlsx("../data/tRNA_anticodons.xlsx")

# loop
lapply(rownames(so_p@assays$tRNA_ac@counts), function(one_ac){
  print(one_ac)
  target_aa <- mapping %>% filter(RNA_anticodon == one_ac) %>% pull(AA)
  all_acs <- mapping %>% filter(AA == target_aa) %>% pull(RNA_anticodon)
  subset_counts <- data.matrix(so_p@assays$tRNA_ac@counts[all_acs,, drop = FALSE])
  data.frame(cell = colnames(subset_counts), one_ac, value = subset_counts[one_ac,]/colSums(subset_counts)*100,
             AA = target_aa,
             target = so_p@meta.data$target, identity_NTC = so_p@meta.data$identity_NTC
  )
}) %>% rbindlist() %>% data.frame() -> long_prop

long_prop$on_target <- str_split_fixed(long_prop[[5]], "-", 3)[,2] == long_prop[,"AA"]
long_prop$target_label <- str_split_fixed(long_prop$target, "-", 4)[, 3]

# change label to AA
long_prop <- long_prop %>% 
  left_join(mapping %>% select(DNA_anticodon, AA), by = c("target_label" = "DNA_anticodon"))
long_prop$AA_label <- paste0(long_prop$AA.y, "_", long_prop$target_label)

# plot for Lys
p <- ggplot(long_prop %>% filter(AA.x == "Lys"), aes(x = AA.y, y = value, color = on_target)) + 
  geom_boxplot(outlier.shape = NA) + facet_wrap(~one_ac, scales = "free_y", ncol = 1) +  pretty_plot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("TRUE" = "violetred4", "FALSE" = "#74ccd4"))
cowplot::ggsave2(p, file = "../output/isodecoder-comp_Lys_label3.pdf", width = 4.2, height = 4, dpi = 400)

# plot for Ile
p <- ggplot(long_prop %>% filter(AA.x == "Ile"), aes(x = AA.y, y = value, color = on_target)) + 
  geom_boxplot(outlier.shape = NA) + facet_wrap(~one_ac, scales = "free_y", ncol = 1) +  pretty_plot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("TRUE" = "#fbb040", "FALSE" = "#74ccd4"))
cowplot::ggsave2(p, file = "../output/isodecoder-comp_Ile_label3.pdf", width = 4.2, height = 5.5, dpi = 400)

# modify labels to only anticodon
long_prop <- long_prop %>% dplyr::mutate(ac_label = str_split_i(target, "-", 3))
long_prop <- long_prop %>%
  mutate(ac_label = if_else(is.na(ac_label), "NTC", ac_label))

# plot for Val -- supplement
p <- ggplot(long_prop %>% dplyr::filter(AA.x == "Val"), aes(x = AA.y, y = value, color = on_target)) + 
  geom_boxplot(outlier.shape = NA) + facet_wrap(~one_ac, scales = "free_y", ncol = 4) +  pretty_plot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("TRUE" = "#fbb040", "FALSE" = "#74ccd4"))
cowplot::ggsave2(p, file = "../output/isodecoder-comp_Val_label3.pdf", width = 20, height = 3, dpi = 400)

# plot for Leu -- supplement
p <- ggplot(long_prop %>% filter(AA.x == "Leu"), aes(x = AA.y, y = value, color = on_target)) + 
  geom_boxplot(outlier.shape = NA) + facet_wrap(~one_ac, scales = "free_y", ncol = 4) +  pretty_plot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("TRUE" = "#fbb040", "FALSE" = "#74ccd4"))
cowplot::ggsave2(p, file = "../output/isodecoder-comp_Leu_label4.pdf", width = 20, height = 5, dpi = 400)
