library(data.table)
library(dplyr)
library(BuenColors)

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_sgRNA_260901.rds")
DefaultAssay(so_p) <- "tRNA_pp_valid"
dim(so_p) #479 89088
so_p[["RNA"]] <- JoinLayers(so_p[["RNA"]])
so_p$target <- str_remove(so_p$identity_NTC, "_[AB]$")

so_p <- JoinLayers(so_p, assay = "RNA")
lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "RNA",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_mRNA
#write.csv(all_diff_out_mRNA, "../output/all_diff_out_mRNA_260901.csv")
all_diff_out <- read.csv("../output/all_diff_out_mRNA_260901.csv")
possible_genes <- rownames(so_p)
possible_gRNAs <- unique(c(so_p$identity_NTC))
candidate_target_vector <- stringr::str_split_fixed(possible_gRNAs, "_", 2)[,1]
unique(candidate_target_vector[!(candidate_target_vector %in% possible_genes)])

# find missing genes
"MARS1" %in% possible_genes
"HARS1" %in% possible_genes
"KARS1" %in% possible_genes
"NARS1" %in% possible_genes
"VARS1" %in% possible_genes
"IARS1" %in% possible_genes
"WARS1" %in% possible_genes
"AARS1" %in% possible_genes
"QARS1" %in% possible_genes
"TARS1" %in% possible_genes
"RARS1" %in% possible_genes

# manually update mismatched guide/gene names
candidate_target_vector <- case_when(
  candidate_target_vector == "MARS" ~ "MARS1",
  candidate_target_vector == "HARS" ~ "HARS1",
  candidate_target_vector == "KARS" ~ "KARS1",
  candidate_target_vector == "NARS" ~ "NARS1",
  candidate_target_vector == "VARS" ~ "VARS1",
  candidate_target_vector == "IARS" ~ "IARS1",
  candidate_target_vector == "WARS" ~ "WARS1",
  candidate_target_vector == "AARS" ~ "AARS1",
  candidate_target_vector == "QARS" ~ "QARS1",
  candidate_target_vector == "TARS" ~ "TARS1",
  candidate_target_vector == "RARS" ~ "RARS1",
  candidate_target_vector %in% c("NTC", "RPL17", "RPL14", "EIF3F", "RPL11",
                                 "GLTSCR1", "RPL13", "TIMM23B", "PNPT1","RPL10") ~ "no_flex_probe",
  TRUE ~ candidate_target_vector
)
names(candidate_target_vector) <- possible_gRNAs

# annotate on-target impact
all_diff_out$target <- candidate_target_vector[all_diff_out$gRNA]
all_diff_out$log2FoldChange <- round(all_diff_out$avg_log2FC, 2)
mean(all_diff_out %>% filter(target == probe) %>% 
       pull(log2FoldChange) < 0)

update_list <- function(probe){
  case_when(
    probe == "AARS" ~ "AARS1",
    probe == "MARS" ~ "MARS1",
    probe == "HARS" ~ "HARS1",
    probe == "KARS" ~ "KARS1",
    probe == "NARS" ~ "NARS1",
    probe == "VARS" ~ "VARS1",
    probe == "IARS" ~ "IARS1",
    probe == "WARS" ~ "WARS1",
    probe == "AARS" ~ "AARS1",
    probe == "QARS" ~ "QARS1",
    probe == "TARS" ~ "TARS1",
    probe == "RARS" ~ "RARS1",
    TRUE ~ probe
  )
}
candidate_genes <- unique(all_diff_out$target)

all_diff_out %>%
  mutate(probe = update_list(probe)) %>%
  mutate(target_gene = update_list(target)) %>%
  filter(probe %in% candidate_genes) %>%
  filter(target_gene != "no_flex_probe") %>%
  mutate(log2FoldChange_plot = ifelse(p_val_adj <= 0.1,
      log2FoldChange, 0)) %>%
  mutate(guideReplicate = stringr::str_extract(gRNA, "[^_]+$")) %>% 
  ggplot(aes(x = target_gene, y = probe, fill = log2FoldChange_plot)) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("solar_flare"), limits = c(-2.75, 2.75), oob = scales::squish) +
  pretty_plot(fontsize = 8) +  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_discrete(limits=rev) +
  labs(x = "gRNA target", y = "mRNA expression") +
  facet_wrap(~guideReplicate)
cowplot::ggsave2(p, file = "../output/heatmap_guideeffect_AB_pvalfiltered3.pdf", width = 15.5, height = 8) # p val adj ≤ 0.1

# get fraction of guides resulting in ≤-0.5 log2fc
guide_summary_rep <- all_diff_out %>%
  mutate(
    probe = update_list(probe),
    target_gene = update_list(target),
    guideReplicate = stringr::str_extract(gRNA, "[^_]+$")
  ) %>%
  filter(probe %in% candidate_genes) %>%
  filter(probe == target_gene) %>%
  distinct(gRNA, .keep_all = TRUE) %>%
  dplyr::summarize(
    n_guides = n(),
    n_success = sum(log2FoldChange <= -0.5, na.rm = TRUE),
    pct_success = 100 * n_success / n_guides,
    .groups = "drop"
  )
guide_summary_rep
