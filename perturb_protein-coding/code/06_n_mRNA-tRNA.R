library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)

filter <- dplyr::filter
select <- dplyr::select

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_sgRNA_260901.rds")

possible_grnas <- unique(so_p$identity_NTC)[!grepl("NTC", unique(so_p$identity_NTC))]

lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "tRNA_pp_valid",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out
#write.csv(all_diff_out, "../output/all_diff_out_tRNApp_260901.csv")
all_diff_out <- read.csv("../output/all_diff_out_tRNApp_260901.csv")

grna_t <- all_diff_out %>% 
  dplyr::group_by(gRNA) %>%
  dplyr::summarize(n_hits = sum(p_val_adj < 0.01), count = n()) %>% arrange(desc(n_hits))

so_p <- JoinLayers(so_p, assay = "RNA")
lapply(possible_grnas, function(one_trna_grna){
  fmout <- FindMarkers(so_p, ident.1 = one_trna_grna, ident.2 = "NTC", assay = "RNA",group.by = "identity_NTC", 
                       logfc.threshold = 0, min.pct = 0)
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout
}) %>% rbindlist() %>% data.frame() -> all_diff_out_mRNA
#write.csv(all_diff_out_mRNA, "../output/all_diff_out_mRNA_260901.csv")
all_diff_out_mRNA <- read.csv("../output/all_diff_out_mRNA_260901.csv")

grna_m <- all_diff_out_mRNA %>% 
  group_by(gRNA) %>%
  dplyr::summarize(n_hits = sum(p_val_adj < 0.01), count = n()) %>% arrange(desc(n_hits))

n_mRNA_df <- all_diff_out_mRNA %>% dplyr::group_by(probe) %>% dplyr::summarize(is_diffential = sum(p_val_adj < 0.01))
often_differential_mRNA <- n_mRNA_df %>% dplyr::filter(is_diffential >= 5) %>% dplyr::pull(probe)

grna_merge <- merge(grna_m, grna_t, by = "gRNA")
grna_merge <- grna_merge %>% 
  dplyr::rename("mRNA"=n_hits.x)
grna_merge <- grna_merge %>% 
  dplyr::rename("tRNA"=n_hits.y)
grna_merge <- grna_merge %>% 
  dplyr::select(gRNA,mRNA,tRNA)

ggplot(grna_merge, aes(x = mRNA, y = tRNA, label = gRNA)) +
  labs(x = "# of differentially expressed genes - mRNA", y = "# of differentially expressed genes - tRNA") +
  geom_point() + pretty_plot() +
  geom_text_repel(data = grna_merge, aes(label = gRNA),
    size = 3, max.overlaps = 5, box.padding = 0.2, point.padding = 0.1) + 
  pretty_plot() + L_border()


# highlight guides belonging to pol. III machinery or chromatin remodellers
unique(grna_merge$gRNA)
guides_trnaps <- read_excel("../data/opool_Twist_Perturb-tRNA_manual_250707 copy.xlsx")
guides_chrom <- guides_trnaps[193:260,1]
guides_chrom_vec <- as.character(guides_chrom[[1]])
guides_chrom_vec <- gsub("-","_",guides_chrom_vec)
guides_chrom_vec2 <- setdiff(guides_chrom_vec, "NTC")
guides_chrom_vec3 <- gsub("_[A/B]", "", guides_chrom_vec2)
guides_chrom_vec4 <- setdiff(guides_chrom_vec3, "BRF2")

gene_sets <- list(
  "Pol III" = c("BDP1","BRF2","GTF3C1","GTF3C6","POLR3C","POLR3D"),
  "guides_chrom" = guides_chrom_vec4
)
gene_anno <- stack(gene_sets) %>% transmute(gene = values, category = as.character(ind))
gene_anno <- tibble::enframe(gene_sets, name = "category", value = "gene") %>%
  tidyr::unnest(gene) %>%
  mutate(gene = as.character(gene))
highlight_sets <- c("Pol III", "guides_chrom")
plot_df <- grna_merge %>%
  mutate(gene = sub("(_P[12])?_[AB]$", "", gRNA)) %>%
  left_join(gene_anno, by = "gene", relationship = "many-to-many") %>%
  mutate(
    color_group = case_when(
      category == "Pol III" ~ "Pol.III",
      category == "guides_chrom" ~ "guides_chrom",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(color_group)

ggplot(plot_df, aes(x = mRNA, y = tRNA, label = gRNA)) +
  labs(x = "# of differentially expressed genes - mRNA",
       y = "# of differentially expressed genes - tRNA") +
  geom_point(aes(color = color_group)) +
  geom_text(data=plot_df %>% filter(tRNA>75 | mRNA>500),
            aes(color = color_group), size = 3, max.overlaps = 5,
            box.padding = 0.2, point.padding = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("Pol.III" = "violetred4", "guides_chrom" = "darkblue", "Other" = "#b4b4b4"), ##5a87c0
                     labels = c("Pol.III" = "Pol III", "guides_chrom" = "Chromatin remodeller guides", "Other" = "Other"),
                     name = NULL) +
  pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/n_mRNA_tRNA_pol3-chrom_noLabel.pdf", width = 5.6, height = 3.5)
