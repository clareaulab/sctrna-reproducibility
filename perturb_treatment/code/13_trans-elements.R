library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(tibble)
library(BuenColors)
library(ggrepel)
library(presto)
library(pheatmap)
library(uwot)
library(ggseqlogo)
filter <- dplyr::filter
select <- dplyr::select

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
so_p2 <- subset(so_p2, subset = target != "TRMT5")
table(so_p2$identity_NTC)

cmp <- read_csv("../output/limma_ControlNTCs-vs-treatmentNTCs.csv")
df_top <- cmp %>% 
  filter(category == "both, same direction" & logFC_aad > 0) %>% 
  #filter(significance == "Significant") %>% 
  pull(trna)
df_toppp <- sub("-.*", "", df_top)
length(df_toppp) #69

# extract sequences
df <- read.csv("../data/tRNA_sequences_Abox-Bbox-custom_quant.csv")
abox_unique <- df %>%
  dplyr::mutate(abox_custom = toupper(abox_custom)) %>%
  dplyr::filter(!is.na(abox_custom)) %>%
  dplyr::group_by(probePairID) %>%
  dplyr::filter(dplyr::n_distinct(abox_custom) == 1) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

fg_seqs <- abox_unique$abox_custom[abox_unique$probePairID %in% df_toppp]
bg_seqs <- abox_unique$abox_custom[!abox_unique$probePairID %in% df_toppp]

# helper functions
count_matrix <- function(seqs, alphabet = c("A", "C", "G", "T")) {
  seqs <- toupper(seqs)
  L <- unique(nchar(seqs))
  if (length(L) != 1) stop("Sequences must all be the same length (align first)")
  sapply(seq_len(L), function(i) table(factor(substr(seqs, i, i), levels = alphabet)))
}

to_freq <- function(cnt, pseudo = 0.5) {
  sweep(cnt + pseudo, 2, colSums(cnt + pseudo), "/")
}

score_pssm <- function(seqs, pssm) {
  seqs <- toupper(seqs)
  L <- ncol(pssm)
  vapply(seqs, function(s) {
    if (is.na(s) || nchar(s) != L) return(NA_real_)
    bases <- strsplit(s, "")[[1]]
    idx <- match(bases, rownames(pssm))
    if (anyNA(idx)) return(NA_real_)
    sum(pssm[cbind(idx, seq_len(L))])
  }, numeric(1), USE.NAMES = FALSE)
}

cnt_fg <- count_matrix(fg_seqs)
cnt_bg <- count_matrix(bg_seqs)
p_fg <- to_freq(cnt_fg)
p_bg <- to_freq(cnt_bg)
pssm <- log2(p_fg / p_bg)
kl <- p_fg * log2(p_fg / p_bg)

scores2 <- abox_unique %>%
  dplyr::transmute(probePairID,
                   abox_score = score_pssm(abox_custom, pssm))

# select factors
so_trans <- subset(so_p2, subset = target %in% c("GTF3C1", "GTF3C6", "NTC", "POLR3C"))
so_trans$target <- str_remove(so_trans$identity_NTC, "_[AB]$")
so_trans <- subset(so_trans, subset = target != "TRMT5")
so_trans$sample2 <- substr(so_trans$sample, 1, 3)
so_trans$target_sample <- paste0(so_trans$identity_NTC, "_", so_trans$sample2)
so_trans$target_sample2 <- paste0(so_trans$target, "_", so_trans$sample2)
table(so_trans$target_sample2)
possible_grnas <- unique(so_trans$target_sample2)[!grepl("^NTC_", unique(so_trans$target_sample2))]
lapply(possible_grnas, function(one_trna_grna) {
  condition <- sub(".*_(Con|Tor|AAd)$", "\\1", one_trna_grna)
  ntc <- paste0("NTC_", condition)
  
  fmout <- FindMarkers(so_trans, ident.1 = one_trna_grna, ident.2 = ntc,
                       assay = "tRNA_pp_valid", group.by = "target_sample2", logfc.threshold = 0, min.pct = 0)
  
  fmout$probe <- rownames(fmout)
  fmout$gRNA <- one_trna_grna
  fmout$NTC <- ntc
  fmout
}) %>% 
  rbindlist() %>% data.frame() -> all_diff_out
#write.csv(all_diff_out, "../output/all_diff_out_tRNApp_261021_trans.csv")
all_diff_out <- read_csv("../output/all_diff_out_tRNApp_261021_trans.csv")

log2FC_cutoff <- 0
pval_cutoff <- 0.05
all_diff_out <- all_diff_out %>%
  mutate(
    significance = case_when(
      abs(avg_log2FC) >= log2FC_cutoff & p_val_adj < pval_cutoff ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )
deg_trans <- all_diff_out %>%
  dplyr::mutate(probePairID = sub("-.*", "", probe)) %>%
  dplyr::left_join(scores2, by = "probePairID") %>%
  dplyr::mutate(regulation = dplyr::case_when(
    significance != "Significant" ~ "Not significant",
    avg_log2FC > 0                ~ "Upregulated",
    TRUE                          ~ "Downregulated")) %>%
  dplyr::arrange(regulation != "Not significant", abox_score)
deg_trans <- deg_trans %>%
  mutate(target    = factor(sub("_[^_]+$", "", gRNA), levels = c("POLR3C", "GTF3C1", "GTF3C6")),
         condition = factor(sub(".*_", "", gRNA), levels = c("Con", "AAd", "Tor")))
p <- ggplot(deg_trans %>% filter(pct.1 > 0.05 & pct.2 > 0.05 & p_val < 0.01), aes(abox_score, avg_log2FC)) +
  geom_point(size = 1) + facet_grid(condition ~ target) +
  geom_smooth(method = "lm", se = TRUE, colour = "violetred4", linewidth = 0.5) +
  stat_cor(method = "spearman", cor.coef.name = "rho",
           label.x.npc = "left", label.y.npc = "top", size = 2.5) +
  labs(x = "A-box motif score (log2 odds)", y = "avg log2FC (targeting guide vs NTC)") +
  pretty_plot() + L_border()
ggsave(p, filename = "../output/corplot_trans-factors_filtered_261021.pdf", width = 5, height = 4, units = "in")
