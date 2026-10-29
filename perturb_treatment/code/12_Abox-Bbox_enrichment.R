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

# extract sequences
df <- read.csv("../data/tRNA_sequences_Abox-Bbox-custom_quant.csv")

# use output from script 08
cmp <- read_csv("../output/limma_ControlNTCs-vs-treatmentNTCs.csv")
df_top <- cmp %>% 
  filter(category == "both, same direction" & logFC_aad > 0) %>% 
  pull(trna)
df_toppp <- sub("-.*", "", df_top)
length(df_toppp) #69

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

# exclude ambiguous PPs
abox_unique <- df %>%
  mutate(abox_custom = toupper(abox_custom)) %>%
  group_by(probePairID) %>%
  filter(n_distinct(abox_custom) == 1) %>%
  slice(1) %>%
  ungroup()

fg_seqs <- abox_unique$abox_custom[abox_unique$probePairID %in% df_toppp]
tested  <- unique(sub("-.*", "", cmp$trna))
bg_seqs <- abox_unique$abox_custom[!abox_unique$probePairID %in% df_toppp]

cnt_fg <- count_matrix(fg_seqs)
cnt_bg <- count_matrix(bg_seqs)
p_fg   <- to_freq(cnt_fg)
p_bg   <- to_freq(cnt_bg)

# enrichment matrix
# Relative-entropy (KL) contributions: letter height = p_fg * log2(p_fg / p_bg)
kl <- p_fg * log2(p_fg / p_bg)

# per-position, per-base test (Fisher's exact, fg vs. bg)
L <- ncol(cnt_fg)
pval <- matrix(NA, 4, L, dimnames = dimnames(cnt_fg))
for (i in seq_len(L)) for (b in rownames(cnt_fg)) {
  tab <- matrix(c(cnt_fg[b, i], sum(cnt_fg[, i]) - cnt_fg[b, i],
                  cnt_bg[b, i], sum(cnt_bg[, i]) - cnt_bg[b, i]), nrow = 2)
  pval[b, i] <- fisher.test(tab)$p.value
}
padj <- matrix(p.adjust(pval, "BH"), nrow = 4, dimnames = dimnames(pval))

kl_sig <- kl
kl_sig[padj > 0.05] <- 0

ggseqlogo(kl, method = "custom", seq_type = "dna") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("G" = "#51bbd6", "A" = "#2e377c",
                               "T" = "#f8b867", "C" = "violetred4")) +
  labs(title = "A box: up vs. reference",
       x = "Position", y = "Enrichment (bits)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 7))

p <- ggseqlogo(kl, method = "custom", seq_type = "dna") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("G" = "#51bbd6", "A" = "#2e377c",
                               "T" = "#f8b867", "C" = "violetred4")) +
  labs(title = "A box: up vs. reference",
       x = "Position", y = "Enrichment (bits)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 7)) +
  ylim(-0.5,1)
ggsave(p, filename = "../output/A-box_enriched_up_lim.pdf", width = 6.5, height = 4.5, units = "in")


# B box
bbox_unique <- df %>%
  mutate(bbox_custom = toupper(bbox_custom)) %>%
  group_by(probePairID) %>%
  filter(n_distinct(bbox_custom) == 1) %>%
  slice(1) %>%
  ungroup()

fg_seqs <- bbox_unique$bbox_custom[bbox_unique$probePairID %in% df_toppp]
tested  <- unique(sub("-.*", "", cmp$trna))
bg_seqs <- bbox_unique$bbox_custom[!bbox_unique$probePairID %in% df_toppp]

cnt_fg <- count_matrix(fg_seqs)
cnt_bg <- count_matrix(bg_seqs)
p_fg   <- to_freq(cnt_fg)
p_bg   <- to_freq(cnt_bg)

kl <- p_fg * log2(p_fg / p_bg)
L <- ncol(cnt_fg)
pval <- matrix(NA, 4, L, dimnames = dimnames(cnt_fg))
for (i in seq_len(L)) for (b in rownames(cnt_fg)) {
  tab <- matrix(c(cnt_fg[b, i], sum(cnt_fg[, i]) - cnt_fg[b, i],
                  cnt_bg[b, i], sum(cnt_bg[, i]) - cnt_bg[b, i]), nrow = 2)
  pval[b, i] <- fisher.test(tab)$p.value
}
padj <- matrix(p.adjust(pval, "BH"), nrow = 4, dimnames = dimnames(pval))

kl_sig <- kl
kl_sig[padj > 0.05] <- 0

ggseqlogo(kl, method = "custom", seq_type = "dna") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("G" = "#51bbd6", "A" = "#2e377c",
                               "T" = "#f8b867", "C" = "violetred4")) +
  labs(title = "B box: up vs. reference",
       x = "Position", y = "Enrichment (bits)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 7))
p_bbox_enr

p <- ggseqlogo(kl, method = "custom", seq_type = "dna") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("G" = "#51bbd6", "A" = "#2e377c",
                               "T" = "#f8b867", "C" = "violetred4")) +
  labs(title = "B box: up vs. reference",
       x = "Position", y = "Enrichment (bits)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 7)) +
  ylim(-0.5,1)
ggsave(p, filename = "../output/B-box_enriched_up_lim.pdf", width = 6.5, height = 4.5, units = "in")


#### volcano for ISR-differential tRNAs with dot size set to A box enrichment score
guides_isr <- c("MTOR", "RPL14", "POLR1B", "MARS", "EIF2B4", "EIF2B3", "HARS", "AARS", 
                "TARS", "POLR3D", "QARS", "HSPA9", "ATF5", "BDP1")
so_p2$sample2 <- substr(so_p2$sample, 1, 3)
table(so_p2$sample2)
so_p2_ctrl <- subset(so_p2, subset = sample2 == "Con")
Idents(so_p2_ctrl) <- so_p2_ctrl$target
DefaultAssay(so_p2_ctrl) <- "tRNA_pp_valid"
deg_isr <- FindMarkers(
  so_p2_ctrl,
  ident.1 = guides_isr,
  ident.2 = "NTC"
)
deg_isr$gene <- rownames(deg_isr)
log2FC_cutoff <- 0.2
pval_cutoff <- 0.05
deg_isr <- deg_isr %>%
  mutate(
    significance = case_when(
      abs(avg_log2FC) >= log2FC_cutoff & p_val_adj < pval_cutoff ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

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

scores2 <- df %>%
  dplyr::transmute(probePairID,
            abox_custom = toupper(abox_custom),
            abox_score  = score_pssm(abox_custom, pssm)) %>%
  dplyr::group_by(probePairID) %>%
  dplyr::summarise(
    n_abox = dplyr::n_distinct(abox_custom, na.rm = TRUE),
    abox_score = if (n_abox == 1) dplyr::first(na.omit(abox_score)) else NA_real_,
    .groups = "drop"
  )
table(scores2$n_abox)

scores2 <- scores2 %>%
  dplyr::transmute(probePairID,
                   bbox_custom = toupper(bbox_custom),
                   bbox_score  = score_pssm(bbox_custom, pssm)) %>%
  dplyr::group_by(probePairID) %>%
  dplyr::summarise(
    n_abox = dplyr::n_distinct(abox_custom, na.rm = TRUE),
    abox_score = if (n_abox == 1) dplyr::first(na.omit(abox_score)) else NA_real_,
    .groups = "drop"
  )
table(scores2$n_bbox)

deg_isr <- deg_isr %>%
  select(-any_of("abox_score")) %>%
  mutate(probePairID = sub("-.*", "", gene)) %>%
  left_join(scores2, by = "probePairID") %>%
  arrange(significance, abox_score)

deg_isr$regulation <- ifelse(
  deg_isr$significance != "Significant",
  "Not significant",
  ifelse(deg_isr$avg_log2FC > 0, "Upregulated", "Downregulated")
)
quantile(
  deg_isr$abox_score,
  probs = c(0, .01, .05, .10, .25, .5, .75, .90, .95, .99, 1),
  na.rm = TRUE
)
lims <- quantile(deg_isr$abox_score, c(0.01, 0.99), na.rm = TRUE)
ggplot(deg_isr, aes(x = avg_log2FC, y = -1*log10(p_val_adj), color = regulation, size = abox_score)) +
  geom_point() +
  scale_color_manual(values = c("Not significant" = "lightgrey",
                                "Upregulated" = "violetred4",
                                "Downregulated" = "darkblue")) +
  theme_minimal() +
  scale_size_continuous(range = c(-6.34, 4.35), # 10th and 90th perc
                        name = "A-box match to\nupregulated set\n(log2 odds)") +
  geom_text_repel(data = filter(deg_isr, significance == "Significant"), aes(label = gene),
                  size = 3.2, max.overlaps = 15, box.padding = 0.3, point.padding = 0.2) + 
  pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/volcano_control_ISR-vs-notISR2.pdf", width = 4, height = 4.5, dpi = 400)


# correlation plot
ggplot(deg_isr, aes(abox_score, avg_log2FC)) +
  geom_point() + pretty_plot() + L_border()

p <- ggplot(deg_isr %>% filter(p_val_adj < 0.05), aes(abox_score, avg_log2FC)) +
  geom_point(size=2) +
  geom_text_repel(data=deg_isr %>% filter(p_val_adj < 0.05), aes(label=gene), max.overlaps = 8) +
  geom_smooth(method = "lm", se = TRUE, colour = "violetred4", linewidth = 0.6) +
  stat_cor(method = "pearson", cor.coef.name = "rho",
           label.x.npc = "left", label.y.npc = "top", size = 3) +
  labs(x = "A-box motif score (log2 odds)", y = "avg log2FC (ISR vs NTC)") +
  pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/cor_volcano_ISR_Abox-score_label.pdf", width = 4.5, height = 4.5, dpi = 400)
