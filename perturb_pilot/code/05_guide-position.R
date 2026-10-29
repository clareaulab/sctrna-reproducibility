library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(tibble)
library(BuenColors)
library(ggrepel)
library(ggstar)
library(purrr)

filter <- dplyr::filter
select <- dplyr::select

RNATtoU <- function(invec) {
  chartr(old = "Tt", new = "UU", invec)
}
DNAtoRNA <- function(dna_vec) {
  sapply(dna_vec, function(dna) {
    dna <- toupper(dna)
    comp <- chartr("ATCG", "TAGC", dna)
    rev_comp <- paste(rev(strsplit(comp, "")[[1]]), collapse = "")
    gsub("T", "U", rev_comp)
  })
}

so_pp <- readRDS(file = "../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA.rds")
dim(so_pp) #18099 10899
unique(so_pp$identity_NTC)
so_pp$target <- str_remove(so_pp$identity_NTC, "_guide[123]$")
unique(so_pp$target)


tRNA_seq <- read_csv("../data/tRNA_sequences_Abox-Bbox-custom_quant.csv")
#target_tRNA <- unique(gsub("tRNA-", "", so_pp$target))
target_tRNA <- unique(so_pp$target)
tRNA_seq$name2 <- paste0("tRNA-", tRNA_seq$name)


tRNA_seq_want <- tRNA_seq %>% filter(name2 %in% target_tRNA)
pool_seq <- read_tsv("../data/pool1_tRNAseqs.tsv")
pool_seq <- pool_seq %>%
  group_by(Symbol) %>%
  dplyr::mutate(Symbol_guide = paste0(Symbol, "_guide", row_number())) %>%
  ungroup()

plot_df <- pool_seq %>% 
  left_join(tRNA_seq_want %>% select(name2, sequences, chr, probePairID, expected_anticodons),
            by = c("Symbol"="name2"))
plot_df <- plot_df[1:62,]
plot_df$guide_start <- mapply(
  function(guide, seq) regexpr(guide, seq, fixed = TRUE)[1],
  plot_df$sgRNA.Sequence,
  plot_df$sequences
)
plot_df$guide_location <- (plot_df$guide_start - 1) /
  (nchar(plot_df$sequences) - 1)

plot_df_filt <- plot_df %>% filter(guide_start != '-1')


# map to fold changes of guide effect
all_diff_out2 <- read_csv(file = "../output/diff_tRNA.csv")
all_diff_out2$name <- gsub("_guide[123]", "", all_diff_out2$gRNA)
gene_want <- intersect(
  unique(all_diff_out2$name),
  unique(plot_df_filt$Symbol)
)
all_diff_out3 <- all_diff_out2 %>% filter(name %in% gene_want)
plot_df2 <- all_diff_out3 %>% 
  left_join(plot_df_filt %>% select(Symbol, Symbol_guide, chr, probePairID, guide_location),
            by = c("name"="Symbol"), relationship = "many-to-many")
plot_df2_filt <- plot_df2[sub("-.*", "", plot_df2$probe) == plot_df2$probePairID,]
plot_df2_filt$AC_want <- sub(
  "^[^-]+-[^-]+-([^-]+).*",
  "\\1",
  plot_df2_filt$name
)
plot_df2_filt$AC_have <- sub(
  "^[^-]+-[^-]+-([^-]+).*",
  "\\1",
  plot_df2_filt$probe
)
plot_df2_filt$AC_want_rc <- DNAtoRNA(plot_df2_filt$AC_want)
plot_df2_filt2 <- plot_df2_filt %>%
  filter(AC_want_rc == AC_have)
plot_df2_filt3 <- plot_df2_filt %>%
  filter(gRNA == Symbol_guide)

ggplot(plot_df2_filt3, aes(x=guide_location, y=avg_log2FC)) +
  geom_point()+
  geom_text_repel(data=plot_df2_filt3, aes(label=name)) +
  pretty_plot() + L_border() +
  stat_cor(aes(label = paste("R ==", signif(after_stat(r), 2), "~`,`~",
                             "p ==", signif(after_stat(p), 2))),
           method = "pearson", label.x.npc = "left", label.y.npc = "top", parse = TRUE)





