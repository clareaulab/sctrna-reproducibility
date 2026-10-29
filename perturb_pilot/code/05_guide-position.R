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
unique(so_pp$identity_NTC)
so_pp$target <- str_remove(so_pp$identity_NTC, "_guide[123]$")
unique(so_pp$target)

grna_df <- read.csv("../data/crispick_output.csv")
grna_df <- grna_df %>%
  group_by(Symbol) %>%
  dplyr::mutate(name = paste0(Symbol, "_guide", row_number())) %>%
  ungroup()
grna_df <- grna_df %>% 
  select(name, Symbol, Locus, Strand.of.Target, Strand.of.sgRNA, sgRNA..Cut..Position, GCprop)
grna_df$start <- sub(".*:(\\d+)-.*", "\\1", grna_df$Locus)
grna_df$end   <- sub(".*-(\\d+)$", "\\1", grna_df$Locus)
grna_df$start <- as.numeric(grna_df$start)
grna_df$end   <- as.numeric(grna_df$end)
grna_df$tRNA_length <- grna_df$end - grna_df$start
grna_df$cut_dist <- grna_df$sgRNA..Cut..Position - grna_df$start

grna_df$cut_position <- ifelse(
  grna_df$Strand.of.Target == "+",
  (grna_df$cut_dist + 1) / grna_df$tRNA_length,
  (grna_df$tRNA_length - grna_df$cut_dist) / grna_df$tRNA_length
)

# map to fold changes of guide effect
all_diff_out2 <- read_csv(file = "../output/diff_tRNA.csv")
all_diff_out2$name <- gsub("_guide[123]", "", all_diff_out2$gRNA)
gene_want <- intersect(
  unique(all_diff_out2$gRNA),
  unique(grna_df$name)
)
all_diff_out3 <- all_diff_out2 %>% filter(gRNA %in% gene_want)
plot_df <- all_diff_out3 %>% 
  left_join(grna_df %>% select(name, cut_position, GCprop),
            by = c("gRNA"="name"), relationship = "many-to-many")

# get probe pair IDs for guide targets
tRNA_seq <- read_csv("../data/tRNA_sequences_Abox-Bbox-custom_quant.csv")
plot_df$name2 <- gsub("tRNA-", "", plot_df$name)
plot_df2 <- plot_df %>% 
  left_join(tRNA_seq %>% select(name, chr, "pp_want"=probePairID, expected_anticodons),
            by = c("name2"="name"))
plot_df2$AC_want <- sub("^[^-]+-[^-]+-([^-]+).*", "\\1", plot_df2$name)
plot_df2$AC_have <- sub("^[^-]+-[^-]+-([^-]+).*", "\\1", plot_df2$probe)
plot_df2$AC_want_rc <- DNAtoRNA(plot_df2$AC_want)
plot_df2_filt <- plot_df2 %>% filter(AC_want_rc == AC_have)

# get only features with matching PP for each guide
plot_df2_filt$pp_have <- sub("-.*", "", plot_df2_filt$probe)
plot_df2_filt2 <- plot_df2_filt %>%
  filter(pp_want == pp_have)

## this plot
ggplot(plot_df2_filt2, aes(x=cut_position, y=avg_log2FC, color = -1*log10(p_val_adj), size = pct.2)) +
  geom_point() + scale_color_viridis_c(option = "viridis") +
  #geom_text_repel(data=plot_df2_filt2, aes(label=gRNA)) +
  pretty_plot(fontsize = 20) + L_border()
#cowplot::ggsave2(p, file = "../output/corplot_GGlocation_FC.pdf", width = 10, height = 4.8)

# get % with log2 FC below -0.5
mean(plot_df2_filt2$avg_log2FC < -0.5, na.rm = TRUE) * 100 # 66.66667
mean(plot_df2_filt2$avg_log2FC < 0, na.rm = TRUE) * 100 # 92.59259

ggplot(plot_df2_filt2, aes(x=GCprop, y=avg_log2FC)) +
  geom_point()+
  geom_text_repel(data=plot_df2_filt2, aes(label=gRNA)) +
  pretty_plot() + L_border() +
  stat_cor(aes(label = paste("R ==", signif(after_stat(r), 2), "~`,`~",
                             "p ==", signif(after_stat(p), 2))),
           method = "pearson", label.x.npc = "left", label.y.npc = "top", parse = TRUE)

plot_df2_filt2$chr <- as.numeric(plot_df2_filt2$chr)
ggplot(plot_df2_filt2, aes(x=chr, y=avg_log2FC)) +
  geom_point()+
  geom_text_repel(data=plot_df2_filt2, aes(label=gRNA)) +
  pretty_plot() + L_border()
