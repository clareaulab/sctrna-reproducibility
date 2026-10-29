library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
filter <- dplyr::filter
select <- dplyr::select

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
table(so_p2$identity_NTC)
table(so_p2$identity)
so_p2$identity2 <- gsub("non-targeting_non-targeting_non-targeting", "NTC", so_p2$identity)
so_p2$sample2 <- substr(so_p2$sample, 1, 3)

# get NTC cells in the three treatments
keep <- grepl("_NTC_[AB]$", so_p2$identity2) & so_p2$sample2 %in% c("Con", "AAd", "Tor")
# keep <- so_p2$target == "NTC" & so_p2$sample2 %in% c("Con", "AAd", "Tor")
cells <- colnames(so_p2)[keep]
grp   <- factor(so_p2$sample2[keep], levels = c("Con", "AAd", "Tor"))
table(grp)
#  Con  AAd  Tor 
#  2975 5564 2975

cnt <- so_p2[["tRNA_aa"]]$counts[, cells]

# mean raw counts per cell
mean_counts <- sapply(split(cells, grp), function(cc) rowMeans(cnt[, cc, drop = FALSE]))

# mean per-cell proportions
tot <- colSums(cnt)
ok <- tot > 0
prop <- t(t(cnt[, ok]) / tot[ok])
mean_prop <- sapply(split(colnames(prop), grp[ok]),
                    function(cc) rowMeans(prop[, cc, drop = FALSE]))

# long format, isoacceptors ordered by overall abundance
to_long <- function(m, value) {
  as.data.frame(m) %>%
    tibble::rownames_to_column("aa") %>%
    pivot_longer(-aa, names_to = "treatment", values_to = value) %>%
    mutate(treatment = factor(treatment, levels = c("Con", "AAd", "Tor")))
}
aa_order <- names(sort(rowMeans(mean_prop), decreasing = TRUE))

df_prop  <- to_long(mean_prop, "prop")   %>% mutate(aa = factor(aa, levels = aa_order))
df_count <- to_long(mean_counts, "count")  %>% mutate(aa = factor(aa, levels = aa_order))

n_label <- paste0(levels(grp), "\n(n = ", as.vector(table(grp)), ")")

top_n <- 11
aa_rank <- names(sort(rowMeans(mean_counts), decreasing = TRUE))
top_aa  <- aa_rank[seq_len(top_n)]
rest_aa <- setdiff(aa_rank, top_aa)
pal_met <- as.character(MetBrewer::met.brewer("Renoir", 11))
aa_cols <- c(
  setNames(pal_met, top_aa),
  setNames(grey.colors(length(rest_aa), start = 0.35, end = 0.9), rest_aa)
)

# colored AAs first, greys after
df_prop$aa  <- factor(df_prop$aa,  levels = c(top_aa, rest_aa))
df_count$aa <- factor(df_count$aa, levels = c(top_aa, rest_aa))

ggplot(df_count, aes(treatment, count, fill = aa)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = n_label) +
  scale_fill_manual(values = aa_cols) +
  labs(x = NULL, y = "Mean tRNA counts per cell", fill = "Amino acid") +
  theme_classic()
#ggsave(filename = "../output/barplot_tRNAaa-counts_treatments2.pdf", plot = p, width = 3.4, height = 2.8, units = "in")

# composition
ggplot(df_prop, aes(treatment, prop, fill = aa)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_discrete(labels = n_lab) +
  scale_fill_manual(values = aa_cols) +
  labs(x = NULL, y = "Mean fraction of tRNA counts per cell", fill = "Amino acid") +
  theme_classic()
#ggsave(filename = "../output/barplot_tRNAaa-counts_treatments_composition.pdf", plot = p, width = 3.4, height = 3.4, units = "in")


# correlation of AAd/Con and Tor/Con
wide <- df_prop %>%
  pivot_wider(names_from = treatment, values_from = prop)

cor.test(wide$AAd, wide$Con, method = "spearman")
cor.test(wide$Tor, wide$Con, method = "spearman")

p <- ggplot(wide, aes(x=log10(Con), y=log10(AAd))) +
  geom_point(size=2.5) + geom_text_repel(data=wide, aes(label=aa), size=4.5) +
  pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/corplot_composition_AAd-Con.pdf", width = 3.8, height = 4)

p <- ggplot(wide, aes(x=log10(Con), y=log10(Tor))) +
  geom_point(size=2.5) + geom_text_repel(data=wide, aes(label=aa), size=4.5) +
  pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/corplot_composition_Tor-Con.pdf", width = 3.8, height = 4)
