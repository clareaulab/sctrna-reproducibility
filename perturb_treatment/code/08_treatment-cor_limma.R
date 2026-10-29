library(dplyr)
library(tidyr)
library(Seurat)
library(tibble)
library(BuenColors)
library(ggrepel)
library(edgeR)
library(limma)
library(Matrix)
filter <- dplyr::filter
select <- dplyr::select

sample_colors <- c(
  Con = "#fdde8f",
  Tor = "#4FBBD6",
  AAd = "#D65147"
)

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
so_p2 <- subset(so_p2, subset = target != "TRMT5")
table(so_p2$identity_NTC)
so_p2$sample2 <- substr(so_p2$sample, 1, 3)
table(so_p2$sample2)
so_p2$identity2 <- gsub("non-targeting_non-targeting_non-targeting", "NTC", so_p2$identity)
table(so_p2$identity2)


# get NTC cells from control and AAd
md <- so_p2@meta.data
md$guide_id <- md[["identity2"]]
md$is_ntc   <- grepl("_NTC_[AB]$", md$guide_id)

keep_cells <- md$is_ntc & md$sample2 %in% c("Con", "AAd")
so_n <- so_p2[, keep_cells]
md   <- md[keep_cells, ]

# pseudobulk per NTC guide x treatment
grp <- factor(paste(md$guide_id, md$sample2, sep = "|"))
pb  <- as.matrix(so_n[["tRNA_pp_valid"]]$counts %*% t(fac2sparse(grp)))

meta <- data.frame(pb_id  = levels(grp),
                   ncells = as.vector(table(grp)),
                   row.names = levels(grp))
meta$guide_id  <- sub("\\|.*$", "", meta$pb_id)
meta$treatment <- sub("^.*\\|", "", meta$pb_id)

# keep guides present in both treatments
meta <- meta[meta$ncells >= 100, ]
paired <- names(which(table(meta$guide_id) == 2))
meta <- meta[meta$guide_id %in% paired, ]
pb   <- pb[, rownames(meta)]
length(paired) # 18
tapply(meta$ncells, meta$treatment, summary)

## paired design: guide as block, control as reference
meta$treatment <- factor(meta$treatment, levels = c("Con", "AAd"))
meta$guide_id  <- factor(meta$guide_id)
design <- model.matrix(~ guide_id + treatment, data = meta)

y <- DGEList(pb)
dim(y) # 429  36
y <- y[filterByExpr(y, group = meta$treatment), , keep.lib.sizes = FALSE]
dim(y) # 216  36
y <- calcNormFactors(y)
v   <- voomWithQualityWeights(y, design, plot = TRUE)
fit <- eBayes(lmFit(v, design))

res_ntc_AAd <- topTable(fit, coef = "treatmentAAd", n = Inf)
summary(decideTests(fit[, "treatmentAAd"]))

## same for Torin-1
md <- so_p2@meta.data
md$guide_id <- md[["identity2"]]
md$is_ntc   <- grepl("_NTC_[AB]$", md$guide_id)
keep_cells <- md$is_ntc & md$sample2 %in% c("Con", "Tor")
so_n <- so_p2[, keep_cells]
md   <- md[keep_cells, ]

# pseudobulk per NTC guide x treatment
grp <- factor(paste(md$guide_id, md$sample2, sep = "|"))
pb  <- as.matrix(so_n[["tRNA_pp_valid"]]$counts %*% t(fac2sparse(grp)))

meta <- data.frame(pb_id  = levels(grp),
                   ncells = as.vector(table(grp)),
                   row.names = levels(grp))
meta$guide_id  <- sub("\\|.*$", "", meta$pb_id)
meta$treatment <- sub("^.*\\|", "", meta$pb_id)

# drop small pseudobulks, keep guides present in both treatments
meta <- meta[meta$ncells >= 10, ]
paired <- names(which(table(meta$guide_id) == 2))
meta <- meta[meta$guide_id %in% paired, ]
pb   <- pb[, rownames(meta)]
length(paired) # number of NTC guide pairs
tapply(meta$ncells, meta$treatment, summary)

meta$treatment <- factor(meta$treatment, levels = c("Con", "Tor"))
meta$guide_id  <- factor(meta$guide_id)
design <- model.matrix(~ guide_id + treatment, data = meta)

y <- DGEList(pb)
dim(y) # 429  36
y <- y[filterByExpr(y, group = meta$treatment), , keep.lib.sizes = FALSE]
dim(y) # 201  36
y <- calcNormFactors(y)
v   <- voomWithQualityWeights(y, design, plot = TRUE)
fit <- eBayes(lmFit(v, design))

res_ntc_Tor <- topTable(fit, coef = "treatmentTor", n = Inf)
summary(decideTests(fit[, "treatmentTor"]))

#plotMDS(v, col = ifelse(meta$treatment == "Tor", "violetred4", "grey40"), pch = 16)
#plotMDS(v, col = ifelse(meta$treatment == "AAd", "violetred4", "grey40"), pch = 16)

res_ntc_AAd <- res_ntc_AAd %>% 
  rownames_to_column("gene")
res_ntc_Tor <- res_ntc_Tor %>% 
  rownames_to_column("gene")

ggplot(res_ntc_AAd, aes(x=logFC, y=adj.P.Val)) +
  geom_point() + geom_text_repel(aes(label=gene)) +
  pretty_plot() + L_border()
ggplot(res_ntc_Tor, aes(x=logFC, y=adj.P.Val)) +
  geom_point() + geom_text_repel(aes(label=gene)) +
  pretty_plot() + L_border()

rownames(res_ntc_AAd) <- res_ntc_AAd$gene
rownames(res_ntc_Tor) <- res_ntc_Tor$gene

prep <- function(df, suffix) {
  if (!"trna" %in% colnames(df)) df$trna <- rownames(df)
  df %>% select(trna, logFC, adj.P.Val) %>%
    rename_with(~ paste0(.x, "_", suffix), -trna)
}

cmp <- inner_join(prep(res_ntc_AAd, "aad"), prep(res_ntc_Tor, "torin"), by = "trna")
nrow(cmp) # tRNAs tested in both: 201

# categories: significant in both/one/neither
alpha <- 0.05
cmp <- cmp %>%
  mutate(
    sig_aad   = adj.P.Val_aad   < alpha,
    sig_torin = adj.P.Val_torin < alpha,
    category = case_when(
      sig_aad & sig_torin & sign(logFC_aad) == sign(logFC_torin) ~ "both, same direction",
      sig_aad & sig_torin ~ "both, opposite",
      sig_aad ~ "AA depletion only",
      sig_torin ~ "Torin-1 only",
      TRUE ~ "neither"
    ),
    diff = logFC_aad - logFC_torin
  )
#write.csv(cmp, "../output/DEGs_treatment-control_limma_261016.csv")
table(cmp$category)

# correlation
r_p <- cor(cmp$logFC_aad, cmp$logFC_torin, method = "pearson")
r_s <- cor(cmp$logFC_aad, cmp$logFC_torin, method = "spearman")

# label the most treatment-specific tRNAs (largest difference, significant in at least one)
to_label <- cmp %>%
  filter(sig_aad | sig_torin) %>%
  slice_max(abs(diff), n = 40)
lim <- max(abs(c(cmp$logFC_aad, cmp$logFC_torin)))
p <- ggplot(cmp, aes(logFC_torin, logFC_aad)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = category, size = -log10(adj.P.Val_aad)), alpha = 0.8) +
  geom_text_repel(data = to_label, aes(label = trna),
                  size = 3, max.overlaps = Inf, min.segment.length = 0) +
  scale_colour_manual(values = c(
    "both, same direction" = "violetred4",
    "both, opposite" = "darkorange",
    "AA depletion only" = "#78b79d",
    "Torin-1 only" = "#5AA9C8",
    "neither" = "grey80")) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  #coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  annotate("text", x = -lim, y = lim, hjust = 0, vjust = 1, size = 3.5,
           label = sprintf("Pearson r = %.2f\nSpearman ρ = %.2f", r_p, r_s)) +
  labs(x = "log2FC Torin-1 vs untreated (NTC)",
       y = "log2FC AA depletion vs untreated (NTC)",
       colour = NULL) + theme_classic()
ggsave(filename = "../output/scatter_treatment-fc_vsNTC4_261016_4.pdf", plot = p, width = 7, height = 4, units = "in")
