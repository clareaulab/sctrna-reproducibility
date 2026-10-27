library(dplyr)
library(Seurat)
library(BuenColors)

filter <- dplyr::filter
select <- dplyr::select

so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_pseudotime_260907.rds")

so_cd34m <- subset(so_cd34, subset = condition %in% c("D0", "D7_Mono"))
DimPlot(so_cd34m, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34m <- subset(so_cd34m, subset = seurat_clusters %in% c(0,2,4,11))
DimPlot(so_cd34m, reduction = "fdl", group.by = c("condition", "seurat_clusters"))

# CD14 expression barplot over time
FeatureScatter(
  object = so_cd34m, 
  feature1 = "palantir_pseudotime", 
  feature2 = "CD14"
) + geom_smooth(method = "gam")

# color gradient bar
df_cd14 <- FetchData(so_cd34m, vars = c("palantir_pseudotime", "CD14"))
p <- ggplot(df_cd14, aes(x = palantir_pseudotime, y = 1, fill = rna_CD14)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("lightgrey", "firebrick3", "violetred4")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) + labs(title = "CD14 Expression Density", fill = "Expression") + ylim(0, 2)
cowplot::ggsave2(p, file = "../output/CD14_bar.pdf", width = 4, height = 4)

# panel for all ACs
DefaultAssay(so_cd34m) <- "tRNA_ac"
trna_features <- rownames(so_cd34m[["tRNA_ac"]])
df_trend_long <- map_dfr(trna_features, function(trna) {
  tRNA_expr <- FetchData(
    so_cd34m,
    vars = trna,
    assay = "tRNA_aa",
    layer = "data"
  )[, 1]
  data.frame(
    Mono = so_cd34m$Mono,
    MNDA = MNDA_expr,
    tRNA = tRNA_expr,
    feature = trna
  )
})
ggplot(df_trend_long, aes(x = Mono)) +
  geom_smooth(aes(y = scale_MNDA(MNDA)), color = "red") +
  geom_smooth(aes(y = tRNA), color = "royalblue") +
  facet_wrap(~ feature, scales = "free_y") +
  labs(
    y = "Expression / Abundance",
    x = "Palantir Mono"
  ) +
  theme_minimal() +
  pretty_plot() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

# plot Asn AAC tRNAs
DefaultAssay(so_cd34m) <- "tRNA_pp_valid"
tRNA1 <- FetchData(so_cd34m, vars = "PP32-Asn-AAC", assay = "tRNA_pp_valid", layer = "data")[,1] # constant/up
tRNA2 <- FetchData(so_cd34m, vars = "PP33-Asn-AAC", assay = "tRNA_pp_valid", layer = "data")[,1] # down
tRNA3 <- FetchData(so_cd34m, vars = "PP35-Asn-AAC", assay = "tRNA_pp_valid", layer = "data")[,1] # down
tRNA4 <- FetchData(so_cd34m, vars = "PP36-Asn-AAC", assay = "tRNA_pp_valid", layer = "data")[,1] # down
tRNA5 <- FetchData(so_cd34m, vars = "PP41-Asn-AAC", assay = "tRNA_pp_valid", layer = "data")[,1] # down
tRNA6 <- FetchData(so_cd34m, vars = "PP42-Asn-AAC", assay = "tRNA_pp_valid", layer = "data")[,1] # constant
df_trend <- data.frame(
  pseudotime = so_cd34m$palantir_pseudotime,
  tRNA1 = tRNA1,
  tRNA3 = tRNA3,
  tRNA4 = tRNA4,
  tRNA5 = tRNA5,
  tRNA6 = tRNA6
)
p <- ggplot(df_trend, aes(x = pseudotime)) +
  geom_smooth(aes(y = tRNA1), color = "violetred4") +
  geom_smooth(aes(y = tRNA2), color = "cornflowerblue") +
  geom_smooth(aes(y = tRNA3), color = "blue4") +
  geom_smooth(aes(y = tRNA4), color = "#4c7f3f") +
  geom_smooth(aes(y = tRNA5), color = "lightblue") +
  geom_smooth(aes(y = tRNA6), color = "mediumorchid") +
  labs(y = "Expression / Abundance", x = "Palantir pseudotime") +
  theme_minimal() + 
  pretty_plot() +
  theme(legend.position = "bottom")
#cowplot::ggsave2(p, file = "../output/trajectory_Mono_pseudotime_Asn.pdf", width = 4, height = 4)
