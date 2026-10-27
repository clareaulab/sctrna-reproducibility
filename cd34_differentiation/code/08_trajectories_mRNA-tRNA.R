library(ggplot2)
library(dplyr)
library(Matrix)
library(mgcv)
filter <- dplyr::filter
select <- dplyr::select

so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_pseudotime_260907.rds")

# subset to stem/ery d7
so_cd34e <- subset(so_cd34, subset = condition %in% c("D0", "D7_Ery"))
DimPlot(so_cd34e, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34e <- subset(so_cd34e, subset = seurat_clusters %in% c(0,3,5,7,10,12))
DimPlot(so_cd34e, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)

# stem/mega d7
so_cd34k <- subset(so_cd34, subset = condition %in% c("D0", "D7_Mega"))
DimPlot(so_cd34k, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34k <- subset(so_cd34k, subset = seurat_clusters %in% c(0,5,8))
DimPlot(so_cd34k, reduction = "fdl", group.by = c("seurat_clusters"))

# stem/mono d7
so_cd34m <- subset(so_cd34, subset = condition %in% c("D0", "D7_Mono"))
DimPlot(so_cd34m, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34m <- subset(so_cd34m, subset = seurat_clusters %in% c(0,1,2,4,11))
DimPlot(so_cd34m, reduction = "fdl", group.by = c("condition", "seurat_clusters"))

# look at lineage markers
DefaultAssay(so_cd34) <- "RNA"
# ery
FeaturePlot(so_cd34, features = c("HBB", "ALAS2", "GATA1", "KLF1"), cols = c("lightgray", "violetred4"), 
            min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
# mono
FeaturePlot(so_cd34, features = c("CD14", "LYZ", "CSF1R", "CD68"), cols = c("lightgray", "violetred4"), 
            min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")
# mega
FeaturePlot(so_cd34, features = c("PPBP", "ITGA2B", "PF4", "GP9"), cols = c("lightgray", "violetred4"), 
            min.cutoff = "q05", max.cutoff = "q95", reduction = "fdl")


## trajectories for mono mRNA and tRNA markers
cd14_expr <- FetchData(so_cd34m, vars = "CD14", assay = "RNA", layer = "data")[, 1]
lyz_expr <- FetchData(so_cd34m, vars = "LYZ", assay = "RNA", layer = "data")[, 1]
DefaultAssay(so_cd34m) <- "tRNA_pp_valid"
trna_features <- rownames(so_cd34m[["tRNA_pp_valid"]])
df_trend_long <- map_dfr(trna_features, function(trna) {
  tRNA_expr <- FetchData(
    so_cd34m,
    vars = trna,
    assay = "tRNA_pp_valid",
    layer = "data"
  )[, 1]
  data.frame(
    cd14 = cd14_expr,
    lyz = lyz_expr,
    tRNA = tRNA_expr,
    pseudotime = so_cd34m$palantir_pseudotime,
    feature = trna
  )
})
df_trend_long <- df_trend_long %>%
  rownames_to_column("cell")
df_trend_summary <- df_trend_long %>%
  group_by(feature) %>%
  reframe(tidy(cor.test(tRNA, pseudotime))) %>%
  select(feature, R = estimate, p = p.value) %>% 
  mutate(p_adj = p.adjust(p, method = "BH"))

ggplot(df_trend_summary, aes(x = R, y = -1*log10(p_adj))) +
  geom_point() +
  theme_minimal() + labs(title = "Correlation over pseudotime, mono") +
  geom_text_repel(data = df_trend_summary,
    aes(label = feature), size = 3, max.overlaps = 15, box.padding = 0.3, point.padding = 0.2) +
  pretty_plot() + L_border()

# look at some hits
selected <- c("PP22-Arg-AGG", "PP45-Asp-GAC")
minmax <- function(x) (x - min(x, na.rm = TRUE)) /
  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
df_sel <- df_trend_long %>% filter(feature %in% selected)
# common x grid to predict on
grid <- data.frame(pseudotime = seq(min(df_sel$pseudotime, na.rm = TRUE),
                                    max(df_sel$pseudotime, na.rm = TRUE),
                                    length.out = 200))

# one fitted, min-max-scaled tRNA curve per feature
fit_trna <- df_sel %>%
  group_by(feature) %>%
  group_modify(~{
    m <- gam(tRNA ~ s(pseudotime, bs = "cs"), data = .x)
    data.frame(pseudotime = grid$pseudotime,
               y = minmax(predict(m, newdata = grid)))
  }) %>% ungroup()

# CD14/LYZ fit once (identical across blocks)
cd14_block <- df_sel %>% filter(feature == selected[1])
m_cd14 <- gam(cd14 ~ s(pseudotime, bs = "cs"), data = cd14_block)
fit_cd14 <- data.frame(pseudotime = grid$pseudotime,
                       y = minmax(predict(m_cd14, newdata = grid)),
                       feature = "CD14")
lyz_block <- df_sel %>% filter(feature == selected[1])
m_lyz <- gam(lyz ~ s(pseudotime, bs = "cs"), data = lyz_block)
fit_lyz <- data.frame(pseudotime = grid$pseudotime,
                      y = minmax(predict(m_lyz, newdata = grid)),
                      feature = "LYZ")

p <- ggplot(mapping = aes(x = pseudotime, y = y, color = feature)) +
  geom_line(data = fit_trna, linewidth = 1) +
  geom_line(data = fit_cd14, linewidth = 1) +
  geom_line(data = fit_lyz, linewidth = 1) +
  scale_color_manual(values = c(
    "PP22-Arg-AGG"  = "#1b9e77",
    "PP45-Asp-GAC"  = "green3",
    "CD14"          = "red",
    "LYZ"          = "violetred4"
  )) +
  labs(y = "Scaled expression", x = "Pseudotime", color = NULL) +
  theme_minimal() + pretty_plot() + theme(strip.text = element_text(size = 8))
cowplot::ggsave2(p, file = "../output/trajectory_mono_mRNA-tRNA2.pdf", width = 5, height = 3.4)


## trajectories for mega mRNA and tRNA markers
m1_expr <- FetchData(so_cd34k, vars = "PF4", assay = "RNA", layer = "data")[, 1]
m2_expr <- FetchData(so_cd34k, vars = "GP9", assay = "RNA", layer = "data")[, 1]
DefaultAssay(so_cd34k) <- "tRNA_pp_valid"
trna_features <- rownames(so_cd34k[["tRNA_pp_valid"]])
df_trend_long <- map_dfr(trna_features, function(trna) {
  tRNA_expr <- FetchData(
    so_cd34k,
    vars = trna,
    assay = "tRNA_pp_valid",
    layer = "data"
  )[, 1]
  data.frame(
    PF4 = m1_expr,
    GP9 = m2_expr,
    pseudotime = so_cd34k$palantir_pseudotime,
    feature = trna
  )
})
df_trend_long <- df_trend_long %>%
  rownames_to_column("cell")
df_trend_summary <- df_trend_long %>%
  group_by(feature) %>%
  reframe(tidy(cor.test(tRNA, pseudotime))) %>%
  select(feature, R = estimate, p = p.value) %>% 
  mutate(p_adj = p.adjust(p, method = "BH"))

ggplot(df_trend_summary, aes(x = R, y = -1*log10(p_adj))) +
  geom_point() +
  theme_minimal() + labs(title = "Correlation over pseudotime, mega") +
  geom_text_repel(data = df_trend_summary,
                  aes(label = feature), size = 3, max.overlaps = 15, box.padding = 0.3, point.padding = 0.2) +
  pretty_plot() + L_border()

# look at some hits
selected <- c("PP79-Ile-AUC", "PP82-Ile-AUC")
minmax <- function(x) (x - min(x, na.rm = TRUE)) /
  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
df_sel <- df_trend_long %>% filter(feature %in% selected)
# common x grid to predict on
grid <- data.frame(pseudotime = seq(min(df_sel$pseudotime, na.rm = TRUE),
                                    max(df_sel$pseudotime, na.rm = TRUE),
                                    length.out = 200))

# one fitted, min-max-scaled tRNA curve per feature
fit_trna <- df_sel %>%
  group_by(feature) %>%
  group_modify(~{
    m <- gam(tRNA ~ s(pseudotime, bs = "cs"), data = .x)
    data.frame(pseudotime = grid$pseudotime,
               y = minmax(predict(m, newdata = grid)))
  }) %>% ungroup()

# PF4/GP9 fit once (identical across blocks)
pf4_block <- df_sel %>% filter(feature == selected[1])
m_pf4 <- gam(PF4 ~ s(pseudotime, bs = "cs"), data = pf4_block)
fit_pf4 <- data.frame(pseudotime = grid$pseudotime,
                       y = minmax(predict(m_pf4, newdata = grid)),
                       feature = "PF4")
gp9_block <- df_sel %>% filter(feature == selected[1])
m_gp9 <- gam(GP9 ~ s(pseudotime, bs = "cs"), data = gp9_block)
fit_gp9 <- data.frame(pseudotime = grid$pseudotime,
                      y = minmax(predict(m_gp9, newdata = grid)),
                      feature = "GP9")

p <- ggplot(mapping = aes(x = pseudotime, y = y, color = feature)) +
  geom_line(data = fit_trna, linewidth = 1) +
  geom_line(data = fit_pf4, linewidth = 1) +
  geom_line(data = fit_gp9, linewidth = 1) +
  scale_color_manual(values = c(
    "PP79-Ile-AUC"  = "#1b9e77",
    "PP82-Ile-AUC"  = "green3",
    "GP9"          = "red",
    "PF4"          = "violetred4"
  )) +
  labs(y = "Scaled expression", x = "Pseudotime", color = NULL) +
  theme_minimal() + pretty_plot() + theme(strip.text = element_text(size = 8))
cowplot::ggsave2(p, file = "../output/trajectory_mega_mRNA-tRNA2.pdf", width = 5, height = 3.4)


## trajectories for ery mRNA and tRNA markers
hbb_expr <- FetchData(so_cd34e, vars = "HBB", assay = "RNA", layer = "data")[, 1]
alas2_expr <- FetchData(so_cd34e, vars = "ALAS2", assay = "RNA", layer = "data")[, 1]
x_expr <- FetchData(so_cd34e, vars = "CD36", assay = "RNA", layer = "data")[, 1]
DefaultAssay(so_cd34e) <- "tRNA_pp_valid"
trna_features <- rownames(so_cd34e[["tRNA_pp_valid"]])
df_trend_long <- map_dfr(trna_features, function(trna) {
  tRNA_expr <- FetchData(
    so_cd34e,
    vars = trna,
    assay = "tRNA_pp_valid",
    layer = "data"
  )[, 1]
  data.frame(
    hbb = hbb_expr,
    alas2 = alas2_expr,
    x_marker = x_expr,
    tRNA = tRNA_expr,
    pseudotime = so_cd34e$palantir_pseudotime,
    feature = trna
  )
})
df_trend_long <- df_trend_long %>%
  rownames_to_column("cell")
df_trend_summary <- df_trend_long %>%
  group_by(feature) %>%
  reframe(tidy(cor.test(tRNA, pseudotime))) %>%
  select(feature, R = estimate, p = p.value) %>% 
  mutate(p_adj = p.adjust(p, method = "BH"))

ggplot(df_trend_summary, aes(x = R, y = -1*log10(p_adj))) +
  geom_point() +
  theme_minimal() + labs(title = "Correlation over pseudotime, mono") +
  pretty_plot() + L_border() +
  geom_text_repel(data = df_trend_summary,
                  aes(label = feature), size = 3, max.overlaps = 15, box.padding = 0.3, point.padding = 0.2)

# look at some hits
selected <- c("PP3-Ala-GCC", "PP177-Arg-AGA")
minmax <- function(x) (x - min(x, na.rm = TRUE)) /
  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
df_sel <- df_trend_long %>% filter(feature %in% selected)
# common x grid to predict on
grid <- data.frame(pseudotime = seq(min(df_sel$pseudotime, na.rm = TRUE),
                                    max(df_sel$pseudotime, na.rm = TRUE),
                                    length.out = 200))

# one fitted, min-max-scaled tRNA curve per feature
fit_trna <- df_sel %>%
  group_by(feature) %>%
  group_modify(~{
    m <- gam(tRNA ~ s(pseudotime, bs = "cs"), data = .x)
    data.frame(pseudotime = grid$pseudotime,
               y = minmax(predict(m, newdata = grid)))
  }) %>% ungroup()

# mRNA marker fit once (identical across blocks)
alas2_block <- df_sel %>% filter(feature == selected[1])
m_alas2 <- gam(alas2 ~ s(pseudotime, bs = "cs"), data = alas2_block)
fit_alas2 <- data.frame(pseudotime = grid$pseudotime,
                       y = minmax(predict(m_alas2, newdata = grid)),
                       feature = "ALAS2")
hbb_block <- df_sel %>% filter(feature == selected[1])
m_hbb <- gam(hbb ~ s(pseudotime, bs = "cs"), data = hbb_block)
fit_hbb <- data.frame(pseudotime = grid$pseudotime,
                      y = minmax(predict(m_hbb, newdata = grid)),
                      feature = "HBB")

p <- ggplot(mapping = aes(x = pseudotime, y = y, color = feature)) +
  geom_line(data = fit_trna, linewidth = 1) +
  geom_line(data = fit_alas2, linewidth = 1) +
  geom_line(data = fit_hbb, linewidth = 1) +
  geom_line(data = fit_x, linewidth = 1) +
  scale_color_manual(values = c(
    "PP177-Arg-AGA"  = "#1b9e77",
    "PP3-Ala-GCC"  = "green3",
    "HBB"          = "red",
    "ALAS2"          = "violetred4"
  )) +
  labs(y = "Scaled expression", x = "Pseudotime", color = NULL) +
  theme_minimal() + pretty_plot() + theme(strip.text = element_text(size = 8))
cowplot::ggsave2(p, file = "../output/trajectory_ery_mRNA-tRNA2.pdf", width = 5, height = 3.4)
