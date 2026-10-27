library(dplyr)
library(tidyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
filter <- dplyr::filter
select <- dplyr::select

so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_pseudotime_260907.rds")

so_cd34m2 <- subset(so_cd34, subset = condition %in% c("D7_Mono"))
DimPlot(so_cd34m2, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
# clean up
so_cd34m2 <- subset(so_cd34m2, subset = seurat_clusters %in% c(0,1,2,4,11))
DimPlot(so_cd34m2, reduction = "fdl", group.by = c("seurat_clusters"), label=T)
dim(so_cd34m2) #414 19066

so_cd34e2 <- subset(so_cd34, subset = condition %in% c("D7_Ery"))
DimPlot(so_cd34e2, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34e2 <- subset(so_cd34e2, subset = seurat_clusters %in% c(3,5,7,10,12))
DimPlot(so_cd34e2, reduction = "fdl", group.by = c("seurat_clusters"), label=T)

so_cd34k2 <- subset(so_cd34, subset = condition %in% c("D7_Mega"))
DimPlot(so_cd34k2, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34k2 <- subset(so_cd34k2, subset = seurat_clusters %in% c(0,5,8))
DimPlot(so_cd34k2, reduction = "fdl", group.by = c("seurat_clusters"))

# get tRNA expression data
DefaultAssay(so_cd34m2) <- "tRNA_pp_valid"
trna_features <- rownames(so_cd34m2[["tRNA_pp_valid"]])
tRNA_expr <- FetchData(so_cd34m2, vars = trna_features, assay = "tRNA_pp_valid", layer = "data")[, 1]
tRNA_data <- data.frame(so_cd34m2@assays$tRNA_pp_valid$data)
tRNA_time <- data.frame(pseudotime = so_cd34m2$palantir_pseudotime)
tRNA_time <- tRNA_time %>% rownames_to_column("cell")
tRNA_time$cell <- gsub("-", ".", tRNA_time$cell)
tRNA_data <- tRNA_data %>% t() %>% data.frame() %>% rownames_to_column("cell")
tRNA_df <- tRNA_data %>% left_join(tRNA_time %>% dplyr::select(cell, pseudotime), by = c("cell" = "cell"))
tRNA_long_M <- tRNA_df %>% pivot_longer(!c("pseudotime","cell"), names_to = "tRNA_name", values_to = "tRNA_expression")
slopes_M <- tRNA_long_M %>%
  dplyr::group_by(tRNA_name) %>%
  dplyr::summarize(slope = coef(lm(tRNA_expression ~ pseudotime))[2])

DefaultAssay(so_cd34e2) <- "tRNA_pp_valid"
trna_features <- rownames(so_cd34e2[["tRNA_pp_valid"]])
tRNA_expr <- FetchData(so_cd34e2, vars = trna_features, assay = "tRNA_pp_valid", layer = "data")[, 1]
tRNA_data <- data.frame(so_cd34e2@assays$tRNA_pp_valid$data)
tRNA_time <- data.frame(pseudotime = so_cd34e2$palantir_pseudotime)
tRNA_time <- tRNA_time %>% rownames_to_column("cell")
tRNA_time$cell <- gsub("-", ".", tRNA_time$cell)
tRNA_data <- tRNA_data %>% t() %>% data.frame() %>% rownames_to_column("cell")
tRNA_df <- tRNA_data %>% left_join(tRNA_time %>% dplyr::select(cell, pseudotime), by = c("cell" = "cell"))
tRNA_long_E <- tRNA_df %>% pivot_longer(!c("pseudotime","cell"), names_to = "tRNA_name", values_to = "tRNA_expression")
slopes_E <- tRNA_long_E %>%
  dplyr::group_by(tRNA_name) %>%
  dplyr::summarize(slope = coef(lm(tRNA_expression ~ pseudotime))[2])

DefaultAssay(so_cd34k2) <- "tRNA_pp_valid"
trna_features <- rownames(so_cd34k2[["tRNA_pp_valid"]])
tRNA_expr <- FetchData(so_cd34k2, vars = trna_features, assay = "tRNA_pp_valid", layer = "data")[, 1]
tRNA_data <- data.frame(so_cd34k2@assays$tRNA_pp_valid$data)
tRNA_time <- data.frame(pseudotime = so_cd34k2$palantir_pseudotime)
tRNA_time <- tRNA_time %>% rownames_to_column("cell")
tRNA_time$cell <- gsub("-", ".", tRNA_time$cell)
tRNA_data <- tRNA_data %>% t() %>% data.frame() %>% rownames_to_column("cell")
tRNA_df <- tRNA_data %>% left_join(tRNA_time %>% dplyr::select(cell, pseudotime), by = c("cell" = "cell"))
tRNA_long_K <- tRNA_df %>% pivot_longer(!c("pseudotime","cell"), names_to = "tRNA_name", values_to = "tRNA_expression")
slopes_K <- tRNA_long_K %>%
  dplyr::group_by(tRNA_name) %>%
  dplyr::summarize(slope = coef(lm(tRNA_expression ~ pseudotime))[2])

# tRNA variability plots
# mono
tRNA_var <- tRNA_long_M %>%
  dplyr::group_by(tRNA_name) %>%
  dplyr::summarize(
    sd_value = sd(tRNA_expression),
    mean_value = mean(tRNA_expression)
  )
label_df2 <- tRNA_var %>% arrange(desc(sd_value)) %>% slice_head(n = 10)
ggplot(tRNA_var, aes(x = sd_value, y = mean_value)) +
  geom_point(size = 2, color = "grey65") +
  geom_point(data = label_df2, size = 2.1, color = "violetred4") +
  geom_text_repel(data = label_df2, aes(label = tRNA_name),
                  size = 4, box.padding = 0.5, force = 7, force_pull = 0.5, 
                  max.overlaps = Inf, direction = "both", hjust = 1) + 
  coord_flip() + labs(y = "tRNA mean value", x = "SD across pseudotime",
                      title = "Fluctuation of tRNA markers in mono differentiation - no stem") +
  pretty_plot() + L_border() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

label_df3 <- tRNA_var %>%
  filter(tRNA_name %in% gsub("-", ".", grep("Asn", rownames(so_cd34m2), value = TRUE)))
label_df4 <- tRNA_var %>%
  dplyr::filter(tRNA_name %in% c("PP32.Asn.AAC", "PP33.Asn.AAC", "PP40.Asn.AAC"))
p <- ggplot(tRNA_var, aes(x = sd_value, y = mean_value)) +
  geom_point(size = 2, color = "grey65") + geom_point(data = label_df3, size = 2.1, color = "violetred4") + 
  geom_point(data = label_df4, size = 2.1, color = "violetred4") + 
  geom_text_repel(data = label_df4, aes(label = tRNA_name),
                  size = 4, box.padding = 0.5, force = 7, force_pull = 0.5,
                  max.overlaps = Inf, direction = "both",hjust = 1) + 
  coord_flip() +
  labs(y = "tRNA mean value", x = "SD across pseudotime", title = "Fluctuation of tRNA markers in myeloid differentiation"
  ) + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/monocyte_fluctuation.pdf", width = 4.2, height = 4)

# ery
tRNA_var <- tRNA_long_E %>%
  dplyr::group_by(tRNA_name) %>%
  dplyr::summarize(
    sd_value = sd(tRNA_expression),
    mean_value = mean(tRNA_expression)
  )
label_df2 <- tRNA_var %>% arrange(desc(sd_value)) %>% slice_head(n = 5)
p <- ggplot(tRNA_var, aes(x = sd_value, y = mean_value)) +
  geom_point(size = 2, color = "grey65") +
  geom_point(data = label_df2, size = 2.1, color = "violetred4") +
  geom_text_repel(data = label_df2, aes(label = tRNA_name),
                  size = 4, box.padding = 0.5, force = 7, force_pull = 0.5, 
                  max.overlaps = Inf, direction = "both", hjust = 1) + 
  coord_flip() + labs(y = "tRNA mean value", x = "SD across pseudotime",
                      title = "erythroid") +
  pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/ery_fluctuation_noStem2.pdf", width = 4.2, height = 2)

# mega
tRNA_var <- tRNA_long_K %>%
  dplyr::group_by(tRNA_name) %>%
  dplyr::summarize(
    sd_value = sd(tRNA_expression),
    mean_value = mean(tRNA_expression)
  )
label_df2 <- tRNA_var %>% arrange(desc(sd_value)) %>% slice_head(n = 5)
p <- ggplot(tRNA_var, aes(x = sd_value, y = mean_value)) +
  geom_point(size = 2, color = "grey65") +
  geom_point(data = label_df2, size = 2.1, color = "violetred4") +
  geom_text_repel(data = label_df2, aes(label = tRNA_name),
    size = 4, box.padding = 0.5, force = 7, force_pull = 0.5, 
    max.overlaps = Inf, direction = "both", hjust = 1) + 
  coord_flip() + labs(y = "tRNA mean value", x = "SD across pseudotime",
                      title = "megakaryocytic") +
  pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/mega_fluctuation_noStem2.pdf", width = 4.2, height = 2)
