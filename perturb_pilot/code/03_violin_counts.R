library(dplyr)
library(Seurat)

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA.rds")

so_p@meta.data[, c("nCount_RNA_log", "nCount_tRNA_pp_valid_log", "nCount_sgRNA_log")] <- 
  log10(so_p@meta.data[, c("nCount_RNA", "nCount_tRNA_pp_valid", "nCount_sgRNA")] + 1)

so_p$AllCells <- "All"
mean(so_p$nCount_tRNA_pp_valid) #5614.843

col1 <- c(
  "nCount_RNA" = "violetred4",
  "nCount_tRNA_pp_valid" = "darkblue",
  "nCount_sgRNA" = "#f7941d"
)
vln <- VlnPlot(object = so_p, features = c("nCount_RNA","nCount_tRNA_pp_valid","nCount_sgRNA"), 
               group.by = "AllCells", cols = col1, pt.size = 0)
vln <- vln + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
vln

df <- so_p@meta.data %>%
  select(nCount_RNA, nCount_tRNA_pp_valid, nCount_sgRNA) %>%
  pivot_longer(
    everything(),
    names_to = "feature",
    values_to = "count"
  )
medians <- so_p@meta.data %>%
  summarise(
    nCount_RNA = median(nCount_RNA),
    nCount_tRNA_pp_valid = median(nCount_tRNA_pp_valid),
    nCount_sgRNA = median(nCount_sgRNA)
  ) %>%
  pivot_longer(everything(), names_to = "assay", values_to = "median_umi") %>%
  mutate(y = log10(median_umi + 1))
label_df <- df %>%
  group_by(feature) %>%
  summarise(
    median = median(count),
    y = max(count),
    .groups = "drop")
p <- ggplot(df, aes(x = feature, y = count, fill = feature)) +
  geom_violin(trim = FALSE, linewidth = 0.8) +
  geom_text(data = medians, aes(x = assay, y = y, label = scales::comma(round(median_umi))),
            hjust = -0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = col1) +
  scale_y_log10() + theme_classic() +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs(x = NULL, y = "UMIs")
ggsave(filename = "../output/vioplot_counts4.pdf", plot = p, width = 3.4, height = 3.4, units = "in")
