library(dplyr)
library(Seurat)
library(BuenColors)
library(ggplot2)

all_diff_out <- read.csv("../output/all_diff_out_tRNApp_260901.csv")
all_diff_out2 <- all_diff_out %>% group_by(gRNA) %>% 
  dplyr::summarize(count_down = sum(avg_log2FC < 0.5 & p_val_adj < 0.05), 
                   count_up = sum(avg_log2FC > .500 & p_val_adj < 0.05)) 

ggplot(all_diff_out2, aes(x = count_down, y = count_up, label = gRNA)) + 
  geom_point() + coord_cartesian(xlim = c(0,120), ylim = c(0,120)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_text_repel(data = filter(all_diff_out4, count_down >= 65), aes(label = gRNA),
    size = 3, max.overlaps = 12, box.padding = 0.2, point.padding = 0.1) + 
  pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/n_tRNA_up_down_abline.pdf", width = 3.5, height = 3.5)

# color by total mean tRNA counts per cell
so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_sgRNA_260901.rds")
so_p[["RNA"]] <- JoinLayers(so_p[["RNA"]])

mean_tRNA <- so_p@meta.data %>%
  group_by(identity_NTC) %>%
  dplyr::summarise(
    mean_nCount_tRNA_pp_valid = mean(nCount_tRNA_pp_valid, na.rm = TRUE),
    .groups = "drop"
  )

all_diff_out2_plot <- all_diff_out2 %>%
  left_join(mean_tRNA, by = c("gRNA" = "identity_NTC"))

ggplot(all_diff_out4_plot, 
       aes(x = count_down, y = count_up, color = mean_nCount_tRNA_pp_valid)) + 
  geom_point() + 
  coord_cartesian(xlim = c(0, 125), ylim = c(0, 125)) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", linewidth = 0.5) +
  geom_text_repel(data = filter(all_diff_out4_plot, count_down >= 65), aes(label = gRNA),
    size = 3, max.overlaps = 12, box.padding = 0.2, point.padding = 0.1) +
  scale_color_viridis_c(name = "Mean tRNA counts") +
  pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/n_tRNA_up_down_abline_color_260904.pdf", width = 5, height = 3.5)
