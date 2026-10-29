library(data.table)
library(tidyverse)
library(Matrix)
library(Seurat)
library(BuenColors)

source("../../global_functions/code/plotting_helpers.R")
source("../../global_functions/code/trna_helpers.R")

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA.rds")
DefaultAssay(so_p) = "tRNA_pp_valid"

unique(so_p$identity_NTC)
so_p$target <- str_remove(so_p$identity_NTC, "_guide[1234567890]*$") %>% gsub("tRNA-", "", .)
so_p$guide_nr <- str_remove(so_p$identity_NTC, ".*_")
unique(so_p$target)
stopifnot(all(unique(so_p$target) %in% c("NTC", trna_gene_info$name_short)))

target_probe <- trna_gene_info %>% 
  mutate(in_perturb_panel = name_short %in% unique(so_p$target)) %>%
  filter(in_perturb_panel, !is.na(specific_probe)) %>%
  dplyr::rename(target = name_short, probe = specific_probe) %>%
  dplyr::select(target, probe) 

target_probe_vector <- target_probe$target
names(target_probe_vector) = target_probe$probe
stopifnot(all(target_probe$target %in% so_p$target))
stopifnot(all(target_probe$probe %in% rownames(so_p)))

target_guides <- so_p@meta.data %>%
  as_tibble() %>%
  dplyr::select(identity_NTC, target, guide_nr) %>%
  distinct() %>%
  inner_join(target_probe)
nrow(target_guides)
target_guides <- arrange(target_guides, identity_NTC)

FC <- list()
for (i in 1:nrow(target_guides)){
  this_target = target_guides[[i,"identity_NTC"]]
  this_probe = target_guides[[i,"probe"]]
  print(this_target)
  print(this_probe)
  Idents(so_p) = "identity_NTC"
  FC[[i]] = 
    FindMarkers(so_p, ident.1 = this_target, ident.2 = "NTC", assay = "tRNA_pp_valid", logfc.threshold = 0, min.pct = 0) %>%
    rownames_to_column(var = "probe") %>%
    mutate(target = this_target) %>%
    as_tibble()
}
FC <- do.call(rbind, FC)

probe_labels <- target_guides %>%
  select(probe, target) %>%
  distinct()
FC %>%
  filter(probe %in% target_guides$probe) %>%
  mutate(target2 = stringr::str_remove(target, "_guide\\d+$")) %>%
  dplyr::group_by(probe, target2) %>%
  dplyr::summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
                   p_val = min(p_val, na.rm = TRUE), .groups = "drop") %>%
  mutate(p_val_adj = p.adjust(p_val, method = "BH"),
         probe = factor(probe, levels = unique(target_guides$probe)),
         target2 = factor(target2, levels = rev(unique(
           stringr::str_remove(target_guides$identity_NTC, "_guide\\d+$")
         )))) %>%
  ggplot(aes(x = target2, y = probe, fill = ifelse(p_val_adj < 0.05, avg_log2FC, NA)
  )) +
  geom_tile() +
  scale_fill_gradientn(
    colors = jdb_palette("brewer_yes"),
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    na.value = "white"
  ) +
  scale_x_discrete(labels = target_probe_vector) +
  ylab("measured tRNA") + xlab("perturbed tRNA") + labs(fill = "log2FC (vs. NTC)")

# add names for tRNA loci
FC %>%
  filter(probe %in% target_guides$probe) %>%
  mutate(target2 = stringr::str_remove(target, "_guide\\d+$")) %>%
  dplyr::group_by(probe, target2) %>%
  dplyr::summarise(
    avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
    p_val = min(p_val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_val_adj = p.adjust(p_val, method = "BH")) %>%
  dplyr::left_join(probe_labels, by = "probe") %>%
  mutate(
    probe = factor(probe, levels = unique(target_guides$probe)),
    target2 = factor(
      target2,
      levels = rev(unique(stringr::str_remove(target_guides$identity_NTC, "_guide\\d+$")))
    )
  ) %>%
  ggplot(aes(x = target2, y = probe,
    fill = ifelse(p_val_adj < 0.05, avg_log2FC, NA))) +
  geom_tile() +
  scale_fill_gradientn(
    colors = jdb_palette("brewer_yes"),
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    na.value = "white"
  ) +
  scale_y_discrete(labels = function(x) {
    probe_labels$target[match(x, probe_labels$probe)]
  }) + ylab("measured tRNA") + xlab("perturbed tRNA") + labs(fill = "log2FC (vs. NTC)")

# sort from A to Z
p <- FC %>%
  filter(probe %in% target_guides$probe) %>%
  mutate(target2 = stringr::str_remove(target, "_guide\\d+$")) %>%
  dplyr::group_by(probe, target2) %>%
  dplyr::summarise(
    avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
    p_val = min(p_val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_val_adj = p.adjust(p_val, method = "BH")) %>%
  dplyr::left_join(probe_labels, by = "probe") %>%
  mutate(
    probe = factor(
      probe,
      levels = probe_labels %>%
        filter(probe %in% target_guides$probe) %>%
        arrange(desc(target)) %>%          # desc → A at the top
        pull(probe)
    ),
    target2 = factor(
      target2,
      levels = sort(unique(stringr::str_remove(target_guides$identity_NTC, "_guide\\d+$")))
    )
  ) %>% 
  ggplot(aes(x = target2, y = probe,
             fill = ifelse(p_val_adj < 0.05, avg_log2FC, NA))) +
  geom_tile() +
  scale_fill_gradientn(
    colors = jdb_palette("brewer_yes"),
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    na.value = "white"
  ) +
  scale_y_discrete(labels = function(x) {
    probe_labels$target[match(x, probe_labels$probe)]
  }) + ylab("measured tRNA") + xlab("perturbed tRNA") + labs(fill = "log2FC (vs. NTC)")
#ggsave("../output/heatmap_target_probe_label.pdf", p, height = 5, width = 6.8)

# plot guides individually
FC %>%
  filter(probe %in% target_guides$probe) %>%
  mutate(p_val_adj = p.adjust(p_val, method = "BH")) %>%
  mutate(
    probe = factor(probe, levels = unique(target_guides$probe)),
    target = factor(target, levels = rev(target_guides$identity_NTC))
  ) %>%
  ggplot(aes(
    x = target,
    y = probe,
    fill = ifelse(p_val_adj < 0.05, avg_log2FC, NA)
  )) +
  geom_tile() +
  scale_fill_gradientn(
    colors = jdb_palette("brewer_yes"),
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    na.value = "white"
  ) +
  scale_y_discrete(labels = target_probe_vector) +
  xlab("perturbed tRNA") +
  ylab("measured tRNA") +
  labs(fill = "log2FC (vs. NTC)")
#ggsave("../output/heatmap_by_guide_pvalfiltered.pdf", height = 4.6, width = 11)
