library(dplyr)
library(tidyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
filter <- dplyr::filter
select <- dplyr::select
summarise <- dplyr::summarise

revComp <- function(invec) {
  chartr(old = "AaTtGgCc", new = "UUAACCGG",  stringi::stri_reverse(invec))
}

so_pp <- readRDS(file = "../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA.rds")
unique(so_pp$identity_NTC)
so_pp$target <- str_remove(so_pp$identity_NTC, "_guide[123]$")
unique(so_pp$target)

so_pp$target_codon <- str_split_fixed(so_pp$target, "-", 5)[,3]
unique(so_pp$target_codon)
so_pp$target_aa <- str_split_fixed(so_pp$target, "-", 5)[,2]
unique(so_pp$target_aa)

so_pp$target_codon_rc <- revComp(so_pp$target_codon)
unique(so_pp$target_codon_rc)

codons <- setdiff(unique(so_pp$target_codon_rc), "")
codon_aas <- setdiff(unique(so_pp$target_aa), "")

plots <- lapply(codons, function(cd) {
  # subset cells: codon vs empty (empty is NTC)
  cells_use <- rownames(so_pp@meta.data)[
    so_pp$target_codon_rc %in% c(cd, "")
  ]
  so_sub <- subset(so_pp, cells = cells_use)
  
  # define grouping variable
  so_sub$group <- ifelse(so_sub$target_codon_rc == cd, cd, "empty")
  Idents(so_sub) <- "group"
  
  # violin plot for that codon feature
  VlnPlot(
    so_sub,
    features = cd,
    assay = "tRNA_ac",
    pt.size = 0
  ) + ggtitle(paste("Codon:", cd))
})
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]
plots[[6]]
plots[[7]]
plots[[8]]

log2FC_cutoff <- 0
pval_cutoff <- 0.05

Idents(so_pp) <- "target_codon_rc"
get_volcano <- function(target_codon, assay){
  diff <- FindMarkers(so_pp, ident.1 = target_codon, ident.2 = "", assay = assay, group.by = "target_codon_rc")
  diff$trna <- rownames(diff)
  diff <- diff %>%
    mutate(
      significance = case_when(
        abs(avg_log2FC) >= log2FC_cutoff & p_val_adj < pval_cutoff ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
  p <- ggplot(diff, aes(x = avg_log2FC, y = -1*log10(p_val_adj), color = significance)) +
    geom_point() +
    scale_color_manual(values = c("lightgrey", "darkblue")) +
    theme_minimal() +
    labs(title = paste(target_codon, assay, sep=", ")) +
    geom_text_repel(
      data = filter(diff, significance == "Significant"), 
      aes(label = trna),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.3,
      point.padding = 0.2
    ) + pretty_plot() + L_border()
  p
}

# get volcano for Ile-AUA
p <- get_volcano("AUA", "tRNA_pp_valid")
cowplot::ggsave2(p, file = "../output/volcano_AUA-KD_260929_unfiltered.pdf", width = 3.5, height = 3, dpi = 400)

# get volcano for Leu-UUG
p <- get_volcano("UUG", "tRNA_pp_valid")
cowplot::ggsave2(p, file = "../output/volcano_UUG-KD_unfiltered.pdf", width = 3.5, height = 3, dpi = 400)
