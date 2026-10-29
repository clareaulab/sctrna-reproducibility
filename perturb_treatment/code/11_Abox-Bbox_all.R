library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(tibble)
library(BuenColors)
library(ggrepel)
library(presto)
library(pheatmap)
library(uwot)
library(ggseqlogo)
filter <- dplyr::filter
select <- dplyr::select

so_p2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_tRNAPS_FlexV2_sgRNA.rds")
so_p2$target <- str_remove(so_p2$identity_NTC, "_[AB]$")
so_p2 <- subset(so_p2, subset = target != "TRMT5")
table(so_p2$identity_NTC)

# extract sequences
df <- read.csv("../data/tRNA_sequences_Abox-Bbox-custom_quant.csv")
df2 <- df
abox_seqs <- df2$abox_custom
bbox_seqs <- df2$bbox_custom

# A box sequence logo
p_abox <- ggseqlogo(abox_seqs, method = "bits",  # information content
                    seq_type = "dna") +
  scale_fill_manual(
    values = c("G" = "#51bbd6", "A" = "#2e377c", "T" = "#f8b867", "C" = "violetred4")) + #f79520 #facb7e #f78c1e
  labs(title = "A box", x = "Position", y = "Information content (bits)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )
p_abox
#ggsave(p_abox, filename = "../output/A-box_all.pdf", width = 6.5, height = 4.5, units = "in")

# B box sequence logo
p_bbox <- ggseqlogo(bbox_seqs,
                    method = "bits",
                    seq_type = "dna") +
  scale_fill_manual(
    values = c("G" = "#51bbd6", "A" = "#2e377c", "T" = "#f8b867", "C" = "violetred4")) +
  labs(title = "B box",
       x = "Position", 
       y = "Information content (bits)") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )
p_bbox
#ggsave(p_bbox, filename = "../output/B-box_all.pdf", width = 6.5, height = 4.5, units = "in")
