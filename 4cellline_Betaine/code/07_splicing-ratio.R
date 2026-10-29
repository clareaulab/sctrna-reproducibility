library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(tibble)
library(BuenColors)
library(ggrepel)

so_4cl <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline_annotated.rds")
so_4clb <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline_Betaine_annotated.rds")

col1 = c("HeLa"="darkseagreen", "HEK293T"="darkblue", "Jurkat"="maroon", "K562"="cyan3")
col2 = c("no"="darkblue", "yes"="cyan3")

exon_ratio <- function(so){
  DefaultAssay(so) <- "tRNA_pp_valid"
  df <- GetAssayData(so, assay = "tRNA_pp_valid")
  df <- data.frame(as.matrix(df))
  
  E1 <- grep("^PP173-Tyr-UAC", rownames(so), value = TRUE)
  I1 <- grep(c("^PP159-Tyr-UAC|^PP153-Tyr-UAC|^PP156-Tyr-UAC|^PP154-Tyr-UAC|^PP155-Tyr-UAC|^PP151|^PP157-Tyr-UAC"), rownames(so), value = TRUE)
  
  E2 <- grep("^PP174-Tyr-UAC", rownames(so), value = TRUE)
  I2 <- grep("^PP158-Tyr-UAC|^PP152-Tyr-UAC|^PP150-Tyr-UAC|^PP160^Tyr-UAC", rownames(so), value = TRUE)
  
  E3 <- grep("^PP175-Ile-AUA", rownames(so), value = TRUE)
  I3 <- grep("^PP89-Ile-AUA", rownames(so), value = TRUE)
  
  E4 <- grep("^PP176-Ile-AUA", rownames(so), value = TRUE)
  I4 <- grep("^PP85-Ile-AUA|^PP86-Ile-AUA|^PP87-Ile-AUA|^PP88-Ile-AUA", rownames(so), value = TRUE)
  
  E5 <- grep("^PP177-Arg-AGA", rownames(so), value = TRUE)
  I5 <- grep("^PP26-Arg-AGA|^PP31-Arg-AGA", rownames(so), value = TRUE)
  
  E6 <- grep("^PP178-Arg-AGA", rownames(so), value = TRUE)
  I6 <- grep("^PP27-Arg-AGA|^PP28-Arg-AGA|^PP29-Arg-AGA", rownames(so), value = TRUE)
  
  E7 <- grep("^PP179-Leu-UUA", rownames(so), value = TRUE)
  I7 <- grep("^PP92-Leu-UUA|^PP93-Leu-UUA", rownames(so), value = TRUE)
  
  E8 <- grep("^PP180-Leu-UUA", rownames(so), value = TRUE)
  I8 <- grep("^PP94-Leu-UUA|^PP95-Leu-UUA|^PP96-Leu-UUA", rownames(so), value = TRUE)

  introns <- c(I1,I2,I3,I4,I5,I6,I7,I8)
  exons <- c(E1,E2,E3,E4,E5,E6,E7,E8)
  
  introns <- gsub("-", ".", introns)
  exons <- gsub("-", ".", exons)
  
  dft <- t(df)
  dft <- data.frame(as.matrix(dft))
  
  dft$intron_count <- rowSums(dft[,introns, drop = FALSE])
  dft$exon_count <- rowSums(dft[,exons, drop = FALSE])
  dft$exon_ratio <- dft$exon_count/(dft$intron_count+dft$exon_count)
  
  dft
}

# subset cell lines
so_4cl_jurkat <- subset(so_4cl, subset = cell_line == "Jurkat")
so_4cl_hela <- subset(so_4cl, subset = cell_line == "HeLa")
so_4cl_hek <- subset(so_4cl, subset = cell_line == "HEK293T")
so_4cl_k562 <- subset(so_4cl, subset = cell_line == "K562")

# now split this up for different cell lines
er_hek <- exon_ratio(so_4cl_hek)
er_hela <- exon_ratio(so_4cl_hela)
er_jurkat <- exon_ratio(so_4cl_jurkat)
er_k562 <- exon_ratio(so_4cl_k562)

er_hek$cell_line <- "HEK293T"
er_hela$cell_line <- "HeLa"
er_jurkat$cell_line <- "Jurkat"
er_k562$cell_line <- "K562"

ratio_4cl <- rbind(er_hek, er_hela, er_jurkat, er_k562)

ggplot(ratio_4cl, aes(x = cell_line, y = exon_ratio, fill = cell_line)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Cell Line", y = "Counts exons/introns+exons") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=col1)
#ggsave(filename = "../output/exon_ratio_4cl.png", plot = p1, width = 5, height = 3, units = "in")


# now add Betaine data
# subset cell lines
so_4clb_jurkat <- subset(so_4clb, subset = cell_line == "Jurkat")
so_4clb_hela <- subset(so_4clb, subset = cell_line == "HeLa")
so_4clb_hek <- subset(so_4clb, subset = cell_line == "HEK293T")
so_4clb_k562 <- subset(so_4clb, subset = cell_line == "K562")

# now split this up for different cell lines
er_hek_b <- exon_ratio(so_4clb_hek)
er_hela_b <- exon_ratio(so_4clb_hela)
er_jurkat_b <- exon_ratio(so_4clb_jurkat)
er_k562_b <- exon_ratio(so_4clb_k562)

er_hek_b$cell_line <- "HEK293T"
er_hela_b$cell_line <- "HeLa"
er_jurkat_b$cell_line <- "Jurkat"
er_k562_b$cell_line <- "K562"

ratio_4clb <- rbind(er_hek_b, er_hela_b, er_jurkat_b, er_k562_b)

ratio_4clb <- ratio_4clb %>%
  dplyr::select(intron_count, exon_count, exon_ratio, cell_line) %>%
  mutate(Betaine = "yes")

ratio_4cl <- ratio_4cl %>%
  dplyr::select(intron_count, exon_count, exon_ratio, cell_line) %>%
  mutate(Betaine = "no")

ratio_comb <- rbind(ratio_4cl, ratio_4clb)

ggplot(ratio_comb, aes(x = cell_line, y = exon_ratio, fill = Betaine)) +
  geom_boxplot() +
  scale_fill_manual(values = col2) +
  theme_void() +
  labs(x = "Cell Line", y = "Counts exon/intron+exon") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(filename = "../output/exon_ratio_4cl_vs_4clb_void.png", plot = p3, width = 6, height = 4, units = "in")

p <- ggplot(ratio_comb, aes(x = cell_line, y = exon_ratio, fill = Betaine)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = col2) +
  pretty_plot() + L_border() +
  labs(x = "Cell Line", y = "Counts exon/intron+exon") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(filename = "../output/exon_ratio_4cl_vs_4clb_noOutliers.pdf", plot = p, width = 6, height = 4, units = "in")
