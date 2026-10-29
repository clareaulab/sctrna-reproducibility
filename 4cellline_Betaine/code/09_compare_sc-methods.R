library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(BuenColors)
library(Seurat)
library(anndataR)

source('../../global_functions/tRNA_import.R')

#so_4clb <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline_Betaine_annotated.rds")
so_4clb2 <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline_Betaine_deepSeq_annotated.rds")
col1 = c("HeLa"="darkseagreen", "HEK293T"="darkblue", "Jurkat"="maroon", "K562"="cyan3", "VASA" = "purple2", "TotalX" = "darkgrey")

# Import totalX
#totalX <- read_h5ad("../../../trna_large_data/public/TotalX_HEK293T.h5ad", as = "Seurat")
totalX <- readRDS("../../../trna_large_data/public/TotalX_HEK293T.h5ad.asSeurat.rds")
totalX_counts <- (totalX@assays$RNA@layers$X)
tx_anno_df <- totalX@assays$RNA@meta.data
tx_anno_df_trna <- tx_anno_df %>% dplyr::mutate(idx = 1:n()) %>% filter(biotype == "tRNA")
tx_trna_counts_mat <- totalX_counts[tx_anno_df_trna %>% pull(idx),]
rownames(tx_trna_counts_mat) <- make.unique(as.character(tx_anno_df_trna$GeneName))
tx_tRNA_umi_count <- colSums(data.matrix(tx_trna_counts_mat))

# Import VASA
possible <- list.files("../../../trna_large_data/public/GSE176588_RAW/", pattern = "total.UFICounts.tsv.gz")[3:47] # remove first two plate-based measures

lapply(possible, function(i){
  dt <- fread(paste0("../../../trna_large_data/public/GSE176588_RAW/", i))
  names <- dt[[1]]
  counts <- data.matrix(data.frame(dt[,-1]))
  data.frame(
    library = gsub("_total.UFICounts.tsv.gz", "", i),
    what = "VASA",
    tRNA_UMIs = colSums(counts[grep("tRNA", names),])
  )
}) %>% rbindlist() -> vasa_df


# put altogether
final_df <- data.frame(
  what = c(so_4clb2@meta.data$cell_line, rep("TotalX", length(tx_tRNA_umi_count)), vasa_df$what), 
  how_many = c(so_4clb2@meta.data$nCount_tRNA_pp, tx_tRNA_umi_count, vasa_df$tRNA_UMIs)
)

final_df %>% group_by(what) %>% 
  dplyr::summarize(median(how_many), mean(how_many), n())
p <- ggplot(final_df, aes(x = what, y = how_many+1, color = what)) + 
  geom_violin(fill = NA) + geom_boxplot(outlier.shape = NA, width = 0.6) + scale_y_log10() + scale_color_manual(values = col1) +
  pretty_plot(fontsize = 11) + L_border() + labs(x = "", y = "# tRNA UMIs (log10 scaled)")
#cowplot::ggsave2(p, file = "../output/violin_scmethods_compare.pdf", width = 4.5, height = 3, dpi = 400)

ggplot(final_df, aes(x = what, y = how_many+1, color = what)) + 
  geom_violin(fill = NA) + geom_boxplot(outlier.shape = NA, width = 0.6) + scale_color_manual(values = col1) +
  pretty_plot(fontsize = 13.5) + L_border() + labs(x = "", y = "#tRNA UMIs") + coord_cartesian(ylim=c(0,24000))
  #ylim(0,10000)
