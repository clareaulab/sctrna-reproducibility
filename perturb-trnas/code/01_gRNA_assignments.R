library(data.table)
library(dplyr)
library(Matrix)
library(matrixStats)
library(Seurat)
library(BuenColors)

import_kb_grna <- function(i){
  genes <- fread(paste0("../data/grna/",i,".genes.txt"), header = FALSE)[[1]]
  barcodes <- fread(paste0("../data/grna/",i,".barcodes.txt"), header = FALSE)[[1]]
  mtx <- fread(paste0("../data/grna/",i,".mtx"), skip = 3)
  
  # create matrix
  sm <- sparseMatrix(
    i = c(mtx[[1]],length(barcodes)),
    j = c(mtx[[2]],length(genes)),
    x = c(mtx[[3]],0)
  )
  rownames(sm) <- paste0(barcodes, "-1")
  colnames(sm) <- genes
  return(data.matrix(sm))
}

gRNA_raw <- import_kb_grna("20250604_CRISPR_72plex_sgRNAoutput")
mRNA_cr <- Read10X_h5("../../../20250524_CRISPR_72plex_WTA.h5")

# Function to assign each barcode
assign_per_barcode <- function(sparse_mat){
  
  second_largest <- apply(sparse_mat, 1, function(x) sort(x, decreasing = TRUE)[2])
  
  summary_df <- data.frame(
    barcode = rownames(sparse_mat),
    max_id = colnames(sparse_mat)[max.col(sparse_mat)],
    n_4 = rowSums(sparse_mat >= 4),
    n_10 = rowSums(sparse_mat >= 10),
    max_umi = rowMaxs(sparse_mat),
    second_umi = second_largest,
    total_umi = rowSums(sparse_mat)
  ) %>%
    mutate(assignable = n_4 == 1 | ((max_umi/total_umi) >= 0.75 & max_umi >= 2))
}

initial_assign_df <- assign_per_barcode(gRNA_raw) %>% mutate(is_cell = barcode %in% colnames(mRNA_cr))

filtered_df <- initial_assign_df %>%
  filter(assignable & is_cell) 

filtered_df[,c("max_umi", "second_umi")] %>%
  rename("x1" = "max_umi") %>%
  rename("x2" = "second_umi") %>%
  reshape2::melt() -> melt_df

p1 <- melt_df %>%
  ggplot(aes(x = variable, y = (value+1))) + 
  geom_boxplot(outlier.shape = NA) + scale_y_log10() +
  labs(x = "Guide rank", y = "#UMIs (log10 scaled)") +
  pretty_plot(fontsize = 7) + L_border()
cowplot::ggsave2(p1, file = "../plots/grna_out.pdf", width = 1.8, height = 1.8)
melt_df %>% group_by(variable) %>% summarize(mean(value))

initial_assign_df %>%
  filter(is_cell) %>% pull(max_id) %>%
  table() %>% sort() %>% data.frame()
