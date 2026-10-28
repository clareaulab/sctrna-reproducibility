library(data.table)
library(dplyr)
library(Matrix)
library(matrixStats)
library(Seurat)
library(BuenColors)

filter <- dplyr::filter
select <- dplyr::select

process_set <- function(grna_id, gex_id, int){
  import_kb_grna <- function(i){
    genes <- fread(paste0(paste0("../data/grna_counts/",i,".genes.txt")), header = FALSE)[[1]]
    barcodes <- fread(paste0("../data/grna_counts/",i,".barcodes.txt"), header = FALSE)[[1]]
    mtx <- fread(paste0("../data/grna_counts/",i,".mtx"), skip = 3)
    
    # create matrix
    sm <- sparseMatrix(
      i = c(mtx[[1]],length(barcodes)),
      j = c(mtx[[2]],length(genes)),
      x = c(mtx[[3]],0)
    )
    rownames(sm) <- paste0(barcodes, "-", int)
    colnames(sm) <- genes
    return(data.matrix(sm))
  }
  
  grna_mat <- import_kb_grna(grna_id)
  
  valid_barcodes <- paste0(substr(colnames(Read10X_h5(paste0("../../../trna_large_data/mRNA/",gex_id,".h5"))), 1, 16), "-", int)
  
  # Function to assign each barcode
  assign_per_barcode <- function(sparse_mat){
    summary_df <- data.frame(
      barcode = rownames(sparse_mat),
      max_id = colnames(sparse_mat)[max.col(sparse_mat)],
      n_4 = rowSums(sparse_mat >= 4),
      n_10 = rowSums(sparse_mat >= 10),
      max_umi = rowMaxs(sparse_mat),
      total_umi = rowSums(sparse_mat)
    ) %>%
      mutate(assignable = n_4 == 1 | ((max_umi/total_umi) >= 0.75 & max_umi >= 2))
  }
  
  grna_df <- assign_per_barcode(grna_mat) %>% mutate(is_cell = barcode %in% valid_barcodes)
  
  cell_df <- grna_df %>% filter(is_cell) %>% mutate(what = gex_id)
  cell_df
}

trna_df  <- rbind(
  process_set("20251207_FlexV2_gRNA_BC001", "tRNAPS_FlexV2_Control.wta", "Control"),
  process_set("20251207_FlexV2_gRNA_BC002", "tRNAPS_FlexV2_AAdepr1.wta", "AAdepr1"),
  process_set("20251207_FlexV2_gRNA_BC003", "tRNAPS_FlexV2_AAdepr2.wta", "AAdepr2"),
  process_set("20251207_FlexV2_gRNA_BC004", "tRNAPS_FlexV2_Torin1.wta", "Torin1")
)

trna_df %>%
  group_by(what) %>%
  summarize(low_UMI = mean(n_4 < 1)*100, doublet = mean(n_4 > 1)*100, assigned = mean(n_4 == 1)*100, n = n())

trna_df %>% filter(assignable) %>% pull(max_id) %>% table() %>% data.frame() %>%
  write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(trna_df, file = "../output/PS2_tRNA_gRNA_calls.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
