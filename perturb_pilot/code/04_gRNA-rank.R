library(ggplot2)
library(data.table)
library(dplyr)
library(Matrix)
library(matrixStats)
library(Seurat)
filter <- dplyr::filter
select <- dplyr::select

# process_set now pulls second highest ranking gRNA too

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
  
  # assign each barcode to guide
  assign_per_barcode <- function(mat){
    nc <- ncol(mat)
    first_umi  <- rowOrderStats(mat, which = nc) # highest ranking guide
    second_umi <- rowOrderStats(mat, which = nc - 1)  # second highest ranking guide
    
    # second-ranking guide ID: mask the top guide, then take the max again
    top_idx <- max.col(mat, ties.method = "first")
    mat2 <- mat
    mat2[cbind(seq_len(nrow(mat2)), top_idx)] <- -Inf
    
    data.frame(
      barcode    = rownames(mat),
      max_id     = colnames(mat)[top_idx],
      second_id  = colnames(mat)[max.col(mat2, ties.method = "first")],
      n_4        = rowSums(mat >= 4),
      n_10       = rowSums(mat >= 10),
      max_umi    = first_umi,
      second_umi = second_umi,
      total_umi  = rowSums(mat)
    ) %>%
      mutate(assignable = n_4 == 1 | ((max_umi / total_umi) >= 0.75 & max_umi >= 2))
  }
  
  grna_df <- assign_per_barcode(grna_mat) %>% mutate(is_cell = barcode %in% valid_barcodes)
  
  cell_df <- grna_df %>% filter(is_cell) %>% mutate(what = gex_id)
  cell_df
}

trna_df  <- rbind(
  process_set("20250604_CRISPR_72plex_sgRNAoutput", "PerturbPilot.wta", "1")
)

trna_df %>%
  dplyr::group_by(what) %>%
  dplyr::summarize(low_UMI = mean(n_4 < 1)*100, doublet = mean(n_4 > 1)*100, assigned = mean(n_4 == 1)*100, n = n())
trna_df %>% filter(assignable) %>% pull(max_id) %>% table() %>% data.frame() %>%
  write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

plot_df <- trna_df %>%
  select(barcode, what, max_umi, second_umi) %>%
  pivot_longer(c(max_umi, second_umi), names_to = "rank", values_to = "umi") %>%
  mutate(
    rank = factor(rank, levels = c("max_umi", "second_umi"),
                  labels = c("1st", "2nd")),
    log10_umi = log10(umi + 1)
  )

cols <- c("1st" = "#f7941d", "2nd" = "#f7941d")
p <- ggplot(plot_df, aes(x = rank, y = log10_umi, fill = rank)) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = expression(log[10]("UMIs + 1"))) +
  theme_classic() +
  theme(legend.position = "none")
ggsave(filename = "../output/boxplot_gRNA_ranks.pdf", plot = p, width = 3.4, height = 2.7, units = "in")
