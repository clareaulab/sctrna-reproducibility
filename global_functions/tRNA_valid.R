library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(stringr)
library(readxl)
library(tibble)
library(ggrepel)
library(coRdon)
library(biomaRt)
library(entropy)
library(lme4)
library(reshape2)

filter <- dplyr::filter
select <- dplyr::select

RNATtoU <- function(invec) {
  chartr(old = "Tt", new = "UU", invec)
}

tRNA_valid <- function(so){
  df <- so@assays$tRNA_pp$counts
  df <- as.data.frame(as.matrix(df))
  df$RowSum <- rowSums(df)
  df$RowMean <- rowMeans(df)
  df$RowMedian <- apply(df, 1, median)
  df$RowSD <- apply(df, 1, sd)
  df_summary <- df %>%
    dplyr::select(RowSum, RowMean, RowMedian, RowSD)
  df_summary <- df_summary %>%
    tibble::rownames_to_column("row_names") %>%
    separate(row_names, into = c("pp", "aa", "ac"), sep = "-")
  rownames(df_summary) <- rownames(df_summary)
  df_summary$ac <- RNATtoU(df_summary$ac)
  
  # Get reference df to find matching and non-encoded anticodons
  ac_map <- read_xlsx("../../global_functions/small_data/tRNA_anticodons.xlsx")
  ac_map2 <- ac_map %>% 
    dplyr::select(map_ac = RNA_anticodon, map_aa = AA)
  
  match_df <- df_summary %>%
    left_join(ac_map2 %>% dplyr::select(map_ac, map_aa),
              by = c("ac" = "map_ac"))
  match_df <- match_df %>% filter(!is.na(map_aa))
  
  # if the amino acids of the observed and the expected/encoded anticodon match,
  # the observed anticodon counts as allowed -> capture non-encoded anticodons
  match_df <- match_df %>%
    mutate(
      allowed_ac = str_detect(aa, map_aa),
      allowed_ac_match = ifelse(allowed_ac, ac, NA_character_)
    )
  
  valid_df <- match_df %>% 
    dplyr::filter(!is.na(allowed_ac_match)) %>% 
    dplyr::select(-allowed_ac, -map_aa)
  
  valid_df <- valid_df %>% 
    left_join(ac_map %>% dplyr::select(RNA_anticodon, encoded = in_gtRNAdb),
              by = c("allowed_ac_match" = "RNA_anticodon"))
  
  valid_df_summary <- valid_df %>%
    dplyr::group_by(ac, encoded) %>%
    dplyr::summarise(sum_counts = sum(as.numeric(RowSum)))
  
  valid_ac_barplot <- valid_df_summary %>% 
    #arrange(desc(amino_acid)) %>%
    ggplot(aes(x = ac, y = sum_counts, fill = encoded)) + 
    geom_bar(stat="identity", position = "dodge", color = "black") +
    #scale_y_log10(limits = c(1,1e6)) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(0.5, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"))
  valid_ac_barplot
  ggsave(filename = "../output/barplot_anticodons.png", plot = valid_ac_barplot, width = 20, height = 3.6, units = "in")
  
  valid_ac_plot <- valid_df_summary %>%  
    arrange(desc(sum_counts)) %>%
    ggplot(aes(x = ac, y = sum_counts, color = encoded)) + 
    geom_point() + #scale_y_log10(limits = c(1,1e6)) +
    theme_bw()
  valid_ac_plot
  ggsave(filename = "../output/valid_ac_plot.png", plot = valid_ac_plot, width = 20, height = 3.6, units = "in")
  
  
  valid_df_summary <- valid_df_summary %>%
    data.frame() %>% 
    arrange(desc(sum_counts)) %>%
    dplyr::mutate(rank = 1:n())
  valid_ac_plot2 <- ggplot(valid_df_summary, aes(x = rank, y = log10(sum_counts), color = encoded)) +
    geom_point(size = 1.75) +
    geom_text_repel(aes(label = ac), size = 4) +
    labs(x = "anticodon rank", y = "log10 raw counts") +
    #theme_minimal()
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black")
    )
  valid_ac_plot2
  ggsave(filename = "../output/valid_ac_plot2.png", plot = valid_ac_plot2, width = 20, height = 3.6, units = "in")
  
  missing_ac <- setdiff(ac_map$RNA_anticodon, valid_df_summary$ac) # TTA, CTA
  # two missing codons make sense -> stop codons: UAA, UAG
  print(paste0("Missing anticodons: ", missing_ac, " Expected: UAA and UAG"))
  
  # Find rows to filter for in single-cell tRNA assay
  valid_df$valid_pp <- paste(valid_df$pp, valid_df$aa, valid_df$ac, sep="-")
  valid_pp <- valid_df$valid_pp
  
  df_valid <- df[rownames(df) %in% valid_pp, ]
  df_valid <- df_valid[, !colnames(df_valid) %in% c("RowSum", "RowMean", "RowMedian", "RowSD")]
  
  return(df_valid)
}
