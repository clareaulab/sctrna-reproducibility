library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(BuenColors)
library(Seurat)
source('../../global_functions/tRNA_import.R')

col1 = c("HeLa"="darkseagreen", "HEK293T"="darkblue", "Jurkat"="violetred4", "K562"="cyan3")
col2 = c("standard"="darkblue", "betaine"="cyan3")

# modified import function to do downsampling
import_giftwrap_dataframe <- function(path_1){
  
  # import data
  av_parse <- fread(paste0("../../../trna_large_data/giftwrap/", path_1, "/probe_reads.tsv.gz"))
  av_anno <- fread(paste0("../../../trna_large_data/giftwrap/", path_1, "/manifest.tsv"), header = TRUE)
  av_anno$index <- as.character(av_anno$index)
  av_anno_vec <- av_anno$name; names(av_anno_vec) <- as.character(av_anno$index)
  
  # more annotations
  probe_pair_anno <- fread("../../global_functions/small_data/tRNA_probes_g3_design_human_pairs_250105.csv", header = TRUE)
  aa_vec <- probe_pair_anno$amino_acid; names(aa_vec) <- probe_pair_anno$probePairID
  
  # annotate
  barcodes_for_av <- fread(paste0("../../../trna_large_data/giftwrap/", path_1, "/barcodes.tsv.gz"))[[1]]
  av_parse$barcode <- barcodes_for_av[av_parse$cell_idx + 1]
  av_parse$pp <- av_anno_vec[as.character(av_parse$probe_idx)]
  av_parse$amino_acid <- aa_vec[av_parse$pp]
  av_parse <- av_parse[complete.cases(av_parse) & nchar(av_parse$gapfill) == 3, ]
  av_parse$feature_name <- paste0(av_parse$pp, "_", av_parse$amino_acid, "_", RNATtoU(av_parse$gapfill))
  av_parse$gapfill <- RNATtoU(av_parse$gapfill)
  return(av_parse)
}

# Function to handle downsampling
create_downsample_curves <- function(gwdf, ds_amounts = c(1:10/10)){
  set.seed(1234)
  
  n_cells <- length(unique(gwdf$barcode))
  
  lapply(ds_amounts, function(ds_amount){
    
    # simulate PCR duplicates observed per molecule at this downsample level
    gwdf$obs_count_in_ds <- rbinom(length(gwdf$pcr_duplicate_count),
                                   gwdf$pcr_duplicate_count, ds_amount)
    
    # total reads actually retained at this downsample level
    reads_this_ds <- sum(gwdf$obs_count_in_ds)
    mean_reads_per_cell <- reads_this_ds / n_cells
    
    # estimate UMIs per cell
    gwdf %>% dplyr::group_by(barcode) %>%
      dplyr::summarize(n_UMIs_per_cell = sum(obs_count_in_ds > 0)) %>%
      dplyr::summarize(mean_UMIs = mean(n_UMIs_per_cell),
                median_UMIs = median(n_UMIs_per_cell)) %>%
      mutate(ds_amount = ds_amount,
             mean_reads_per_cell = mean_reads_per_cell)
    
  }) %>% rbindlist() %>% data.frame() -> ds_df
  return(ds_df)
}

# import cells
cells_standard <- Read10X_h5("../../../trna_large_data/mRNA/4CL.wta.h5") %>% colnames()
cells_betaine <- Read10X_h5("../../../trna_large_data/mRNA/4CLB.wta.h5")  %>% colnames()

# import giftwrap probe reads
gwdf_standard <- import_giftwrap_dataframe("4CL")
gwdf_betaine <- import_giftwrap_dataframe("4CLB")

# remove some columns
gwdf_standard_ss <- gwdf_standard[gwdf_standard$barcode %in% cells_standard,c("gapfill", "pcr_duplicate_count", "barcode", "feature_name")]
gwdf_betaine_ss <- gwdf_betaine[gwdf_betaine$barcode %in% cells_betaine,c("gapfill", "pcr_duplicate_count", "barcode", "feature_name")]

# compute DS curves
dsdf_standard <- gwdf_standard_ss %>% create_downsample_curves() %>% mutate(experiment = "standard")
dsdf_betaine <- gwdf_betaine_ss %>% create_downsample_curves() %>% mutate(experiment = "betaine")

# plot
p <- rbind(dsdf_standard, dsdf_betaine) %>%
  ggplot(aes(x = mean_reads_per_cell, y = mean_UMIs, color = experiment)) +
  scale_color_manual(values = col2) +
  geom_point() + geom_line() +
  labs(x = "Mean reads per cell", y = expression(10^3~"mean UMIs"), title = "Read saturation") +
  theme_classic()
#ggsave(filename = "../output/saturation_all_updated_260916.pdf", plot = p, width = 3.1, height = 2.5, units = "in")
