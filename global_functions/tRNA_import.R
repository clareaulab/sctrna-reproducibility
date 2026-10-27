# functions to import mRNA and tRNA data and combine into one seurat object

filter <- dplyr::filter
select <- dplyr::select

# global mini funcitons
RNATtoU <- function(invec) {
  chartr(old = "Tt", new = "UU", invec)
}

# global files
ac_map <- read_xlsx("../../global_functions/small_data/tRNA_anticodons.xlsx")


# import mRNA data and create seurat object
createSeurat <- function(h5_files, sample_names) {
  if (length(h5_files) != length(sample_names)) {
    stop("Number of files and sample names must match.")
  }
  
  seurat_list <- mapply(function(file, sample) {
    # Read 10x h5 file
    data <- Read10X_h5(paste0("../../../trna_large_data/mRNA/", file, ".wta.h5"))
    
    # if multiple assays, select only "Gene Expression"
    if (is.list(data)) data <- data[["Gene Expression"]]
    
    # create Seurat object
    so <- CreateSeuratObject(counts = data)
    
    # add sample identifier to cell barcodes
    colnames(so) <- paste0(substr(colnames(so),1,16), "-", sample)
    
    return(so)
  }, h5_files, sample_names, SIMPLIFY = FALSE)
  
  # merge seurat objects
  so_merged <- Reduce(function(x, y) merge(x, y), seurat_list)
  return(so_merged)
}

# import tRNA data
# the only difference compared to standard import is this: match_df <- match_df %>% filter(!is.na(AA))
import_tRNA <- function(so, paths, sample_names, feature_min = 25){
  
  if (length(paths) != length(sample_names)) {
    stop("Number of files, barcode sets and sample names must match.")
  }
  
  ac_map <- read_xlsx("../../global_functions/small_data/tRNA_anticodons.xlsx")
  
  av_list <- lapply(seq_along(paths), function(i) {
    path_1 <- paths[[i]]
    sample_1 <- sample_names[[i]]
    
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
    
    # add sample identifier to cell barcodes
    av_parse$barcode <- paste0(substr(av_parse$barcode, 1, 16), "-", sample_1)
    
    return(av_parse)
  })
  
  # merge seurat objects
  av_merged <- data.table::rbindlist(av_list, use.names = TRUE, fill = TRUE)
  
  valid_barcodes <- colnames(so)
  
  # filter probe pair features
  for_pp <- av_merged %>% 
    group_by(feature_name) %>% 
    dplyr::mutate(total_feature = n()) %>%
    dplyr::filter(total_feature >= feature_min) %>%
    dplyr::filter(barcode %in% valid_barcodes) %>%
    group_by(barcode, feature_name) %>%
    dplyr::summarize( x = n())
  df_pp <- sparseMatrix(i = as.factor(for_pp$feature_name), j = as.factor(for_pp$barcode), x = for_pp$x)
  rownames(df_pp) <- levels(as.factor(for_pp$feature_name))
  colnames(df_pp) <- levels(as.factor(for_pp$barcode))
  
  # ac level: filter valid tRNA data
  match_df <- av_merged %>%
    left_join(ac_map %>% dplyr::select(RNA_anticodon, AA),
              by = c("gapfill" = "RNA_anticodon"))
  # if the amino acids of the observed and the expected/encoded anticodon match,
  # the observed anticodon counts as allowed -> capture non-encoded anticodons
  dim(match_df) #82635932  12
  match_df <- match_df %>%
    filter(!is.na(AA))
  dim(match_df) #82634920  12
  match_df <- match_df %>%
    mutate(
      allowed_ac = str_detect(amino_acid, AA),
      allowed_ac_match = ifelse(allowed_ac, gapfill, NA_character_)
    )
  valid_df <- match_df %>% 
    dplyr::filter(!is.na(allowed_ac_match)) %>% 
    dplyr::select(-allowed_ac, -AA)
  # filter features that are rare
  for_ac <- valid_df %>% 
    #group_by(gapfill) %>% 
    #mutate(total_feature = n()) %>%
    #filter(total_feature >= feature_min) %>%
    dplyr::filter(barcode %in% valid_barcodes) %>%
    group_by(barcode, gapfill) %>%
    dplyr::summarize(x = n())
  df_ac <- sparseMatrix(i = as.factor(for_ac$gapfill), j = as.factor(for_ac$barcode), x = for_ac$x)
  rownames(df_ac) <- levels(as.factor(for_ac$gapfill))
  colnames(df_ac) <- levels(as.factor(for_ac$barcode))
  
  # filter AA level
  for_aa <- valid_df %>% 
    #group_by(amino_acid) %>% 
    #mutate(total_feature = n()) %>%
    #filter(total_feature >= 25) %>%
    dplyr::filter(barcode %in% valid_barcodes) %>%
    group_by(barcode, amino_acid) %>%
    dplyr::summarize(x = n())
  df_aa <- sparseMatrix(i = as.factor(for_aa$amino_acid), j = as.factor(for_aa$barcode), x = for_aa$x)
  rownames(df_aa) <- levels(as.factor(for_aa$amino_acid))
  colnames(df_aa) <- levels(as.factor(for_aa$barcode))
  
  return(list(df_pp, df_ac, df_aa))
  
}

# add tRNA data to seurat object
append_tRNA <- function(so, tRNA_dfs){
  df_pp <- tRNA_dfs[[1]]
  df_ac <- tRNA_dfs[[2]]
  df_aa <- tRNA_dfs[[3]]

  df_pp <- df_pp[,colSums(df_pp) >= 50] # eventually make more stringent
  cells_pp <- intersect(colnames(df_pp), colnames(so))
  so <- so[,cells_pp]
  so[["tRNA_pp"]] <- CreateAssayObject(df_pp[,cells_pp])
  
  df_ac <- df_ac[,colSums(df_ac) >= 50] # eventually make more stringent
  cells_ac <- intersect(colnames(df_ac), colnames(so))
  so <- so[,cells_ac]
  so[["tRNA_ac"]] <- CreateAssayObject(df_ac[,cells_ac])
  
  df_aa <- df_aa[,colSums(df_aa) >= 50] # eventually make more stringent
  cells_aa <- intersect(colnames(df_aa), colnames(so))
  so <- so[,cells_aa]
  so[["tRNA_aa"]] <- CreateAssayObject(df_aa[,cells_aa])
  
  return(so)

}

# example, setwd = git_repo/a549/code
#so_a549 <- createSeurat(c("A549_Dex_3hr", "A549_EtOH_3hr"), c("dex", "etoh"))
#trna_raw_data <- import_tRNA(so_a549, c("A549_Dex", "A549_EtOH"), c("dex", "etoh"))
#so_a549 <- append_tRNA(so_a549, trna_raw_data)
