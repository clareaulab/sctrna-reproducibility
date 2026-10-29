library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(tibble)
library(BuenColors)
library(ggrepel)
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

so_4clb <- readRDS(file = "../../../trna_large_data/output/seurat_object_4cellline_Betaine_deepSeq_annotated.rds")
table(so_4clb$cell_line)

# codon demand
#load reference data
codon_demand <- readRDS("../data/human_hg38_codon_usage.rds")

# compute PCAs for all modalities
so_4clb <- NormalizeData(so_4clb, assay = "RNA") 
so_4clb <- FindVariableFeatures(so_4clb, assay = "RNA")
so_4clb <- ScaleData(so_4clb, assay = "RNA")
so_4clb <- RunPCA(so_4clb, assay = "RNA", reduction.name = "pca")

so_4clb <- NormalizeData(so_4clb, assay = "tRNA_pp_valid") 
so_4clb <- FindVariableFeatures(so_4clb, assay = "tRNA_pp_valid")
so_4clb <- ScaleData(so_4clb, assay = "tRNA_pp_valid")
so_4clb <- RunPCA(so_4clb, assay = "tRNA_pp_valid", reduction.name = "tpca")

so_4clb <- FindMultiModalNeighbors(so_4clb, 
                                    reduction.list = list("pca", "tpca"), 
                                    dims.list = list(1:30, 1:10))
so_4clb@graphs
so_4clb <- RunUMAP(so_4clb, nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap", 
                    reduction.key = "wnnUMAP_")
so_4clb@reductions
features_use <- intersect(rownames(codon_demand), rownames(so_4clb@reductions$pca@feature.loadings))
missing_features <- features_use[is.na(match(features_use, dimnames(so_4clb@assays$RNA@features)[[1]]))]
features_use_clean <- features_use[!is.na(match(features_use, dimnames(so_4clb@assays$RNA@features)[[1]]))]
codon_cell_mat <- t(codon_demand[features_use_clean,]) %*% so_4clb@assays$RNA@layers$data[match(features_use_clean, dimnames(so_4clb@assays$RNA@features)[[1]]),]
cd <- codon_cell_mat; colnames(cd) <- dimnames(so_4clb@assays$RNA@cells)[[1]]
cd_mat <- as.matrix(cd)
so_4clb[["codondemand"]] <- CreateAssayObject(counts = cd_mat)


# amino acid demand
cd2 <- data.frame(as.matrix(so_4clb[["codondemand"]]@counts))
cd_rn <- data.frame(rownames(cd2))

# Get reference df to find matching and non-encoded anticodons
ac_map <- read_xlsx("../../global_functions/small_data/tRNA_anticodons.xlsx")
ac_map2 <- ac_map %>% 
  dplyr::select(map_ac = RNA_anticodon, map_aa = AA)
cd_rn$rownames.cd2. <- gsub("T", "U", cd_rn$rownames.cd2.)
rn_aa <- cd_rn %>%
  left_join(ac_map2 %>% dplyr::select(map_ac, map_aa),
            by = c("rownames.cd2." = "map_ac"))
rn_aa$aa_ac <- paste(rn_aa$map_aa, rn_aa$rownames.cd2., sep="-")
cd_aa <- cd2
rownames(cd_aa) <- rn_aa$aa_ac
cd_aa_filt <- cd_aa %>% 
  filter(rownames(cd_aa) != "Stop-UAA")
cd_aa_filt <- cd_aa_filt %>% 
  filter(rownames(cd_aa_filt) != "Stop-UAG")

# aggregate codon demand matrix to AAs
groups <- sub("-.*", "", rownames(cd_aa_filt))
counts_grouped <- rowsum(cd_aa_filt, group = groups)
groups <- sub("-.*", "", rownames(cd_aa_filt))
counts_grouped <- as.matrix(
  t(sapply(split(seq_len(nrow(cd_aa_filt)), groups),
           function(idx) Matrix::colSums(cd_aa_filt[idx, , drop = FALSE])))
)
aa_demand <- data.frame(counts_grouped)
colnames(aa_demand) <- gsub("\\.", "-", colnames(aa_demand))
#aa_demand_mat <- as.matrix(aa_demand)
so_4clb[["aa_demand3"]] <- CreateAssayObject(counts = aa_demand)

cells_keep <- readLines("../data/4cellline_filtered_barcodes_keep.txt")
dim(so_4clb) # 18129 22682
so_4clb_filt <- so_4clb[,cells_keep]
dim(so_4clb_filt) # 18129 18528

# correlate tRNAac and codondemand
aa_repr_tRNA_norm <- function(so){
  mat <- as.matrix(so[["tRNA_aa"]]@counts)
  
  # CPM
  cell_sums <- colSums(mat)
  mat_cpm <- t(t(mat) / cell_sums) * 1e6
  df_norm <- as.data.frame(mat_cpm)
  df_norm$RowMean <- rowMeans(df_norm)
  
  df_summary <- df_norm %>%
    dplyr::select(RowMean) %>%
    tibble::rownames_to_column("AA") %>%
    dplyr::arrange(desc(RowMean)) %>%
    #mutate(rank = seq_len(n()))
    dplyr::mutate(rank = dplyr::row_number()) #%>% 
    #mutate(celltype = ct)
  
  return(df_summary)
}

aa_repr_cd_norm <- function(so){
  mat <- as.matrix(so[["aa_demand3"]]@counts)
  
  # CPM
  cell_sums <- colSums(mat)
  mat_cpm <- t(t(mat) / cell_sums) * 1e6
  df_norm <- as.data.frame(mat_cpm)
  df_norm$RowMean <- rowMeans(df_norm)
  
  df_summary <- df_norm %>%
    dplyr::select(RowMean) %>%
    tibble::rownames_to_column("AA") %>%
    dplyr::arrange(desc(RowMean)) %>%
    dplyr::mutate(rank = dplyr::row_number()) #%>%
    #mutate(celltype = ct)
  
  return(df_summary)
}

celltype_list <- SplitObject(so_4clb_filt, split.by = "cell_line")

# apply over the list and 
lapply(1:length(celltype_list), function(idx){
  aa_repr_cd_norm(celltype_list[[1]], names(celltype_list)[idx]) 
}) %>% bind_rows -> ac_cd_df
lapply(1:length(celltype_list), function(idx){
  aa_repr_tRNA_norm(celltype_list[[1]], names(celltype_list)[idx]) 
}) %>% bind_rows -> ac_tRNA_df

ac_tRNA_results <- purrr::map(celltype_list, ~ aa_repr_tRNA_norm(.x))
ac_cd_results <- purrr:::map(celltype_list, ~ aa_repr_cd_norm(.x))

ac_tRNA_df <- bind_rows(
  lapply(names(ac_tRNA_results), function(ct) {
    df <- ac_tRNA_results[[ct]]
    df$celltype <- ct
    return(df)
  })
)
##ac_tRNA_df$RowMean <- ac_tRNA_df$RowMean+0.0001
##ac_tRNA_df$RowMean <- ac_tRNA_df$RowMean

ac_cd_df <- bind_rows(
  lapply(names(ac_cd_results), function(ct) {
    df <- ac_cd_results[[ct]]
    df$celltype <- ct
    return(df)
  })
)
##ac_cd_df$RowMean <- ac_cd_df$RowMean+0.0001

ac_tRNA_count <- ac_tRNA_df %>% 
  dplyr::select(AA, tRNA_supply = RowMean, rank, celltype)
ac_cd_count <- ac_cd_df %>% 
  dplyr::select(AA, tRNA_demand = RowMean, rank, celltype)

#ac_tRNA_count$AA <- gsub("U", "T", ac_tRNA_count$AC)

ac_tRNA_count$AA_celltype <- paste(ac_tRNA_count$AA, ac_tRNA_count$celltype, sep="-")
ac_cd_count$AA_celltype <- paste(ac_cd_count$AA, ac_cd_count$celltype, sep="-")

ac_count <- ac_cd_count %>% 
  left_join(ac_tRNA_count %>% dplyr::select(AA_celltype, tRNA_supply),
            by = c("AA_celltype" = "AA_celltype"))

# remove SeC bc can't differentiate between stop codon and amino acid in codon demand assay
ac_count <- ac_count %>% 
  filter(AA != "SeC")

ggplot(ac_count, aes(x=log(tRNA_demand), y=log(tRNA_supply))) +
  geom_star(aes(starshape = celltype, fill=celltype), size = 2.5) + pretty_plot() + L_border() + 
  stat_cor(
    aes(label = paste("R ==", signif(after_stat(r), 2), "~`,`~",
                      "p ==", signif(after_stat(p), 2))),
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    parse = TRUE
  )

cor_stats <- ac_count %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    r = cor(log(tRNA_demand), log(tRNA_supply), method = "spearman"), #spearman
    p = cor.test(log(tRNA_demand), log(tRNA_supply), method = "spearman")$p.value,
    .groups = "drop"
  )
cor_stats
cor_stats <- cor_stats %>%
  mutate(label = paste0("R = ", signif(r, 2), ", p = ", signif(p, 2)))
cor_stats <- cor_stats %>%
  dplyr::arrange(celltype) %>%
  dplyr::mutate(y_pos = Inf, 
                x_pos = -Inf,
                offset = row_number())

celltype_cols = c("HeLa"="darkseagreen", "HEK293T"="darkblue", "Jurkat"="violetred4", "K562"="cyan3")

ggplot(ac_count, aes(x = log(tRNA_demand), y = log(tRNA_supply))) +
  geom_star(aes(starshape = celltype, fill = celltype), size = 2.5) +
  pretty_plot() +
  L_border() +
  geom_text(data = cor_stats, aes(x = x_pos, y = y_pos, label = label, color = celltype),
    hjust = -0.1, vjust = 1.2 + cor_stats$offset * 1.2) +
  geom_text(data=ac_count, aes(label=AA)) +
  scale_fill_manual(values = celltype_cols) +
  scale_color_manual(values = celltype_cols)

write.csv(ac_count, "../output/ac_count.csv")
cor_stats2 <- ac_count %>%
  #dplyr::group_by(celltype) %>%
  dplyr::summarise(
    r = cor(log10(tRNA_demand), log10(tRNA_supply), method = "spearman"),
    p = cor.test(log10(tRNA_demand), log10(tRNA_supply), method = "spearman")$p.value,
    .groups = "drop"
  )
cor_stats2
format.pval(
  cor.test(log10(ac_count$tRNA_demand), log10(ac_count$tRNA_supply), method = "spearman")$p.value,
  digits = 22) 
#"< 2.2204460492503130808e-16"
# pearson: 6.41181454572862765475e-12

# this plot
p <- ggplot(ac_count, aes(x = log10(tRNA_demand), y = log10(tRNA_supply), color=celltype)) +
  geom_point(size = 2.5) +
  pretty_plot() +
  L_border() +
  geom_text(data = cor_stats, aes(x = x_pos, y = y_pos, label = label, color = celltype),
            hjust = -0.1, vjust = 1.2 + cor_stats$offset * 1.2) +
  #geom_text(data=ac_count, aes(label=AA)) +
  scale_fill_manual(values = celltype_cols) +
  scale_color_manual(values = celltype_cols)
cowplot::ggsave2(p, file = "../output/codondemand-supply_AA2.pdf", width = 3.5, height = 3, dpi = 400)
cowplot::ggsave2(p, file = "../output/codondemand-supply_AA3.pdf", width = 3.5, height = 3.5, dpi = 400)
