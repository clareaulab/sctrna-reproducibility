library(dplyr)
library(Seurat)
library(BuenColors)
library(ggrepel)
library(Biostrings)
filter <- dplyr::filter
select <- dplyr::select

so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_pseudotime_260907.rds")

so_cd34m <- subset(so_cd34, subset = condition == "D7_Mono")
so_cd34e <- subset(so_cd34, subset = condition == "D7_Ery")

aa_repr <- function(so){
  mat <- so@assays$tRNA_aa$counts
  df <- as.data.frame(t(t(mat) / colSums(mat) * 1e6))
  df$RowSum <- rowSums(df)
  df$RowMean <- rowMeans(df)
  df$RowMedian <- apply(df, 1, median)
  df$RowSD <- apply(df, 1, sd)
  df_summary <- df %>%
    dplyr::select(RowSum, RowMean, RowMedian, RowSD)

  df_summary <- df_summary %>%
    data.frame() %>%
    rownames_to_column("AA") %>% 
    arrange(RowMean) %>%
    dplyr::mutate(rank = 1:n())
  df_summary
}

aa_m <- aa_repr(so_cd34m)
aa_e <- aa_repr(so_cd34e)
aa_m <- aa_m %>% 
  dplyr::rename("counts_Mono"=RowMean, "rank_Mono"=rank)
aa_e <- aa_e %>% 
  dplyr::rename("counts_Ery"=RowMean, "rank_Ery"=rank)
aa_df <- aa_m %>%
  left_join(aa_e %>% dplyr::select(AA, counts_Ery, rank_Ery),
            by = c("AA" = "AA"))

aa_df$AA_diff <- aa_df$counts_Ery - aa_df$counts_Mono
aa_df$AA_fc <- aa_df$counts_Ery/aa_df$counts_Mono
aa_df$AA_fc_rev <- aa_df$counts_Mono/aa_df$counts_Ery
aa_df <- aa_df %>%
  arrange(desc(AA_fc)) %>%
  dplyr::mutate(rank_diff = 1:n())
ggplot(aa_df, aes(x = rank_diff, y = AA_fc)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = AA), size = 4.5, point.padding = 1) +
  labs(x = "AA rank", y = "AA difference Ery/Mono") + pretty_plot() + L_border()

# color by %AA in Hb
# Hb AA sequences
alpha <- AAString("ZVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR")
beta <- AAString("ZVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH")
hb <- paste0(
  as.character(alpha),
  as.character(alpha),
  as.character(beta),
  as.character(beta)
)
aa <- strsplit(hb, "")[[1]]
hb_freq <- as.data.frame(table(aa))
hb_freq$frequency <- hb_freq$Freq / sum(hb_freq$Freq)
aa_map <- data.frame(
  AA = c("iMet",
    "Ala","Arg","Asn","Asp","Cys",
    "Gln","Glu","Gly","His","Ile",
    "Leu","Lys","Met","Phe","Pro",
    "Ser","Thr","Trp","Tyr","Val",
    "SeC"
  ),
  aa = c("Z",
    "A","R","N","D","C",
    "Q","E","G","H","I",
    "L","K","M","F","P",
    "S","T","W","Y","V",
    "U"
  ),
  stringsAsFactors = FALSE
)
hb_freq <- hb_freq %>%
  left_join(aa_map, by = "aa")

aa_df2 <- aa_df %>% 
  left_join(hb_freq %>% select("frequency", "AA"), by = "AA")
p <- ggplot(aa_df2, aes(x = rank_diff, y = log2(AA_fc), color = (frequency)*100)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = AA), size = 4.5, point.padding = 1) +
  labs(x = "AA rank", y = "AA difference Ery/Mono") + pretty_plot() + L_border() + 
  scale_color_gradientn(colors = c("grey",jdb_palette("solar_rojos")[2:9]))
ggsave(filename = "../output/aa_representation_Ery_Mono.pdf", plot = p, width = 4.8, height = 3.9, units = "in")

cor.test(aa_df2$frequency, log2(aa_df2$AA_fc), method = "spearman", use = "pairwise.complete.obs")
