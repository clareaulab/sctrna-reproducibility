library(Seurat)

so_p <- readRDS(file = "../../../trna_large_data/output/seurat_object_Perturb_pilot_sgRNA.rds")
so_p$target <- str_remove(so_p$identity_NTC, "_guide[1234567890]*$") %>% gsub("tRNA-", "", .)

expr <- GetAssayData(so_p, assay = "tRNA_aa", layer = "counts")["SeC", ]
df <- data.frame(SeC = expr, target = so_p$target)
df <- df[df$target %in% c("SeC-TCA-1-1", "NTC"), ]

m_ntc <- mean(df$SeC[df$target == "NTC"])
m_kd  <- mean(df$SeC[df$target == "SeC-TCA-1-1"])

pct_kd <- (1 - m_kd / m_ntc) * 100
pct_kd #44.7197

# cell counts and significance test
table(df$target)
# NTC  SeC-TCA-1-1 
# 1303         343 
wilcox.test(SeC ~ target, data = df) # W = 336776, p-value < 2.2e-16
str(wilcox.test(SeC ~ target, data = df))
