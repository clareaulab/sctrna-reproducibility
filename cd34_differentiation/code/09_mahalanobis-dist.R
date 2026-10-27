library(ggplot2)
library(dplyr)
library(Matrix)
filter <- dplyr::filter
select <- dplyr::select

# import data
so_cd34 <- readRDS(file = "../../../trna_large_data/output/seurat_object_CD34_diff_pseudotime_260907.rds")
FeaturePlot(so_cd34, features = c("palantir_pseudotime"), cols = c("lightgrey", "darkblue"), reduction="fdl")
DimPlot(so_cd34, group.by = c("condition"), reduction="fdl")

# subset to stem/ery d7
so_cd34e <- subset(so_cd34, subset = condition %in% c("D0", "D7_Ery"))
DimPlot(so_cd34e, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34e <- subset(so_cd34e, subset = seurat_clusters %in% c(0,3,5,7,10,12))
DimPlot(so_cd34e, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)

# stem/mega d7
so_cd34k <- subset(so_cd34, subset = condition %in% c("D0", "D7_Mega"))
DimPlot(so_cd34k, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34k <- subset(so_cd34k, subset = seurat_clusters %in% c(0,5,8))
DimPlot(so_cd34k, reduction = "fdl", group.by = c("seurat_clusters"))

# stem/mono d7
so_cd34m <- subset(so_cd34, subset = condition %in% c("D0", "D7_Mono"))
DimPlot(so_cd34m, reduction = "fdl", group.by = c("condition", "seurat_clusters"), label=T)
so_cd34m <- subset(so_cd34m, subset = seurat_clusters %in% c(0,1,2,4,11))
DimPlot(so_cd34m, reduction = "fdl", group.by = c("seurat_clusters"))


compute_pseudotime_simliarities_SE2 <- function(so, n_bins = 100, n_tRNAs = 100,
                                                pt_col = "palantir_pseudotime",
                                                n_cells_extreme = 200,
                                                n_pcs = 30, n_pcs_trna = 15, PERMUTE = FALSE) {
  
  # build PCA embeddings
  so <- so %>% FindVariableFeatures(assay = "tRNA_pp_valid", nfeatures = n_tRNAs) %>%
    RunPCA(npcs = n_pcs, reduction.name = "pca_trna", assay = "tRNA_pp_valid")
  so <- so %>% FindVariableFeatures(assay = "RNA") %>%
    RunPCA(npcs = n_pcs, reduction.name = "pca_rna")
  
  # eigenvalue-weighted/Mahalanobis-style scaling
  scale_embedding <- function(emb, stdev, n_pcs = NULL) {
    if (!is.null(n_pcs)) {
      emb   <- emb[, seq_len(n_pcs), drop = FALSE]
      stdev <- stdev[seq_len(n_pcs)]
    }
    centered <- scale(emb, center = TRUE, scale = FALSE)
    sweep(centered, 2, stdev, "/")
  }
  
  codon_demand <- t(scale_embedding(so@reductions$pca_rna@cell.embeddings,
                                    so@reductions$pca_rna@stdev,  n_pcs))
  codon_supply <- t(scale_embedding(so@reductions$pca_trna@cell.embeddings,
                                    so@reductions$pca_trna@stdev, n_pcs_trna))
  colnames(codon_demand) <- rownames(so@meta.data)
  colnames(codon_supply) <- rownames(so@meta.data)
  
  ordered_cells <- so@meta.data %>% mutate(barcode = rownames(so@meta.data)) %>%
    arrange(palantir_pseudotime) %>% dplyr::mutate(palantir_pseudotime_rank = 1:n())
  n_cells_extreme <- round(dim(ordered_cells)[1]*0.05)
  print(n_cells_extreme)
  
  # pull out the extreme early and late cells
  cells_early <- ordered_cells %>% head(n_cells_extreme) %>% pull(barcode)
  cells_late <- ordered_cells %>% tail(n_cells_extreme) %>% pull(barcode)
  
  if (PERMUTE) {
    set.seed(1234)
    ordered_cells$palantir_pseudotime <- sample(ordered_cells$palantir_pseudotime)
    ordered_cells <- ordered_cells %>% arrange(palantir_pseudotime) %>%
      dplyr::mutate(palantir_pseudotime_rank = 1:n())
  }
  
  # see: https://github.com/GreenleafLab/ArchR/blob/5f855d7b7ff3f57eb0d28f312c5ea5d373d27ebd/R/Trajectory.R
  n_bins <- n_bins
  groups_ps <- lapply(1:n_bins, function(idx){
    multiplier <- 1/(n_bins/dim(ordered_cells)[1])
    bcs <- ordered_cells %>% dplyr::filter(palantir_pseudotime_rank >= (idx-1)*multiplier & palantir_pseudotime_rank < (idx*multiplier)) %>% pull(barcode)
    print(length(bcs))
    bcs
  })
  
  # centroid Euclidean distance in the already-scaled PC space
  # Because codon_demand / codon_supply are pre-scaled by per-PC SD, plain
  # Euclidean distance here is the Mahalanobis distance.
  centroid_dist <- function(mat, reference_cells, query_cells) {
    ref_centroid   <- rowMeans(mat[, reference_cells, drop = FALSE])
    query_centroid <- rowMeans(mat[, query_cells,     drop = FALSE])
    sqrt(sum((ref_centroid - query_centroid)^2))
  }
  
  # compute simliarities
  similarity_df <- lapply(1:length(groups_ps), function(idx){
    barcodes <- groups_ps[[idx]]
    main_estimate_df <- data.frame(
      idx = idx,
      mean_pt = ordered_cells %>%
        filter(barcode %in% barcodes) %>%
        pull(palantir_pseudotime) %>%
        mean(),
      dist_early_demand = centroid_dist(codon_demand, cells_early, barcodes),
      dist_early_supply = centroid_dist(codon_supply, cells_early, barcodes),
      dist_late_demand  = centroid_dist(codon_demand, cells_late,  barcodes),
      dist_late_supply  = centroid_dist(codon_supply, cells_late,  barcodes)
    )
    
    # compute bootstrap standard errors
    boostrap_mat <- sapply(1:100, function(seed){
      set.seed(seed)
      sampled_bcs <- sample(barcodes, replace = TRUE)
      data.frame(
        dist_early_demand = centroid_dist(codon_demand, cells_early, sampled_bcs),
        dist_early_supply = centroid_dist(codon_supply, cells_early, sampled_bcs),
        dist_late_demand  = centroid_dist(codon_demand, cells_late,  sampled_bcs),
        dist_late_supply  = centroid_dist(codon_supply, cells_late,  sampled_bcs)
      ) %>% data.matrix()
    })
    
    # save standard errors
    se <- matrixStats::rowSds(boostrap_mat)
    names_columns_basic <- paste0("dist_", c("early_demand", "early_supply", "late_demand", "late_supply"))
    theta_hat <- unlist(main_estimate_df[, names_columns_basic])
    se_df  <- setNames(as.data.frame(t(se)), paste0(names_columns_basic, "_se"))
    
    cbind(
      main_estimate_df,
      se_df
    )
    
  }) %>% rbindlist() %>% data.frame()
  
  # plot
  pf <- function(feature, z = qnorm(0.975)) {
    ggplot(similarity_df, aes(x = idx)) +
      geom_ribbon(
        aes(ymin = .data[[feature]] - z * .data[[paste0(feature, "_se")]],
            ymax = .data[[feature]] + z * .data[[paste0(feature, "_se")]]),
        fill = "grey40", alpha = 0.25
      ) +
      geom_line(aes(y = .data[[feature]])) +
      labs(title = feature, y = feature) +
      theme_bw()
  }
  cowplot::plot_grid(pf("dist_early_demand"), pf("dist_early_supply"),
                     pf("dist_late_demand"),  pf("dist_late_supply"))
  
  # try smoothing things
  pf_smooth <- function(feature, span = 0.8, z = qnorm(0.975)) {
    y  <- similarity_df[[feature]]
    se <- similarity_df[[paste0(feature, "_se")]]
    
    k  <- max(3, 2 * floor(span * length(y) / 2) + 1)   # odd window
    ys <- runmed(y,  k, endrule = "median")
    ss <- runmed(se, k, endrule = "median")
    
    sm <- data.frame(idx = similarity_df$idx, y = ys,
                     ymin = ys - z * ss, ymax = ys + z * ss)
    
    ggplot(sm, aes(x = idx, y = y)) +
      geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey40", alpha = 0.25) +
      geom_line(linewidth = 0.7) +
      labs(title = feature, y = feature) +
      theme_bw()
  }
  
  cowplot::plot_grid(pf_smooth("dist_early_demand"), pf_smooth("dist_early_supply"),
                     pf_smooth("dist_late_demand"),  pf_smooth("dist_late_supply"))
  
  #z value for a 95% confidence interval
  z <- qnorm(0.975)
  
  #z value for a 99% confidence interval
  #z <- qnorm(0.995)
  
  norm_diff <- function(df, metric, label, span = 0.5) {
    early <- paste0("dist_early_", metric)
    late  <- paste0("dist_late_",  metric)
    
    y  <- df[[early]] - df[[late]]
    se <- sqrt(df[[paste0(early, "_se")]]^2 + df[[paste0(late, "_se")]]^2)  # conservative
    
    k  <- max(3, 2 * floor(span * length(y) / 2) + 1)
    ys <- runmed(y,  k, endrule = "median")
    ss <- runmed(se, k, endrule = "median")
    
    m <- min(ys); M <- max(ys)
    rs <- function(v) (v - m) / (M - m)
    
    data.frame(idx = df$idx, mean_pt = df$mean_pt,
               values = rs(ys), ymin = rs(ys - z * ss), ymax = rs(ys + z * ss),
               what = label)
  }
  
  plot_df <- rbind(
    norm_diff(similarity_df, "supply", "tRNA"),
    norm_diff(similarity_df, "demand", "mRNA")
  )
  
  ggplot(plot_df, aes(x = idx, y = values, color = what, fill = what)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    labs(x = "pseudotime", y = "normalized (early - late)",
         color = NULL, fill = NULL) +
    theme_bw() + pretty_plot() + L_border()

  plot_df
}


# helper function to scale pseudotime from 0 to 1
rng <- range(so_cd34e$palantir_pseudotime)

plot_df_ex <- compute_pseudotime_simliarities_SE2(so_cd34e)
plot_df_ex$pt_rel <- (plot_df_ex$mean_pt - rng[1]) / diff(rng)
ggplot(plot_df_ex, aes(x = pt_rel, y = values, color = what, fill = what)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"), #6E91CB #violetred3
                     labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  scale_fill_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                    labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  labs(x = "pseudotime", y = "normalized (early - late)", color = NULL, fill = NULL) +
  theme_bw() + pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/mahalanobis-dist_ery.pdf", width = 4.5, height = 4, dpi = 400)
#cowplot::ggsave2(p, file = "../output/mahalanobis-dist_ery2.pdf", width = 2.7, height = 1.5, dpi = 400)

plot_df_mx <- compute_pseudotime_simliarities_SE2(so_cd34m)
plot_df_mx$pt_rel <- (plot_df_mx$mean_pt - rng[1]) / diff(rng)
ggplot(plot_df_mx, aes(x = pt_rel, y = values, color = what, fill = what)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                     labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  scale_fill_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                    labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  labs(x = "pseudotime", y = "normalized (early - late)", color = NULL, fill = NULL) +
  theme_bw() + pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/mahalanobis-dist_mono.pdf", width = 4.5, height = 4, dpi = 400)
#cowplot::ggsave2(p, file = "../output/mahalanobis-dist_mono2.pdf", width = 2.7, height = 1.5, dpi = 400)

plot_df_kx <- compute_pseudotime_simliarities_SE2(so_cd34k)
plot_df_kx$pt_rel <- (plot_df_kx$mean_pt - rng[1]) / diff(rng)
ggplot(plot_df_kx, aes(x = pt_rel, y = values, color = what, fill = what)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                     labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  scale_fill_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                    labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  labs(x = "pseudotime", y = "normalized (early - late)", color = NULL, fill = NULL) +
  theme_bw() + pretty_plot() + L_border()
#cowplot::ggsave2(p, file = "../output/mahalanobis-dist_mega.pdf", width = 4.5, height = 4, dpi = 400)
#cowplot::ggsave2(p, file = "../output/mahalanobis-dist_mega2.pdf", width = 2.7, height = 1.5, dpi = 400)


# sanity check: permute = TRUE
plot_df_exp <- compute_pseudotime_simliarities_SE2(so_cd34e, PERMUTE = TRUE)
p <- ggplot(plot_df_exp, aes(x = idx, y = values, color = what, fill = what)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                     labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  scale_fill_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                    labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  labs(x = "pseudotime", y = "normalized (early - late)",
       color = NULL, fill = NULL) +
  theme_bw() + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/mahalanobis-dist_ery_permuted.pdf", width = 2.7, height = 1.5, dpi = 400)

plot_df_mxp <- compute_pseudotime_simliarities_SE2(so_cd34m, PERMUTE = TRUE)
p <- ggplot(plot_df_mxp, aes(x = idx, y = values, color = what, fill = what)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                     labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  scale_fill_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                    labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  labs(x = "pseudotime", y = "normalized (early - late)",
       color = NULL, fill = NULL) +
  theme_bw() + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/mahalanobis-dist_mono_permuted.pdf", width = 2.7, height = 1.5, dpi = 400)

plot_df_kxp <- compute_pseudotime_simliarities_SE2(so_cd34k, PERMUTE = TRUE)
p <- ggplot(plot_df_kxp, aes(x = idx, y = values, color = what, fill = what)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                     labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  scale_fill_manual(values = c("tRNA" = "#2f3a73", "mRNA" = "#b92f54"),
                    labels = c("tRNA" = "tRNAome", "mRNA" = "Transcriptome")) +
  labs(x = "pseudotime", y = "normalized (early - late)",
       color = NULL, fill = NULL) +
  theme_bw() + pretty_plot() + L_border()
cowplot::ggsave2(p, file = "../output/mahalanobis-dist_mega_permuted.pdf", width = 2.7, height = 1.5, dpi = 400)
