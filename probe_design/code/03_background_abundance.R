library(stringdist)
library(reshape2)
library(dplyr)
library(data.table)
library(readxl)
library(BuenColors)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load reference data
# -----------------------------------------------------------------------------

full_reference <- fread("../output/trna_reference.csv")
colnames(full_reference) <- c("V1", "name_long", "name_short", "aa", "ac", "RHS", "LHS", "rhs_probe_g3", "lhs_probe_g3", "gapfill_g3")
g3_cols <- c("V1", "name_long", "name_short", "aa", "ac",
             "RHS", "LHS", "rhs_probe_g3", "lhs_probe_g3", "gapfill_g3")
simplified_g3 <- full_reference[, ..g3_cols]

# -----------------------------------------------------------------------------
# 2. Build exon-only rows for intron-containing tRNAs
# -----------------------------------------------------------------------------

reverse_complement <- function(seq) {
  complement_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  paste(rev(complement_map[unlist(strsplit(seq, ""))]), collapse = "")
}

g3_exon_ref <- read_excel("../output/tRNA_probes_intron.xlsx") %>%
  dplyr::rename(name_short = names_short) %>%
  dplyr::mutate(
    LHS_exon     = mapply(function(lhs, intron) gsub(intron, "", lhs), LHS, intron),
    name_long    = paste0(name_long, " exon"),
    name_short   = paste0(name_short, "-exon"),
    lhs_probe_g3 = substr(LHS_exon, 1, 25),
    gapfill_g3   = sapply(codon, reverse_complement)
  )

g3_exon_ref_reduced <- g3_exon_ref %>%
  dplyr::select(V1 = X, name_long, name_short, aa = AA, ac = codon, RHS,
                LHS = LHS_exon, rhs_probe_g3 = RHS_25mer, lhs_probe_g3, gapfill_g3)

simplified_g3_exon <- rbind(simplified_g3, g3_exon_ref_reduced)

full_reference_exon <- rbind(
  full_reference %>% dplyr::select(name_short, rhs_probe_g3, lhs_probe_g3),
  g3_exon_ref   %>% dplyr::select(name_short, rhs_probe_g3 = RHS_25mer, lhs_probe_g3)
)

# -----------------------------------------------------------------------------
# 3. Pairwise edit distances (lower triangle, long format)
# -----------------------------------------------------------------------------

pairwise_distances_long <- function(probes, names) {
  dist_mat <- data.matrix(stringdistmatrix(probes, method = "lv"))
  rownames(dist_mat) <- names
  colnames(dist_mat) <- names
  reshape2::melt(dist_mat) %>%
    dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
    dplyr::filter(Var1 >= Var2)
}

melted_lhs_df <- pairwise_distances_long(simplified_g3_exon$lhs_probe_g3,
                                         simplified_g3_exon$name_short)
melted_rhs_df <- pairwise_distances_long(simplified_g3_exon$rhs_probe_g3,
                                         simplified_g3_exon$name_short)

# -----------------------------------------------------------------------------
# 4. Greedy random clustering: assign each tRNA a representative probe
#
#    Among tRNAs still unassigned, find all near-duplicate pairs (LV dist < 2),
#    pick one at random as the representative, and assign its probe to every
#    tRNA linked to it. set.seed() fires before every draw, so LHS and RHS
#    passes are fully RNG-isolated regardless of call order.
# -----------------------------------------------------------------------------

design_probes <- function(items, melted_dist, seed, verbose = FALSE) {
  design <- rep("missing", length(items))
  names(design) <- items

  while (any(design == "missing")) {
    if (verbose) print(sum(design == "missing"))

    need <- names(design)[design == "missing"]

    candidate_pairs <- melted_dist %>%
      dplyr::filter(value < 2, Var1 %in% need, Var2 %in% need)

    set.seed(seed)
    chosen <- sample(unique(c(candidate_pairs$Var1, candidate_pairs$Var2)), 1)

    bag <- candidate_pairs %>%
      dplyr::filter(Var1 %in% chosen | Var2 %in% chosen) %>%
      { unique(c(.$Var1, .$Var2)) }

    design[names(design) %in% bag] <- chosen
  }

  unname(design)
}

# -----------------------------------------------------------------------------
# 5. Per-seed run: assign LHS and RHS together, then count unique pairs
#
#    A "probe pair" is the (designed_lhs, designed_rhs) combination for a
#    tRNA. n_pairs counts distinct combinations -- two tRNAs share a pair only
#    if they were assigned the same representative on *both* sides.
# -----------------------------------------------------------------------------

run_seed <- function(seed, items, melted_lhs, melted_rhs, verbose = FALSE) {
  message("Running with seed = ", seed)
  lhs <- design_probes(items, melted_lhs, seed, verbose = verbose)
  rhs <- design_probes(items, melted_rhs, seed, verbose = verbose)
  list(
    lhs     = lhs,
    rhs     = rhs,
    n_lhs   = length(unique(lhs)),
    n_rhs   = length(unique(rhs)),
    n_pairs = nrow(unique(data.frame(lhs = lhs, rhs = rhs)))
  )
}

# -----------------------------------------------------------------------------
# 6. Run across seeds 1-1000
# -----------------------------------------------------------------------------

seeds <- 1:1000
items <- simplified_g3_exon$name_short

results_list <- lapply(seeds, run_seed,
                       items      = items,
                       melted_lhs = melted_lhs_df,
                       melted_rhs = melted_rhs_df)
names(results_list) <- as.character(seeds)

results_summary <- data.frame(
  seed    = seeds,
  n_lhs   = sapply(results_list, `[[`, "n_lhs"),
  n_rhs   = sapply(results_list, `[[`, "n_rhs"),
  n_pairs = sapply(results_list, `[[`, "n_pairs")
)

print(results_summary)

# -----------------------------------------------------------------------------
# 7. Z-scores / p-values
#
#    Mean and SD are computed from the random-seed distribution only
#    The rationally designed counts are compared against that null --
#    not included in the mean/SD calculation.
# -----------------------------------------------------------------------------

# Rationally designed probe counts -- update n_pairs when known
designed_lhs   <- 145
designed_rhs   <- 89
designed_pairs <- 181 # verify these values

stats <- results_summary %>%
  dplyr::summarize(
    mean_lhs   = mean(n_lhs,   na.rm = TRUE),
    sd_lhs     = sd(n_lhs,     na.rm = TRUE),
    mean_rhs   = mean(n_rhs,   na.rm = TRUE),
    sd_rhs     = sd(n_rhs,     na.rm = TRUE),
    mean_pairs = mean(n_pairs, na.rm = TRUE),
    sd_pairs   = sd(n_pairs,   na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    z_score_lhs   = (designed_lhs   - mean_lhs)   / sd_lhs,
    p_value_lhs   = 2 * pnorm(-abs(z_score_lhs)),
    z_score_rhs   = (designed_rhs   - mean_rhs)   / sd_rhs,
    p_value_rhs   = 2 * pnorm(-abs(z_score_rhs)),
    z_score_pairs = (designed_pairs - mean_pairs) / sd_pairs,
    p_value_pairs = 2 * pnorm(-abs(z_score_pairs))
  )

cat("\n=== Stats (z-scores relative to random-seed distribution) ===\n")
print(stats[, c("mean_lhs", "sd_lhs", "z_score_lhs", "p_value_lhs",
                 "mean_rhs", "sd_rhs", "z_score_rhs", "p_value_rhs",
                 "mean_pairs", "sd_pairs", "z_score_pairs", "p_value_pairs")])

# -----------------------------------------------------------------------------
# 8. Histograms: random-seed distribution with median and designed-set lines
# -----------------------------------------------------------------------------

make_hist <- function(counts, designed_val, x_label) {
  freq_df <- as.data.frame(table(counts))
  freq_df$counts <- as.numeric(as.character(freq_df$counts))
  colnames(freq_df) <- c("count", "n_freq")

  med_val  <- median(counts, na.rm = TRUE)
  bar_min  <- min(counts, na.rm = TRUE)
  bar_max  <- max(counts, na.rm = TRUE)
  bar_span <- bar_max - bar_min
  print(med_val)
  # Extend x range to show designed value; add 5% padding on each side
  x_lo <- min(bar_min, designed_val, na.rm = TRUE) 
  x_hi <- max(c(max(bar_max, designed_val, na.rm = TRUE) , designed_val))

  # Place labels to left of their line if at the far right, else to the right
  med_hjust      <- if (med_val > (x_lo + x_hi) / 2) 1.15 else -0.15
  designed_hjust <- if (!is.na(designed_val) && designed_val > med_val) 1.15 else -0.15

  p <- ggplot(freq_df, aes(x = count, y = n_freq)) +
    geom_bar(stat = "identity", fill = "lightgrey", color = "black", width = 0.5) +
    geom_vline(xintercept = med_val,
               color = "grey30", linetype = "dashed", linewidth = 0.5) +
    scale_x_continuous(limits = c(x_lo, x_hi)) +
    labs(x = x_label, y = "Frequency") +
    pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0))

  if (!is.na(designed_val)) {
    p <- p +
      geom_vline(xintercept = designed_val,
                 color = "darkblue", linetype = "dashed", linewidth = 0.5)
  }

  p
}

p_lhs   <- make_hist(results_summary$n_lhs,   designed_lhs,   "Number of LHS probes")
p_rhs   <- make_hist(results_summary$n_rhs,   designed_rhs,   "Number of RHS probes")
p_pairs <- make_hist(results_summary$n_pairs, designed_pairs, "Number of probe pairs")

panel <- p_lhs + p_rhs + p_pairs +
  plot_layout(ncol = 3) 

cowplot::ggsave2("../plots/probe_design_histograms.pdf", panel, width = 5.5, height = 1.6)
