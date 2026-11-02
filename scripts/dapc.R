#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# DAPC (Discriminant Analysis of Principal Components)
#
# Reads a VCF + popmap (Stacks format: SampleID<TAB>PopulationName) and runs
# the same analysis pipeline twice:
#   1) On all individuals.
#   2) On a subset with RomeroCanyon, PenaBlanca, UniversityArizonaPond, and
#      WildhorsePool individuals excluded.
#
# For each run, the pipeline:
#   1) Infers genetic clusters de novo using adegenet::find.clusters()
#      (K-means on PCA-transformed genotypes, K chosen automatically by
#      lowest BIC across K = 1..max_k).
#   2) Runs DAPC on the inferred clusters (adegenet::dapc) and plots
#      individuals on the first two discriminant functions, colored by
#      inferred cluster, plus a labeled version and a text summary.
#
# Requirements: vcfR, adegenet, ggplot2, dplyr, ggrepel, RColorBrewer, scales
# (This is a refactor of dapc_plot.R: the two near-identical analysis blocks
# there have been consolidated into run_dapc_pipeline() below, called once
# per dataset.)
# ---------------------------------------------------------------------------

suppressMessages({
  library(vcfR)
  library(adegenet)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

set.seed(42)

# ---- paths -------------------------------------------------------------
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))))

vcf_file    <- "populations.snps.vcf"
popmap_file <- "llf_popmap.txt"
max_k       <- 20

# ---- load data -----------------------------------------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)

cat("Reading popmap...\n")
popmap <- read.table(popmap_file, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("sample", "pop"))

gl <- vcfR2genlight(vcf)

sample_order <- indNames(gl)
popmap <- popmap[match(sample_order, popmap$sample), ]

if (any(is.na(popmap$pop))) {
  missing <- sample_order[is.na(popmap$pop)]
  stop(paste("The following samples in the VCF are missing from the popmap:",
             paste(missing, collapse = ", ")))
}

pop(gl) <- popmap$pop
ploidy(gl) <- 2

# ---------------------------------------------------------------------------
# run_dapc_pipeline(): shared find.clusters + DAPC + plotting + summary logic
#
# Arguments:
#   gl                    genlight object (already subset, if applicable)
#   popmap                data.frame with columns sample/pop, matching gl
#   suffix                string appended to output filenames (e.g. "" or
#                         "_noROM_PB_UA_WHP")
#   exclude_pops          character vector of population names to drop from
#                         gl/popmap before clustering (NULL = use as-is)
#   max_k                 maximum K to test in find.clusters()
#   cluster_name_map      optional named character vector, e.g.
#                         c("1"="GSR","2"="ROM",...), used to relabel
#                         numeric cluster IDs with meaningful names. Only
#                         applied when the observed K exactly matches
#                         length(cluster_name_map); otherwise raw numeric
#                         cluster labels are kept (avoids partial relabeling
#                         if K differs from what was expected).
#   cluster_colors_fixed  optional named color vector matching the (possibly
#                         relabeled) cluster levels. Used only when its
#                         length matches the observed K; otherwise colors
#                         are chosen automatically (RColorBrewer -> scales
#                         -> rainbow fallback).
#   summary_header        header line written at the top of the summary .txt
#
# Returns (invisibly): list(grp, dapc_res, dapc_df, best_k)
# ---------------------------------------------------------------------------
run_dapc_pipeline <- function(gl, popmap, suffix = "",
                               exclude_pops = NULL,
                               max_k = 20,
                               cluster_name_map = NULL,
                               cluster_colors_fixed = NULL,
                               summary_header = "DAPC Summary") {

  out_dapc_png         <- paste0("DAPC_plot", suffix, ".png")
  out_dapc_labeled_png <- paste0("DAPC_plot_labeled", suffix, ".png")
  out_bic_png          <- paste0("DAPC_BIC", suffix, ".png")
  out_dapc_csv         <- paste0("DAPC_scores", suffix, ".csv")
  out_summary          <- paste0("DAPC_summary", suffix, ".txt")

  # ---- Optionally exclude populations -------------------------------------
  if (!is.null(exclude_pops)) {
    keep_idx <- !(pop(gl) %in% exclude_pops)
    cat("Removing", sum(!keep_idx), "individual(s) from population(s):",
        paste(exclude_pops, collapse = ", "), "\n")
    gl <- gl[keep_idx, ]
    popmap <- popmap[keep_idx, ]

    # Drop loci that are entirely missing (all NA) in this subset of
    # individuals, since glPca cannot compute a column mean for such loci.
    geno_check <- as.matrix(gl)
    loci_all_na <- colSums(!is.na(geno_check)) == 0
    if (any(loci_all_na)) {
      cat("Dropping", sum(loci_all_na), "loci with no data in this subset...\n")
      gl <- gl[, !loci_all_na]
    }
    rm(geno_check)
  }

  # ---- Identify clusters de novo (find.clusters) --------------------------
  cat("Running find.clusters() to identify genetic clusters de novo...\n")
  cat("Testing K = 1 to", max_k, ", selecting K with lowest BIC (fully automated)...\n")

  n_pca_axes <- min(100, nInd(gl) - 1)

  grp <- find.clusters(
    gl,
    n.pca       = n_pca_axes,
    max.n.clust = max_k,
    choose.n.clust = FALSE,   # fully automated, no user prompt
    criterion   = "diffNgroup", # pick K at the elbow / diminishing returns in BIC
    stat        = "BIC"
  )

  best_k <- length(levels(grp$grp))
  cat("Selected K =", best_k, "clusters (lowest BIC)\n")

  # Save BIC curve plot
  bic_df <- data.frame(K = seq_along(grp$Kstat), BIC = grp$Kstat)
  p_bic <- ggplot(bic_df, aes(x = K, y = BIC)) +
    geom_line(linewidth = 0.8, color = "steelblue") +
    geom_point(size = 2, color = "steelblue") +
    geom_vline(xintercept = best_k, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = bic_df$K) +
    labs(x = "Number of clusters (K)", y = "BIC",
         title = sprintf("BIC vs. K (selected K = %d)", best_k)) +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"))

  ggsave(out_bic_png, plot = p_bic, width = 8, height = 6, dpi = 300)
  cat("Saved BIC plot to", out_bic_png, "\n")

  # ---- Run DAPC on inferred clusters ---------------------------------------
  cat("Running DAPC on inferred clusters...\n")

  # Choose number of PCs to retain for DAPC (avoid overfitting): a
  # conservative rule of thumb of min(n.pca axes used above, nInd/3)
  n_pcs_dapc <- min(n_pca_axes, floor(nInd(gl) / 3))
  n_da <- min(best_k - 1, 10)  # number of discriminant functions to retain

  dapc_res <- dapc(gl, grp$grp, n.pca = n_pcs_dapc, n.da = n_da)

  # Percent variance explained by each discriminant function
  var_explained <- 100 * dapc_res$eig / sum(dapc_res$eig)

  dapc_df <- data.frame(
    sample  = indNames(gl),
    pop     = pop(gl),
    cluster = factor(grp$grp),
    LD1     = dapc_res$ind.coord[, 1]
  )
  dapc_df$LD2 <- if (n_da >= 2) dapc_res$ind.coord[, 2] else 0

  # Optionally relabel numeric cluster IDs with meaningful region names.
  # Only applied when the observed K exactly matches length(cluster_name_map)
  # to avoid partially relabeling clusters if K differs from what's expected.
  if (!is.null(cluster_name_map) && best_k == length(cluster_name_map)) {
    dapc_df$cluster <- factor(
      cluster_name_map[as.character(dapc_df$cluster)],
      levels = cluster_name_map
    )
  }

  write.csv(dapc_df, out_dapc_csv, row.names = FALSE)
  cat("Saved DAPC scores to", out_dapc_csv, "\n")

  # ---- Choose cluster colors -----------------------------------------------
  n_clusters <- nlevels(dapc_df$cluster)

  if (!is.null(cluster_colors_fixed) && length(cluster_colors_fixed) == n_clusters) {
    cluster_colors <- cluster_colors_fixed
  } else if (n_clusters <= 8) {
    safe_cols <- tryCatch(RColorBrewer::brewer.pal(max(3, n_clusters), "Set1"),
                           error = function(e) NULL)
    cluster_colors <- if (!is.null(safe_cols)) safe_cols[seq_len(n_clusters)] else scales::hue_pal()(n_clusters)
  } else {
    cluster_colors <- grDevices::rainbow(n_clusters)
  }

  # ---- Plot DAPC (LD1 vs LD2, colored by inferred cluster) -----------------
  if (n_da >= 2) {
    p_dapc <- ggplot(dapc_df, aes(x = LD1, y = LD2, color = cluster)) +
      geom_point(size = 2.2, alpha = 0.85) +
      stat_ellipse(level = 0.95, linewidth = 0.6, na.rm = TRUE) +
      scale_color_manual(values = cluster_colors, name = "Inferred\nCluster") +
      labs(
        x = sprintf("DAPC1 (%.1f%%)", var_explained[1]),
        y = sprintf("DAPC2 (%.1f%%)", var_explained[2])
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)
      )
  } else {
    # Only one discriminant function available (K = 2): plot density instead
    p_dapc <- ggplot(dapc_df, aes(x = LD1, fill = cluster)) +
      geom_density(alpha = 0.6) +
      scale_fill_manual(values = cluster_colors, name = "Inferred\nCluster") +
      labs(
        x = sprintf("DAPC1 (%.1f%%)", var_explained[1]),
        y = "Density"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)
      )
  }

  ggsave(out_dapc_png, plot = p_dapc, width = 10, height = 8, dpi = 300)
  cat("Saved DAPC plot to", out_dapc_png, "\n")

  # ---- Labeled version: same DAPC plot with sample ID labels on each point --
  if (n_da >= 2) {
    p_dapc_labeled <- ggplot(dapc_df, aes(x = LD1, y = LD2, color = cluster)) +
      geom_point(size = 2.2, alpha = 0.85) +
      stat_ellipse(level = 0.95, linewidth = 0.6, na.rm = TRUE) +
      ggrepel::geom_text_repel(aes(label = sample), size = 2.5,
                                max.overlaps = Inf, show.legend = FALSE,
                                segment.size = 0.2, seed = 42) +
      scale_color_manual(values = cluster_colors, name = "Inferred\nCluster") +
      labs(
        x = sprintf("DAPC1 (%.1f%%)", var_explained[1]),
        y = sprintf("DAPC2 (%.1f%%)", var_explained[2])
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)
      )
  } else {
    p_dapc_labeled <- ggplot(dapc_df, aes(x = LD1, y = 0, color = cluster)) +
      geom_point(size = 2.2, alpha = 0.85) +
      ggrepel::geom_text_repel(aes(label = sample), size = 2.5,
                                max.overlaps = Inf, show.legend = FALSE,
                                segment.size = 0.2, seed = 42) +
      scale_color_manual(values = cluster_colors, name = "Inferred\nCluster") +
      labs(
        x = sprintf("DAPC1 (%.1f%%)", var_explained[1]),
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  ggsave(out_dapc_labeled_png, plot = p_dapc_labeled, width = 14, height = 11, dpi = 300)
  cat("Saved labeled DAPC plot to", out_dapc_labeled_png, "\n")

  # ---- Cross-tabulation: inferred cluster vs. known population -------------
  cross_tab <- table(Population = dapc_df$pop, Cluster = dapc_df$cluster)

  sink(out_summary)
  cat(summary_header, "\n")
  cat(strrep("=", nchar(summary_header)), "\n\n")
  cat("Method: adegenet::find.clusters() (K-means on PCA of SNP genotypes)\n")
  cat("K search range: 1 to", max_k, "\n")
  cat("Selection criterion: diffNgroup (elbow / diminishing returns in BIC)\n")
  cat("Selected K:", best_k, "\n")
  cat("PCA axes retained for find.clusters():", n_pca_axes, "\n")
  cat("PCA axes retained for DAPC:", n_pcs_dapc, "\n")
  cat("Discriminant functions retained:", n_da, "\n\n")
  cat("Percent variance explained by discriminant functions:\n")
  print(round(var_explained, 2))
  cat("\nCross-tabulation of known population vs. inferred cluster:\n\n")
  print(cross_tab)
  sink()

  cat("Saved DAPC summary to", out_summary, "\n")

  invisible(list(grp = grp, dapc_res = dapc_res, dapc_df = dapc_df, best_k = best_k))
}

# ---------------------------------------------------------------------------
# Run 1: all individuals
#
# Cluster 1 -> GSR (Galiuro Mtn / San Pedro River Valley / Reddington Pass)
# Cluster 2 -> ROM (RomeroCanyon)
# Cluster 3 -> WSC (Whetstone-Empire Complex)
# Cluster 4 -> RM  (Rincon Mtn)
# (relabeling is only applied if the observed K is exactly 4; see
# run_dapc_pipeline()'s cluster_name_map argument for details)
# ---------------------------------------------------------------------------
cat("==== Running DAPC on all individuals ====\n")

cluster_name_map_all <- c("1" = "GSR", "2" = "ROM", "3" = "WSC", "4" = "RM")
cluster_colors_all   <- c(GSR = "#E41A1C", ROM = "#377EB8",
                           WSC = "#984EA3", RM = "#FF7F00")

run_dapc_pipeline(
  gl, popmap,
  suffix               = "",
  max_k                = max_k,
  cluster_name_map     = cluster_name_map_all,
  cluster_colors_fixed = cluster_colors_all,
  summary_header       = "DAPC Summary"
)

# ---------------------------------------------------------------------------
# Run 2: RomeroCanyon, PenaBlanca, UniversityArizonaPond, and WildhorsePool
# individuals excluded
# ---------------------------------------------------------------------------
cat("\n\n==== Re-running DAPC with RomeroCanyon, PenaBlanca,",
    "UniversityArizonaPond, and WildhorsePool excluded ====\n")

exclude_pops <- c("RomeroCanyon", "PenaBlanca", "UniversityArizonaPond", "WildhorsePool")

run_dapc_pipeline(
  gl, popmap,
  suffix         = "_noROM_PB_UA_WHP",
  exclude_pops   = exclude_pops,
  max_k          = max_k,
  summary_header = "DAPC Summary (RomeroCanyon, PenaBlanca, UniversityArizonaPond, WildhorsePool excluded)"
)

cat("\nDone.\n")
