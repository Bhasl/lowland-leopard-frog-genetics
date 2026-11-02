#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Pairwise Fst heatmap (AlambreTank excluded - only 1 individual)
# Reads a VCF + popmap (Stacks format: SampleID<TAB>PopulationName)
# Computes pairwise Weir & Cockerham Fst (via StAMPP) with bootstrapped
# p-values and plots a heatmap with significance markers.
# ---------------------------------------------------------------------------

suppressMessages({
  library(vcfR)
  library(adegenet)
  library(StAMPP)
  library(reshape2)
  library(ggplot2)
})

# ---- paths -----------------------------------------------------------------
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))))

vcf_file    <- "populations.snps.vcf"
popmap_file <- "llf_popmap.txt"
out_csv     <- "pairwise_fst_matrix_noAlambre.csv"
out_pval_csv <- "pairwise_fst_pvalues_noAlambre.csv"
out_png     <- "pairwise_fst_heatmap_noAlambre.png"

# Population to exclude
exclude_pops <- c("AlambreTank")

# Number of bootstrap replicates for significance testing
n_bootstraps <- 1000

# ---- load data ---------------------------------------------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)

cat("Reading popmap...\n")
popmap <- read.table(popmap_file, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("sample", "pop"))

# Convert VCF -> genlight
gl <- vcfR2genlight(vcf)

# Match sample order between genlight object and popmap
sample_order <- indNames(gl)
popmap <- popmap[match(sample_order, popmap$sample), ]

if (any(is.na(popmap$pop))) {
  missing <- sample_order[is.na(popmap$pop)]
  stop(paste("The following samples in the VCF are missing from the popmap:",
             paste(missing, collapse = ", ")))
}

pop(gl) <- popmap$pop
ploidy(gl) <- 2

# ---- remove excluded population(s) -----------------------------------------
keep_idx <- !(pop(gl) %in% exclude_pops)
cat("Removing", sum(!keep_idx), "individual(s) from population(s):",
    paste(exclude_pops, collapse = ", "), "\n")
gl <- gl[keep_idx, ]

# ---- compute pairwise Fst (Weir & Cockerham 1984) with bootstrapped p-values ----
cat("Calculating pairwise Fst with", n_bootstraps, "bootstraps (this may take a while)...\n")
fst_result <- stamppFst(gl, nboots = n_bootstraps, percent = 95, nclusters = 1)

fst_mat  <- fst_result$Fsts
pval_mat <- fst_result$Pvalues

# stamppFst returns a lower-triangular matrix; mirror it to full symmetric matrix
fst_full <- fst_mat
fst_full[upper.tri(fst_full)] <- t(fst_mat)[upper.tri(fst_mat)]
diag(fst_full) <- 0

pval_full <- pval_mat
pval_full[upper.tri(pval_full)] <- t(pval_mat)[upper.tri(pval_mat)]
diag(pval_full) <- NA

write.csv(fst_full, out_csv)
write.csv(pval_full, out_pval_csv)
cat("Saved pairwise Fst matrix to", out_csv, "\n")
cat("Saved pairwise Fst p-values to", out_pval_csv, "\n")

# ---- plot heatmap (upper triangle only) ------------------------------------
# Custom population ordering (AlambreTank removed)
custom_order <- c(
  "AravaipaEast", "AravaipaHorsecamp", "AravaipaWest", "Dudleyville",
  "LowerSanPedro", "7BRanch", "WhitlockWash", "EdgarCanyon", "BuehmanCanyon",
  "BullockCanyon", "YoutcyCanyon", "EspirituCanyon",
  "ChivaFalls", "JoaquinCanyon", "TanqueVerde", "3LinkFarm", "HookerHotSpring", "MintSpring",
  "BassSpring", "MadronaCanyon", "ChimeneaCanyon", "WolfHouse", "RinconCreek",
  "BoxCanyon", "CienegaDelLago", "DavidsonCanyon", "CienegaHorseshoeBend",
  "CienegaCreekPreserve", "Narrows", "PumpCanyon", "Coldsprings",
  "LowerEmpireGulch", "CatalinaRegionalPonds", "RomeroCanyon",
  "UniversityArizonaPond", "WildhorsePool", "PenaBlanca"
)

available_pops <- rownames(fst_full)
missing_from_custom <- setdiff(available_pops, custom_order)
missing_from_data   <- setdiff(custom_order, available_pops)

if (length(missing_from_data) > 0) {
  warning(paste("The following populations from the custom order were not found in the data and will be skipped:",
                paste(missing_from_data, collapse = ", ")))
}
if (length(missing_from_custom) > 0) {
  warning(paste("The following populations in the data were not in the custom order list and will be appended at the end:",
                paste(missing_from_custom, collapse = ", ")))
}

pop_order <- c(custom_order[custom_order %in% available_pops], missing_from_custom)

fst_full  <- fst_full[pop_order, pop_order]
pval_full <- pval_full[pop_order, pop_order]

# Keep only the upper triangle (including diagonal); mask lower triangle as NA
# This puts AravaipaEast (first in pop_order) on the far LEFT of the x-axis
# and at the TOP of the y-axis, with the legend able to sit in the empty
# lower-left region of the plot.
fst_upper <- fst_full
fst_upper[lower.tri(fst_upper)] <- NA

pval_upper <- pval_full
pval_upper[lower.tri(pval_upper)] <- NA

fst_df <- melt(fst_upper, na.rm = FALSE)
colnames(fst_df) <- c("Pop1", "Pop2", "Fst")

pval_df <- melt(pval_upper, na.rm = FALSE)
colnames(pval_df) <- c("Pop1", "Pop2", "Pvalue")

fst_df$Pvalue <- pval_df$Pvalue

# Significance stars based on bootstrapped p-values
fst_df$sig <- ""
fst_df$sig[!is.na(fst_df$Pvalue) & fst_df$Pvalue < 0.05]  <- "*"
fst_df$sig[!is.na(fst_df$Pvalue) & fst_df$Pvalue < 0.01]  <- "**"
fst_df$sig[!is.na(fst_df$Pvalue) & fst_df$Pvalue < 0.001] <- "***"
# Diagonal (Fst = 0, comparing pop to itself) should have no stars
fst_df$sig[fst_df$Pop1 == fst_df$Pop2] <- ""

# Preserve population order (Pop1 = rows top->bottom, Pop2 = columns left->right)
fst_df$Pop1 <- factor(fst_df$Pop1, levels = pop_order)
fst_df$Pop2 <- factor(fst_df$Pop2, levels = pop_order)


n_pops <- length(pop_order)
plot_size <- max(8, n_pops * 0.35)

p <- ggplot(fst_df, aes(x = Pop2, y = Pop1, fill = Fst)) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_text(aes(label = sig), color = "black", size = 3, na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = mean(fst_df$Fst, na.rm = TRUE),
                        na.value = NA, name = "Fst") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.position = c(0.15, 0.85),
    legend.justification = c(0, 1),

    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.2, "cm")
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed()


ggsave(out_png, plot = p, width = plot_size, height = plot_size, dpi = 300, limitsize = FALSE)
cat("Saved heatmap to", out_png, "\n")
cat("Significance legend: * p<0.05, ** p<0.01, *** p<0.001 (bootstrapped, n =", n_bootstraps, ")\n")
