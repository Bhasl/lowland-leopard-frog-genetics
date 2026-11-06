library(vcfR)
library(adegenet)
library(vegan)
library(ggplot2)
library(parallel)

vcf_file      <- "populations.snps.vcf"
popmap_file   <- "llf_popmap.txt"
out_pca_png   <- "PCA_plot_noRomero_noMisc.png"
out_pca_csv   <- "PCA_scores_noRomero_noMisc.csv"
out_permanova <- "PERMANOVA_results_noRomero_noMisc.txt"

exclude_pops <- c("RomeroCanyon", "UniversityArizonaPond", "WildhorsePool", "PenaBlanca")

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

keep_idx <- !(pop(gl) %in% exclude_pops)
cat("Removing", sum(!keep_idx), "individual(s) from population(s):",
    paste(exclude_pops, collapse = ", "), "\n")
gl <- gl[keep_idx, ]
popmap <- popmap[keep_idx, ]

# Drop loci that became fully missing after excluding these populations,
# then mean-impute any remaining scattered missing genotypes
# (glPca requires no NAs when computing means across the remaining individuals)
gl_mat <- as.matrix(gl)
fully_na_loci <- colSums(is.na(gl_mat)) == nrow(gl_mat)
cat("Dropping", sum(fully_na_loci), "of", ncol(gl_mat),
    "loci that are fully missing after exclusion\n")
gl_mat <- gl_mat[, !fully_na_loci]

for (j in seq_len(ncol(gl_mat))) {
  col <- gl_mat[, j]
  if (any(is.na(col))) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    gl_mat[, j] <- col
  }
}
gl_imputed <- new("genlight", gl_mat, ploidy = rep(2, nrow(gl_mat)))

indNames(gl_imputed) <- indNames(gl)
pop(gl_imputed) <- pop(gl)


cat("Running PCA...\n")
pca <- glPca(gl_imputed, nf = 10, parallel = TRUE, n.cores = n_cores)

var_explained <- 100 * pca$eig / sum(pca$eig)


pca_df <- data.frame(
  sample = indNames(gl),
  pop    = pop(gl),
  PC1    = pca$scores[, 1],
  PC2    = pca$scores[, 2]
)

write.csv(pca_df, out_pca_csv, row.names = FALSE)
cat("Saved PCA scores to", out_pca_csv, "\n")

write.csv(data.frame(PC = paste0("PC", seq_along(var_explained)),
                      pct_variance = var_explained),
          "PCA_variance_explained_noRomero_noMisc.csv", row.names = FALSE)
cat("Saved PCA variance explained to PCA_variance_explained_noRomero_noMisc.csv\n")

# Regional grouping for coloring
region_map <- list(
  "Galiuro Mtn" = c("AravaipaEast", "AravaipaHorsecamp", "AravaipaWest",
                     "HookerHotSpring", "MintSpring", "BassSpring"),
  "San Pedro River Valley" = c("Dudleyville", "LowerSanPedro", "7BRanch",
                                "WhitlockWash", "3LinkFarm"),
  "Reddington Pass" = c("EdgarCanyon", "BuehmanCanyon", "BullockCanyon",
                         "YoutcyCanyon", "EspirituCanyon", "AlambreTank",
                         "ChivaFalls", "JoaquinCanyon", "TanqueVerde"),
  "Rincon Mtn" = c("MadronaCanyon", "ChimeneaCanyon", "WolfHouse",
                    "RinconCreek", "BoxCanyon"),
  "Catalina Mtn" = c("CatalinaRegionalPonds"),  # RomeroCanyon removed
  "Whetstone-Empire Complex" = c("CienegaDelLago", "DavidsonCanyon",
                                   "CienegaHorseshoeBend", "CienegaCreekPreserve",
                                   "Narrows", "PumpCanyon", "Coldsprings",
                                   "LowerEmpireGulch")
  # "Other" region removed entirely (misc sites excluded)
)

region_colors <- c(
  "Galiuro Mtn" = "green",
  "San Pedro River Valley" = "yellow",
  "Reddington Pass" = "#D9534F",
  "Rincon Mtn" = "orange",
  "Catalina Mtn" = "blue",
  "Whetstone-Empire Complex" = "purple"
)

pop_to_region <- setNames(
  rep(names(region_map), lengths(region_map)),
  unlist(region_map)
)

pca_df$region <- pop_to_region[as.character(pca_df$pop)]

unmatched <- unique(pca_df$pop[is.na(pca_df$region)])
if (length(unmatched) > 0) {
  warning(paste("The following populations were not found in the region map:",
                paste(unmatched, collapse = ", ")))
}

pca_df$region <- factor(pca_df$region, levels = names(region_map))

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = region_colors, name = "Region") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )

ggsave(out_pca_png, plot = p, width = 10, height = 8, dpi = 300)
cat("Saved PCA plot to", out_pca_png, "\n")


# ---------------------------------------------------------------------------
# Same PCA plot, but with 95% confidence ellipses drawn around each region
# ---------------------------------------------------------------------------
out_pca_png_ellipse <- "PCA_plot_noRomero_noMisc_ellipses.png"

p_ellipse <- ggplot(pca_df, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 0.7, na.rm = TRUE) +
  scale_color_manual(values = region_colors, name = "Region") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )

ggsave(out_pca_png_ellipse, plot = p_ellipse, width = 10, height = 8, dpi = 300)

cat("Saved PCA plot with 95% CI ellipses to", out_pca_png_ellipse, "\n")

cat("Running PERMANOVA...\n")


geno_mat <- as.matrix(gl)
for (j in seq_len(ncol(geno_mat))) {
  col <- geno_mat[, j]
  if (any(is.na(col))) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    geno_mat[, j] <- col
  }
}

dist_mat <- dist(geno_mat, method = "euclidean")

permanova_result <- adonis2(dist_mat ~ pop, data = popmap, permutations = 999)

print(permanova_result)

sink(out_permanova)
cat("PERMANOVA (vegan::adonis2)\n")
cat("RomeroCanyon and misc/'Other' sites (UniversityArizonaPond, WildhorsePool, PenaBlanca) excluded\n")
cat("Distance: Euclidean (SNP allele dosage)\n")
cat("Model: distance ~ pop\n")
cat("Permutations: 999\n\n")
print(permanova_result)
sink()

cat("Saved PERMANOVA results to", out_permanova, "\n")
