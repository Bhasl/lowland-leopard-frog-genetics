library(vcfR)
library(adegenet)
library(vegan)
library(ggplot2)
library(parallel)

vcf_file      <- "populations.snps.vcf"
popmap_file   <- "llf_popmap.txt"
out_pca_png   <- "PCA_plot.png"
out_pca_csv   <- "PCA_scores.csv"
out_permanova <- "PERMANOVA_results.txt"

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

# ---- PCA --------------------------------------------------------------------
cat("Running PCA...\n")
# glPca handles missing data internally (uses mean allele frequency imputation)
pca <- glPca(gl, nf = 10, parallel = TRUE, n.cores = n_cores)


# Percent variance explained
var_explained <- 100 * pca$eig / sum(pca$eig)

pca_df <- data.frame(
  sample = indNames(gl),
  pop    = pop(gl),
  PC1    = pca$scores[, 1],
  PC2    = pca$scores[, 2]
)

write.csv(pca_df, out_pca_csv, row.names = FALSE)
cat("Saved PCA scores to", out_pca_csv, "\n")

# Save variance explained so it can be reused for re-plotting without rerunning PCA
write.csv(data.frame(PC = paste0("PC", seq_along(var_explained)),
                      pct_variance = var_explained),
          "PCA_variance_explained.csv", row.names = FALSE)
cat("Saved PCA variance explained to PCA_variance_explained.csv\n")


# ---- Regional grouping for coloring ----------------------------------------
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
  "Catalina Mtn" = c("CatalinaRegionalPonds", "RomeroCanyon"),
  "Whetstone-Empire Complex" = c("CienegaDelLago", "DavidsonCanyon",
                                   "CienegaHorseshoeBend", "CienegaCreekPreserve",
                                   "Narrows", "PumpCanyon", "Coldsprings",
                                   "LowerEmpireGulch"),
  "Misc. Site" = c("UniversityArizonaPond", "WildhorsePool", "PenaBlanca")
)

region_colors <- c(
  "Galiuro Mtn" = "green",
  "San Pedro River Valley" = "yellow",
  "Reddington Pass" = "#D9534F",

  "Rincon Mtn" = "orange",
  "Catalina Mtn" = "blue",
  "Whetstone-Empire Complex" = "purple",
  "Misc. Site" = "grey"
)


# Build a lookup table: population -> region
pop_to_region <- setNames(
  rep(names(region_map), lengths(region_map)),
  unlist(region_map)
)

pca_df$region <- pop_to_region[as.character(pca_df$pop)]

unmatched <- unique(pca_df$pop[is.na(pca_df$region)])
if (length(unmatched) > 0) {
  warning(paste("The following populations were not found in the region map and will be grouped as 'Other':",
                paste(unmatched, collapse = ", ")))
  pca_df$region[is.na(pca_df$region)] <- "Other"
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

# ---- PERMANOVA ----------------------------------------------------------------
cat("Running PERMANOVA (this may take a while for many SNPs/samples)...\n")

# Build a genotype matrix (individuals x loci), impute missing with mean allele freq
geno_mat <- as.matrix(gl)
# Impute missing values (NA) with column (locus) mean
for (j in seq_len(ncol(geno_mat))) {
  col <- geno_mat[, j]
  if (any(is.na(col))) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    geno_mat[, j] <- col
  }
}

# Euclidean distance among individuals based on SNP genotypes
dist_mat <- dist(geno_mat, method = "euclidean")

permanova_result <- adonis2(dist_mat ~ pop, data = popmap, permutations = 999)

print(permanova_result)

sink(out_permanova)
cat("PERMANOVA (vegan::adonis2)\n")
cat("Distance: Euclidean (SNP allele dosage)\n")
cat("Model: distance ~ pop\n")
cat("Permutations: 999\n\n")
print(permanova_result)
sink()

cat("Saved PERMANOVA results to", out_permanova, "\n")

# ---------------------------------------------------------------------------
# ---- Repeat PCA + PERMANOVA with RomeroCanyon excluded ---------------------
# ---------------------------------------------------------------------------
cat("\n\n==== Re-running PCA + PERMANOVA with RomeroCanyon excluded ====\n")

exclude_pops_2 <- c("RomeroCanyon")

out_pca_png_2   <- "PCA_plot_noRomero.png"
out_pca_csv_2   <- "PCA_scores_noRomero.csv"
out_permanova_2 <- "PERMANOVA_results_noRomero.txt"

keep_idx_2 <- !(pop(gl) %in% exclude_pops_2)
cat("Removing", sum(!keep_idx_2), "individual(s) from population(s):",
    paste(exclude_pops_2, collapse = ", "), "\n")

gl2 <- gl[keep_idx_2, ]
popmap2 <- popmap[keep_idx_2, ]

cat("Running PCA (RomeroCanyon excluded)...\n")
pca2 <- glPca(gl2, nf = 10, parallel = TRUE, n.cores = n_cores)

var_explained_2 <- 100 * pca2$eig / sum(pca2$eig)

pca_df2 <- data.frame(
  sample = indNames(gl2),
  pop    = pop(gl2),
  PC1    = pca2$scores[, 1],
  PC2    = pca2$scores[, 2]
)

write.csv(pca_df2, out_pca_csv_2, row.names = FALSE)
cat("Saved PCA scores to", out_pca_csv_2, "\n")

write.csv(data.frame(PC = paste0("PC", seq_along(var_explained_2)),
                      pct_variance = var_explained_2),
          "PCA_variance_explained_noRomero.csv", row.names = FALSE)
cat("Saved PCA variance explained to PCA_variance_explained_noRomero.csv\n")

# Regional grouping for coloring (RomeroCanyon removed from Catalina Mtn)
region_map_2 <- region_map
region_map_2[["Catalina Mtn"]] <- setdiff(region_map_2[["Catalina Mtn"]], "RomeroCanyon")

pop_to_region_2 <- setNames(
  rep(names(region_map_2), lengths(region_map_2)),
  unlist(region_map_2)
)

pca_df2$region <- pop_to_region_2[as.character(pca_df2$pop)]

unmatched_2 <- unique(pca_df2$pop[is.na(pca_df2$region)])
if (length(unmatched_2) > 0) {
  warning(paste("The following populations were not found in the region map and will be grouped as 'Other':",
                paste(unmatched_2, collapse = ", ")))
  pca_df2$region[is.na(pca_df2$region)] <- "Other"
}

pca_df2$region <- factor(pca_df2$region, levels = names(region_map_2))

p2 <- ggplot(pca_df2, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = region_colors, name = "Region") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained_2[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained_2[2])
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )

ggsave(out_pca_png_2, plot = p2, width = 10, height = 8, dpi = 300)

cat("Saved PCA plot to", out_pca_png_2, "\n")

cat("Running PERMANOVA (RomeroCanyon excluded)...\n")

geno_mat2 <- as.matrix(gl2)
for (j in seq_len(ncol(geno_mat2))) {
  col <- geno_mat2[, j]
  if (any(is.na(col))) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    geno_mat2[, j] <- col
  }
}

dist_mat2 <- dist(geno_mat2, method = "euclidean")

permanova_result_2 <- adonis2(dist_mat2 ~ pop, data = popmap2, permutations = 999)

print(permanova_result_2)

sink(out_permanova_2)
cat("PERMANOVA (vegan::adonis2)\n")
cat("RomeroCanyon excluded\n")
cat("Distance: Euclidean (SNP allele dosage)\n")
cat("Model: distance ~ pop\n")
cat("Permutations: 999\n\n")
print(permanova_result_2)
sink()

cat("Saved PERMANOVA results to", out_permanova_2, "\n")

# ---------------------------------------------------------------------------
# ---- Repeat PCA + PERMANOVA with RomeroCanyon AND misc/"Other" sites excluded
# ---------------------------------------------------------------------------
cat("\n\n==== Re-running PCA + PERMANOVA with RomeroCanyon and misc sites excluded ====\n")

exclude_pops_3 <- c("RomeroCanyon", "UniversityArizonaPond", "WildhorsePool", "PenaBlanca")

out_pca_png_3   <- "PCA_plot_noRomero_noMisc.png"
out_pca_csv_3   <- "PCA_scores_noRomero_noMisc.csv"
out_permanova_3 <- "PERMANOVA_results_noRomero_noMisc.txt"

keep_idx_3 <- !(pop(gl) %in% exclude_pops_3)
cat("Removing", sum(!keep_idx_3), "individual(s) from population(s):",
    paste(exclude_pops_3, collapse = ", "), "\n")

gl3 <- gl[keep_idx_3, ]
popmap3 <- popmap[keep_idx_3, ]

cat("Running PCA (RomeroCanyon and misc sites excluded)...\n")
pca3 <- glPca(gl3, nf = 10, parallel = TRUE, n.cores = n_cores)

var_explained_3 <- 100 * pca3$eig / sum(pca3$eig)

pca_df3 <- data.frame(
  sample = indNames(gl3),
  pop    = pop(gl3),
  PC1    = pca3$scores[, 1],
  PC2    = pca3$scores[, 2]
)

write.csv(pca_df3, out_pca_csv_3, row.names = FALSE)
cat("Saved PCA scores to", out_pca_csv_3, "\n")

write.csv(data.frame(PC = paste0("PC", seq_along(var_explained_3)),
                      pct_variance = var_explained_3),
          "PCA_variance_explained_noRomero_noMisc.csv", row.names = FALSE)
cat("Saved PCA variance explained to PCA_variance_explained_noRomero_noMisc.csv\n")

# Regional grouping for coloring (RomeroCanyon removed from Catalina Mtn, "Other" region dropped entirely)
region_map_3 <- region_map
region_map_3[["Catalina Mtn"]] <- setdiff(region_map_3[["Catalina Mtn"]], "RomeroCanyon")
region_map_3[["Other"]] <- NULL

pop_to_region_3 <- setNames(
  rep(names(region_map_3), lengths(region_map_3)),
  unlist(region_map_3)
)

pca_df3$region <- pop_to_region_3[as.character(pca_df3$pop)]

unmatched_3 <- unique(pca_df3$pop[is.na(pca_df3$region)])
if (length(unmatched_3) > 0) {
  warning(paste("The following populations were not found in the region map:",
                paste(unmatched_3, collapse = ", ")))
}

pca_df3$region <- factor(pca_df3$region, levels = names(region_map_3))

p3 <- ggplot(pca_df3, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = region_colors, name = "Region") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained_3[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained_3[2])
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )

ggsave(out_pca_png_3, plot = p3, width = 10, height = 8, dpi = 300)

cat("Saved PCA plot to", out_pca_png_3, "\n")

cat("Running PERMANOVA (RomeroCanyon and misc sites excluded)...\n")

geno_mat3 <- as.matrix(gl3)
for (j in seq_len(ncol(geno_mat3))) {
  col <- geno_mat3[, j]
  if (any(is.na(col))) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    geno_mat3[, j] <- col
  }
}

dist_mat3 <- dist(geno_mat3, method = "euclidean")

permanova_result_3 <- adonis2(dist_mat3 ~ pop, data = popmap3, permutations = 999)

print(permanova_result_3)

sink(out_permanova_3)
cat("PERMANOVA (vegan::adonis2)\n")
cat("RomeroCanyon and misc/'Other' sites (UniversityArizonaPond, WildhorsePool, PenaBlanca) excluded\n")
cat("Distance: Euclidean (SNP allele dosage)\n")
cat("Model: distance ~ pop\n")
cat("Permutations: 999\n\n")
print(permanova_result_3)
sink()

cat("Saved PERMANOVA results to", out_permanova_3, "\n")


