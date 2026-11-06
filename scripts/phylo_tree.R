library(vcfR)
library(ape)
library(adegenet)
library(StAMPP)
library(ggplot2)
library(ggtree)

# Regions
regions <- c("HookerHotSpring" = "Galiuro Mtn", "MintSpring" = "Galiuro Mtn", "BassSpring" = "Galiuro Mtn",
             "BuehmanCanyon" = "Reddington Pass", "BullockCanyon" = "Reddington Pass", "EdgarCanyon" = "Reddington Pass",
             "YoutcyCanyon" = "Reddington Pass", "EspirituCanyon" = "Reddington Pass", "CienegaCreekPreserve" = "Whetstone-Empire Complex",
             "DavidsonCanyon" = "Whetstone-Empire Complex", "Coldsprings" = "Whetstone-Empire Complex", "Narrows" = "Whetstone-Empire Complex",
             "LowerEmpireGulch" = "Whetstone-Empire Complex", "PumpCanyon" = "Whetstone-Empire Complex", "CienegaDelLago" = "Whetstone-Empire Complex",
             "CienegaHorseshoeBend" = "Whetstone-Empire Complex", "RinconCreek" = "Rincon Mtn", "WolfHouse" = "Rincon Mtn",
             "ChimeneaCanyon" = "Rincon Mtn", "BoxCanyon" = "Rincon Mtn", "MadronaCanyon" = "Rincon Mtn",
             "TanqueVerde" = "Reddington Pass", "3LinkFarm" = "San Pedro River Valley", "Dudleyville" = "San Pedro River Valley",
             "LowerSanPedro" = "San Pedro River Valley", "WhitlockWash" = "San Pedro River Valley", "AravaipaEast" = "Galiuro Mtn",
             "AravaipaHorsecamp" = "Galiuro Mtn", "AravaipaWest" = "Galiuro Mtn", "UniversityArizonaPond" = "Misc. Site",
             "RomeroCanyon" = "Catalina Mtn", "7BRanch" = "San Pedro River Valley", "PenaBlanca" = "Misc. Site",
             "CatalinaRegionalPonds" = "Catalina Mtn", "ChivaFalls" = "Reddington Pass", "JoaquinCanyon" = "Reddington Pass",
             "AlambreTank" = "Reddington Pass", "WildhorsePool" = "Misc. Site")

# Colors
colors <- c("Galiuro Mtn" = "#00BA38", "Catalina Mtn" = "#619CFF", "Rincon Mtn" = "#FF8C00",
            "Whetstone-Empire Complex" = "#C77CFF", "San Pedro River Valley" = "#FFD700", "Reddington Pass" = "#F8766D",
            "Misc. Site" = "#999999", "PenaBlanca" = "gray50")


vcf_file    <- "populations.snps.vcf"
popmap_file <- "llf_popmap.txt"

vcf <- read.vcfR(vcf_file, verbose = FALSE)
gl <- vcfR2genlight(vcf)

popmap <- read.table(popmap_file, header = FALSE, stringsAsFactors = FALSE)
colnames(popmap) <- c("Individual", "Population")

# Filter genlight object to only include individuals present in popmap
gl <- gl[indNames(gl) %in% popmap$Individual]

# Assign populations directly from the popmap (population names already given)
pop_map_vec <- setNames(popmap$Population, popmap$Individual)
pop(gl) <- pop_map_vec[indNames(gl)]
ploidy(gl) <- 2

cat("Calculating population-level Nei's genetic distance...\n")
pop_dist <- as.dist(stamppNeisD(gl, pop = TRUE))
tree <- nj(pop_dist)

tree$edge.length[tree$edge.length < 0] <- 0

if ("PenaBlanca" %in% tree$tip.label) {
  tree <- root(tree, outgroup = "PenaBlanca", resolve.root = TRUE)
}


# Prepare plot data
tree_data <- data.frame(
  label = tree$tip.label,
  Region = sapply(tree$tip.label, function(x) {
    r <- regions[x]
    if (is.na(r)) "Misc. Site" else r
  })
)

# ---- Create vertical tree plot ------------------------------------------------
out_pdf_vertical <- "population_tree_vertical.pdf"
pdf(out_pdf_vertical, width = 14, height = 14)
p_vertical <- ggtree(tree, layout = "rectangular", ladderize = TRUE)
print(
  p_vertical %<+% tree_data +
    geom_tiplab(aes(color = Region), size = 6, hjust = -0.05, fontface = "bold", show.legend = FALSE) +
    geom_point(aes(color = Region), size = 0, alpha = 0) +
    scale_color_manual(values = colors, na.value = "gray50", na.translate = FALSE,
                       breaks = c("Whetstone-Empire Complex", "Catalina Mtn", "Galiuro Mtn",
                                 "Rincon Mtn", "Reddington Pass", "San Pedro River Valley", "Misc. Site")) +
    theme_tree() +
    theme(legend.position = "right",
          legend.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 18),
          plot.margin = margin(1, 3, 1, 1, "cm")) +
    guides(color = guide_legend(title = "Region", override.aes = list(size = 6, alpha = 1, shape = 16, linetype = 0))) +
    xlim(NA, max(tree$edge.length) * 1.3)
)
dev.off()
cat("Saved", out_pdf_vertical, "\n")


# ---- Create circular tree plot -------------------------------------------------
out_pdf_circular <- "population_tree_circular.pdf"
pdf(out_pdf_circular, width = 14, height = 14)
p_circular <- ggtree(tree, layout = "inward_circular", ladderize = TRUE)
print(
  p_circular %<+% tree_data +
    geom_tiplab(aes(color = Region), size = 6, hjust = 0.5, show.legend = FALSE) +
    geom_point(aes(color = Region), size = 0, alpha = 0) +
    scale_color_manual(values = colors, na.value = "gray50", na.translate = FALSE,
                       breaks = c("Whetstone-Empire Complex", "Catalina Mtn", "Galiuro Mtn",
                                 "Rincon Mtn", "Reddington Pass", "San Pedro River Valley", "Misc. Site")) +
    theme_tree() +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    guides(color = guide_legend(title = "Region", override.aes = list(size = 6, alpha = 1)))
)
dev.off()
cat("Saved", out_pdf_circular, "\n")

# ---- Create radial/fan tree plot ------------------------------------------------
out_pdf_fan <- "population_tree_fan.pdf"
pdf(out_pdf_fan, width = 14, height = 14)
p_fan <- ggtree(tree, layout = "fan", ladderize = TRUE, open.angle = 10)
print(
  p_fan %<+% tree_data +
    geom_tiplab(aes(color = Region), size = 6, hjust = -0.05, show.legend = FALSE) +
    geom_point(aes(color = Region), size = 0, alpha = 0) +
    scale_color_manual(values = colors, na.value = "gray50", na.translate = FALSE,
                       breaks = c("Whetstone-Empire Complex", "Catalina Mtn", "Galiuro Mtn",
                                 "Rincon Mtn", "Reddington Pass", "San Pedro River Valley", "Misc. Site")) +
    theme_tree() +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    guides(color = guide_legend(title = "Region", override.aes = list(size = 6, alpha = 1)))
)
dev.off()
cat("Saved", out_pdf_fan, "\n")

# ---- Create vertical tree with colored edges by region -------------------------
out_pdf_colored_edges <- "population_tree_vertical_colored_edges.pdf"
pdf(out_pdf_colored_edges, width = 14, height = 14)
p_colored_edges <- ggtree(tree, layout = "rectangular", ladderize = TRUE, aes(color = Region))
print(
  p_colored_edges %<+% tree_data +
    geom_tiplab(aes(color = Region), size = 6, hjust = -0.05, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = colors, na.value = "gray50", na.translate = FALSE,
                       breaks = c("Whetstone-Empire Complex", "Catalina Mtn", "Galiuro Mtn",
                                 "Rincon Mtn", "Reddington Pass", "San Pedro River Valley", "Misc. Site")) +
    theme_tree() +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14),
          plot.margin = margin(1, 3, 1, 1, "cm")) +
    guides(color = guide_legend(title = "Region", override.aes = list(size = 4, alpha = 1))) +
    xlim(NA, max(tree$edge.length) * 1.3)
)
dev.off()
cat("Saved", out_pdf_colored_edges, "\n")

write.tree(tree, "population_tree.nwk")
cat("Done! Trees saved:\n")
cat("  - population_tree_vertical.pdf\n")
cat("  - population_tree_circular.pdf\n")
cat("  - population_tree_fan.pdf\n")
cat("  - population_tree_vertical_colored_edges.pdf\n")
cat("  - population_tree.nwk\n")
