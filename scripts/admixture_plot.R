#!/usr/bin/env Rscript
# ADMIXTURE barplots for K = 4, 7, 8, 9, 14, 15, 19

# Individuals ordered and grouped by population (from llf_popmap.txt)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(viridis)
})


setwd("/Users/bethhasl/Desktop/new_figs")

# ---- Read popmap ----
popmap <- read.table("llf_popmap.txt", header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE,
                      col.names = c("Sample", "Population"))

# ---- Custom population order ----
custom_order <- c(
  "AravaipaEast", "AravaipaHorsecamp", "AravaipaWest", "Dudleyville",
  "LowerSanPedro", "7BRanch", "WhitlockWash", "EdgarCanyon", "BuehmanCanyon",
  "BullockCanyon", "YoutcyCanyon", "EspirituCanyon", "AlambreTank",
  "ChivaFalls", "JoaquinCanyon", "TanqueVerde", "3LinkFarm", "HookerHotSpring", "MintSpring",
  "BassSpring", "MadronaCanyon", "ChimeneaCanyon", "WolfHouse", "RinconCreek",
  "BoxCanyon", "CienegaDelLago", "DavidsonCanyon", "CienegaHorseshoeBend",
  "CienegaCreekPreserve", "Narrows", "PumpCanyon", "Coldsprings",
  "LowerEmpireGulch", "CatalinaRegionalPonds", "RomeroCanyon",
  "UniversityArizonaPond", "WildhorsePool", "PenaBlanca"
)

# Sanity check: make sure every population in popmap is present in custom_order
missing_pops <- setdiff(unique(popmap$Population), custom_order)
if (length(missing_pops) > 0) {
  stop(paste0("Populations in popmap not found in custom_order: ",
              paste(missing_pops, collapse = ", ")))
}


# ---- Function to build one ADMIXTURE plot for a given K ----
make_admixture_plot <- function(qfile, k, outfile) {

  q <- read.table(qfile, header = FALSE, stringsAsFactors = FALSE)
  colnames(q) <- paste0("Q", seq_len(ncol(q)))

  if (nrow(q) != nrow(popmap)) {
    stop(paste0("Row count mismatch between ", qfile,
                " (", nrow(q), ") and popmap (", nrow(popmap), ")"))
  }

  q$Sample <- popmap$Sample
  q$Population <- popmap$Population

  # Order individuals: grouped by population using custom_order,
  # keep original sample order within each population.
  # Reverse the population order so that, after coord_flip(),
  # AravaipaEast (first in custom_order) appears at the TOP of the plot.
  q$Population <- factor(q$Population, levels = rev(custom_order))
  q <- q %>% arrange(Population, Sample)

  q$Sample <- factor(q$Sample, levels = q$Sample)  # preserve plotting order

  # Assign a plotting position for each individual, adding an extra gap
  # between WildhorsePool and PenaBlanca blocks for visual separation.
  q <- q %>% mutate(row_order = row_number())
  gap_after_pops <- c("WildhorsePool", "PenaBlanca")  # gap inserted after these populations' blocks
  gap_size <- 2  # number of blank "slots" to insert as spacing per gap



  q <- q %>% mutate(x_pos = row_order)
  for (gp in gap_after_pops) {
    gap_start_row <- max(q$row_order[q$Population == gp])
    q <- q %>% mutate(x_pos = ifelse(row_order > gap_start_row,
                                      x_pos + gap_size,
                                      x_pos))
  }


  # Long format for stacked bar plot
  q_long <- q %>%
    pivot_longer(cols = starts_with("Q"), names_to = "Cluster", values_to = "Ancestry")

  # Compute population label positions (center of each population block)
  pop_positions <- q %>%
    group_by(Population) %>%
    summarise(mid = mean(range(x_pos)), xmin = min(x_pos) - 0.5, xmax = max(x_pos) + 0.5) %>%
    ungroup()


  n_clusters <- k
  cluster_colors <- viridis(n_clusters, option = "turbo")



  p <- ggplot(q_long, aes(x = x_pos, y = Ancestry, fill = Cluster, group = Sample)) +
    geom_col(width = 1, color = NA) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_vline(xintercept = pop_positions$xmax[-nrow(pop_positions)],
               color = "black", linewidth = 0.8) +
    labs(x = NULL, y = "Ancestry Proportion") +
    coord_flip(ylim = c(0, 1), clip = "off") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 16),
      panel.grid = element_blank(),

      plot.margin = margin(30, 100, 30, 240)

    ) +
    annotate("text",
             x = pop_positions$mid,
             y = -0.03,
             label = pop_positions$Population,
             angle = 0, hjust = 1, vjust = 0.5, size = 6)



  ggsave(outfile, p, width = 12, height = max(17, nrow(q) * 0.055),
         dpi = 300, limitsize = FALSE)




  message("Saved: ", outfile)
}

# ---- Generate plots for each K ----
make_admixture_plot("recode_reference_plink.4.Q",  4,  "ADMIXTURE_K4_byPopulation.png")
make_admixture_plot("recode_reference_plink.7.Q",  7,  "ADMIXTURE_K7_byPopulation.png")
make_admixture_plot("recode_reference_plink.8.Q",  8,  "ADMIXTURE_K8_byPopulation.png")
make_admixture_plot("recode_reference_plink.9.Q",  9,  "ADMIXTURE_K9_byPopulation.png")
make_admixture_plot("recode_reference_plink.14.Q", 14, "ADMIXTURE_K14_byPopulation.png")
make_admixture_plot("recode_reference_plink.15.Q", 15, "ADMIXTURE_K15_byPopulation.png")
make_admixture_plot("recode_reference_plink.19.Q", 19, "ADMIXTURE_K19_byPopulation.png")


