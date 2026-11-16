# ============================================================
# 06 - Temporal beta diversity (betapart)
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Measure how bee composition changes among months in Forest and Scrubland and split it into turnover and nestedness (Sorensen, incidence) and balanced variation and abundance gradient (Bray-Curtis, abundance).
# Input:   data/database.csv
# Output:  figures/baselga/S4_Fig.tif, Fig6.tif, S5_Fig.tif, Fig7.tif
# ============================================================

# Packages
library(betapart)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Helper for figure labels: drop the leading zero (0.37 -> .37),
# and write 0.0 as 0 and 1.0 as 1, so 8-pt labels fit in cells and bars
strip0 <- function(x) {
  x <- sub("^0\\.0+$", "0", x)
  x <- sub("^1\\.0+$", "1", x)
  sub("^0\\.", ".", x)
}

# 1. Paths and data
# Project root: the folder that contains data/ and figures/.
# Run this script from the project root (or change this path).
base_dir <- "."

data <- read.csv(
  file.path(base_dir, "data", "database.csv")
)

out_dir <- file.path(base_dir, "figures", "baselga")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Sampling months in chronological order
all_months <- c(
  "Feb","Mar","Apr","May","Jun","Jul",
  "Aug","Sep","Oct","Nov","Dec","Jan"
)

# Forest September had no bees, so it is excluded from the Forest matrices
months_forest <- c(
  "Feb","Mar","Apr","May","Jun","Jul",
  "Aug","Oct","Nov","Dec","Jan"
)

months_scrub <- all_months
all_species <- unique(data$Genus_species)

# 2. Subsets: Forest and Scrubland (called "Scrubs" in the database),
# each as all bees, solitary (non-social) and social (eusocial) bees
forest_all <- data[data$Area == "Forest", ]
forest_non <- data[data$Area == "Forest" &
                     data$Sociality == "Non-social", ]
forest_eus <- data[data$Area == "Forest" &
                     data$Sociality == "Eusocial", ]

scrub_all <- data[data$Area == "Scrubs", ]
scrub_non <- data[data$Area == "Scrubs" &
                    data$Sociality == "Non-social", ]
scrub_eus <- data[data$Area == "Scrubs" &
                    data$Sociality == "Eusocial", ]

# 3. Month x species abundance matrices
forest_all_t <- table(
  factor(forest_all$Month, levels = all_months),
  factor(forest_all$Genus_species, levels = all_species)
)

forest_non_t <- table(
  factor(forest_non$Month, levels = all_months),
  factor(forest_non$Genus_species, levels = all_species)
)

forest_eus_t <- table(
  factor(forest_eus$Month, levels = all_months),
  factor(forest_eus$Genus_species, levels = all_species)
)

scrub_all_t <- table(
  factor(scrub_all$Month, levels = all_months),
  factor(scrub_all$Genus_species, levels = all_species)
)

scrub_non_t <- table(
  factor(scrub_non$Month, levels = all_months),
  factor(scrub_non$Genus_species, levels = all_species)
)

scrub_eus_t <- table(
  factor(scrub_eus$Month, levels = all_months),
  factor(scrub_eus$Genus_species, levels = all_species)
)

# Remove September from Forest
forest_all_t <- forest_all_t[rownames(forest_all_t) != "Sep", ]
forest_non_t <- forest_non_t[rownames(forest_non_t) != "Sep", ]
forest_eus_t <- forest_eus_t[rownames(forest_eus_t) != "Sep", ]

# 4. Presence/absence matrices
forest_all_pa <- (forest_all_t > 0) * 1
forest_non_pa <- (forest_non_t > 0) * 1
forest_eus_pa <- (forest_eus_t > 0) * 1

scrub_all_pa <- (scrub_all_t > 0) * 1
scrub_non_pa <- (scrub_non_t > 0) * 1
scrub_eus_pa <- (scrub_eus_t > 0) * 1

# Optional check: richness per month
rowSums(forest_all_pa)
rowSums(forest_non_pa)
rowSums(forest_eus_pa)

rowSums(scrub_all_pa)
rowSums(scrub_non_pa)
rowSums(scrub_eus_pa)

# PART I. PRESENCE/ABSENCE
# Sorensen dissimilarity = turnover + nestedness
#
# 5. Sorensen beta diversity among months
beta_forest_all <- beta.pair(forest_all_pa, index.family = "sorensen")
beta_forest_non <- beta.pair(forest_non_pa, index.family = "sorensen")
beta_forest_eus <- beta.pair(forest_eus_pa, index.family = "sorensen")

beta_scrub_all <- beta.pair(scrub_all_pa, index.family = "sorensen")
beta_scrub_non <- beta.pair(scrub_non_pa, index.family = "sorensen")
beta_scrub_eus <- beta.pair(scrub_eus_pa, index.family = "sorensen")

# 6. Extract the Sorensen matrices (total, turnover = sim, nestedness = sne)
forest_all_sor <- as.matrix(beta_forest_all$beta.sor)
forest_all_sim <- as.matrix(beta_forest_all$beta.sim)
forest_all_sne <- as.matrix(beta_forest_all$beta.sne)

forest_non_sor <- as.matrix(beta_forest_non$beta.sor)
forest_non_sim <- as.matrix(beta_forest_non$beta.sim)
forest_non_sne <- as.matrix(beta_forest_non$beta.sne)

forest_eus_sor <- as.matrix(beta_forest_eus$beta.sor)
forest_eus_sim <- as.matrix(beta_forest_eus$beta.sim)
forest_eus_sne <- as.matrix(beta_forest_eus$beta.sne)

scrub_all_sor <- as.matrix(beta_scrub_all$beta.sor)
scrub_all_sim <- as.matrix(beta_scrub_all$beta.sim)
scrub_all_sne <- as.matrix(beta_scrub_all$beta.sne)

scrub_non_sor <- as.matrix(beta_scrub_non$beta.sor)
scrub_non_sim <- as.matrix(beta_scrub_non$beta.sim)
scrub_non_sne <- as.matrix(beta_scrub_non$beta.sne)

scrub_eus_sor <- as.matrix(beta_scrub_eus$beta.sor)
scrub_eus_sim <- as.matrix(beta_scrub_eus$beta.sim)
scrub_eus_sne <- as.matrix(beta_scrub_eus$beta.sne)

# 7. Function that converts a matrix into a table for heatmaps
mat_to_df <- function(mat, months){
  df <- as.data.frame(as.table(mat))
  names(df) <- c("Month1","Month2","Beta")
  df$Month1 <- factor(df$Month1, levels = months)
  df$Month2 <- factor(df$Month2, levels = rev(months))
  df
}

# 8. Sorensen heatmaps (S4 Fig)
forest_all_sor_df <- mat_to_df(forest_all_sor, months_forest)
forest_non_sor_df <- mat_to_df(forest_non_sor, months_forest)
forest_eus_sor_df <- mat_to_df(forest_eus_sor, months_forest)

scrub_all_sor_df <- mat_to_df(scrub_all_sor, months_scrub)
scrub_non_sor_df <- mat_to_df(scrub_non_sor, months_scrub)
scrub_eus_sor_df <- mat_to_df(scrub_eus_sor, months_scrub)

# Function for one heatmap
plot_heat_sor <- function(df, title){
  
  ggplot(df, aes(Month1, Month2, fill = Beta)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = strip0(sprintf("%.1f", Beta))),
      size = 8 / .pt
    ) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      limits = c(0,1)
    ) +
    coord_equal() +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = expression(beta[SOR])
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0
      )
    )
}

# Panels a-c: Forest (all bees, solitary, social); d-f: Scrubland
h_sor_a <- plot_heat_sor(forest_all_sor_df, "a")#"a  Forest - Overall"
h_sor_b <- plot_heat_sor(forest_non_sor_df, "b")#"b  Forest - Non-social"
h_sor_c <- plot_heat_sor(forest_eus_sor_df, "c")#"c  Forest - Eusocial"
h_sor_d <- plot_heat_sor(scrub_all_sor_df, "d")#"d  Scrubs - Overall"
h_sor_e <- plot_heat_sor(scrub_non_sor_df, "e")#"e  Scrubs - Non-social"
h_sor_f <- plot_heat_sor(scrub_eus_sor_df, "f")#"f  Scrubs - Eusocial"

# Combine the six panels (3 columns x 2 rows) with one legend
heatmaps_sor_3x2 <- ggarrange(
  h_sor_a, h_sor_b, h_sor_c,
  h_sor_d, h_sor_e, h_sor_f,
  ncol = 3,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)


# Figures
heatmaps_sor_3x2

# Save S4 Fig
ggsave(
  file.path(out_dir, "S4_Fig.tif"),
  heatmaps_sor_3x2,
  width = 19,
  height = 14,
  dpi = 300,
  compression = "lzw",
  units = "cm",
  bg = "white"
  
)

# 9. Dissimilarity between consecutive months, partitioned into
# turnover and nestedness (proportion of the total Sorensen value)
get_consecutive_baselga <- function(sim_mat, sne_mat, sor_mat, months){
  
  from <- months[-length(months)]
  to <- months[-1]
  
  out <- do.call(rbind, lapply(seq_along(from), function(i){
    
    bsim <- sim_mat[from[i], to[i]]
    bsne <- sne_mat[from[i], to[i]]
    bsor <- sor_mat[from[i], to[i]]
    
    data.frame(
      Transition = paste(from[i], to[i], sep = "-"),
      BetaSOR = bsor,
      Component = c("Turnover","Nestedness"),
      Value = c(bsim, bsne),
      Proportion = if(bsor == 0) c(0,0)
      else c(bsim/bsor, bsne/bsor)
    )
  }))
  
  out$Transition <- factor(
    out$Transition,
    levels = paste(from, to, sep = "-")
  )
  
  out
}

# Apply to each group
forest_all_consec <- get_consecutive_baselga(
  forest_all_sim, forest_all_sne, forest_all_sor, months_forest
)

forest_non_consec <- get_consecutive_baselga(
  forest_non_sim, forest_non_sne, forest_non_sor, months_forest
)

forest_eus_consec <- get_consecutive_baselga(
  forest_eus_sim, forest_eus_sne, forest_eus_sor, months_forest
)

scrub_all_consec <- get_consecutive_baselga(
  scrub_all_sim, scrub_all_sne, scrub_all_sor, months_scrub
)

scrub_non_consec <- get_consecutive_baselga(
  scrub_non_sim, scrub_non_sne, scrub_non_sor, months_scrub
)

scrub_eus_consec <- get_consecutive_baselga(
  scrub_eus_sim, scrub_eus_sne, scrub_eus_sor, months_scrub
)

# 10. Sorensen bar plots (Fig 6)
# Number above each bar = total Sorensen dissimilarity of that transition
plot_sor_bar <- function(df, title){
  
  labels <- df %>% distinct(Transition, BetaSOR)
  
  ggplot(df, aes(Transition, Proportion, fill = Component)) +
    geom_col(width = 0.8) +
    geom_text(
      data = labels,
      aes(Transition, 1.05, label = strip0(as.character(round(BetaSOR, 2)))),
      inherit.aes = FALSE,
      size = 8 / .pt
    ) +
    scale_fill_manual(
      values = c(
        "Turnover" = "#3B7191",
        "Nestedness" = "#BDBDBD"
      )
    ) +
    scale_y_continuous(
      limits = c(0,1.12),
      breaks = seq(0,1,0.25),
      expand = c(0,0)
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold")
    )
}

# Panels a-f
p_sor_a <- plot_sor_bar(forest_all_consec, "a")#"a  Forest - Overall"
p_sor_b <- plot_sor_bar(forest_non_consec, "b")#"b  Forest - Non-social"
p_sor_c <- plot_sor_bar(forest_eus_consec, "c")#"c  Forest - Eusocial"
p_sor_d <- plot_sor_bar(scrub_all_consec, "d")#"d  Scrubs - Overall"
p_sor_e <- plot_sor_bar(scrub_non_consec, "e")#"e  Scrubs - Non-social"
p_sor_f <- plot_sor_bar(scrub_eus_consec, "f")#"f  Scrubs - Eusocial"

# Combine the six panels with one legend and one common y-axis title
bars_sor_3x2 <- ggarrange(
  p_sor_a, p_sor_b, p_sor_c,
  p_sor_d, p_sor_e, p_sor_f,
  ncol = 3,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

bars_sor_3x2 <- annotate_figure(
  bars_sor_3x2,
  left = text_grob(
    "Proportion of beta diversity",
    rot = 90,
    face = "bold",
    size = 12
  )
)

bars_sor_3x2

# Save Fig 6
ggsave(
  file.path(out_dir, "Fig6.tif"),
  bars_sor_3x2,
  width = 19,
  height = 14,
  units = "cm",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

# PART II. ABUNDANCE
# Bray-Curtis dissimilarity = balanced variation + abundance gradient
#
# 11. Bray-Curtis beta diversity among months
beta_ab_forest_all <- beta.pair.abund(
  forest_all_t,
  index.family = "bray"
)

beta_ab_forest_non <- beta.pair.abund(
  forest_non_t,
  index.family = "bray"
)

beta_ab_forest_eus <- beta.pair.abund(
  forest_eus_t,
  index.family = "bray"
)

beta_ab_scrub_all <- beta.pair.abund(
  scrub_all_t,
  index.family = "bray"
)

beta_ab_scrub_non <- beta.pair.abund(
  scrub_non_t,
  index.family = "bray"
)

beta_ab_scrub_eus <- beta.pair.abund(
  scrub_eus_t,
  index.family = "bray"
)

# 12. Extract the Bray-Curtis matrices (total, balanced variation = bal,
# abundance gradient = gra)
forest_all_bray <- as.matrix(beta_ab_forest_all$beta.bray)
forest_all_bal <- as.matrix(beta_ab_forest_all$beta.bray.bal)
forest_all_gra <- as.matrix(beta_ab_forest_all$beta.bray.gra)

forest_non_bray <- as.matrix(beta_ab_forest_non$beta.bray)
forest_non_bal <- as.matrix(beta_ab_forest_non$beta.bray.bal)
forest_non_gra <- as.matrix(beta_ab_forest_non$beta.bray.gra)

forest_eus_bray <- as.matrix(beta_ab_forest_eus$beta.bray)
forest_eus_bal <- as.matrix(beta_ab_forest_eus$beta.bray.bal)
forest_eus_gra <- as.matrix(beta_ab_forest_eus$beta.bray.gra)

scrub_all_bray <- as.matrix(beta_ab_scrub_all$beta.bray)
scrub_all_bal <- as.matrix(beta_ab_scrub_all$beta.bray.bal)
scrub_all_gra <- as.matrix(beta_ab_scrub_all$beta.bray.gra)

scrub_non_bray <- as.matrix(beta_ab_scrub_non$beta.bray)
scrub_non_bal <- as.matrix(beta_ab_scrub_non$beta.bray.bal)
scrub_non_gra <- as.matrix(beta_ab_scrub_non$beta.bray.gra)

scrub_eus_bray <- as.matrix(beta_ab_scrub_eus$beta.bray)
scrub_eus_bal <- as.matrix(beta_ab_scrub_eus$beta.bray.bal)
scrub_eus_gra <- as.matrix(beta_ab_scrub_eus$beta.bray.gra)

# 13. Bray-Curtis heatmaps (S5 Fig)
forest_all_bray_df <- mat_to_df(forest_all_bray, months_forest)
forest_non_bray_df <- mat_to_df(forest_non_bray, months_forest)
forest_eus_bray_df <- mat_to_df(forest_eus_bray, months_forest)

scrub_all_bray_df <- mat_to_df(scrub_all_bray, months_scrub)
scrub_non_bray_df <- mat_to_df(scrub_non_bray, months_scrub)
scrub_eus_bray_df <- mat_to_df(scrub_eus_bray, months_scrub)

# Function for one heatmap
plot_heat_bray <- function(df, title){
  
  ggplot(df, aes(Month1, Month2, fill = Beta)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = strip0(sprintf("%.1f", Beta))),
      size = 8 / .pt
    ) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      limits = c(0,1)
    ) +
    coord_equal() +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = expression(beta[Bray])
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0
      )
    )
}

# Panels a-c: Forest (all bees, solitary, social); d-f: Scrubland
h_bray_a <- plot_heat_bray(forest_all_bray_df, "a") #"a  Forest - Overall"
h_bray_b <- plot_heat_bray(forest_non_bray_df, "b")#"b  Forest - Non-social"
h_bray_c <- plot_heat_bray(forest_eus_bray_df, "c")#"c  Forest - Eusocial"
h_bray_d <- plot_heat_bray(scrub_all_bray_df, "d")#"d  Scrubs - Overall"
h_bray_e <- plot_heat_bray(scrub_non_bray_df, "e")#"e  Scrubs - Non-social"
h_bray_f <- plot_heat_bray(scrub_eus_bray_df, "f")#"f  Scrubs - Eusocial"

# Combine the six panels with one legend
heatmaps_bray_3x2 <- ggarrange(
  h_bray_a, h_bray_b, h_bray_c,
  h_bray_d, h_bray_e, h_bray_f,
  ncol = 3,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

heatmaps_bray_3x2

# Save S5 Fig
ggsave(
  filename = file.path(out_dir, "S5_Fig.tif"),
  plot = heatmaps_bray_3x2,
  width = 19,
  height = 14,
  units = "cm",
  dpi = 300,
  compression = "lzw",
  bg = "white",
  device = "tiff"
)

# 14. Dissimilarity between consecutive months, partitioned into
# balanced variation and abundance gradient
get_consecutive_bray <- function(bal_mat, gra_mat, bray_mat, months){
  
  from <- months[-length(months)]
  to <- months[-1]
  
  out <- do.call(rbind, lapply(seq_along(from), function(i){
    
    bbal <- bal_mat[from[i], to[i]]
    bgra <- gra_mat[from[i], to[i]]
    bbray <- bray_mat[from[i], to[i]]
    
    data.frame(
      Transition = paste(from[i], to[i], sep = "-"),
      BetaBray = bbray,
      Component = c(
        "Balanced variation",
        "Abundance gradient"
      ),
      Value = c(bbal, bgra),
      Proportion = if(bbray == 0) c(0,0)
      else c(bbal/bbray, bgra/bbray)
    )
  }))
  
  out$Transition <- factor(
    out$Transition,
    levels = paste(from, to, sep = "-")
  )
  
  out
}

# Apply to each group
forest_all_bray_consec <- get_consecutive_bray(
  forest_all_bal, forest_all_gra,
  forest_all_bray, months_forest
)

forest_non_bray_consec <- get_consecutive_bray(
  forest_non_bal, forest_non_gra,
  forest_non_bray, months_forest
)

forest_eus_bray_consec <- get_consecutive_bray(
  forest_eus_bal, forest_eus_gra,
  forest_eus_bray, months_forest
)

scrub_all_bray_consec <- get_consecutive_bray(
  scrub_all_bal, scrub_all_gra,
  scrub_all_bray, months_scrub
)

scrub_non_bray_consec <- get_consecutive_bray(
  scrub_non_bal, scrub_non_gra,
  scrub_non_bray, months_scrub
)

scrub_eus_bray_consec <- get_consecutive_bray(
  scrub_eus_bal, scrub_eus_gra,
  scrub_eus_bray, months_scrub
)

# 15. Bray-Curtis bar plots (Fig 7)
# Number above each bar = total Bray-Curtis dissimilarity of that transition
plot_bray_bar <- function(df, title){
  
  labels <- df %>% distinct(Transition, BetaBray)
  
  ggplot(df, aes(Transition, Proportion, fill = Component)) +
    geom_col(width = 0.8) +
    geom_text(
      data = labels,
      aes(Transition, 1.05, label = strip0(as.character(round(BetaBray, 2)))),
      inherit.aes = FALSE,
      size = 8 / .pt
    ) +
    scale_fill_manual(
      values = c(
        "Balanced variation" = "#756BB1",
        "Abundance gradient" = "#BCBDDC"
      )
    ) +
    scale_y_continuous(
      limits = c(0,1.12),
      breaks = seq(0,1,0.25),
      expand = c(0,0)
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold")
    )
}

# Panels a-f
p_bray_a <- plot_bray_bar(forest_all_bray_consec,"a")#"a  Forest - Overall"
p_bray_b <- plot_bray_bar(forest_non_bray_consec,"b")#"b  Forest - Non-social"
p_bray_c <- plot_bray_bar(forest_eus_bray_consec,"c")#"c  Forest - Eusocial"
p_bray_d <- plot_bray_bar(scrub_all_bray_consec,"d")#"d  Scrubs - Overall"
p_bray_e <- plot_bray_bar(scrub_non_bray_consec,"e")#"e  Scrubs - Non-social"
p_bray_f <- plot_bray_bar(scrub_eus_bray_consec,"f")#"f  Scrubs - Eusocial"

# Combine the six panels with one legend and one common y-axis title
bars_bray_3x2 <- ggarrange(
  p_bray_a, p_bray_b, p_bray_c,
  p_bray_d, p_bray_e, p_bray_f,
  ncol = 3,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

bars_bray_3x2 <- annotate_figure(
  bars_bray_3x2,
  left = text_grob(
    "Proportion of beta diversity",
    rot = 90,
    face = "bold",
    size = 12
  )
)

bars_bray_3x2

# Save Fig 7
ggsave(
  file.path(out_dir, "Fig7.tif"),
  bars_bray_3x2,
  width = 19,
  height = 14,
  units = "cm",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)


# 16. View results
# Sorensen: consecutive months
forest_all_consec
forest_non_consec
forest_eus_consec

scrub_all_consec
scrub_non_consec
scrub_eus_consec


# Bray-Curtis: consecutive months
forest_all_bray_consec
forest_non_bray_consec
forest_eus_bray_consec

scrub_all_bray_consec
scrub_non_bray_consec
scrub_eus_bray_consec

heatmaps_sor_3x2
bars_sor_3x2
heatmaps_bray_3x2
bars_bray_3x2
