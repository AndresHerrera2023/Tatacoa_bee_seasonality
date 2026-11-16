# ============================================================
# 02 - Rarefaction and extrapolation (iNEXT)
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Compare species richness and diversity (Hill numbers q = 0, 1, 2) between Forest and Scrubland using incidence data.
# Input:   data/database.csv
# Output:  figures/inext/S2_Fig.tiff, figures/inext/Fig2.tiff, tables/inext/*.csv
# ============================================================

# Packages
library(iNEXT)
library(ggplot2)

# 1. Paths and data
# Project root: the folder that contains data/, figures/ and tables/.
# Run this script from the project root (or change this path).
base_dir <- "."

fig_dir <- file.path(base_dir, "figures", "inext")
tab_dir <- file.path(base_dir, "tables", "inext")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

database <- read.csv(
  file.path(base_dir, "data", "database.csv"),
  stringsAsFactors = FALSE
)

# 2. The 12 sampling months
# (only used to count the 12 sampling units, so the order does not matter here)
months <- c(
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

# 3. Function to prepare incidence-frequency data
# Output: number of sampling units (12 months) followed by the number of months
# in which each species was recorded
make_incidence_freq <- function(data) {
  
  data$Month <- factor(data$Month, levels = months)
  
  incidence_table <- table(
    data$Genus_species,
    data$Month
  )
  
  incidence_matrix <- incidence_table > 0
  
  species_frequency <- rowSums(incidence_matrix)
  
  c(ncol(incidence_matrix), species_frequency)
}

# 4. Separate the two vegetation types
forest <- subset(database, Area == "Forest")
scrubland <- subset(database, Area == "Scrubs")

# 5. Incidence-frequency vectors
freq_all <- make_incidence_freq(database)
freq_forest <- make_incidence_freq(forest)
freq_scrubland <- make_incidence_freq(scrubland)

freq_all[1]
freq_forest[1]
freq_scrubland[1]

# 6. OVERALL ASSEMBLAGE (S2 Fig)
# Rarefaction and extrapolation of Hill numbers q = 0, 1, 2
all_data <- list(
  "Overall" = freq_all
)

set.seed(123)

inext_all <- iNEXT(
  all_data,
  q = c(0, 1, 2),
  datatype = "incidence_freq",
  nboot = 200
)

inext_all

# Save asymptotic estimates
write.csv(
  inext_all$AsyEst,
  file.path(tab_dir, "all_inext.csv"),
  row.names = FALSE
)

# Plot
gg_all <- ggiNEXT(
  inext_all,
  type = 1,
  facet.var = "Order.q",
  color.var = "Assemblage"
) +
  theme_bw() +
  facet_wrap(
    ~ factor(
      Order.q,
      levels = c(0, 1, 2),
      labels = c("q0", "q1", "q2")
    )
  ) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.spacing = unit(0.5, "cm"),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = "darkred") +
  scale_colour_manual(values = "darkred")

gg_all

# Save figure (S2 Fig)
ggsave(
  filename = "S2_Fig.tiff",
  plot = gg_all,
  path = fig_dir,
  width = 15,
  height = 10,
  units = "cm",
  compression = "lzw",
  dpi = 300
)

# 7. FOREST vs SCRUBLAND (Fig 2)
forest_scrubland_data <- list(
  "Forest" = freq_forest,
  "Scrubland" = freq_scrubland
)

set.seed(123)

inext_fst_scb <- iNEXT(
  forest_scrubland_data,
  q = c(0, 1, 2),
  datatype = "incidence_freq",
  nboot = 200
)

inext_fst_scb

# Save asymptotic estimates
write.csv(
  inext_fst_scb$AsyEst,
  file.path(tab_dir, "forest_scrubland_inext.csv"),
  row.names = FALSE
)

# Plot
gg_fst_scb <- ggiNEXT(
  inext_fst_scb,
  type = 1,
  facet.var = "Order.q",
  color.var = "Assemblage"
) +
  theme_bw() +
  facet_wrap(
    ~ factor(
      Order.q,
      levels = c(0, 1, 2),
      labels = c("q0", "q1", "q2")
    )
  ) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.spacing = unit(0.5, "cm"),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(
    values = c(
      "Forest" = "green4",
      "Scrubland" = "orange"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Forest" = "green4",
      "Scrubland" = "orange"
    )
  )

gg_fst_scb

# Save figure (Fig 2)
ggsave(
  filename = "Fig2.tiff",
  plot = gg_fst_scb,
  path = fig_dir,
  width = 15,
  height = 10,
  units = "cm",
  compression = "lzw",
  dpi = 300
)
