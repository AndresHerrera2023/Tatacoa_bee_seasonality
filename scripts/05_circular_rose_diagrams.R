# ============================================================
# 05 - Circular statistics and rose diagrams
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Test whether bee richness and abundance are concentrated at particular times of the year (Rayleigh test) and draw rose diagrams.
# Input:   data/database.csv
# Output:  figures/circular/Fig4.tiff, Fig5.tiff, tables/phenology_richness_abundance.csv
# ============================================================

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(circular)
library(ggpubr)
library(grid)
# 1. Paths and data
# Project root: the folder that contains data/, figures/ and tables/.
# Run this script from the project root (or change this path).
base_dir <- "."

fig_dir <- file.path(base_dir, "figures", "circular")
tab_dir <- file.path(base_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(
  file.path(base_dir, "data", "database.csv"),
  stringsAsFactors = FALSE
)

# 2. Sampling-year order: February 2018 = first month, January 2019 = last month
month_levels <- c(
  "Feb", "Mar", "Apr", "May", "Jun", "Jul",
  "Aug", "Sep", "Oct", "Nov", "Dec", "Jan"
)

# Angles used for the STATISTICS: Feb = 0, Mar = 30, ..., Jan = 330 degrees
month_angles <- seq(
  0,
  330,
  by = 30
)

# Angles used ONLY for plotting: bars are shifted by 15 degrees so that
# February is not cut. This does not change the statistics.
plot_angles <- month_angles + 15

# 3. Groups: Forest and Scrubland, each as assemblage, solitary (non-social)
# and social (eusocial) bees
groups <- list(
  fst_assemblage = data %>%
    filter(
      Area == "Forest"
    ),
  fst_solitary = data %>%
    filter(
      Area == "Forest",
      Sociality == "Non-social"
    ),
  fst_social = data %>%
    filter(
      Area == "Forest",
      Sociality == "Eusocial"
    ),
  scb_assemblage = data %>%
    filter(
      Area == "Scrubs"
    ),
  scb_solitary = data %>%
    filter(
      Area == "Scrubs",
      Sociality == "Non-social"
    ),
  scb_social = data %>%
    filter(
      Area == "Scrubs",
      Sociality == "Eusocial"
    )
)

# 4. Monthly richness and abundance of each group
monthly_summary <- function(df, group_name) {

  # Count species and individuals per month
  temp <- df %>%
    group_by(Month) %>%
    summarise(
      Richness = n_distinct(
        Genus_species,
        na.rm = TRUE
      ),
      Abundance = n(),
      .groups = "drop"
    )

  # Complete 12-month calendar (months without bees get zero)
  result <- data.frame(
    Month = month_levels
  ) %>%
    left_join(
      temp,
      by = "Month"
    ) %>%
    mutate(
      Richness = replace_na(
        Richness,
        0
      ),
      Abundance = replace_na(
        Abundance,
        0
      ),
      Month = factor(
        Month,
        levels = month_levels,
        ordered = TRUE
      ),
      Group = group_name
    ) %>%
    arrange(Month)
  return(result)
}

# Monthly table for all groups
monthly_data <- bind_rows(
  lapply(
    names(groups),
    function(x) {
      monthly_summary(
        groups[[x]],
        x
      )
    }
  )
)

# 5. Checks
# Every group must have 12 months
monthly_data %>%
  count(Group)

# Forest September must have richness = 0 and abundance = 0
monthly_data %>%
  filter(
    Group == "fst_assemblage",
    Month == "Sep"
  )

# 6. Function that makes one rose diagram and its circular statistics
make_rose <- function(
    monthly_data,
    group_name,
    metric,
    panel_title,
    fill_color
) {

  # Select the group
  df <- monthly_data %>%
    filter(
      Group == group_name
    ) %>%
    arrange(Month)

  # Monthly values of richness or abundance
  values <- df[[metric]]

  # Angles in radians
  angles_rad <- month_angles * pi / 180

  # Sample size: richness = total species-month records; abundance = total individuals
  n_value <- sum(values)

  # Weighted circular components
  C <- sum(
    values * cos(angles_rad)
  )
  S <- sum(
    values * sin(angles_rad)
  )

  # Mean resultant length R (0 = evenly spread over the year, 1 = fully concentrated)
  R_value <- sqrt(
    C^2 + S^2
  ) / n_value

  # Mean direction in degrees (Feb = 0, ..., Jan = 330)
  mean_deg <- (
    atan2(S, C) *
      180 / pi
  ) %% 360

  # Rayleigh statistic
  Z_value <- n_value * R_value^2

  # Rayleigh test of uniformity (H0: no seasonality)
  weighted_angles <- rep(
    angles_rad,
    times = values
  )
  theta <- circular(
    weighted_angles,
    units = "radians"
  )
  rayleigh_result <- rayleigh.test(
    theta
  )
  p_value <- rayleigh_result$p.value

  # Month closest to the mean direction
  month_index <- (round(mean_deg / 30) %% 12) + 1
  
  mean_month <- month_levels[
    month_index
  ]

  # Month with the highest observed value
  peak_index <- which.max(
    values
  )
  peak_month <- month_levels[
    peak_index
  ]

  # Data for plotting
  plot_df <- data.frame(
    Month = month_levels,
    Value = values,
    Angle = plot_angles
  )

  # Each plot has its own radial scale (25% extra room for labels)
  radial_max <- max(values) * 1.25
  if (radial_max == 0) {
    radial_max <- 1
  }

  # Small offset so value labels sit just above the bars
  label_offset <- radial_max * 0.025
  plot_df <- plot_df %>%
    mutate(
      Label_y = Value + label_offset
    )

  # Mean vector: its length represents R
  vector_length <- radial_max * R_value

  # Plot angle of the vector (+15 degrees, plotting only)
  plot_mean_angle <- (
    mean_deg + 15
  ) %% 360
  vector_df <- data.frame(
    Mean_angle = plot_mean_angle,
    Vector_length = vector_length
  )

  # Rose diagram
  p <- ggplot(
    plot_df,
    aes(
      x = Angle,
      y = Value
    )
  ) +
  # Monthly bars
  geom_col(
    width = 30,
    fill = fill_color,
    color = fill_color,
    alpha = 0.80
  ) +
  # Monthly values
  geom_text(
    aes(
      label = Value,
      y = Label_y
    ),
    size = 8 / .pt,
    color = "gray20"
  ) +
  # Mean vector (red arrow)
  geom_segment(
    data = vector_df,
    aes(
      x = Mean_angle,
      xend = Mean_angle,
      y = 0,
      yend = Vector_length
    ),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 1.1,
    arrow = arrow(
      length = unit(
        0.20,
        "cm"
      ),
      type = "closed"
    )
  ) +
  # Polar coordinates (-15 degrees compensates for the +15 degree plotting shift)
  coord_polar(
    start = -pi / 2 - pi / 12,
    clip = "off"
  ) +
  # Month positions
  scale_x_continuous(
    breaks = plot_angles,
    labels = month_levels,
    limits = c(
      0,
      360
    ),
    expand = c(
      0,
      0
    )
  ) +
  # Independent radial scale
  scale_y_continuous(
    limits = c(
      0,
      radial_max
    ),
    expand = c(
      0,
      0
    )
  ) +
  labs(
    title = panel_title,
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
    theme(
      axis.text.x = element_text(
        size = 12,
        color = "gray20",
        margin = margin(t = 5)
      ),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_line(
        linewidth = 0.20,
        color = "gray70",
        linetype = "solid"
      ),
      plot.title = element_text(
        hjust = 0,
        face = "bold",
        size = 12
      ),
      plot.margin = margin(t = 2, r = 8, b = 2, l = 8, unit = "mm")
    )

  # Statistics of this group
  stats <- data.frame(
    Group = group_name,
    Metric = metric,
    n = n_value,
    R = R_value,
    Z = Z_value,
    p = p_value,
    Mean_direction = mean_deg,
    Mean_month = mean_month,
    Peak_month = peak_month
  )
  return(
    list(
      plot = p,
      stats = stats
    )
  )
}

# 7. Richness rose diagrams (a-c: Forest; d-f: Scrubland)
r1 <- make_rose(
  monthly_data,
  "fst_assemblage",
  "Richness",
  "a",
  "green4" 
)
r2 <- make_rose(
  monthly_data,
  "fst_solitary",
  "Richness",
  "b",
  "green4"
)
r3 <- make_rose(
  monthly_data,
  "fst_social",
  "Richness",
  "c",
  "green4"
)
r4 <- make_rose(
  monthly_data,
  "scb_assemblage",
  "Richness",
  "d",
  "orange"
)
r5 <- make_rose(
  monthly_data,
  "scb_solitary",
  "Richness",
  "e",
  "orange"
)
r6 <- make_rose(
  monthly_data,
  "scb_social",
  "Richness",
  "f",
  "orange"
)

# 8. Abundance rose diagrams (a-c: Forest; d-f: Scrubland)
a1 <- make_rose(
  monthly_data,
  "fst_assemblage",
  "Abundance",
  "a",
  "green4"
)
a2 <- make_rose(
  monthly_data,
  "fst_solitary",
  "Abundance",
  "b",
  "green4"
)
a3 <- make_rose(
  monthly_data,
  "fst_social",
  "Abundance",
  "c",
  "green4"
)
a4 <- make_rose(
  monthly_data,
  "scb_assemblage",
  "Abundance",
  "d",
  "orange"
)
a5 <- make_rose(
  monthly_data,
  "scb_solitary",
  "Abundance",
  "e",
  "orange"
)
a6 <- make_rose(
  monthly_data,
  "scb_social",
  "Abundance",
  "f",
  "orange"
)

# 9. Combine panels (2 rows x 3 columns)
richness_panel <- ggarrange(
  r1$plot,
  r2$plot,
  r3$plot,
  r4$plot,
  r5$plot,
  r6$plot,
  ncol = 3,
  nrow = 2,
  widths = c(
    1,
    1,
    1
  ),
  heights = c(
    1,
    1
  ),
  align = "hv"
)
richness_panel

# Abundance panel
abundance_panel <- ggarrange(
  a1$plot,
  a2$plot,
  a3$plot,
  a4$plot,
  a5$plot,
  a6$plot,
  ncol = 3,
  nrow = 2,
  widths = c(
    1,
    1,
    1
  ),
  heights = c(
    1,
    1
  ),
  align = "hv"
)
abundance_panel

# 10. Save figures
# Fig 4: richness
ggsave(
  filename =
    "Fig4.tiff",
  plot =
    richness_panel,
  path =
    fig_dir,
  width = 19,
  height = 12,
  units = "cm",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

# Fig 5: abundance
ggsave(
  filename =
    "Fig5.tiff",
  plot =
    abundance_panel,
  path =
    fig_dir,
  width = 19,
  height = 12,
  units = "cm",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

# 11. Statistics table (Table 3)
phenology_stats <- bind_rows(
  r1$stats,
  r2$stats,
  r3$stats,
  r4$stats,
  r5$stats,
  r6$stats,
  a1$stats,
  a2$stats,
  a3$stats,
  a4$stats,
  a5$stats,
  a6$stats
)

# Round values for the manuscript
phenology_stats_table <- phenology_stats %>%
  mutate(
    R = round(
      R,
      3
    ),
    Z = round(
      Z,
      3
    ),
    Mean_direction = round(
      Mean_direction,
      3
    ),
    p = ifelse(
      p < 0.001,
      "<0.001",
      sprintf(
        "%.3f",
        p
      )
    )
  )
phenology_stats_table

# Save the statistics
write.csv(
  phenology_stats_table,
  file.path(tab_dir, "phenology_richness_abundance.csv"),
  row.names = FALSE
)
