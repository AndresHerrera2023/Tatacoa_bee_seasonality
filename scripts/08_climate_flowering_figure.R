# ============================================================
# 08 - Climate and flowering figure
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Plot monthly precipitation, temperature and flowering plant richness for Forest and Scrubland.
# Input:   data/variables.csv
# Output:  figures/climate/S1_Fig.tiff
# ============================================================

# Packages
library(ggplot2)
library(gridExtra)
library(cowplot)    # align_plots() and plot_grid()

# 1. Paths and data
# Project root: the folder that contains data/ and figures/.
# Run this script from the project root (or change this path).
base_dir <- "."

fig_dir <- file.path(base_dir, "figures", "climate")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(
  file.path(base_dir, "data", "variables.csv"),
  stringsAsFactors = FALSE
)

# 2. Order the months: February 2018 -> January 2019
month_order <- c(
  "Feb", "Mar", "Apr", "May", "Jun", "Jul",
  "Aug", "Sep", "Oct", "Nov", "Dec", "Jan"
)

data$Month <- factor(
  data$Month,
  levels = month_order,
  ordered = TRUE
)

# 3. Separate the two vegetation types (called "Forest" and "Scrub" in the file)
forest <- subset(data, Site == "Forest")
scrub <- subset(data, Site == "Scrub")

# 4. Scales for the climographs
# Temperature is drawn on the precipitation axis and shown on a secondary axis (28-46 C)
pp_max <- 300
temp_min <- 28
temp_max <- 46

temp_to_pp <- function(x) {
  (x - temp_min) * pp_max / (temp_max - temp_min)
}

pp_to_temp <- function(x) {
  x * (temp_max - temp_min) / pp_max + temp_min
}

# 5. Panel a: Forest climograph (bars = precipitation, line = temperature)
pA <- ggplot(forest, aes(x = Month)) +
  geom_col(
    aes(y = PP),
    fill = "#4C9BE8",
    width = 0.7
  ) +
  geom_line(
    aes(y = temp_to_pp(TEMP), group = 1),
    color = "#D1495B",
    linewidth = 1
  ) +
  geom_point(
    aes(y = temp_to_pp(TEMP)),
    color = "#D1495B",
    size = 2.5
  ) +
  scale_y_continuous(
    name = "Precipitation (mm)",
    limits = c(0, 300),
    breaks = seq(0, 300, 50),
    sec.axis = sec_axis(
      ~ pp_to_temp(.),
      name = "Temperature (°C)",
      breaks = seq(28, 46, 2)
    )
  ) +
  labs(title = "a", x = NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold", size = 12, hjust = 0
    ),
    axis.title.y.left = element_text(
      face = "bold", size = 12
    ),
    axis.title.y.right = element_text(
      face = "bold", size = 12
    ),
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5,
      face = "plain", size = 10, color = "black"
    ),
    axis.text.y.left = element_text(
      size = 10, color = "black"
    ),
    axis.text.y.right = element_text(
      size = 10, color = "black"
    ),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# 6. Panel b: Scrubland climograph
pB <- ggplot(scrub, aes(x = Month)) +
  geom_col(
    aes(y = PP),
    fill = "#4C9BE8",
    width = 0.7
  ) +
  geom_line(
    aes(y = temp_to_pp(TEMP), group = 1),
    color = "#D1495B",
    linewidth = 1
  ) +
  geom_point(
    aes(y = temp_to_pp(TEMP)),
    color = "#D1495B",
    size = 2.5
  ) +
  scale_y_continuous(
    name = "Precipitation (mm)",
    limits = c(0, 300),
    breaks = seq(0, 300, 50),
    sec.axis = sec_axis(
      ~ pp_to_temp(.),
      name = "Temperature (°C)",
      breaks = seq(28, 46, 2)
    )
  ) +
  labs(title = "b", x = NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold", size = 12, hjust = 0
    ),
    axis.title.y.left = element_text(
      face = "bold", size = 12
    ),
    axis.title.y.right = element_text(
      face = "bold", size = 12
    ),
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5,
      face = "plain", size = 10, color = "black"
    ),
    axis.text.y.left = element_text(
      size = 10, color = "black"
    ),
    axis.text.y.right = element_text(
      size = 10, color = "black"
    ),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# 7. Panel c: Forest, flowering plant richness per month
pC <- ggplot(forest, aes(x = Month, y = FLW)) +
  geom_col(
    fill = "#5AA469",
    width = 0.7
  ) +
  scale_y_continuous(
    limits = c(0, 16),
    breaks = seq(0, 16, 2)
  ) +
  labs(
    title = "c",
    x = NULL,
    y = "Flowering plant richness"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold", size = 12, hjust = 0
    ),
    axis.title.y = element_text(
      face = "bold", size = 12
    ),
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5,
      face = "plain", size = 10, color = "black"
    ),
    axis.text.y = element_text(
      size = 10, color = "black"
    ),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# 8. Panel d: Scrubland, flowering plant richness per month
pD <- ggplot(scrub, aes(x = Month, y = FLW)) +
  geom_col(
    fill = "#5AA469",
    width = 0.7
  ) +
  scale_y_continuous(
    limits = c(0, 16),
    breaks = seq(0, 16, 2)
  ) +
  labs(
    title = "d",
    x = NULL,
    y = "Flowering plant richness"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold", size = 12, hjust = 0
    ),
    axis.title.y = element_text(
      face = "bold", size = 12
    ),
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5,
      face = "plain", size = 10, color = "black"
    ),
    axis.text.y = element_text(
      size = 10, color = "black"
    ),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )


# 9. Align and combine the four panels
# Align plot areas
aligned_plots <- align_plots(
  pA, pB, pC, pD,
  align = "hv",
  axis = "tblr"
)

# Combine (2 columns x 2 rows)
final_figure <- plot_grid(
  plotlist = aligned_plots,
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

final_figure

# 10. Save figure (S1 Fig)
ggsave(
  filename = file.path(fig_dir, "S1_Fig.tiff"),
  plot = final_figure,
  width = 19,
  height = 15,
  units = "cm",
  compression = "lzw",
  dpi = 300
)
