# ============================================================
# 03 - Rank-abundance curves
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Show dominance and evenness of Forest and Scrubland assemblages, overall and for dry and wet seasons.
# Input:   data/database.csv
# Output:  figures/rank/S3_Fig.tif
# ============================================================

# Packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Helper: clean species names for labels
# (underscores -> spaces, subgenus in parentheses removed,
# e.g. "Ceratina (Zadontomerus) nautlana" -> "Ceratina nautlana")
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\s*\\([^)]*\\)", "", x)
  trimws(gsub("\\s+", " ", x))
}

# 1. Paths and data
# Project root: the folder that contains data/ and figures/.
# Run this script from the project root (or change this path).
base_dir <- "."

fig_dir <- file.path(base_dir, "figures", "rank")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(
  file.path(base_dir, "data", "database.csv")
)

# 2. Overall rank-abundance: Forest vs Scrubland
# Rank species from the most to the least abundant within each vegetation type
rank_overall <- data %>%
  group_by(Area, Genus_species) %>%
  summarise(Abundance = n(), .groups = "drop") %>%
  group_by(Area) %>%
  arrange(desc(Abundance), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# Keep the 3 most abundant species of each vegetation type for labels
top3_overall <- rank_overall %>%
  group_by(Area) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "italic('",
      clean_name(Genus_species),
      "')"
    )
  )

# 3. Seasonal rank-abundance (dry vs wet months)
wet_months <- c("Mar", "Apr", "May", "Oct", "Nov")
dry_months <- c("Feb", "Jun", "Jul", "Aug", "Sep", "Dec", "Jan")

data_season <- data %>%
  mutate(
    Season = case_when(
      Month %in% wet_months ~ "Wet",
      Month %in% dry_months ~ "Dry",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Season))

# Rank species within each vegetation type and season
rank_season <- data_season %>%
  group_by(Area, Season, Genus_species) %>%
  summarise(Abundance = n(), .groups = "drop") %>%
  group_by(Area, Season) %>%
  arrange(desc(Abundance), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# Keep the 3 most abundant species of each type and season for labels
top3_labels <- rank_season %>%
  group_by(Area, Season) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "italic('",
      clean_name(Genus_species),
      "')"
    )
  )

# Split by vegetation type
rank_forest_season <- rank_season %>%
  filter(Area == "Forest")

rank_scrubland_season <- rank_season %>%
  filter(Area == "Scrubs")

top3_forest <- top3_labels %>%
  filter(Area == "Forest")

top3_scrubland <- top3_labels %>%
  filter(Area == "Scrubs")

# 4. Common axis limits for all panels
x_max_all <- max(
  c(rank_overall$Rank, rank_season$Rank),
  na.rm = TRUE
)

y_max_all <- max(
  c(rank_overall$Abundance, rank_season$Abundance),
  na.rm = TRUE
)

# 5. Panel a: Forest vs Scrubland
p_rank_overall <- ggplot(
  rank_overall,
  aes(x = Rank, y = Abundance, color = Area)
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_text_repel(
    data = top3_overall,
    aes(label = label),
    parse = TRUE,
    size = 8 / .pt,
    color = "black",
    hjust = 0,
    nudge_x = 14,
    direction = "y",
    segment.size = 0.3,
    segment.color = "grey40",
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Forest" = "darkgreen",
      "Scrubs" = "orange3"
    ),
    labels = c(
      "Forest" = "Forest",
      "Scrubs" = "Scrubland"
    )
  ) +
  scale_x_continuous(
    limits = c(0, x_max_all)
  ) +
  scale_y_continuous(
    limits = c(0, y_max_all)
  ) +
  labs(
    title = "a",
    x = "Species rank",
    y = "Abundance",
    color = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0,
      color = "black"
    ),
    legend.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "bottom",
    legend.margin = margin(t = -4)
  )

# 6. Panel b: Forest, dry vs wet
p_rank_forest_season <- ggplot(
  rank_forest_season,
  aes(
    x = Rank,
    y = Abundance,
    color = Season,
    linetype = Season
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_text_repel(
    data = top3_forest,
    aes(label = label),
    parse = TRUE,
    size = 8 / .pt,
    color = "black",
    hjust = 0,
    nudge_x = 14,
    direction = "y",
    segment.size = 0.3,
    segment.color = "grey40",
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Dry" = "red3",
      "Wet" = "steelblue"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Dry" = "dashed",
      "Wet" = "solid"
    )
  ) +
  scale_x_continuous(
    limits = c(0, x_max_all)
  ) +
  scale_y_continuous(
    limits = c(0, y_max_all)
  ) +
  labs(
    title = "b",
    x = "Species rank",
    y = "Abundance",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0,
      color = "black"
    ),
    legend.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "bottom",
    legend.margin = margin(t = -4)
  )

# 7. Panel c: Scrubland, dry vs wet
p_rank_scrubland_season <- ggplot(
  rank_scrubland_season,
  aes(
    x = Rank,
    y = Abundance,
    color = Season,
    linetype = Season
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_text_repel(
    data = top3_scrubland,
    aes(label = label),
    parse = TRUE,
    size = 8 / .pt,
    color = "black",
    hjust = 0,
    nudge_x = 14,
    direction = "y",
    segment.size = 0.3,
    segment.color = "grey40",
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Dry" = "red3",
      "Wet" = "steelblue"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Dry" = "dashed",
      "Wet" = "solid"
    )
  ) +
  scale_x_continuous(
    limits = c(0, x_max_all)
  ) +
  scale_y_continuous(
    limits = c(0, y_max_all)
  ) +
  labs(
    title = "c",
    x = "Species rank",
    y = "Abundance",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0,
      color = "black"
    ),
    legend.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "bottom",
    legend.margin = margin(t = -4)
  )

# 8. Combine panels (each panel has its own legend)
final_rank_plot <- ggarrange(
  p_rank_overall,
  p_rank_forest_season,
  p_rank_scrubland_season,
  ncol = 3,
  nrow = 1,
  align = "hv"
)

final_rank_plot

# 9. Save figure (S3 Fig)
ggsave(
  filename = file.path(fig_dir, "S3_Fig.tif"),
  plot = final_rank_plot,
  width = 19,
  height = 9,
  units = "cm",
  dpi = 300,
  bg = "white",
  device = "tiff",
  compression = "lzw"
)
