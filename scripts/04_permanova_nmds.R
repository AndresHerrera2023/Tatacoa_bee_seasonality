# ============================================================
# 04 - Community composition: PERMANOVA, PERMDISP and NMDS
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Test whether bee assemblage composition differs between vegetation types (Forest vs Scrubland) and seasons (dry vs wet), and visualize it with NMDS. Hellinger-transformed abundances and Euclidean distances.
# Input:   data/database.csv
# Output:  figures/NMDS/Fig3.tiff, tables/permanova/*.csv
# ============================================================

# Packages
library(vegan)
library(permute)
library(ggplot2)
library(ggpubr)
library(patchwork)   # install.packages("patchwork") if missing

# Paths
# Project root: the folder that contains data/, figures/ and tables/.
# Run this script from the project root (or change this path).
base_dir <- "."
fig_dir  <- file.path(base_dir, "figures", "NMDS")
tab_dir  <- file.path(base_dir, "tables", "permanova")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# Package versions (to report in the Methods)
packageVersion("vegan")
packageVersion("permute")

# 1. Load data and build the community matrix
# (rows = monthly sample of each site, columns = species, values = abundance)
data <- read.csv(
  file.path(base_dir, "data", "database.csv"),
  stringsAsFactors = FALSE
)

# Sampling months in chronological order
all_months <- c(
  "Feb","Mar","Apr","May","Jun","Jul",
  "Aug","Sep","Oct","Nov","Dec","Jan"
)

all_species <- sort(unique(na.omit(data$Genus_species)))

# Separate the two vegetation types (called "Scrubs" in the database)
forest <- data[data$Area == "Forest", ]
scrub  <- data[data$Area == "Scrubs", ]

forest_t <- table(
  factor(forest$Month, levels = all_months),
  factor(forest$Genus_species, levels = all_species)
)

scrub_t <- table(
  factor(scrub$Month, levels = all_months),
  factor(scrub$Genus_species, levels = all_species)
)

# Combine Forest and Scrubland into one matrix
assemblage_matrix <- as.matrix(rbind(forest_t, scrub_t))

# Metadata in exactly the same row order as the matrix
metadata <- data.frame(
  Site  = rep(c("Forest", "Scrub"), each = 12),
  Month = rep(all_months, times = 2)
)

# Wet months: March-May and October-November; all other months are dry
metadata$Season <- ifelse(
  metadata$Month %in% c("Mar","Apr","May","Oct","Nov"),
  "Wet", "Dry"
)

metadata$Site   <- factor(metadata$Site,   levels = c("Forest", "Scrub"))
metadata$Month  <- factor(metadata$Month,  levels = all_months)
metadata$Season <- factor(metadata$Season, levels = c("Dry", "Wet"))

# Unique row names (site_month)
rownames(assemblage_matrix) <- paste(metadata$Site, metadata$Month, sep = "_")
rownames(metadata)          <- rownames(assemblage_matrix)

# 2. Checks: individuals and richness per month
dim(assemblage_matrix)
nrow(metadata)

# Individuals (N) and richness per sample
ab <- data.frame(
  metadata,
  N        = rowSums(assemblage_matrix),
  Richness = rowSums(assemblage_matrix > 0)
)

# Months x sites table
ab_w <- reshape(
  ab[, c("Site", "Month", "N")],
  idvar = "Month", timevar = "Site", direction = "wide"
)

reshape(
  ab[, c("Site", "Month", "N", "Richness")],
  idvar = "Month", timevar = "Site", direction = "wide"
)

# Forest September had zero bees. It is valid for richness and abundance,
# but a sample without individuals has no relative composition, so it is
# excluded only from the composition analyses.
assemblage_matrix[metadata$Site == "Forest" & metadata$Month == "Sep", ]

# Keep only samples with at least one individual
keep <- rowSums(assemblage_matrix) > 0

assemblage_comp <- assemblage_matrix[keep, ]
metadata_comp   <- droplevels(metadata[keep, ])

# Check: 11 Forest + 12 Scrubland samples
table(metadata_comp$Site)                    # 11 Forest + 12 Scrub
table(metadata_comp$Season, metadata_comp$Site)

# 3. Helper functions
# Months with at least min_n individuals in BOTH sites
months_min <- function(min_n) {
  as.character(ab_w$Month[ab_w$N.Forest >= min_n & ab_w$N.Scrub >= min_n])
}

# Number of permutations actually used by an adonis2 result
n_perm <- function(a) {
  tryCatch(permute::getNperm(attr(a, "control")), error = function(e) NA_integer_)
}

# Q1 and Q2 on a balanced set of months (each month present in both sites)
paired_tests <- function(months) {
  
  # Keep the selected months
  idx <- metadata$Month %in% months
  m   <- droplevels(metadata[idx, ])
  cm  <- assemblage_matrix[idx, ]
  
  # Order by month, then site (each month is one contiguous pair)
  o  <- order(m$Month, m$Site)
  m  <- m[o, ]
  cm <- cm[o, ]
  stopifnot(all(table(m$Month, m$Site) == 1))
  
  # Hellinger transformation + Euclidean distance
  d <- vegdist(decostand(cm, "hellinger"), "euclidean")
  
  # Site test: permute site labels only within months (paired design)
  perm_site <- how(nperm = 999, blocks = m$Month)
  
  # Season test: permute whole months (both sites together)
  perm_season <- how(
    nperm  = 999,
    plots  = Plots(strata = m$Month, type = "free"),
    within = Within(type = "none")
  )
  
  # PERMANOVA: vegetation type (site)
  set.seed(123)
  a_site <- adonis2(d ~ Site, data = m, permutations = perm_site)
  
  # PERMDISP: are dispersions similar between sites?
  set.seed(123)
  pd_site <- permutest(betadisper(d, m$Site, type = "median"),
                       permutations = perm_site)
  
  # PERMANOVA: season, adjusted for site (marginal effect)
  set.seed(123)
  a_season <- adonis2(d ~ Site + Season, data = m,
                      permutations = perm_season, by = "margin")
  
  # PERMDISP: are dispersions similar between seasons?
  set.seed(123)
  pd_season <- permutest(betadisper(d, m$Season, type = "median"),
                         permutations = perm_season)
  
  list(months = levels(m$Month), meta = m, dist = d,
       site = a_site, permdisp_site = pd_site,
       season = a_season, permdisp_season = pd_season)
}

# Q3: season within one vegetation type (all its months with individuals)
within_site <- function(site) {
  
  # Select the samples of this site
  k <- metadata_comp$Site == site
  m <- metadata_comp[k, ]
  d <- vegdist(decostand(assemblage_comp[k, ], "hellinger"), "euclidean")
  
  # Free permutations
  set.seed(123)
  a <- adonis2(d ~ Season, data = m, permutations = 999)
  
  # PERMDISP
  set.seed(123)
  pd <- permutest(betadisper(d, m$Season, type = "median"), permutations = 999)
  
  # Series permutations (cyclic shifts of the monthly series).
  # There are few distinct permutations, so this is a conservative check
  # for temporal autocorrelation.
  set.seed(123)
  a_ser <- adonis2(d ~ Season, data = m,
                   permutations = how(within = Within(type = "series"),
                                      nperm = 999))
  
  list(meta = m, dist = d, free = a, permdisp = pd, series = a_ser)
}

# 4. Q1 and Q2: Forest vs Scrubland, dry vs wet (all months present in both sites)
# September is removed in both sites because Forest had no bees.
# Q1: read the Site row of res_main$site (and PERMDISP).
# Q2: read ONLY the Season row of res_main$season;
#     the Site row is not valid under whole-month permutations.
res_main <- paired_tests(months_min(1))
res_main$months

res_main$site
res_main$permdisp_site

res_main$season
res_main$permdisp_season

# 5. Q3: dry vs wet within each vegetation type
res_forest <- within_site("Forest")
res_scrub  <- within_site("Scrub")

res_forest$free;  res_forest$permdisp;  res_forest$series
res_scrub$free;   res_scrub$permdisp;   res_scrub$series

# 6. Sensitivity: exclude months with few individuals
# (>= 15 individuals: 10 months; >= 20 individuals: 7 months, in both sites)
res_15 <- paired_tests(months_min(15))
res_20 <- paired_tests(months_min(20))

res_15$months; res_15$site; res_15$season
res_20$months; res_20$site; res_20$season

# 7. Supplementary table (S4 Table), built directly from the results
tab_row <- function(a, term, analysis, factor_name, test, n, pd = NULL) {
  r <- if (term %in% rownames(a)) a[term, ] else a[1, ]
  data.frame(
    Analysis     = analysis,
    Factor       = factor_name,
    Test         = test,
    N_samples    = n,
    Df           = r$Df,
    SumOfSqs     = round(r$SumOfSqs, 4),
    MeanSqs      = round(r$SumOfSqs / r$Df, 4),
    F            = round(r$F, 4),
    R2           = round(r$R2, 4),
    p            = round(r$`Pr(>F)`, 3),
    PERMDISP_F   = if (is.null(pd)) NA else round(pd$tab[1, "F"], 4),
    PERMDISP_p   = if (is.null(pd)) NA else round(pd$tab[1, "Pr(>F)"], 3),
    Permutations = n_perm(a),
    row.names    = NULL
  )
}

# One row per test
permanova_supp <- rbind(
  tab_row(res_main$site, "Site", "Overall", "Site",
          "Paired by month", nrow(res_main$meta), res_main$permdisp_site),
  tab_row(res_main$season, "Season", "Overall", "Season",
          "Whole-month permutation", nrow(res_main$meta), res_main$permdisp_season),
  tab_row(res_forest$free, "Season", "Forest", "Season",
          "Free permutation", nrow(res_forest$meta), res_forest$permdisp),
  tab_row(res_scrub$free, "Season", "Scrubland", "Season",
          "Free permutation", nrow(res_scrub$meta), res_scrub$permdisp),
  tab_row(res_forest$series, "Season", "Forest", "Season",
          "Series permutation", nrow(res_forest$meta)),
  tab_row(res_scrub$series, "Season", "Scrubland", "Season",
          "Series permutation", nrow(res_scrub$meta)),
  tab_row(res_15$site, "Site", "Sensitivity (>= 15 ind.)", "Site",
          "Paired by month", nrow(res_15$meta)),
  tab_row(res_15$season, "Season", "Sensitivity (>= 15 ind.)", "Season",
          "Whole-month permutation", nrow(res_15$meta)),
  tab_row(res_20$site, "Site", "Sensitivity (>= 20 ind.)", "Site",
          "Paired by month", nrow(res_20$meta)),
  tab_row(res_20$season, "Season", "Sensitivity (>= 20 ind.)", "Season",
          "Whole-month permutation", nrow(res_20$meta))
)

permanova_supp

# Save the table
write.csv(
  permanova_supp,
  file.path(tab_dir, "supplementary_PERMANOVA_PERMDISP.csv"),
  row.names = FALSE
)

# 8. NMDS (k = 2). NMDS is only a visual aid; inference comes from PERMANOVA.
# Panel a: both sites (same samples as Q1 and Q2)
set.seed(123)
nmds_sites <- metaMDS(res_main$dist, k = 2, trymax = 100, trace = FALSE)

# Panels b and c: Forest and Scrubland separately
set.seed(123)
nmds_forest <- metaMDS(res_forest$dist, k = 2, trymax = 100, trace = FALSE)

set.seed(123)
nmds_scrub <- metaMDS(res_scrub$dist, k = 2, trymax = 100, trace = FALSE)

# Stress values (fit of each ordination)
c(sites  = nmds_sites$stress,
  forest = nmds_forest$stress,
  scrub  = nmds_scrub$stress)

# Site scores for plotting
scores_sites <- as.data.frame(scores(nmds_sites, display = "sites"))
scores_sites$Site   <- factor(ifelse(res_main$meta$Site == "Scrub", "Scrubland", "Forest"),
                              levels = c("Forest", "Scrubland"))
scores_sites$Month  <- as.character(res_main$meta$Month)
scores_sites$Season <- res_main$meta$Season

scores_forest <- as.data.frame(scores(nmds_forest, display = "sites"))
scores_forest$Month  <- as.character(res_forest$meta$Month)
scores_forest$Season <- res_forest$meta$Season

scores_scrub <- as.data.frame(scores(nmds_scrub, display = "sites"))
scores_scrub$Month  <- as.character(res_scrub$meta$Month)
scores_scrub$Season <- res_scrub$meta$Season

# Convex hulls (polygons around each group)
find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

sites_hull  <- do.call(rbind, lapply(split(scores_sites,  scores_sites$Site),    find_hull))
forest_hull <- do.call(rbind, lapply(split(scores_forest, scores_forest$Season), find_hull))
scrub_hull  <- do.call(rbind, lapply(split(scores_scrub,  scores_scrub$Season),  find_hull))

# Common theme and colours
theme_nmds <- theme_bw() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_text(size = 10, color = "black"),
    axis.title   = element_text(size = 12, color = "black", face = "bold"),
    plot.title   = element_text(hjust = 0, color = "black", size = 12, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.title = element_blank()
  )

season_cols <- c("Dry" = "#e41a1c", "Wet" = "#377eb8")
site_cols   <- c("Forest" = "#1b7837", "Scrubland" = "#d95f02")

# Extra room around the points so month labels are not cut at the panel edge
axes_room <- list(
  scale_x_continuous(expand = expansion(mult = 0.08),
                     breaks = scales::breaks_width(0.5),
                     labels = scales::label_number(accuracy = 0.1)),
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.14)),
                     breaks = scales::breaks_width(0.2),
                     labels = scales::label_number(accuracy = 0.1))
)

# Panel a: both sites (colour and polygon = site)
NMDS_sites <- ggplot(scores_sites, aes(NMDS1, NMDS2)) +
  geom_polygon(data = sites_hull, aes(fill = Site, group = Site), alpha = 0.2) +
  geom_point(aes(color = Site), shape = 16, size = 2) +
  geom_text(aes(label = Month, color = Site),
            vjust = -0.7, size = 8 / .pt, show.legend = FALSE) +
  scale_fill_manual(values  = site_cols) +
  scale_color_manual(values = site_cols) +
  axes_room +
  theme_nmds +
  ggtitle("a")

# Panel b: Forest (colour and polygon = season)
NMDS_forest <- ggplot(scores_forest, aes(NMDS1, NMDS2)) +
  geom_polygon(data = forest_hull, aes(fill = Season, group = Season), alpha = 0.3) +
  geom_point(aes(color = Season), size = 1.5) +
  geom_text(aes(label = Month, color = Season),
            vjust = -0.6, size = 8 / .pt, show.legend = FALSE) +
  scale_fill_manual(values  = season_cols) +
  scale_color_manual(values = season_cols) +
  axes_room +
  theme_nmds +
  ggtitle("b")

# Panel c: Scrubland (colour and polygon = season)
NMDS_scrub <- ggplot(scores_scrub, aes(NMDS1, NMDS2)) +
  geom_polygon(data = scrub_hull, aes(fill = Season, group = Season), alpha = 0.3) +
  geom_point(aes(color = Season), size = 1.5) +
  geom_text(aes(label = Month, color = Season),
            vjust = -0.6, size = 8 / .pt, show.legend = FALSE) +
  scale_fill_manual(values  = season_cols) +
  scale_color_manual(values = season_cols) +
  axes_room +
  theme_nmds +
  ggtitle("c")

# Combine panels a | b | c with one common legend at the bottom (Fig 3)
one_row <- guides(fill  = guide_legend(nrow = 1),
                  color = guide_legend(nrow = 1))

NMDS_fig3 <- wrap_plots(NMDS_sites + one_row,
                        NMDS_forest + one_row,
                        NMDS_scrub + one_row,
                        nrow = 1, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box = "horizontal",
                                legend.spacing.x = unit(0.6, "cm")))

NMDS_fig3

# Save figure (Fig 3)
ggsave(
  filename = "Fig3.tiff",
  plot = NMDS_fig3,
  path = fig_dir,
  width = 19, height = 10.5, units = "cm",
  dpi = 300, compression = "lzw", bg = "white"
)

# 9. Results to report
res_main$site;   res_main$permdisp_site
res_main$season; res_main$permdisp_season
res_forest$free; res_forest$permdisp; res_forest$series
res_scrub$free;  res_scrub$permdisp;  res_scrub$series
permanova_supp

c(sites  = nmds_sites$stress,
  forest = nmds_forest$stress,
  scrub  = nmds_scrub$stress)

# Save session information
writeLines(capture.output(sessionInfo()),
           file.path(tab_dir, "sessionInfo.txt"))

# 10. Short version of the table (as in the manuscript)
tab <- permanova_supp
tab$Factor   <- sub("^Site$", "Veg. type", tab$Factor)
tab$Test     <- sub("Paired by month", "Paired", tab$Test)
tab$Test     <- sub("Whole-month permutation", "Whole-month", tab$Test)
tab$Test     <- sub("Free permutation", "Free", tab$Test)
tab$Test     <- sub("Series permutation", "Series", tab$Test)
tab$Analysis <- sub("Sensitivity \\(>= (\\d+) ind\\.\\)", "Sens. >=\\1", tab$Analysis)

num <- sapply(tab, is.numeric)
tab[num] <- lapply(tab[num], round, 3)
names(tab) <- c("Analysis","Factor","Test","n","Df","SS","MS","F","R2","p",
                "F_disp","p_disp","Perm")

tab
write.csv(tab, file.path(tab_dir, "supplementary_PERMANOVA_PERMDISP_short.csv"),
          row.names = FALSE, na = "-")
