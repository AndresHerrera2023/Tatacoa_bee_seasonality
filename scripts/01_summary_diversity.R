# ============================================================
# 01 - Diversity summary
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Count tribes, genera, species and individuals (overall, Forest, Scrubland) and monthly richness, abundance and Pielou evenness.
# Input:   data/database.csv (one row per individual bee)
# Output:  tables/summary_diversity/*.csv
# ============================================================

# Packages
library(dplyr)
library(vegan)   # diversity() for Shannon index

# 1. Paths and data
# Project root: the folder that contains data/ and tables/.
# Run this script from the project root (or change this path).
base_dir <- "."

out_dir <- file.path(base_dir, "tables", "summary_diversity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(file.path(base_dir, "data", "database.csv"))
head(data)

# 2. Totals: tribes, genera, species and individuals (all data)
summary_all <- data %>%
  summarise(Tribes = n_distinct(Tribe), 
            Genera = n_distinct(Genus), 
            Species = n_distinct(Genus_species),
            Individuals = n())
print(summary_all) 

# Totals by family (all data)
summary_family_all <- data %>%
  group_by(Family) %>% 
  summarise(Tribes = n_distinct(Tribe), 
            Genera = n_distinct(Genus), 
            Species = n_distinct(Genus_species),
            Individuals = n()) %>%
  add_row(Family = "Total", summary_all)
print(summary_family_all)

# Forest only: totals
summary_forest <- data[data$Area == "Forest",] %>%
  summarise(Tribes = n_distinct(Tribe), 
            Genera = n_distinct(Genus), 
            Species = n_distinct(Genus_species),
            Individuals = n())
print(summary_forest) 

# Forest only: totals by family
summary_family_forest <- data[data$Area == "Forest",] %>%
  group_by(Family) %>% 
  summarise(Tribes = n_distinct(Tribe), 
            Genera = n_distinct(Genus), 
            Species = n_distinct(Genus_species),
            Individuals = n())

summary_family_forest <- summary_family_forest %>% 

  # Colletidae was not recorded in Forest: add it as a row of zeros
  add_row(Family = "Colletidae",Tribes = 0, Genera = 0, Species = 0, .before = 2) 

# Replace NA values with 0
summary_family_forest[is.na(summary_family_forest)] <- 0

summary_family_forest <- summary_family_forest %>%

  # Add the total row
  add_row(Family = "Total", summary_forest)
print(summary_family_forest)

# Scrubland only: totals and totals by family
summary_scrub <- data[data$Area == "Scrubs",] %>%
  summarise(Tribes = n_distinct(Tribe), 
            Genera = n_distinct(Genus), 
            Species = n_distinct(Genus_species),
            Individuals = n())
print(summary_scrub) 

summary_family_scrub <- data[data$Area == "Scrubs",] %>%
  group_by(Family) %>% 
  summarise(Tribes = n_distinct(Tribe), 
            Genera = n_distinct(Genus), 
            Species = n_distinct(Genus_species),
            Individuals = n()) %>%
  add_row(Family = "Total", summary_scrub)
print(summary_family_scrub)

# 3. Final table: all data, Forest and Scrubland side by side
summarized_table <- cbind(summary_family_all, 
                          summary_family_forest[,2:5],
                          summary_family_scrub[,2:5])
print(summarized_table)
write.csv(summarized_table, file.path(out_dir, "summary_diversity.csv"), row.names = FALSE)


# 4. Monthly richness, abundance and Pielou evenness
all_months <- c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan")

# Forest
summary_month_forest <- data %>%
  filter(Area == "Forest") %>%
  mutate(Month = factor(Month, levels = all_months)) %>%
  group_by(Month) %>%
  summarise(
    Richness = n_distinct(Genus_species),
    Abundance = n(),
    Pielou = {
      abund <- table(Genus_species)
      S <- length(abund)
      ifelse(S > 1, diversity(abund, index = "shannon") / log(S), NA)
    },
    .groups = "drop"
  ) %>%
  mutate(Pielou = round(Pielou, 3))

print(summary_month_forest)

# Scrubland
summary_month_scrub <- data %>%
  filter(Area == "Scrubs") %>%
  mutate(Month = factor(Month, levels = all_months)) %>%
  group_by(Month) %>%
  summarise(
    Richness = n_distinct(Genus_species),
    Abundance = n(),
    Pielou = {
      abund <- table(Genus_species)
      S <- length(abund)
      ifelse(S > 1, diversity(abund, index = "shannon") / log(S), NA)
    },
    .groups = "drop"
  ) %>%
  mutate(Pielou = round(Pielou, 3))

print(summary_month_scrub)

# Forest and Scrubland in one table
summary_monthly <- full_join(
  summary_month_forest,
  summary_month_scrub,
  by = "Month",
  suffix = c("_Forest", "_Scrubs")
)

print(summary_monthly)

write.csv(
  summary_monthly,
  file.path(out_dir, "summary_monthly_diversity.csv"),
  row.names = FALSE
)
