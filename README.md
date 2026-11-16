# Seasonal and monthly dynamics of bee assemblages in a Colombian tropical dry forest

## Overview

This repository contains the data and R scripts associated with a study of seasonal and monthly bee assemblage dynamics in two vegetation types (Forest and Scrubland) in La Tatacoa, Huila, Colombia. Monthly sampling covered February 2018 through January 2019. Analyses examine bee diversity, community composition, phenology, temporal beta diversity, and relationships with flowering plant richness and climatic variables.

**Associated manuscript:** *Seasonal and monthly dynamics of bee assemblages are associated with flowering plant richness in an inter-Andean tropical dry forest in Colombia*.

## Repository structure

```text
repository/
├── README.md
├── data/
│   ├── database.csv
│   └── variables.csv
└── scripts/
    ├── 01_summary_diversity.R
    ├── 02_rarefaction_extrapolation.R
    ├── 03_rank_abundance.R
    ├── 04_permanova_nmds.R
    ├── 05_circular_rose_diagrams.R
    ├── 06_temporal_beta_diversity.R
    ├── 07_glm_environmental_drivers.R
    └── 08_climate_flowering_figure.R
```

## Data files

- **`data/database.csv`**: individual bee collection records, including sampling date, vegetation type, taxonomic identification, sociality, and associated collection information.
- **`data/variables.csv`**: monthly observations by vegetation type, including temperature (`TEMP`), relative humidity (`RH`), flowering plant richness (`FLW`), precipitation (`PP`), bee richness and abundance, and corresponding solitary/social bee summaries.

**Before making the repository public:** `database.csv` contains geographic coordinates and collection locality information. Review whether those records can be shared openly, especially for sensitive taxa or sites, and check any applicable permissions or restrictions.

## R scripts

Run the scripts from the **repository root** (the folder containing `data/` and `scripts/`). The scripts set `base_dir <- "."` and create output directories when needed.

| Script | Main task |
|---|---|
| `01_summary_diversity.R` | Summary of taxonomic diversity, monthly richness, abundance, and Pielou evenness |
| `02_rarefaction_extrapolation.R` | Incidence-based rarefaction and extrapolation using iNEXT |
| `03_rank_abundance.R` | Rank-abundance curves for vegetation types and seasons |
| `04_permanova_nmds.R` | PERMANOVA, PERMDISP, and NMDS community-composition analyses |
| `05_circular_rose_diagrams.R` | Circular phenology statistics and rose diagrams |
| `06_temporal_beta_diversity.R` | Temporal beta-diversity analyses and decomposition |
| `07_glm_environmental_drivers.R` | GLMs relating bee richness/abundance to flowering, precipitation, and temperature |
| `08_climate_flowering_figure.R` | Figure summarizing monthly precipitation, temperature, and flowering |

## Requirements

R and the following packages are used across the scripts:

`AER`, `betapart`, `circular`, `cowplot`, `dplyr`, `ggplot2`, `ggpubr`, `ggrepel`, `gridExtra`, `iNEXT`, `MASS`, `patchwork`, `permute`, `tidyr`, `vegan`.

Install packages that are not already present, for example:

```r
packages <- c("AER", "betapart", "circular", "cowplot", "dplyr", "ggplot2",
              "ggpubr", "ggrepel", "gridExtra", "iNEXT", "MASS",
              "patchwork", "permute", "tidyr", "vegan")
install.packages(setdiff(packages, rownames(installed.packages())))
```

`grid` is included with base R and does not need separate installation.

## Running the analyses

1. Download or clone the repository.
2. Open R or RStudio and set the working directory to the repository root.
3. Install the required packages.
4. Execute the scripts in numeric order, for example:

```r
source("scripts/01_summary_diversity.R")
source("scripts/02_rarefaction_extrapolation.R")
source("scripts/03_rank_abundance.R")
source("scripts/04_permanova_nmds.R")
source("scripts/05_circular_rose_diagrams.R")
source("scripts/06_temporal_beta_diversity.R")
source("scripts/07_glm_environmental_drivers.R")
source("scripts/08_climate_flowering_figure.R")
```

Generated figures and tables are written to the `figures/` and `tables/` directories, respectively. These are created by the scripts as necessary. Some scripts use `set.seed(123)` for stochastic analyses.


