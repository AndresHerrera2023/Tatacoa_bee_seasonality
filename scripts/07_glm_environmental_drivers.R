# ============================================================
# 07 - GLMs: bee richness and abundance vs. flowering, precipitation and temperature
# Project: Bee assemblages in a tropical dry forest, La Tatacoa, Colombia
#          (sampling: February 2018 - January 2019)
#
# Purpose: Test whether monthly bee richness and abundance are associated with flowering plant richness (FLW), precipitation (PP) and temperature (TEMP). 12 models = 2 vegetation types x (overall, social, solitary) x (richness, abundance); n = 12 months per model.
# Input:   data/variables.csv
# Output:  tables/glm/*.csv and *.txt (Tables 4, 5 and S5)
#
# Rule for the error distribution (fixed in advance, same for all models):
#   1. Fit Poisson and negative binomial (NB).
#   2. Dispersion test (AER) p < 0.05                  -> NB
#   3. p < 0.10 OR Pearson ratio > 1.5 (equivocal)     -> model with lower AIC
#   4. otherwise                                       -> Poisson
# Extra checks: VIF, lag-1 autocorrelation of residuals, NB theta and convergence,
# Holm-adjusted p-values over all 36 coefficient tests.
# ============================================================

# Packages
# (MASS: negative binomial models; AER: dispersion test)
library(MASS)
library(AER)

# Paths
# Project root: the folder that contains data/ and tables/.
# Run this script from the project root (or change this path).
base_dir <- "."
out_dir  <- file.path(base_dir, "tables", "glm")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Data
# One row per month (chronological order, Feb ... Jan) and vegetation type
# (called "Scrub" in the file).
vars <- read.csv(file.path(base_dir, "data", "variables.csv"),
                 stringsAsFactors = FALSE)   # 'vars', not 'var' (var is a base function)

# Split by vegetation type (12 months each)
dat_list <- list(Forest    = subset(vars, Site == "Forest"),
                 Scrubland = subset(vars, Site == "Scrub"))
stopifnot(all(sapply(dat_list, nrow) == 12))

# Responses: overall, social and solitary bees; richness and abundance
responses <- data.frame(
  Group    = c("Overall", "Overall", "Social", "Social", "Solitary", "Solitary"),
  Response = c("Richness", "Abundance", "Richness", "Abundance", "Richness", "Abundance"),
  Column   = c("Richness", "Abundance", "Social_richness", "Social_abundance",
               "Solitary_richness", "Solitary_abundance"),
  stringsAsFactors = FALSE
)

# Predictors: flowering plant richness (FLW), precipitation (PP) and temperature (TEMP)
predictors <- c("FLW", "PP", "TEMP")

# 2. Helper functions
# Pearson dispersion ratio (about 1 = no overdispersion)
pearson_ratio <- function(m) sum(residuals(m, type = "pearson")^2) / df.residual(m)

# Variance inflation factors (collinearity among predictors)
vif_manual <- function(dat, preds) {
  sapply(preds, function(p) {
    r2 <- summary(lm(reformulate(setdiff(preds, p), response = p), data = dat))$r.squared
    1 / (1 - r2)
  })
}

# Negative binomial fit that also keeps convergence warnings
fit_nb <- function(f, dat) {
  warn <- character(0)
  m <- withCallingHandlers(
    tryCatch(glm.nb(f, data = dat, control = glm.control(maxit = 100)),
             error = function(e) NULL),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(model = m, warn = paste(unique(warn), collapse = " | "))
}

# 3. Fit all models
# For each vegetation type, group and response: fit Poisson and negative binomial (NB),
# then choose the final error distribution with the rule above.
coef_tab <- list(); diag_tab <- list(); final_models <- list()

for (veg in names(dat_list)) {
  dat <- dat_list[[veg]]
  dat_s <- dat

  # Standardized copy of the data (coefficients per 1 SD, comparable among predictors)
  dat_s[predictors] <- scale(dat[predictors])      # standardized copy (beta per 1 SD)
  
  for (i in seq_len(nrow(responses))) {
    r  <- responses[i, ]
    f  <- reformulate(predictors, response = r$Column)
    id <- paste(veg, r$Group, r$Response)
    
    # Poisson model and overdispersion tests
    pois <- glm(f, family = poisson, data = dat)
    dt   <- AER::dispersiontest(pois)
    pr_p <- pearson_ratio(pois)

    # Negative binomial model
    nbf  <- fit_nb(f, dat)
    nb   <- nbf$model
    
    aic_p  <- AIC(pois)
    aic_nb <- if (is.null(nb)) NA else AIC(nb)
    
    # Choose the final model (rule in the header)
    final <- if (is.null(nb)) "Poisson"
    else if (dt$p.value < 0.05) "NB"
    else if (dt$p.value < 0.10 || pr_p > 1.5) { if (aic_nb < aic_p) "NB" else "Poisson" }
    else "Poisson"
    
    m_final <- if (final == "NB") nb else pois
    m_std   <- if (final == "NB") fit_nb(f, dat_s)$model else glm(f, family = poisson, data = dat_s)
    final_models[[id]] <- m_final
    
    # Residual autocorrelation (lag 1, Pearson residuals, chronological order)
    res <- residuals(m_final, type = "pearson")
    ac1 <- cor(res[-1], res[-length(res)])
    lb  <- Box.test(res, lag = 1, type = "Ljung-Box")$p.value
    
    # Collinearity
    vif <- vif_manual(dat, predictors)
    
    # Store diagnostics
    diag_tab[[id]] <- data.frame(
      Vegetation = veg, Group = r$Group, Response = r$Response,
      Disp_p = dt$p.value, Disp_est = unname(dt$estimate), Pearson_Pois = pr_p,
      AIC_Pois = aic_p, AIC_NB = aic_nb,
      Final = final, Pearson_final = pearson_ratio(m_final), AIC_final = AIC(m_final),
      theta = if (final == "NB") m_final$theta else NA,
      NB_warning = nbf$warn,
      AC1_resid = ac1, AC1_p = lb,
      VIF_FLW = vif["FLW"], VIF_PP = vif["PP"], VIF_TEMP = vif["TEMP"],
      row.names = NULL)
    
    # Store coefficients (raw and standardized) with 95% confidence intervals
    sm  <- summary(m_final)$coefficients[predictors, , drop = FALSE]
    sms <- summary(m_std)$coefficients[predictors, , drop = FALSE]
    ci  <- suppressMessages(confint.default(m_final))[predictors, , drop = FALSE]
    coef_tab[[id]] <- data.frame(
      Vegetation = veg, Group = r$Group, Response = r$Response, Model = final,
      Predictor = predictors,
      beta = sm[, 1], SE = sm[, 2], z = sm[, 3], p = sm[, 4],
      CI_low = ci[, 1], CI_high = ci[, 2],
      beta_std = sms[, 1],          # per 1 SD: comparable among predictors
      row.names = NULL)
  }
}

# Combine results; Holm-adjusted p-values over all 36 coefficient tests
coef_all <- do.call(rbind, coef_tab)
diag_all <- do.call(rbind, diag_tab)
coef_all$p_holm <- p.adjust(coef_all$p, method = "holm")   # 36 tests

# 4. Checks to read before writing
# Dispersion, final model, theta, warnings and autocorrelation
print(diag_all[, c("Vegetation","Group","Response","Disp_p","Pearson_Pois","Final",
                   "Pearson_final","theta","NB_warning","AC1_resid","AC1_p")], digits = 3)

# Collinearity: VIF > 5 (or > 3 if strict) means predictors cannot be separated with n = 12
print(diag_all[, c("Vegetation","Group","Response","VIF_FLW","VIF_PP","VIF_TEMP")], digits = 3)

# Correlations among predictors (one table per vegetation type)
lapply(dat_list, function(d) round(cor(d[predictors]), 2))

# Significant terms before and after Holm adjustment
print(subset(coef_all, p < 0.05)[, c("Vegetation","Group","Response","Predictor",
                                     "beta","SE","p","p_holm","beta_std")], digits = 3)

# 5. Export raw tables
fmt <- function(x, d = 3) round(x, d)
tab45 <- coef_all
tab45[c("beta","SE","z","beta_std","CI_low","CI_high")] <-
  lapply(tab45[c("beta","SE","z","beta_std","CI_low","CI_high")], fmt)
tab45$p      <- ifelse(coef_all$p < 0.001, "<0.001", fmt(coef_all$p))
tab45$p_holm <- ifelse(coef_all$p_holm < 0.001, "<0.001", fmt(coef_all$p_holm))
write.csv(tab45, file.path(out_dir, "Table4_5_GLM_coefficients.csv"), row.names = FALSE)

tabS5 <- diag_all
num <- sapply(tabS5, is.numeric); tabS5[num] <- lapply(tabS5[num], round, 3)
write.csv(tabS5, file.path(out_dir, "TableS5_GLM_diagnostics.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo_GLM.txt"))

# 6. Tables 4, 5 and S5 in manuscript format
# use_minus: TRUE = typographic minus sign (U+2212); FALSE = plain hyphen
use_minus <- TRUE   # TRUE: typographic minus (U+2212); FALSE: plain hyphen

# Formatting helpers
fmt_p  <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_b  <- function(x) ifelse(abs(x) < 0.01, sprintf("%.4f", x), sprintf("%.3f", x))  # small PP betas keep 4 decimals
fmt_3  <- function(x) sprintf("%.3f", x)
fmt_2  <- function(x) sprintf("%.2f", x)
fix_minus <- function(x) if (use_minus) sub("^-", "\u2212", x) else x

# Blank repeated labels (a label is blanked only if it and all labels to its left
# equal the previous row)
blank_repeats <- function(df, cols) {
  key <- rep("", nrow(df)); drop <- list()
  for (cl in cols) {
    key <- paste(key, df[[cl]], sep = "|")
    drop[[cl]] <- c(FALSE, key[-1] == key[-length(key)])
  }
  for (cl in cols) df[[cl]][drop[[cl]]] <- ""
  df
}

# Save a table as CSV (UTF-8 with BOM so Excel reads symbols correctly)
# and as tab-separated .txt (to paste into Word: Insert > Table > Convert Text to Table)
save_tab <- function(df, name) {
  f_csv <- file.path(out_dir, paste0(name, ".csv"))
  con <- file(f_csv, open = "wb")                       # write the UTF-8 BOM by hand
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  close(con)
  write.table(df, f_csv, sep = ",", qmethod = "double", row.names = FALSE,
              append = TRUE, fileEncoding = "UTF-8")
  write.table(df, file.path(out_dir, paste0(name, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
}

# Short labels
short_resp <- c(Richness = "Rich.", Abundance = "Abun.")

# Common ordering of the coefficient table
ct <- coef_all
ct$Vegetation <- factor(ct$Vegetation, levels = c("Forest", "Scrubland"))
ct$Group      <- factor(ct$Group,      levels = c("Overall", "Solitary", "Social"))
ct$Response   <- factor(ct$Response,   levels = c("Richness", "Abundance"))
ct$Predictor  <- factor(ct$Predictor,  levels = c("FLW", "PP", "TEMP"))
ct <- ct[order(ct$Vegetation, ct$Group, ct$Response, ct$Predictor), ]

# Build a table (one row per predictor)
make_tab <- function(d, with_group) {
  out <- data.frame(
    Vegetation = as.character(d$Vegetation),
    Group      = as.character(d$Group),
    Resp.      = unname(short_resp[as.character(d$Response)]),
    Model      = d$Model,
    Predictor  = as.character(d$Predictor),
    Beta       = fix_minus(fmt_b(d$beta)),
    SE         = fmt_b(d$SE),
    z          = fix_minus(fmt_3(d$z)),
    p          = fmt_p(d$p),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "Beta"] <- "\u03b2"
  cols <- c("Vegetation", if (with_group) "Group", "Resp.", "Model")
  if (!with_group) out$Group <- NULL
  out <- blank_repeats(out, cols)
  out
}

# Table 4: complete assemblage
tab4 <- make_tab(subset(ct, Group == "Overall"), with_group = FALSE)
print(tab4, row.names = FALSE)
save_tab(tab4, "Table4_GLM_overall")

# Table 5: solitary and social bees
tab5 <- make_tab(subset(ct, Group %in% c("Solitary", "Social")), with_group = TRUE)
print(tab5, row.names = FALSE)
save_tab(tab5, "Table5_GLM_sociality")

# Table S5: dispersion diagnostics
ds <- diag_all
ds$Vegetation <- factor(ds$Vegetation, levels = c("Forest", "Scrubland"))
ds$Group      <- factor(ds$Group,      levels = c("Overall", "Solitary", "Social"))
ds$Response   <- factor(ds$Response,   levels = c("Richness", "Abundance"))
ds <- ds[order(ds$Vegetation, ds$Group, ds$Response), ]

tabS5 <- data.frame(
  Vegetation = as.character(ds$Vegetation),
  Group      = as.character(ds$Group),
  Resp.      = unname(short_resp[as.character(ds$Response)]),
  Disp.p     = fmt_p(ds$Disp_p),
  Disp.est   = fmt_3(ds$Disp_est),
  PR.init    = fmt_3(ds$Pearson_Pois),
  AIC.Pois   = fmt_2(ds$AIC_Pois),
  AIC.NB     = fmt_2(ds$AIC_NB),
  Final      = ifelse(ds$Final == "NB", "NB", "Poisson"),
  PR.final   = fmt_3(ds$Pearson_final),
  theta      = ifelse(is.na(ds$theta), "\u2013", fmt_2(ds$theta)),
  Lag1.r     = fix_minus(fmt_3(ds$AC1_resid)),
  stringsAsFactors = FALSE
)
names(tabS5) <- c("Vegetation", "Group", "Resp.", "Disp. p", "Disp. est.", "PR init.",
                  "AIC Pois.", "AIC NB", "Final model", "PR final", "\u03b8", "Lag-1 r")
tabS5 <- blank_repeats(tabS5, c("Vegetation", "Group"))
print(tabS5, row.names = FALSE)
save_tab(tabS5, "TableS5_GLM_diagnostics_manuscript")
