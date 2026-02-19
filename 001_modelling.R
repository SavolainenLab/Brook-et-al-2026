################################################################################
# SCRIPTS FOR ALL THE DROUGHT MODELLING ANALYSIS IN BROOK ET AL. (2026)
#
#
# Article title: Genomic architecture and breeding trade-offs of coconut
#                drought tolerance
#
# Authors (R script): Theodore S. Brook, Tom Versluys (script 1)
#
# Authors (manuscript): Theodore S. Brook, Tom Versluys, Bi Tra Serges Doubi,
#          Jahcub Trew, Sam Acors, Kyle Macleod, Louis Bell-Roberts, 
#          Katarina Pipoini, Jean Louis Konan, Konan Engueran Djaha, Hala N'klo,
#          Benjamin J Roberts, CM (Tilly) Collins, Matteo Fumagalli,
#          Julien Godwin, Vincent Savolainen
################################################################################


################################################################################
# SCRIPT 1: ANNUAL CLIMATE PROCESSING
#
# Description:
#   Processes monthly climate data to derive drought indices (SPI, SPEI, RDI),
#   rolling climate metrics, and joins processed
#   climate data to filtered yield records and produces an annualised model
#   dataset for downstream mixed-model analyses.
#
# Input files required:
#   - yield_data.csv
#   - climate_data.csv
#
# Output:
#   - annualised_model_data.csv
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(ggthemes)
  library(SPEI)
  library(zoo)
})

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

yield_data   <- fread("yield_data.csv")
climate_data <- fread("climate_data.csv")

cat("Yield data loaded:", nrow(yield_data), "rows,",
    length(unique(yield_data$new_tree_id)), "unique trees\n")
cat("Unique crosses:", length(unique(yield_data$cross)), "\n\n")


# ==============================================================================
# 3. ANNUAL CLIMATE SUMMARIES
# ==============================================================================

climate_annual <- climate_data %>%
  group_by(year) %>%
  mutate(
    annual_temp           = mean(mean_monthly_temperature, na.rm = TRUE),
    total_annual_rainfall = sum(total_monthly_rainfall,   na.rm = TRUE)
  ) %>%
  data.table()


# ==============================================================================
# 4. POTENTIAL EVAPOTRANSPIRATION (PET)
# ==============================================================================

# PET estimated via Hargreaves and Thornthwaite methods (latitude = 3.23 N)
# Climatic water balance = Precipitation - PET (Hargreaves)

climate_annual <- climate_annual %>%
  mutate(
    PET_hargreaves         = hargreaves(min_monthly_temperature, max_monthly_temperature,
                                        lat = 3.23, na.rm = TRUE),
    PET_thornthwaite       = thornthwaite(mean_monthly_temperature, lat = 3.23, na.rm = TRUE),
    climatic_water_balance = total_monthly_rainfall - PET_hargreaves
  )


# ==============================================================================
# 5. DROUGHT INDICES
# ==============================================================================

# --- 5a. Standardised Precipitation Index (SPI) ---
# Fitted to monthly rainfall series at 6-, 12-, and 18-month time scales

spi_6  <- spi(climate_annual$total_monthly_rainfall, 6)
spi_12 <- spi(climate_annual$total_monthly_rainfall, 12)
spi_18 <- spi(climate_annual$total_monthly_rainfall, 18)

# --- 5b. Standardised Precipitation-Evapotranspiration Index (SPEI) ---
# Fitted to climatic water balance series at the same time scales

spei_6  <- spei(climate_annual$climatic_water_balance, 6)
spei_12 <- spei(climate_annual$climatic_water_balance, 12)
spei_18 <- spei(climate_annual$climatic_water_balance, 18)

# Append fitted values to data frame
climate_annual <- climate_annual %>%
  mutate(
    spi_6   = as.numeric(spi_6$fitted),
    spi_12  = as.numeric(spi_12$fitted),
    spi_18  = as.numeric(spi_18$fitted),
    spei_6  = as.numeric(spei_6$fitted),
    spei_12 = as.numeric(spei_12$fitted),
    spei_18 = as.numeric(spei_18$fitted)
  )


# ==============================================================================
# 6. ROLLING CLIMATE METRICS
# ==============================================================================

# Helper: rolling sum (returns NA where window evaluates to 0)
compute_rolling_totals <- function(data, n, var = "total_monthly_rainfall") {
  x      <- data[[var]]
  totals <- numeric(length(x))
  for (i in n:length(x)) totals[i] <- sum(x[(i - n + 1):i])
  ifelse(totals == 0, NA, totals)
}

# Helper: rolling mean
compute_rolling_means <- function(data, n, var = "mean_monthly_temperature") {
  x     <- data[[var]]
  means <- numeric(length(x))
  for (i in n:length(x)) means[i] <- mean(x[(i - n + 1):i])
  ifelse(means == 0, NA, means)
}

# Helper: rolling count of dry months (rainfall < 150 mm)
compute_rolling_dry_months <- function(data, n, threshold = 150,
                                       var = "total_monthly_rainfall") {
  x      <- data[[var]]
  counts <- numeric(length(x))
  for (i in n:length(x)) counts[i] <- sum(x[(i - n + 1):i] < threshold)
  counts
}

# Apply rolling metrics at 6-, 12-, and 18-month windows
for (n in c(6, 12, 18)) {
  climate_annual[[paste0("total_rainfall_",   n)]] <- compute_rolling_totals(climate_annual, n)
  climate_annual[[paste0("mean_temp_",        n)]] <- compute_rolling_means(climate_annual, n)
  climate_annual[[paste0("months_under_150_", n)]] <- compute_rolling_dry_months(climate_annual, n)
}


# ==============================================================================
# 7. RECONNAISSANCE DROUGHT INDEX (RDI)
# ==============================================================================

# Rolling PET sums (Thornthwaite) via zoo::rollapply
for (n in c(6, 12, 18)) {
  climate_annual[[paste0("pet_", n)]] <- zoo::rollapply(
    climate_annual$PET_thornthwaite,
    width = n, FUN = sum, align = "right", fill = NA
  )
}

# Standardised RDI = (alpha - mean_alpha) / sd_alpha, where alpha = P / PET
climate_annual <- climate_annual %>%
  mutate(
    ratio_6  = total_rainfall_6  / pet_6,
    rdi_6    = (ratio_6  - mean(ratio_6,  na.rm = TRUE)) / sd(ratio_6,  na.rm = TRUE),
    ratio_12 = total_rainfall_12 / pet_12,
    rdi_12   = (ratio_12 - mean(ratio_12, na.rm = TRUE)) / sd(ratio_12, na.rm = TRUE),
    ratio_18 = total_rainfall_18 / pet_18,
    rdi_18   = (ratio_18 - mean(ratio_18, na.rm = TRUE)) / sd(ratio_18, na.rm = TRUE)
  )


# ==============================================================================
# 8. DIAGNOSTIC PLOTS: DISTRIBUTIONS AND MISSINGNESS
# ==============================================================================

plot_distributions <- function(data, cols, title) {
  data %>%
    pivot_longer(cols = all_of(cols), names_to = "variable", values_to = "value") %>%
    ggplot(aes(value)) +
    geom_histogram(bins = 30, colour = "white") +
    facet_wrap(~variable, scales = "free") +
    theme_few() +
    labs(title = title, x = NULL, y = "Count")
}

print(plot_distributions(climate_annual,
                         c("spi_6", "spi_12", "spi_18"),
                         "SPI Distributions"))
print(plot_distributions(climate_annual,
                         c("spei_6", "spei_12", "spei_18"),
                         "SPEI Distributions"))
print(plot_distributions(climate_annual,
                         c("total_rainfall_6", "total_rainfall_12",
                           "total_rainfall_18"),
                         "Rolling Total Rainfall Distributions"))
print(plot_distributions(climate_annual,
                         c("rdi_6", "rdi_12", "rdi_18"),
                         "RDI Distributions"))

plot_missingness <- function(data, cols, title) {
  data %>%
    pivot_longer(cols = all_of(cols), names_to = "variable", values_to = "value") %>%
    mutate(status = factor(ifelse(is.na(value), "Missing", "Present"))) %>%
    group_by(variable, month, status) %>%
    summarise(n = n(), .groups = "drop") %>%
    ggplot(aes(x = status, y = n, fill = status)) +
    geom_bar(stat = "identity") +
    facet_grid(variable ~ month) +
    theme_few() +
    labs(title = title, x = NULL, y = "Count", fill = NULL)
}

print(plot_missingness(climate_annual,
                       c("spi_6", "spi_12", "spi_18"),
                       "SPI Missingness by Month"))
print(plot_missingness(climate_annual,
                       c("spei_6", "spei_12", "spei_18"),
                       "SPEI Missingness by Month"))


# ==============================================================================
# 9. JOIN YIELD AND CLIMATE DATA
# ==============================================================================

all_data <- inner_join(yield_data, climate_annual, by = c("year", "month"))

cat("After joining yield and climate data:", nrow(all_data), "rows\n")


# ==============================================================================
# 10. CONSTRUCT MODEL DATASET
# ==============================================================================

all_data <- all_data %>%
  mutate(
    new_tree_id   = as.factor(new_tree_id),
    num_fruit_log = as.numeric(log(num_fruit + 1)),
    tree_age      = as.numeric(year - year_planted),
    year_planted  = as.factor(year_planted),
    year          = as.factor(year),
    field         = as.factor(field),
    month         = as.factor(month)
  ) %>%
  group_by(new_tree_id) %>%
  mutate(
    n    = n(),
    time = seq_len(n())
  ) %>%
  ungroup()

model_vars <- c(
  "new_tree_id", "row", "tree", "mother", "father", "cross",
  "field", "tree_age", "year_planted", "num_fruit", "num_fruit_log",
  "n", "month", "time", "year",
  "spi_6",  "spi_12",  "spi_18",
  "rdi_6",  "rdi_12",  "rdi_18",
  "spei_6", "spei_12", "spei_18",
  "mean_temp_6",      "mean_temp_12",      "mean_temp_18",
  "total_rainfall_6", "total_rainfall_12", "total_rainfall_18"
)

all_data <- all_data %>%
  dplyr::select(all_of(model_vars)) %>%
  data.table()

cat("Unique trees in model dataset:", length(unique(all_data$new_tree_id)), "\n")


# ==============================================================================
# 11. ANNUALISE DATA
# ==============================================================================

# Average monthly climate metrics within each tree-year, then collapse to one
# row per tree-year with annual nut count summaries

annualised <- all_data %>%
  pivot_longer(
    cols = c(spi_6, spi_12, spi_18,
             rdi_6, rdi_12, rdi_18,
             spei_6, spei_12, spei_18,
             total_rainfall_6, total_rainfall_12, total_rainfall_18,
             mean_temp_6, mean_temp_12, mean_temp_18),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  group_by(new_tree_id, year, variable) %>%
  mutate(value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  data.table()

annualised <- annualised %>%
  group_by(new_tree_id, year) %>%
  mutate(
    yearly_nuts       = mean(num_fruit, na.rm = TRUE),
    yearly_nuts_log   = mean(num_fruit_log, na.rm = TRUE),
    yearly_nuts_log_1 = log(yearly_nuts + 1),
    yearly_nuts_count = as.integer(round(yearly_nuts, 0))
  ) %>%
  slice(1) %>%
  group_by(new_tree_id) %>%
  mutate(
    time = seq_len(n()),
    n    = n()
  ) %>%
  # Uncomment to restrict to trees with >= 12 years of data:
  # filter(n >= 12) %>%
  drop_na() %>%
  ungroup() %>%
  data.table()


# ==============================================================================
# 12. FINAL SUMMARY AND EXPORT
# ==============================================================================

cat("\n=== FINAL DATASET SUMMARY ===\n")
cat("Rows:           ", nrow(annualised), "\n")
cat("Unique trees:   ", length(unique(annualised$new_tree_id)), "\n")
cat("Unique crosses: ", length(unique(annualised$cross)), "\n")
cat("Years covered:  ",
    paste(range(as.numeric(as.character(annualised$year))), collapse = " to "), "\n")
cat("\nCross distribution:\n")
print(table(annualised$cross))

fwrite(annualised, "annualised_model_data.csv")
cat("\nSaved: annualised_model_data.csv\n")
cat("=== Script 1 complete ===\n")


################################################################################
# SCRIPT 2: ANNUAL DROUGHT METRIC SELECTION
#
# Description:
#   Screens candidate drought metrics by fitting linear and quadratic
#   negative-binomial mixed models (glmmTMB, nbinom2) to the annualised
#   yield dataset and comparing them by AIC. A coverage filter retains only
#   trees with sufficient observations across both dry and wet conditions
#   (base thresholds; not stringent). The best-fitting metric and
#   evidence for nonlinearity (ΔAIC = AIC_lin - AIC_quad) are reported.
#
# Filtering applied:  BASE coverage filter only
#   min_obs_per_tree : 6 finite observations (yield + metric_z + tree_age)
#   min_dry          : 3 observations with metric_z < 0
#   min_wet          : 3 observations with metric_z >= 0
#   (No overall-range threshold at this stage)
#
# Random effects (annual):
#   (1 | new_tree_id) + (1 | cross) + (1 | year_f)
#
# Outputs:
#   script_2/
#     metric_qc_annual.csv
#     metric_correlations_annual.csv
#     plot_correlations_annual.png
#     coverage_filter_log_annual.csv
#     metric_model_comparison_annual_lin_vs_quad.csv
#     plot_annual_metric_ranking_by_AICquad.png
#     plot_annual_deltaAIC_lin_vs_quad.png
#     annual_screening_summary.txt
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(glmmTMB)
})

# --- User settings ---

annual_file <- "annualised_model_data.csv"
yield_col   <- "yearly_nuts_count"
tree_id_col <- "new_tree_id"
year_col    <- "year"
cross_col   <- "cross"

candidate_metrics <- c(
  "rdi_6",  "rdi_12",  "rdi_18",
  "spei_6", "spei_12", "spei_18",
  "spi_6",  "spi_12",  "spi_18",
  "total_rainfall_6", "total_rainfall_12", "total_rainfall_18"
)

# Base coverage filter thresholds
min_obs_per_tree <- 6   # minimum finite complete cases per tree
min_dry          <- 3   # minimum observations with metric_z < 0
min_wet          <- 3   # minimum observations with metric_z >= 0

# Random effects switches
include_year_RE  <- TRUE
include_cross_RE <- TRUE

# Nonlinearity support threshold (ΔAIC = AIC_lin - AIC_quad)
deltaAIC_threshold <- 2

output_dir <- "script_2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

# --- 2a. Metric quality-control summary ---
metric_qc <- function(df, metrics) {
  tibble(metric = metrics) %>%
    rowwise() %>%
    mutate(
      n            = sum(!is.na(df[[metric]])),
      prop_missing = mean(is.na(df[[metric]])),
      mean         = mean(df[[metric]], na.rm = TRUE),
      sd           = sd(df[[metric]],   na.rm = TRUE),
      p01          = quantile(df[[metric]], 0.01, na.rm = TRUE),
      p99          = quantile(df[[metric]], 0.99, na.rm = TRUE)
    ) %>%
    ungroup()
}

# --- 2b. Safe z-score (returns NA vector if sd == 0 or < 2 finite values) ---
zscore_safe <- function(x) {
  x  <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) < 2) return(rep(NA_real_, length(x)))
  sx <- sd(x[ok])
  if (!is.finite(sx) || sx == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

# --- 2c. Base coverage filter (applied on standardised metric_z) ---
# Retains trees with >= min_obs finite complete cases, >= min_dry observations
# with metric_z < 0, and >= min_wet observations with metric_z >= 0
filter_trees_base <- function(df, metric_z_col, yield_col, age_col, tree_id_col,
                              min_obs   = min_obs_per_tree,
                              min_dry_n = min_dry,
                              min_wet_n = min_wet) {
  tree_stats <- df %>%
    group_by(.data[[tree_id_col]]) %>%
    summarise(
      n_obs = sum(is.finite(.data[[yield_col]]) &
                    is.finite(.data[[metric_z_col]]) &
                    is.finite(.data[[age_col]])),
      n_dry = sum(.data[[metric_z_col]] < 0,  na.rm = TRUE),
      n_wet = sum(.data[[metric_z_col]] >= 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_obs >= min_obs, n_dry >= min_dry_n, n_wet >= min_wet_n)
  
  df %>% filter(.data[[tree_id_col]] %in% tree_stats[[tree_id_col]])
}

# --- 2d. Build random-effects string ---
build_re_string <- function(tree_id_col, cross_col,
                            include_cross = TRUE, include_year = TRUE) {
  terms <- paste0("(1|", tree_id_col, ")")
  if (include_cross) terms <- c(terms, paste0("(1|", cross_col, ")"))
  if (include_year)  terms <- c(terms, "(1|year_f)")
  paste(terms, collapse = " + ")
}

# --- 2e. Fit linear and quadratic annual models; return AIC comparison row ---
fit_lin_quad_annual <- function(df, yield_col, re_str) {
  
  f_lin  <- as.formula(paste0(yield_col, " ~ metric_z + tree_age + ", re_str))
  f_quad <- as.formula(paste0(yield_col, " ~ metric_z + I(metric_z^2) + tree_age + ", re_str))
  
  m_lin  <- tryCatch(glmmTMB(f_lin,  family = nbinom2, data = df), error = function(e) NULL)
  m_quad <- tryCatch(glmmTMB(f_quad, family = nbinom2, data = df), error = function(e) NULL)
  if (is.null(m_lin) || is.null(m_quad)) return(NULL)
  
  s_lin  <- summary(m_lin)$coefficients$cond
  s_quad <- summary(m_quad)$coefficients$cond
  
  tibble(
    AIC_lin          = AIC(m_lin),
    AIC_quad         = AIC(m_quad),
    delta_AIC        = AIC(m_lin) - AIC(m_quad),
    beta_lin         = s_lin["metric_z",       "Estimate"],
    se_lin           = s_lin["metric_z",       "Std. Error"],
    p_lin            = s_lin["metric_z",       "Pr(>|z|)"],
    beta_lin_in_quad = s_quad["metric_z",      "Estimate"],
    se_lin_in_quad   = s_quad["metric_z",      "Std. Error"],
    p_lin_in_quad    = s_quad["metric_z",      "Pr(>|z|)"],
    beta_quad        = s_quad["I(metric_z^2)", "Estimate"],
    se_quad          = s_quad["I(metric_z^2)", "Std. Error"],
    p_quad           = s_quad["I(metric_z^2)", "Pr(>|z|)"],
    converged_lin    = isTRUE(m_lin$sdr$pdHess),
    converged_quad   = isTRUE(m_quad$sdr$pdHess)
  )
}


# ==============================================================================
# 3. LOAD ANNUAL DATA
# ==============================================================================

cat("=============================================================================\n")
cat("SCRIPT 2: ANNUAL DROUGHT METRIC SELECTION\n")
cat("  Base coverage filter | Linear vs Quadratic | glmmTMB nbinom2\n")
cat("=============================================================================\n\n")

stopifnot(file.exists(annual_file))
a <- fread(annual_file) %>%
  mutate(
    !!tree_id_col := as.factor(.data[[tree_id_col]]),
    !!cross_col   := as.factor(.data[[cross_col]]),
    year_f        = as.factor(.data[[year_col]])
  )

stopifnot("tree_age" %in% names(a))
stopifnot(yield_col  %in% names(a))

cat("Annual dataset:", nrow(a), "rows |",
    length(unique(a[[tree_id_col]])), "trees |",
    length(unique(a[[cross_col]])), "crosses\n\n")

candidate_metrics <- candidate_metrics[candidate_metrics %in% names(a)]
cat("Candidate metrics available:", length(candidate_metrics), "\n\n")
if (length(candidate_metrics) == 0) stop("No candidate metrics found in annual data.")


# ==============================================================================
# 4. QUALITY CONTROL AND CORRELATIONS
# ==============================================================================

qc <- metric_qc(a, candidate_metrics)
fwrite(qc, file.path(output_dir, "metric_qc_annual.csv"))

cor_mat <- cor(a[, ..candidate_metrics], use = "pairwise.complete.obs")
cor_df  <- as.data.frame(as.table(cor_mat)) %>%
  rename(metric1 = Var1, metric2 = Var2, r = Freq)
fwrite(cor_df, file.path(output_dir, "metric_correlations_annual.csv"))

p_cor <- ggplot(cor_df, aes(metric1, metric2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Annual drought metric correlations (raw values)",
       x = NULL, y = NULL, fill = "r")
ggsave(file.path(output_dir, "plot_correlations_annual.png"),
       p_cor, width = 10, height = 8, dpi = 150)


# ==============================================================================
# 5. MODEL FITTING: STANDARDISE -> BASE COVERAGE FILTER -> LINEAR VS QUADRATIC
# ==============================================================================

cat("Fitting annual linear vs quadratic models...\n")
cat("  Order: z-score metric -> base coverage filter -> glmmTMB (nbinom2)\n")
cat("  Base coverage filter thresholds: min_obs =", min_obs_per_tree,
    "| min_dry =", min_dry, "| min_wet =", min_wet, "\n\n")

re_str <- build_re_string(tree_id_col, cross_col,
                          include_cross = include_cross_RE,
                          include_year  = include_year_RE)

results_list   <- list()
coverage_filter_log <- list()

for (met in candidate_metrics) {
  
  cat("  Metric:", met)
  
  # Step 1: Standardise
  a_std <- a %>% mutate(metric_z = zscore_safe(.data[[met]]))
  if (all(!is.finite(a_std$metric_z))) {
    cat("  !! Skipped — zero variance or insufficient values\n")
    next
  }
  
  # Step 2: Base coverage filter on metric_z
  a_met   <- filter_trees_base(a_std, "metric_z", yield_col, "tree_age", tree_id_col)
  n_trees <- length(unique(a_met[[tree_id_col]]))
  cat("  [trees post-filter:", n_trees, "| rows:", nrow(a_met), "]\n")
  
  coverage_filter_log[[met]] <- tibble(metric = met, n_trees_filtered = n_trees, n_rows_filtered = nrow(a_met))
  
  if (n_trees < 10) {
    cat("    !! Skipped — fewer than 10 trees after coverage filter\n")
    next
  }
  
  # Step 3: Fit models
  result <- fit_lin_quad_annual(a_met, yield_col, re_str)
  
  if (!is.null(result)) {
    results_list[[met]] <- result %>%
      mutate(metric = met, scale = "annual", n = nrow(a_met), n_trees = n_trees)
  } else {
    cat("    !! Skipped — model failed to converge\n")
  }
}

coverage_log <- bind_rows(coverage_filter_log)
fwrite(coverage_log, file.path(output_dir, "coverage_filter_log_annual.csv"))

results <- bind_rows(results_list) %>%
  mutate(nonlinear_supported = delta_AIC >= deltaAIC_threshold) %>%
  arrange(AIC_quad)

fwrite(results, file.path(output_dir, "metric_model_comparison_annual_lin_vs_quad.csv"))
cat("\nResults saved:", file.path(output_dir, "metric_model_comparison_annual_lin_vs_quad.csv"), "\n\n")


# ==============================================================================
# 6. DIAGNOSTIC PLOTS
# ==============================================================================

if (nrow(results) > 0) {
  
  p_rank <- results %>%
    mutate(metric = reorder(metric, -AIC_quad)) %>%
    ggplot(aes(metric, AIC_quad)) +
    geom_point(size = 2) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Annual: Metric ranking by quadratic model AIC",
         x = NULL, y = "AIC (quadratic model)")
  ggsave(file.path(output_dir, "plot_annual_metric_ranking_by_AICquad.png"),
         p_rank, width = 9, height = 6, dpi = 150)
  
  p_delta <- results %>%
    mutate(metric = reorder(metric, delta_AIC)) %>%
    ggplot(aes(metric, delta_AIC, colour = nonlinear_supported)) +
    geom_point(size = 2) +
    geom_hline(yintercept = deltaAIC_threshold, linetype = "dashed", colour = "grey40") +
    coord_flip() +
    scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#D6604D"),
                        labels = c("FALSE" = "Not supported", "TRUE" = "Supported")) +
    theme_minimal() +
    labs(
      title    = "Annual: Support for nonlinearity (ΔAIC = AIC_lin − AIC_quad)",
      subtitle = paste0("Dashed line = ΔAIC threshold of ", deltaAIC_threshold),
      x        = NULL,
      y        = "ΔAIC (positive = quadratic preferred)",
      colour   = paste0("Nonlinear (ΔAIC >= ", deltaAIC_threshold, ")")
    )
  ggsave(file.path(output_dir, "plot_annual_deltaAIC_lin_vs_quad.png"),
         p_delta, width = 9, height = 6, dpi = 150)
}


# ==============================================================================
# 7. SCREENING SUMMARY
# ==============================================================================

if (nrow(results) > 0) {
  
  top    <- results %>% slice(1)
  runner <- results %>% slice(min(2, n()))
  
  summary_lines <- c(
    "ANNUAL DROUGHT METRIC SCREENING SUMMARY",
    "========================================",
    paste0("Date: ", Sys.Date()),
    "",
    paste0("Top metric by AIC_quad : ", top$metric,
           "  (AIC_quad = ", round(top$AIC_quad, 1),
           ", ΔAIC = ", round(top$delta_AIC, 1), ")"),
    paste0("Runner-up              : ", runner$metric,
           "  (AIC_quad = ", round(runner$AIC_quad, 1),
           ", ΔAIC = ", round(runner$delta_AIC, 1), ")"),
    "",
    "Random effects:",
    paste0("  ", re_str),
    "",
    "Filtering and standardisation order:",
    "  1. metric_z = z-score(metric) across full annual dataset",
    "  2. Base coverage filter applied on metric_z:",
    paste0("       min_obs = ", min_obs_per_tree,
           " | min_dry (metric_z < 0) = ", min_dry,
           " | min_wet (metric_z >= 0) = ", min_wet),
    "  3. Linear and quadratic models fitted using metric_z",
    "",
    "Decision checklist:",
    "  - Prefer metrics with lower AIC_quad and confirmed convergence.",
    "  - Metrics within ~2 AIC units are broadly equivalent; use biological",
    "    reasoning to choose among them.",
    paste0("  - Nonlinearity supported when ΔAIC >= ", deltaAIC_threshold, ".")
  )
  
  writeLines(summary_lines, con = file.path(output_dir, "annual_screening_summary.txt"))
  cat("Wrote:", file.path(output_dir, "annual_screening_summary.txt"), "\n")
  
} else {
  cat("No annual models were successfully fitted.\n")
}

cat("\n=== Script 2 complete ===\n")


################################################################################
# SCRIPT 3: REACTION NORMS (MIXED-MODEL ADJUSTMENT)
#
# Description:
#   Estimates per-tree drought response curves (reaction norms) using the
#   mixed-model adjustment strategy:
#
#     Stage 1 — Population adjustment
#       lmer(yield_raw ~ age_c + (1|new_tree_id) + (1|year))
#       yield_adj = residuals + grand_mean
#       This removes the shared effects of age and year from yield before
#       fitting per-tree curves, accounting for repeated measures on trees.
#
#     Stage 2 — Per-tree quadratic regression
#       yield_adj ~ drought_c + drought_c^2
#       where drought_c = drought_metric - tree_mean(drought_metric)
#
#   A linearised slope is computed over [x1_linearise, x2_linearise] as the
#   average derivative of the quadratic across that interval. This is used as
#   the primary GWAS phenotype.
#
#   Two filtered datasets are produced:
#     main_dataset        : base filter — minimum data quality thresholds
#     supporting_dataset_1: stringent filter — additionally requires coverage
#                           of both the dry and wet tails of the drought metric
#
# Inputs:
#   - annualised_model_data.csv         (from Script 1)
#   - individually_sequenced_crosses.csv
#
# Outputs:  script_3/
#   main_dataset/
#     reaction_norms_all_trees.csv
#     reaction_norms_main_dataset.csv
#     phenotypes_sequenced_main_dataset.csv
#     pheno_slope_q_lin_<x1>_<x2>.txt
#     pheno_intercept.txt  |  pheno_beta1.txt  |  pheno_beta2.txt
#     pipeline_summary.txt
#   supporting_dataset_1/
#     reaction_norms_all_trees_trimmed.csv
#     reaction_norms_supporting_dataset_1.csv
#     phenotypes_sequenced_supporting_dataset_1.csv
#     pheno_slope_q_lin_<x1>_<x2>_stringent.txt  (and intercept/beta variants)
#     pipeline_summary_stringent.txt
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lme4)
  library(ggplot2)
})

# --- User settings ---

annualised_file <- "annualised_model_data.csv"
sequenced_file  <- "individually_sequenced_crosses.csv"

drought_metric <- "rdi_18"
yield_col      <- "yearly_nuts_count"
age_col        <- "tree_age"
cross_col      <- "cross"

# Stage 1 mixed model: allow per-tree age slope in addition to intercept?
# FALSE = yield_raw ~ age_c + (1|new_tree_id) + (1|year)   [used in the paper]
# TRUE  = yield_raw ~ age_c + (age_c|new_tree_id) + (1|year) [heavier and thus not used]
mixed_use_random_age_slope       <- FALSE
mixed_fallback_to_random_intercept <- TRUE   # fall back if random-slope model fails

# Linearised slope interval [x1, x2]:
# slope = average derivative of the quadratic curve between these two drought values
x1_linearise <- -1
x2_linearise <-  0

# --- Base filter (main_dataset) ---
min_obs_per_tree  <- 6     # minimum finite (yield_adj, drought, age) observations
min_dry           <- 2     # minimum years with drought_metric < 0
min_wet           <- 2     # minimum years with drought_metric >= 0
min_drought_range <- 0.5   # minimum range of drought_metric across observed years

# --- Stringent filter (supporting_dataset_1) ---
# Applied on top of min_obs_per_tree; also trims extreme drought years
stringent_min_dry       <- 3      # minimum dry observations
stringent_min_wet       <- 3      # minimum wet observations
stringent_min_neg_range <- 0.5    # minimum range within dry years (drought_metric < 0)
stringent_min_pos_range <- 0.5    # minimum range within wet years (drought_metric >= 0)
stringent_rdi_max       <- 1000   # remove years with drought_metric > this before fitting

output_dir <- "script_3"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

stopifnot(file.exists(annualised_file), file.exists(sequenced_file))

annualised_raw <- fread(annualised_file) %>%
  mutate(
    new_tree_id = as.character(new_tree_id),
    year        = as.factor(year),
    yield_raw   = .data[[yield_col]],
    age_raw     = .data[[age_col]]
  )

sequenced <- fread(sequenced_file) %>%
  mutate(
    new_tree_id   = as.character(new_tree_id),
    sequencing_id = as.character(sequencing_id)
  )

stopifnot(all(c("new_tree_id", "year", drought_metric, yield_col, age_col) %in% names(annualised_raw)))
stopifnot(all(c("new_tree_id", "sequencing_id") %in% names(sequenced)))

grand_mean <- mean(annualised_raw$yield_raw, na.rm = TRUE)
cat("Grand mean yield:", round(grand_mean, 2), "\n\n")


# ==============================================================================
# 3. STAGE 1: MIXED-MODEL POPULATION ADJUSTMENT
# ==============================================================================

cat("=== Stage 1: Mixed-model adjustment ===\n")

mix_dat <- annualised_raw %>%
  filter(is.finite(yield_raw), is.finite(age_raw), !is.na(year), !is.na(new_tree_id)) %>%
  mutate(age_c = as.numeric(scale(age_raw, center = TRUE, scale = FALSE)))

mix_form <- if (mixed_use_random_age_slope) {
  yield_raw ~ age_c + (age_c | new_tree_id) + (1 | year)
} else {
  yield_raw ~ age_c + (1 | new_tree_id) + (1 | year)
}

cat("Stage 1 model formula:\n")
print(mix_form)

m_mixed <- tryCatch(
  lmer(mix_form, data = mix_dat, REML = TRUE),
  error = function(e) e
)

if (inherits(m_mixed, "error") && mixed_use_random_age_slope && mixed_fallback_to_random_intercept) {
  cat("WARNING: random-slope model failed — falling back to random-intercept model.\n")
  mix_form <- yield_raw ~ age_c + (1 | new_tree_id) + (1 | year)
  m_mixed  <- lmer(mix_form, data = mix_dat, REML = TRUE)
}

cat("Fixed age_c effect:", round(fixef(m_mixed)["age_c"], 4), "\n")
cat("isSingular:", isSingular(m_mixed, tol = 1e-4), "\n")

mix_dat <- mix_dat %>%
  mutate(yield_adj = resid(m_mixed) + grand_mean) %>%
  select(new_tree_id, year, yield_adj)

annualised_mixed <- annualised_raw %>%
  left_join(mix_dat, by = c("new_tree_id", "year"))

cat("Rows with yield_adj:", sum(is.finite(annualised_mixed$yield_adj)),
    "/", nrow(annualised_mixed), "\n\n")


# ==============================================================================
# 4. HELPER FUNCTIONS
# ==============================================================================

# --- 4a. Per-tree quadratic regression (no age term; age removed in Stage 1) ---
# Returns named vector of coefficients and derived quantities.
# drought_c = drought - tree_mean(drought) to improve numerical stability.
# Linearised slope = average derivative of quadratic over [x1, x2].
calc_quad <- function(y, x, x1, x2) {
  
  ok <- is.finite(y) & is.finite(x)
  y  <- y[ok]; x <- x[ok]
  
  na_out <- c(beta0 = NA_real_, beta1 = NA_real_, beta2 = NA_real_,
              slope_lin = NA_real_, slope_lin_se = NA_real_,
              x_mean = NA_real_, baseline_at0 = NA_real_)
  
  if (length(y) < 4)   return(na_out)
  if (sd(x) == 0) { na_out["x_mean"] <- mean(x); return(na_out) }
  
  x_mean <- mean(x)
  x_c    <- x - x_mean
  
  fit <- tryCatch(lm(y ~ x_c + I(x_c^2)), error = function(e) NULL)
  if (is.null(fit)) { na_out["x_mean"] <- x_mean; return(na_out) }
  
  cf <- coef(fit);  V <- vcov(fit)
  b0 <- unname(cf["(Intercept)"])
  b1 <- unname(cf["x_c"])
  b2 <- unname(cf["I(x_c^2)"])
  
  # Linearised slope over [x1, x2] (integral of derivative / interval width)
  x1c <- x1 - x_mean;  x2c <- x2 - x_mean
  A1  <- (x2c - x1c)         / (x2 - x1)
  A2  <- (x2c^2 - x1c^2)     / (x2 - x1)
  slope_lin    <- b1 * A1 + b2 * A2
  slope_lin_se <- sqrt(A1^2 * V["x_c", "x_c"] +
                         A2^2 * V["I(x_c^2)", "I(x_c^2)"] +
                         2 * A1 * A2 * V["x_c", "I(x_c^2)"])
  
  # Predicted yield at drought_metric = 0 (population mean drought)
  x0c <- 0 - x_mean
  baseline_at0 <- b0 + b1 * x0c + b2 * x0c^2
  
  c(beta0 = b0, beta1 = b1, beta2 = b2,
    slope_lin = slope_lin, slope_lin_se = slope_lin_se,
    x_mean = x_mean, baseline_at0 = baseline_at0)
}

# --- 4b. Fit per-tree reaction norms across all trees in a dataset ---
fit_trees <- function(dat) {
  dat %>%
    group_by(new_tree_id) %>%
    summarise(
      n_obs         = sum(is.finite(yield_adj) & is.finite(.data[[drought_metric]])),
      n_dry         = sum(.data[[drought_metric]] <  0, na.rm = TRUE),
      n_wet         = sum(.data[[drought_metric]] >= 0, na.rm = TRUE),
      drought_min   = suppressWarnings(min(.data[[drought_metric]], na.rm = TRUE)),
      drought_max   = suppressWarnings(max(.data[[drought_metric]], na.rm = TRUE)),
      drought_range = drought_max - drought_min,
      neg_rdi_min   = suppressWarnings(min(.data[[drought_metric]][.data[[drought_metric]] <  0], na.rm = TRUE)),
      neg_rdi_max   = suppressWarnings(max(.data[[drought_metric]][.data[[drought_metric]] <  0], na.rm = TRUE)),
      neg_rdi_range = ifelse(is.finite(neg_rdi_min) & is.finite(neg_rdi_max),
                             neg_rdi_max - neg_rdi_min, 0),
      pos_rdi_min   = suppressWarnings(min(.data[[drought_metric]][.data[[drought_metric]] >= 0], na.rm = TRUE)),
      pos_rdi_max   = suppressWarnings(max(.data[[drought_metric]][.data[[drought_metric]] >= 0], na.rm = TRUE)),
      pos_rdi_range = ifelse(is.finite(pos_rdi_min) & is.finite(pos_rdi_max),
                             pos_rdi_max - pos_rdi_min, 0),
      tmp           = list(calc_quad(
        y  = yield_adj, x = .data[[drought_metric]],
        x1 = x1_linearise, x2 = x2_linearise
      )),
      intercept      = tmp[[1]]["beta0"],
      beta1          = tmp[[1]]["beta1"],
      beta2          = tmp[[1]]["beta2"],
      slope_q_lin    = tmp[[1]]["slope_lin"],
      slope_q_lin_se = tmp[[1]]["slope_lin_se"],
      x_mean         = tmp[[1]]["x_mean"],
      baseline_at0   = tmp[[1]]["baseline_at0"],
      .groups = "drop"
    ) %>%
    select(-tmp) %>%
    mutate(slope_t = abs(slope_q_lin) / slope_q_lin_se)
}

# --- 4c. Write PLINK-format phenotype file ---
write_pheno <- function(df, trait_col, path) {
  df %>%
    transmute(
      FID   = sequencing_id,
      IID   = sequencing_id,
      PHENO = as.numeric(.data[[trait_col]])
    ) %>%
    filter(is.finite(PHENO)) %>%
    arrange(FID) %>%
    fwrite(path, sep = "\t", quote = FALSE, col.names = TRUE)
}

# --- 4d. Print and log reaction norm descriptive summary ---
summarise_rn <- function(tree_all, tree_filtered, annualised, label) {
  
  all_ok <- tree_all %>% filter(is.finite(slope_q_lin))
  
  years_per_tree <- annualised %>%
    semi_join(all_ok %>% select(new_tree_id), by = "new_tree_id") %>%
    group_by(new_tree_id) %>%
    summarise(n_years = n_distinct(year), .groups = "drop")
  
  filtered_ids    <- unique(tree_filtered$new_tree_id)
  raw_obs         <- annualised_raw %>%
    filter(as.character(new_tree_id) %in% filtered_ids, is.finite(yield_raw))
  per_tree_means  <- raw_obs %>%
    group_by(new_tree_id) %>%
    summarise(tree_mean = mean(yield_raw, na.rm = TRUE), .groups = "drop")
  
  cat("\n---", label, "---\n")
  cat("Trees with fitted slopes (all)    :", nrow(all_ok), "\n")
  cat("Mean years per tree               :", round(mean(years_per_tree$n_years), 2),
      "±", round(sd(years_per_tree$n_years), 2), "\n")
  cat("Trees passing filter              :", nrow(tree_filtered), "\n")
  cat("Slope (linearised): mean =",
      round(mean(tree_filtered$slope_q_lin, na.rm = TRUE), 6),
      "SD =", round(sd(tree_filtered$slope_q_lin, na.rm = TRUE), 6), "\n")
  cat("Proportion positive slopes        :",
      round(mean(tree_filtered$slope_q_lin > 0, na.rm = TRUE), 3), "\n")
  cat("Proportion negative quadratic     :",
      round(mean(tree_filtered$beta2 < 0, na.rm = TRUE), 3), "\n")
  cat("Baseline yield @ drought = 0: mean =",
      round(mean(tree_filtered$baseline_at0, na.rm = TRUE), 2),
      "SD =", round(sd(tree_filtered$baseline_at0, na.rm = TRUE), 2), "\n")
  cat("Raw nut count (filtered trees): mean =",
      round(mean(raw_obs$yield_raw), 2),
      "SD =", round(sd(raw_obs$yield_raw), 2),
      "(N tree-years =", nrow(raw_obs), ")\n")
  cat("Per-tree mean yield: mean =",
      round(mean(per_tree_means$tree_mean), 2),
      "SD =", round(sd(per_tree_means$tree_mean), 2), "\n")
}


# ==============================================================================
# 5. STAGE 2: PER-TREE REACTION NORMS (ALL TREES)
# ==============================================================================

cat("=== Stage 2: Fitting per-tree reaction norms ===\n\n")

tree_all <- fit_trees(annualised_mixed)

cat("Trees with finite slopes:", sum(is.finite(tree_all$slope_q_lin)), "\n")


# ==============================================================================
# 6. MAIN DATASET (BASE FILTER)
# ==============================================================================

cat("\n=== Main dataset (base filter) ===\n")
cat("Thresholds: min_obs =", min_obs_per_tree,
    "| min_dry =", min_dry,
    "| min_wet =", min_wet,
    "| min_range =", min_drought_range, "\n\n")

main_dir <- file.path(output_dir, "main_dataset")
dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

log_main <- file.path(main_dir, "pipeline_summary.txt")
sink(log_main, split = TRUE)

tree_main <- tree_all %>%
  filter(
    n_obs         >= min_obs_per_tree,
    n_dry         >= min_dry,
    n_wet         >= min_wet,
    drought_range >= min_drought_range,
    is.finite(slope_q_lin)
  )

pheno_main <- sequenced %>% inner_join(tree_main, by = "new_tree_id")

cat("Trees passing base filter    :", nrow(tree_main), "\n")
cat("Sequenced trees (main_dataset):", nrow(pheno_main), "\n")
cat("Linearised slope interval    :", x1_linearise, "to", x2_linearise, "\n")

summarise_rn(tree_all, tree_main, annualised_mixed, "MAIN DATASET")

# Save outputs
fwrite(tree_all,  file.path(main_dir, "reaction_norms_all_trees.csv"))
fwrite(tree_main, file.path(main_dir, "reaction_norms_main_dataset.csv"))
fwrite(pheno_main, file.path(main_dir, "phenotypes_sequenced_main_dataset.csv"))

write_pheno(pheno_main, "slope_q_lin",
            file.path(main_dir, sprintf("pheno_slope_q_lin_%s_%s.txt", x1_linearise, x2_linearise)))
write_pheno(pheno_main, "intercept", file.path(main_dir, "pheno_intercept.txt"))
write_pheno(pheno_main, "beta1",     file.path(main_dir, "pheno_beta1.txt"))
write_pheno(pheno_main, "beta2",     file.path(main_dir, "pheno_beta2.txt"))

sink()
cat("Outputs saved to:", main_dir, "\n")
cat("Summary log     :", log_main, "\n\n")


# ==============================================================================
# 7. SUPPORTING DATASET 1 (STRINGENT FILTER)
# ==============================================================================

cat("=== Supporting dataset 1 (stringent filter) ===\n")
cat("Additional thresholds: min_dry =", stringent_min_dry,
    "| min_wet =", stringent_min_wet,
    "| neg_range >= ", stringent_min_neg_range,
    "| pos_range >= ", stringent_min_pos_range,
    "| rdi_max =", stringent_rdi_max, "\n\n")

supp1_dir <- file.path(output_dir, "supporting_dataset_1")
dir.create(supp1_dir, showWarnings = FALSE, recursive = TRUE)

log_supp1 <- file.path(supp1_dir, "pipeline_summary_stringent.txt")
sink(log_supp1, split = TRUE)

# Trim extreme drought years before fitting
n_before <- nrow(annualised_mixed)
annualised_trimmed <- annualised_mixed %>%
  filter(.data[[drought_metric]] <= stringent_rdi_max)
n_after <- nrow(annualised_trimmed)

cat("Rows before RDI trim (>", stringent_rdi_max, "):", n_before, "\n")
cat("Rows after  RDI trim                           :", n_after, "\n")
cat("Rows removed                                   :", n_before - n_after, "\n\n")

tree_stringent_all <- fit_trees(annualised_trimmed)

tree_supp1 <- tree_stringent_all %>%
  filter(
    n_obs         >= min_obs_per_tree,
    n_dry         >= stringent_min_dry,
    n_wet         >= stringent_min_wet,
    neg_rdi_range >= stringent_min_neg_range,
    pos_rdi_range >= stringent_min_pos_range,
    is.finite(slope_q_lin)
  )

pheno_supp1 <- sequenced %>% inner_join(tree_supp1, by = "new_tree_id")

cat("Trees passing stringent filter         :", nrow(tree_supp1), "\n")
cat("Sequenced trees (supporting_dataset_1) :", nrow(pheno_supp1), "\n")

summarise_rn(tree_stringent_all, tree_supp1, annualised_trimmed, "SUPPORTING DATASET 1")

# Save outputs
fwrite(tree_stringent_all, file.path(supp1_dir, "reaction_norms_all_trees_trimmed.csv"))
fwrite(tree_supp1,         file.path(supp1_dir, "reaction_norms_supporting_dataset_1.csv"))
fwrite(pheno_supp1,        file.path(supp1_dir, "phenotypes_sequenced_supporting_dataset_1.csv"))

write_pheno(pheno_supp1, "slope_q_lin",
            file.path(supp1_dir, sprintf("pheno_slope_q_lin_%s_%s_stringent.txt", x1_linearise, x2_linearise)))
write_pheno(pheno_supp1, "intercept", file.path(supp1_dir, "pheno_intercept_stringent.txt"))
write_pheno(pheno_supp1, "beta1",     file.path(supp1_dir, "pheno_beta1_stringent.txt"))
write_pheno(pheno_supp1, "beta2",     file.path(supp1_dir, "pheno_beta2_stringent.txt"))

sink()
cat("Outputs saved to:", supp1_dir, "\n")
cat("Summary log     :", log_supp1, "\n\n")


# ==============================================================================
# 8. FINAL SUMMARY
# ==============================================================================

cat("\n=============================================================================\n")
cat("SCRIPT 3 COMPLETE\n")
cat("=============================================================================\n")
cat("Stage 1 model  : lmer(yield_raw ~ age_c + (1|new_tree_id) + (1|year))\n")
cat("Stage 2 model  : per-tree lm(yield_adj ~ drought_c + drought_c^2)\n")
cat("Drought metric :", drought_metric, "\n")
cat("Slope interval :", x1_linearise, "to", x2_linearise, "\n")
cat("\nDatasets written:\n")
cat("  main_dataset/          (base filter)      :", nrow(tree_main),  "trees |",
    nrow(pheno_main),  "sequenced\n")
cat("  supporting_dataset_1/  (stringent filter) :", nrow(tree_supp1), "trees |",
    nrow(pheno_supp1), "sequenced\n")
cat("=== Script 3 complete ===\n")


################################################################################
# SCRIPT 4: DROUGHT METRIC REACTION NORM CORRELATIONS
#
# Description:
#   For each candidate drought metric, fits per-tree quadratic reaction norms
#   using the mixed-adjusted phenotype (Stage 1, same approach as Script 3):
#     Stage 1: lmer(yield_raw ~ age_c + (1|new_tree_id) + (1|year))
#              yield_adj = residuals + grand_mean
#   Each metric is z-scored before fitting so that dry/wet classification
#   (metric_z < 0 vs >= 0) and filtering thresholds are comparable across
#   all metric types, including total rainfall.
#   Base filtering is applied (n_obs >= 6, n_dry >= 2, n_wet >= 2,
#   drought_range >= 0.5) and slopes are retained for sequenced trees only.
#
#   Then computes:
#     - Pairwise Pearson and Spearman correlation matrices of per-tree
#       linearised slopes across all metrics
#     - Mean pairwise correlations (off-diagonal)
#     - Heatmaps of both correlation matrices
#
# Inputs:
#   - annualised_model_data.csv         (from Script 1)
#   - individually_sequenced_crosses.csv
#
# Outputs:  script_4/
#   filter_log.csv
#   per_tree_slopes_all_metrics.csv
#   correlation_matrix_pearson.csv
#   correlation_matrix_spearman.csv
#   correlation_matrix_n.csv
#   correlation_summary.csv
#   heatmap_pearson.png
#   heatmap_spearman.png
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(tidyr)
})

# --- User settings ---

annualised_file <- "annualised_model_data.csv"
sequenced_file  <- "individually_sequenced_crosses.csv"

id_col    <- "new_tree_id"
year_col  <- "year"
yield_col <- "yearly_nuts_count"
age_col   <- "tree_age"

candidate_metrics <- c(
  "rdi_6",  "rdi_12",  "rdi_18",
  "spei_6", "spei_12", "spei_18",
  "spi_6",  "spi_12",  "spi_18",
  "total_rainfall_6",  "total_rainfall_12",
  "total_rainfall_18"
)

# Base filter settings
min_obs_per_tree  <- 6
min_dry           <- 2
min_wet           <- 2
min_drought_range <- 0.5

# Linearise quadratic slope between these values (on z-scored scale)
x1_linearise <- -1
x2_linearise <-  0

# Mixed-adjustment options
mixed_use_random_age_slope         <- FALSE
mixed_fallback_to_random_intercept <- TRUE

output_dir <- "script_4"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

stopifnot(file.exists(annualised_file), file.exists(sequenced_file))

annualised <- fread(annualised_file) %>%
  mutate(
    new_tree_id = as.character(.data[[id_col]]),
    year        = as.factor(.data[[year_col]]),
    yield_raw   = .data[[yield_col]],
    age_raw     = .data[[age_col]]
  )

sequenced <- fread(sequenced_file) %>%
  transmute(
    new_tree_id   = as.character(.data[[id_col]]),
    sequencing_id = as.character(sequencing_id)
  )

available_metrics <- candidate_metrics[candidate_metrics %in% names(annualised)]
missing_metrics   <- setdiff(candidate_metrics, available_metrics)

if (length(missing_metrics) > 0)
  warning("These metrics are not in the data and will be skipped: ",
          paste(missing_metrics, collapse = ", "))
if (length(available_metrics) == 0) stop("No candidate metrics found in the data.")

cat("Script 4: ", length(available_metrics), "metrics available\n")


# ==============================================================================
# 3. STAGE 1: MIXED-MODEL POPULATION ADJUSTMENT
# ==============================================================================

cat("=== Stage 1: Mixed-model adjustment ===\n")

grand_mean <- mean(annualised$yield_raw, na.rm = TRUE)

annualised <- annualised %>%
  mutate(age_c = age_raw - mean(age_raw, na.rm = TRUE))

m_mixed <- NULL
if (mixed_use_random_age_slope) {
  m_mixed <- tryCatch(
    lmer(yield_raw ~ age_c + (1|year) + (age_c|new_tree_id),
         data = annualised, REML = TRUE),
    error = function(e) NULL
  )
  if (!is.null(m_mixed) && isSingular(m_mixed, tol = 1e-4)) {
    if (mixed_fallback_to_random_intercept) m_mixed <- NULL
  }
}
if (is.null(m_mixed)) {
  m_mixed <- lmer(yield_raw ~ age_c + (1|year) + (1|new_tree_id),
                  data = annualised, REML = TRUE)
}

annualised <- annualised %>%
  mutate(yield_adj = resid(m_mixed) + grand_mean)

cat("Stage 1 complete. Grand mean =", round(grand_mean, 2), "\n\n")


# ==============================================================================
# 4. HELPER FUNCTIONS
# ==============================================================================

# --- 4a. Safe z-score ---
zscore_safe <- function(x) {
  x  <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) < 2) return(rep(NA_real_, length(x)))
  sx <- sd(x[ok])
  if (!is.finite(sx) || sx == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

# --- 4b. Per-tree quadratic reaction norm: linearised slope only ---
calc_quad_slope <- function(y, x, x1, x2) {
  ok <- is.finite(y) & is.finite(x)
  y  <- y[ok]; x <- x[ok]
  if (length(y) < 6 || sd(x) == 0) return(NA_real_)
  xm  <- mean(x); xc <- x - xm
  fit <- tryCatch(lm(y ~ xc + I(xc^2)), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  cf  <- coef(fit)
  b1  <- unname(cf["xc"]); b2 <- unname(cf["I(xc^2)"])
  x1c <- x1 - xm; x2c <- x2 - xm
  (b1 * (x2c - x1c) + b2 * (x2c^2 - x1c^2)) / (x2 - x1)
}


# ==============================================================================
# 5. FIT PER-TREE SLOPES FOR EACH METRIC
# ==============================================================================

cat("=== Fitting per-tree slopes for", length(available_metrics), "drought metrics ===\n")
cat("  (each metric is z-scored before dry/wet classification and filtering)\n\n")

slope_list  <- list()
filter_log  <- list()

for (metric in available_metrics) {
  
  cat("  ", metric, "...")
  
  # Z-score the metric across the full dataset
  annualised <- annualised %>%
    mutate(metric_z = zscore_safe(.data[[metric]]))
  
  if (sum(is.finite(annualised$metric_z)) == 0) {
    cat("  SKIPPED (z-scoring produced no finite values)\n")
    next
  }
  
  # Per-tree summaries and slope fitting on z-scored metric
  tree_slopes <- annualised %>%
    group_by(new_tree_id) %>%
    summarise(
      n_obs         = sum(is.finite(yield_adj) & is.finite(metric_z)),
      n_dry         = sum(metric_z <  0, na.rm = TRUE),
      n_wet         = sum(metric_z >= 0, na.rm = TRUE),
      drought_range = suppressWarnings(
        max(metric_z, na.rm = TRUE) - min(metric_z, na.rm = TRUE)),
      slope = calc_quad_slope(yield_adj, metric_z, x1_linearise, x2_linearise),
      .groups = "drop"
    ) %>%
    filter(
      n_obs >= min_obs_per_tree,
      n_dry >= min_dry,
      n_wet >= min_wet,
      drought_range >= min_drought_range,
      is.finite(slope)
    )
  
  # Restrict to sequenced trees
  tree_seq <- sequenced %>%
    inner_join(tree_slopes %>% select(new_tree_id, slope), by = "new_tree_id") %>%
    rename(!!metric := slope)
  
  slope_list[[metric]] <- tree_seq
  filter_log[[metric]] <- tibble(
    metric          = metric,
    trees_filtered  = nrow(tree_slopes),
    sequenced_trees = nrow(tree_seq)
  )
  
  cat("  trees (filtered):", nrow(tree_slopes),
      "| sequenced:", nrow(tree_seq), "\n")
}

fwrite(bind_rows(filter_log), file.path(output_dir, "filter_log.csv"))


# ==============================================================================
# 6. BUILD WIDE SLOPE TABLE (sequenced trees x metrics)
# ==============================================================================

valid_metrics <- names(slope_list)
if (length(valid_metrics) == 0) stop("No metrics produced any slopes.")

slope_wide <- slope_list[[1]]
for (i in seq_along(slope_list)[-1])
  slope_wide <- full_join(slope_wide, slope_list[[i]],
                          by = c("new_tree_id", "sequencing_id"))

metric_cols <- intersect(valid_metrics, names(slope_wide))

cat("\nCombined slope table:", nrow(slope_wide), "trees x",
    length(metric_cols), "metrics\n")
cat("Trees with slopes for ALL metrics:",
    sum(complete.cases(slope_wide[, ..metric_cols])), "\n\n")

fwrite(slope_wide, file.path(output_dir, "per_tree_slopes_all_metrics.csv"))


# ==============================================================================
# 7. PAIRWISE CORRELATION MATRICES (PEARSON + SPEARMAN)
# ==============================================================================

cat("=== Computing pairwise correlation matrices ===\n")

n_metrics    <- length(metric_cols)
cor_pearson  <- matrix(NA_real_,    n_metrics, n_metrics,
                       dimnames = list(metric_cols, metric_cols))
cor_spearman <- matrix(NA_real_,    n_metrics, n_metrics,
                       dimnames = list(metric_cols, metric_cols))
cor_n        <- matrix(NA_integer_, n_metrics, n_metrics,
                       dimnames = list(metric_cols, metric_cols))

for (i in seq_len(n_metrics)) {
  for (j in seq_len(n_metrics)) {
    x    <- slope_wide[[metric_cols[i]]]
    y    <- slope_wide[[metric_cols[j]]]
    ok   <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)
    cor_n[i, j] <- n_ok
    if (n_ok >= 3) {
      cor_pearson[i, j]  <- cor(x[ok], y[ok], method = "pearson")
      cor_spearman[i, j] <- cor(x[ok], y[ok], method = "spearman")
    }
  }
}

fwrite(as.data.frame(cor_pearson,  check.names = FALSE) %>%
         tibble::rownames_to_column("metric"),
       file.path(output_dir, "correlation_matrix_pearson.csv"))
fwrite(as.data.frame(cor_spearman, check.names = FALSE) %>%
         tibble::rownames_to_column("metric"),
       file.path(output_dir, "correlation_matrix_spearman.csv"))
fwrite(as.data.frame(cor_n, check.names = FALSE) %>%
         tibble::rownames_to_column("metric"),
       file.path(output_dir, "correlation_matrix_n.csv"))


# ==============================================================================
# 8. MEAN PAIRWISE CORRELATIONS (OFF-DIAGONAL)
# ==============================================================================

off_diag   <- function(mat) mat[lower.tri(mat)]

mean_r     <- mean(off_diag(cor_pearson),  na.rm = TRUE)
mean_rho   <- mean(off_diag(cor_spearman), na.rm = TRUE)
min_r      <- min(off_diag(cor_pearson),   na.rm = TRUE)
max_r      <- max(off_diag(cor_pearson),   na.rm = TRUE)
min_rho    <- min(off_diag(cor_spearman),  na.rm = TRUE)
max_rho    <- max(off_diag(cor_spearman),  na.rm = TRUE)

fwrite(
  tibble(
    statistic = c("mean_pearson_r",    "min_pearson_r",    "max_pearson_r",
                  "mean_spearman_rho", "min_spearman_rho", "max_spearman_rho",
                  "n_metrics",         "n_pairwise_comparisons"),
    value     = c(round(mean_r,   3), round(min_r,   3), round(max_r,   3),
                  round(mean_rho, 3), round(min_rho, 3), round(max_rho, 3),
                  n_metrics, choose(n_metrics, 2))
  ),
  file.path(output_dir, "correlation_summary.csv")
)

cat("\n=== CORRELATION SUMMARY ===\n")
cat("  Number of metrics  :", n_metrics, "\n")
cat("  Mean Pearson r     :", round(mean_r,   3),
    "(range:", round(min_r,   3), "to", round(max_r,   3), ")\n")
cat("  Mean Spearman rho  :", round(mean_rho, 3),
    "(range:", round(min_rho, 3), "to", round(max_rho, 3), ")\n\n")


# ==============================================================================
# 9. HEATMAPS
# ==============================================================================

cat("Generating correlation heatmaps...\n")

make_heatmap <- function(mat, stat_name, mean_val, metric_cols, outfile) {
  
  lab <- if (grepl("Pearson", stat_name, ignore.case = TRUE)) {
    bquote("Mean" ~ italic(r) ~ "=" ~ .(round(mean_val, 3)))
  } else {
    bquote("Mean" ~ rho ~ "=" ~ .(round(mean_val, 3)))
  }
  
  as.data.frame(mat, check.names = FALSE) %>%
    tibble::rownames_to_column("metric_x") %>%
    pivot_longer(-metric_x, names_to = "metric_y", values_to = "val") %>%
    mutate(
      metric_x = factor(metric_x, levels = metric_cols),
      metric_y = factor(metric_y, levels = rev(metric_cols))
    ) %>%
    ggplot(aes(x = metric_x, y = metric_y, fill = val)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", val)), size = 2.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = stat_name) +
    labs(x = NULL, y = NULL, subtitle = lab) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7)) +
    coord_fixed() -> p
  
  ggsave(outfile, p, width = 10, height = 9, dpi = 400)
}

make_heatmap(cor_pearson,  "Pearson r",    mean_r,
             metric_cols, file.path(output_dir, "heatmap_pearson.png"))
make_heatmap(cor_spearman, "Spearman rho", mean_rho,
             metric_cols, file.path(output_dir, "heatmap_spearman.png"))

cat("=== Script 4 complete. Outputs in:", output_dir, "===\n")


################################################################################
# SCRIPT 5: LINEAR VS QUADRATIC DROUGHT SLOPE COMPARISON
#
# Description:
#   Compares per-tree drought slope estimates derived from:
#     A) Linear model:    yield_adj ~ drought
#     B) Quadratic model: yield_adj ~ drought + drought^2
#                         with linearised slope computed between x1 and x2
#   using the mixed-adjusted Stage 1 phenotype and base filtering on rdi_18.
#
#   Stage 1: lmer(yield_raw ~ age_c + (1|new_tree_id) + (1|year))
#            yield_adj = residuals + grand_mean
#
#   Both slope estimates are computed for all trees passing the base filter.
#   Pearson and Spearman correlations between the two estimates are computed
#   and a scatter plot is produced using all filtered trees.
#
# Inputs:
#   - annualised_model_data.csv         (from Script 1)
#   - individually_sequenced_crosses.csv
#
# Outputs:  script_5/
#   per_tree_slopes_linear_and_quad.csv
#   per_tree_slopes_sequenced.csv
#   linear_vs_quad_correlation.csv
#   linear_vs_quad_slope_correlation.png
#   linear_vs_quad_slope_correlation.jpg
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(lme4)
})

# --- User settings ---

annualised_file <- "annualised_model_data.csv"
sequenced_file  <- "individually_sequenced_crosses.csv"

drought_metric <- "rdi_18"
yield_col      <- "yearly_nuts_count"
age_col        <- "tree_age"

# Base filter thresholds
min_obs_per_tree  <- 6
min_dry           <- 2
min_wet           <- 2
min_drought_range <- 0.5

# Linearisation interval for the quadratic slope
x1_linearise <- -1
x2_linearise <-  0

# Mixed-model Stage 1 options
mixed_use_random_age_slope         <- FALSE
mixed_fallback_to_random_intercept <- TRUE

output_dir <- "script_5"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

stopifnot(file.exists(annualised_file), file.exists(sequenced_file))

annualised_raw <- fread(annualised_file) %>%
  mutate(
    new_tree_id = as.character(new_tree_id),
    year        = as.factor(year),
    yield_raw   = .data[[yield_col]],
    age_raw     = .data[[age_col]]
  )

sequenced <- fread(sequenced_file) %>%
  mutate(
    new_tree_id   = as.character(new_tree_id),
    sequencing_id = as.character(sequencing_id)
  )

stopifnot(all(c("new_tree_id", "year", drought_metric, yield_col, age_col) %in%
                names(annualised_raw)))
stopifnot(all(c("new_tree_id", "sequencing_id") %in% names(sequenced)))

grand_mean <- mean(annualised_raw$yield_raw, na.rm = TRUE)
cat("Grand mean yield:", round(grand_mean, 2), "\n\n")


# ==============================================================================
# 3. STAGE 1: MIXED-MODEL POPULATION ADJUSTMENT
# ==============================================================================

cat("=== Stage 1: Mixed-model adjustment ===\n")

annualised_raw <- annualised_raw %>%
  mutate(age_c = age_raw - mean(age_raw, na.rm = TRUE))

m_mixed <- NULL
if (mixed_use_random_age_slope) {
  m_mixed <- tryCatch(
    lmer(yield_raw ~ age_c + (age_c | new_tree_id) + (1 | year),
         data = annualised_raw, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(m_mixed) || isSingular(m_mixed)) {
    cat("  Random-slope model singular/failed; falling back to random-intercept.\n")
    m_mixed <- NULL
  }
}
if (is.null(m_mixed)) {
  m_mixed <- lmer(yield_raw ~ age_c + (1 | new_tree_id) + (1 | year),
                  data = annualised_raw, REML = TRUE)
}

annualised <- annualised_raw %>%
  mutate(yield_adj = resid(m_mixed) + grand_mean)

cat("Stage 1 complete. Grand mean =", round(grand_mean, 2), "\n\n")


# ==============================================================================
# 4. HELPER FUNCTIONS
# ==============================================================================

# --- 4a. Per-tree linear slope ---
calc_linear <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  y  <- y[ok]; x <- x[ok]
  
  na_out <- c(slope_linear = NA_real_, slope_linear_se = NA_real_,
              n = NA_integer_)
  
  if (length(y) < 3 || sd(x) == 0) return(na_out)
  
  fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
  if (is.null(fit)) return(na_out)
  
  cf <- coef(fit); se <- sqrt(diag(vcov(fit)))
  c(slope_linear    = unname(cf["x"]),
    slope_linear_se = unname(se["x"]),
    n               = length(y))
}

# --- 4b. Per-tree quadratic linearised slope ---
calc_quad <- function(y, x, x1, x2) {
  ok <- is.finite(y) & is.finite(x)
  y  <- y[ok]; x <- x[ok]
  
  na_out <- c(beta1 = NA_real_, beta2 = NA_real_,
              slope_quad = NA_real_, slope_quad_se = NA_real_,
              baseline_at0 = NA_real_, x_mean = NA_real_)
  
  if (length(y) < 6 || sd(x) == 0) return(na_out)
  
  xm <- mean(x); xc <- x - xm
  fit <- tryCatch(lm(y ~ xc + I(xc^2)), error = function(e) NULL)
  if (is.null(fit)) return(na_out)
  
  cf <- coef(fit); V <- vcov(fit)
  b1 <- unname(cf["xc"]); b2 <- unname(cf["I(xc^2)"])
  
  x1c <- x1 - xm; x2c <- x2 - xm
  A1  <- (x2c - x1c) / (x2 - x1)
  A2  <- (x2c^2 - x1c^2) / (x2 - x1)
  sq  <- b1 * A1 + b2 * A2
  sq_se <- sqrt(A1^2 * V["xc", "xc"] +
                  A2^2 * V["I(xc^2)", "I(xc^2)"] +
                  2 * A1 * A2 * V["xc", "I(xc^2)"])
  
  x0c <- 0 - xm
  c(beta1        = b1,
    beta2        = b2,
    slope_quad    = sq,
    slope_quad_se = sq_se,
    baseline_at0  = unname(cf["(Intercept)"]) + b1 * x0c + b2 * x0c^2,
    x_mean        = xm)
}


# ==============================================================================
# 5. STAGE 2: FIT PER-TREE LINEAR AND QUADRATIC SLOPES
# ==============================================================================

cat("=== Stage 2: Fitting per-tree linear and quadratic slopes ===\n")

tree_slopes <- annualised %>%
  group_by(new_tree_id) %>%
  summarise(
    n_obs         = sum(is.finite(yield_adj) & is.finite(.data[[drought_metric]])),
    n_dry         = sum(.data[[drought_metric]] <  0, na.rm = TRUE),
    n_wet         = sum(.data[[drought_metric]] >= 0, na.rm = TRUE),
    drought_range = suppressWarnings(
      max(.data[[drought_metric]], na.rm = TRUE) -
        min(.data[[drought_metric]], na.rm = TRUE)),
    lin_tmp          = list(calc_linear(yield_adj, .data[[drought_metric]])),
    slope_linear     = lin_tmp[[1]]["slope_linear"],
    slope_linear_se  = lin_tmp[[1]]["slope_linear_se"],
    quad_tmp         = list(calc_quad(yield_adj, .data[[drought_metric]],
                                      x1_linearise, x2_linearise)),
    slope_quad       = quad_tmp[[1]]["slope_quad"],
    slope_quad_se    = quad_tmp[[1]]["slope_quad_se"],
    beta2            = quad_tmp[[1]]["beta2"],
    baseline_at0     = quad_tmp[[1]]["baseline_at0"],
    .groups = "drop"
  ) %>%
  select(-lin_tmp, -quad_tmp)


# ==============================================================================
# 6. BASE FILTERING
# ==============================================================================

tree_filtered <- tree_slopes %>%
  filter(
    n_obs         >= min_obs_per_tree,
    n_dry         >= min_dry,
    n_wet         >= min_wet,
    drought_range >= min_drought_range,
    is.finite(slope_linear),
    is.finite(slope_quad)
  )

cat("Trees passing base filter:", nrow(tree_filtered), "\n")
fwrite(tree_filtered, file.path(output_dir, "per_tree_slopes_linear_and_quad.csv"))

# Also save sequenced-tree subset for downstream GWAS use
pheno <- sequenced %>%
  inner_join(tree_filtered, by = "new_tree_id")

cat("Sequenced trees            :", nrow(pheno), "\n\n")
fwrite(pheno, file.path(output_dir, "per_tree_slopes_sequenced.csv"))


# ==============================================================================
# 7. CORRELATIONS (computed on ALL filtered trees)
# ==============================================================================

cat("=== Correlations: linear vs quadratic slope (all filtered trees) ===\n")

make_cor_label <- function(x, y) {
  ct    <- cor.test(x, y, method = "pearson")
  r_val <- round(ct$estimate, 3)
  p_raw <- ct$p.value
  if (p_raw == 0 || !is.finite(log10(p_raw))) {
    p_label <- bquote(italic(r) ~ "=" ~ .(r_val) * "," ~
                        italic(p) ~ "< 2.2" %*% 10^-16)
  } else {
    pe <- floor(log10(p_raw)); pm <- round(p_raw / 10^pe, 1)
    p_label <- bquote(italic(r) ~ "=" ~ .(r_val) * "," ~
                        italic(p) ~ "=" ~ .(pm) %*% 10^.(pe))
  }
  list(r = r_val, p_label = p_label)
}

cor_lq   <- make_cor_label(tree_filtered$slope_linear, tree_filtered$slope_quad)
rho_lq   <- cor.test(tree_filtered$slope_linear, tree_filtered$slope_quad, method = "spearman")

fwrite(
  tibble(
    comparison   = "linear vs quadratic linearised slope (rdi_18, all filtered trees)",
    pearson_r    = cor_lq$r,
    spearman_rho = round(rho_lq$estimate, 3),
    pearson_p    = round(cor.test(tree_filtered$slope_linear, tree_filtered$slope_quad)$p.value, 6),
    n            = nrow(tree_filtered)
  ),
  file.path(output_dir, "linear_vs_quad_correlation.csv")
)

cat("  Pearson r    :", cor_lq$r, "\n")
cat("  Spearman rho :", round(rho_lq$estimate, 3), "\n\n")


# ==============================================================================
# 8. PLOT (all filtered trees)
# ==============================================================================

p <- ggplot(tree_filtered, aes(x = slope_quad, y = slope_linear)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "red") +
  labs(
    x        = "Drought slope estimate (quadratic)",
    y        = "Drought slope estimate (linear)",
    subtitle = cor_lq$p_label
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(file.path(output_dir, "FigS10.png"),
       p, width = 6, height = 5, dpi = 400)

cat("=== Script 5 complete. Outputs in:", output_dir, "===\n")


################################################################################
# SCRIPT 6: TRADE-OFF BETWEEN DROUGHT SLOPE AND MEAN YIELD
#
# Description:
#   Tests and plots the trade-off between per-tree drought slope estimates
#   (derived from Script 3: mixed-model adjusted, quadratic reaction norm,
#   linearised between x1 and x2) and two per-tree mean yield estimates
#   computed directly from the raw annualised data:
#     Panel A) Mean observed yield across years where rdi_18 >= 0
#              (non-drought years only)
#     Panel B) Mean observed yield across all rdi_18 years
#
#   Drought slope estimates are read directly from the Script 3
#   main_dataset output (reaction_norms_main_dataset.csv) rather than
#   being re-derived here, so that the slope phenotype is identical to
#   that used in downstream GWAS.
#
#   Pearson correlations between the drought slope and each mean yield
#   estimate are reported, and a two-panel figure is produced.
#
# Inputs:
#   - annualised_model_data.csv                              (from Script 1)
#   - script_3/main_dataset/reaction_norms_main_dataset.csv  (from Script 3)
#
# Outputs:  script_6/
#   trade_off_slope_vs_yield_wet.csv
#   trade_off_slope_vs_yield_all.csv
#   trade_off_correlations.csv
#   FigS1.png
#   FigS1.jpg
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# --- User settings ---

annualised_file  <- "annualised_model_data.csv"
rn_file          <- "script_3/main_dataset/reaction_norms_main_dataset.csv"

drought_metric   <- "rdi_18"
yield_col        <- "yearly_nuts_count"

output_dir <- "script_6"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

stopifnot(file.exists(annualised_file))
stopifnot(file.exists(rn_file))

annualised_raw <- fread(annualised_file) %>%
  mutate(new_tree_id = as.character(new_tree_id))

rn <- fread(rn_file) %>%
  mutate(new_tree_id = as.character(new_tree_id))

stopifnot(all(c("new_tree_id", drought_metric, yield_col) %in% names(annualised_raw)))
stopifnot(all(c("new_tree_id", "slope_q_lin") %in% names(rn)))

# Retain only trees with a finite slope estimate
rn_valid <- rn %>% filter(is.finite(slope_q_lin))

cat("Trees with finite slope estimates:", nrow(rn_valid), "\n\n")


# ==============================================================================
# 3. HELPER: FORMAT P-VALUE IN SCIENTIFIC NOTATION FOR bquote
# ==============================================================================

fmt_p <- function(p) {
  if (!is.finite(log10(p)) || p == 0) return(list(man = 2.2, exp = -16, lt = TRUE))
  list(man = round(p / 10^floor(log10(p)), 1),
       exp = floor(log10(p)),
       lt  = FALSE)
}

make_subtitle <- function(cor_test) {
  r_val <- round(cor_test$estimate, 3)
  pf    <- fmt_p(cor_test$p.value)
  if (isTRUE(pf$lt)) {
    bquote(italic(r) ~ "=" ~ .(r_val) * "," ~
             italic(p) ~ "< 2.2" %*% 10^-16)
  } else {
    bquote(italic(r) ~ "=" ~ .(r_val) * "," ~
             italic(p) ~ "=" ~ .(pf$man) %*% 10^.(pf$exp))
  }
}


# ==============================================================================
# 4. PANEL A: NON-DROUGHT MEAN YIELD (RDI-18 >= 0) VS DROUGHT SLOPE
# ==============================================================================

cat("=== Panel A: Non-drought mean yield (RDI-18 >= 0) vs drought slope ===\n")

mean_wet <- annualised_raw %>%
  filter(new_tree_id %in% rn_valid$new_tree_id,
         .data[[drought_metric]] >= 0) %>%
  group_by(new_tree_id) %>%
  summarise(mean_yield = mean(.data[[yield_col]], na.rm = TRUE), .groups = "drop")

trade_wet <- mean_wet %>%
  inner_join(rn_valid %>% select(new_tree_id, slope_q_lin), by = "new_tree_id") %>%
  filter(is.finite(mean_yield), is.finite(slope_q_lin))

cat("Trees in Panel A:", nrow(trade_wet), "\n")

cor_A <- cor.test(trade_wet$mean_yield, trade_wet$slope_q_lin, method = "pearson")
rho_A <- cor.test(trade_wet$mean_yield, trade_wet$slope_q_lin, method = "spearman")

cat("  Pearson r    :", round(cor_A$estimate, 3), "\n")
cat("  Spearman rho :", round(rho_A$estimate, 3), "\n\n")

fwrite(trade_wet, file.path(output_dir, "trade_off_slope_vs_yield_wet.csv"))

pA <- ggplot(trade_wet, aes(x = slope_q_lin, y = mean_yield)) +
  geom_point(size = 2, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "red") +
  labs(
    x        = "Drought slope estimate",
    y        = bquote("Mean observed yield (RDI-18" ~ "\u2265" ~ "0)"),
    subtitle = make_subtitle(cor_A)
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())


# ==============================================================================
# 5. PANEL B: OVERALL MEAN YIELD (ALL RDI-18 YEARS) VS DROUGHT SLOPE
# ==============================================================================

cat("=== Panel B: Overall mean yield (all years) vs drought slope ===\n")

mean_all <- annualised_raw %>%
  filter(new_tree_id %in% rn_valid$new_tree_id) %>%
  group_by(new_tree_id) %>%
  summarise(mean_yield = mean(.data[[yield_col]], na.rm = TRUE), .groups = "drop")

trade_all <- mean_all %>%
  inner_join(rn_valid %>% select(new_tree_id, slope_q_lin), by = "new_tree_id") %>%
  filter(is.finite(mean_yield), is.finite(slope_q_lin))

cat("Trees in Panel B:", nrow(trade_all), "\n")

cor_B <- cor.test(trade_all$mean_yield, trade_all$slope_q_lin, method = "pearson")
rho_B <- cor.test(trade_all$mean_yield, trade_all$slope_q_lin, method = "spearman")

cat("  Pearson r    :", round(cor_B$estimate, 3), "\n")
cat("  Spearman rho :", round(rho_B$estimate, 3), "\n\n")

fwrite(trade_all, file.path(output_dir, "trade_off_slope_vs_yield_all.csv"))

pB <- ggplot(trade_all, aes(x = slope_q_lin, y = mean_yield)) +
  geom_point(size = 2, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "red") +
  labs(
    x        = "Drought slope estimate",
    y        = "Mean observed yield (all years)",
    subtitle = make_subtitle(cor_B)
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())


# ==============================================================================
# 6. SAVE CORRELATION SUMMARY
# ==============================================================================

fwrite(
  tibble(
    panel        = c("A: wet years (RDI-18 >= 0)", "B: all years"),
    n            = c(nrow(trade_wet), nrow(trade_all)),
    pearson_r    = c(round(cor_A$estimate, 3), round(cor_B$estimate, 3)),
    pearson_p    = c(signif(cor_A$p.value, 3), signif(cor_B$p.value, 3)),
    spearman_rho = c(round(rho_A$estimate, 3), round(rho_B$estimate, 3)),
    spearman_p   = c(signif(rho_A$p.value, 3), signif(rho_B$p.value, 3))
  ),
  file.path(output_dir, "trade_off_correlations.csv")
)


# ==============================================================================
# 7. COMBINE PANELS AND SAVE FIGURE
# ==============================================================================

combined_plot <- pA + pB +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(output_dir, "FigS1.png"),
       plot = combined_plot, width = 12, height = 5, dpi = 400)
ggsave(file.path(output_dir, "FigS1.jpg"),
       plot = combined_plot, width = 12, height = 5, dpi = 400)

cat("Figure saved: FigS1.png / FigS1.jpg\n")
cat("=== Script 6 complete. Outputs in:", output_dir, "===\n")

