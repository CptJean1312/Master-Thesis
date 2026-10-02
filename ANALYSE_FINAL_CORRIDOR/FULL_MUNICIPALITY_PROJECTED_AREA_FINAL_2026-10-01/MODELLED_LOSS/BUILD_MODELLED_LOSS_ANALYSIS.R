#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(tidyr)
})

options(scipen = 999)

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"
external_root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS"
module_root <- file.path(
  root,
  "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/MODELLED_LOSS"
)
table_dir <- file.path(module_root, "outputs/tables")
gpkg_dir <- file.path(module_root, "outputs/gpkg")
log_dir <- file.path(module_root, "outputs/logs")

for (path in c(table_dir, gpkg_dir, log_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

paths <- list(
  analysis_csv = file.path(
    root,
    "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/outputs/tables/corridor_analysis_rp100.csv"
  ),
  analysis_gpkg = file.path(
    root,
    "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/outputs/gpkg/corridor_wide_pca_rp100_analysis.gpkg"
  ),
  landuse_csv = file.path(
    root,
    "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/LANDUSE/outputs/tables/corridor_landuse_exposure_wide.csv"
  ),
  curve_csv = file.path(
    root,
    "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/EXPOSURE_CURVES/outputs/tables/corridor_exposure_curve_metrics.csv"
  ),
  loss_csv = file.path(
    external_root,
    "PROTECTION/elbe_protection_level_mun.csv"
  )
)

missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) {
  stop("Missing required inputs:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

safe_cor <- function(x, y, method = "pearson") {
  suppressWarnings(cor(x, y, use = "complete.obs", method = method))
}

analysis <- read_csv(
  paths$analysis_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

landuse <- read_csv(
  paths$landuse_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
) %>%
  select(
    AGS,
    rp100_artificial_flooded_area_m2,
    rp100_artificial_flood_share_of_group
  )

curves <- read_csv(
  paths$curve_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
) %>%
  select(
    AGS,
    exposure_curve_type,
    realized_by_rp100,
    late_growth_share,
    normalized_auc
  )

loss_raw <- read_csv(
  paths$loss_csv,
  col_types = cols(
    ags = col_character(),
    municipality_name = col_character(),
    n_nonzero_years = col_double(),
    annual_loss_probability = col_double(),
    protection_return_period = col_double()
  ),
  show_col_types = FALSE
) %>%
  rename(
    AGS = ags,
    loss_dataset_municipality_name = municipality_name,
    finite_loss_occurrence_return_period = protection_return_period
  )

if (any(loss_raw$n_nonzero_years < 1) || anyNA(loss_raw$n_nonzero_years)) {
  stop("The modelled-loss extract contains an invalid positive-loss row.", call. = FALSE)
}

formula_check <- loss_raw %>%
  summarise(
    rows = n(),
    all_probability_values_match = all(
      abs(annual_loss_probability - n_nonzero_years / 5000) < 1e-12
    ),
    all_return_period_values_match = all(
      abs(finite_loss_occurrence_return_period - 5000 / n_nonzero_years) < 1e-9
    ),
    minimum_nonzero_years = min(n_nonzero_years),
    maximum_nonzero_years = max(n_nonzero_years)
  )

joined <- analysis %>%
  left_join(landuse, by = "AGS") %>%
  left_join(curves, by = "AGS") %>%
  left_join(loss_raw, by = "AGS") %>%
  mutate(
    positive_modelled_loss = !is.na(n_nonzero_years),
    modelled_loss_status = if_else(
      positive_modelled_loss,
      "At least one modelled loss year",
      "No modelled loss year in the 5,000-year catalogue"
    ),
    catalogue_nonzero_loss_years = coalesce(n_nonzero_years, 0),
    catalogue_annual_loss_probability = coalesce(annual_loss_probability, 0),
    finite_loss_occurrence_return_period = if_else(
      positive_modelled_loss,
      finite_loss_occurrence_return_period,
      NA_real_
    )
  )

if (nrow(joined) != 835 || sum(joined$positive_modelled_loss) != 280) {
  stop("Unexpected corridor or positive-loss sample size.", call. = FALSE)
}

write_csv(joined, file.path(table_dir, "corridor_modelled_loss_analysis.csv"))
write_csv(
  joined %>% filter(positive_modelled_loss),
  file.path(table_dir, "corridor_positive_modelled_loss_municipalities.csv")
)
write_csv(
  joined %>% filter(!positive_modelled_loss),
  file.path(table_dir, "corridor_no_event_municipalities.csv")
)
write_csv(formula_check, file.path(table_dir, "modelled_loss_formula_check.csv"))

overall_summary <- joined %>%
  summarise(
    corridor_municipalities = n(),
    positive_loss_municipalities = sum(positive_modelled_loss),
    no_event_municipalities = sum(!positive_modelled_loss),
    positive_loss_share = mean(positive_modelled_loss),
    minimum_positive_annual_probability = min(
      catalogue_annual_loss_probability[positive_modelled_loss]
    ),
    maximum_positive_annual_probability = max(
      catalogue_annual_loss_probability[positive_modelled_loss]
    ),
    minimum_finite_return_period = min(
      finite_loss_occurrence_return_period,
      na.rm = TRUE
    ),
    maximum_finite_return_period = max(
      finite_loss_occurrence_return_period,
      na.rm = TRUE
    )
  )
write_csv(overall_summary, file.path(table_dir, "modelled_loss_overall_summary.csv"))

correlation_targets <- tibble::tribble(
  ~x, ~y,
  "vuln_index_main_z", "positive_modelled_loss",
  "vuln_index_main_z", "catalogue_annual_loss_probability",
  "vuln_index_main_z", "finite_loss_occurrence_return_period",
  "flood_share_rp100", "positive_modelled_loss",
  "flood_share_rp100", "catalogue_annual_loss_probability",
  "flood_share_rp100", "finite_loss_occurrence_return_period",
  "rp100_artificial_flood_share_of_group", "positive_modelled_loss",
  "rp100_artificial_flood_share_of_group", "catalogue_annual_loss_probability",
  "rp100_artificial_flood_share_of_group", "finite_loss_occurrence_return_period"
)

correlations <- bind_rows(lapply(seq_len(nrow(correlation_targets)), function(i) {
  x_name <- correlation_targets$x[i]
  y_name <- correlation_targets$y[i]
  data_pair <- joined %>%
    transmute(
      x = as.numeric(.data[[x_name]]),
      y = as.numeric(.data[[y_name]]),
      positive_modelled_loss = positive_modelled_loss
    )
  if (y_name == "finite_loss_occurrence_return_period") {
    data_pair <- data_pair %>% filter(positive_modelled_loss)
  }
  tibble(
    x = x_name,
    y = y_name,
    complete_cases = sum(complete.cases(data_pair)),
    pearson = safe_cor(data_pair$x, data_pair$y, "pearson"),
    spearman = safe_cor(data_pair$x, data_pair$y, "spearman")
  )
}))
write_csv(correlations, file.path(table_dir, "modelled_loss_correlations.csv"))

quintile_data <- joined %>%
  filter(!is.na(vuln_index_main_z)) %>%
  mutate(vulnerability_quintile = ntile(vuln_index_main_z, 5))

quintile_summary <- quintile_data %>%
  group_by(vulnerability_quintile) %>%
  summarise(
    municipalities = n(),
    positive_loss_municipalities = sum(positive_modelled_loss),
    positive_loss_share = mean(positive_modelled_loss),
    mean_annual_loss_probability = mean(catalogue_annual_loss_probability),
    median_annual_loss_probability = median(catalogue_annual_loss_probability),
    median_finite_return_period_positive_only = median(
      finite_loss_occurrence_return_period[positive_modelled_loss],
      na.rm = TRUE
    ),
    .groups = "drop"
  )
write_csv(quintile_summary, file.path(table_dir, "modelled_loss_by_vulnerability_quintile.csv"))

data_dictionary <- tibble::tribble(
  ~variable, ~definition,
  "positive_modelled_loss", "TRUE when at least one non-zero modelled loss year occurs in the 5,000-year catalogue.",
  "catalogue_nonzero_loss_years", "Number of non-zero modelled loss years; zero for corridor municipalities absent from the positive-loss extract.",
  "catalogue_annual_loss_probability", "Catalogue frequency n/5,000; zero for no-event corridor municipalities.",
  "finite_loss_occurrence_return_period", "5,000/n for positive-loss municipalities only; undefined for no-event cases.",
  "modelled_loss_status", "Two-category catalogue outcome. It is not an observation of an engineering protection standard or zero real-world flood risk."
)
write_csv(data_dictionary, file.path(table_dir, "modelled_loss_data_dictionary.csv"))

geometry <- st_read(paths$analysis_gpkg, quiet = TRUE)
geometry_output <- geometry %>%
  select(AGS) %>%
  left_join(joined, by = "AGS")

st_write(
  geometry_output,
  file.path(gpkg_dir, "corridor_modelled_loss_analysis.gpkg"),
  layer = "corridor_modelled_loss_analysis",
  delete_dsn = TRUE,
  quiet = TRUE
)

log_lines <- c(
  paste("Raw positive-loss rows:", nrow(loss_raw)),
  paste("Corridor positive-loss municipalities:", sum(joined$positive_modelled_loss)),
  paste("Corridor no-event municipalities:", sum(!joined$positive_modelled_loss)),
  paste("Probability formula verified:", formula_check$all_probability_values_match),
  paste("Return-period formula verified:", formula_check$all_return_period_values_match)
)
writeLines(log_lines, file.path(log_dir, "build_modelled_loss_analysis.log"))
cat(paste(log_lines, collapse = "\n"), "\n")
