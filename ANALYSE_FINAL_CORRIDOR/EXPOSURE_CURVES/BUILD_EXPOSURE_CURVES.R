#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

options(scipen = 999)

# ------------------------------------------------------------
# Multi-RP exposure curves for corridor municipalities
# ------------------------------------------------------------
# Goal:
# - prepare municipality-level exposure curves across RP10-RP500
# - derive simple curve metrics that can later inform protection logic
# - create a first typology:
#   early exposure / gradual increase / delayed jump
# ------------------------------------------------------------

paths <- list(
  exposure_csv = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv",
  quality_csv = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/exposure_quality_checks.csv",
  corridor_gpkg = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg",
  output_dir = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/EXPOSURE_CURVES/outputs"
)

dir.create(file.path(paths$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "gpkg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(paths$output_dir, "logs", "build_exposure_curves_log.txt")
if (file.exists(log_file)) unlink(log_file)

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

write_gpkg_layer <- function(x, dsn, layer) {
  if (file.exists(dsn)) unlink(dsn)
  write_sf(x, dsn = dsn, layer = layer, quiet = TRUE)
}

standardize_ags <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  str_pad(x, width = 8, side = "left", pad = "0")
}

first_rp_above <- function(values, rps, threshold) {
  idx <- which(values >= threshold)
  if (length(idx) == 0) return(NA_integer_)
  as.integer(rps[min(idx)])
}

trapz_scaled <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2) return(NA_real_)
  area <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  area / diff(range(x))
}

classify_curve <- function(onset_1pct, realized_by_rp100, late_growth_share, max_jump_interval) {
  if (!is.na(realized_by_rp100) && realized_by_rp100 >= 0.75) {
    return("early exposure")
  }

  if (
    !is.na(realized_by_rp100) &&
    realized_by_rp100 < 0.50 &&
    !is.na(late_growth_share) &&
    late_growth_share >= 0.40 &&
    max_jump_interval %in% c("RP100-RP200", "RP200-RP500")
  ) {
    return("delayed jump")
  }

  if (!is.na(onset_1pct) && onset_1pct <= 50 && !is.na(realized_by_rp100) && realized_by_rp100 >= 0.60) {
    return("early exposure")
  }

  "gradual increase"
}

rps <- c(10, 20, 50, 100, 200, 500)
rp_labels <- paste0("RP", rps)
interval_labels <- c("RP10-RP20", "RP20-RP50", "RP50-RP100", "RP100-RP200", "RP200-RP500")
share_cols <- paste0("flood_share_rp", rps)
area_cols <- paste0("flood_area_rp", rps, "_m2")
log_rps <- log10(rps)

log_message("Loading exposure table, quality checks, and corridor geometry.")

exposure <- read_csv(paths$exposure_csv, show_col_types = FALSE)
quality <- read_csv(paths$quality_csv, show_col_types = FALSE)
corridor <- st_read(paths$corridor_gpkg, quiet = TRUE)

corridor <- corridor %>%
  mutate(
    AGS_join = dplyr::coalesce(
      standardize_ags(AGS),
      standardize_ags(Gemeindeschlüssel_AGS)
    ),
    mun_name_join = dplyr::coalesce(
      as.character(mun_name),
      as.character(Gemeinde),
      as.character(GeografischerName_GEN)
    )
  )

if (!"municipality_area_m2" %in% names(corridor)) {
  stop("Corridor geometry has no municipality_area_m2 field.")
}

exposure <- exposure %>%
  mutate(
    AGS = standardize_ags(AGS),
    mun_name = as.character(mun_name)
  )

quality <- quality %>%
  mutate(AGS = standardize_ags(AGS))

if (nrow(exposure) == nrow(corridor)) {
  same_area <- isTRUE(all.equal(as.numeric(exposure$municipality_area_m2), as.numeric(corridor$municipality_area_m2)))

  if (same_area) {
    missing_idx <- which(is.na(exposure$AGS) | exposure$AGS == "")
    if (length(missing_idx) > 0) {
      exposure$AGS[missing_idx] <- corridor$AGS_join[missing_idx]
      exposure$mun_name[missing_idx] <- corridor$mun_name_join[missing_idx]
      log_message("Filled missing AGS/name in exposure table from corridor geometry by matched row order: ", length(missing_idx))
    }
  }
}

if (nrow(quality) == nrow(corridor) && "municipality_area_m2" %in% names(quality)) {
  same_area_quality <- isTRUE(all.equal(as.numeric(quality$municipality_area_m2), as.numeric(corridor$municipality_area_m2)))

  if (same_area_quality) {
    missing_idx_quality <- which(is.na(quality$AGS) | quality$AGS == "")
    if (length(missing_idx_quality) > 0) {
      quality$AGS[missing_idx_quality] <- corridor$AGS_join[missing_idx_quality]
      log_message("Filled missing AGS in quality table from corridor geometry by matched row order: ", length(missing_idx_quality))
    }
  }
}

if (any(is.na(exposure$AGS) | exposure$AGS == "")) {
  stop("Exposure table still contains missing AGS after key repair.")
}

exposure <- exposure %>%
  left_join(
    quality %>% select(AGS, monotonicity_all_ok, suspicious_any),
    by = "AGS"
  )

compute_curve_metrics <- function(ags, mun_name, municipality_area_m2, shares_obs, areas_obs, monotonicity_all_ok, suspicious_any) {
  shares_obs <- as.numeric(shares_obs)
  areas_obs <- as.numeric(areas_obs)
  shares_adj <- cummax(shares_obs)
  areas_adj <- shares_adj * municipality_area_m2
  deltas_adj <- diff(shares_adj)
  normalized <- if (shares_adj[length(shares_adj)] > 0) shares_adj / shares_adj[length(shares_adj)] else rep(0, length(shares_adj))

  slope_logrp <- if (stats::sd(shares_adj) == 0) 0 else unname(coef(lm(shares_adj ~ log_rps))[2])
  slope_norm_logrp <- if (stats::sd(normalized) == 0) 0 else unname(coef(lm(normalized ~ log_rps))[2])

  max_jump_idx <- if (all(deltas_adj == 0)) NA_integer_ else which.max(deltas_adj)
  max_jump_interval <- if (is.na(max_jump_idx)) NA_character_ else interval_labels[max_jump_idx]
  max_jump_abs <- if (is.na(max_jump_idx)) 0 else deltas_adj[max_jump_idx]
  max_jump_rel <- if (shares_adj[length(shares_adj)] > 0) max_jump_abs / shares_adj[length(shares_adj)] else 0

  onset_any <- first_rp_above(shares_adj, rps, threshold = 1e-10)
  onset_1pct <- first_rp_above(shares_adj, rps, threshold = 0.01)
  onset_5pct <- first_rp_above(shares_adj, rps, threshold = 0.05)

  realized_by_rp100 <- if (shares_adj[6] > 0) shares_adj[4] / shares_adj[6] else NA_real_
  realized_by_rp50 <- if (shares_adj[6] > 0) shares_adj[3] / shares_adj[6] else NA_real_
  late_growth_share <- if (shares_adj[6] > 0) (shares_adj[6] - shares_adj[4]) / shares_adj[6] else NA_real_
  normalized_auc <- trapz_scaled(log_rps, normalized)
  low_signal_flag <- shares_adj[6] < 0.02

  curve_type <- classify_curve(
    onset_1pct = onset_1pct,
    realized_by_rp100 = realized_by_rp100,
    late_growth_share = late_growth_share,
    max_jump_interval = max_jump_interval
  )

  metrics <- tibble(
    AGS = ags,
    mun_name = mun_name,
    municipality_area_m2 = municipality_area_m2,
    monotonicity_observed_ok = monotonicity_all_ok,
    suspicious_any = suspicious_any,
    monotone_adjustment_applied = !isTRUE(monotonicity_all_ok),
    flood_share_rp10_obs = shares_obs[1],
    flood_share_rp20_obs = shares_obs[2],
    flood_share_rp50_obs = shares_obs[3],
    flood_share_rp100_obs = shares_obs[4],
    flood_share_rp200_obs = shares_obs[5],
    flood_share_rp500_obs = shares_obs[6],
    flood_share_rp10_adj = shares_adj[1],
    flood_share_rp20_adj = shares_adj[2],
    flood_share_rp50_adj = shares_adj[3],
    flood_share_rp100_adj = shares_adj[4],
    flood_share_rp200_adj = shares_adj[5],
    flood_share_rp500_adj = shares_adj[6],
    onset_rp_any = onset_any,
    onset_rp_1pct = onset_1pct,
    onset_rp_5pct = onset_5pct,
    slope_logrp = slope_logrp,
    slope_norm_logrp = slope_norm_logrp,
    normalized_auc = normalized_auc,
    realized_by_rp50 = realized_by_rp50,
    realized_by_rp100 = realized_by_rp100,
    late_growth_share = late_growth_share,
    jump_rp10_rp20 = deltas_adj[1],
    jump_rp20_rp50 = deltas_adj[2],
    jump_rp50_rp100 = deltas_adj[3],
    jump_rp100_rp200 = deltas_adj[4],
    jump_rp200_rp500 = deltas_adj[5],
    max_jump_interval = max_jump_interval,
    max_jump_abs = max_jump_abs,
    max_jump_rel = max_jump_rel,
    low_signal_flag = low_signal_flag,
    exposure_curve_type = curve_type
  )

  long <- tibble(
    AGS = ags,
    mun_name = mun_name,
    return_period = rps,
    rp_label = rp_labels,
    flood_area_m2_obs = areas_obs,
    flood_share_obs = shares_obs,
    flood_area_m2_adj = areas_adj,
    flood_share_adj = shares_adj,
    normalized_share = normalized,
    exposure_curve_type = curve_type,
    monotone_adjustment_applied = !isTRUE(monotonicity_all_ok),
    low_signal_flag = low_signal_flag
  )

  list(metrics = metrics, long = long)
}

log_message("Computing curve metrics for each municipality.")

curve_results <- lapply(seq_len(nrow(exposure)), function(i) {
  compute_curve_metrics(
    ags = exposure$AGS[i],
    mun_name = exposure$mun_name[i],
    municipality_area_m2 = exposure$municipality_area_m2[i],
    shares_obs = exposure[i, share_cols, drop = TRUE],
    areas_obs = exposure[i, area_cols, drop = TRUE],
    monotonicity_all_ok = exposure$monotonicity_all_ok[i],
    suspicious_any = exposure$suspicious_any[i]
  )
})

curve_metrics <- bind_rows(lapply(curve_results, `[[`, "metrics"))
curve_long <- bind_rows(lapply(curve_results, `[[`, "long"))

type_summary <- curve_metrics %>%
  count(exposure_curve_type, sort = TRUE) %>%
  mutate(share_pct = round(100 * n / sum(n), 2))

onset_summary <- curve_metrics %>%
  count(onset_rp_1pct, sort = FALSE, name = "municipalities") %>%
  mutate(onset_rp_1pct = if_else(is.na(onset_rp_1pct), "not reached", as.character(onset_rp_1pct)))

metric_summary <- curve_metrics %>%
  group_by(exposure_curve_type) %>%
  summarise(
    municipalities = n(),
    mean_rp500_share = mean(flood_share_rp500_adj, na.rm = TRUE),
    mean_realized_by_rp100 = mean(realized_by_rp100, na.rm = TRUE),
    mean_late_growth_share = mean(late_growth_share, na.rm = TRUE),
    mean_normalized_auc = mean(normalized_auc, na.rm = TRUE),
    mean_max_jump_rel = mean(max_jump_rel, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(curve_metrics, file.path(paths$output_dir, "tables", "corridor_exposure_curve_metrics.csv"))
write_csv(curve_long, file.path(paths$output_dir, "tables", "corridor_exposure_curves_long.csv"))
write_csv(type_summary, file.path(paths$output_dir, "tables", "corridor_exposure_curve_type_summary.csv"))
write_csv(onset_summary, file.path(paths$output_dir, "tables", "corridor_onset_1pct_summary.csv"))
write_csv(metric_summary, file.path(paths$output_dir, "tables", "corridor_exposure_curve_metric_summary_by_type.csv"))

log_message("Joining metrics back to geometry.")

corridor_metrics_gpkg <- corridor %>%
  mutate(
    AGS = AGS_join,
    mun_name = mun_name_join
  ) %>%
  select(AGS, mun_name) %>%
  left_join(curve_metrics, by = c("AGS", "mun_name"))

write_gpkg_layer(
  corridor_metrics_gpkg,
  file.path(paths$output_dir, "gpkg", "corridor_exposure_curve_metrics.gpkg"),
  "corridor_exposure_curve_metrics"
)

log_message("Creating plots.")

curve_type_levels <- c("early exposure", "gradual increase", "delayed jump")
curve_long$exposure_curve_type <- factor(curve_long$exposure_curve_type, levels = curve_type_levels)
curve_metrics$exposure_curve_type <- factor(curve_metrics$exposure_curve_type, levels = curve_type_levels)

type_colors <- c(
  "early exposure" = "#b33b24",
  "gradual increase" = "#2f7e8a",
  "delayed jump" = "#d8a227"
)

mean_curve_plot <- curve_long %>%
  group_by(exposure_curve_type, return_period) %>%
  summarise(
    mean_normalized_share = mean(normalized_share, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = return_period, y = mean_normalized_share, color = exposure_curve_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = type_colors, drop = FALSE) +
  scale_x_log10(breaks = rps, labels = rp_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Mean normalized exposure curves by municipality type",
    x = "Return period",
    y = "Share of final RP500 exposure already realized",
    color = "Curve type"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(paths$output_dir, "plots", "plot_mean_normalized_curves_by_type.png"),
  mean_curve_plot,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

type_count_plot <- type_summary %>%
  mutate(exposure_curve_type = factor(exposure_curve_type, levels = curve_type_levels)) %>%
  ggplot(aes(x = exposure_curve_type, y = n, fill = exposure_curve_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = type_colors, drop = FALSE) +
  labs(
    title = "Municipality counts by exposure-curve type",
    x = NULL,
    y = "Municipalities"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave(
  file.path(paths$output_dir, "plots", "plot_exposure_curve_type_counts.png"),
  type_count_plot,
  width = 7,
  height = 5,
  dpi = 300
)

onset_plot <- curve_metrics %>%
  mutate(onset_rp_1pct_label = if_else(is.na(onset_rp_1pct), "not reached", paste0("RP", onset_rp_1pct))) %>%
  count(onset_rp_1pct_label) %>%
  ggplot(aes(x = onset_rp_1pct_label, y = n)) +
  geom_col(fill = "#2f7e8a", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  labs(
    title = "Onset of exposure at 1% municipal flood share",
    x = "First return period with share >= 1%",
    y = "Municipalities"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(paths$output_dir, "plots", "plot_onset_1pct_counts.png"),
  onset_plot,
  width = 7,
  height = 5,
  dpi = 300
)

type_map <- ggplot(corridor_metrics_gpkg) +
  geom_sf(aes(fill = exposure_curve_type), color = NA) +
  scale_fill_manual(values = type_colors, drop = FALSE) +
  labs(
    title = "Municipality exposure-curve types in the HQ500 corridor",
    fill = "Curve type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  file.path(paths$output_dir, "plots", "map_exposure_curve_types.png"),
  type_map,
  width = 7.5,
  height = 8,
  dpi = 300
)

adjusted_n <- sum(curve_metrics$monotone_adjustment_applied, na.rm = TRUE)

overview_lines <- c(
  "# Exposure Curves Overview",
  "",
  sprintf("Created: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Purpose",
  "",
  "This output prepares municipality-level exposure curves across return periods so that later protection-related interpretations can build on an already structured set of curve metrics.",
  "",
  "## Input",
  "",
  "- Exposure table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv`",
  "- Quality checks: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/exposure_quality_checks.csv`",
  "- Corridor geometry: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg`",
  "",
  "## Metric definitions",
  "",
  "- `onset_rp_any`: first return period with any positive flood share",
  "- `onset_rp_1pct`: first return period with flood share >= 1% of municipality area",
  "- `onset_rp_5pct`: first return period with flood share >= 5% of municipality area",
  "- `slope_logrp`: slope of monotone-adjusted flood share over log10(return period)",
  "- `normalized_auc`: area under the normalized exposure curve, scaled to `[0,1]`; higher values indicate earlier realization of final exposure",
  "- `realized_by_rp100`: share of final RP500 exposure already realized by RP100",
  "- `late_growth_share`: share of final RP500 exposure added after RP100",
  "- `max_jump_interval`: interval with the largest increase in flood share",
  "- `max_jump_rel`: largest interval increase relative to final RP500 exposure",
  "",
  "## Treatment of non-monotonic cases",
  "",
  sprintf("- Observed non-monotonic municipalities repaired with a simple `cummax` monotone adjustment: `%s`", adjusted_n),
  "- The adjusted curve is used for all derived curve metrics.",
  "",
  "## Typology rules",
  "",
  "- `early exposure`: at least 75% of final RP500 exposure already realized by RP100, or onset <= RP50 together with substantial early realization",
  "- `delayed jump`: less than 50% of final RP500 exposure realized by RP100, with at least 40% of final exposure added after RP100 and the largest jump in RP100-RP200 or RP200-RP500",
  "- `gradual increase`: all remaining municipalities",
  "",
  "## Output files",
  "",
  "- Metrics table: `outputs/tables/corridor_exposure_curve_metrics.csv`",
  "- Long curve table: `outputs/tables/corridor_exposure_curves_long.csv`",
  "- Type summary: `outputs/tables/corridor_exposure_curve_type_summary.csv`",
  "- Onset summary: `outputs/tables/corridor_onset_1pct_summary.csv`",
  "- Geometry output: `outputs/gpkg/corridor_exposure_curve_metrics.gpkg`",
  "- Mean normalized curves plot: `outputs/plots/plot_mean_normalized_curves_by_type.png`",
  "- Type-count plot: `outputs/plots/plot_exposure_curve_type_counts.png`",
  "- Type map: `outputs/plots/map_exposure_curve_types.png`"
)

writeLines(
  overview_lines,
  "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/EXPOSURE_CURVES/EXPOSURE_CURVES_OVERVIEW.md"
)

log_message("Done.")
