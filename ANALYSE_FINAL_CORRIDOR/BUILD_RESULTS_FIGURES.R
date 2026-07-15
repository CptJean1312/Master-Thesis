#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

options(scipen = 999)

# ============================================================
# Results figures for the final corridor thesis workflow
# ------------------------------------------------------------
# This script builds the result-facing maps and plots that connect
# exposure, the final vulnerability index, exposure curves, and the
# protection/loss portfolio.
# ============================================================

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"

paths <- list(
  analysis_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/outputs/gpkg/corridor_wide_pca_rp100_analysis.gpkg"),
  analysis_layer = "corridor_wide_pca_rp100",
  protection_csv = file.path(root, "ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/corridor_protection_level_all_corridor_municipalities.csv"),
  elbe_gpkg = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg",
  elbe_layer = "elbe",
  output_dir = file.path(root, "ANALYSE_FINAL_CORRIDOR/outputs/results_figures")
)

dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)

save_result_plot <- function(plot, filename, width = 9, height = 6) {
  ggsave(
    filename = file.path(paths$output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white",
    limitsize = FALSE
  )
}

map_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_line(color = "grey90", linewidth = 0.2),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 0),
      strip.text = element_text(face = "bold")
    )
}

plot_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 0),
      legend.position = "bottom"
    )
}

message("Loading final analysis layer ...")
corridor <- st_read(paths$analysis_gpkg, layer = paths$analysis_layer, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

message("Loading protection table ...")
protection <- read_csv(
  paths$protection_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

message("Loading Elbe geometry ...")
elbe <- st_read(paths$elbe_gpkg, layer = paths$elbe_layer, quiet = TRUE)
if (st_crs(elbe) != st_crs(corridor)) {
  elbe <- st_transform(elbe, st_crs(corridor))
}

analysis <- corridor %>%
  left_join(
    protection %>%
      select(
        AGS,
        protection_status,
        loss_portfolio_status,
        no_simulated_loss_event,
        annual_loss_probability_model,
        protection_return_period
      ),
    by = "AGS"
  ) %>%
  mutate(
    vuln_quintile_num = if_else(
      is.na(vuln_index_main_z),
      NA_integer_,
      ntile(vuln_index_main_z, 5)
    ),
    vuln_quintile = factor(
      vuln_quintile_num,
      levels = 1:5,
      labels = c("Q1 lowest", "Q2", "Q3", "Q4", "Q5 highest")
    ),
    vulnerability_tercile = factor(
      ntile(vuln_index_main_z, 3),
      levels = 1:3,
      labels = c("Low vulnerability", "Medium vulnerability", "High vulnerability")
    ),
    exposure_tercile = factor(
      ntile(flood_share_rp100, 3),
      levels = 1:3,
      labels = c("Low exposure", "Medium exposure", "High exposure")
    ),
    bivariate_class = paste(exposure_tercile, vulnerability_tercile, sep = " / "),
    bivariate_class = if_else(
      is.na(vulnerability_tercile),
      "No vulnerability data",
      bivariate_class
    ),
    annual_loss_class = case_when(
      annual_loss_probability_model == 0 ~ "0 no simulated loss",
      annual_loss_probability_model <= 0.01 ~ ">0-1%",
      annual_loss_probability_model <= 0.02 ~ "1-2%",
      annual_loss_probability_model <= 0.05 ~ "2-5%",
      annual_loss_probability_model <= 0.10 ~ "5-10%",
      annual_loss_probability_model > 0.10 ~ ">10%",
      TRUE ~ NA_character_
    ),
    annual_loss_class = factor(
      annual_loss_class,
      levels = c("0 no simulated loss", ">0-1%", "1-2%", "2-5%", "5-10%", ">10%")
    )
  )

analysis_tbl <- analysis %>% st_drop_geometry()

# ---------------------------
# 1) Vulnerability dimension map
# ---------------------------

dimension_long <- analysis %>%
  select(
    AGS,
    mun_name,
    vuln_dim_demographic_household_z,
    vuln_dim_deprivation_labour_z,
    vuln_dim_access_adaptive_capacity_z
  ) %>%
  pivot_longer(
    cols = starts_with("vuln_dim_"),
    names_to = "dimension",
    values_to = "score_z"
  ) %>%
  mutate(
    dimension = recode(
      dimension,
      vuln_dim_demographic_household_z = "Demographic and household structure",
      vuln_dim_deprivation_labour_z = "Deprivation and labour market",
      vuln_dim_access_adaptive_capacity_z = "Access and adaptive capacity"
    )
  )

map_vulnerability_dimensions <- ggplot(dimension_long) +
  geom_sf(aes(fill = score_z), color = NA) +
  geom_sf(data = elbe, color = "black", linewidth = 0.25, inherit.aes = FALSE) +
  facet_wrap(~dimension, ncol = 1) +
  scale_fill_viridis_c(option = "C", name = "Score (z)") +
  labs(
    title = "Dimensions of the socio-economic vulnerability index",
    subtitle = "Direction-coded dimension scores before final equal-weight aggregation",
    caption = "Data: INKAR (BBSR), EFAS/JRC-derived RP500 corridor. Own processing."
  ) +
  map_theme()

save_result_plot(map_vulnerability_dimensions, "map_vulnerability_dimensions.png", width = 8.5, height = 12)

# ---------------------------
# 2) Exposure-vulnerability plots
# ---------------------------

plot_vuln_vs_rp100 <- ggplot(analysis_tbl, aes(vuln_index_main_z, flood_share_rp100)) +
  geom_point(alpha = 0.45, size = 1.6, color = "#31572c") +
  geom_smooth(method = "lm", se = TRUE, color = "#1f2937", linewidth = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "RP100 flood exposure and socio-economic vulnerability",
    subtitle = "Municipality-level comparison within the RP500 Elbe corridor",
    x = "Socio-economic vulnerability index (z)",
    y = "RP100 flood share",
    caption = "Each point represents one corridor municipality. One municipality lacks INKAR data and is omitted."
  ) +
  plot_theme()

save_result_plot(plot_vuln_vs_rp100, "plot_vulnerability_vs_rp100_exposure.png", width = 8.5, height = 6)

plot_rp100_by_quintile <- ggplot(analysis_tbl, aes(vuln_quintile, flood_share_rp100)) +
  geom_boxplot(fill = "#d8e2dc", color = "#1f2937", outlier.alpha = 0.35) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "#b91c1c", color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "RP100 exposure by vulnerability quintile",
    subtitle = "Red dots show quintile means; boxes show municipal distribution",
    x = "Vulnerability quintile",
    y = "RP100 flood share",
    caption = "Data: final corridor exposure table and final vulnerability index. Own processing."
  ) +
  plot_theme()

save_result_plot(plot_rp100_by_quintile, "plot_rp100_exposure_by_vulnerability_quintile.png", width = 8.5, height = 6)

# ---------------------------
# 3) Bivariate exposure-vulnerability map
# ---------------------------

bivar_levels <- c(
  "Low exposure / Low vulnerability",
  "Medium exposure / Low vulnerability",
  "High exposure / Low vulnerability",
  "Low exposure / Medium vulnerability",
  "Medium exposure / Medium vulnerability",
  "High exposure / Medium vulnerability",
  "Low exposure / High vulnerability",
  "Medium exposure / High vulnerability",
  "High exposure / High vulnerability",
  "No vulnerability data"
)

bivar_colors <- c(
  "Low exposure / Low vulnerability" = "#e8e8e8",
  "Medium exposure / Low vulnerability" = "#ace4e4",
  "High exposure / Low vulnerability" = "#5ac8c8",
  "Low exposure / Medium vulnerability" = "#dfb0d6",
  "Medium exposure / Medium vulnerability" = "#a5add3",
  "High exposure / Medium vulnerability" = "#5698b9",
  "Low exposure / High vulnerability" = "#be64ac",
  "Medium exposure / High vulnerability" = "#8c62aa",
  "High exposure / High vulnerability" = "#3b4994",
  "No vulnerability data" = "#808080"
)

analysis <- analysis %>%
  mutate(bivariate_class = factor(bivariate_class, levels = bivar_levels))

map_exposure_vulnerability_bivariate <- ggplot(analysis) +
  geom_sf(aes(fill = bivariate_class), color = NA) +
  geom_sf(data = elbe, color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_manual(values = bivar_colors, drop = FALSE, name = "RP100 exposure / vulnerability") +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE)) +
  labs(
    title = "Bivariate pattern of RP100 exposure and socio-economic vulnerability",
    subtitle = "Tercile classification within the RP500 Elbe corridor",
    caption = "High-high municipalities combine high RP100 flood share and high socio-economic vulnerability. Data: EFAS/JRC, INKAR. Own processing."
  ) +
  map_theme() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.margin = margin(10, 20, 12, 20)
  )

save_result_plot(map_exposure_vulnerability_bivariate, "map_exposure_vulnerability_bivariate.png", width = 11.5, height = 8.7)

# ---------------------------
# 4) Multi-return-period exposure curves by vulnerability
# ---------------------------

curve_long <- analysis_tbl %>%
  select(
    AGS,
    vuln_quintile,
    flood_share_rp10,
    flood_share_rp20,
    flood_share_rp50,
    flood_share_rp100,
    flood_share_rp200,
    flood_share_rp500
  ) %>%
  filter(!is.na(vuln_quintile)) %>%
  pivot_longer(
    cols = starts_with("flood_share_rp"),
    names_to = "return_period",
    values_to = "flood_share"
  ) %>%
  mutate(
    return_period = as.numeric(gsub("flood_share_rp", "", return_period))
  ) %>%
  group_by(vuln_quintile, return_period) %>%
  summarise(
    mean_flood_share = mean(flood_share, na.rm = TRUE),
    median_flood_share = median(flood_share, na.rm = TRUE),
    .groups = "drop"
  )

plot_exposure_curves_by_vulnerability <- ggplot(
  curve_long,
  aes(return_period, mean_flood_share, color = vuln_quintile, group = vuln_quintile)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(10, 20, 50, 100, 200, 500)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_viridis_d(option = "C", end = 0.9, name = "Vulnerability") +
  labs(
    title = "Mean flood-exposure curves by vulnerability quintile",
    subtitle = "Average municipal flood share across return periods",
    x = "Return period",
    y = "Mean flood share",
    caption = "Data: multi-return-period exposure table and final vulnerability index. Own processing."
  ) +
  plot_theme()

save_result_plot(plot_exposure_curves_by_vulnerability, "plot_exposure_curves_by_vulnerability_quintile.png", width = 8.5, height = 6)

# ---------------------------
# 5) Protection/loss portfolio patterns
# ---------------------------

portfolio_summary <- analysis_tbl %>%
  filter(!is.na(vuln_quintile), !is.na(loss_portfolio_status)) %>%
  count(vuln_quintile, loss_portfolio_status) %>%
  group_by(vuln_quintile) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

plot_portfolio_status <- ggplot(portfolio_summary, aes(vuln_quintile, share, fill = loss_portfolio_status)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      "positive_modeled_loss" = "#1967a3",
      "zero_modeled_loss_probability" = "#d8d3c7"
    ),
    labels = c(
      "positive_modeled_loss" = "Positive modeled loss / finite RP",
      "zero_modeled_loss_probability" = "No simulated loss event"
    ),
    name = NULL
  ) +
  labs(
    title = "Provider loss-portfolio status by vulnerability quintile",
    subtitle = "No-event municipalities are interpreted as zero modeled annual loss probability",
    x = "Vulnerability quintile",
    y = "Share of municipalities",
    caption = "Data: final vulnerability index and provider protection/loss table. Own processing."
  ) +
  plot_theme()

save_result_plot(plot_portfolio_status, "plot_protection_status_by_vulnerability_quintile.png", width = 8.5, height = 6)

plot_annual_loss_by_quintile <- ggplot(
  analysis_tbl %>% filter(!is.na(vuln_quintile), !is.na(annual_loss_probability_model)),
  aes(vuln_quintile, annual_loss_probability_model)
) +
  geom_boxplot(fill = "#fee8c8", color = "#1f2937", outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "#b91c1c", color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Modeled annual loss probability by vulnerability quintile",
    subtitle = "No-event municipalities are retained as zero modeled annual loss probability",
    x = "Vulnerability quintile",
    y = "Annual loss probability",
    caption = "Data: provider loss/protection table and final vulnerability index. Own processing."
  ) +
  plot_theme()

save_result_plot(plot_annual_loss_by_quintile, "plot_annual_loss_probability_by_vulnerability_quintile.png", width = 8.5, height = 6)

plot_rp_by_quintile <- ggplot(
  analysis_tbl %>% filter(!is.na(vuln_quintile), !is.na(protection_return_period)),
  aes(vuln_quintile, protection_return_period)
) +
  geom_boxplot(fill = "#dbeafe", color = "#1f2937", outlier.alpha = 0.35) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.6, fill = "#b91c1c", color = "white") +
  scale_y_log10(labels = comma_format()) +
  labs(
    title = "Finite protection return periods by vulnerability quintile",
    subtitle = "Positive-loss municipalities only; red diamonds show medians",
    x = "Vulnerability quintile",
    y = "Protection return period (log scale)",
    caption = "Municipalities without simulated loss events are not assigned an artificial finite return period."
  ) +
  plot_theme()

save_result_plot(plot_rp_by_quintile, "plot_protection_return_period_by_vulnerability_quintile.png", width = 8.5, height = 6)

map_annual_loss_probability <- ggplot(analysis) +
  geom_sf(aes(fill = annual_loss_class), color = NA) +
  geom_sf(data = elbe, color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(
      "0 no simulated loss" = "#e5e7eb",
      ">0-1%" = "#ffffcc",
      "1-2%" = "#ffeda0",
      "2-5%" = "#feb24c",
      "5-10%" = "#f03b20",
      ">10%" = "#bd0026"
    ),
    drop = FALSE,
    name = "Annual loss probability"
  ) +
  labs(
    title = "Modeled annual loss probability in the RP500 Elbe corridor",
    subtitle = "No-event municipalities are mapped as zero modeled annual loss probability",
    caption = "Data: provider loss/protection table joined to corridor municipalities. Own processing."
  ) +
  map_theme()

save_result_plot(map_annual_loss_probability, "map_modeled_annual_loss_probability.png", width = 10.5, height = 8)

# ---------------------------
# 6) Summary tables for Results writing
# ---------------------------

write_csv(
  portfolio_summary,
  file.path(paths$output_dir, "table_portfolio_status_by_vulnerability_quintile.csv")
)

write_csv(
  curve_long,
  file.path(paths$output_dir, "table_mean_exposure_curves_by_vulnerability_quintile.csv")
)

result_summary <- tibble(
  metric = c(
    "municipalities",
    "municipalities_with_vulnerability_index",
    "corr_vulnerability_rp100",
    "corr_main_index_pca23_70pct",
    "corr_pca23_70pct_pca23_kaiser",
    "positive_loss_municipalities",
    "zero_modeled_loss_probability_municipalities"
  ),
  value = c(
    nrow(analysis_tbl),
    sum(!is.na(analysis_tbl$vuln_index_main_z)),
    cor(analysis_tbl$vuln_index_main_z, analysis_tbl$flood_share_rp100, use = "complete.obs"),
    cor(analysis_tbl$vuln_index_main_z, analysis_tbl$vuln_index_pca23_70pct_z, use = "complete.obs"),
    cor(analysis_tbl$vuln_index_pca23_70pct_z, analysis_tbl$vuln_index_pca23_kaiser_z, use = "complete.obs"),
    sum(analysis_tbl$loss_portfolio_status == "positive_modeled_loss", na.rm = TRUE),
    sum(analysis_tbl$loss_portfolio_status == "zero_modeled_loss_probability", na.rm = TRUE)
  )
)

write_csv(
  result_summary,
  file.path(paths$output_dir, "table_results_key_summary.csv")
)

message("Results figures written to: ", paths$output_dir)
