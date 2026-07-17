#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

options(scipen = 999)

# Result-facing diagnostics for the DLR land-use exposure extension.
# This combines land-use exposure, the final vulnerability index, and
# protection/loss outputs so Chapter 5 can use one coherent analysis table.

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"

paths <- list(
  landuse_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/gpkg/corridor_landuse_exposure.gpkg"),
  landuse_layer = "corridor_landuse_exposure",
  protection_csv = file.path(root, "ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/corridor_protection_level_all_corridor_municipalities.csv"),
  large_city_csv = file.path(root, "ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/large_municipalities_portfolio_status.csv"),
  elbe_gpkg = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg",
  elbe_layer = "elbe",
  table_dir = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/tables"),
  gpkg_dir = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/gpkg"),
  plot_dir = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/plots")
)

for (dir_path in c(paths$table_dir, paths$gpkg_dir, paths$plot_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

plot_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 0)
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
      plot.caption = element_text(size = base_size - 3, hjust = 0)
    )
}

save_plot <- function(plot, filename, width = 8.5, height = 6) {
  ggsave(
    filename = file.path(paths$plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white",
    limitsize = FALSE
  )
}

make_tercile <- function(x, labels) {
  q <- quantile(x, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
  factor(
    case_when(
      is.na(x) ~ NA_character_,
      x <= q[1] ~ labels[1],
      x <= q[2] ~ labels[2],
      TRUE ~ labels[3]
    ),
    levels = labels
  )
}

safe_cor <- function(x, y, method) {
  if (sum(complete.cases(x, y)) < 3) {
    return(NA_real_)
  }
  cor(x, y, use = "complete.obs", method = method)
}

message("Loading land-use exposure layer ...")
landuse <- st_read(paths$landuse_gpkg, layer = paths$landuse_layer, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

if (st_crs(landuse)$epsg != 25832) {
  landuse <- st_transform(landuse, 25832)
}

message("Loading protection/loss table ...")
protection <- read_csv(paths$protection_csv, col_types = cols(AGS = col_character()), show_col_types = FALSE) %>%
  select(
    AGS,
    pop_total,
    protection_available,
    protection_status,
    loss_portfolio_status,
    no_simulated_loss_event,
    annual_loss_probability_model,
    protection_return_period,
    n_nonzero_years
  )

large_cities <- tibble::tibble(AGS = character(), city_label = character())
if (file.exists(paths$large_city_csv)) {
  large_cities <- read_csv(paths$large_city_csv, col_types = cols(AGS = col_character()), show_col_types = FALSE) %>%
    select(AGS, city_label)
}

analysis <- landuse %>%
  left_join(protection, by = "AGS", suffix = c("", "_protection")) %>%
  left_join(large_cities, by = "AGS") %>%
  mutate(
    vuln_quintile_num = if_else(is.na(vuln_index_main_z), NA_integer_, ntile(vuln_index_main_z, 5)),
    vuln_quintile = factor(
      vuln_quintile_num,
      levels = 1:5,
      labels = c("Q1 lowest", "Q2", "Q3", "Q4", "Q5 highest")
    ),
    vulnerability_tercile = make_tercile(
      vuln_index_main_z,
      c("Low vulnerability", "Medium vulnerability", "High vulnerability")
    ),
    total_area_exposure_tercile = make_tercile(
      flood_share_rp100,
      c("Low total-area exposure", "Medium total-area exposure", "High total-area exposure")
    ),
    artificial_exposure_tercile = make_tercile(
      rp100_artificial_flood_share_of_group,
      c("Low artificial-land exposure", "Medium artificial-land exposure", "High artificial-land exposure")
    ),
    high_total_area_exposure = total_area_exposure_tercile == "High total-area exposure",
    high_artificial_exposure = artificial_exposure_tercile == "High artificial-land exposure",
    high_vulnerability = vulnerability_tercile == "High vulnerability",
    high_high_artificial_vulnerability = high_artificial_exposure & high_vulnerability,
    bivariate_artificial_vulnerability = case_when(
      is.na(vulnerability_tercile) | is.na(artificial_exposure_tercile) ~ "Missing",
      TRUE ~ paste(artificial_exposure_tercile, vulnerability_tercile, sep = " / ")
    )
  )

analysis_tbl <- analysis %>% st_drop_geometry()

write_csv(
  analysis_tbl,
  file.path(paths$table_dir, "corridor_landuse_protection_analysis_rp100.csv")
)

out_gpkg <- file.path(paths$gpkg_dir, "corridor_landuse_protection_analysis_rp100.gpkg")
if (file.exists(out_gpkg)) {
  unlink(out_gpkg)
}
st_write(
  analysis,
  dsn = out_gpkg,
  layer = "corridor_landuse_protection_analysis_rp100",
  quiet = TRUE
)

correlation_comparison <- tibble::tibble(
  relationship = c(
    "vulnerability_vs_rp100_total_area_flood_share",
    "vulnerability_vs_rp100_artificial_land_flood_share",
    "access_adaptive_capacity_vs_rp100_artificial_land_flood_share",
    "demographic_household_vs_rp100_artificial_land_flood_share",
    "deprivation_labour_vs_rp100_artificial_land_flood_share",
    "annual_loss_probability_vs_rp100_total_area_flood_share",
    "annual_loss_probability_vs_rp100_artificial_land_flood_share",
    "finite_protection_rp_vs_rp100_artificial_land_flood_share_positive_loss_only"
  ),
  pearson = c(
    safe_cor(analysis_tbl$vuln_index_main_z, analysis_tbl$flood_share_rp100, "pearson"),
    safe_cor(analysis_tbl$vuln_index_main_z, analysis_tbl$rp100_artificial_flood_share_of_group, "pearson"),
    safe_cor(analysis_tbl$vuln_dim_access_adaptive_capacity_z, analysis_tbl$rp100_artificial_flood_share_of_group, "pearson"),
    safe_cor(analysis_tbl$vuln_dim_demographic_household_z, analysis_tbl$rp100_artificial_flood_share_of_group, "pearson"),
    safe_cor(analysis_tbl$vuln_dim_deprivation_labour_z, analysis_tbl$rp100_artificial_flood_share_of_group, "pearson"),
    safe_cor(analysis_tbl$annual_loss_probability_model, analysis_tbl$flood_share_rp100, "pearson"),
    safe_cor(analysis_tbl$annual_loss_probability_model, analysis_tbl$rp100_artificial_flood_share_of_group, "pearson"),
    safe_cor(analysis_tbl$protection_return_period[analysis_tbl$protection_available], analysis_tbl$rp100_artificial_flood_share_of_group[analysis_tbl$protection_available], "pearson")
  ),
  spearman = c(
    safe_cor(analysis_tbl$vuln_index_main_z, analysis_tbl$flood_share_rp100, "spearman"),
    safe_cor(analysis_tbl$vuln_index_main_z, analysis_tbl$rp100_artificial_flood_share_of_group, "spearman"),
    safe_cor(analysis_tbl$vuln_dim_access_adaptive_capacity_z, analysis_tbl$rp100_artificial_flood_share_of_group, "spearman"),
    safe_cor(analysis_tbl$vuln_dim_demographic_household_z, analysis_tbl$rp100_artificial_flood_share_of_group, "spearman"),
    safe_cor(analysis_tbl$vuln_dim_deprivation_labour_z, analysis_tbl$rp100_artificial_flood_share_of_group, "spearman"),
    safe_cor(analysis_tbl$annual_loss_probability_model, analysis_tbl$flood_share_rp100, "spearman"),
    safe_cor(analysis_tbl$annual_loss_probability_model, analysis_tbl$rp100_artificial_flood_share_of_group, "spearman"),
    safe_cor(analysis_tbl$protection_return_period[analysis_tbl$protection_available], analysis_tbl$rp100_artificial_flood_share_of_group[analysis_tbl$protection_available], "spearman")
  ),
  n_complete = c(
    sum(complete.cases(analysis_tbl$vuln_index_main_z, analysis_tbl$flood_share_rp100)),
    sum(complete.cases(analysis_tbl$vuln_index_main_z, analysis_tbl$rp100_artificial_flood_share_of_group)),
    sum(complete.cases(analysis_tbl$vuln_dim_access_adaptive_capacity_z, analysis_tbl$rp100_artificial_flood_share_of_group)),
    sum(complete.cases(analysis_tbl$vuln_dim_demographic_household_z, analysis_tbl$rp100_artificial_flood_share_of_group)),
    sum(complete.cases(analysis_tbl$vuln_dim_deprivation_labour_z, analysis_tbl$rp100_artificial_flood_share_of_group)),
    sum(complete.cases(analysis_tbl$annual_loss_probability_model, analysis_tbl$flood_share_rp100)),
    sum(complete.cases(analysis_tbl$annual_loss_probability_model, analysis_tbl$rp100_artificial_flood_share_of_group)),
    sum(complete.cases(analysis_tbl$protection_return_period[analysis_tbl$protection_available], analysis_tbl$rp100_artificial_flood_share_of_group[analysis_tbl$protection_available]))
  )
)

write_csv(
  correlation_comparison,
  file.path(paths$table_dir, "table_landuse_exposure_correlation_comparison.csv")
)

quintile_summary <- analysis_tbl %>%
  group_by(vuln_quintile) %>%
  summarise(
    municipalities = n(),
    mean_rp100_total_area_flood_share = mean(flood_share_rp100, na.rm = TRUE),
    median_rp100_total_area_flood_share = median(flood_share_rp100, na.rm = TRUE),
    mean_rp100_artificial_land_flood_share = mean(rp100_artificial_flood_share_of_group, na.rm = TRUE),
    median_rp100_artificial_land_flood_share = median(rp100_artificial_flood_share_of_group, na.rm = TRUE),
    mean_annual_loss_probability_model = mean(annual_loss_probability_model, na.rm = TRUE),
    positive_loss_share = mean(protection_available, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  quintile_summary,
  file.path(paths$table_dir, "table_artificial_exposure_by_vulnerability_quintile.csv")
)

high_exposure_overlap <- analysis_tbl %>%
  mutate(
    exposure_overlap_class = case_when(
      high_total_area_exposure & high_artificial_exposure ~ "High total-area and high artificial-land exposure",
      high_total_area_exposure & !high_artificial_exposure ~ "High total-area only",
      !high_total_area_exposure & high_artificial_exposure ~ "High artificial-land only",
      TRUE ~ "Neither high"
    )
  ) %>%
  count(exposure_overlap_class, name = "municipalities") %>%
  mutate(share_of_corridor = municipalities / sum(municipalities))

write_csv(
  high_exposure_overlap,
  file.path(paths$table_dir, "table_high_exposure_overlap.csv")
)

protection_by_artificial_tercile <- analysis_tbl %>%
  group_by(artificial_exposure_tercile) %>%
  summarise(
    municipalities = n(),
    positive_loss_municipalities = sum(protection_available, na.rm = TRUE),
    positive_loss_share = mean(protection_available, na.rm = TRUE),
    mean_annual_loss_probability_model = mean(annual_loss_probability_model, na.rm = TRUE),
    median_finite_protection_return_period = median(protection_return_period[protection_available], na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  protection_by_artificial_tercile,
  file.path(paths$table_dir, "table_protection_by_artificial_exposure_tercile.csv")
)

data_dictionary <- tibble::tribble(
  ~field, ~meaning,
  "rp100_artificial_flood_share_of_group", "Share of DLR Artificial Land within a municipality intersecting RP100 modeled flood cells.",
  "rp100_artificial_flooded_area_m2", "Area of DLR Artificial Land intersecting RP100 modeled flood cells.",
  "rp100_artificial_share_of_lc_flooded_area", "Share of all flooded land-cover area within a municipality that is DLR Artificial Land.",
  "vuln_quintile", "Quintile of the final socio-economic vulnerability index.",
  "artificial_exposure_tercile", "Tercile classification of RP100 artificial-land flood share.",
  "high_high_artificial_vulnerability", "TRUE if municipality is in the high artificial-land exposure tercile and high vulnerability tercile.",
  "protection_available", "TRUE if the provider loss/protection table contains a positive modeled loss and finite protection return-period estimate.",
  "annual_loss_probability_model", "Annual loss probability from the provider table; no-event municipalities retained as modeled zero following provider clarification."
)

write_csv(
  data_dictionary,
  file.path(paths$table_dir, "corridor_landuse_protection_data_dictionary.csv")
)

boxplot_artificial <- ggplot(
  analysis_tbl %>% filter(!is.na(vuln_quintile), !is.na(rp100_artificial_flood_share_of_group)),
  aes(vuln_quintile, rp100_artificial_flood_share_of_group)
) +
  geom_boxplot(fill = "#d8e2dc", color = "#1f2937", outlier.alpha = 0.35) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "#b91c1c", color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "RP100 artificial-land exposure by vulnerability quintile",
    subtitle = "Red dots show quintile means; boxes show municipal distribution",
    x = "Vulnerability quintile",
    y = "Flooded share of DLR Artificial Land",
    caption = "Data: DLR Land Cover DE 2015, EFAS/JRC RP100 flood raster, INKAR. Own processing."
  ) +
  plot_theme()

save_plot(boxplot_artificial, "plot_rp100_artificial_land_exposure_by_vulnerability_quintile.png")

message("Loading Elbe geometry for bivariate map ...")
elbe <- st_read(paths$elbe_gpkg, layer = paths$elbe_layer, quiet = TRUE)
if (st_crs(elbe)$epsg != 25832) {
  elbe <- st_transform(elbe, 25832)
}

bivar_colors <- c(
  "Low artificial-land exposure / Low vulnerability" = "#e8e8e8",
  "Medium artificial-land exposure / Low vulnerability" = "#ace4e4",
  "High artificial-land exposure / Low vulnerability" = "#5ac8c8",
  "Low artificial-land exposure / Medium vulnerability" = "#dfb0d6",
  "Medium artificial-land exposure / Medium vulnerability" = "#a5add3",
  "High artificial-land exposure / Medium vulnerability" = "#5698b9",
  "Low artificial-land exposure / High vulnerability" = "#be64ac",
  "Medium artificial-land exposure / High vulnerability" = "#8c62aa",
  "High artificial-land exposure / High vulnerability" = "#3b4994",
  "Missing" = "#eeeeee"
)

map_bbox <- st_bbox(analysis)

bivariate_map <- ggplot(analysis) +
  geom_sf(aes(fill = bivariate_artificial_vulnerability), color = "white", linewidth = 0.035) +
  geom_sf(data = elbe, color = "#0b4f8a", linewidth = 0.55, inherit.aes = FALSE) +
  scale_fill_manual(values = bivar_colors, drop = FALSE, name = "Artificial-land exposure / vulnerability") +
  coord_sf(
    xlim = c(unname(map_bbox[["xmin"]]), unname(map_bbox[["xmax"]])),
    ylim = c(unname(map_bbox[["ymin"]]), unname(map_bbox[["ymax"]])),
    expand = FALSE
  ) +
  labs(
    title = "RP100 artificial-land exposure and socio-economic vulnerability",
    subtitle = "Tercile-based bivariate classification for corridor municipalities",
    caption = "Data: DLR Land Cover DE 2015, EFAS/JRC RP100 flood raster, INKAR. Blue line: Elbe. Own processing."
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  map_theme(base_size = 10)

save_plot(bivariate_map, "map_artificial_exposure_vulnerability_bivariate.png", width = 10.5, height = 8)

message("Done.")
message("Analysis table: ", file.path(paths$table_dir, "corridor_landuse_protection_analysis_rp100.csv"))
message("Analysis GPKG: ", out_gpkg)
message("Correlation comparison: ", file.path(paths$table_dir, "table_landuse_exposure_correlation_comparison.csv"))
