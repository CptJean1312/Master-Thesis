#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(ggspatial)
})

options(scipen = 999)
set.seed(42)

# ============================================================
# ANALYSE FINAL CORRIDOR
# ------------------------------------------------------------
# Wide PCA on the EFAS-defined RP500 corridor only
# - same wide PCA logic as the original clean analysis
# - sample restricted to corridor municipalities
# - hazard focus for now: RP100 only
# ============================================================

# ---------------------------
# 1) User-editable input paths
# ---------------------------

corridor_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg"
corridor_layer <- "municipalities_corridor"

exposure_csv <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv"

basin_municipalities_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/gemeinden_elbe_landonly_basin_inkar.gpkg"
basin_municipalities_layer <- "gemeinden_elbe_landonly_basin_inkar"

germany_municipalities_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/vg250_gemeinden_landonly.gpkg"
germany_municipalities_layer <- "vg250_gemeinden_landonly"

elbe_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"
elbe_layer <- "elbe"

# ---------------------------
# 2) Output structure
# ---------------------------

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else file.path(getwd(), "ANALYSE_FINAL_CORRIDOR.R")
script_dir <- normalizePath(dirname(script_file), winslash = "/")
output_root <- file.path(script_dir, "outputs")
plot_dir <- file.path(output_root, "plots")
map_dir <- file.path(output_root, "maps")
table_dir <- file.path(output_root, "tables")
gpkg_dir <- file.path(output_root, "gpkg")
log_dir <- file.path(output_root, "logs")

for (dir_path in c(output_root, plot_dir, map_dir, table_dir, gpkg_dir, log_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "processing_log.txt")
if (file.exists(log_file)) unlink(log_file)

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

save_plot <- function(plot_obj, filename, width = 8, height = 5, subdir = "plots") {
  ggsave(
    filename = file.path(output_root, subdir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 320,
    bg = "white",
    limitsize = FALSE
  )
}

save_table <- function(data_obj, filename) {
  write_csv(data_obj, file.path(table_dir, filename))
}

map_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.caption = element_text(size = base_size - 3, hjust = 1)
    )
}

map_annotations <- function() {
  list(
    annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = north_arrow_fancy_orienteering,
      height = unit(1.2, "cm"),
      width = unit(1.2, "cm")
    ),
    annotation_scale(
      location = "bl",
      width_hint = 0.25,
      line_width = 0.6,
      text_cex = 0.7
    )
  )
}

impute_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

build_weighted_index <- function(spatial_data, n_components, variance_table) {
  weights <- variance_table$Variance[seq_len(n_components)]
  names(weights) <- paste0("PC", seq_len(n_components))

  pc_matrix <- spatial_data %>%
    st_drop_geometry() %>%
    select(all_of(names(weights))) %>%
    as.matrix()

  raw_index <- as.numeric(pc_matrix %*% weights)
  z_index <- as.numeric(scale(raw_index))

  list(raw = raw_index, z = z_index)
}

write_gpkg_layer <- function(x, dsn, layer, replace_dsn = FALSE) {
  if (replace_dsn && file.exists(dsn)) {
    unlink(dsn)
  }

  if (!file.exists(dsn)) {
    write_sf(x, dsn = dsn, layer = layer, quiet = TRUE)
  } else {
    write_sf(x, dsn = dsn, layer = layer, append = TRUE, quiet = TRUE)
  }
}

# ---------------------------
# 3) Load corridor + exposure
# ---------------------------

log_message("Loading corridor municipalities ...")
corridor_sf <- st_read(corridor_path, layer = corridor_layer, quiet = TRUE)

if (st_crs(corridor_sf)$epsg != 25832) {
  corridor_sf <- st_transform(corridor_sf, 25832)
}

log_message("Corridor municipalities loaded: ", nrow(corridor_sf))

log_message("Loading municipal exposure table ...")
exposure_tbl <- read_csv(
  exposure_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

log_message("Exposure rows loaded: ", nrow(exposure_tbl))

corridor_sf <- corridor_sf %>%
  mutate(AGS = as.character(AGS)) %>%
  left_join(select(exposure_tbl, -mun_name), by = "AGS")

if (all(c("municipality_area_m2.x", "municipality_area_m2.y") %in% names(corridor_sf))) {
  corridor_sf <- corridor_sf %>%
    mutate(municipality_area_m2 = coalesce(municipality_area_m2.y, municipality_area_m2.x)) %>%
    select(-municipality_area_m2.x, -municipality_area_m2.y)
}

missing_rp100 <- sum(is.na(corridor_sf$flood_share_rp100))
if (missing_rp100 > 0) {
  stop("Missing RP100 exposure after join for ", missing_rp100, " corridor municipalities.")
}

log_message("RP100 exposure joined successfully.")

# ---------------------------
# 4) Wide PCA inputs
# ---------------------------

pca_columns <- c(
  "share_alg2_sgb2",
  "share_bg_single_parent",
  "share_bg_5plus",
  "share_bg_with_children",
  "share_sgb2_with_housing_costs",
  "share_longterm_unemp",
  "unemp_u25_per_1000",
  "unemp_55plus_per_1000",
  "income_tax_per_capita",
  "purchasing_power",
  "trade_tax_per_capita",
  "tax_revenue_total",
  "share_age_0_3",
  "share_age_3_6",
  "share_age_6_18",
  "share_age_18_25",
  "share_age_25_30",
  "share_age_30_50",
  "share_age_50_65",
  "share_age_65_75",
  "share_age_65plus",
  "share_age_75plus",
  "natural_pop_change",
  "migration_balance",
  "old_age_dependency",
  "youth_dependency",
  "young_old_ratio",
  "share_single_households",
  "share_households_with_children",
  "share_hh_income_high",
  "share_hh_income_medium",
  "share_hh_income_low",
  "students_total_per_1000",
  "students_18_25_per_1000",
  "students_fh_per_1000",
  "gp_general_per_1000",
  "gp_primary_per_1000",
  "internists_per_1000",
  "pediatricians_per_1000_children",
  "doctors_total_per_1000",
  "share_bb_1000mbit",
  "share_bb_100mbit",
  "share_bb_50mbit",
  "share_4g",
  "dist_supermarket_m",
  "dist_pharmacy_m",
  "dist_gp_m",
  "dist_public_transport_m",
  "dist_primary_school_m",
  "pop_density_per_km2",
  "employment_density_per_km2"
)

missing_pca_columns <- setdiff(pca_columns, names(corridor_sf))
if (length(missing_pca_columns) > 0) {
  stop("Missing PCA columns: ", paste(missing_pca_columns, collapse = ", "))
}

wide_pca_data <- corridor_sf %>%
  st_drop_geometry() %>%
  select(all_of(pca_columns)) %>%
  mutate(across(everything(), as.numeric))

missingness_tbl <- tibble(
  variable = names(wide_pca_data),
  n_missing = colSums(is.na(wide_pca_data)),
  share_missing = colMeans(is.na(wide_pca_data))
) %>%
  arrange(desc(share_missing), variable)

save_table(missingness_tbl, "corridor_pca_missingness.csv")

wide_pca_imputed <- wide_pca_data %>%
  mutate(across(everything(), impute_median))

wide_pca_scaled <- scale(wide_pca_imputed)
wide_pca <- prcomp(wide_pca_scaled)

component_scores <- as.data.frame(wide_pca$x)
names(component_scores) <- paste0("PC", seq_len(ncol(component_scores)))
corridor_sf <- bind_cols(corridor_sf, component_scores)

eigenvalues <- wide_pca$sdev^2
variance_share <- eigenvalues / sum(eigenvalues)

scree_table <- tibble(
  PC = seq_along(eigenvalues),
  Eigenvalue = eigenvalues,
  Variance = variance_share,
  Cumulative = cumsum(variance_share)
)

save_table(scree_table, "corridor_scree_table.csv")

scree_plot <- ggplot(scree_table, aes(x = PC, y = Eigenvalue)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Scree Plot (Kaiser)", x = "PC", y = "Eigenvalue") +
  theme_classic()

cumulative_variance_plot <- ggplot(scree_table, aes(x = PC, y = Cumulative)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  labs(title = "Cumulative Variance", x = "PC", y = "Cumulative") +
  theme_classic()

save_plot(scree_plot, "scree_kaiser_corridor.png", width = 8, height = 5, subdir = "plots")
save_plot(cumulative_variance_plot, "cumulative_variance_corridor.png", width = 8, height = 5, subdir = "plots")

loading_table <- as.data.frame(wide_pca$rotation)
loading_table$variable <- rownames(loading_table)

top_loadings <- loading_table %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") %>%
  mutate(abs_loading = abs(loading)) %>%
  group_by(PC) %>%
  arrange(desc(abs_loading)) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  arrange(as.integer(str_remove(PC, "PC")), desc(abs_loading))

save_table(top_loadings, "corridor_pca_top_loadings_top8_per_pc.csv")
save_table(tibble(variable = pca_columns), "corridor_wide_pca_variables.csv")

# ---------------------------
# 5) Vulnerability index
# ---------------------------

main_index <- build_weighted_index(corridor_sf, n_components = 8, variance_table = scree_table)
corridor_sf$vuln_index_main <- main_index$raw
corridor_sf$vuln_index_main_z <- main_index$z

sensitivity_index <- build_weighted_index(corridor_sf, n_components = 12, variance_table = scree_table)
corridor_sf$vuln_index_sens12 <- sensitivity_index$raw
corridor_sf$vuln_index_sens12_z <- sensitivity_index$z

anchor_correlation <- suppressWarnings(
  cor(corridor_sf$vuln_index_main_z, wide_pca_imputed$share_alg2_sgb2, use = "complete.obs")
)

if (!is.na(anchor_correlation) && anchor_correlation < 0) {
  corridor_sf$vuln_index_main <- -corridor_sf$vuln_index_main
  corridor_sf$vuln_index_main_z <- -corridor_sf$vuln_index_main_z
  corridor_sf$vuln_index_sens12 <- -corridor_sf$vuln_index_sens12
  corridor_sf$vuln_index_sens12_z <- -corridor_sf$vuln_index_sens12_z
  log_message("Index flipped so that higher values indicate higher vulnerability.")
} else {
  log_message("Index direction OK.")
}

index_summary <- corridor_sf %>%
  st_drop_geometry() %>%
  summarise(
    municipalities = n(),
    rp100_mean = mean(flood_share_rp100, na.rm = TRUE),
    rp100_median = median(flood_share_rp100, na.rm = TRUE),
    vuln_mean = mean(vuln_index_main_z, na.rm = TRUE),
    vuln_sd = sd(vuln_index_main_z, na.rm = TRUE),
    corr_vuln_rp100 = cor(vuln_index_main_z, flood_share_rp100, use = "complete.obs")
  )

save_table(index_summary, "corridor_rp100_index_summary.csv")

# ---------------------------
# 6) Final analysis layer
# ---------------------------

analysis_table <- corridor_sf %>%
  st_drop_geometry() %>%
  select(
    AGS,
    mun_name,
    municipality_area_m2,
    starts_with("flood_area_rp"),
    starts_with("flood_share_rp"),
    starts_with("PC"),
    vuln_index_main,
    vuln_index_main_z,
    vuln_index_sens12,
    vuln_index_sens12_z
  )

save_table(analysis_table, "corridor_analysis_rp100.csv")

final_gpkg <- file.path(gpkg_dir, "corridor_wide_pca_rp100_analysis.gpkg")
write_gpkg_layer(corridor_sf, final_gpkg, layer = "corridor_wide_pca_rp100", replace_dsn = TRUE)

saveRDS(corridor_sf, file.path(output_root, "corridor_wide_pca_rp100.rds"))

# ---------------------------
# 7) Maps
# ---------------------------

log_message("Building study area map ...")

germany_vect <- vect(germany_municipalities_path, layer = germany_municipalities_layer)
germany_outline_sf <- st_as_sf(aggregate(germany_vect))

elbe_sf <- st_read(elbe_path, layer = elbe_layer, quiet = TRUE)
elbe_sf <- st_zm(elbe_sf, drop = TRUE, what = "ZM")
elbe_sf <- st_make_valid(elbe_sf)
if (st_crs(elbe_sf) != st_crs(corridor_sf)) {
  elbe_sf <- st_transform(elbe_sf, st_crs(corridor_sf))
}

basin_sf <- st_read(basin_municipalities_path, layer = basin_municipalities_layer, quiet = TRUE)
if (st_crs(basin_sf) != st_crs(corridor_sf)) {
  basin_sf <- st_transform(basin_sf, st_crs(corridor_sf))
}

study_area_map <- ggplot() +
  geom_sf(data = germany_outline_sf, fill = "#f4f1eb", color = "grey70", linewidth = 0.25) +
  geom_sf(data = basin_sf, fill = "#d6dde2", color = NA, alpha = 0.55) +
  geom_sf(data = corridor_sf, fill = "#0f766e", color = NA, alpha = 0.88) +
  geom_sf(data = elbe_sf, color = "#0b3954", linewidth = 0.45) +
  coord_sf(xlim = c(350000, 950000), ylim = c(5220000, 6105000), expand = FALSE) +
  map_annotations() +
  labs(
    title = "Study area: municipalities intersecting the RP500 flood corridor",
    subtitle = "EFAS-derived corridor municipalities within the Elbe basin",
    caption = "Background: Germany outline. Light grey: Elbe-basin municipalities. Green: corridor municipalities. Line: Elbe."
  ) +
  map_theme()

save_plot(study_area_map, "map_study_area_corridor.png", width = 9.5, height = 8, subdir = "maps")

map_vulnerability <- ggplot(corridor_sf) +
  geom_sf(aes(fill = vuln_index_main_z), color = NA) +
  geom_sf(data = elbe_sf, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  labs(
    title = "Socio-economic vulnerability index in the EFAS corridor",
    subtitle = "Wide PCA recalculated on the RP500 corridor municipalities only",
    fill = "Index (z)",
    caption = "Source: INKAR (BBSR), EFAS corridor definition. Own processing."
  ) +
  scale_fill_viridis_c(option = "C") +
  map_theme()

save_plot(map_vulnerability, "map_vulnerability_index_corridor.png", width = 9, height = 7, subdir = "maps")

map_rp100 <- ggplot(corridor_sf) +
  geom_sf(aes(fill = flood_share_rp100), color = NA) +
  geom_sf(data = elbe_sf, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  labs(
    title = "Flood exposure in the EFAS corridor (RP100)",
    subtitle = "Share of municipality area with valid RP100 flood-depth pixels",
    fill = "Share",
    caption = "Source: EFAS flood rasters. Own processing."
  ) +
  scale_fill_viridis_c(option = "C") +
  map_theme()

save_plot(map_rp100, "map_rp100_exposure_corridor.png", width = 9, height = 7, subdir = "maps")

pc_titles <- c(
  PC1 = "PC1: Socio-economic disadvantage and welfare dependency",
  PC2 = "PC2: Urbanisation, density and accessibility",
  PC3 = "PC3: Demographic ageing and dependency",
  PC4 = "PC4: Education, students and human capital"
)

pc_subtitles <- c(
  PC1 = "High loadings: ALG II, low income, unemployment, single households",
  PC2 = "High loadings: population density, transport access, broadband",
  PC3 = "High loadings: 65+, old-age dependency, low youth share",
  PC4 = "High loadings: students, higher income, education indicators"
)

for (pc_name in paste0("PC", 1:4)) {
  pc_map <- ggplot(corridor_sf) +
    geom_sf(aes(fill = .data[[pc_name]]), color = NA) +
    geom_sf(data = elbe_sf, inherit.aes = FALSE, color = "black", linewidth = 0.3) +
    coord_sf(crs = 25832) +
    map_annotations() +
    labs(
      title = pc_titles[[pc_name]],
      subtitle = pc_subtitles[[pc_name]],
      fill = pc_name,
      caption = "Source: INKAR (BBSR). Corridor-specific PCA. Own processing."
    ) +
    scale_fill_viridis_c(option = "C") +
    map_theme()

  save_plot(pc_map, paste0("map_", pc_name, "_corridor.png"), width = 9, height = 7, subdir = "maps")
}

log_message("Corridor PCA pipeline finished successfully.")
log_message("Final analysis GPKG: ", final_gpkg)
