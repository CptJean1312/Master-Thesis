#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

# Land-cover exposure module for the final Elbe corridor thesis workflow.
#
# Purpose:
# - add a land-use sensitive exposure metric to the municipality-level analysis;
# - distinguish area-based flood exposure from artificial/built land-cover exposure;
# - keep the large DLR Land Cover DE raster outside the Git repository.
#
# Method:
# - DLR Land Cover DE 2015 is a 10 m categorical raster in EPSG:3035.
# - EFAS/JRC flood rasters are 100 m rasters in EPSG:25832.
# - Each DLR class is converted to a binary mask and projected to the EFAS grid
#   with method = "average". The resulting raster is the fraction of each
#   100 m EFAS cell covered by that DLR class.
# - Class area per EFAS cell is fraction * cell area.
# - Flooded class area is class area masked by valid flood-depth cells.
# - Municipality-level sums are extracted with exact polygon extraction.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

paths <- list(
  landcover_tif = Sys.getenv(
    "DLR_LAND_COVER_TIF",
    unset = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/LANDUSE/data_raw/Land_Cover_DE_2015.tif"
  ),
  corridor_gpkg = Sys.getenv(
    "CORRIDOR_GPKG",
    unset = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg"
  ),
  exposure_csv = Sys.getenv(
    "EXPOSURE_CSV",
    unset = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv"
  ),
  analysis_csv = file.path(root, "ANALYSE_FINAL_CORRIDOR/outputs/tables/corridor_analysis_rp100.csv"),
  output_root = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs"),
  flood_raster_dir = Sys.getenv(
    "FLOOD_RASTER_DIR",
    unset = file.path(root, "outputs_eu_flood_25832")
  )
)

paths$table_dir <- file.path(paths$output_root, "tables")
paths$gpkg_dir <- file.path(paths$output_root, "gpkg")
paths$plot_dir <- file.path(paths$output_root, "plots")
paths$log_dir <- file.path(paths$output_root, "logs")
paths$cache_dir <- file.path(paths$output_root, "cache", "corridor_template_v1")

for (dir_path in c(paths$output_root, paths$table_dir, paths$gpkg_dir, paths$plot_dir, paths$log_dir, paths$cache_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(paths$log_dir, "build_landuse_exposure_log.txt")
if (file.exists(log_path)) {
  file.remove(log_path)
}

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_path, append = TRUE)
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

clean_slug <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

extract_sum <- function(raster, polygons_vect) {
  extracted <- terra::extract(
    raster,
    polygons_vect,
    fun = sum,
    na.rm = TRUE,
    exact = TRUE
  )
  values <- extracted[[2]]
  values[is.na(values) | is.nan(values)] <- 0
  as.numeric(values)
}

project_class_fraction <- function(landcover_crop, class_code, template, temp_dir) {
  tmp <- file.path(temp_dir, paste0("dlr_class_", class_code, "_fraction_", Sys.getpid(), ".tif"))
  class_mask <- terra::ifel(landcover_crop == class_code, 1, 0)
  fraction <- terra::project(
    class_mask,
    template,
    method = "average",
    threads = TRUE,
    filename = tmp,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
  fraction <- terra::clamp(fraction, lower = 0, upper = 1, values = TRUE)
  names(fraction) <- paste0("class_", class_code, "_fraction")
  fraction
}

rp_requested <- Sys.getenv("LANDUSE_RPS", unset = "rp100")
rp_names <- trimws(unlist(strsplit(rp_requested, ",")))
rp_names <- rp_names[nzchar(rp_names)]
use_cache <- toupper(Sys.getenv("LANDUSE_USE_CACHE", unset = "TRUE")) != "FALSE"

flood_rasters <- c(
  rp10 = file.path(paths$flood_raster_dir, "floodmap_EFAS_RP010_C_25832_elbe_basin.tif"),
  rp20 = file.path(paths$flood_raster_dir, "floodmap_EFAS_RP020_C_25832_elbe_basin.tif"),
  rp50 = file.path(paths$flood_raster_dir, "floodmap_EFAS_RP050_C_25832_elbe_basin.tif"),
  rp100 = file.path(paths$flood_raster_dir, "floodmap_EFAS_RP100_C_25832_elbe_basin.tif"),
  rp200 = file.path(paths$flood_raster_dir, "floodmap_EFAS_RP200_C_25832_elbe_basin.tif"),
  rp500 = file.path(paths$flood_raster_dir, "floodmap_EFAS_RP500_C_25832_elbe_basin.tif")
)

unknown_rps <- setdiff(rp_names, names(flood_rasters))
if (length(unknown_rps) > 0) {
  stop("Unknown return period(s) requested: ", paste(unknown_rps, collapse = ", "), call. = FALSE)
}

for (path_name in c("landcover_tif", "corridor_gpkg", "exposure_csv")) {
  stop_if_missing(paths[[path_name]], path_name)
}
for (rp_name in rp_names) {
  stop_if_missing(flood_rasters[[rp_name]], paste0("flood raster ", rp_name))
}

class_lookup <- tibble::tibble(
  class_code = 1:7,
  class_name = c(
    "Artificial Land",
    "Open Soil",
    "High Seasonal Vegetation",
    "High Perennial Vegetation",
    "Low Seasonal Vegetation",
    "Low Perennial Vegetation",
    "Water"
  ),
  class_slug = clean_slug(class_name),
  thesis_group = c(
    "artificial",
    "open_or_seasonal",
    "open_or_seasonal",
    "perennial_vegetation",
    "open_or_seasonal",
    "perennial_vegetation",
    "water"
  )
)

write_csv(class_lookup, file.path(paths$table_dir, "dlr_land_cover_class_lookup.csv"))

log_message("Loading corridor municipalities.")
corridor <- st_read(paths$corridor_gpkg, quiet = TRUE)
if (st_crs(corridor)$epsg != 25832) {
  corridor <- st_transform(corridor, 25832)
}

coalesce_character_columns <- function(data, cols) {
  out <- rep(NA_character_, nrow(data))
  for (col in cols) {
    values <- as.character(data[[col]])
    fill <- is.na(out) | out == ""
    out[fill] <- values[fill]
  }
  out
}

ags_candidates <- intersect(c("AGS_final", "AGS", "Gemeindeschlüssel_AGS"), names(corridor))
name_candidates <- intersect(c("mun_name_final", "mun_name", "GeografischerName_GEN"), names(corridor))
area_col <- intersect(c("municipality_area_m2", "mun_area_basin_m2"), names(corridor))[1]

if (length(ags_candidates) == 0 || length(name_candidates) == 0) {
  stop("Could not identify AGS/name columns in corridor layer.", call. = FALSE)
}

corridor$AGS <- coalesce_character_columns(corridor, ags_candidates)
corridor$mun_name <- coalesce_character_columns(corridor, name_candidates)
corridor$municipality_area_m2 <- if (!is.na(area_col)) {
  as.numeric(corridor[[area_col]])
} else {
  as.numeric(st_area(corridor))
}

corridor <- corridor %>%
  select(AGS, mun_name, municipality_area_m2)

if (anyNA(corridor$AGS)) {
  stop("At least one corridor municipality has no AGS.", call. = FALSE)
}

corridor_vect <- terra::vect(corridor)
template_full <- terra::rast(flood_rasters[[rp_names[1]]])
template_epsg <- as.character(terra::crs(template_full, describe = TRUE)$code)
if (!template_epsg %in% c("25832", "EPSG:25832")) {
  stop("Flood template is not EPSG:25832 as expected.", call. = FALSE)
}
template <- terra::crop(template_full, corridor_vect, snap = "out")
if (terra::ncell(template) == 0) {
  stop("Corridor crop of the flood template has no cells.", call. = FALSE)
}
log_message("Analysis template CRS: EPSG:", terra::crs(template, describe = TRUE)$code)
log_message("Analysis template resolution: ", paste(terra::res(template), collapse = " x "))
log_message("Analysis template dimensions: ", paste(dim(template), collapse = " x "))
log_message("Analysis template extent: ", paste(round(as.vector(terra::ext(template)), 2), collapse = ", "))

log_message("Loading DLR Land Cover DE raster.")
landcover <- terra::rast(paths$landcover_tif)
log_message("DLR raster CRS: ", terra::crs(landcover, describe = TRUE)$code)
log_message("DLR raster resolution: ", paste(terra::res(landcover), collapse = " x "))

expected_cache_paths <- unlist(lapply(class_lookup$class_slug, function(class_slug) {
  c(
    file.path(paths$cache_dir, paste0("class_total_", class_slug, ".rds")),
    file.path(paths$cache_dir, paste0("class_flood_", rp_names, "_", class_slug, ".rds"))
  )
}))
all_requested_results_cached <- use_cache && all(file.exists(expected_cache_paths))

if (all_requested_results_cached) {
  log_message("All requested class results are cached; skipping DLR crop and reprojection.")
  landcover_crop <- NULL
} else {
  template_poly <- terra::as.polygons(terra::ext(template), crs = terra::crs(template))
  template_poly_lc <- terra::project(template_poly, terra::crs(landcover))
  landcover_crop <- terra::crop(landcover, template_poly_lc, snap = "out")
  log_message("Cropped DLR raster dimensions: ", paste(dim(landcover_crop), collapse = " x "))
}

cell_area <- terra::cellSize(template, unit = "m")
names(cell_area) <- "cell_area_m2"

temp_dir <- file.path(tempdir(), "landuse_exposure_terra")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

exposure_table <- read_csv(paths$exposure_csv, show_col_types = FALSE) %>%
  mutate(AGS = as.character(AGS)) %>%
  select(AGS, starts_with("flood_area_"), starts_with("flood_share_"))

exposure_area_long <- exposure_table %>%
  select(AGS, starts_with("flood_area_")) %>%
  pivot_longer(
    cols = -AGS,
    names_to = "flood_area_metric",
    values_to = "original_flood_area_m2"
  ) %>%
  mutate(
    rp = sub("^flood_area_(rp[0-9]+)_m2$", "\\1", flood_area_metric)
  ) %>%
  select(AGS, rp, original_flood_area_m2)

log_message("Corridor municipalities: ", nrow(corridor))
log_message("Return periods requested: ", paste(rp_names, collapse = ", "))

class_total_rows <- list()
class_flood_rows <- list()

for (i in seq_len(nrow(class_lookup))) {
  class_code <- class_lookup$class_code[i]
  class_name <- class_lookup$class_name[i]
  class_slug <- class_lookup$class_slug[i]
  thesis_group <- class_lookup$thesis_group[i]
  total_cache_path <- file.path(paths$cache_dir, paste0("class_total_", class_slug, ".rds"))
  flood_cache_paths <- setNames(
    file.path(paths$cache_dir, paste0("class_flood_", rp_names, "_", class_slug, ".rds")),
    rp_names
  )

  if (use_cache && file.exists(total_cache_path) && all(file.exists(flood_cache_paths))) {
    log_message("Reading cached DLR class ", class_code, " (", class_name, ") results.")
    class_total_rows[[class_slug]] <- readRDS(total_cache_path)
    for (rp_name in rp_names) {
      class_flood_rows[[paste(rp_name, class_slug, sep = "_")]] <- readRDS(flood_cache_paths[[rp_name]])
    }
    next
  }

  log_message("Projecting DLR class ", class_code, " (", class_name, ") to EFAS grid as area fraction.")
  class_fraction <- project_class_fraction(landcover_crop, class_code, template, temp_dir)
  class_area_raster <- class_fraction * cell_area
  names(class_area_raster) <- paste0(class_slug, "_area_m2")

  total_area <- extract_sum(class_area_raster, corridor_vect)

  class_total_rows[[class_slug]] <- tibble::tibble(
    AGS = corridor$AGS,
    mun_name = corridor$mun_name,
    municipality_area_m2 = corridor$municipality_area_m2,
    class_code = class_code,
    class_name = class_name,
    class_slug = class_slug,
    thesis_group = thesis_group,
    total_area_m2 = total_area,
    total_share_of_municipality = total_area / corridor$municipality_area_m2
  )
  saveRDS(class_total_rows[[class_slug]], total_cache_path)

  for (rp_name in rp_names) {
    log_message("  Extracting flooded class area for ", rp_name, ".")
    flood_raster <- terra::crop(terra::rast(flood_rasters[[rp_name]]), template, snap = "near")
    if (!terra::compareGeom(flood_raster, template, stopOnError = FALSE)) {
      stop("Flood raster geometry does not match template for ", rp_name, call. = FALSE)
    }
    flood_presence <- terra::ifel(!is.na(flood_raster), 1, NA)
    flooded_class_area <- class_area_raster * flood_presence
    names(flooded_class_area) <- paste0(rp_name, "_", class_slug, "_flooded_m2")
    flooded_area <- extract_sum(flooded_class_area, corridor_vect)

    class_flood_rows[[paste(rp_name, class_slug, sep = "_")]] <- tibble::tibble(
      AGS = corridor$AGS,
      mun_name = corridor$mun_name,
      municipality_area_m2 = corridor$municipality_area_m2,
      rp = rp_name,
      class_code = class_code,
      class_name = class_name,
      class_slug = class_slug,
      thesis_group = thesis_group,
      total_area_m2 = total_area,
      flooded_area_m2 = flooded_area,
      flood_share_of_class = dplyr::if_else(total_area > 0, flooded_area / total_area, NA_real_)
    )
    saveRDS(class_flood_rows[[paste(rp_name, class_slug, sep = "_")]], flood_cache_paths[[rp_name]])

    rm(flood_raster, flood_presence, flooded_class_area)
    gc(verbose = FALSE)
  }

  rm(class_fraction, class_area_raster)
  gc(verbose = FALSE)
}

class_area_long <- bind_rows(class_total_rows)
class_flood_long <- bind_rows(class_flood_rows)

landcover_totals <- class_area_long %>%
  group_by(AGS) %>%
  summarise(
    landcover_total_area_m2 = sum(total_area_m2, na.rm = TRUE),
    .groups = "drop"
  )

flooded_totals <- class_flood_long %>%
  group_by(AGS, rp) %>%
  summarise(
    flooded_landcover_total_area_m2 = sum(flooded_area_m2, na.rm = TRUE),
    .groups = "drop"
  )

class_flood_long <- class_flood_long %>%
  left_join(flooded_totals, by = c("AGS", "rp")) %>%
  left_join(exposure_area_long, by = c("AGS", "rp")) %>%
  left_join(exposure_table, by = "AGS") %>%
  mutate(
    share_of_lc_flooded_area = dplyr::if_else(
      flooded_landcover_total_area_m2 > 0,
      flooded_area_m2 / flooded_landcover_total_area_m2,
      NA_real_
    ),
    share_of_original_flood_area = dplyr::if_else(
      original_flood_area_m2 > 0,
      flooded_area_m2 / original_flood_area_m2,
      NA_real_
    )
  )

group_area_long <- class_area_long %>%
  group_by(AGS, mun_name, municipality_area_m2, thesis_group) %>%
  summarise(
    total_area_m2 = sum(total_area_m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_share_of_municipality = dplyr::if_else(
      municipality_area_m2 > 0,
      total_area_m2 / municipality_area_m2,
      NA_real_
    )
  )

group_flood_long <- class_flood_long %>%
  group_by(AGS, mun_name, municipality_area_m2, rp, thesis_group) %>%
  summarise(
    total_area_m2 = sum(total_area_m2, na.rm = TRUE),
    flooded_area_m2 = sum(flooded_area_m2, na.rm = TRUE),
    flooded_landcover_total_area_m2 = first(flooded_landcover_total_area_m2),
    original_flood_area_m2 = first(original_flood_area_m2),
    flood_share_of_group = dplyr::if_else(total_area_m2 > 0, flooded_area_m2 / total_area_m2, NA_real_),
    share_of_lc_flooded_area = dplyr::if_else(
      flooded_landcover_total_area_m2 > 0,
      flooded_area_m2 / flooded_landcover_total_area_m2,
      NA_real_
    ),
    share_of_original_flood_area = dplyr::if_else(
      original_flood_area_m2 > 0,
      flooded_area_m2 / original_flood_area_m2,
      NA_real_
    ),
    .groups = "drop"
  )

qa_table <- corridor %>%
  st_drop_geometry() %>%
  left_join(landcover_totals, by = "AGS") %>%
  mutate(
    landcover_total_area_m2 = coalesce(landcover_total_area_m2, 0),
    landcover_area_ratio_to_municipality = landcover_total_area_m2 / municipality_area_m2
  ) %>%
  left_join(
    class_area_long %>%
      filter(class_slug == "artificial_land") %>%
      select(AGS, artificial_total_area_m2 = total_area_m2),
    by = "AGS"
  ) %>%
  mutate(
    artificial_total_area_m2 = coalesce(artificial_total_area_m2, 0),
    has_artificial_land = artificial_total_area_m2 > 0
  )

class_area_wide <- class_area_long %>%
  select(AGS, class_slug, total_area_m2, total_share_of_municipality) %>%
  pivot_wider(
    names_from = class_slug,
    values_from = c(total_area_m2, total_share_of_municipality),
    names_glue = "lc_{class_slug}_{.value}"
  )

group_area_wide <- group_area_long %>%
  select(AGS, thesis_group, total_area_m2, total_share_of_municipality) %>%
  pivot_wider(
    names_from = thesis_group,
    values_from = c(total_area_m2, total_share_of_municipality),
    names_glue = "lc_group_{thesis_group}_{.value}"
  )

class_flood_wide <- class_flood_long %>%
  select(
    AGS,
    rp,
    class_slug,
    flooded_area_m2,
    flood_share_of_class,
    share_of_lc_flooded_area,
    share_of_original_flood_area
  ) %>%
  pivot_wider(
    names_from = c(rp, class_slug),
    values_from = c(
      flooded_area_m2,
      flood_share_of_class,
      share_of_lc_flooded_area,
      share_of_original_flood_area
    ),
    names_glue = "{rp}_lc_class_{class_slug}_{.value}"
  )

group_flood_wide <- group_flood_long %>%
  select(
    AGS,
    rp,
    thesis_group,
    flooded_area_m2,
    flood_share_of_group,
    share_of_lc_flooded_area,
    share_of_original_flood_area
  ) %>%
  pivot_wider(
    names_from = c(rp, thesis_group),
    values_from = c(
      flooded_area_m2,
      flood_share_of_group,
      share_of_lc_flooded_area,
      share_of_original_flood_area
    ),
    names_glue = "{rp}_lc_group_{thesis_group}_{.value}"
  )

landuse_wide <- corridor %>%
  st_drop_geometry() %>%
  left_join(landcover_totals, by = "AGS") %>%
  left_join(class_area_wide, by = "AGS") %>%
  left_join(group_area_wide, by = "AGS") %>%
  left_join(class_flood_wide, by = "AGS") %>%
  left_join(group_flood_wide, by = "AGS") %>%
  left_join(flooded_totals %>% pivot_wider(names_from = rp, values_from = flooded_landcover_total_area_m2, names_prefix = "lc_flooded_total_area_m2_"), by = "AGS") %>%
  left_join(exposure_table, by = "AGS")

if ("rp100_lc_group_artificial_flooded_area_m2" %in% names(landuse_wide)) {
  landuse_wide$rp100_artificial_flooded_area_m2 <- landuse_wide$rp100_lc_group_artificial_flooded_area_m2
}
if ("rp100_lc_group_artificial_flood_share_of_group" %in% names(landuse_wide)) {
  landuse_wide$rp100_artificial_flood_share_of_group <- landuse_wide$rp100_lc_group_artificial_flood_share_of_group
}
if ("rp100_lc_group_artificial_share_of_lc_flooded_area" %in% names(landuse_wide)) {
  landuse_wide$rp100_artificial_share_of_lc_flooded_area <- landuse_wide$rp100_lc_group_artificial_share_of_lc_flooded_area
}
if ("rp100_lc_group_artificial_share_of_original_flood_area" %in% names(landuse_wide)) {
  landuse_wide$rp100_artificial_share_of_original_flood_area <- landuse_wide$rp100_lc_group_artificial_share_of_original_flood_area
}

if (file.exists(paths$analysis_csv)) {
  log_message("Joining final vulnerability analysis table for diagnostic summaries.")
  analysis <- read_csv(paths$analysis_csv, show_col_types = FALSE) %>%
    mutate(AGS = as.character(AGS)) %>%
    select(
      AGS,
      vuln_index_main_z,
      vuln_dim_access_adaptive_capacity_z,
      vuln_dim_demographic_household_z,
      vuln_dim_deprivation_labour_z,
      flood_share_rp100_original = flood_share_rp100
    )
  landuse_wide <- landuse_wide %>%
    left_join(analysis, by = "AGS")
}

summary_by_group <- group_flood_long %>%
  group_by(rp, thesis_group) %>%
  summarise(
    municipalities = n(),
    total_area_m2 = sum(total_area_m2, na.rm = TRUE),
    flooded_area_m2 = sum(flooded_area_m2, na.rm = TRUE),
    mean_flood_share_of_group = mean(flood_share_of_group, na.rm = TRUE),
    median_flood_share_of_group = median(flood_share_of_group, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(rp) %>%
  mutate(
    rp_flooded_landcover_area_m2 = sum(flooded_area_m2, na.rm = TRUE),
    share_of_all_flooded_landcover_area = dplyr::if_else(
      rp_flooded_landcover_area_m2 > 0,
      flooded_area_m2 / rp_flooded_landcover_area_m2,
      NA_real_
    )
  ) %>%
  ungroup()

summary_key_metrics <- tibble::tibble(
  metric = c(
    "corridor_municipalities",
    "municipalities_with_artificial_land",
    "municipalities_without_artificial_land",
    "mean_landcover_area_ratio_to_municipality",
    "median_landcover_area_ratio_to_municipality"
  ),
  value = c(
    nrow(corridor),
    sum(qa_table$has_artificial_land, na.rm = TRUE),
    sum(!qa_table$has_artificial_land, na.rm = TRUE),
    mean(qa_table$landcover_area_ratio_to_municipality, na.rm = TRUE),
    median(qa_table$landcover_area_ratio_to_municipality, na.rm = TRUE)
  )
)

diagnostic_correlations <- tibble::tibble()
if ("vuln_index_main_z" %in% names(landuse_wide) && "rp100_artificial_flood_share_of_group" %in% names(landuse_wide)) {
  diagnostic_correlations <- tibble::tibble(
    relationship = c(
      "vulnerability_vs_rp100_total_area_flood_share",
      "vulnerability_vs_rp100_artificial_group_flood_share",
      "access_adaptive_capacity_vs_rp100_artificial_group_flood_share",
      "demographic_household_vs_rp100_artificial_group_flood_share",
      "deprivation_labour_vs_rp100_artificial_group_flood_share"
    ),
    pearson = c(
      cor(landuse_wide$vuln_index_main_z, landuse_wide$flood_share_rp100, use = "complete.obs", method = "pearson"),
      cor(landuse_wide$vuln_index_main_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "pearson"),
      cor(landuse_wide$vuln_dim_access_adaptive_capacity_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "pearson"),
      cor(landuse_wide$vuln_dim_demographic_household_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "pearson"),
      cor(landuse_wide$vuln_dim_deprivation_labour_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "pearson")
    ),
    spearman = c(
      cor(landuse_wide$vuln_index_main_z, landuse_wide$flood_share_rp100, use = "complete.obs", method = "spearman"),
      cor(landuse_wide$vuln_index_main_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "spearman"),
      cor(landuse_wide$vuln_dim_access_adaptive_capacity_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "spearman"),
      cor(landuse_wide$vuln_dim_demographic_household_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "spearman"),
      cor(landuse_wide$vuln_dim_deprivation_labour_z, landuse_wide$rp100_artificial_flood_share_of_group, use = "complete.obs", method = "spearman")
    )
  )
}

write_csv(class_area_long, file.path(paths$table_dir, "corridor_landcover_area_by_class_long.csv"))
write_csv(class_flood_long, file.path(paths$table_dir, "corridor_flooded_landcover_by_class_long.csv"))
write_csv(group_area_long, file.path(paths$table_dir, "corridor_landcover_area_by_group_long.csv"))
write_csv(group_flood_long, file.path(paths$table_dir, "corridor_flooded_landcover_by_group_long.csv"))
write_csv(landuse_wide, file.path(paths$table_dir, "corridor_landuse_exposure_wide.csv"))
write_csv(qa_table, file.path(paths$table_dir, "corridor_landuse_exposure_qa.csv"))
write_csv(summary_by_group, file.path(paths$table_dir, "corridor_landuse_exposure_summary_by_group.csv"))
write_csv(summary_key_metrics, file.path(paths$table_dir, "corridor_landuse_exposure_key_metrics.csv"))
write_csv(diagnostic_correlations, file.path(paths$table_dir, "corridor_landuse_exposure_diagnostic_correlations.csv"))

landuse_gpkg <- corridor %>%
  left_join(landuse_wide %>% select(-mun_name, -municipality_area_m2), by = "AGS")

st_write(
  landuse_gpkg,
  file.path(paths$gpkg_dir, "corridor_landuse_exposure.gpkg"),
  layer = "corridor_landuse_exposure",
  delete_dsn = TRUE,
  quiet = TRUE
)

if ("rp100_artificial_flood_share_of_group" %in% names(landuse_wide)) {
  plot_tbl <- landuse_wide %>%
    filter(!is.na(vuln_index_main_z), !is.na(rp100_artificial_flood_share_of_group))

  p_scatter <- ggplot(plot_tbl, aes(vuln_index_main_z, rp100_artificial_flood_share_of_group)) +
    geom_point(alpha = 0.65, size = 1.7, color = "#2f6f73") +
    geom_smooth(method = "lm", se = FALSE, color = "#202020", linewidth = 0.5) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "RP100 artificial-land flood share and socio-economic vulnerability",
      subtitle = "DLR Land Cover DE class 1, aggregated to corridor municipalities",
      x = "Socio-economic vulnerability index (z-score)",
      y = "Flooded share of artificial land",
      caption = "Data: DLR Land Cover DE 2015, EFAS/JRC RP100 flood raster, INKAR. Own processing."
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    file.path(paths$plot_dir, "plot_vulnerability_vs_rp100_artificial_land_flood_share.png"),
    p_scatter,
    width = 8.5,
    height = 6,
    dpi = 320
  )

  map_data <- landuse_gpkg
  p_map <- ggplot(map_data) +
    geom_sf(aes(fill = rp100_artificial_flood_share_of_group), color = "white", linewidth = 0.03) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      labels = percent_format(accuracy = 1),
      na.value = "#efefef",
      name = "Flooded artificial land"
    ) +
    labs(
      title = "RP100 artificial-land exposure in the Elbe corridor",
      subtitle = "Share of DLR artificial land class intersecting modeled RP100 flooding",
      caption = "Data: DLR Land Cover DE 2015, EFAS/JRC RP100 flood raster. Own processing."
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )

  ggsave(
    file.path(paths$plot_dir, "map_rp100_artificial_land_flood_share.png"),
    p_map,
    width = 10.5,
    height = 8,
    dpi = 320
  )
}

overview <- c(
  "# DLR Land-Use Exposure Module",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Processing summary",
  "",
  "- DLR Land Cover DE is read as a 10 m categorical raster in EPSG:3035.",
  "- The analysis template is the corridor-cropped EFAS/JRC flood grid in EPSG:25832.",
  paste0("- Analysis template resolution: `", paste(terra::res(template), collapse = " x "), " m`."),
  paste0("- Analysis template dimensions: `", paste(dim(template), collapse = " x "), "`."),
  paste0("- Analysis template extent: `", paste(round(as.vector(terra::ext(template)), 2), collapse = ", "), "`."),
  "- Each DLR class is converted to a binary raster and aggregated to the 100 m EFAS grid with `terra::project(..., method = \"average\")`.",
  "- Class fractions are multiplied by the EPSG:25832 cell area and summed within corridor municipality polygons.",
  "- Flooded class areas use the same class-area rasters masked by valid RP100 flood-depth cells.",
  "",
  "## Inputs",
  "",
  paste0("- DLR Land Cover DE raster: `", paths$landcover_tif, "`"),
  paste0("- Corridor municipalities: `", paths$corridor_gpkg, "`"),
  paste0("- Exposure CSV: `", paths$exposure_csv, "`"),
  paste0("- Flood raster directory: `", paths$flood_raster_dir, "`"),
  paste0("- Return periods processed: `", paste(rp_names, collapse = ", "), "`"),
  "",
  "## DLR classes",
  "",
  paste0("- `", class_lookup$class_code, "`: ", class_lookup$class_name, " -> `", class_lookup$thesis_group, "`"),
  "",
  "## Main outputs",
  "",
  "- `outputs/tables/corridor_landcover_area_by_class_long.csv`",
  "- `outputs/tables/corridor_flooded_landcover_by_class_long.csv`",
  "- `outputs/tables/corridor_landcover_area_by_group_long.csv`",
  "- `outputs/tables/corridor_flooded_landcover_by_group_long.csv`",
  "- `outputs/tables/corridor_landuse_exposure_wide.csv`",
  "- `outputs/tables/corridor_landuse_exposure_diagnostic_correlations.csv`",
  "- `outputs/gpkg/corridor_landuse_exposure.gpkg`",
  "",
  "## RP100 key results",
  "",
  paste0(
    "- `", summary_by_group$thesis_group, "`: ",
    round(summary_by_group$flooded_area_m2 / 1e6, 1), " km² flooded; ",
    round(summary_by_group$share_of_all_flooded_landcover_area * 100, 1),
    "% of all flooded land-cover area; median group flood share ",
    round(summary_by_group$median_flood_share_of_group * 100, 1), "%."
  ),
  "",
  "## Diagnostic correlations",
  "",
  if (nrow(diagnostic_correlations) > 0) {
    paste0(
      "- `", diagnostic_correlations$relationship, "`: Pearson ",
      round(diagnostic_correlations$pearson, 3), ", Spearman ",
      round(diagnostic_correlations$spearman, 3), "."
    )
  } else {
    "- Diagnostic correlations were not calculated because the vulnerability table or RP100 artificial-land metric was unavailable."
  },
  "",
  "## Interpretation note",
  "",
  "The `artificial` group is based on DLR class 1, `Artificial Land`. It is used as a built/artificial land-cover exposure proxy, not as a cadastral building-footprint or residential-population exposure measure.",
  "The `open_or_seasonal` and `perennial_vegetation` groups are analytical land-cover proxies and should not be overinterpreted as exact agricultural or natural land-use classes."
)

writeLines(overview, file.path(paths$output_root, "LANDUSE_EXPOSURE_OVERVIEW.md"))

log_message("Land-use exposure module complete.")
log_message("Wide table: ", file.path(paths$table_dir, "corridor_landuse_exposure_wide.csv"))
log_message("GPKG: ", file.path(paths$gpkg_dir, "corridor_landuse_exposure.gpkg"))
