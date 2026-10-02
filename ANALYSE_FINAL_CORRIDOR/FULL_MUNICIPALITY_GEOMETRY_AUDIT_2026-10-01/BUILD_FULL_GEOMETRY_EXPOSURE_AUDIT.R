#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(scipen = 999)

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"
external_root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS"

audit_root <- file.path(
  root,
  "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_GEOMETRY_AUDIT_2026-10-01"
)
table_dir <- file.path(audit_root, "outputs/tables")
gpkg_dir <- file.path(audit_root, "outputs/gpkg")
log_dir <- file.path(audit_root, "outputs/logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gpkg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  old_corridor = file.path(
    external_root,
    "analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg"
  ),
  old_exposure = file.path(
    external_root,
    "analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv"
  ),
  full_municipalities = file.path(
    external_root,
    "ANALYSIS.nosync/vg250_gemeinden_landonly.gpkg"
  ),
  current_analysis = file.path(
    root,
    "ANALYSE_FINAL_CORRIDOR/outputs/tables/corridor_analysis_rp100.csv"
  )
)

rp_files <- c(
  rp10 = file.path(root, "outputs_eu_flood_25832/floodmap_EFAS_RP010_C_25832_elbe_basin.tif"),
  rp20 = file.path(root, "outputs_eu_flood_25832/floodmap_EFAS_RP020_C_25832_elbe_basin.tif"),
  rp50 = file.path(root, "outputs_eu_flood_25832/floodmap_EFAS_RP050_C_25832_elbe_basin.tif"),
  rp100 = file.path(root, "outputs_eu_flood_25832/floodmap_EFAS_RP100_C_25832_elbe_basin.tif"),
  rp200 = file.path(root, "outputs_eu_flood_25832/floodmap_EFAS_RP200_C_25832_elbe_basin.tif"),
  rp500 = file.path(root, "outputs_eu_flood_25832/floodmap_EFAS_RP500_C_25832_elbe_basin.tif")
)

required_files <- c(unlist(paths), rp_files)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
}

log_file <- file.path(log_dir, "audit_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_message <- function(...) {
  message_text <- paste(..., collapse = "")
  cat(message_text, "\n")
  cat(message_text, "\n", file = log_file, append = TRUE)
}

safe_cor <- function(x, y, method) {
  suppressWarnings(cor(x, y, use = "complete.obs", method = method))
}

extract_flooded_area <- function(depth_raster, municipalities, cell_area_raster) {
  flooded_cell_area <- terra::mask(cell_area_raster, depth_raster)
  extracted <- terra::extract(
    flooded_cell_area,
    terra::vect(municipalities),
    fun = sum,
    na.rm = TRUE,
    exact = TRUE
  )[[2]]
  extracted[is.na(extracted) | is.nan(extracted)] <- 0
  as.numeric(extracted)
}

log_message("Loading existing basin-clipped corridor geometry.")
old_corridor <- st_read(paths$old_corridor, quiet = TRUE) %>%
  st_transform(25832) %>%
  mutate(
    AGS = coalesce(as.character(AGS), as.character(Gemeindeschlüssel_AGS)),
    mun_name = coalesce(as.character(mun_name), as.character(GeografischerName_GEN)),
    old_geometry_area_m2 = as.numeric(st_area(.))
  )

if (anyNA(old_corridor$AGS) || anyDuplicated(old_corridor$AGS)) {
  stop("Existing corridor AGS values are missing or duplicated after repair.", call. = FALSE)
}

log_message("Loading complete BKG VG250 municipality geometry.")
full_municipalities <- st_read(paths$full_municipalities, quiet = TRUE) %>%
  st_transform(25832) %>%
  transmute(
    AGS = as.character(Gemeindeschlüssel_AGS),
    mun_name_full = as.character(GeografischerName_GEN)
  )

full_corridor <- full_municipalities %>%
  filter(AGS %in% old_corridor$AGS) %>%
  arrange(match(AGS, old_corridor$AGS)) %>%
  mutate(municipality_area_m2 = as.numeric(st_area(.)))

if (nrow(full_corridor) != nrow(old_corridor)) {
  stop(
    "Expected ", nrow(old_corridor), " full municipality geometries but found ",
    nrow(full_corridor), ".",
    call. = FALSE
  )
}

if (!identical(full_corridor$AGS, old_corridor$AGS)) {
  stop("Full municipality geometries do not align with the existing corridor AGS order.", call. = FALSE)
}

st_write(
  full_corridor,
  file.path(gpkg_dir, "full_municipality_corridor_geometry.gpkg"),
  layer = "full_municipality_corridor_geometry",
  delete_dsn = TRUE,
  quiet = TRUE
)

geometry_comparison <- old_corridor %>%
  st_drop_geometry() %>%
  select(AGS, mun_name, old_geometry_area_m2) %>%
  left_join(
    full_corridor %>%
      st_drop_geometry() %>%
      select(AGS, mun_name_full, full_geometry_area_m2 = municipality_area_m2),
    by = "AGS"
  ) %>%
  mutate(
    old_to_full_area_ratio = old_geometry_area_m2 / full_geometry_area_m2,
    omitted_area_m2 = full_geometry_area_m2 - old_geometry_area_m2,
    geometry_was_clipped = old_to_full_area_ratio < 0.999
  ) %>%
  arrange(old_to_full_area_ratio)

write_csv(geometry_comparison, file.path(table_dir, "municipality_geometry_area_comparison.csv"))

log_message(
  "Municipalities with old/full area ratio below 0.999: ",
  sum(geometry_comparison$geometry_was_clipped)
)

log_message("Recalculating exposure with complete municipality geometries.")
rasters <- lapply(rp_files, terra::rast)
template <- rasters[[1]]
cell_area <- terra::cellSize(template, unit = "m")

exposure <- full_corridor %>%
  st_drop_geometry() %>%
  select(AGS, mun_name = mun_name_full, municipality_area_m2)

for (rp_name in names(rasters)) {
  log_message("Extracting ", toupper(rp_name), ".")
  flooded_area <- extract_flooded_area(rasters[[rp_name]], full_corridor, cell_area)
  exposure[[paste0("flood_area_", rp_name, "_m2")]] <- flooded_area
  exposure[[paste0("flood_share_", rp_name)]] <- flooded_area / exposure$municipality_area_m2
}

write_csv(exposure, file.path(table_dir, "full_municipality_exposure_all_RPs.csv"))

full_corridor_output <- full_corridor %>%
  left_join(exposure %>% select(-mun_name, -municipality_area_m2), by = "AGS")

st_write(
  full_corridor_output,
  file.path(gpkg_dir, "full_municipality_corridor_exposure.gpkg"),
  layer = "full_municipality_corridor_exposure",
  delete_dsn = TRUE,
  quiet = TRUE
)

old_exposure <- read_csv(
  paths$old_exposure,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

if (sum(is.na(old_exposure$AGS)) == 1) {
  old_exposure$AGS[is.na(old_exposure$AGS)] <- "16076094"
}

comparison <- old_exposure %>%
  select(
    AGS,
    old_municipality_area_m2 = municipality_area_m2,
    old_flood_area_rp100_m2 = flood_area_rp100_m2,
    old_flood_share_rp100 = flood_share_rp100,
    old_flood_share_rp500 = flood_share_rp500
  ) %>%
  left_join(
    exposure %>%
      select(
        AGS,
        new_municipality_area_m2 = municipality_area_m2,
        new_flood_area_rp100_m2 = flood_area_rp100_m2,
        new_flood_share_rp100 = flood_share_rp100,
        new_flood_share_rp500 = flood_share_rp500
      ),
    by = "AGS"
  ) %>%
  mutate(
    area_ratio_old_to_new = old_municipality_area_m2 / new_municipality_area_m2,
    rp100_share_change_pp = 100 * (new_flood_share_rp100 - old_flood_share_rp100),
    rp500_share_change_pp = 100 * (new_flood_share_rp500 - old_flood_share_rp500),
    rp100_flood_area_change_m2 = new_flood_area_rp100_m2 - old_flood_area_rp100_m2
  ) %>%
  arrange(rp100_share_change_pp)

write_csv(comparison, file.path(table_dir, "old_vs_full_geometry_exposure_comparison.csv"))

analysis <- read_csv(
  paths$current_analysis,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

analysis_comparison <- analysis %>%
  select(AGS, vuln_index_main_z) %>%
  left_join(
    comparison %>%
      select(AGS, old_flood_share_rp100, new_flood_share_rp100),
    by = "AGS"
  )

impact_summary <- tibble(
  metric = c(
    "municipalities",
    "municipalities_old_area_ratio_below_0.999",
    "municipalities_old_area_ratio_below_0.95",
    "old_rp100_mean",
    "new_rp100_mean",
    "old_rp100_median",
    "new_rp100_median",
    "old_vulnerability_rp100_pearson",
    "new_vulnerability_rp100_pearson",
    "old_vulnerability_rp100_spearman",
    "new_vulnerability_rp100_spearman",
    "maximum_absolute_rp100_share_change_pp"
  ),
  value = c(
    nrow(comparison),
    sum(geometry_comparison$old_to_full_area_ratio < 0.999),
    sum(geometry_comparison$old_to_full_area_ratio < 0.95),
    mean(comparison$old_flood_share_rp100),
    mean(comparison$new_flood_share_rp100),
    median(comparison$old_flood_share_rp100),
    median(comparison$new_flood_share_rp100),
    safe_cor(
      analysis_comparison$vuln_index_main_z,
      analysis_comparison$old_flood_share_rp100,
      "pearson"
    ),
    safe_cor(
      analysis_comparison$vuln_index_main_z,
      analysis_comparison$new_flood_share_rp100,
      "pearson"
    ),
    safe_cor(
      analysis_comparison$vuln_index_main_z,
      analysis_comparison$old_flood_share_rp100,
      "spearman"
    ),
    safe_cor(
      analysis_comparison$vuln_index_main_z,
      analysis_comparison$new_flood_share_rp100,
      "spearman"
    ),
    max(abs(comparison$rp100_share_change_pp))
  )
)

write_csv(impact_summary, file.path(table_dir, "full_geometry_impact_summary.csv"))

top_changes <- comparison %>%
  mutate(abs_rp100_share_change_pp = abs(rp100_share_change_pp)) %>%
  arrange(desc(abs_rp100_share_change_pp)) %>%
  slice_head(n = 50)

write_csv(top_changes, file.path(table_dir, "largest_rp100_share_changes.csv"))

log_message("Audit completed without modifying existing outputs.")
print(impact_summary, n = Inf)
