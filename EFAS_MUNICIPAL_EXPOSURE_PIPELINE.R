#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(readr)
})

options(scipen = 999)

# ============================================================
# Municipal flood-exposure pipeline across EFAS return periods
# ------------------------------------------------------------
# Logic
# - corridor is defined empirically by RP500 flood-depth cells
# - flooded area is derived from valid depth pixels (NoData = not flooded)
# - no raster polygonization is required for the main workflow
# - exposure is calculated for RP10, RP20, RP50, RP100, RP200, RP500
# ============================================================

# ---------------------------
# 1) User-editable paths
# ---------------------------

municipality_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/gemeinden_elbe_landonly_basin_inkar.gpkg"
municipality_layer <- "gemeinden_elbe_landonly_basin_inkar"

rp_files <- c(
  rp10 = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_eu_flood_25832/floodmap_EFAS_RP010_C_25832_elbe_basin.tif",
  rp20 = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_eu_flood_25832/floodmap_EFAS_RP020_C_25832_elbe_basin.tif",
  rp50 = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_eu_flood_25832/floodmap_EFAS_RP050_C_25832_elbe_basin.tif",
  rp100 = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_eu_flood_25832/floodmap_EFAS_RP100_C_25832_elbe_basin.tif",
  rp200 = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_eu_flood_25832/floodmap_EFAS_RP200_C_25832_elbe_basin.tif",
  rp500 = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_eu_flood_25832/floodmap_EFAS_RP500_C_25832_elbe_basin.tif"
)

base_output_dir <- Sys.getenv(
  "EXPOSURE_PIPELINE_BASE_OUT",
  unset = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2"
)

# Optional:
# polygonizing the RP500 extent is disabled by default because the workflow
# should remain raster-based and polygonization can become expensive.
export_rp500_extent_vector <- FALSE
use_exact_extract <- !identical(Sys.getenv("EXPOSURE_PIPELINE_EXACT", unset = "TRUE"), "FALSE")

# ---------------------------
# 2) Output structure
# ---------------------------

output_root <- file.path(base_output_dir, "outputs_exposure_pipeline")
corridor_dir <- file.path(output_root, "corridor")
tables_dir <- file.path(output_root, "tables")
gpkg_dir <- file.path(output_root, "gpkg")
logs_dir <- file.path(output_root, "logs")

for (dir_path in c(output_root, corridor_dir, tables_dir, gpkg_dir, logs_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(logs_dir, "processing_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

save_csv <- function(x, filename) {
  write_csv(x, file.path(tables_dir, filename))
}

write_gpkg_layer <- function(x, dsn, layer, replace_dsn = FALSE) {
  if (replace_dsn && file.exists(dsn)) {
    file.remove(dsn)
  }

  if (!file.exists(dsn)) {
    write_sf(x, dsn = dsn, layer = layer, quiet = TRUE)
  } else {
    write_sf(x, dsn = dsn, layer = layer, append = TRUE, quiet = TRUE)
  }
}

pick_column <- function(candidates, available, label) {
  hit <- intersect(candidates, available)
  if (length(hit) == 0) {
    stop("Could not find a column for ", label, ". Checked: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

extract_flooded_area <- function(depth_raster, municipalities_vect, cell_area_raster, exact = TRUE) {
  target_extent <- ext(municipalities_vect)
  depth_crop <- crop(depth_raster, target_extent, snap = "out")
  area_crop <- crop(cell_area_raster, target_extent, snap = "out")
  flooded_cell_area <- mask(area_crop, depth_crop)
  extraction <- extract(
    flooded_cell_area,
    municipalities_vect,
    fun = sum,
    na.rm = TRUE,
    exact = exact
  )
  flooded_area <- extraction[[2]]
  flooded_area[is.na(flooded_area) | is.nan(flooded_area)] <- 0
  as.numeric(flooded_area)
}

# ---------------------------
# 3) Load municipalities
# ---------------------------

log_message("Loading municipality layer ...")

municipalities_sf <- st_read(municipality_path, layer = municipality_layer, quiet = TRUE)

if (is.na(st_crs(municipalities_sf)$epsg) || st_crs(municipalities_sf)$epsg != 25832) {
  log_message("Municipality layer is not in EPSG:25832. Reprojecting now.")
  municipalities_sf <- st_transform(municipalities_sf, 25832)
}

ags_col <- pick_column(
  c("AGS", "Gemeindeschlüssel_AGS", "GemeindeschlüsselAufgefüllt"),
  names(municipalities_sf),
  "municipality AGS"
)

name_col <- pick_column(
  c("mun_name", "Gemeinde", "GeografischerName_GEN", "Bezeichnung"),
  names(municipalities_sf),
  "municipality name"
)

municipality_area_values <- as.numeric(expanse(vect(municipalities_sf), unit = "m"))

municipalities_sf <- municipalities_sf %>%
  mutate(
    zone_id = seq_len(n()),
    AGS_final = as.character(.data[[ags_col]]),
    mun_name_final = as.character(.data[[name_col]]),
    municipality_area_m2 = municipality_area_values
  )

if ("mun_area_basin_m2" %in% names(municipalities_sf)) {
  area_diff <- max(abs(municipalities_sf$municipality_area_m2 - municipalities_sf$mun_area_basin_m2), na.rm = TRUE)
  log_message("Max difference between geometry area and existing mun_area_basin_m2: ", round(area_diff, 3), " m2")
}

log_message("Municipalities loaded: ", nrow(municipalities_sf))
log_message("AGS column used: ", ags_col)
log_message("Name column used: ", name_col)
log_message("Exact extraction enabled: ", use_exact_extract)

municipalities_vect <- vect(municipalities_sf)

# ---------------------------
# 4) Load and validate rasters
# ---------------------------

log_message("Loading flood rasters ...")

for (file_path in rp_files) {
  if (!file.exists(file_path)) stop("Missing raster: ", file_path)
}

rasters <- lapply(rp_files, rast)

template_raster <- rasters[[1]]
template_crs <- crs(template_raster, proj = TRUE)
template_extent <- ext(template_raster)
template_resolution <- res(template_raster)

for (rp_name in names(rasters)) {
  r <- rasters[[rp_name]]
  if (crs(r, proj = TRUE) != template_crs) {
    stop("Raster CRS mismatch for ", rp_name)
  }
  if (!all.equal(as.vector(ext(r)), as.vector(template_extent))) {
    stop("Raster extent mismatch for ", rp_name)
  }
  if (!all.equal(as.numeric(res(r)), as.numeric(template_resolution))) {
    stop("Raster resolution mismatch for ", rp_name)
  }
}

log_message("All rasters share the same CRS, extent, and resolution.")
log_message("Raster resolution: ", paste(template_resolution, collapse = " x "), " m")

cell_area_raster <- cellSize(template_raster, unit = "m")

# ---------------------------
# 5) Define RP500 corridor
# ---------------------------

log_message("Deriving municipality corridor from RP500 valid flood-depth pixels ...")
log_message("Step 1/2: fast raster preselection of potentially intersecting municipalities ...")

zone_raster <- rasterize(
  vect(municipalities_sf[, "zone_id"]),
  template_raster,
  field = "zone_id",
  touches = TRUE
)

rp500_presence <- ifel(!is.na(rasters[["rp500"]]), 1, NA)
rp500_zonal <- as.data.frame(zonal(rp500_presence, zone_raster, fun = "sum", na.rm = TRUE))
names(rp500_zonal) <- c("zone_id", "rp500_flooded_cells")
rp500_zonal$rp500_flooded_cells[is.na(rp500_zonal$rp500_flooded_cells)] <- 0

municipalities_sf <- municipalities_sf %>%
  left_join(rp500_zonal, by = "zone_id") %>%
  mutate(rp500_flooded_cells = if_else(is.na(rp500_flooded_cells), 0, rp500_flooded_cells))

pre_corridor_sf <- municipalities_sf %>%
  filter(rp500_flooded_cells > 0)

if (nrow(pre_corridor_sf) == 0) {
  stop("No municipalities intersect the RP500 flood extent in the raster preselection step.")
}

log_message("Potential corridor municipalities after raster preselection: ", nrow(pre_corridor_sf))
  log_message("Step 2/2: exact flooded-area extract for the RP500 preselection ...")
  if (!use_exact_extract) {
    log_message("Exact extraction is disabled for this run; using approximate raster extraction instead.")
  }

pre_corridor_vect <- vect(pre_corridor_sf)
rp500_area_exact <- extract_flooded_area(
  depth_raster = rasters[["rp500"]],
  municipalities_vect = pre_corridor_vect,
  cell_area_raster = cell_area_raster,
  exact = use_exact_extract
)

pre_corridor_sf$rp500_flood_area_m2_exact <- rp500_area_exact

corridor_sf <- pre_corridor_sf %>%
  filter(rp500_flood_area_m2_exact > 0) %>%
  select(-rp500_flooded_cells)

if (nrow(corridor_sf) == 0) {
  stop("No municipalities intersect the RP500 flood extent.")
}

log_message("Corridor municipalities retained: ", nrow(corridor_sf))

corridor_summary <- tibble(
  municipalities_total = nrow(municipalities_sf),
  municipalities_in_corridor = nrow(corridor_sf),
  municipalities_outside_corridor = nrow(municipalities_sf) - nrow(corridor_sf)
)

save_csv(corridor_summary, "corridor_summary.csv")

corridor_gpkg <- file.path(corridor_dir, "municipalities_corridor.gpkg")
write_gpkg_layer(corridor_sf, corridor_gpkg, layer = "municipalities_corridor", replace_dsn = TRUE)

if (isTRUE(export_rp500_extent_vector)) {
  log_message("Exporting RP500 flood extent polygon (optional heavy step) ...")
  rp500_mask <- ifel(!is.na(rasters[["rp500"]]), 1, NA)
  rp500_poly <- as.polygons(rp500_mask, dissolve = TRUE, na.rm = TRUE)
  writeVector(rp500_poly, file.path(corridor_dir, "RP500_flood_extent.gpkg"), overwrite = TRUE)
}

corridor_vect <- vect(corridor_sf)

# ---------------------------
# 6) Exposure per return period
# ---------------------------

log_message("Calculating municipal exposure for all return periods ...")

rp_results <- list()
rp_gpkg <- file.path(gpkg_dir, "municipalities_corridor_exposure_by_rp.gpkg")
if (file.exists(rp_gpkg)) file.remove(rp_gpkg)

for (i in seq_along(rasters)) {
  rp_name <- names(rasters)[i]
  rp_raster <- rasters[[i]]

  log_message("Processing ", rp_name, " ...")

  flooded_area <- extract_flooded_area(
    depth_raster = rp_raster,
    municipalities_vect = corridor_vect,
      cell_area_raster = cell_area_raster,
      exact = use_exact_extract
  )

  flooded_share <- flooded_area / corridor_sf$municipality_area_m2

  rp_table <- tibble(
    AGS = corridor_sf$AGS_final,
    mun_name = corridor_sf$mun_name_final,
    municipality_area_m2 = corridor_sf$municipality_area_m2,
    !!paste0("flood_area_", rp_name, "_m2") := flooded_area,
    !!paste0("flood_share_", rp_name) := flooded_share
  )

  rp_results[[rp_name]] <- rp_table

  save_csv(rp_table, paste0("municipality_flood_exposure_", rp_name, ".csv"))

  rp_layer_sf <- corridor_sf %>%
    select(AGS_final, mun_name_final, municipality_area_m2) %>%
    rename(AGS = AGS_final, mun_name = mun_name_final) %>%
    mutate(
      !!paste0("flood_area_", rp_name, "_m2") := flooded_area,
      !!paste0("flood_share_", rp_name) := flooded_share
    )

  write_gpkg_layer(
    rp_layer_sf,
    dsn = rp_gpkg,
    layer = paste0("exposure_", rp_name)
  )
}

# ---------------------------
# 7) Merge final wide dataset
# ---------------------------

log_message("Merging wide municipality exposure table ...")

final_table <- corridor_sf %>%
  st_drop_geometry() %>%
  transmute(
    AGS = AGS_final,
    mun_name = mun_name_final,
    municipality_area_m2 = municipality_area_m2
  )

for (rp_name in names(rp_results)) {
  final_table <- final_table %>%
    left_join(
      rp_results[[rp_name]] %>% select(-mun_name, -municipality_area_m2),
      by = "AGS"
    )
}

final_gpkg_sf <- corridor_sf %>%
  select(AGS_final, mun_name_final, municipality_area_m2) %>%
  rename(AGS = AGS_final, mun_name = mun_name_final)

for (rp_name in names(rp_results)) {
  rp_cols <- rp_results[[rp_name]] %>%
    select(AGS, starts_with(paste0("flood_area_", rp_name)), starts_with(paste0("flood_share_", rp_name)))
  final_gpkg_sf <- final_gpkg_sf %>% left_join(rp_cols, by = "AGS")
}

save_csv(final_table, "municipality_flood_exposure_all_RPs.csv")

final_gpkg <- file.path(gpkg_dir, "municipalities_corridor_exposure_all_RPs.gpkg")
write_gpkg_layer(final_gpkg_sf, final_gpkg, layer = "municipalities_corridor_exposure_all_rps", replace_dsn = TRUE)

# ---------------------------
# 8) Quality checks
# ---------------------------

log_message("Running quality checks ...")

share_cols <- paste0("flood_share_", names(rp_results))
area_cols <- paste0("flood_area_", names(rp_results), "_m2")
tol <- 1e-9

quality_checks <- final_table

for (rp_name in names(rp_results)) {
  share_col <- paste0("flood_share_", rp_name)
  area_col <- paste0("flood_area_", rp_name, "_m2")

  quality_checks[[paste0("share_bounds_ok_", rp_name)]] <-
    !is.na(quality_checks[[share_col]]) &
    quality_checks[[share_col]] >= -tol &
    quality_checks[[share_col]] <= 1 + tol

  quality_checks[[paste0("area_bounds_ok_", rp_name)]] <-
    !is.na(quality_checks[[area_col]]) &
    quality_checks[[area_col]] >= -tol &
    quality_checks[[area_col]] <= quality_checks$municipality_area_m2 + tol
}

monotonic_pairs <- list(
  c("rp10", "rp20"),
  c("rp20", "rp50"),
  c("rp50", "rp100"),
  c("rp100", "rp200"),
  c("rp200", "rp500")
)

for (pair in monotonic_pairs) {
  left <- paste0("flood_share_", pair[1])
  right <- paste0("flood_share_", pair[2])
  check_name <- paste0("mono_", pair[1], "_le_", pair[2])
  diff_name <- paste0("delta_", pair[2], "_minus_", pair[1])

  quality_checks[[diff_name]] <- quality_checks[[right]] - quality_checks[[left]]
  quality_checks[[check_name]] <- quality_checks[[diff_name]] >= -tol
}

mono_check_cols <- grep("^mono_", names(quality_checks), value = TRUE)
share_check_cols <- grep("^share_bounds_ok_", names(quality_checks), value = TRUE)
area_check_cols <- grep("^area_bounds_ok_", names(quality_checks), value = TRUE)

quality_checks <- quality_checks %>%
  mutate(
    share_bounds_all_ok = if_all(all_of(share_check_cols), identity),
    area_bounds_all_ok = if_all(all_of(area_check_cols), identity),
    monotonicity_all_ok = if_all(all_of(mono_check_cols), identity),
    suspicious_any = !(share_bounds_all_ok & area_bounds_all_ok & monotonicity_all_ok)
  )

save_csv(quality_checks, "exposure_quality_checks.csv")

suspicious_municipalities <- quality_checks %>%
  filter(suspicious_any)

save_csv(suspicious_municipalities, "exposure_quality_checks_suspicious_only.csv")

log_message("Municipalities flagged as suspicious: ", nrow(suspicious_municipalities))

summary_checks <- tibble(
  municipalities_in_corridor = nrow(final_table),
  suspicious_municipalities = nrow(suspicious_municipalities),
  all_share_bounds_ok = all(quality_checks$share_bounds_all_ok),
  all_area_bounds_ok = all(quality_checks$area_bounds_all_ok),
  all_monotonicity_ok = all(quality_checks$monotonicity_all_ok)
)

save_csv(summary_checks, "quality_check_summary.csv")

# ---------------------------
# 9) Final notes
# ---------------------------

log_message("Pipeline finished successfully.")
log_message("Final wide CSV: ", file.path(tables_dir, "municipality_flood_exposure_all_RPs.csv"))
log_message("Final exposure GPKG: ", final_gpkg)
log_message("Quality checks CSV: ", file.path(tables_dir, "exposure_quality_checks.csv"))
