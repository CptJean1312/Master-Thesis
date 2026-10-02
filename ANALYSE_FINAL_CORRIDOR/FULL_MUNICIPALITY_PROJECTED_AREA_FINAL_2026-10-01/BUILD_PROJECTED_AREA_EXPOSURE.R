#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(readr)
})

options(scipen = 999)

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"
output_root <- file.path(
  root,
  "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01"
)
table_dir <- file.path(output_root, "outputs/tables")
gpkg_dir <- file.path(output_root, "outputs/gpkg")
log_dir <- file.path(output_root, "outputs/logs")

for (path in c(table_dir, gpkg_dir, log_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

paths <- list(
  full_corridor = file.path(
    root,
    "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_GEOMETRY_AUDIT_2026-10-01/outputs/gpkg/full_municipality_corridor_geometry.gpkg"
  ),
  previous_exposure = file.path(
    root,
    "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_GEOMETRY_AUDIT_2026-10-01/outputs/tables/full_municipality_exposure_all_RPs.csv"
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

required <- c(unlist(paths), rp_files)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required inputs:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

log_file <- file.path(log_dir, "build_projected_area_exposure.log")
if (file.exists(log_file)) file.remove(log_file)

log_message <- function(...) {
  message_text <- paste(..., collapse = "")
  cat(message_text, "\n")
  cat(message_text, "\n", file = log_file, append = TRUE)
}

extract_projected_flood_area <- function(depth_raster, municipalities) {
  projected_cell_area_m2 <- prod(terra::res(depth_raster))
  flooded_area_raster <- terra::ifel(
    !is.na(depth_raster),
    projected_cell_area_m2,
    NA
  )
  extracted <- terra::extract(
    flooded_area_raster,
    terra::vect(municipalities),
    fun = sum,
    na.rm = TRUE,
    exact = TRUE
  )[[2]]
  extracted[is.na(extracted) | is.nan(extracted)] <- 0
  as.numeric(extracted)
}

log_message("Loading complete corridor municipality geometries.")
corridor <- st_read(paths$full_corridor, quiet = TRUE) %>%
  st_transform(25832) %>%
  mutate(
    AGS = as.character(AGS),
    mun_name = as.character(mun_name_full),
    municipality_area_m2 = as.numeric(st_area(.))
  ) %>%
  select(AGS, mun_name, municipality_area_m2)

if (nrow(corridor) != 835 || anyNA(corridor$AGS) || anyDuplicated(corridor$AGS)) {
  stop("Corridor geometry failed the row or AGS checks.", call. = FALSE)
}

exposure <- corridor %>%
  st_drop_geometry() %>%
  select(AGS, mun_name, municipality_area_m2)

log_message("Calculating flooded area with projected 100 m by 100 m cells.")
for (rp_name in names(rp_files)) {
  depth <- terra::rast(rp_files[[rp_name]])
  if (!identical(as.numeric(terra::res(depth)), c(100, 100))) {
    stop("Unexpected raster resolution for ", rp_name, call. = FALSE)
  }
  flooded_area <- extract_projected_flood_area(depth, corridor)
  area_col <- paste0("flood_area_", rp_name, "_m2")
  share_col <- paste0("flood_share_", rp_name)
  exposure[[area_col]] <- flooded_area
  exposure[[share_col]] <- flooded_area / exposure$municipality_area_m2
  log_message("Finished ", toupper(rp_name), ".")
}

write_csv(
  exposure,
  file.path(table_dir, "full_municipality_exposure_all_RPs_projected_area.csv")
)

corridor_output <- corridor %>%
  left_join(exposure %>% select(-mun_name, -municipality_area_m2), by = "AGS")

st_write(
  corridor_output,
  file.path(gpkg_dir, "full_municipality_corridor_exposure_projected_area.gpkg"),
  layer = "full_municipality_exposure",
  delete_dsn = TRUE,
  quiet = TRUE
)

rp_names <- names(rp_files)
tol <- 1e-9
quality <- exposure

for (rp_name in rp_names) {
  area_col <- paste0("flood_area_", rp_name, "_m2")
  share_col <- paste0("flood_share_", rp_name)
  quality[[paste0("share_bounds_ok_", rp_name)]] <-
    quality[[share_col]] >= -tol & quality[[share_col]] <= 1 + tol
  quality[[paste0("area_bounds_ok_", rp_name)]] <-
    quality[[area_col]] >= -tol &
    quality[[area_col]] <= quality$municipality_area_m2 + tol
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
  difference <- paste0("delta_", pair[2], "_minus_", pair[1])
  check <- paste0("mono_", pair[1], "_le_", pair[2])
  quality[[difference]] <- quality[[right]] - quality[[left]]
  quality[[check]] <- quality[[difference]] >= -tol
}

share_checks <- grep("^share_bounds_ok_", names(quality), value = TRUE)
area_checks <- grep("^area_bounds_ok_", names(quality), value = TRUE)
monotonic_checks <- grep("^mono_", names(quality), value = TRUE)

quality <- quality %>%
  mutate(
    share_bounds_all_ok = if_all(all_of(share_checks), identity),
    area_bounds_all_ok = if_all(all_of(area_checks), identity),
    monotonicity_all_ok = if_all(all_of(monotonic_checks), identity),
    suspicious_any = !(share_bounds_all_ok & area_bounds_all_ok & monotonicity_all_ok)
  )

write_csv(quality, file.path(table_dir, "exposure_quality_checks_projected_area.csv"))
write_csv(
  quality %>% filter(suspicious_any),
  file.path(table_dir, "exposure_quality_checks_projected_area_suspicious_only.csv")
)

quality_summary <- tibble(
  municipalities = nrow(quality),
  projected_cell_area_m2 = prod(terra::res(terra::rast(rp_files[[1]]))),
  all_share_bounds_ok = all(quality$share_bounds_all_ok),
  all_area_bounds_ok = all(quality$area_bounds_all_ok),
  monotonicity_violations = sum(!quality$monotonicity_all_ok),
  suspicious_municipalities = sum(quality$suspicious_any)
)
write_csv(quality_summary, file.path(table_dir, "exposure_quality_summary_projected_area.csv"))

previous <- read_csv(
  paths$previous_exposure,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

comparison <- previous %>%
  select(
    AGS,
    previous_rp100_area_m2 = flood_area_rp100_m2,
    previous_rp100_share = flood_share_rp100
  ) %>%
  left_join(
    exposure %>%
      select(
        AGS,
        projected_rp100_area_m2 = flood_area_rp100_m2,
        projected_rp100_share = flood_share_rp100
      ),
    by = "AGS"
  ) %>%
  mutate(
    rp100_area_change_pct = 100 *
      (projected_rp100_area_m2 - previous_rp100_area_m2) /
      previous_rp100_area_m2,
    rp100_share_change_pp = 100 *
      (projected_rp100_share - previous_rp100_share)
  )

write_csv(
  comparison,
  file.path(table_dir, "geodesic_vs_projected_cell_area_comparison.csv")
)

log_message("Quality summary:")
capture.output(print(quality_summary), file = log_file, append = TRUE)
print(quality_summary)
