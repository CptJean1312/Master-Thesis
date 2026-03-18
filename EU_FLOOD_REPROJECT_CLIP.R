#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(terra)
})

source_dir <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/EU_Flood_Maps/Flood_Hazard"
basin_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/ELBE BASIN PRO.gpkg"
basin_layer <- "elbe_basin"

out_dir <- file.path(getwd(), "outputs_eu_flood_25832")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_crs <- "EPSG:25832"
target_res <- 100
overwrite_outputs <- identical(Sys.getenv("EU_FLOOD_OVERWRITE", unset = "FALSE"), "TRUE")

wopt_float <- list(
  datatype = "FLT4S",
  gdal = c(
    "COMPRESS=DEFLATE",
    "PREDICTOR=2",
    "TILED=YES",
    "BIGTIFF=IF_SAFER"
  )
)

if (!dir.exists(source_dir)) stop("Missing source directory: ", source_dir)
if (!file.exists(basin_file)) stop("Missing basin layer: ", basin_file)

source_files <- sort(list.files(
  source_dir,
  pattern = "^floodmap_EFAS_RP[0-9]{3}_C\\.tif$",
  full.names = TRUE
))

if (length(source_files) == 0) {
  stop("No EFAS RP rasters found in: ", source_dir)
}

message("Loading Elbe basin polygon ...")
basin_25832 <- vect(basin_file, layer = basin_layer)

if (is.na(crs(basin_25832))) {
  stop("Elbe basin layer has no CRS.")
}

if (!grepl("25832", crs(basin_25832))) {
  basin_25832 <- project(basin_25832, target_crs)
}

aligned_extent <- align(ext(basin_25832), target_res)
template_25832 <- rast(aligned_extent, resolution = target_res, crs = target_crs)

process_raster <- function(raster_file, basin_25832, template_25832, out_dir, overwrite_outputs) {
  raster_name <- tools::file_path_sans_ext(basename(raster_file))
  out_file <- file.path(out_dir, paste0(raster_name, "_25832_elbe_basin.tif"))

  if (file.exists(out_file) && !overwrite_outputs) {
    message("")
    message("Skipping existing output: ", basename(out_file))
    return(data.frame(
      raster_name = raster_name,
      output_file = out_file,
      status = "skipped_existing",
      stringsAsFactors = FALSE
    ))
  }

  message("")
  message("Processing: ", raster_name)
  src <- rast(raster_file)

  message("Source CRS: ", crs(src, proj = TRUE))
  message("Source resolution: ", paste(res(src), collapse = " x "))

  basin_in_src_crs <- project(basin_25832, crs(src))

  message("Cropping and masking in source CRS to reduce runtime ...")
  src_crop <- crop(src, basin_in_src_crs, snap = "out")
  src_mask <- mask(src_crop, basin_in_src_crs)

  message("Projecting to EPSG:25832 ...")
  projected <- project(src_mask, template_25832, method = "bilinear")

  message("Applying final Elbe basin mask in EPSG:25832 ...")
  projected_masked <- mask(projected, basin_25832)

  message("Writing: ", out_file)
  writeRaster(projected_masked, out_file, overwrite = TRUE, wopt = wopt_float)

  data.frame(
    raster_name = raster_name,
    output_file = out_file,
    status = "written",
    stringsAsFactors = FALSE
  )
}

results <- do.call(
  rbind,
  lapply(
    source_files,
    process_raster,
    basin_25832 = basin_25832,
    template_25832 = template_25832,
    out_dir = out_dir,
    overwrite_outputs = overwrite_outputs
  )
)

write.csv(
  results,
  file = file.path(out_dir, "reprojection_run_log.csv"),
  row.names = FALSE
)

message("")
message("Done.")
message("Total source rasters found: ", length(source_files))
message("Written: ", sum(results$status == "written"))
message("Skipped existing: ", sum(results$status == "skipped_existing"))
message("Run log: ", file.path(out_dir, "reprojection_run_log.csv"))
