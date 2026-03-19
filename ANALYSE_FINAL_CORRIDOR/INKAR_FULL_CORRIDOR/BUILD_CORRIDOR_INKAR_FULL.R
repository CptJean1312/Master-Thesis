#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(sf)
  library(stringr)
  library(readr)
})

# -------------------------------------------------------------------
# Build a full INKAR corridor layer directly from the original INKAR CSV.
# This keeps the original INKAR codes and official indicator names.
# -------------------------------------------------------------------

paths <- list(
  corridor_gpkg = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg",
  inkar_raw_csv = "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/SOCIOECONOMIC.nosync/Downloads/inkar_2025/inkar_2025.csv",
  inkar_meta_csv = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_inkar_overview/inkar_indicator_availability_all_states.csv",
  output_dir = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs"
)

dir.create(file.path(paths$output_dir, "gpkg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(paths$output_dir, "logs", "build_corridor_inkar_full_log.txt")
log_lines <- character()

log_message <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  log_lines <<- c(log_lines, line)
  message(line)
}

safe_write_gpkg <- function(x, path, layer) {
  if (file.exists(path)) unlink(path)
  write_sf(x, path, layer = layer, quiet = TRUE)
}

standardize_ags <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  str_pad(x, width = 8, side = "left", pad = "0")
}

log_message("Loading corridor municipalities.")
corridor <- st_read(paths$corridor_gpkg, quiet = TRUE)

corridor <- corridor %>%
  mutate(
    AGS_corridor = dplyr::coalesce(
      standardize_ags(AGS),
      standardize_ags(Gemeindeschlüssel_AGS)
    ),
    mun_name_corridor = dplyr::coalesce(
      as.character(mun_name),
      as.character(Gemeinde),
      as.character(GeografischerName_GEN)
    )
  )

corridor_keys <- corridor %>%
  st_drop_geometry() %>%
  transmute(
    AGS = AGS_corridor,
    mun_name_corridor
  ) %>%
  distinct()

if (anyNA(corridor_keys$AGS)) {
  stop("At least one corridor municipality still has no valid AGS after key consolidation.")
}

log_message("Corridor municipalities loaded: ", nrow(corridor_keys))

log_message("Reading original INKAR raw CSV.")
inkar_raw <- fread(
  paths$inkar_raw_csv,
  select = c("Bereich", "Raumbezug", "Kennziffer", "Name", "Kuerzel", "Indikator", "Zeitbezug", "Wert"),
  showProgress = TRUE
)

inkar_raw[, Kennziffer := standardize_ags(Kennziffer)]
inkar_raw[, Zeitbezug := as.integer(Zeitbezug)]

inkar_corridor_long <- inkar_raw[
  Bereich == "LRB" &
    Raumbezug == "Gemeinden" &
    Kennziffer %in% corridor_keys$AGS
]

rm(inkar_raw)
gc()

log_message("Filtered INKAR long rows for corridor: ", nrow(inkar_corridor_long))

setorder(inkar_corridor_long, Kennziffer, Kuerzel, -Zeitbezug)
inkar_latest <- inkar_corridor_long[, .SD[1], by = .(Kennziffer, Kuerzel)]

log_message(
  "Latest municipality x indicator rows retained: ",
  nrow(inkar_latest),
  " across ",
  uniqueN(inkar_latest$Kennziffer),
  " municipalities and ",
  uniqueN(inkar_latest$Kuerzel),
  " indicators."
)

indicator_meta <- unique(inkar_latest[, .(Kuerzel, Indikator)])
setorder(indicator_meta, Kuerzel)

values_wide <- dcast(
  inkar_latest,
  Kennziffer + Name ~ Kuerzel,
  value.var = "Wert"
)

years_wide <- dcast(
  inkar_latest,
  Kennziffer + Name ~ Kuerzel,
  value.var = "Zeitbezug"
)

year_cols <- setdiff(names(years_wide), c("Kennziffer", "Name"))
setnames(years_wide, year_cols, paste0(year_cols, "_year"))

inkar_wide <- merge(values_wide, years_wide, by = c("Kennziffer", "Name"), all = TRUE)
setnames(inkar_wide, c("Kennziffer", "Name"), c("AGS", "inkar_name"))

for (col in setdiff(names(inkar_wide), c("AGS", "inkar_name", grep("_year$", names(inkar_wide), value = TRUE)))) {
  suppressWarnings(set(inkar_wide, j = col, value = type.convert(inkar_wide[[col]], as.is = TRUE)))
}

for (col in grep("_year$", names(inkar_wide), value = TRUE)) {
  suppressWarnings(set(inkar_wide, j = col, value = as.integer(inkar_wide[[col]])))
}

inventory <- merge(
  data.table(AGS = corridor_keys$AGS),
  inkar_wide,
  by = "AGS",
  all.x = TRUE
)

indicator_codes <- indicator_meta$Kuerzel

coverage_summary <- rbindlist(lapply(indicator_codes, function(code) {
  year_col <- paste0(code, "_year")
  values <- inventory[[code]]
  years <- inventory[[year_col]]

  data.table(
    Kuerzel = code,
    Indikator = indicator_meta[Kuerzel == code, Indikator][1],
    n_non_missing = sum(!is.na(values)),
    n_missing = sum(is.na(values)),
    coverage_pct = round(mean(!is.na(values)) * 100, 2),
    min_year = if (all(is.na(years))) NA_integer_ else min(years, na.rm = TRUE),
    max_year = if (all(is.na(years))) NA_integer_ else max(years, na.rm = TRUE)
  )
}))

setorder(coverage_summary, -coverage_pct, Kuerzel)

pair_summary <- inventory[, .(
  n_complete_pairs = sum(complete.cases(a_ALGII_SGBII, a_Unterkunft_SGBII)),
  correlation = suppressWarnings(cor(a_ALGII_SGBII, a_Unterkunft_SGBII, use = "pairwise.complete.obs"))
)]

missing_ags <- setdiff(corridor_keys$AGS, inkar_wide$AGS)

corridor_full <- corridor %>%
  transmute(
    AGS = AGS_corridor,
    mun_name_corridor = mun_name_corridor
  ) %>%
  left_join(as.data.frame(inkar_wide), by = "AGS")

output_gpkg <- file.path(paths$output_dir, "gpkg", "corridor_full_inkar_original_codes.gpkg")
safe_write_gpkg(corridor_full, output_gpkg, "corridor_full_inkar")

write_csv(as.data.frame(coverage_summary), file.path(paths$output_dir, "tables", "corridor_inkar_indicator_inventory.csv"))
write_csv(as.data.frame(indicator_meta), file.path(paths$output_dir, "tables", "corridor_inkar_indicator_meta.csv"))
write_csv(as.data.frame(pair_summary), file.path(paths$output_dir, "tables", "corridor_sgb2_pair_check.csv"))
write_csv(as.data.frame(inkar_wide), file.path(paths$output_dir, "tables", "corridor_full_inkar_original_codes.csv"))
write_csv(data.frame(AGS = missing_ags), file.path(paths$output_dir, "tables", "corridor_inkar_missing_ags.csv"))

if (file.exists(paths$inkar_meta_csv)) {
  meta_overview <- read_csv(paths$inkar_meta_csv, show_col_types = FALSE)

  coverage_dimensioned <- coverage_summary %>%
    left_join(
      meta_overview %>%
        select(Kuerzel, official_name, main_dimension, sub_dimension, broad_dimension),
      by = "Kuerzel"
    )

  dimension_summary <- coverage_dimensioned %>%
    group_by(broad_dimension) %>%
    summarise(
      indicators = n(),
      median_coverage_pct = median(coverage_pct, na.rm = TRUE),
      min_coverage_pct = min(coverage_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(indicators))

  write_csv(
    coverage_dimensioned,
    file.path(paths$output_dir, "tables", "corridor_inkar_indicator_inventory_dimensioned.csv")
  )
  write_csv(
    dimension_summary,
    file.path(paths$output_dir, "tables", "corridor_inkar_dimension_summary.csv")
  )
}

headline <- data.table(
  corridor_municipalities = nrow(corridor_keys),
  matched_municipalities_with_any_inkar = sum(inventory[, rowSums(!is.na(.SD)) > 0, .SDcols = indicator_codes]),
  municipalities_missing_any_inkar_row = length(missing_ags),
  indicators_available = length(indicator_codes),
  median_indicator_coverage_pct = median(coverage_summary$coverage_pct),
  min_indicator_coverage_pct = min(coverage_summary$coverage_pct),
  max_indicator_coverage_pct = max(coverage_summary$coverage_pct)
)

write_csv(as.data.frame(headline), file.path(paths$output_dir, "tables", "corridor_inkar_headline.csv"))

log_message("Wrote GPKG: ", output_gpkg)
log_message("Wrote indicator inventory and headline tables.")
log_message(
  "Pair check a_ALGII_SGBII vs a_Unterkunft_SGBII: n_complete=",
  pair_summary$n_complete_pairs,
  ", correlation=",
  round(pair_summary$correlation, 6)
)

writeLines(log_lines, log_path)
log_message("Done.")
