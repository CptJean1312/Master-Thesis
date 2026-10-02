#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"
input_csv <- file.path(
  root,
  "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_GEOMETRY_AUDIT_2026-10-01/outputs/tables/full_municipality_exposure_all_RPs.csv"
)
output_dir <- file.path(
  root,
  "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_CORRECTED_2026-10-01/exposure_quality_outputs"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

exposure <- read_csv(
  input_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

rp_names <- c("rp10", "rp20", "rp50", "rp100", "rp200", "rp500")
tol <- 1e-9
quality <- exposure

for (rp_name in rp_names) {
  share_col <- paste0("flood_share_", rp_name)
  area_col <- paste0("flood_area_", rp_name, "_m2")
  quality[[paste0("share_bounds_ok_", rp_name)]] <-
    !is.na(quality[[share_col]]) &
    quality[[share_col]] >= -tol &
    quality[[share_col]] <= 1 + tol
  quality[[paste0("area_bounds_ok_", rp_name)]] <-
    !is.na(quality[[area_col]]) &
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

write_csv(quality, file.path(output_dir, "exposure_quality_checks.csv"))
write_csv(
  quality %>% filter(suspicious_any),
  file.path(output_dir, "exposure_quality_checks_suspicious_only.csv")
)

summary <- tibble(
  municipalities = nrow(quality),
  all_share_bounds_ok = all(quality$share_bounds_all_ok),
  all_area_bounds_ok = all(quality$area_bounds_all_ok),
  monotonicity_violations = sum(!quality$monotonicity_all_ok),
  suspicious_municipalities = sum(quality$suspicious_any)
)

write_csv(summary, file.path(output_dir, "exposure_quality_summary.csv"))
print(summary)
