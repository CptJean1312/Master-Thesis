#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"
table_dir <- file.path(
  root,
  "ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/outputs/tables"
)
input_csv <- file.path(table_dir, "full_municipality_exposure_all_RPs_projected_area.csv")

exposure <- read_csv(
  input_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

rp_names <- c("rp10", "rp20", "rp50", "rp100", "rp200", "rp500")
share_tolerance <- 1e-8
area_tolerance_m2 <- 1
monotonicity_tolerance <- 1e-9
quality <- exposure

for (rp_name in rp_names) {
  share_col <- paste0("flood_share_", rp_name)
  area_col <- paste0("flood_area_", rp_name, "_m2")
  quality[[paste0("share_bounds_ok_", rp_name)]] <-
    quality[[share_col]] >= -share_tolerance &
    quality[[share_col]] <= 1 + share_tolerance
  quality[[paste0("area_bounds_ok_", rp_name)]] <-
    quality[[area_col]] >= -area_tolerance_m2 &
    quality[[area_col]] <= quality$municipality_area_m2 + area_tolerance_m2
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
  quality[[check]] <- quality[[difference]] >= -monotonicity_tolerance
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

write_csv(
  quality,
  file.path(table_dir, "exposure_quality_checks_projected_area_with_numeric_tolerance.csv")
)
write_csv(
  quality %>% filter(suspicious_any),
  file.path(table_dir, "exposure_quality_checks_projected_area_suspicious_with_numeric_tolerance.csv")
)

summary <- tibble(
  municipalities = nrow(quality),
  share_tolerance = share_tolerance,
  area_tolerance_m2 = area_tolerance_m2,
  all_share_bounds_ok = all(quality$share_bounds_all_ok),
  all_area_bounds_ok = all(quality$area_bounds_all_ok),
  monotonicity_violations = sum(!quality$monotonicity_all_ok),
  suspicious_municipalities = sum(quality$suspicious_any)
)

write_csv(
  summary,
  file.path(table_dir, "exposure_quality_summary_projected_area_with_numeric_tolerance.csv")
)
print(summary)
