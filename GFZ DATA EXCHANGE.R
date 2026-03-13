#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(sf)
  library(ggplot2)
})

# ============================================================
# GFZ DATA EXCHANGE
# ------------------------------------------------------------
# Purpose:
# - reproduce the original 52-variable wide PCA input set
# - link each project variable back to the original INKAR code
# - calculate Germany-wide municipality coverage for these 52 indicators
# - export a clean variable list and municipality coverage tables
# - create a Germany map showing municipal coverage across all 52 variables
# ============================================================

raw_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/SOCIOECONOMIC.nosync/Downloads/inkar_2025/inkar_2025.csv"
meta_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/SOCIOECONOMIC.nosync/INKAR Uebersicht der Indikatoren.xlsx"
geom_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/vg250_gemeinden_landonly.gpkg"
out_dir <- file.path(getwd(), "outputs_gfz_data_exchange")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_file)) stop("Raw INKAR CSV not found: ", raw_file)
if (!file.exists(meta_file)) stop("INKAR metadata workbook not found: ", meta_file)
if (!file.exists(geom_file)) stop("Municipality geometry not found: ", geom_file)

save_csv <- function(x, name) {
  fwrite(x, file.path(out_dir, name))
}

save_md <- function(lines, name) {
  writeLines(lines, con = file.path(out_dir, name))
}

safe_empty_to_na <- function(x) {
  if (!is.character(x)) return(x)
  y <- trimws(x)
  y[y == ""] <- NA_character_
  y
}

all_states <- sprintf("%02d", 1:16)

wide_pca_map <- data.table(
  inkar_code = c(
    "a_ALGII_SGBII",
    "a_BG1P",
    "a_BG5um",
    "a_BGKind",
    "a_Unterkunft_SGBII",
    "a_aloLang",
    "q_alo_u25_einw",
    "q_alo_ü55_einw",
    "q_einkst_bev",
    "q_kaufkraft",
    "q_gewst_bev",
    "d_steuereinnahme",
    "a_bev0003",
    "a_bev0306",
    "a_bev0618",
    "a_bev1825",
    "a_bev2530",
    "a_bev3050",
    "a_bev5065",
    "a_bev6575",
    "a_bev65um",
    "a_bev75um",
    "i_saldo_nat",
    "i_wans",
    "q_abhg_alt",
    "q_abhg_jung",
    "r_ewf_jungalt",
    "q_HH1",
    "a_hh_kind",
    "a_hheink_hoch",
    "a_hheink_mittel",
    "a_hheink_niedrig",
    "q_stud",
    "q_stud_1825",
    "q_stud_fh",
    "q_allgemeinärzte_bev",
    "q_hausarzt_bev",
    "q_internist_bev",
    "q_kinderarzt_kinder",
    "q_ärzte_bev",
    "a_bb_1000Mbits",
    "a_bb_100Mbits",
    "a_bb_50Mbits",
    "a_bb_4G",
    "m_G02_SUP_DIST",
    "m_Q01_APO_DIST",
    "m_Q07_HA_DIST",
    "m_OEV20_DIST",
    "m_P01_PRIM_DIST",
    "q_bev_fl",
    "q_bevsva_qkm"
  ),
  project_name = c(
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
)

wide_pca_map[, analysis_name := paste0("exposure_", project_name)]

if (nrow(wide_pca_map) != 51L) {
  stop("Expected 51 wide PCA variables from Analyse_CLEAN.R, found: ", nrow(wide_pca_map))
}

message("1) Reading INKAR metadata workbook ...")

meta_raw <- read_excel(
  meta_file,
  sheet = "Raumbeobachtung DE",
  skip = 1,
  col_names = c(
    "Kurzname", "Name", "Algorithmus", "M_ID", "Kuerzel",
    "Anmerkungen", "Grundlagen", "Gemeinden", "Kreise"
  )
)

meta_raw[] <- lapply(meta_raw, safe_empty_to_na)
meta_raw$row_id <- seq_len(nrow(meta_raw))

current_main <- NA_character_
current_sub <- NA_character_
main_vec <- character(nrow(meta_raw))
sub_vec <- character(nrow(meta_raw))

for (i in seq_len(nrow(meta_raw))) {
  short_label <- meta_raw$Kurzname[i]
  kuerzel <- meta_raw$Kuerzel[i]

  if (is.na(kuerzel) && !is.na(short_label) && short_label != "Kurzname") {
    if (grepl("\\s[-\u2013]\\s", short_label)) {
      current_sub <- short_label
    } else {
      current_main <- short_label
      current_sub <- NA_character_
    }
  }

  main_vec[i] <- current_main
  sub_vec[i] <- current_sub
}

meta <- as.data.table(meta_raw)[
  !is.na(Kuerzel) & Kuerzel != "K\u00fcrzel",
  .(
    inkar_code = Kuerzel,
    official_short = Kurzname,
    official_name = Name,
    main_dimension = main_vec[row_id],
    sub_dimension = sub_vec[row_id]
  )
]

message("2) Reading raw INKAR CSV and filtering the 51 original wide PCA indicators ...")

inkar <- fread(
  raw_file,
  sep = ";",
  quote = "",
  select = c("Bereich", "Raumbezug", "Kennziffer", "Kuerzel", "Wert", "Zeitbezug"),
  encoding = "UTF-8",
  showProgress = FALSE
)

inkar <- inkar[
  Bereich == "LRB" &
    Raumbezug == "Gemeinden" &
    substr(Kennziffer, 1, 2) %in% all_states &
    Kuerzel %in% wide_pca_map$inkar_code
]

inkar[, state := substr(Kennziffer, 1, 2)]

message("3) Keeping latest observation per municipality x indicator ...")

setorder(inkar, Kennziffer, Kuerzel, -Zeitbezug)
latest <- inkar[, .SD[1], by = .(Kennziffer, state, Kuerzel)]

all_municipalities <- sort(unique(latest$Kennziffer))
indicator_template <- CJ(Kennziffer = all_municipalities, inkar_code = wide_pca_map$inkar_code)

latest[, value_chr := trimws(as.character(Wert))]
latest[, has_value := !is.na(Wert) & !(value_chr %in% c("", ".", "..", "-"))]

latest_small <- latest[, .(
  Kennziffer,
  state,
  inkar_code = Kuerzel,
  has_value,
  Zeitbezug
)]

coverage_long <- merge(
  indicator_template,
  latest_small,
  by = c("Kennziffer", "inkar_code"),
  all.x = TRUE
)

coverage_long[is.na(state), state := substr(Kennziffer, 1, 2)]
coverage_long[is.na(has_value), has_value := FALSE]

message("4) Calculating variable coverage and municipality coverage ...")

total_municipalities <- uniqueN(coverage_long$Kennziffer)

variable_coverage <- coverage_long[
  ,
  .(
    municipalities_with_value = sum(has_value),
    municipalities_total = .N,
    coverage_pct = 100 * mean(has_value),
    missing_municipalities = sum(!has_value)
  ),
  by = inkar_code
]

state_coverage <- coverage_long[
  ,
  .(
    municipalities_with_value = sum(has_value),
    municipalities_total = .N,
    state_coverage_pct = 100 * mean(has_value)
  ),
  by = .(inkar_code, state)
]

state_summary <- state_coverage[
  ,
  .(
    states_with_any_data = sum(municipalities_with_value > 0),
    present_in_all_16_states = all(municipalities_with_value > 0),
    full_coverage_all_16_states = all(state_coverage_pct == 100),
    min_state_coverage_pct = min(state_coverage_pct),
    mean_state_coverage_pct = mean(state_coverage_pct),
    max_state_coverage_pct = max(state_coverage_pct)
  ),
  by = inkar_code
]

variable_coverage <- merge(variable_coverage, state_summary, by = "inkar_code", all.x = TRUE)
variable_coverage <- merge(variable_coverage, wide_pca_map, by = "inkar_code", all.x = TRUE)
variable_coverage <- merge(variable_coverage, meta, by = "inkar_code", all.x = TRUE)

variable_coverage[inkar_code == "a_bb_4G", `:=`(
  official_short = fifelse(
    is.na(official_short) | official_short == "",
    "Bandbreitenverfugbarkeit mindestens 4G",
    official_short
  ),
  official_name = fifelse(
    is.na(official_name) | official_name == "",
    "Anteil der Haushalte mit einer 4G-Mobilfunkversorgung in %",
    official_name
  ),
  main_dimension = fifelse(
    is.na(main_dimension) | main_dimension == "",
    "Erreichbarkeit",
    main_dimension
  ),
  sub_dimension = fifelse(
    is.na(sub_dimension) | sub_dimension == "",
    "Digitale Erreichbarkeit",
    sub_dimension
  )
)]

variable_coverage[is.na(main_dimension), main_dimension := "(missing in metadata)"]
variable_coverage[is.na(sub_dimension), sub_dimension := "(no sub-dimension)"]

municipality_coverage <- coverage_long[
  ,
  .(
    n_variables_available = sum(has_value),
    n_variables_total = .N,
    coverage_pct = 100 * mean(has_value)
  ),
  by = Kennziffer
]

municipality_coverage[, state := substr(Kennziffer, 1, 2)]

message("5) Reading municipality geometry and building Germany map ...")

mun_geom <- st_read(geom_file, quiet = TRUE)
mun_geom <- mun_geom[, c("Gemeindeschlüssel_AGS", "Gemeinde")]
names(mun_geom) <- c("AGS", "mun_name")

map_idx <- match(mun_geom$AGS, municipality_coverage$Kennziffer)
map_attrs <- municipality_coverage[map_idx, .(
  n_variables_available,
  n_variables_total,
  coverage_pct,
  state
)]
map_sf <- cbind(mun_geom, as.data.frame(map_attrs))

join_summary <- data.table(
  municipalities_in_inkar = total_municipalities,
  polygons_in_geometry = nrow(mun_geom),
  polygons_with_coverage = sum(!is.na(map_sf$coverage_pct)),
  polygons_without_coverage = sum(is.na(map_sf$coverage_pct)),
  inkar_ags_missing_in_geometry = total_municipalities - sum(municipality_coverage$Kennziffer %in% mun_geom$AGS)
)

map_plot <- ggplot(map_sf) +
  geom_sf(aes(fill = coverage_pct), color = NA, size = 0) +
  scale_fill_gradientn(
    colours = c("#efe3d3", "#d6c0a9", "#8ca79d", "#2b6f77"),
    limits = c(0, 100),
    na.value = "grey90",
    name = "Coverage (%)"
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Municipality-level availability across the original wide PCA variables",
    subtitle = "Share of the 51 original INKAR indicators with a valid latest value per municipality",
    caption = "Source: INKAR 2025, municipality level (Bereich = LRB, Raumbezug = Gemeinden). Own processing."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8),
    legend.position = c(0.12, 0.18),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA)
  )

ggsave(
  filename = file.path(out_dir, "map_germany_wide_pca_coverage.png"),
  plot = map_plot,
  width = 11,
  height = 8.5,
  dpi = 320
)

bar_plot <- ggplot(
  variable_coverage[order(coverage_pct)],
  aes(x = reorder(project_name, coverage_pct), y = coverage_pct)
) +
  geom_col(fill = "#2b6f77") +
  coord_flip() +
  labs(
    title = "Coverage of the original wide PCA variables",
    subtitle = "Germany-wide municipality coverage based on the latest available INKAR value",
    x = NULL,
    y = "Coverage (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "plot_wide_pca_variable_coverage.png"),
  plot = bar_plot,
  width = 9.5,
  height = 11,
  dpi = 320
)

message("6) Writing output files ...")

setorder(variable_coverage, coverage_pct, project_name)
setorder(municipality_coverage, coverage_pct, Kennziffer)
setorder(state_coverage, inkar_code, state)

variables_used_with_coverage <- merge(
  wide_pca_map,
  variable_coverage,
  by = c("inkar_code", "project_name", "analysis_name"),
  all.x = TRUE,
  sort = FALSE
)

save_csv(variable_coverage, "wide_pca_variable_coverage_germany.csv")
save_csv(municipality_coverage, "wide_pca_municipality_coverage_germany.csv")
save_csv(state_coverage, "wide_pca_state_coverage_germany.csv")
save_csv(join_summary, "wide_pca_map_join_summary.csv")
save_csv(wide_pca_map, "wide_pca_variables_used.csv")
save_csv(variables_used_with_coverage, "wide_pca_variables_used_with_coverage.csv")

md_lines <- c(
  "# Original Wide PCA Variables Used in Analyse_CLEAN.R",
  "",
  sprintf("Total variables in the original wide PCA input set: %d", nrow(wide_pca_map)),
  "",
  "Each bullet shows the original INKAR code, the project rename used in the thesis code, and the Germany-wide municipality coverage based on the latest available INKAR value.",
  ""
)

for (i in seq_len(nrow(variables_used_with_coverage))) {
  row <- variables_used_with_coverage[i]
  md_lines <- c(
    md_lines,
    sprintf("- `%s` -> `%s` (`%s`)", row$inkar_code, row$project_name, row$analysis_name),
    sprintf("  Official INKAR name: %s", row$official_short),
    sprintf("  Germany-wide municipality coverage: %.2f%% (%d / %d municipalities)", row$coverage_pct, row$municipalities_with_value, row$municipalities_total),
    sprintf("  State coverage range: %.2f%% to %.2f%%", row$min_state_coverage_pct, row$max_state_coverage_pct),
    ""
  )
}

save_md(md_lines, "wide_pca_variables_used.md")

message("")
message("GFZ DATA EXCHANGE completed.")
message("- Original wide PCA variables documented: ", nrow(wide_pca_map))
message("- Germany municipalities in INKAR subset: ", total_municipalities)
message("- Coverage table written to: ", file.path(out_dir, "wide_pca_variable_coverage_germany.csv"))
message("- Municipality coverage map written to: ", file.path(out_dir, "map_germany_wide_pca_coverage.png"))
