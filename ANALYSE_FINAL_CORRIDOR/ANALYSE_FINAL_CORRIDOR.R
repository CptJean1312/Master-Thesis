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

corridor_inkar_original_codes_csv <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_full_inkar_original_codes.csv"

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

scale_to_z <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

build_direction_coded_inputs <- function(data_imputed, reverse_variables) {
  scaled <- as.data.frame(scale(data_imputed))
  direction_coded <- scaled
  direction_coded[reverse_variables] <- lapply(direction_coded[reverse_variables], function(x) -x)
  direction_coded
}

build_dimension_balanced_index <- function(direction_coded_data, dimensions) {
  dimension_scores <- setNames(lapply(names(dimensions), function(dimension_name) {
    vars <- dimensions[[dimension_name]]
    raw_score <- rowMeans(direction_coded_data[, vars, drop = FALSE], na.rm = TRUE)
    scale_to_z(raw_score)
  }), names(dimensions))

  dimension_scores <- as_tibble(dimension_scores)
  names(dimension_scores) <- paste0("vuln_dim_", names(dimensions), "_z")

  raw_index <- rowMeans(dimension_scores, na.rm = TRUE)

  list(
    raw = raw_index,
    z = scale_to_z(raw_index),
    dimension_scores = dimension_scores
  )
}

build_aligned_pca_index <- function(direction_coded_data, anchor_index, n_components) {
  pca <- prcomp(direction_coded_data)
  eigenvalues <- pca$sdev^2
  variance_share <- eigenvalues / sum(eigenvalues)

  component_scores <- as.data.frame(pca$x)
  names(component_scores) <- paste0("PC", seq_len(ncol(component_scores)))

  retained_pcs <- paste0("PC", seq_len(n_components))
  pc_alignment <- sapply(retained_pcs, function(pc) {
    pc_correlation <- suppressWarnings(cor(component_scores[[pc]], anchor_index, use = "complete.obs"))
    ifelse(!is.na(pc_correlation) && pc_correlation < 0, -1, 1)
  })

  aligned_scores <- sweep(
    as.matrix(component_scores[, retained_pcs, drop = FALSE]),
    2,
    pc_alignment,
    `*`
  )

  weights <- variance_share[seq_len(n_components)]
  raw_index <- as.numeric(aligned_scores %*% (weights / sum(weights)))

  list(
    pca = pca,
    scores = component_scores,
    aligned_scores = aligned_scores,
    eigenvalues = eigenvalues,
    variance_share = variance_share,
    alignment = pc_alignment,
    raw = raw_index,
    z = scale_to_z(raw_index)
  )
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

corridor_sf <- corridor_sf %>%
  mutate(
    AGS = coalesce(as.character(AGS), as.character(Gemeindeschlüssel_AGS)),
    mun_name = coalesce(mun_name, GeografischerName_GEN)
  )

log_message("Corridor municipalities loaded: ", nrow(corridor_sf))

log_message("Loading municipal exposure table ...")
exposure_tbl <- read_csv(
  exposure_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
)

log_message("Exposure rows loaded: ", nrow(exposure_tbl))

if (
  sum(is.na(exposure_tbl$AGS)) == 1 &&
    "16076094" %in% corridor_sf$AGS &&
    !"16076094" %in% exposure_tbl$AGS
) {
  exposure_tbl$AGS[is.na(exposure_tbl$AGS)] <- "16076094"
  exposure_tbl$mun_name[exposure_tbl$AGS == "16076094" & is.na(exposure_tbl$mun_name)] <- "Berga-Wünschendorf"
  log_message("Repaired one missing exposure-table AGS as 16076094 (Berga-Wünschendorf).")
}

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

vulnerability_variables <- tribble(
  ~variable, ~dimension, ~indicator_label, ~expected_direction, ~reverse_for_index,
  "a_bev65um", "demographic_household", "Residents aged 65 years and older", "higher vulnerability", FALSE,
  "a_hh_kind", "demographic_household", "Households with children", "higher vulnerability", FALSE,
  "q_HH1", "demographic_household", "Single-person households", "higher vulnerability", FALSE,
  "a_bev_0006", "demographic_household", "Residents under 6 years", "higher vulnerability", FALSE,
  "a_ewfBG_allein", "demographic_household", "Single-parent employable SGB-II recipients", "higher vulnerability", FALSE,
  "a_ALGII_SGBII", "deprivation_labour", "ALG-II share within SGB-II benefits", "higher vulnerability", FALSE,
  "a_hheink_niedrig", "deprivation_labour", "Low-income households", "higher vulnerability", FALSE,
  "q_einkst_bev", "deprivation_labour", "Income tax per capita", "higher capacity", TRUE,
  "q_kaufkraft", "deprivation_labour", "Purchasing power per capita", "higher capacity", TRUE,
  "q_newfBGu15_bev", "deprivation_labour", "Child poverty", "higher vulnerability", FALSE,
  "a_aloLang", "deprivation_labour", "Long-term unemployed", "higher vulnerability", FALSE,
  "q_alo_u25_einw", "deprivation_labour", "Unemployed under 25", "higher vulnerability", FALSE,
  "a_Minijobs", "deprivation_labour", "Marginal employment", "higher vulnerability", FALSE,
  "q_alo_ü55_einw", "deprivation_labour", "Unemployed aged 55 and older", "higher vulnerability", FALSE,
  "q_svw", "deprivation_labour", "Employment rate", "higher capacity", TRUE,
  "m_G02_SUP_DIST", "access_adaptive_capacity", "Distance to supermarket/discounter", "higher vulnerability", FALSE,
  "m_OEV20_DIST", "access_adaptive_capacity", "Distance to public transport stop", "higher vulnerability", FALSE,
  "m_Q01_APO_DIST", "access_adaptive_capacity", "Distance to pharmacy", "higher vulnerability", FALSE,
  "m_Q07_HA_DIST", "access_adaptive_capacity", "Distance to general practitioner", "higher vulnerability", FALSE,
  "q_ärzte_bev", "access_adaptive_capacity", "Doctors per 10,000 residents", "higher capacity", TRUE,
  "a_bb_4G", "access_adaptive_capacity", "4G mobile broadband availability", "higher capacity", TRUE,
  "a_bb_50Mbits", "access_adaptive_capacity", "Fixed broadband availability at least 50 Mbit/s", "higher capacity", TRUE,
  "m_P01_PRIM_DIST", "access_adaptive_capacity", "Distance to primary school", "higher vulnerability", FALSE
)

pca_columns <- vulnerability_variables$variable
reverse_index_columns <- vulnerability_variables %>%
  filter(reverse_for_index) %>%
  pull(variable)

vulnerability_dimensions <- split(vulnerability_variables$variable, vulnerability_variables$dimension)

log_message("Loading final corridor INKAR table with original indicator codes ...")
inkar_original_tbl <- read_csv(
  corridor_inkar_original_codes_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
) %>%
  select(AGS, all_of(pca_columns))

corridor_sf <- corridor_sf %>%
  select(-any_of(pca_columns)) %>%
  left_join(inkar_original_tbl, by = "AGS")

missing_pca_columns <- setdiff(pca_columns, names(corridor_sf))
if (length(missing_pca_columns) > 0) {
  stop("Missing PCA columns: ", paste(missing_pca_columns, collapse = ", "))
}

vulnerability_raw_data <- corridor_sf %>%
  st_drop_geometry() %>%
  select(all_of(pca_columns)) %>%
  mutate(across(everything(), as.numeric))

has_inkar_record <- rowSums(!is.na(vulnerability_raw_data)) > 0
corridor_sf$has_inkar_vulnerability_record <- has_inkar_record

log_message("Municipalities with INKAR records for vulnerability index: ", sum(has_inkar_record))
log_message("Municipalities without INKAR records for vulnerability index: ", sum(!has_inkar_record))

missingness_tbl <- tibble(
  variable = names(vulnerability_raw_data),
  n_missing_all_corridor = colSums(is.na(vulnerability_raw_data)),
  share_missing_all_corridor = colMeans(is.na(vulnerability_raw_data)),
  n_missing_inkar_matched = colSums(is.na(vulnerability_raw_data[has_inkar_record, , drop = FALSE])),
  share_missing_inkar_matched = colMeans(is.na(vulnerability_raw_data[has_inkar_record, , drop = FALSE]))
) %>%
  left_join(vulnerability_variables, by = "variable") %>%
  arrange(dimension, variable)

save_table(missingness_tbl, "corridor_pca_missingness.csv")
save_table(vulnerability_variables, "corridor_vulnerability_index_variables.csv")

vulnerability_imputed <- vulnerability_raw_data[has_inkar_record, , drop = FALSE] %>%
  mutate(across(everything(), impute_median))

vulnerability_direction_coded <- build_direction_coded_inputs(
  data_imputed = vulnerability_imputed,
  reverse_variables = reverse_index_columns
)

main_index <- build_dimension_balanced_index(
  direction_coded_data = vulnerability_direction_coded,
  dimensions = vulnerability_dimensions
)

pca_preliminary <- prcomp(vulnerability_direction_coded)
prelim_eigenvalues <- pca_preliminary$sdev^2
prelim_variance_share <- prelim_eigenvalues / sum(prelim_eigenvalues)
n_components_70pct <- which(cumsum(prelim_variance_share) >= 0.7)[1]
n_components_kaiser <- sum(prelim_eigenvalues > 1)

pca_sensitivity_70pct <- build_aligned_pca_index(
  direction_coded_data = vulnerability_direction_coded,
  anchor_index = main_index$z,
  n_components = n_components_70pct
)

pca_sensitivity_kaiser <- build_aligned_pca_index(
  direction_coded_data = vulnerability_direction_coded,
  anchor_index = main_index$z,
  n_components = n_components_kaiser
)

component_scores <- pca_sensitivity_70pct$scores
component_cols <- paste0("PC", seq_len(ncol(component_scores)))
for (pc_name in component_cols) {
  corridor_sf[[pc_name]] <- NA_real_
  corridor_sf[[pc_name]][has_inkar_record] <- component_scores[[pc_name]]
}

dimension_score_cols <- names(main_index$dimension_scores)
for (dimension_col in dimension_score_cols) {
  corridor_sf[[dimension_col]] <- NA_real_
  corridor_sf[[dimension_col]][has_inkar_record] <- main_index$dimension_scores[[dimension_col]]
}

corridor_sf$vuln_index_main <- NA_real_
corridor_sf$vuln_index_main_z <- NA_real_
corridor_sf$vuln_index_pca23_70pct <- NA_real_
corridor_sf$vuln_index_pca23_70pct_z <- NA_real_
corridor_sf$vuln_index_pca23_kaiser <- NA_real_
corridor_sf$vuln_index_pca23_kaiser_z <- NA_real_

corridor_sf$vuln_index_main[has_inkar_record] <- main_index$raw
corridor_sf$vuln_index_main_z[has_inkar_record] <- main_index$z
corridor_sf$vuln_index_pca23_70pct[has_inkar_record] <- pca_sensitivity_70pct$raw
corridor_sf$vuln_index_pca23_70pct_z[has_inkar_record] <- pca_sensitivity_70pct$z
corridor_sf$vuln_index_pca23_kaiser[has_inkar_record] <- pca_sensitivity_kaiser$raw
corridor_sf$vuln_index_pca23_kaiser_z[has_inkar_record] <- pca_sensitivity_kaiser$z

eigenvalues <- pca_sensitivity_70pct$eigenvalues
variance_share <- pca_sensitivity_70pct$variance_share

scree_table <- tibble(
  PC = seq_along(eigenvalues),
  Eigenvalue = eigenvalues,
  Variance = variance_share,
  Cumulative = cumsum(variance_share),
  retained_70pct = PC <= n_components_70pct,
  retained_kaiser = Eigenvalue > 1
)

save_table(scree_table, "corridor_scree_table.csv")

scree_plot <- ggplot(scree_table, aes(x = PC, y = Eigenvalue)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Scree Plot: final 23-variable vulnerability PCA sensitivity", x = "PC", y = "Eigenvalue") +
  theme_classic()

cumulative_variance_plot <- ggplot(scree_table, aes(x = PC, y = Cumulative)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  labs(title = "Cumulative variance: final 23-variable PCA sensitivity", x = "PC", y = "Cumulative") +
  theme_classic()

save_plot(scree_plot, "scree_kaiser_corridor.png", width = 8, height = 5, subdir = "plots")
save_plot(cumulative_variance_plot, "cumulative_variance_corridor.png", width = 8, height = 5, subdir = "plots")

loading_table <- as.data.frame(pca_sensitivity_70pct$pca$rotation)
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

pc_alignment_table <- tibble(
  PC = paste0("PC", seq_len(n_components_70pct)),
  alignment = as.numeric(pca_sensitivity_70pct$alignment),
  variance = variance_share[seq_len(n_components_70pct)],
  cumulative = cumsum(variance_share)[seq_len(n_components_70pct)]
)
save_table(pc_alignment_table, "corridor_pca23_component_alignment.csv")

dimension_summary <- main_index$dimension_scores %>%
  summarise(across(everything(), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value")

save_table(dimension_summary, "corridor_vulnerability_dimension_summary.csv")

dimension_correlation <- cor(main_index$dimension_scores, use = "complete.obs") %>%
  as.data.frame() %>%
  rownames_to_column("dimension")

save_table(dimension_correlation, "corridor_vulnerability_dimension_correlations.csv")

index_variable_correlations <- tibble(
  variable = pca_columns,
  corr_with_main_index = sapply(
    pca_columns,
    function(variable_name) cor(
      main_index$z,
      vulnerability_imputed[[variable_name]],
      use = "complete.obs"
    )
  ),
  corr_with_pca23_70pct_index = sapply(
    pca_columns,
    function(variable_name) cor(
      pca_sensitivity_70pct$z,
      vulnerability_imputed[[variable_name]],
      use = "complete.obs"
    )
  )
) %>%
  left_join(vulnerability_variables, by = "variable")

save_table(index_variable_correlations, "corridor_vulnerability_index_variable_correlations.csv")

pca_sensitivity_summary <- tibble(
  comparison = c(
    "main_vs_pca23_70pct_pearson",
    "main_vs_pca23_70pct_spearman",
    "pca23_70pct_vs_kaiser_pearson",
    "pca23_70pct_vs_kaiser_spearman"
  ),
  value = c(
    cor(corridor_sf$vuln_index_main_z, corridor_sf$vuln_index_pca23_70pct_z, use = "complete.obs"),
    cor(corridor_sf$vuln_index_main_z, corridor_sf$vuln_index_pca23_70pct_z, use = "complete.obs", method = "spearman"),
    cor(corridor_sf$vuln_index_pca23_70pct_z, corridor_sf$vuln_index_pca23_kaiser_z, use = "complete.obs"),
    cor(corridor_sf$vuln_index_pca23_70pct_z, corridor_sf$vuln_index_pca23_kaiser_z, use = "complete.obs", method = "spearman")
  )
)

save_table(pca_sensitivity_summary, "corridor_vulnerability_index_sensitivity_summary.csv")

log_message("Main vulnerability index: direction-coded, dimension-balanced 23-variable composite.")
log_message("PCA sensitivity components retained by 70% cumulative variance: ", n_components_70pct)
log_message("PCA sensitivity components retained by Kaiser criterion: ", n_components_kaiser)

index_summary <- corridor_sf %>%
  st_drop_geometry() %>%
  summarise(
    municipalities = n(),
    municipalities_with_vulnerability_index = sum(!is.na(vuln_index_main_z)),
    rp100_mean = mean(flood_share_rp100, na.rm = TRUE),
    rp100_median = median(flood_share_rp100, na.rm = TRUE),
    vuln_mean = mean(vuln_index_main_z, na.rm = TRUE),
    vuln_sd = sd(vuln_index_main_z, na.rm = TRUE),
    corr_vuln_rp100 = cor(vuln_index_main_z, flood_share_rp100, use = "complete.obs"),
    corr_main_pca23_70pct = cor(vuln_index_main_z, vuln_index_pca23_70pct_z, use = "complete.obs"),
    corr_pca23_70pct_pca23_kaiser = cor(vuln_index_pca23_70pct_z, vuln_index_pca23_kaiser_z, use = "complete.obs")
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
    has_inkar_vulnerability_record,
    starts_with("vuln_dim_"),
    starts_with("PC"),
    vuln_index_main,
    vuln_index_main_z,
    vuln_index_pca23_70pct,
    vuln_index_pca23_70pct_z,
    vuln_index_pca23_kaiser,
    vuln_index_pca23_kaiser_z
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
    subtitle = "Direction-coded, dimension-balanced composite from 23 screened INKAR indicators",
    fill = "Index (z)",
    caption = "Source: INKAR (BBSR), EFAS corridor definition. Own processing. One corridor municipality lacks INKAR data."
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
  PC1 = "PCA sensitivity PC1: household and accessibility contrast",
  PC2 = "PCA sensitivity PC2: income and service-access gradient",
  PC3 = "PCA sensitivity PC3: labour and deprivation structure",
  PC4 = "PCA sensitivity PC4: age and household composition"
)

pc_subtitles <- c(
  PC1 = "Direction-coded 23-variable PCA; see loading table for exact signs",
  PC2 = "Direction-coded 23-variable PCA; see loading table for exact signs",
  PC3 = "Direction-coded 23-variable PCA; see loading table for exact signs",
  PC4 = "Direction-coded 23-variable PCA; see loading table for exact signs"
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
