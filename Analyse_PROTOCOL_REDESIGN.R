# ============================================================
# Elbe Flood Justice Analysis - Protocol Redesign
# Author: Codex
# Date: 2026-03-10
# Purpose:
#   Rebuild the analysis around the updated method protocol:
#   - keep BfG HQ100 Zone 1 as main unprotected exposure outcome
#   - replace the old protection proxy in the models with EU NUTS3 protection data
#   - build a cleaner vulnerability index from theory-led domain PCAs
#   - add depth-weighted hazard sensitivity from municipality x zone x depth aggregates
#   - keep spatial diagnostics and publication-ready outputs in one script
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(broom)
  library(spdep)
  library(spatialreg)
  library(tibble)
})

options(scipen = 999)
set.seed(42)

# ============================================================
# 0) Paths
# ============================================================

mun_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/gemeinden_elbe_final_full.gpkg"
elbe_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"
eu_prot_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/EU_Flood_Maps/Flood_protection/peseta4_protection_nuts3.shp"
agg1_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/AGG1_flood_by_zone_depth.gpkg"

out_dir <- "outputs_method_redesign"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "maps"), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) Helpers
# ============================================================

save_plot <- function(p, filename, w = 8, h = 5, subdir = "plots") {
  ggsave(
    filename = file.path(out_dir, subdir, filename),
    plot = p,
    width = w,
    height = h,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

save_table <- function(x, filename) {
  readr::write_csv(x, file.path(out_dir, "tables", filename))
}

save_note <- function(lines, filename) {
  writeLines(lines, con = file.path(out_dir, "tables", filename), useBytes = TRUE)
}

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

z_num <- function(x) {
  as.numeric(scale(x))
}

safe_cor <- function(x, y) {
  suppressWarnings(cor(x, y, use = "complete.obs"))
}

coalesce_numeric <- function(...) {
  xs <- list(...)
  out <- xs[[1]]
  if (length(xs) == 1) return(out)
  for (i in 2:length(xs)) {
    out <- dplyr::coalesce(out, xs[[i]])
  }
  out
}

median_impute_df <- function(df) {
  df %>%
    mutate(across(
      everything(),
      ~ {
        x <- as.numeric(.x)
        if (all(is.na(x))) {
          x
        } else {
          x[is.na(x)] <- median(x, na.rm = TRUE)
          x
        }
      }
    ))
}

map_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 1)
    )
}

map_annotations <- function() {
  if (!requireNamespace("ggspatial", quietly = TRUE)) {
    return(list())
  }

  list(
    ggspatial::annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering,
      height = grid::unit(1.2, "cm"),
      width = grid::unit(1.2, "cm")
    ),
    ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.25,
      line_width = 0.6,
      text_cex = 0.7
    )
  )
}

make_listw <- function(sf_polys, queen = TRUE, snap = 50, k = 4) {
  nb1 <- spdep::poly2nb(sf_polys, queen = queen, snap = snap)
  has_islands <- any(spdep::card(nb1) == 0)
  n_subgraphs <- tryCatch(
    length(unique(spdep::n.comp.nb(nb1)$comp.id)),
    error = function(e) NA_integer_
  )

  if (isTRUE(has_islands) || (!is.na(n_subgraphs) && n_subgraphs > 1)) {
    message("Neighbour warning (islands/subgraphs). Using k-nearest neighbours (k=", k, ").")
    pts <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(sf_polys)))
    nb_knn <- spdep::knn2nb(spdep::knearneigh(pts, k = k))
    return(spdep::nb2listw(nb_knn, style = "W", zero.policy = TRUE))
  }

  spdep::nb2listw(nb1, style = "W", zero.policy = TRUE)
}

run_domain_pca <- function(sf_obj, vars, domain_name) {
  raw_df <- sf_obj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(all_of(vars))

  keep <- names(raw_df)[vapply(
    raw_df,
    function(x) {
      x <- as.numeric(x)
      sum(!is.na(x)) > 30 && isTRUE(sd(x, na.rm = TRUE) > 0)
    },
    logical(1)
  )]

  if (length(keep) == 0) {
    stop(paste("No usable variables left for domain:", domain_name))
  }

  dat_imp <- raw_df %>%
    dplyr::select(all_of(keep)) %>%
    median_impute_df()

  dat_scaled <- scale(dat_imp)
  anchor <- rowMeans(dat_scaled, na.rm = TRUE)

  if (ncol(dat_scaled) == 1) {
    score <- as.numeric(dat_scaled[, 1])
    out <- list(
      score = z_num(score),
      summary = tibble(
        domain = domain_name,
        n_vars = 1,
        method = "single_variable_z",
        pc1_variance = 1,
        kept_variables = paste(keep, collapse = "; ")
      ),
      loadings = tibble(
        domain = domain_name,
        variable = keep,
        loading = 1,
        abs_loading = 1
      )
    )
    return(out)
  }

  pca <- prcomp(dat_scaled)
  score <- pca$x[, 1]
  load <- pca$rotation[, 1]

  if (!is.na(safe_cor(score, anchor)) && safe_cor(score, anchor) < 0) {
    score <- -score
    load <- -load
  }

  list(
    score = z_num(score),
    summary = tibble(
      domain = domain_name,
      n_vars = ncol(dat_scaled),
      method = "PC1",
      pc1_variance = (pca$sdev[1]^2) / sum(pca$sdev^2),
      kept_variables = paste(keep, collapse = "; ")
    ),
    loadings = tibble(
      domain = domain_name,
      variable = names(load),
      loading = as.numeric(load),
      abs_loading = abs(as.numeric(load))
    ) %>%
      arrange(desc(abs_loading))
  )
}

# ============================================================
# 2) Load and harmonise base layers
# ============================================================

MUN <- sf::st_read(mun_path, quiet = TRUE)
ELBE <- sf::st_read(elbe_path, quiet = TRUE)

ELBE <- sf::st_zm(ELBE, drop = TRUE, what = "ZM")
ELBE <- sf::st_make_valid(ELBE)
if (is.na(sf::st_crs(ELBE))) stop("ELBE layer has no CRS.")
if (sf::st_crs(ELBE) != sf::st_crs(MUN)) ELBE <- sf::st_transform(ELBE, sf::st_crs(MUN))
ELBE <- suppressWarnings(sf::st_cast(ELBE, "LINESTRING"))

char_cols <- names(MUN)[vapply(sf::st_drop_geometry(MUN), is.character, logical(1))]
MUN <- MUN %>%
  mutate(across(all_of(char_cols), ~ na_if(str_trim(.x), "")))

year_cols <- names(MUN)[str_detect(names(MUN), "_year$")]
if (length(year_cols) > 0) {
  MUN <- MUN %>%
    mutate(across(
      all_of(year_cols),
      ~ as.integer(readr::parse_number(as.character(.x)))
    ))
}

protected_text_cols <- names(MUN)[str_detect(
  names(MUN),
  regex(
    "AGS|Gemeinde|Kreis|Land|Objekt|Identifikator|Bezeichnung|Name|mun_name|LKZ|ARS|GEN|Geofaktor|BSG|GF|Regierungsbezirk|Verwaltungssitz|layer|path",
    ignore_case = TRUE
  )
)]
parse_cols <- setdiff(char_cols, protected_text_cols)
if (length(parse_cols) > 0) {
  MUN <- MUN %>%
    mutate(across(all_of(parse_cols), readr::parse_number))
}

MUN <- MUN %>%
  mutate(
    AGS_final = coalesce(
      str_pad(as.character(exposure_AGS_final), width = 8, side = "left", pad = "0"),
      str_pad(as.character(exposure_AGS), width = 8, side = "left", pad = "0"),
      str_pad(as.character(Gemeindeschlüssel_AGS), width = 8, side = "left", pad = "0"),
      str_pad(as.character(exposure_Gemeindeschlüssel_AGS_2), width = 8, side = "left", pad = "0")
    ),
    state_code = substr(AGS_final, 1, 2),
    state_factor = factor(state_code),
    mun_area_m2 = coalesce_numeric(exposure_mun_area_m2, mun_area_basin_m2, exposure_mun_area_basin_m2)
  )

state_counts <- table(MUN$state_code)
small_states <- names(state_counts[state_counts < 10])
MUN <- MUN %>%
  mutate(
    state_model = factor(if_else(state_code %in% small_states, "small_state", state_code))
  )

# ============================================================
# 3) Hazard metrics (BfG HQ100)
# ============================================================

req_cols <- c(
  "exposure_share_flood_total",
  "exposure_share_flood_zone1",
  "exposure_share_flood_zone2",
  "exposure_share_flood_zone3"
)
missing_hazard <- setdiff(req_cols, names(MUN))
if (length(missing_hazard) > 0) {
  stop(paste("Missing hazard columns:", paste(missing_hazard, collapse = ", ")))
}

MUN <- MUN %>%
  mutate(
    exposure_flood_area_total_m2 = coalesce_numeric(exposure_flood_area_total_m2, 0),
    exposure_flood_area_zone1_m2 = coalesce_numeric(exposure_flood_area_zone1_m2, 0),
    exposure_flood_area_zone2_m2 = coalesce_numeric(exposure_flood_area_zone2_m2, 0),
    exposure_flood_area_zone3_m2 = coalesce_numeric(exposure_flood_area_zone3_m2, 0),

    exposure_share_flood_total = coalesce_numeric(exposure_share_flood_total, 0),
    exposure_share_flood_zone1 = coalesce_numeric(exposure_share_flood_zone1, 0),
    exposure_share_flood_zone2 = coalesce_numeric(exposure_share_flood_zone2, 0),
    exposure_share_flood_zone3 = coalesce_numeric(exposure_share_flood_zone3, 0),

    share_total_hazard = exposure_share_flood_total,
    share_zone1 = exposure_share_flood_zone1,
    share_zone2 = exposure_share_flood_zone2,
    share_zone3_bfg = exposure_share_flood_zone3,

    risk_zone1_share = share_zone1,
    hazard_any = as.integer(share_total_hazard > 0),
    zone1_any = as.integer(share_zone1 > 0)
  )

hazard_balance <- MUN %>%
  sf::st_drop_geometry() %>%
  transmute(
    diff = share_total_hazard - (share_zone1 + share_zone2 + share_zone3_bfg)
  )

save_table(
  tibble(
    metric = c("min_diff", "median_diff", "mean_diff", "max_diff"),
    value = c(
      min(hazard_balance$diff, na.rm = TRUE),
      median(hazard_balance$diff, na.rm = TRUE),
      mean(hazard_balance$diff, na.rm = TRUE),
      max(hazard_balance$diff, na.rm = TRUE)
    )
  ),
  "hazard_balance_check.csv"
)

# ============================================================
# 4) EU flood protection integration (NUTS3)
# ============================================================

EU_PROT <- sf::st_read(eu_prot_path, quiet = TRUE) %>%
  sf::st_make_valid() %>%
  filter(str_starts(nuts_id, "DE"))

if (sf::st_crs(EU_PROT) != sf::st_crs(MUN)) {
  EU_PROT <- sf::st_transform(EU_PROT, sf::st_crs(MUN))
}

mun_points <- sf::st_as_sf(
  MUN %>%
    sf::st_drop_geometry() %>%
    dplyr::select(AGS_final),
  geometry = sf::st_point_on_surface(sf::st_geometry(MUN)),
  crs = sf::st_crs(MUN)
)

prot_join <- sf::st_join(
  mun_points,
  EU_PROT %>% dplyr::select(nuts_id, eu_protection_years = prtctn_),
  join = sf::st_within,
  left = TRUE
)

missing_prot <- is.na(prot_join$eu_protection_years)
if (any(missing_prot)) {
  nearest_idx <- sf::st_nearest_feature(prot_join[missing_prot, ], EU_PROT)
  prot_join$eu_protection_years[missing_prot] <- EU_PROT$prtctn_[nearest_idx]
  prot_join$nuts_id[missing_prot] <- EU_PROT$nuts_id[nearest_idx]
}

prot_lookup <- sf::st_drop_geometry(prot_join)

MUN$eu_nuts3 <- prot_lookup$nuts_id
MUN$eu_protection_years <- as.numeric(prot_lookup$eu_protection_years)
MUN$eu_protection_years_z <- z_num(MUN$eu_protection_years)
MUN$eu_protection_class <- case_when(
  is.na(MUN$eu_protection_years) ~ "Missing",
  MUN$eu_protection_years < 100 ~ "<100y",
  MUN$eu_protection_years < 150 ~ "100-149y",
  TRUE ~ ">=150y"
)
MUN$eu_protection_class <- factor(
  MUN$eu_protection_class,
  levels = c("<100y", "100-149y", ">=150y", "Missing")
)

# ============================================================
# 5) Depth-weighted sensitivity from municipality x zone x depth
# ============================================================

depth_lookup <- tibble(
  depth_class = 1:5,
  depth_mid_m = c(0.25, 0.75, 1.5, 3.0, 5.0)
)

AGG1 <- sf::st_read(
  agg1_path,
  query = 'SELECT "Gemeindeschlüssel_AGS", Zone, depth_class, flood_area_m2_sum FROM AGG1_flood_by_zone_depth',
  quiet = TRUE
)

AGG1 <- as_tibble(AGG1) %>%
  mutate(
    AGS_final = str_pad(as.character(`Gemeindeschlüssel_AGS`), width = 8, side = "left", pad = "0")
  ) %>%
  left_join(depth_lookup, by = "depth_class") %>%
  mutate(
    known_depth_area_m2 = if_else(!is.na(depth_mid_m), flood_area_m2_sum, 0),
    unknown_depth_area_m2 = if_else(is.na(depth_mid_m), flood_area_m2_sum, 0),
    weighted_depth_area = flood_area_m2_sum * coalesce(depth_mid_m, 0)
  )

depth_by_zone <- AGG1 %>%
  mutate(
    zone_name = case_when(
      Zone == 1 ~ "zone1",
      Zone == 2 ~ "zone2",
      Zone == 3 ~ "zone3",
      TRUE ~ paste0("zone", Zone)
    )
  ) %>%
  group_by(AGS_final, zone_name) %>%
  summarise(
    flooded_area_m2 = sum(flood_area_m2_sum, na.rm = TRUE),
    known_depth_area_m2 = sum(known_depth_area_m2, na.rm = TRUE),
    unknown_depth_area_m2 = sum(unknown_depth_area_m2, na.rm = TRUE),
    weighted_depth_area = sum(weighted_depth_area, na.rm = TRUE),
    mean_depth_m = {
      kd <- sum(known_depth_area_m2, na.rm = TRUE)
      if (kd > 0) sum(weighted_depth_area, na.rm = TRUE) / kd else NA_real_
    },
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = zone_name,
    values_from = c(flooded_area_m2, known_depth_area_m2, unknown_depth_area_m2, weighted_depth_area, mean_depth_m),
    names_glue = "{zone_name}_{.value}"
  )

MUN <- MUN %>%
  left_join(depth_by_zone, by = "AGS_final") %>%
  mutate(
    zone1_depth_weighted_exposure = zone1_weighted_depth_area / mun_area_m2,
    total_depth_weighted_exposure = (
      coalesce_numeric(zone1_weighted_depth_area, 0) +
        coalesce_numeric(zone2_weighted_depth_area, 0) +
        coalesce_numeric(zone3_weighted_depth_area, 0)
    ) / mun_area_m2,
    zone1_unknown_depth_share = zone1_unknown_depth_area_m2 / mun_area_m2,
    total_unknown_depth_share = (
      coalesce_numeric(zone1_unknown_depth_area_m2, 0) +
        coalesce_numeric(zone2_unknown_depth_area_m2, 0) +
        coalesce_numeric(zone3_unknown_depth_area_m2, 0)
    ) / mun_area_m2
  )

depth_diag <- MUN %>%
  sf::st_drop_geometry() %>%
  summarise(
    mean_zone1_unknown_depth_share = mean(zone1_unknown_depth_share, na.rm = TRUE),
    max_zone1_unknown_depth_share = max(zone1_unknown_depth_share, na.rm = TRUE),
    mean_zone1_depth_weighted_exposure = mean(zone1_depth_weighted_exposure, na.rm = TRUE),
    max_zone1_depth_weighted_exposure = max(zone1_depth_weighted_exposure, na.rm = TRUE)
  )

save_table(depth_diag, "depth_metric_diagnostics.csv")

# ============================================================
# 6) Indicator engineering for revised vulnerability modelling
# ============================================================

MUN <- MUN %>%
  mutate(
    v_dep_alg2 = exposure_share_alg2_sgb2,
    v_dep_sgb2_housing = exposure_share_sgb2_with_housing_costs,
    v_dep_longterm_unemp = exposure_share_longterm_unemp,
    v_dep_unemp_u25 = exposure_unemp_u25_per_1000,
    v_dep_purchasing_power = -exposure_purchasing_power,
    v_dep_income_tax = -exposure_income_tax_per_capita,
    v_dep_income_low = exposure_share_hh_income_low,

    v_age_65plus = exposure_share_age_65plus,
    v_age_75plus = exposure_share_age_75plus,
    v_age_old_dep = exposure_old_age_dependency,

    v_house_single_parent = exposure_share_bg_single_parent,
    v_house_single = exposure_share_single_households,

    v_access_supermarket = exposure_dist_supermarket_m,
    v_access_pharmacy = exposure_dist_pharmacy_m,
    v_access_gp = exposure_dist_gp_m,
    v_access_pt = exposure_dist_public_transport_m,
    v_access_doctors = -exposure_doctors_total_per_1000,
    v_access_broadband = -exposure_share_bb_100mbit,

    ctl_density_log = log1p(exposure_pop_density_per_km2),
    ctl_pop_log = log1p(exposure_pop_total),
    ctl_area_log = log1p(mun_area_m2)
  )

indicator_map <- tribble(
  ~analysis_var, ~source_var, ~source_year_var,
  "v_dep_alg2", "exposure_share_alg2_sgb2", "exposure_share_alg2_sgb2_year",
  "v_dep_sgb2_housing", "exposure_share_sgb2_with_housing_costs", "exposure_share_sgb2_with_housing_costs_year",
  "v_dep_longterm_unemp", "exposure_share_longterm_unemp", "exposure_share_longterm_unemp_year",
  "v_dep_unemp_u25", "exposure_unemp_u25_per_1000", "exposure_unemp_u25_per_1000_year",
  "v_dep_purchasing_power", "exposure_purchasing_power", "exposure_purchasing_power_year",
  "v_dep_income_tax", "exposure_income_tax_per_capita", "exposure_income_tax_per_capita_year",
  "v_dep_income_low", "exposure_share_hh_income_low", "exposure_share_hh_income_low_year",
  "v_age_65plus", "exposure_share_age_65plus", "exposure_share_age_65plus_year",
  "v_age_75plus", "exposure_share_age_75plus", "exposure_share_age_75plus_year",
  "v_age_old_dep", "exposure_old_age_dependency", "exposure_old_age_dependency_year",
  "v_house_single_parent", "exposure_share_bg_single_parent", "exposure_share_bg_single_parent_year",
  "v_house_single", "exposure_share_single_households", "exposure_share_single_households_year",
  "v_access_supermarket", "exposure_dist_supermarket_m", "exposure_dist_supermarket_m_year",
  "v_access_pharmacy", "exposure_dist_pharmacy_m", "exposure_dist_pharmacy_m_year",
  "v_access_gp", "exposure_dist_gp_m", "exposure_dist_gp_m_year",
  "v_access_pt", "exposure_dist_public_transport_m", "exposure_dist_public_transport_m_year",
  "v_access_doctors", "exposure_doctors_total_per_1000", "exposure_doctors_total_per_1000_year",
  "v_access_broadband", "exposure_share_bb_100mbit", "exposure_share_bb_100mbit_year"
)

year_audit <- indicator_map %>%
  rowwise() %>%
  mutate(
    distinct_years = length(unique(stats::na.omit(MUN[[source_year_var]]))),
    min_year = suppressWarnings(min(MUN[[source_year_var]], na.rm = TRUE)),
    max_year = suppressWarnings(max(MUN[[source_year_var]], na.rm = TRUE))
  ) %>%
  ungroup()

save_table(year_audit, "indicator_year_audit.csv")

# ============================================================
# 7) Domain PCA vulnerability index
# ============================================================

domain_specs <- list(
  deprivation = c(
    "v_dep_alg2",
    "v_dep_sgb2_housing",
    "v_dep_longterm_unemp",
    "v_dep_unemp_u25",
    "v_dep_purchasing_power",
    "v_dep_income_tax",
    "v_dep_income_low"
  ),
  age = c(
    "v_age_65plus",
    "v_age_75plus",
    "v_age_old_dep"
  ),
  household = c(
    "v_house_single_parent",
    "v_house_single"
  ),
  access = c(
    "v_access_supermarket",
    "v_access_pharmacy",
    "v_access_gp",
    "v_access_pt",
    "v_access_doctors",
    "v_access_broadband"
  )
)

domain_results <- lapply(names(domain_specs), function(dom) {
  run_domain_pca(MUN, domain_specs[[dom]], dom)
})
names(domain_results) <- names(domain_specs)

for (dom in names(domain_results)) {
  MUN[[paste0("dom_", dom, "_z")]] <- domain_results[[dom]]$score
}

domain_summary <- bind_rows(lapply(domain_results, `[[`, "summary"))
domain_loadings <- bind_rows(lapply(domain_results, `[[`, "loadings"))

save_table(domain_summary, "domain_pca_summary.csv")
save_table(domain_loadings, "domain_pca_loadings.csv")

domain_cols <- paste0("dom_", names(domain_results), "_z")

MUN <- MUN %>%
  mutate(
    vuln_index_domain_raw = rowMeans(across(all_of(domain_cols)), na.rm = TRUE),
    vuln_index_domain_z = z_num(vuln_index_domain_raw),
    vuln_index_domain_01 = rescale01(vuln_index_domain_z)
  )

# ============================================================
# 8) Full PCA sensitivity with selected relevant components
# ============================================================

candidate_vars <- unique(unlist(domain_specs))
candidate_df <- MUN %>%
  sf::st_drop_geometry() %>%
  dplyr::select(all_of(candidate_vars)) %>%
  median_impute_df()

candidate_scaled <- scale(candidate_df)
full_pca <- prcomp(candidate_scaled)

anchor_components <- intersect(
  c("v_dep_alg2", "v_dep_longterm_unemp", "v_dep_purchasing_power", "v_dep_income_tax"),
  colnames(candidate_scaled)
)
anchor_score <- rowMeans(candidate_scaled[, anchor_components, drop = FALSE], na.rm = TRUE)

eig <- full_pca$sdev^2
var_prop <- eig / sum(eig)
pc_tbl <- tibble(
  PC = paste0("PC", seq_along(eig)),
  eigenvalue = eig,
  variance = var_prop,
  cumulative = cumsum(var_prop),
  anchor_cor = vapply(seq_along(eig), function(i) safe_cor(full_pca$x[, i], anchor_score), numeric(1))
) %>%
  mutate(
    retain_kaiser = eigenvalue > 1,
    retain_relevant = retain_kaiser & coalesce(abs(anchor_cor) >= 0.25, FALSE),
    sign = if_else(!is.na(anchor_cor) & anchor_cor < 0, -1, 1)
  )

selected_pc_ids <- which(pc_tbl$retain_relevant)
if (length(selected_pc_ids) == 0) {
  selected_pc_ids <- which(pc_tbl$retain_kaiser)
}
if (length(selected_pc_ids) == 0) {
  selected_pc_ids <- 1
}

selected_scores <- full_pca$x[, selected_pc_ids, drop = FALSE]
selected_signs <- pc_tbl$sign[selected_pc_ids]
selected_scores <- sweep(selected_scores, 2, selected_signs, `*`)

selected_weights <- pc_tbl$variance[selected_pc_ids]
selected_weights <- selected_weights / sum(selected_weights)

MUN$vuln_index_selected_pca_raw <- as.numeric(selected_scores %*% selected_weights)
MUN$vuln_index_selected_pca_z <- z_num(MUN$vuln_index_selected_pca_raw)
MUN$vuln_index_selected_pca_01 <- rescale01(MUN$vuln_index_selected_pca_z)

selected_loadings <- as.data.frame(full_pca$rotation[, selected_pc_ids, drop = FALSE])
for (i in seq_along(selected_pc_ids)) {
  selected_loadings[, i] <- selected_loadings[, i] * selected_signs[i]
}
selected_loadings$variable <- rownames(selected_loadings)
selected_loadings <- selected_loadings %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") %>%
  mutate(abs_loading = abs(loading)) %>%
  group_by(PC) %>%
  arrange(desc(abs_loading)) %>%
  slice_head(n = 10) %>%
  ungroup()

save_table(pc_tbl, "selected_pca_component_table.csv")
save_table(selected_loadings, "selected_pca_top_loadings.csv")

# ============================================================
# 9) Validation and descriptive justice tables
# ============================================================

validation_tbl <- tibble(
  index = c("domain_main", "selected_pca"),
  cor_with_alg2 = c(
    safe_cor(MUN$vuln_index_domain_z, MUN$exposure_share_alg2_sgb2),
    safe_cor(MUN$vuln_index_selected_pca_z, MUN$exposure_share_alg2_sgb2)
  ),
  cor_with_longterm_unemp = c(
    safe_cor(MUN$vuln_index_domain_z, MUN$exposure_share_longterm_unemp),
    safe_cor(MUN$vuln_index_selected_pca_z, MUN$exposure_share_longterm_unemp)
  ),
  cor_with_zone1_share = c(
    safe_cor(MUN$vuln_index_domain_z, MUN$risk_zone1_share),
    safe_cor(MUN$vuln_index_selected_pca_z, MUN$risk_zone1_share)
  )
)

save_table(validation_tbl, "vulnerability_validation.csv")

descriptive_df <- MUN %>%
  sf::st_drop_geometry() %>%
  mutate(
    vuln_q = ntile(vuln_index_domain_z, 5)
  )

by_vuln_q <- descriptive_df %>%
  group_by(vuln_q) %>%
  summarise(
    n = n(),
    mean_zone1_share = mean(risk_zone1_share, na.rm = TRUE),
    median_zone1_share = median(risk_zone1_share, na.rm = TRUE),
    mean_total_share = mean(share_total_hazard, na.rm = TRUE),
    mean_zone1_depth_weighted = mean(zone1_depth_weighted_exposure, na.rm = TRUE),
    mean_eu_protection_years = mean(eu_protection_years, na.rm = TRUE),
    .groups = "drop"
  )

by_protection_class <- descriptive_df %>%
  group_by(eu_protection_class) %>%
  summarise(
    n = n(),
    mean_vulnerability = mean(vuln_index_domain_z, na.rm = TRUE),
    mean_zone1_share = mean(risk_zone1_share, na.rm = TRUE),
    mean_zone1_depth_weighted = mean(zone1_depth_weighted_exposure, na.rm = TRUE),
    .groups = "drop"
  )

save_table(by_vuln_q, "justice_summary_by_vulnerability_quintile.csv")
save_table(by_protection_class, "justice_summary_by_eu_protection_class.csv")

bfg_vs_eu <- descriptive_df %>%
  summarise(
    cor_bfg_zone3_vs_eu_protection = safe_cor(share_zone3_bfg, eu_protection_years),
    mean_bfg_zone3_share = mean(share_zone3_bfg, na.rm = TRUE),
    mean_eu_protection_years = mean(eu_protection_years, na.rm = TRUE)
  )
save_table(bfg_vs_eu, "bfg_zone3_vs_eu_protection.csv")

# ============================================================
# 10) Two-part models + spatial sensitivity
# ============================================================

analysis_df <- descriptive_df %>%
  transmute(
    AGS_final,
    state_model,
    risk_zone1_share,
    zone1_depth_weighted_exposure,
    exposed = as.integer(risk_zone1_share > 0),
    vuln = vuln_index_domain_z,
    vuln_pca = vuln_index_selected_pca_z,
    eu_protection = eu_protection_years_z,
    density = z_num(ctl_density_log),
    population = z_num(ctl_pop_log),
    mun_area = z_num(ctl_area_log)
  ) %>%
  filter(
    !is.na(vuln),
    !is.na(eu_protection),
    !is.na(density),
    !is.na(population),
    !is.na(mun_area)
  ) %>%
  mutate(
    log_zone1_share = log(pmax(risk_zone1_share, 0.000001)),
    log_zone1_depth_weighted = log(pmax(zone1_depth_weighted_exposure, 0.000001))
  )

m_presence <- glm(
  exposed ~ vuln + eu_protection + density + population + mun_area + state_model,
  data = analysis_df,
  family = binomial()
)

m_amount <- lm(
  log_zone1_share ~ vuln + eu_protection + density + population + mun_area + state_model,
  data = analysis_df %>% filter(exposed == 1)
)

m_depth <- lm(
  log_zone1_depth_weighted ~ vuln + eu_protection + density + population + mun_area + state_model,
  data = analysis_df %>% filter(zone1_depth_weighted_exposure > 0)
)

save_table(broom::tidy(m_presence, conf.int = TRUE), "model_presence_logit_domain.csv")
save_table(broom::glance(m_presence), "model_presence_logit_domain_fit.csv")
save_table(broom::tidy(m_amount, conf.int = TRUE), "model_amount_lm_domain.csv")
save_table(broom::glance(m_amount), "model_amount_lm_domain_fit.csv")
save_table(broom::tidy(m_depth, conf.int = TRUE), "model_depth_lm_domain.csv")
save_table(broom::glance(m_depth), "model_depth_lm_domain_fit.csv")

spatial_sf <- MUN %>%
  mutate(
    vuln = vuln_index_domain_z,
    eu_protection = eu_protection_years_z,
    density = z_num(ctl_density_log),
    population = z_num(ctl_pop_log),
    mun_area_ctl = z_num(ctl_area_log),
    risk = risk_zone1_share
  ) %>%
  filter(
    !is.na(vuln),
    !is.na(eu_protection),
    !is.na(density),
    !is.na(population),
    !is.na(mun_area_ctl),
    !is.na(risk)
  )

lw <- make_listw(spatial_sf, queen = TRUE, snap = 50, k = 4)

m_spatial_ols <- lm(
  risk ~ vuln + eu_protection + density + population + mun_area_ctl,
  data = sf::st_drop_geometry(spatial_sf)
)

moran_res <- spdep::moran.test(resid(m_spatial_ols), lw, zero.policy = TRUE)
moran_tbl <- tibble(
  moran_i = unname(moran_res$estimate[[1]]),
  expectation = unname(moran_res$estimate[[2]]),
  variance = unname(moran_res$estimate[[3]]),
  p_value = moran_res$p.value
)
save_table(moran_tbl, "spatial_moran_residuals.csv")

sar <- tryCatch(
  spatialreg::lagsarlm(
    risk ~ vuln + eu_protection + density + population + mun_area_ctl,
    data = spatial_sf,
    listw = lw,
    zero.policy = TRUE
  ),
  error = function(e) NULL
)

sem <- tryCatch(
  spatialreg::errorsarlm(
    risk ~ vuln + eu_protection + density + population + mun_area_ctl,
    data = spatial_sf,
    listw = lw,
    zero.policy = TRUE
  ),
  error = function(e) NULL
)

aic_tbl <- tibble(model = "OLS", AIC = AIC(m_spatial_ols))
if (!is.null(sar)) aic_tbl <- bind_rows(aic_tbl, tibble(model = "SAR_lag", AIC = AIC(sar)))
if (!is.null(sem)) aic_tbl <- bind_rows(aic_tbl, tibble(model = "SEM_error", AIC = AIC(sem)))
save_table(aic_tbl, "model_aic_comparison_redesign.csv")

# ============================================================
# 11) LISA clusters
# ============================================================

x <- spatial_sf$risk
x_lag <- spdep::lag.listw(lw, x, zero.policy = TRUE)
lisa <- spdep::localmoran(x, lw, zero.policy = TRUE)

p_col <- grep("^Pr\\(", colnames(lisa), value = TRUE)
if (length(p_col) == 0) p_col <- grep("^Pr", colnames(lisa), value = TRUE)
if (length(p_col) == 0) stop("Could not find p-value column in localmoran output.")

spatial_sf$lisa_I <- lisa[, "Ii"]
spatial_sf$lisa_p <- lisa[, p_col[1]]

x_c <- x - mean(x, na.rm = TRUE)
lag_c <- x_lag - mean(x_lag, na.rm = TRUE)

spatial_sf$lisa_class <- "Not significant"
sig <- !is.na(spatial_sf$lisa_p) & spatial_sf$lisa_p < 0.05
spatial_sf$lisa_class[sig & x_c > 0 & lag_c > 0] <- "High-High"
spatial_sf$lisa_class[sig & x_c < 0 & lag_c < 0] <- "Low-Low"
spatial_sf$lisa_class[sig & x_c > 0 & lag_c < 0] <- "High-Low"
spatial_sf$lisa_class[sig & x_c < 0 & lag_c > 0] <- "Low-High"
spatial_sf$lisa_class <- factor(
  spatial_sf$lisa_class,
  levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
)

save_table(
  spatial_sf %>%
    sf::st_drop_geometry() %>%
    dplyr::select(AGS_final, risk, lisa_I, lisa_p, lisa_class),
  "lisa_zone1_clusters_redesign.csv"
)

# ============================================================
# 12) Plots and maps
# ============================================================

scree_plot <- ggplot(pc_tbl, aes(x = seq_along(PC), y = eigenvalue)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    title = "Selected-PCA scree plot",
    x = "Component",
    y = "Eigenvalue"
  ) +
  theme_classic()

anchor_plot <- ggplot(pc_tbl, aes(x = seq_along(PC), y = anchor_cor, color = retain_relevant)) +
  geom_hline(yintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("FALSE" = "grey55", "TRUE" = "#b2182b")) +
  labs(
    title = "Component relevance to deprivation anchor",
    x = "Component",
    y = "Correlation with deprivation anchor",
    color = "Selected"
  ) +
  theme_classic()

save_plot(scree_plot, "selected_pca_scree.png")
save_plot(anchor_plot, "selected_pca_anchor_relevance.png")

quintile_plot <- descriptive_df %>%
  ggplot(aes(x = factor(vuln_q), y = risk_zone1_share)) +
  geom_boxplot(fill = "#d9e6f2", color = "#254b6d", outlier.alpha = 0.25) +
  labs(
    title = "Zone 1 exposure by vulnerability quintile",
    x = "Vulnerability quintile",
    y = "Zone 1 share"
  ) +
  theme_classic()

save_plot(quintile_plot, "zone1_by_vulnerability_quintile.png")

map_vuln <- ggplot(MUN) +
  geom_sf(aes(fill = vuln_index_domain_z), color = NA) +
  geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  scale_fill_viridis_c(option = "C") +
  labs(
    title = "Revised vulnerability index",
    subtitle = "Theory-led domain PCA: deprivation, age, household, access",
    fill = "Index (z)",
    caption = "Source: INKAR, own processing."
  ) +
  map_theme()

map_risk <- ggplot(MUN) +
  geom_sf(aes(fill = risk_zone1_share), color = NA) +
  geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  scale_fill_viridis_c(option = "C") +
  labs(
    title = "Unprotected flood exposure (HQ100, Zone 1)",
    subtitle = "Share of municipal area in official floodplain",
    fill = "Share",
    caption = "Source: BfG flood hazard data. Own processing."
  ) +
  map_theme()

map_depth <- ggplot(MUN) +
  geom_sf(aes(fill = zone1_depth_weighted_exposure), color = NA) +
  geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  scale_fill_viridis_c(option = "C") +
  labs(
    title = "Depth-weighted Zone 1 exposure",
    subtitle = "Weighted by known depth classes 1-5; classes 6-7 excluded from depth score",
    fill = "Depth-weighted",
    caption = "Source: BfG AGG1 zone-depth aggregate. Own processing."
  ) +
  map_theme()

map_protection <- ggplot(MUN) +
  geom_sf(aes(fill = eu_protection_class), color = NA) +
  geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  scale_fill_manual(
    values = c(
      "<100y" = "#d73027",
      "100-149y" = "#fdae61",
      ">=150y" = "#1a9850",
      "Missing" = "grey80"
    ),
    drop = FALSE
  ) +
  labs(
    title = "EU flood protection context",
    subtitle = "NUTS3 protection standard joined to municipalities by point-on-surface",
    fill = "Protection",
    caption = "Source: EU Flood Maps / PESETA4. Own processing."
  ) +
  map_theme()

lisa_map <- ggplot(spatial_sf) +
  geom_sf(aes(fill = lisa_class), color = NA) +
  geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  scale_fill_manual(
    values = c(
      "High-High" = "#d73027",
      "Low-Low" = "#4575b4",
      "High-Low" = "#fdae61",
      "Low-High" = "#74add1",
      "Not significant" = "grey85"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Local Moran clusters of Zone 1 exposure",
    subtitle = "LISA on municipal Zone 1 share",
    fill = "Cluster",
    caption = "Source: BfG HQ100 Zone 1. Own processing."
  ) +
  map_theme()

save_plot(map_vuln, "map_vulnerability_index_redesign.png", w = 9, h = 7, subdir = "maps")
save_plot(map_risk, "map_zone1_exposure_redesign.png", w = 9, h = 7, subdir = "maps")
save_plot(map_depth, "map_zone1_depth_weighted_exposure.png", w = 9, h = 7, subdir = "maps")
save_plot(map_protection, "map_eu_protection_class.png", w = 9, h = 7, subdir = "maps")
save_plot(lisa_map, "map_lisa_zone1_clusters_redesign.png", w = 9, h = 7, subdir = "maps")

for (dom in names(domain_results)) {
  dom_col <- paste0("dom_", dom, "_z")

  p <- ggplot(MUN) +
    geom_sf(aes(fill = .data[[dom_col]]), color = NA) +
    geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
    coord_sf(crs = 25832) +
    map_annotations() +
    scale_fill_viridis_c(option = "C") +
    labs(
      title = paste("Domain score:", str_to_title(dom)),
      subtitle = paste("PC1-based domain index for", dom),
      fill = "Score (z)",
      caption = "Source: INKAR, own processing."
    ) +
    map_theme()

  save_plot(p, paste0("map_domain_", dom, ".png"), w = 9, h = 7, subdir = "maps")
}

# ============================================================
# 13) Optional land-use detection
# ============================================================

landuse_candidates <- list.files(
  "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS",
  recursive = TRUE,
  full.names = TRUE,
  pattern = "landuse|land_use|corine|clc|urban.?atlas|nutz",
  ignore.case = TRUE
)

if (length(landuse_candidates) == 0) {
  save_note(
    c(
      "No land-use dataset detected under DATEN + GIS.",
      "Land-use flood overlay was skipped in this redesign script."
    ),
    "landuse_status.txt"
  )
} else {
  save_note(
    c(
      "Land-use datasets detected, but not integrated automatically because class fields are unknown.",
      "Detected paths:",
      landuse_candidates
    ),
    "landuse_status.txt"
  )
}

# ============================================================
# 14) Save enriched dataset + method notes
# ============================================================

saveRDS(MUN, file.path(out_dir, "FULL_protocol_redesign.rds"))

save_note(
  c(
    "Main vulnerability index: vuln_index_domain_z",
    "Sensitivity vulnerability index: vuln_index_selected_pca_z",
    "Main exposure outcome: risk_zone1_share (BfG HQ100 Zone 1 share)",
    "Protection context in models: EU NUTS3 protection standard (eu_protection_years)",
    "GISD is intentionally excluded from this script and reserved for a later validation/sensitivity stage.",
    "BfG Zone 3 retained only as diagnostic because of known inconsistency.",
    "Depth-weighted exposure excludes depth classes 6 and 7 from the depth score and reports them separately as unknown-depth share.",
    "Two-part models are the primary non-spatial inference layer because Zone 1 share is strongly zero-inflated.",
    "OLS + Moran + SAR/SEM are retained as spatial sensitivity checks."
  ),
  "method_notes.txt"
)

message("Protocol redesign analysis finished. Outputs written to: ", out_dir)
