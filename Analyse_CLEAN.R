# ============================================================
# Elbe Flood Risk & Justice Analysis — CLEAN PIPELINE
# Author: Maxi Ahl
# Date: 2026-02-19
# Purpose: Reproducible PCA-based vulnerability index + regression + spatial diagnostics
# ============================================================

# ---- 0) Libraries ----
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
  library(ggspatial)
})

options(scipen = 999)
set.seed(42)

# ============================================================
# 1) Paths + Output helpers  (MUST be defined BEFORE first use)
# ============================================================

gpkg_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/gemeinden_elbe_final_full.gpkg"
elbe_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"

out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "maps"), showWarnings = FALSE, recursive = TRUE)

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

map_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),

      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),

      legend.position  = "right",
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA),

      plot.title    = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.caption  = element_text(size = base_size - 3, hjust = 1)
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

# ============================================================
# 2) Load data (FULL + ELBE)
# ============================================================

FULL <- sf::st_read(gpkg_path, quiet = TRUE)

ELBE <- sf::st_read(elbe_path, quiet = TRUE)
ELBE <- sf::st_zm(ELBE, drop = TRUE, what = "ZM")
ELBE <- sf::st_make_valid(ELBE)

if (is.na(sf::st_crs(ELBE))) stop("ELBE layer has no CRS.")
if (sf::st_crs(ELBE) != sf::st_crs(FULL)) ELBE <- sf::st_transform(ELBE, sf::st_crs(FULL))
ELBE <- suppressWarnings(sf::st_cast(ELBE, "LINESTRING"))

# ============================================================
# 3) Define hazard outcomes (HQ100 zones)
# ============================================================

req_cols <- c(
  "exposure_share_flood_total",
  "exposure_share_flood_zone1",
  "exposure_share_flood_zone2",
  "exposure_share_flood_zone3",
  "share_protected"
)
missing_cols <- setdiff(req_cols, names(FULL))
if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse=", ")))

FULL <- FULL %>%
  mutate(
    exposure_share_flood_total = if_else(is.na(exposure_share_flood_total), 0, exposure_share_flood_total),
    exposure_share_flood_zone1 = if_else(is.na(exposure_share_flood_zone1), 0, exposure_share_flood_zone1),
    exposure_share_flood_zone2 = if_else(is.na(exposure_share_flood_zone2), 0, exposure_share_flood_zone2),
    exposure_share_flood_zone3 = if_else(is.na(exposure_share_flood_zone3), 0, exposure_share_flood_zone3)
  ) %>%
  mutate(
    share_total_hazard      = exposure_share_flood_total,
    share_unprotected       = exposure_share_flood_zone1,
    share_informational     = exposure_share_flood_zone2,
    share_protected_hazard  = exposure_share_flood_zone3,

    # LOCKED definition (your supervisor-consistent choice):
    residual_risk = share_unprotected,

    protection_share_of_hazard = if_else(
      share_total_hazard > 0,
      share_protected_hazard / share_total_hazard,
      NA_real_
    )
  )

check_diff <- FULL %>%
  sf::st_drop_geometry() %>%
  transmute(diff = share_total_hazard - (share_unprotected + share_informational + share_protected_hazard))
print(summary(check_diff$diff))

# ============================================================
# 4) Wide PCA inputs
# ============================================================

wide_pca_vars <- FULL %>%
  sf::st_drop_geometry() %>%
  dplyr::select(
    exposure_share_alg2_sgb2,
    exposure_share_bg_single_parent,
    exposure_share_bg_5plus,
    exposure_share_bg_with_children,
    exposure_share_sgb2_with_housing_costs,

    exposure_share_longterm_unemp,
    exposure_unemp_u25_per_1000,
    exposure_unemp_55plus_per_1000,

    exposure_income_tax_per_capita,
    exposure_purchasing_power,
    exposure_trade_tax_per_capita,
    exposure_tax_revenue_total,

    exposure_share_age_0_3,
    exposure_share_age_3_6,
    exposure_share_age_6_18,
    exposure_share_age_18_25,
    exposure_share_age_25_30,
    exposure_share_age_30_50,
    exposure_share_age_50_65,
    exposure_share_age_65_75,
    exposure_share_age_65plus,
    exposure_share_age_75plus,

    exposure_natural_pop_change,
    exposure_migration_balance,

    exposure_old_age_dependency,
    exposure_youth_dependency,
    exposure_young_old_ratio,

    exposure_share_single_households,
    exposure_share_households_with_children,

    exposure_share_hh_income_high,
    exposure_share_hh_income_medium,
    exposure_share_hh_income_low,
    exposure_students_total_per_1000,
    exposure_students_18_25_per_1000,
    exposure_students_fh_per_1000,

    exposure_gp_general_per_1000,
    exposure_gp_primary_per_1000,
    exposure_internists_per_1000,
    exposure_pediatricians_per_1000_children,
    exposure_doctors_total_per_1000,

    exposure_share_bb_1000mbit,
    exposure_share_bb_100mbit,
    exposure_share_bb_50mbit,
    exposure_share_4g,

    exposure_dist_supermarket_m,
    exposure_dist_pharmacy_m,
    exposure_dist_gp_m,
    exposure_dist_public_transport_m,
    exposure_dist_primary_school_m,

    exposure_pop_density_per_km2,
    exposure_employment_density_per_km2
  ) %>%
  mutate(across(everything(), as.numeric)) %>%
  rename_with(~ stringr::str_remove(.x, "^exposure_"))

# ============================================================
# 5) PCA (median imputation)
# ============================================================

impute_median <- function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }

wide_pca_imp <- wide_pca_vars %>%
  mutate(across(everything(), impute_median))

wide_pca_scaled <- scale(wide_pca_imp)
wide_pca <- prcomp(wide_pca_scaled)

pc_scores <- as.data.frame(wide_pca$x)
names(pc_scores) <- paste0("PC", seq_len(ncol(pc_scores)))
FULL <- bind_cols(FULL, pc_scores)

pca_var <- wide_pca$sdev^2
pca_var_prop <- pca_var / sum(pca_var)

scree_df <- tibble(
  PC = seq_along(pca_var),
  Eigenvalue = pca_var,
  Variance = pca_var_prop,
  Cumulative = cumsum(pca_var_prop)
)

# ============================================================
# 6) Vulnerability index (PC-weighted) + direction lock
# ============================================================

build_index <- function(dat, k, scree_tbl) {
  w <- scree_tbl$Variance[1:k]
  names(w) <- paste0("PC", 1:k)

  pc_mat <- dat %>%
    sf::st_drop_geometry() %>%
    dplyr::select(all_of(names(w))) %>%
    as.matrix()

  idx <- as.numeric(pc_mat %*% w)
  z <- as.numeric(scale(idx))
  list(raw = idx, z = z)
}

idx8 <- build_index(FULL, k = 8, scree_tbl = scree_df)
FULL$vuln_index_main <- idx8$raw
FULL$vuln_index_main_z <- idx8$z

idx12 <- build_index(FULL, k = 12, scree_tbl = scree_df)
FULL$vuln_index_sens12 <- idx12$raw
FULL$vuln_index_sens12_z <- idx12$z

# FIX: correct anchor column name
anchor_cor <- suppressWarnings(cor(FULL$vuln_index_main_z, FULL$exposure_share_alg2_sgb2, use = "complete.obs"))
if (!is.na(anchor_cor) && anchor_cor < 0) {
  FULL$vuln_index_main_z <- -FULL$vuln_index_main_z
  FULL$vuln_index_main <- -FULL$vuln_index_main
  FULL$vuln_index_sens12_z <- -FULL$vuln_index_sens12_z
  FULL$vuln_index_sens12 <- -FULL$vuln_index_sens12
  message("Index flipped: higher = more vulnerable.")
} else {
  message("Index direction OK.")
}

# ============================================================
# 7) Regression block
# ============================================================

df <- FULL %>%
  sf::st_drop_geometry() %>%
  transmute(
    vuln = vuln_index_main_z,
    prot = share_protected,
    risk = residual_risk
  ) %>%
  filter(!is.na(vuln), !is.na(prot), !is.na(risk))

m_ols <- lm(risk ~ vuln * prot, data = df)
ols_tidy <- broom::tidy(m_ols, conf.int = TRUE)
ols_fit  <- broom::glance(m_ols)

print(summary(m_ols))
print(ols_tidy)

# ============================================================
# 8) Spatial diagnostics
# ============================================================

sf_dat <- FULL %>%
  mutate(vuln = vuln_index_main_z, prot = share_protected, risk = residual_risk) %>%
  filter(!is.na(vuln), !is.na(prot), !is.na(risk))

# ---- weights: try contiguity, fall back to k-nearest neighbours if needed ----
make_listw <- function(sf_polys, queen = TRUE, snap = 50, k = 4) {
  nb1 <- spdep::poly2nb(sf_polys, queen = queen, snap = snap)
  has_islands <- any(spdep::card(nb1) == 0)
  n_subgraphs <- tryCatch({
    length(spdep::n.comp.nb(nb1)$comp.id)
  }, error = function(e) NA_integer_)

  if (isTRUE(has_islands) || (!is.na(n_subgraphs) && n_subgraphs > 1)) {
    message("Neighbour warning (islands/subgraphs). Using k-nearest neighbours (k=", k, ") on centroids.")
    ctr <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(sf_polys)))
    nb_knn <- spdep::knn2nb(spdep::knearneigh(ctr, k = k))
    return(spdep::nb2listw(nb_knn, style = "W", zero.policy = TRUE))
  }

  spdep::nb2listw(nb1, style = "W", zero.policy = TRUE)
}

lw <- make_listw(sf_dat, queen = TRUE, snap = 50, k = 4)

e <- resid(lm(risk ~ vuln * prot, data = sf::st_drop_geometry(sf_dat)))
moran_res <- spdep::moran.test(e, lw, zero.policy = TRUE)
print(moran_res)

sar <- spatialreg::lagsarlm(risk ~ vuln * prot, data = sf_dat, listw = lw, zero.policy = TRUE)
sem <- spatialreg::errorsarlm(risk ~ vuln * prot, data = sf_dat, listw = lw, zero.policy = TRUE)
print(AIC(m_ols, sar, sem))

# ============================================================
# 9) OUTPUTS: tables + plots + maps
# ============================================================

# Scree plots
scree_plot <- ggplot(scree_df, aes(x = PC, y = Eigenvalue)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Scree Plot (Kaiser)", x = "PC", y = "Eigenvalue") +
  theme_classic()

cumvar_plot <- ggplot(scree_df, aes(x = PC, y = Cumulative)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  labs(title = "Cumulative Variance", x = "PC", y = "Cumulative") +
  theme_classic()

save_plot(scree_plot, "scree_kaiser.png")
save_plot(cumvar_plot, "cumulative_variance.png")
save_table(scree_df, "scree_table.csv")

# Top loadings
loadings <- as.data.frame(wide_pca$rotation)
loadings$variable <- rownames(loadings)

top_loadings <- loadings %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") %>%
  mutate(abs_loading = abs(loading)) %>%
  group_by(PC) %>%
  arrange(desc(abs_loading)) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  arrange(as.integer(str_remove(PC, "PC")), desc(abs_loading))

save_table(top_loadings, "pca_top_loadings_top8_per_pc.csv")

# Regression exports
save_table(ols_tidy, "regression_ols_interaction_tidy.csv")
save_table(as.data.frame(ols_fit), "regression_ols_interaction_fit.csv")

#
# Maps (with Elbe overlay)
map_vuln <- ggplot(FULL) +
  geom_sf(aes(fill = vuln_index_main_z), color = NA) +
  geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  labs(
    title = "Socio-economic vulnerability index (PCA-weighted)",
    subtitle = "Composite index based on income, unemployment, demography, accessibility",
    fill = "Index (z)",
    caption = "Source: BfG, BBSR (INKAR), RKI (GISD). Own processing."
  ) +
  scale_fill_viridis_c(option = "C") +
  map_theme()

map_risk <- ggplot(FULL) +
  geom_sf(aes(fill = residual_risk), color = NA) +
  geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  labs(
    title = "Unprotected flood exposure (HQ100)",
    subtitle = "Share of municipality area in official floodplain (Zone 1)",
    fill = "Share",
    caption = "Source: BfG flood hazard data. Own processing."
  ) +
  scale_fill_viridis_c(option = "C") +
  map_theme()

save_plot(map_vuln, "map_vulnerability_index.png", w = 9, h = 7, subdir = "maps")
save_plot(map_risk, "map_residual_risk.png", w = 9, h = 7, subdir = "maps")

for (pc in paste0("PC", 1:4)) {

  pc_titles <- c(
    PC1 = "PC1: Socio-economic disadvantage & welfare dependency",
    PC2 = "PC2: Urbanisation, density & accessibility",
    PC3 = "PC3: Demographic ageing & dependency",
    PC4 = "PC4: Education, students & human capital"
  )

  pc_subtitles <- c(
    PC1 = "High loadings: ALG II, low income, unemployment, single households",
    PC2 = "High loadings: population density, transport access, broadband",
    PC3 = "High loadings: 65+, old-age dependency, low youth share",
    PC4 = "High loadings: students, higher income, education indicators"
  )

  p <- ggplot(FULL) +
    geom_sf(aes(fill = .data[[pc]]), color = NA) +
    geom_sf(data = ELBE, inherit.aes = FALSE, color = "black", linewidth = 0.3) +
    coord_sf(crs = 25832) +
    map_annotations() +
    labs(
      title = pc_titles[[pc]],
      subtitle = pc_subtitles[[pc]],
      fill = pc,
      caption = "Source: INKAR (BBSR). PCA-based component scores. Own processing."
    ) +
    scale_fill_viridis_c(option = "C") +
    map_theme()

  save_plot(p, paste0("map_", pc, ".png"), w = 9, h = 7, subdir = "maps")
}

# LISA clusters
x <- sf_dat$risk
x_lag <- spdep::lag.listw(lw, x, zero.policy = TRUE)
lisa <- spdep::localmoran(x, lw, zero.policy = TRUE)

p_col <- grep("^Pr\\(", colnames(lisa), value = TRUE)
if (length(p_col) == 0) p_col <- grep("^Pr", colnames(lisa), value = TRUE)
if (length(p_col) == 0) stop("Could not find p-value column in localmoran output.")

sf_dat$lisa_I <- lisa[, "Ii"]
sf_dat$lisa_p <- lisa[, p_col[1]]

x_c <- x - mean(x, na.rm = TRUE)
lag_c <- x_lag - mean(x_lag, na.rm = TRUE)

sf_dat$hotspot <- !is.na(sf_dat$lisa_p) & (sf_dat$lisa_p < 0.05) & (sf_dat$lisa_I > 0)

sf_dat$lisa_class <- "Not significant"
sig <- !is.na(sf_dat$lisa_p) & (sf_dat$lisa_p < 0.05)
sf_dat$lisa_class[sig & x_c > 0 & lag_c > 0] <- "High-High"
sf_dat$lisa_class[sig & x_c < 0 & lag_c < 0] <- "Low-Low"
sf_dat$lisa_class[sig & x_c > 0 & lag_c < 0] <- "High-Low"
sf_dat$lisa_class[sig & x_c < 0 & lag_c > 0] <- "Low-High"

sf_dat$lisa_class <- factor(sf_dat$lisa_class,
  levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
)

#
lisa_map <- ggplot(sf_dat) +
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
    title = "Local Moran’s I clusters of unprotected flood exposure",
    subtitle = "Spatial clusters of Zone 1 exposure (p < 0.05)",
    fill = "Cluster",
    caption = "Source: BfG flood hazard data. Spatial analysis in R. Own processing."
  ) +
  map_theme()

save_plot(lisa_map, "map_lisa_clusters.png", w = 9, h = 7, subdir = "maps")

save_table(sf::st_drop_geometry(sf_dat) %>% select(lisa_I, lisa_p, hotspot, lisa_class),
           "lisa_clusters.csv")

saveRDS(FULL, file.path(out_dir, "FULL_with_PCs_index.rds"))
message("DONE. Outputs saved to: ", normalizePath(out_dir))