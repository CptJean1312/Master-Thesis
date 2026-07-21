#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(ggspatial)
})

options(scipen = 999)

# ============================================================
# Final thesis figure builder
# ------------------------------------------------------------
# Rebuilds every map and chart referenced in the active thesis draft
# with one design system and exports both PNG and vector PDF files.
# ============================================================

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"
external_root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS"

paths <- list(
  analysis_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/outputs/gpkg/corridor_wide_pca_rp100_analysis.gpkg"),
  curve_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/EXPOSURE_CURVES/outputs/gpkg/corridor_exposure_curve_metrics.gpkg"),
  landuse_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/gpkg/corridor_landuse_protection_analysis_rp100.gpkg"),
  landuse_summary_csv = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/tables/corridor_landuse_exposure_summary_by_group.csv"),
  protection_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/gpkg/corridor_protection_coverage.gpkg"),
  stream_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/gpkg/corridor_stream_context.gpkg"),
  rivers_gpkg = file.path(external_root, "CLIPPED + PROJECTED/am_riverwaterbody_de-basin.gpkg"),
  basin_gpkg = file.path(external_root, "PHYSISCH.nosync/RIGHT PROJECTION/ELBE BASIN PRO.gpkg"),
  germany_gpkg = file.path(external_root, "ANALYSIS.nosync/vg250_gemeinden_landonly.gpkg"),
  elbe_gpkg = file.path(external_root, "PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"),
  state_dir = file.path(external_root, "PHYSISCH.nosync/RIGHT PROJECTION"),
  city_dir = file.path(external_root, "MAP FEATURES/CITIES"),
  png_dir = file.path(root, "THESIS_DRAFTS/FIGURES/FINAL/png"),
  pdf_dir = file.path(root, "THESIS_DRAFTS/FIGURES/FINAL/pdf")
)

dir.create(paths$png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$pdf_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- unlist(paths[c(
  "analysis_gpkg", "curve_gpkg", "landuse_gpkg", "landuse_summary_csv",
  "protection_gpkg", "stream_gpkg", "rivers_gpkg", "basin_gpkg",
  "germany_gpkg", "elbe_gpkg"
)])

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required inputs:\n", paste(missing_files, collapse = "\n"))
}

# ---------------------------
# 1) Design system
# ---------------------------

font_family <- "sans"

colors <- list(
  ink = "#20252B",
  muted = "#626A73",
  grid = "#D9DEE3",
  state = "#8E979F",
  river = "#236B8E",
  teal = "#1B6A73",
  teal_light = "#B8D9D7",
  rust = "#B55432",
  no_event = "#D8DCE0",
  context = "#EEF1F3"
)

theme_thesis <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size, base_family = font_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = colors$grid, linewidth = 0.25),
      axis.title = element_text(color = colors$ink, size = base_size),
      axis.text = element_text(color = colors$muted, size = base_size - 1),
      legend.position = "bottom",
      legend.justification = "left",
      legend.title = element_text(face = "bold", color = colors$ink, size = base_size - 0.5),
      legend.text = element_text(color = colors$ink, size = base_size - 1.2),
      plot.title = element_text(face = "bold", color = colors$ink, size = base_size + 3, margin = margin(b = 4)),
      plot.subtitle = element_text(color = colors$muted, size = base_size, margin = margin(b = 8)),
      plot.caption = element_text(color = colors$muted, size = base_size - 2.4, hjust = 0, margin = margin(t = 7)),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = margin(10, 13, 10, 13)
    )
}

theme_thesis_map <- function(base_size = 10.5) {
  theme_thesis(base_size) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_line(color = "#E5E9EC", linewidth = 0.2)
    )
}

map_annotations <- function() {
  list(
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = unit(0.25, "cm"),
      pad_y = unit(0.25, "cm"),
      height = unit(0.8, "cm"),
      width = unit(0.8, "cm"),
      style = north_arrow_orienteering(text_family = font_family)
    ),
    annotation_scale(
      location = "bl",
      width_hint = 0.19,
      pad_x = unit(0.25, "cm"),
      pad_y = unit(0.25, "cm"),
      line_width = 0.45,
      text_cex = 0.62,
      text_family = font_family
    )
  )
}

save_figure <- function(plot, filename, width = 8.2, height = 5.8) {
  ggsave(
    file.path(paths$png_dir, paste0(filename, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 360,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(paths$pdf_dir, paste0(filename, ".pdf")),
    plot,
    width = width,
    height = height,
    device = "pdf",
    bg = "white",
    limitsize = FALSE
  )
}

# ---------------------------
# 2) Read analytical data
# ---------------------------

analysis <- st_read(paths$analysis_gpkg, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

curves <- st_read(paths$curve_gpkg, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

landuse <- st_read(paths$landuse_gpkg, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

protection <- st_read(paths$protection_gpkg, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

stream_context <- st_read(paths$stream_gpkg, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

landuse_summary <- read_csv(paths$landuse_summary_csv, show_col_types = FALSE)

target_crs <- st_crs(analysis)

to_target_crs <- function(x) {
  if (st_crs(x) != target_crs) st_transform(x, target_crs) else x
}

curves <- to_target_crs(curves)
landuse <- to_target_crs(landuse)
protection <- to_target_crs(protection)
stream_context <- to_target_crs(stream_context)

elbe <- st_read(paths$elbe_gpkg, quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_make_valid() %>%
  to_target_crs()

basin <- st_read(paths$basin_gpkg, quiet = TRUE) %>%
  st_make_valid() %>%
  to_target_crs()

germany_outline <- st_read(paths$germany_gpkg, quiet = TRUE) %>%
  summarise() %>%
  st_make_valid() %>%
  to_target_crs()

state_names <- c(
  "Niedersachsen", "Brandenburg", "Sachsen-Anhalt", "Sachsen",
  "Schleswig-Holstein", "Thüringen", "Hamburg", "Berlin",
  "Mecklenburg-Vorpommern", "Bremen"
)

states <- bind_rows(lapply(state_names, function(state_name) {
  state_path <- file.path(paths$state_dir, paste0(state_name, ".gpkg"))
  x <- st_read(state_path, quiet = TRUE) %>% st_make_valid() %>% to_target_crs()
  st_sf(state_name = state_name, geometry = st_geometry(x))
}))

city_files <- c(
  Hamburg = "Hamburg.gpkg",
  Magdeburg = "magdeburg.gpkg",
  Leipzig = "Leipzig.gpkg",
  Dresden = "Dresden.gpkg",
  Berlin = "Berlin.gpkg"
)

cities <- bind_rows(lapply(names(city_files), function(city_name) {
  x <- st_read(file.path(paths$city_dir, city_files[[city_name]]), quiet = TRUE) %>%
    st_make_valid() %>%
    to_target_crs() %>%
    summarise() %>%
    st_point_on_surface()
  st_sf(city = city_name, geometry = st_geometry(x))
}))

city_xy <- cbind(st_drop_geometry(cities), st_coordinates(cities)) %>%
  as_tibble() %>%
  mutate(
    label_x = X,
    label_y = Y + case_when(
      city %in% c("Hamburg", "Magdeburg", "Berlin") ~ 13000,
      city %in% c("Leipzig", "Dresden") ~ -15000,
      TRUE ~ 0
    )
  )

corridor_bbox <- st_bbox(st_union(st_geometry(analysis), st_geometry(cities)))
corridor_x <- c(corridor_bbox[["xmin"]] - 15000, corridor_bbox[["xmax"]] + 15000)
corridor_y <- c(corridor_bbox[["ymin"]] - 12000, corridor_bbox[["ymax"]] + 12000)

context_layers <- function(show_cities = TRUE) {
  layers <- list(
    geom_sf(data = states, fill = NA, color = colors$state, linewidth = 0.28, inherit.aes = FALSE),
    geom_sf(data = elbe, color = colors$river, linewidth = 0.55, inherit.aes = FALSE)
  )
  if (show_cities) {
    layers <- c(
      layers,
      list(
        geom_sf(data = cities, color = colors$ink, fill = "white", shape = 21, size = 1.8, stroke = 0.45, inherit.aes = FALSE),
        geom_label(
          data = city_xy,
          aes(label_x, label_y, label = city),
          inherit.aes = FALSE,
          family = font_family,
          size = 2.55,
          label.padding = unit(0.08, "lines"),
          label.r = unit(0.08, "lines"),
          linewidth = 0.15,
          color = colors$ink,
          fill = alpha("white", 0.88)
        )
      )
    )
  }
  layers
}

coord_corridor <- function() {
  coord_sf(xlim = corridor_x, ylim = corridor_y, expand = FALSE, datum = target_crs)
}

# ---------------------------
# 3) Study area and exposure maps
# ---------------------------

study_area_map <- ggplot() +
  geom_sf(data = germany_outline, fill = "#F7F8F9", color = "#A7AFB6", linewidth = 0.3) +
  geom_sf(data = basin, aes(fill = "Elbe basin"), color = "#8DA5AF", linewidth = 0.35, alpha = 0.72) +
  geom_sf(data = analysis, aes(fill = "RP500 corridor"), color = NA, alpha = 0.88) +
  context_layers() +
  coord_sf(xlim = c(390000, 940000), ylim = c(5390000, 6055000), expand = FALSE, datum = target_crs) +
  map_annotations() +
  scale_fill_manual(
    values = c("Elbe basin" = "#DDE8EC", "RP500 corridor" = colors$teal),
    name = "Study-area layer"
  ) +
  labs(
    title = "Study area and RP500 flood corridor",
    subtitle = "835 municipalities within the wider Elbe basin context",
    caption = str_wrap("Blue line: Elbe. Data: BKG VG250, BfG, and EFAS/JRC. Own processing.", 110)
  ) +
  theme_thesis_map()

save_figure(study_area_map, "map_study_area_corridor", height = 6.5)

map_rp100 <- ggplot(analysis) +
  geom_sf(aes(fill = flood_share_rp100), color = NA) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    name = "Flooded share"
  ) +
  labs(
    title = "Municipal RP100 flood exposure",
    subtitle = "Flooded share of total municipal area",
    caption = "Data: EFAS/JRC flood hazard rasters and BKG VG250 municipalities. Own processing."
  ) +
  theme_thesis_map()

save_figure(map_rp100, "map_rp100_exposure_corridor", height = 6.35)

vuln_limit <- max(abs(analysis$vuln_index_main_z), na.rm = TRUE)
map_vulnerability <- ggplot(analysis) +
  geom_sf(aes(fill = vuln_index_main_z), color = NA) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "#F7F7F7",
    high = "#D95F0E",
    midpoint = 0,
    limits = c(-vuln_limit, vuln_limit),
    name = "Vulnerability (z)"
  ) +
  labs(
    title = "Socio-economic vulnerability",
    subtitle = "Direction-coded, dimension-balanced index",
    caption = "Data: INKAR 07/2025 and the RP500 corridor definition. Own processing. One municipality has no INKAR value."
  ) +
  theme_thesis_map()

save_figure(map_vulnerability, "map_vulnerability_index_corridor", height = 6.35)

# ---------------------------
# 4) Vulnerability dimensions and RP100 comparisons
# ---------------------------

dimension_long <- analysis %>%
  select(
    AGS,
    vuln_dim_demographic_household_z,
    vuln_dim_deprivation_labour_z,
    vuln_dim_access_adaptive_capacity_z
  ) %>%
  pivot_longer(
    cols = starts_with("vuln_dim_"),
    names_to = "dimension",
    values_to = "score_z"
  ) %>%
  mutate(
    dimension = recode(
      dimension,
      vuln_dim_demographic_household_z = "Demographic and\nhousehold structure",
      vuln_dim_deprivation_labour_z = "Deprivation and\nlabour market",
      vuln_dim_access_adaptive_capacity_z = "Access and adaptive\ncapacity"
    )
  )

dimension_limit <- max(abs(dimension_long$score_z), na.rm = TRUE)

map_dimensions <- ggplot(dimension_long) +
  geom_sf(aes(fill = score_z), color = NA) +
  geom_sf(data = states, fill = NA, color = colors$state, linewidth = 0.22, inherit.aes = FALSE) +
  geom_sf(data = elbe, color = colors$river, linewidth = 0.45, inherit.aes = FALSE) +
  geom_sf(data = cities, color = colors$ink, fill = "white", shape = 21, size = 1.2, stroke = 0.35, inherit.aes = FALSE) +
  geom_text(
    data = city_xy,
    aes(label_x, label_y, label = city),
    inherit.aes = FALSE,
    family = font_family,
    size = 2.0,
    color = colors$ink,
    check_overlap = TRUE
  ) +
  facet_wrap(~dimension, ncol = 1) +
  coord_corridor() +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "#F7F7F7",
    high = "#D95F0E",
    midpoint = 0,
    limits = c(-dimension_limit, dimension_limit),
    name = "Dimension score (z)"
  ) +
  labs(
    title = "Dimensions of socio-economic vulnerability",
    subtitle = "Higher values indicate greater vulnerability or lower adaptive capacity",
    caption = str_wrap("Blue line: Elbe. Data: INKAR 07/2025 and the RP500 corridor definition. Own processing.", 105)
  ) +
  theme_thesis_map() +
  theme(strip.text = element_text(face = "bold", color = colors$ink, size = 10.5))

save_figure(map_dimensions, "map_vulnerability_dimensions", height = 10.2)

analysis_tbl <- analysis %>%
  st_drop_geometry() %>%
  mutate(
    vuln_quintile_num = if_else(is.na(vuln_index_main_z), NA_integer_, ntile(vuln_index_main_z, 5)),
    vuln_quintile = factor(
      vuln_quintile_num,
      levels = 1:5,
      labels = c("Q1 lowest", "Q2", "Q3", "Q4", "Q5 highest")
    ),
    exposure_tercile = factor(
      ntile(flood_share_rp100, 3),
      levels = 1:3,
      labels = c("Low exposure", "Medium exposure", "High exposure")
    ),
    vulnerability_tercile = factor(
      ntile(vuln_index_main_z, 3),
      levels = 1:3,
      labels = c("Low vulnerability", "Medium vulnerability", "High vulnerability")
    )
  )

pearson_rp100 <- cor(analysis_tbl$vuln_index_main_z, analysis_tbl$flood_share_rp100, use = "complete.obs")
spearman_rp100 <- cor(
  analysis_tbl$vuln_index_main_z,
  analysis_tbl$flood_share_rp100,
  use = "complete.obs",
  method = "spearman"
)

scatter_rp100 <- ggplot(analysis_tbl, aes(vuln_index_main_z, flood_share_rp100)) +
  geom_point(alpha = 0.42, size = 1.5, color = colors$teal) +
  geom_smooth(method = "lm", se = TRUE, color = colors$rust, fill = alpha(colors$rust, 0.16), linewidth = 0.75) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "RP100 exposure and socio-economic vulnerability",
    subtitle = sprintf("Pearson r = %.3f; Spearman rho = %.3f", pearson_rp100, spearman_rp100),
    x = "Socio-economic vulnerability index (z)",
    y = "Flooded share of municipal area",
    caption = "Each point represents one municipality. Berga-Wünschendorf is omitted because no INKAR value is available."
  ) +
  theme_thesis()

save_figure(scatter_rp100, "plot_vulnerability_vs_rp100_exposure")

box_rp100 <- ggplot(
  analysis_tbl %>% filter(!is.na(vuln_quintile)),
  aes(vuln_quintile, flood_share_rp100)
) +
  geom_boxplot(fill = colors$teal_light, color = colors$ink, linewidth = 0.45, outlier.alpha = 0.26, outlier.size = 1.1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = colors$rust, color = "white", stroke = 0.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "RP100 exposure by vulnerability quintile",
    subtitle = "Rust circles show group means; boxes show municipal distributions",
    x = "Vulnerability quintile",
    y = "Flooded share of municipal area",
    caption = "Data: EFAS/JRC flood hazard rasters and final INKAR vulnerability index. Own processing."
  ) +
  theme_thesis()

save_figure(box_rp100, "plot_rp100_exposure_by_vulnerability_quintile")

bivar_levels <- c(
  "Low exposure / Low vulnerability",
  "Medium exposure / Low vulnerability",
  "High exposure / Low vulnerability",
  "Low exposure / Medium vulnerability",
  "Medium exposure / Medium vulnerability",
  "High exposure / Medium vulnerability",
  "Low exposure / High vulnerability",
  "Medium exposure / High vulnerability",
  "High exposure / High vulnerability",
  "No vulnerability data"
)

bivar_colors <- c(
  "Low exposure / Low vulnerability" = "#E8E8E8",
  "Medium exposure / Low vulnerability" = "#ACE4E4",
  "High exposure / Low vulnerability" = "#5AC8C8",
  "Low exposure / Medium vulnerability" = "#DFB0D6",
  "Medium exposure / Medium vulnerability" = "#A5ADD3",
  "High exposure / Medium vulnerability" = "#5698B9",
  "Low exposure / High vulnerability" = "#BE64AC",
  "Medium exposure / High vulnerability" = "#8C62AA",
  "High exposure / High vulnerability" = "#3B4994",
  "No vulnerability data" = "#777777"
)

bivar_labels <- c(
  "Low exposure / Low vulnerability" = "Low exp. / Low vuln.",
  "Medium exposure / Low vulnerability" = "Medium exp. / Low vuln.",
  "High exposure / Low vulnerability" = "High exp. / Low vuln.",
  "Low exposure / Medium vulnerability" = "Low exp. / Medium vuln.",
  "Medium exposure / Medium vulnerability" = "Medium exp. / Medium vuln.",
  "High exposure / Medium vulnerability" = "High exp. / Medium vuln.",
  "Low exposure / High vulnerability" = "Low exp. / High vuln.",
  "Medium exposure / High vulnerability" = "Medium exp. / High vuln.",
  "High exposure / High vulnerability" = "High exp. / High vuln.",
  "No vulnerability data" = "No vulnerability data"
)

bivar_map_data <- analysis %>%
  left_join(
    analysis_tbl %>%
      select(AGS, exposure_tercile, vulnerability_tercile),
    by = "AGS"
  ) %>%
  mutate(
    bivariate_class = if_else(
      is.na(vulnerability_tercile),
      "No vulnerability data",
      paste(exposure_tercile, vulnerability_tercile, sep = " / ")
    ),
    bivariate_class = factor(bivariate_class, levels = bivar_levels)
  )

bivar_map <- ggplot(bivar_map_data) +
  geom_sf(aes(fill = bivariate_class), color = NA) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_manual(
    values = bivar_colors,
    labels = bivar_labels,
    drop = FALSE,
    name = "RP100 exposure / vulnerability"
  ) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE)) +
  labs(
    title = "Combined RP100 exposure and vulnerability",
    subtitle = "Terciles within the RP500 corridor; dark purple marks high-high municipalities",
    caption = "Data: EFAS/JRC flood hazard rasters and final INKAR vulnerability index. Own processing."
  ) +
  theme_thesis_map() +
  theme(
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 8.4),
    legend.key.width = unit(0.42, "cm"),
    legend.spacing.x = unit(0.12, "cm")
  )

save_figure(bivar_map, "map_exposure_vulnerability_bivariate", height = 6.65)

# ---------------------------
# 5) Land-cover refinement
# ---------------------------

landuse_summary_plot <- landuse_summary %>%
  mutate(
    thesis_group = factor(
      thesis_group,
      levels = c("water", "perennial_vegetation", "open_or_seasonal", "artificial"),
      labels = c("Water", "Perennial vegetation", "Open or seasonal", "Artificial land")
    ),
    flooded_area_km2 = flooded_area_m2 / 1e6,
    label = paste0(
      number(flooded_area_km2, accuracy = 0.1, big.mark = ","),
      " km²  (", percent(share_of_all_flooded_landcover_area, accuracy = 0.1), ")"
    )
  )

composition_plot <- ggplot(landuse_summary_plot, aes(thesis_group, flooded_area_km2, fill = thesis_group)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = label), hjust = -0.05, family = font_family, size = 3.15, color = colors$ink) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Artificial land" = colors$rust,
      "Open or seasonal" = "#D6A84B",
      "Perennial vegetation" = "#5A8F62",
      "Water" = colors$river
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 1, big.mark = ","),
    expand = expansion(mult = c(0, 0.30))
  ) +
  labs(
    title = "Land-cover composition of modeled RP100 flooding",
    subtitle = "Artificial land accounts for 6.0% of flooded land-cover area",
    x = NULL,
    y = "Flooded area (km²)",
    caption = "Data: DLR Land Cover DE 2015 and EFAS/JRC RP100 flood hazard raster. Own processing."
  ) +
  theme_thesis()

save_figure(composition_plot, "plot_rp100_flooded_landcover_composition", height = 5.3)

artificial_map <- ggplot(landuse) +
  geom_sf(aes(fill = rp100_artificial_flood_share_of_group), color = NA) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    na.value = "#C8CDD1",
    name = "Flooded artificial land"
  ) +
  labs(
    title = "RP100 artificial-land exposure",
    subtitle = "Share of DLR Artificial Land flooded in each municipality",
    caption = "Data: DLR Land Cover DE 2015 and EFAS/JRC RP100 flood hazard raster. Own processing. Grey: no Artificial Land."
  ) +
  theme_thesis_map()

save_figure(artificial_map, "map_rp100_artificial_land_flood_share", height = 6.35)

landuse_tbl <- landuse %>%
  st_drop_geometry() %>%
  mutate(
    vuln_quintile = factor(
      vuln_quintile_num,
      levels = 1:5,
      labels = c("Q1 lowest", "Q2", "Q3", "Q4", "Q5 highest")
    )
  )

box_artificial <- ggplot(
  landuse_tbl %>% filter(!is.na(vuln_quintile), !is.na(rp100_artificial_flood_share_of_group)),
  aes(vuln_quintile, rp100_artificial_flood_share_of_group)
) +
  geom_boxplot(fill = colors$teal_light, color = colors$ink, linewidth = 0.45, outlier.alpha = 0.26, outlier.size = 1.1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = colors$rust, color = "white", stroke = 0.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Artificial-land exposure by vulnerability quintile",
    subtitle = "Rust circles show group means; boxes show municipal distributions",
    x = "Vulnerability quintile",
    y = "Flooded share of artificial land",
    caption = "Data: DLR Land Cover DE 2015, EFAS/JRC RP100 flood hazard raster, and final INKAR vulnerability index. Own processing."
  ) +
  theme_thesis()

save_figure(box_artificial, "plot_rp100_artificial_land_exposure_by_vulnerability_quintile")

# ---------------------------
# 6) Exposure curves
# ---------------------------

curve_levels <- c("early exposure", "gradual increase", "delayed jump")
curve_labels <- c(
  "early exposure" = "Early exposure",
  "gradual increase" = "Gradual increase",
  "delayed jump" = "Delayed jump"
)
curve_colors <- c(
  "early exposure" = colors$teal,
  "gradual increase" = "#D6A84B",
  "delayed jump" = colors$rust
)

curve_map_data <- curves %>%
  mutate(exposure_curve_type = factor(exposure_curve_type, levels = curve_levels))

curve_type_map <- ggplot(curve_map_data) +
  geom_sf(aes(fill = exposure_curve_type), color = NA) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_manual(values = curve_colors, labels = curve_labels, drop = FALSE, name = "Exposure-curve type") +
  labs(
    title = "Exposure dynamics across return periods",
    subtitle = "Classification from RP10, RP20, RP50, RP100, RP200, and RP500 flood shares",
    caption = "Data: EFAS/JRC multi-return-period flood hazard rasters. Own processing."
  ) +
  theme_thesis_map()

save_figure(curve_type_map, "map_exposure_curve_types", height = 6.35)

curve_long <- analysis_tbl %>%
  select(
    AGS,
    vuln_quintile,
    flood_share_rp10,
    flood_share_rp20,
    flood_share_rp50,
    flood_share_rp100,
    flood_share_rp200,
    flood_share_rp500
  ) %>%
  filter(!is.na(vuln_quintile)) %>%
  pivot_longer(
    starts_with("flood_share_rp"),
    names_to = "return_period",
    values_to = "flood_share"
  ) %>%
  mutate(return_period = as.numeric(str_remove(return_period, "flood_share_rp"))) %>%
  group_by(vuln_quintile, return_period) %>%
  summarise(mean_flood_share = mean(flood_share, na.rm = TRUE), .groups = "drop")

quintile_colors <- c("#2C7BB6", "#6FB1B3", "#B8C76B", "#E4A44A", "#C34F34")

curve_plot <- ggplot(
  curve_long,
  aes(return_period, mean_flood_share, color = vuln_quintile, group = vuln_quintile)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  scale_x_continuous(breaks = c(10, 20, 50, 100, 200, 500)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.06))) +
  scale_color_manual(values = quintile_colors, name = "Vulnerability quintile") +
  labs(
    title = "Mean flood-exposure curves by vulnerability quintile",
    subtitle = "Average municipal flood share across six return periods",
    x = "Return period (years)",
    y = "Mean flooded share of municipal area",
    caption = "Data: EFAS/JRC multi-return-period flood hazard rasters and final INKAR vulnerability index. Own processing."
  ) +
  theme_thesis()

save_figure(curve_plot, "plot_exposure_curves_by_vulnerability_quintile")

# ---------------------------
# 7) Protection and modeled loss
# ---------------------------

protection <- protection %>%
  mutate(
    protection_status_display = factor(
      if_else(protection_available, "Finite protection RP", "No simulated loss event"),
      levels = c("No simulated loss event", "Finite protection RP")
    )
  )

coverage_map <- ggplot(protection) +
  geom_sf(aes(fill = protection_status_display), color = "white", linewidth = 0.04) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_manual(
    values = c("No simulated loss event" = colors$no_event, "Finite protection RP" = colors$teal),
    drop = FALSE,
    name = "Portfolio status"
  ) +
  labs(
    title = "Modeled loss/protection portfolio coverage",
    subtitle = "280 municipalities have positive modeled losses; 555 have no simulated loss event",
    caption = "Data: provider loss/protection portfolio and RP500 corridor municipalities. Own processing."
  ) +
  theme_thesis_map()

save_figure(coverage_map, "map_corridor_protection_coverage", height = 6.35)

rp_levels <- c("<10", "10-20", "20-50", "50-100", "100-200", "200-500", "500-1000", ">=1000")
protection_rp <- protection %>%
  mutate(
    protection_class = cut(
      protection_return_period,
      breaks = c(-Inf, 10, 20, 50, 100, 200, 500, 1000, Inf),
      labels = rp_levels,
      right = FALSE
    )
  )

rp_map <- ggplot() +
  geom_sf(data = protection, fill = colors$no_event, color = "white", linewidth = 0.035) +
  geom_sf(data = protection_rp %>% filter(protection_available), aes(fill = protection_class), color = "white", linewidth = 0.04) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_manual(
    values = c(
      "<10" = "#9E3D22", "10-20" = "#C75C2D", "20-50" = "#E18A42",
      "50-100" = "#E8C267", "100-200" = "#B8D59A", "200-500" = "#75B8A5",
      "500-1000" = "#3E8FA3", ">=1000" = "#285A8C"
    ),
    drop = FALSE,
    name = "Finite protection RP (years)"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "Finite protection return periods",
    subtitle = "Positive-loss municipalities only; grey indicates no simulated loss event",
    caption = "Higher return periods imply less frequent modeled loss occurrence. Data: provider loss/protection portfolio. Own processing."
  ) +
  theme_thesis_map()

save_figure(rp_map, "map_corridor_protection_return_period", height = 6.45)

analysis_protection <- analysis %>%
  left_join(
    protection %>%
      st_drop_geometry() %>%
      select(AGS, protection_available, annual_loss_probability_model, protection_return_period),
    by = "AGS"
  ) %>%
  left_join(
    analysis_tbl %>% select(AGS, vuln_quintile),
    by = "AGS"
  ) %>%
  mutate(
    annual_loss_class = factor(
      case_when(
        annual_loss_probability_model == 0 ~ "0: no simulated loss",
        annual_loss_probability_model <= 0.01 ~ ">0-1%",
        annual_loss_probability_model <= 0.02 ~ "1-2%",
        annual_loss_probability_model <= 0.05 ~ "2-5%",
        annual_loss_probability_model <= 0.10 ~ "5-10%",
        annual_loss_probability_model > 0.10 ~ ">10%",
        TRUE ~ NA_character_
      ),
      levels = c("0: no simulated loss", ">0-1%", "1-2%", "2-5%", "5-10%", ">10%")
    )
  )

annual_loss_map <- ggplot(analysis_protection) +
  geom_sf(aes(fill = annual_loss_class), color = NA) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_manual(
    values = c(
      "0: no simulated loss" = colors$no_event,
      ">0-1%" = "#FFF2B2",
      "1-2%" = "#F9D56E",
      "2-5%" = "#EE994A",
      "5-10%" = "#D65A35",
      ">10%" = "#9E2A2B"
    ),
    drop = FALSE,
    name = "Annual loss probability"
  ) +
  labs(
    title = "Modeled annual loss probability",
    subtitle = "No-event municipalities are retained as zero within the provider portfolio",
    caption = "Data: provider loss/protection portfolio and RP500 corridor municipalities. Own processing."
  ) +
  theme_thesis_map()

save_figure(annual_loss_map, "map_modeled_annual_loss_probability", height = 6.35)

portfolio_summary <- analysis_protection %>%
  st_drop_geometry() %>%
  filter(!is.na(vuln_quintile), !is.na(protection_available)) %>%
  count(vuln_quintile, protection_available) %>%
  group_by(vuln_quintile) %>%
  mutate(share = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    status = factor(
      if_else(protection_available, "Positive modeled loss", "No simulated loss event"),
      levels = c("No simulated loss event", "Positive modeled loss")
    )
  )

portfolio_plot <- ggplot(portfolio_summary, aes(vuln_quintile, share, fill = status)) +
  geom_col(color = "white", linewidth = 0.25, width = 0.76) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(
    values = c("No simulated loss event" = colors$no_event, "Positive modeled loss" = colors$teal),
    name = "Portfolio status"
  ) +
  labs(
    title = "Modeled loss occurrence by vulnerability quintile",
    subtitle = "Positive-loss status is more common in the upper vulnerability groups",
    x = "Vulnerability quintile",
    y = "Share of municipalities",
    caption = "Data: provider loss/protection portfolio and final INKAR vulnerability index. Own processing."
  ) +
  theme_thesis()

save_figure(portfolio_plot, "plot_protection_status_by_vulnerability_quintile")

annual_loss_plot <- ggplot(
  analysis_protection %>% st_drop_geometry() %>% filter(!is.na(vuln_quintile)),
  aes(vuln_quintile, annual_loss_probability_model)
) +
  geom_boxplot(fill = "#F2D5C8", color = colors$ink, linewidth = 0.45, outlier.alpha = 0.25, outlier.size = 1.1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = colors$rust, color = "white", stroke = 0.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "Annual loss probability by vulnerability quintile",
    subtitle = "Rust circles show means; no-event municipalities are retained as zero",
    x = "Vulnerability quintile",
    y = "Modeled annual loss probability",
    caption = "Data: provider loss/protection portfolio and final INKAR vulnerability index. Own processing."
  ) +
  theme_thesis()

save_figure(annual_loss_plot, "plot_annual_loss_probability_by_vulnerability_quintile")

rp_quintile_plot <- ggplot(
  analysis_protection %>%
    st_drop_geometry() %>%
    filter(!is.na(vuln_quintile), !is.na(protection_return_period)),
  aes(vuln_quintile, protection_return_period)
) +
  geom_boxplot(fill = "#C9DEE8", color = colors$ink, linewidth = 0.45, outlier.alpha = 0.28, outlier.size = 1.1) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.7, fill = colors$rust, color = "white", stroke = 0.5) +
  scale_y_log10(labels = label_number(big.mark = ",")) +
  labs(
    title = "Finite protection return period by vulnerability quintile",
    subtitle = "Positive-loss municipalities only; rust diamonds show medians",
    x = "Vulnerability quintile",
    y = "Finite protection return period (years, log scale)",
    caption = "No-event municipalities are not assigned an artificial finite return period. Data: provider portfolio and INKAR."
  ) +
  theme_thesis()

save_figure(rp_quintile_plot, "plot_protection_return_period_by_vulnerability_quintile")

# ---------------------------
# 8) Supporting stream-context diagnostic
# ---------------------------

rivers <- st_read(paths$rivers_gpkg, quiet = TRUE) %>% to_target_crs()
stream_display_levels <- c(
  "directly on/within 1 km of Elbe main stem",
  "intersects major WFD river/tributary",
  "only minor/other WFD water bodies",
  "no WFD river body intersects municipality"
)

stream_context <- stream_context %>%
  mutate(river_context_class = factor(river_context_class, levels = stream_display_levels))

stream_map <- ggplot() +
  geom_sf(data = stream_context, fill = colors$context, color = "white", linewidth = 0.035) +
  geom_sf(data = rivers, color = "#9CC7DB", linewidth = 0.16, alpha = 0.45) +
  geom_sf(
    data = stream_context %>% filter(no_simulated_loss_event),
    aes(fill = river_context_class),
    color = "white",
    linewidth = 0.035
  ) +
  context_layers() +
  coord_corridor() +
  map_annotations() +
  scale_fill_manual(
    values = c(
      "directly on/within 1 km of Elbe main stem" = "#2C7BB6",
      "intersects major WFD river/tributary" = "#72B7B2",
      "only minor/other WFD water bodies" = "#E3B34B",
      "no WFD river body intersects municipality" = "#B55432"
    ),
    labels = c(
      "directly on/within 1 km of Elbe main stem" = "Elbe main stem (<=1 km)",
      "intersects major WFD river/tributary" = "Major WFD river or tributary",
      "only minor/other WFD water bodies" = "Only minor/other WFD water bodies",
      "no WFD river body intersects municipality" = "No WFD river-body intersection"
    ),
    drop = FALSE,
    name = "No-event municipality context"
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(
    title = "No-event municipalities by river context",
    subtitle = "A GIS diagnostic of main-stem, tributary, and smaller-waterbody settings",
    caption = "Pale blue: WFD river water bodies. Dark blue: Elbe. Smaller streams are not fully represented in the provider simulation."
  ) +
  theme_thesis_map() +
  theme(legend.text = element_text(size = 7.8))

save_figure(stream_map, "map_no_loss_stream_context", height = 6.55)

# ---------------------------
# 9) Export manifest
# ---------------------------

manifest <- tribble(
  ~figure, ~filename, ~description,
  "2.1", "map_study_area_corridor", "Study area and RP500 corridor",
  "5.1", "map_rp100_exposure_corridor", "RP100 total-area exposure",
  "5.2", "map_vulnerability_index_corridor", "Final vulnerability index",
  "5.3", "map_vulnerability_dimensions", "Three vulnerability dimensions",
  "5.4", "plot_vulnerability_vs_rp100_exposure", "Exposure-vulnerability scatterplot",
  "5.5", "plot_rp100_exposure_by_vulnerability_quintile", "RP100 exposure by vulnerability quintile",
  "5.6", "map_exposure_vulnerability_bivariate", "Bivariate RP100 exposure-vulnerability map",
  "5.7", "plot_rp100_flooded_landcover_composition", "Flooded land-cover composition",
  "5.8", "map_rp100_artificial_land_flood_share", "Artificial-land exposure map",
  "5.9", "plot_rp100_artificial_land_exposure_by_vulnerability_quintile", "Artificial-land exposure by vulnerability quintile",
  "5.10", "map_exposure_curve_types", "Exposure-curve types",
  "5.11", "plot_exposure_curves_by_vulnerability_quintile", "Mean exposure curves by vulnerability quintile",
  "5.12", "map_corridor_protection_coverage", "Provider portfolio coverage",
  "5.13", "map_corridor_protection_return_period", "Finite protection return periods",
  "5.14", "map_modeled_annual_loss_probability", "Modeled annual loss probability",
  "5.15", "plot_protection_status_by_vulnerability_quintile", "Portfolio status by vulnerability quintile",
  "5.16", "plot_annual_loss_probability_by_vulnerability_quintile", "Annual loss probability by vulnerability quintile",
  "5.17", "plot_protection_return_period_by_vulnerability_quintile", "Finite protection return period by vulnerability quintile",
  "support", "map_no_loss_stream_context", "No-event stream-context diagnostic"
) %>%
  mutate(
    png_path = file.path(paths$png_dir, paste0(filename, ".png")),
    pdf_path = file.path(paths$pdf_dir, paste0(filename, ".pdf")),
    png_exists = file.exists(png_path),
    pdf_exists = file.exists(pdf_path)
  )

write_csv(manifest, file.path(root, "THESIS_DRAFTS/FIGURES/FINAL/FIGURE_MANIFEST.csv"))

message("Final figures written to: ", dirname(paths$png_dir))
message("Manifest: ", file.path(root, "THESIS_DRAFTS/FIGURES/FINAL/FIGURE_MANIFEST.csv"))
