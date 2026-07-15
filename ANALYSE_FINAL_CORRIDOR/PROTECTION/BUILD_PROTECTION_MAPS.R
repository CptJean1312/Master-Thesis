#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

options(scipen = 999)

# ============================================================
# Protection coverage maps for the EFAS RP500 Elbe corridor
# ------------------------------------------------------------
# Outputs:
# - coverage map: municipalities with/without protection data
# - protection return-period map for matched municipalities
# - CSV tables of matched and unmatched corridor municipalities
# - GeoPackage with joined protection attributes
# ============================================================

# ---------------------------
# 1) Input paths
# ---------------------------

corridor_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg"
corridor_layer <- "municipalities_corridor"

protection_csv <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PROTECTION/elbe_protection_level_mun.csv"

elbe_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"
elbe_layer <- "elbe"

# ---------------------------
# 2) Output paths
# ---------------------------

output_root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs"
map_dir <- file.path(output_root, "maps")
table_dir <- file.path(output_root, "tables")
gpkg_dir <- file.path(output_root, "gpkg")

for (dir_path in c(output_root, map_dir, table_dir, gpkg_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# ---------------------------
# 3) Load and join data
# ---------------------------

message("Loading RP500 corridor municipalities ...")
corridor <- st_read(corridor_path, layer = corridor_layer, quiet = TRUE) %>%
  mutate(
    AGS = coalesce(as.character(AGS), as.character(Gemeindeschlüssel_AGS)),
    mun_name = coalesce(mun_name, GeografischerName_GEN)
  )

if (st_crs(corridor)$epsg != 25832) {
  corridor <- st_transform(corridor, 25832)
}

message("Loading protection table ...")
protection <- read_csv(
  protection_csv,
  col_types = cols(
    ags = col_character(),
    municipality_name = col_character(),
    n_nonzero_years = col_double(),
    annual_loss_probability = col_double(),
    protection_return_period = col_double()
  ),
  show_col_types = FALSE
)

message("Loading Elbe geometry ...")
elbe <- st_read(elbe_path, layer = elbe_layer, quiet = TRUE)

if (st_crs(elbe)$epsg != 25832) {
  elbe <- st_transform(elbe, 25832)
}

map_data <- corridor %>%
  left_join(protection, by = c("AGS" = "ags")) %>%
  mutate(
    protection_available = !is.na(protection_return_period),
    protection_status = if_else(
      protection_available,
      "Finite protection RP available",
      "No simulated loss event"
    ),
    loss_portfolio_status = if_else(
      protection_available,
      "positive_modeled_loss",
      "zero_modeled_loss_probability"
    ),
    provider_interpretation = if_else(
      protection_available,
      "Municipality has positive modeled loss in the provider output and an estimated finite protection return period.",
      "Municipality has no simulated loss event in the provider output. Following provider clarification, modeled annual loss probability is treated as zero; small streams are not represented and should be noted as a limitation."
    ),
    municipality_name_final = coalesce(municipality_name, mun_name),
    pop_total_num = suppressWarnings(as.numeric(pop_total)),
    no_simulated_loss_event = !protection_available,
    annual_loss_probability_model = if_else(
      protection_available,
      annual_loss_probability,
      0
    )
  )

matched_n <- sum(map_data$protection_available, na.rm = TRUE)
total_n <- nrow(map_data)
missing_n <- total_n - matched_n
coverage_share <- matched_n / total_n

message("Matched corridor municipalities: ", matched_n, " of ", total_n)
message("Unmatched corridor municipalities: ", missing_n)

# ---------------------------
# 4) Export joined data tables
# ---------------------------

matched_table <- map_data %>%
  st_drop_geometry() %>%
  filter(protection_available) %>%
  transmute(
    AGS,
    municipality_name = municipality_name_final,
    protection_return_period,
    annual_loss_probability,
    n_nonzero_years
  ) %>%
  arrange(AGS)

all_corridor_table <- map_data %>%
  st_drop_geometry() %>%
  transmute(
    AGS,
    municipality_name = municipality_name_final,
    pop_total = pop_total_num,
    protection_available,
    protection_status,
    loss_portfolio_status,
    provider_interpretation,
    no_simulated_loss_event,
    annual_loss_probability_model,
    protection_return_period,
    annual_loss_probability,
    n_nonzero_years
  ) %>%
  arrange(desc(protection_available), AGS)

missing_table <- map_data %>%
  st_drop_geometry() %>%
  filter(!protection_available) %>%
  transmute(
    AGS,
    municipality_name = municipality_name_final
  ) %>%
  arrange(AGS)

write_csv(
  matched_table,
  file.path(table_dir, "corridor_municipalities_with_protection.csv")
)

write_csv(
  all_corridor_table,
  file.path(table_dir, "corridor_protection_level_all_corridor_municipalities.csv")
)

write_csv(
  missing_table,
  file.path(table_dir, "corridor_municipalities_without_protection.csv")
)

joined_gpkg <- file.path(gpkg_dir, "corridor_protection_coverage.gpkg")
if (file.exists(joined_gpkg)) {
  unlink(joined_gpkg)
}

st_write(
  map_data,
  dsn = joined_gpkg,
  layer = "corridor_protection_coverage",
  quiet = TRUE
)

# ---------------------------
# 5) Map helpers
# ---------------------------

map_bbox <- st_bbox(map_data)

base_map_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_line(color = "grey90", linewidth = 0.2),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      plot.title = element_text(face = "bold", size = base_size + 4),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 0)
    )
}

elbe_layer_geom <- geom_sf(
  data = elbe,
  color = "#08306b",
  linewidth = 0.55,
  alpha = 0.95,
  inherit.aes = FALSE
)

make_city_label <- function(x) {
  x %>%
    str_replace(", Freie und Hansestadt", "") %>%
    str_replace(", Landeshauptstadt", "") %>%
    str_replace(", Hansestadt", "") %>%
    str_replace(", Lutherstadt", "") %>%
    str_replace(", Hochschulstadt", "") %>%
    str_replace(", Stadt", "") %>%
    str_replace(" / .*", "") %>%
    str_trim()
}

elbe_union <- st_union(elbe)
distance_to_elbe_m <- as.numeric(st_distance(map_data, elbe_union))

city_labels <- map_data %>%
  mutate(
    distance_to_elbe_m = distance_to_elbe_m,
    city_label = make_city_label(municipality_name_final),
    city_portfolio_status = if_else(
      protection_available,
      "Large municipality with finite protection RP",
      "Large municipality without simulated loss event"
    )
  ) %>%
  filter(
    !is.na(pop_total_num),
    pop_total_num >= 50000 | (distance_to_elbe_m <= 5000 & pop_total_num >= 25000)
  ) %>%
  arrange(desc(pop_total_num)) %>%
  st_point_on_surface()

large_city_table <- city_labels %>%
  st_drop_geometry() %>%
  transmute(
    AGS,
    municipality_name = municipality_name_final,
    city_label,
    pop_total = pop_total_num,
    protection_available,
    loss_portfolio_status,
    protection_return_period,
    annual_loss_probability_model,
    annual_loss_probability,
    distance_to_elbe_m
  ) %>%
  arrange(desc(pop_total))

write_csv(
  large_city_table,
  file.path(table_dir, "large_municipalities_portfolio_status.csv")
)

city_point_layer <- geom_sf(
  data = city_labels,
  aes(color = city_portfolio_status),
  size = 2.0,
  alpha = 0.95,
  inherit.aes = FALSE
)

city_label_layer <- geom_sf_label(
  data = city_labels,
  aes(label = city_label, color = city_portfolio_status),
  fill = "white",
  linewidth = 0.12,
  label.padding = unit(0.08, "lines"),
  size = 2.45,
  fontface = "bold",
  alpha = 0.9,
  show.legend = FALSE,
  inherit.aes = FALSE
)

city_color_scale <- scale_color_manual(
  values = c(
    "Large municipality with finite protection RP" = "#08306b",
    "Large municipality without simulated loss event" = "#b2182b"
  ),
  labels = c(
    "Large municipality with finite protection RP" = "City with finite RP",
    "Large municipality without simulated loss event" = "City no event"
  ),
  name = "Large municipalities"
)

# ---------------------------
# 6) Coverage map
# ---------------------------

coverage_map <- ggplot(map_data) +
  geom_sf(aes(fill = protection_status), color = "white", linewidth = 0.06) +
  elbe_layer_geom +
  city_point_layer +
  city_label_layer +
  scale_fill_manual(
    values = c(
      "Finite protection RP available" = "#1967a3",
      "No simulated loss event" = "#d8d3c7"
    ),
    labels = c(
      "Finite protection RP available" = "Finite protection RP",
      "No simulated loss event" = "No simulated loss event"
    ),
    name = NULL
  ) +
  city_color_scale +
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE),
    color = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  coord_sf(
    xlim = c(map_bbox$xmin, map_bbox$xmax),
    ylim = c(map_bbox$ymin, map_bbox$ymax),
    expand = FALSE
  ) +
  labs(
    title = "Protection data coverage in the RP500 Elbe corridor",
    subtitle = paste0(
      matched_n, " of ", total_n,
      " corridor municipalities matched to elbe_protection_level_mun.csv (",
      round(100 * coverage_share, 1), "%)."
    ),
    caption = "Data: EFAS/JRC-derived RP500 corridor; elbe_protection_level_mun.csv. Dark blue line: Elbe. Red city labels: large municipalities with no simulated loss event."
  ) +
  base_map_theme()

ggsave(
  filename = file.path(map_dir, "map_corridor_protection_coverage.png"),
  plot = coverage_map,
  width = 11,
  height = 8.5,
  dpi = 320,
  bg = "white"
)

# ---------------------------
# 7) Protection return-period map
# ---------------------------

rp_map_data <- map_data %>%
  filter(protection_available) %>%
  mutate(
    protection_class = cut(
      protection_return_period,
      breaks = c(-Inf, 10, 20, 50, 100, 200, 500, 1000, Inf),
      labels = c("<10", "10-20", "20-50", "50-100", "100-200", "200-500", "500-1000", ">=1000"),
      right = FALSE
    )
  )

rp_map <- ggplot() +
  geom_sf(data = map_data, fill = "#eeeeea", color = "white", linewidth = 0.04) +
  geom_sf(data = rp_map_data, aes(fill = protection_class), color = "white", linewidth = 0.05) +
  elbe_layer_geom +
  city_point_layer +
  city_label_layer +
  scale_fill_manual(
    values = c(
      "<10" = "#8c2d04",
      "10-20" = "#cc4c02",
      "20-50" = "#ec7014",
      "50-100" = "#fec44f",
      "100-200" = "#c7e9b4",
      "200-500" = "#7fcdbb",
      "500-1000" = "#41b6c4",
      ">=1000" = "#225ea8"
    ),
    drop = FALSE,
    name = "Protection RP"
  ) +
  city_color_scale +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    color = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  coord_sf(
    xlim = c(map_bbox$xmin, map_bbox$xmax),
    ylim = c(map_bbox$ymin, map_bbox$ymax),
    expand = FALSE
  ) +
  labs(
    title = "Protection return period for positive-loss corridor municipalities",
    subtitle = paste0(
      matched_n,
      " municipalities with finite RP estimates; grey = no simulated loss event."
    ),
    caption = "Higher return periods indicate lower annual exceedance/loss probability among positive-loss municipalities. Dark blue line: Elbe. Red city labels: no simulated loss event."
  ) +
  base_map_theme()

ggsave(
  filename = file.path(map_dir, "map_corridor_protection_return_period.png"),
  plot = rp_map,
  width = 11,
  height = 8.5,
  dpi = 320,
  bg = "white"
)

message("Done.")
message("Coverage map: ", file.path(map_dir, "map_corridor_protection_coverage.png"))
message("Return-period map: ", file.path(map_dir, "map_corridor_protection_return_period.png"))
message("Matched table: ", file.path(table_dir, "corridor_municipalities_with_protection.csv"))
message("All corridor table: ", file.path(table_dir, "corridor_protection_level_all_corridor_municipalities.csv"))
message("Large municipality table: ", file.path(table_dir, "large_municipalities_portfolio_status.csv"))
message("Unmatched table: ", file.path(table_dir, "corridor_municipalities_without_protection.csv"))
message("Joined GeoPackage: ", joined_gpkg)
