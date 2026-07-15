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
# Check stream context of corridor municipalities without
# simulated loss events in the provider flood portfolio
# ------------------------------------------------------------
# Purpose:
# Test whether municipalities absent from elbe_protection_level_mun.csv
# are mainly located away from the Elbe main stem and associated with
# smaller / other WFD river water bodies.
# ============================================================

# ---------------------------
# 1) Input paths
# ---------------------------

corridor_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg"
corridor_layer <- "municipalities_corridor"

exposure_csv <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv"

protection_csv <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PROTECTION/elbe_protection_level_mun.csv"

elbe_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"
elbe_layer <- "elbe"

riverwaterbody_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/CLIPPED + PROJECTED/am_riverwaterbody_de-basin.gpkg"
riverwaterbody_layer <- "am_riverwaterbody_de-basin"

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
# 3) Load data
# ---------------------------

message("Loading corridor municipalities ...")
corridor <- st_read(corridor_path, layer = corridor_layer, quiet = TRUE) %>%
  mutate(
    AGS = coalesce(as.character(AGS), as.character(Gemeindeschlüssel_AGS)),
    mun_name = coalesce(mun_name, GeografischerName_GEN)
  )

if (st_crs(corridor)$epsg != 25832) {
  corridor <- st_transform(corridor, 25832)
}

message("Loading exposure table ...")
exposure <- read_csv(
  exposure_csv,
  col_types = cols(AGS = col_character()),
  show_col_types = FALSE
) %>%
  select(
    AGS,
    starts_with("flood_area_"),
    starts_with("flood_share_")
  )

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

message("Loading Elbe main-stem geometry ...")
elbe <- st_read(elbe_path, layer = elbe_layer, quiet = TRUE)
if (st_crs(elbe)$epsg != 25832) {
  elbe <- st_transform(elbe, 25832)
}

message("Loading WFD river water bodies ...")
rivers <- st_read(riverwaterbody_path, layer = riverwaterbody_layer, quiet = TRUE)
if (st_crs(rivers)$epsg != 25832) {
  rivers <- st_transform(rivers, 25832)
}

# ---------------------------
# 4) Prepare joined municipality dataset
# ---------------------------

mun <- corridor %>%
  left_join(exposure, by = "AGS") %>%
  left_join(protection, by = c("AGS" = "ags")) %>%
  mutate(
    protection_available = !is.na(protection_return_period),
    loss_portfolio_status = if_else(
      protection_available,
      "positive_modeled_loss",
      "zero_modeled_loss_probability"
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

message("Corridor municipalities: ", nrow(mun))
message("Protection values available: ", sum(mun$protection_available))
message("No simulated loss event / zero modeled loss probability: ", sum(!mun$protection_available))

# ---------------------------
# 5) Classify WFD river-body context
# ---------------------------

# Conservative operationalization:
# - Elbe main stem: municipality lies directly on or within 1 km of the
#   explicit Elbe line layer.
# - Major WFD river/tributary: WFD river body has RIVER_CAT 1 or 100,
#   or map drawing category 1/2. This captures Elbe-like main channels
#   and larger tributaries such as Saale, Havel, Mulde, Bode etc.
# - Minor/other WFD: municipality intersects only other WFD river bodies.

rivers <- rivers %>%
  mutate(
    river_name = coalesce(S_NAME, NAMETEXT, TEXTINTERN, as.character(RIVER_CD)),
    is_major_wfd = RIVER_CAT %in% c(1, 100) |
      DRAW_CAT %in% c("N1", "M1", "A1", "N2", "M2", "A2")
  )

major_rivers <- rivers %>%
  filter(is_major_wfd)

elbe_union <- st_union(elbe)
distance_to_elbe_m <- as.numeric(st_distance(mun, elbe_union))

any_river_idx <- st_intersects(mun, rivers)
major_river_idx <- st_intersects(mun, major_rivers)

collapse_names <- function(idx, source, max_names = 8) {
  if (length(idx) == 0) {
    return(NA_character_)
  }
  values <- unique(na.omit(source$river_name[idx]))
  values <- values[values != ""]
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(head(values, max_names), collapse = "; ")
}

mun <- mun %>%
  mutate(
    distance_to_elbe_m = distance_to_elbe_m,
    intersects_elbe_main = distance_to_elbe_m < 1e-6,
    near_elbe_1km = distance_to_elbe_m <= 1000,
    near_elbe_5km = distance_to_elbe_m <= 5000,
    intersects_any_wfd_river = lengths(any_river_idx) > 0,
    intersects_major_wfd_river = lengths(major_river_idx) > 0,
    river_context_class = case_when(
      near_elbe_1km ~ "directly on/within 1 km of Elbe main stem",
      intersects_major_wfd_river ~ "intersects major WFD river/tributary",
      intersects_any_wfd_river ~ "only minor/other WFD water bodies",
      TRUE ~ "no WFD river body intersects municipality"
    ),
    major_wfd_names = vapply(
      major_river_idx,
      collapse_names,
      source = major_rivers,
      FUN.VALUE = character(1)
    ),
    any_wfd_names = vapply(
      any_river_idx,
      collapse_names,
      source = rivers,
      FUN.VALUE = character(1)
    )
  )

# ---------------------------
# 6) Export tables
# ---------------------------

stream_context_table <- mun %>%
  st_drop_geometry() %>%
  transmute(
    AGS,
    municipality_name = municipality_name_final,
    protection_available,
    loss_portfolio_status,
    river_context_class,
    distance_to_elbe_m,
    intersects_elbe_main,
    near_elbe_1km,
    near_elbe_5km,
    intersects_major_wfd_river,
    intersects_any_wfd_river,
    major_wfd_names,
    any_wfd_names,
    flood_share_rp10,
    flood_share_rp20,
    flood_share_rp50,
    flood_share_rp100,
    flood_share_rp200,
    flood_share_rp500,
    protection_return_period,
    annual_loss_probability,
    annual_loss_probability_model,
    no_simulated_loss_event,
    n_nonzero_years
  ) %>%
  arrange(protection_available, river_context_class, AGS)

summary_by_status <- stream_context_table %>%
  count(loss_portfolio_status, protection_available, river_context_class, name = "n") %>%
  group_by(loss_portfolio_status) %>%
  mutate(share_within_status = n / sum(n)) %>%
  ungroup() %>%
  arrange(loss_portfolio_status, desc(n))

no_loss_summary <- stream_context_table %>%
  filter(!protection_available) %>%
  group_by(river_context_class) %>%
  summarise(
    n = n(),
    share = n / nrow(filter(stream_context_table, !protection_available)),
    mean_distance_to_elbe_m = mean(distance_to_elbe_m, na.rm = TRUE),
    median_distance_to_elbe_m = median(distance_to_elbe_m, na.rm = TRUE),
    mean_flood_share_rp100 = mean(flood_share_rp100, na.rm = TRUE),
    median_flood_share_rp100 = median(flood_share_rp100, na.rm = TRUE),
    mean_flood_share_rp500 = mean(flood_share_rp500, na.rm = TRUE),
    median_flood_share_rp500 = median(flood_share_rp500, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

write_csv(
  stream_context_table,
  file.path(table_dir, "corridor_stream_context_by_municipality.csv")
)

write_csv(
  summary_by_status,
  file.path(table_dir, "corridor_stream_context_summary_by_portfolio_status.csv")
)

write_csv(
  no_loss_summary,
  file.path(table_dir, "no_loss_stream_context_summary.csv")
)

stream_gpkg <- file.path(gpkg_dir, "corridor_stream_context.gpkg")
if (file.exists(stream_gpkg)) {
  unlink(stream_gpkg)
}

st_write(
  mun,
  dsn = stream_gpkg,
  layer = "corridor_stream_context",
  quiet = TRUE
)

# ---------------------------
# 7) Map
# ---------------------------

map_bbox <- st_bbox(mun)

major_river_for_map <- major_rivers %>%
  st_crop(st_as_sfc(map_bbox))

no_loss_map_data <- mun %>%
  filter(!protection_available)

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

city_labels <- mun %>%
  mutate(
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

context_colors <- c(
  "directly on/within 1 km of Elbe main stem" = "#54278f",
  "intersects major WFD river/tributary" = "#238b45",
  "only minor/other WFD water bodies" = "#fdae6b",
  "no WFD river body intersects municipality" = "#969696"
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

stream_context_map <- ggplot() +
  geom_sf(data = mun, fill = "#eeeeea", color = "white", linewidth = 0.04) +
  geom_sf(
    data = major_river_for_map,
    color = "#9ecae1",
    linewidth = 0.18,
    alpha = 0.45
  ) +
  geom_sf(
    data = no_loss_map_data,
    aes(fill = river_context_class),
    color = "white",
    linewidth = 0.05
  ) +
  geom_sf(
    data = elbe,
    color = "#08306b",
    linewidth = 0.65,
    alpha = 0.95
  ) +
  geom_sf(
    data = city_labels,
    aes(color = city_portfolio_status),
    size = 2.0,
    alpha = 0.95,
    inherit.aes = FALSE
  ) +
  geom_sf_label(
    data = city_labels,
    aes(label = city_label, color = city_portfolio_status),
    fill = "white",
    linewidth = 0.12,
    label.padding = unit(0.08, "lines"),
    size = 2.35,
    fontface = "bold",
    alpha = 0.9,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = context_colors,
    labels = c(
      "directly on/within 1 km of Elbe main stem" = "Elbe main stem (<=1 km)",
      "intersects major WFD river/tributary" = "Major WFD river/tributary",
      "only minor/other WFD water bodies" = "Minor/other WFD water bodies",
      "no WFD river body intersects municipality" = "No WFD river body"
    ),
    drop = FALSE,
    name = "Context"
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
    title = "No-event corridor municipalities by stream context",
    subtitle = "Colored municipalities have no simulated loss event in the provider output. Pale blue: major WFD rivers; dark blue: Elbe.",
    caption = "Classification is a GIS diagnostic. Smaller streams are not represented in the provider simulation."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave(
  filename = file.path(map_dir, "map_no_loss_stream_context.png"),
  plot = stream_context_map,
  width = 12,
  height = 9,
  dpi = 320,
  bg = "white"
)

# ---------------------------
# 8) Console summary
# ---------------------------

message("Done.")
message("No-event summary:")
print(no_loss_summary)
message("Summary by portfolio status:")
print(summary_by_status)
message("Municipality table: ", file.path(table_dir, "corridor_stream_context_by_municipality.csv"))
message("No-event summary table: ", file.path(table_dir, "no_loss_stream_context_summary.csv"))
message("Map: ", file.path(map_dir, "map_no_loss_stream_context.png"))
message("GeoPackage: ", stream_gpkg)
