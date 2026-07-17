#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

options(scipen = 999)

# Thesis-facing maps and summary graphics for the DLR land-use exposure module.
# The analytical land-use exposure data are built in BUILD_LANDUSE_EXPOSURE.R.
# This script only handles presentation maps with orientation layers.

root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis"
data_root <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS"

paths <- list(
  landuse_gpkg = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/gpkg/corridor_landuse_exposure.gpkg"),
  landuse_layer = "corridor_landuse_exposure",
  summary_by_group = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/tables/corridor_landuse_exposure_summary_by_group.csv"),
  elbe_gpkg = file.path(data_root, "PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"),
  elbe_layer = "elbe",
  state_dir = file.path(data_root, "PHYSISCH.nosync/RIGHT PROJECTION"),
  city_dir = file.path(data_root, "MAP FEATURES/CITIES"),
  output_dir = file.path(root, "ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/maps")
)

dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)

read_first_layer <- function(path) {
  st_read(path, layer = st_layers(path)$name[1], quiet = TRUE)
}

read_optional_gpkg <- function(path, label) {
  if (!file.exists(path)) {
    warning(label, " not found: ", path, call. = FALSE)
    return(NULL)
  }
  read_first_layer(path)
}

save_map <- function(plot, filename, width = 11, height = 8.5) {
  ggsave(
    filename = file.path(paths$output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white",
    limitsize = FALSE
  )
}

map_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_line(color = "grey90", linewidth = 0.2),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = base_size + 4),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 0),
      strip.text = element_text(face = "bold")
    )
}

plot_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = base_size + 4),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 0)
    )
}

message("Loading land-use exposure layer ...")
landuse <- st_read(paths$landuse_gpkg, layer = paths$landuse_layer, quiet = TRUE) %>%
  mutate(AGS = as.character(AGS))

if (st_crs(landuse)$epsg != 25832) {
  landuse <- st_transform(landuse, 25832)
}

message("Loading Elbe geometry ...")
elbe <- st_read(paths$elbe_gpkg, layer = paths$elbe_layer, quiet = TRUE)
if (st_crs(elbe)$epsg != 25832) {
  elbe <- st_transform(elbe, 25832)
}

message("Loading federal-state orientation layers ...")
state_files <- c(
  "Schleswig-Holstein.gpkg",
  "Niedersachsen.gpkg",
  "Hamburg.gpkg",
  "Bremen.gpkg",
  "Mecklenburg-Vorpommern.gpkg",
  "Brandenburg.gpkg",
  "Berlin.gpkg",
  "Sachsen-Anhalt.gpkg",
  "Sachsen.gpkg",
  "Thüringen.gpkg"
)

states <- bind_rows(lapply(state_files, function(file_name) {
  path <- file.path(paths$state_dir, file_name)
  state <- read_optional_gpkg(path, file_name)
  if (is.null(state)) {
    return(NULL)
  }
  state %>%
    st_transform(25832) %>%
    mutate(state_name = sub("\\.gpkg$", "", file_name))
}))

message("Loading city orientation layers ...")
city_files <- c(
  Hamburg = "Hamburg.gpkg",
  Magdeburg = "magdeburg.gpkg",
  Leipzig = "Leipzig.gpkg",
  Dresden = "Dresden.gpkg",
  Berlin = "Berlin.gpkg"
)

cities <- bind_rows(lapply(names(city_files), function(city_name) {
  path <- file.path(paths$city_dir, city_files[[city_name]])
  city <- read_optional_gpkg(path, city_name)
  if (is.null(city)) {
    return(NULL)
  }
  city %>%
    st_transform(25832) %>%
    st_make_valid() %>%
    mutate(city = city_name)
}))

city_points <- suppressWarnings(st_point_on_surface(cities))
city_xy <- st_coordinates(city_points)
city_points <- city_points %>%
  mutate(
    x = city_xy[, 1],
    y = city_xy[, 2],
    x_label = case_when(
      city == "Hamburg" ~ x - 19000,
      city == "Magdeburg" ~ x + 17000,
      city == "Leipzig" ~ x + 15000,
      city == "Dresden" ~ x + 15000,
      city == "Berlin" ~ x + 17000,
      TRUE ~ x + 12000
    ),
    y_label = case_when(
      city == "Hamburg" ~ y + 13000,
      city == "Magdeburg" ~ y - 7000,
      city == "Leipzig" ~ y - 10000,
      city == "Dresden" ~ y - 8000,
      city == "Berlin" ~ y + 11000,
      TRUE ~ y
    )
  )

state_points <- suppressWarnings(st_point_on_surface(states))
state_xy <- st_coordinates(state_points)
state_points <- state_points %>%
  mutate(
    x = state_xy[, 1],
    y = state_xy[, 2],
    label = recode(
      state_name,
      "Schleswig-Holstein" = "Schleswig-\nHolstein",
      "Mecklenburg-Vorpommern" = "Mecklenburg-\nVorpommern",
      "Sachsen-Anhalt" = "Sachsen-\nAnhalt",
      .default = state_name
    )
  )

summary_by_group <- read_csv(paths$summary_by_group, show_col_types = FALSE) %>%
  mutate(
    thesis_group_label = recode(
      thesis_group,
      artificial = "Artificial land",
      open_or_seasonal = "Open/seasonal land cover",
      perennial_vegetation = "Perennial vegetation",
      water = "Water"
    ),
    share_pct = share_of_all_flooded_landcover_area * 100,
    flooded_area_km2 = flooded_area_m2 / 1e6
  )

map_bbox <- st_bbox(landuse)
xmin <- unname(map_bbox[["xmin"]])
ymin <- unname(map_bbox[["ymin"]])
xmax <- unname(map_bbox[["xmax"]])
ymax <- unname(map_bbox[["ymax"]])
pad_x <- as.numeric(xmax - xmin) * 0.035
pad_y <- as.numeric(ymax - ymin) * 0.035

states_context <- states %>%
  st_crop(st_as_sfc(st_bbox(c(
    xmin = xmin - 80000,
    ymin = ymin - 80000,
    xmax = xmax + 80000,
    ymax = ymax + 80000
  ), crs = st_crs(landuse))))

elbe_context <- st_crop(elbe, st_as_sfc(st_bbox(c(
  xmin = xmin - 50000,
  ymin = ymin - 50000,
  xmax = xmax + 50000,
  ymax = ymax + 50000
), crs = st_crs(landuse))))

city_points_context <- city_points %>%
  filter(
    x >= xmin - 50000,
    x <= xmax + 50000,
    y >= ymin - 50000,
    y <= ymax + 50000
  )

state_points_context <- state_points %>%
  filter(
    x >= xmin - 90000,
    x <= xmax + 90000,
    y >= ymin - 90000,
    y <= ymax + 90000
  )

map_artificial <- ggplot() +
  geom_sf(data = states_context, fill = "#f6f5f0", color = "#b9b4a7", linewidth = 0.25) +
  geom_sf(data = landuse, aes(fill = rp100_artificial_flood_share_of_group), color = "white", linewidth = 0.035) +
  geom_sf(data = elbe_context, color = "#0b4f8a", linewidth = 0.6, inherit.aes = FALSE) +
  geom_sf(data = city_points_context, shape = 21, fill = "#202020", color = "white", linewidth = 0.18, size = 2.2, inherit.aes = FALSE) +
  geom_segment(
    data = st_drop_geometry(city_points_context),
    aes(x = x, y = y, xend = x_label, yend = y_label),
    linewidth = 0.18,
    color = "#333333"
  ) +
  geom_label(
    data = st_drop_geometry(city_points_context),
    aes(x = x_label, y = y_label, label = city),
    size = 2.65,
    linewidth = 0.12,
    label.padding = unit(0.11, "lines"),
    fill = "white",
    color = "#202020"
  ) +
  geom_text(
    data = st_drop_geometry(state_points_context),
    aes(x = x, y = y, label = label),
    size = 2.35,
    color = "#6b655b",
    alpha = 0.78,
    lineheight = 0.92
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    labels = percent_format(accuracy = 1),
    na.value = "#e5e5e5",
    name = "Flooded share\nof artificial land"
  ) +
  coord_sf(
    xlim = c(xmin - pad_x, xmax + pad_x),
    ylim = c(ymin - pad_y, ymax + pad_y),
    expand = FALSE
  ) +
  labs(
    title = "RP100 artificial-land exposure in the Elbe corridor",
    subtitle = "Share of DLR Artificial Land intersecting modeled RP100 flooding, by corridor municipality",
    caption = "Data: DLR Land Cover DE 2015, EFAS/JRC RP100 flood raster, VG250 municipalities. Blue line: Elbe. Own processing."
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(7.8, "cm"), barheight = unit(0.35, "cm"))) +
  map_theme()

save_map(map_artificial, "map_rp100_artificial_land_flood_share_context.png")

composition_plot <- ggplot(summary_by_group, aes(x = reorder(thesis_group_label, share_pct), y = share_pct)) +
  geom_col(fill = "#2f6f73", width = 0.72) +
  geom_text(
    aes(label = paste0(round(share_pct, 1), "%")),
    hjust = -0.12,
    size = 3.4,
    color = "#202020"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, max(summary_by_group$share_pct) * 1.18),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Land-cover composition of modeled RP100 flooded area",
    subtitle = "Municipality-level DLR class areas aggregated to four thesis groups",
    x = NULL,
    y = "Share of all flooded land-cover area",
    caption = "Data: DLR Land Cover DE 2015 and EFAS/JRC RP100 flood raster. Own processing."
  ) +
  plot_theme()

save_map(composition_plot, "plot_rp100_flooded_landcover_composition.png", width = 8.5, height = 5.5)

message("Done.")
message("Context map: ", file.path(paths$output_dir, "map_rp100_artificial_land_flood_share_context.png"))
message("Composition plot: ", file.path(paths$output_dir, "plot_rp100_flooded_landcover_composition.png"))
