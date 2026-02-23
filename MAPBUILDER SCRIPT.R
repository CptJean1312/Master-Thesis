# =========================
# Elbe Basin Overview Map
# (Basin polygon + rivers + major cities)
# =========================

# Packages
library(sf)
library(ggplot2)
library(dplyr)
library(units)
library(ggspatial)
library(prettymapr)
library(ggplot2)


# ---- Paths (from you) ----
p_berlin    <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Berlin.gpkg"
p_dresden   <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Dresden.gpkg"
p_hamburg   <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Hamburg.gpkg"
p_magdeburg <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/magdeburg.gpkg"

p_basin     <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/ELBE BASIN PRO.gpkg"
p_elbe      <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"
p_rivers    <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/CLIPPED + PROJECTED/am_riverwaterbody_de-basin.gpkg"

# ---- Helper: read first layer of a gpkg safely ----
read_gpkg_first_layer <- function(path) {
  layers <- st_layers(path)$name
  st_read(path, layer = layers[1], quiet = TRUE)
}

# ---- Read data ----
basin <- read_gpkg_first_layer(p_basin)
elbe  <- read_gpkg_first_layer(p_elbe)
rivs  <- read_gpkg_first_layer(p_rivers)

berlin    <- read_gpkg_first_layer(p_berlin)    %>% mutate(city = "Berlin")
hamburg   <- read_gpkg_first_layer(p_hamburg)   %>% mutate(city = "Hamburg")
dresden   <- read_gpkg_first_layer(p_dresden)   %>% mutate(city = "Dresden")
magdeburg <- read_gpkg_first_layer(p_magdeburg) %>% mutate(city = "Magdeburg")

cities <- bind_rows(berlin, hamburg, dresden, magdeburg)

# ---- Make sure everything is in the same CRS ----
target_crs <- st_crs(basin)

elbe   <- st_transform(elbe, target_crs)
rivs   <- st_transform(rivs, target_crs)
cities <- st_transform(cities, target_crs)

# ---- Geometry hygiene (optional but helps if something is weird) ----
basin  <- st_make_valid(basin)
cities <- st_make_valid(cities)

# ---- Create label points (nice: point on surface, works for polygons) ----
# ---- Create label points (robust) ----
city_labels <- st_point_on_surface(cities)

coords <- st_coordinates(st_geometry(city_labels))
city_labels <- city_labels %>%
  mutate(x = coords[, 1],
         y = coords[, 2])

# ---- Define map extent with a little padding around basin ----
bb <- st_bbox(basin)
pad_x <- (bb$xmax - bb$xmin) * 0.08
pad_y <- (bb$ymax - bb$ymin) * 0.08
xlim <- c(bb$xmin - pad_x, bb$xmax + pad_x)
ylim <- c(bb$ymin - pad_y, bb$ymax + pad_y)

# ---- Styling parameters (easy to tweak) ----
line_main   <- 0.8   # Elbe main stem width
line_other  <- 0.35  # tributaries/other river bodies width
alpha_other <- 0.55

# ---- Plot ----
p <- ggplot() +
  # Basin outline + light fill
  geom_sf(data = basin, fill = "grey95", color = "grey40", linewidth = 0.6) +
  
  # River bodies / tributaries
  geom_sf(data = rivs, color = "steelblue4", linewidth = line_other, alpha = alpha_other) +
  
  # Main Elbe highlighted
  geom_sf(data = elbe, color = "steelblue4", linewidth = line_main) +
  
  # Cities polygons
  geom_sf(data = cities, fill = "grey20", color = NA, alpha = 0.75) +
  
  # City labels
  geom_text(
    data = city_labels,
    aes(x = x, y = y, label = city),
    size = 3.6,
    fontface = "bold",
    color = "black",
    nudge_y = (ylim[2] - ylim[1]) * 0.01
  ) +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  
  # North arrow + scale bar
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1.1, "cm"), width = unit(1.1, "cm")) +
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7) +
  
  labs(
    title = "Elbe Basin (Germany/Czechia) – Overview",
    subtitle = "Catchment boundary, Elbe River and major cities",
    caption = "Data: own processing (GPKG layers); Map: R (sf, ggplot2)"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9, color = "grey40"),
    axis.title = element_blank()
  )

p


########################## VARIANT 2

# =========================
# Elbe Basin Overview Map
# Basemap via maptiles (no prettymapr needed)
# =========================

library(sf)
library(ggplot2)
library(dplyr)
library(units)
library(ggspatial)
library(maptiles)   # <-- install if missing
library(terra)      # maptiles uses this
library(tidyterra)



# ---- Paths ----
p_berlin    <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Berlin.gpkg"
p_dresden   <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Dresden.gpkg"
p_hamburg   <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Hamburg.gpkg"
p_magdeburg <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/magdeburg.gpkg"
p_erfurt <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Erfurt.gpkg" 
p_leipzig <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/MAP FEATURES/CITIES/Leipzig.gpkg"

p_basin     <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/ELBE BASIN PRO.gpkg"
p_elbe      <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"
p_rivers    <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/CLIPPED + PROJECTED/am_riverwaterbody_de-basin.gpkg"

read_gpkg_first_layer <- function(path) {
  layers <- st_layers(path)$name
  st_read(path, layer = layers[1], quiet = TRUE)
}

# ---- Read ----
basin <- read_gpkg_first_layer(p_basin)
elbe  <- read_gpkg_first_layer(p_elbe)
rivs  <- read_gpkg_first_layer(p_rivers)

berlin    <- read_gpkg_first_layer(p_berlin)    %>% mutate(city = "Berlin")
hamburg   <- read_gpkg_first_layer(p_hamburg)   %>% mutate(city = "Hamburg")
dresden   <- read_gpkg_first_layer(p_dresden)   %>% mutate(city = "Dresden")
magdeburg <- read_gpkg_first_layer(p_magdeburg) %>% mutate(city = "Magdeburg")
erfurt <- read_gpkg_first_layer(p_erfurt) %>% mutate(city = "Erfurt")
leipzig <- read_gpkg_first_layer(p_leipzig) %>% mutate(city = "Leipzig")


cities <- bind_rows(berlin, hamburg, dresden, magdeburg, erfurt, leipzig)

# ---- Valid + CRS for basemap ----
basin  <- st_make_valid(basin)
cities <- st_make_valid(cities)

target_crs <- 3857
basin  <- st_transform(basin, target_crs)
elbe   <- st_transform(elbe, target_crs)
rivs   <- st_transform(rivs, target_crs)
cities <- st_transform(cities, target_crs)

# ---- Labels ----
city_labels <- suppressWarnings(st_point_on_surface(cities))
coords <- st_coordinates(st_geometry(city_labels))
city_labels <- city_labels %>% mutate(x = coords[,1], y = coords[,2])

# ---- Extent + mask ----
basin_union <- st_union(st_make_valid(basin))
bb <- st_bbox(basin_union)

xrange <- as.numeric(c(bb["xmin"], bb["xmax"]))
yrange <- as.numeric(c(bb["ymin"], bb["ymax"]))
pad_x <- diff(xrange) * 0.08
pad_y <- diff(yrange) * 0.08

xlim <- c(xrange[1] - pad_x, xrange[2] + pad_x)
ylim <- c(yrange[1] - pad_y, yrange[2] + pad_y)

bbox_poly <- st_as_sfc(st_bbox(c(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]),
                               crs = st_crs(basin)))
outside <- st_difference(st_make_valid(bbox_poly), st_make_valid(basin_union)) %>% st_as_sf()

# ---- Basemap tiles (choose provider) ----
# Nice clean B/W: "CartoDB.Positron"
# Topo-ish:       "Esri.WorldTopoMap"
basemap <- tryCatch(
  {
    get_tiles(bbox_poly, provider = "CartoDB.Positron", crop = TRUE)
  },
  error = function(e) {
    message("Basemap download failed -> plotting without basemap. Reason: ", e$message)
    NULL
  }
)
# ---- Style ----
line_main   <- 0.9
line_other  <- 0.35
alpha_other <- 0.55

city_fill  <- "#d95f02"
label_text <- "#1f3b73"
label_bg   <- "white"

# ---- Plot ----
if (!is.null(basemap)) {
  p <- p + layer_spatial(basemap)
}

p2 <- ggplot() +
  layer_spatial(basemap) +
  
  # wash out outside basin
  geom_sf(data = outside, fill = "white", color = NA, alpha = 0.55) +
  
  # basin outline
  geom_sf(data = basin, fill = "white", alpha = 0.06, color = "grey20", linewidth = 0.8) +
  
  # rivers
  geom_sf(data = rivs, color = "steelblue4", linewidth = line_other, alpha = alpha_other) +
  geom_sf(data = elbe, color = "steelblue4", linewidth = line_main) +
  
  # cities
  geom_sf(data = cities, fill = city_fill, color = NA, alpha = 0.75) +
  
  # labels
  geom_label(
    data = city_labels,
    aes(x = x, y = y, label = city),
    size = 3.6, fontface = "bold",
    color = label_text, fill = label_bg,
    label.size = 0.2,
    label.r = unit(0.15, "lines"),
    alpha = 0.95
  ) +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1.1, "cm"), width = unit(1.1, "cm")) +
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.75) +
  
  labs(
    title = "Elbe Basin – Overview",
    subtitle = "Catchment boundary, Elbe river network and major cities",
    caption = "Basemap: OpenStreetMap; Data: own processing (GPKG layers)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9, color = "grey35")
  )


p2

