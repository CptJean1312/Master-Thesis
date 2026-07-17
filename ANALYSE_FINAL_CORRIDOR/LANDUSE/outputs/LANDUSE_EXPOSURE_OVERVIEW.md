# DLR Land-Use Exposure Module

Generated: 2026-07-17 12:02:47

## Processing summary

- DLR Land Cover DE is read as a 10 m categorical raster in EPSG:3035.
- The analysis template is the corridor-cropped EFAS/JRC flood grid in EPSG:25832.
- Analysis template resolution: `100 x 100 m`.
- Analysis template dimensions: `4568 x 4373 x 1`.
- Analysis template extent: `461900, 899200, 5572300, 6029100`.
- Each DLR class is converted to a binary raster and aggregated to the 100 m EFAS grid with `terra::project(..., method = "average")`.
- Class fractions are multiplied by the EPSG:25832 cell area and summed within corridor municipality polygons.
- Flooded class areas use the same class-area rasters masked by valid RP100 flood-depth cells.

## Inputs

- DLR Land Cover DE raster: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/LANDUSE/data_raw/Land_Cover_DE_2015.tif`
- Corridor municipalities: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg`
- Exposure CSV: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv`
- Flood raster directory: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_eu_flood_25832`
- Return periods processed: `rp100`

## DLR classes

- `1`: Artificial Land -> `artificial`
- `2`: Open Soil -> `open_or_seasonal`
- `3`: High Seasonal Vegetation -> `open_or_seasonal`
- `4`: High Perennial Vegetation -> `perennial_vegetation`
- `5`: Low Seasonal Vegetation -> `open_or_seasonal`
- `6`: Low Perennial Vegetation -> `perennial_vegetation`
- `7`: Water -> `water`

## Main outputs

- `outputs/tables/corridor_landcover_area_by_class_long.csv`
- `outputs/tables/corridor_flooded_landcover_by_class_long.csv`
- `outputs/tables/corridor_landcover_area_by_group_long.csv`
- `outputs/tables/corridor_flooded_landcover_by_group_long.csv`
- `outputs/tables/corridor_landuse_exposure_wide.csv`
- `outputs/tables/corridor_landuse_exposure_diagnostic_correlations.csv`
- `outputs/gpkg/corridor_landuse_exposure.gpkg`

## RP100 key results

- `artificial`: 499.4 km² flooded; 6% of all flooded land-cover area; median group flood share 6.2%.
- `open_or_seasonal`: 3883.8 km² flooded; 46.5% of all flooded land-cover area; median group flood share 8.7%.
- `perennial_vegetation`: 3289.5 km² flooded; 39.4% of all flooded land-cover area; median group flood share 14.7%.
- `water`: 682.4 km² flooded; 8.2% of all flooded land-cover area; median group flood share 52.3%.

## Diagnostic correlations

- `vulnerability_vs_rp100_total_area_flood_share`: Pearson -0.05, Spearman -0.05.
- `vulnerability_vs_rp100_artificial_group_flood_share`: Pearson -0.088, Spearman -0.063.
- `access_adaptive_capacity_vs_rp100_artificial_group_flood_share`: Pearson 0.039, Spearman -0.088.
- `demographic_household_vs_rp100_artificial_group_flood_share`: Pearson -0.015, Spearman 0.046.
- `deprivation_labour_vs_rp100_artificial_group_flood_share`: Pearson -0.165, Spearman -0.016.

## Interpretation note

The `artificial` group is based on DLR class 1, `Artificial Land`. It is used as a built/artificial land-cover exposure proxy, not as a cadastral building-footprint or residential-population exposure measure.
The `open_or_seasonal` and `perennial_vegetation` groups are analytical land-cover proxies and should not be overinterpreted as exact agricultural or natural land-use classes.
