# DLR Land-Use Exposure Module

Generated: 2026-10-01 18:29:09

## Processing summary

- DLR Land Cover DE is read as a 10 m categorical raster in EPSG:3035.
- The analysis template is the corridor-cropped EFAS/JRC flood grid in EPSG:25832.
- Analysis template resolution: `100 x 100 m`.
- Analysis template dimensions: `4601 x 4376 x 1`.
- Analysis template extent: `461900, 899500, 5572300, 6032400`.
- Each DLR class is converted to a binary raster and aggregated to the 100 m EFAS grid with `terra::project(..., method = "average")`.
- Class fractions are multiplied by the EPSG:25832 cell area and summed within corridor municipality polygons.
- Flooded class areas use the same class-area rasters masked by valid RP100 flood-depth cells.

## Inputs

- DLR Land Cover DE raster: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/LANDUSE/data_raw/Land_Cover_DE_2015.tif`
- Corridor municipalities: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/outputs/gpkg/full_municipality_corridor_exposure_projected_area.gpkg`
- Exposure CSV: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/FULL_MUNICIPALITY_PROJECTED_AREA_FINAL_2026-10-01/outputs/tables/full_municipality_exposure_all_RPs_projected_area.csv`
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

- `artificial`: 499.5 km² flooded; 6% of all flooded land-cover area; median group flood share 6.1%.
- `open_or_seasonal`: 3884.8 km² flooded; 46.5% of all flooded land-cover area; median group flood share 8.7%.
- `perennial_vegetation`: 3290.3 km² flooded; 39.4% of all flooded land-cover area; median group flood share 14.4%.
- `water`: 682.5 km² flooded; 8.2% of all flooded land-cover area; median group flood share 51.6%.

## Diagnostic correlations

- `vulnerability_vs_rp100_total_area_flood_share`: Pearson -0.053, Spearman -0.053.
- `vulnerability_vs_rp100_artificial_group_flood_share`: Pearson -0.094, Spearman -0.068.
- `access_adaptive_capacity_vs_rp100_artificial_group_flood_share`: Pearson 0.033, Spearman -0.088.
- `demographic_household_vs_rp100_artificial_group_flood_share`: Pearson -0.022, Spearman 0.041.
- `deprivation_labour_vs_rp100_artificial_group_flood_share`: Pearson -0.162, Spearman -0.016.

## Interpretation note

The `artificial` group is based on DLR class 1, `Artificial Land`. It is used as a built/artificial land-cover exposure proxy, not as a cadastral building-footprint or residential-population exposure measure.
The `open_or_seasonal` and `perennial_vegetation` groups are analytical land-cover proxies and should not be overinterpreted as exact agricultural or natural land-use classes.
