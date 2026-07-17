# DLR Land-Use Exposure QA Report

Generated: 2026-07-17

## Status

The RP100 land-use exposure extension is complete and internally consistent.

## Processing Scope

- Study population: `835` RP500 corridor municipalities.
- Analysis CRS: `EPSG:25832` (ETRS89 / UTM Zone 32N).
- DLR source CRS: `EPSG:3035`.
- DLR raster resolution: `10 x 10 m`.
- Analysis template: corridor-cropped EFAS/JRC flood grid.
- Analysis template resolution: `100 x 100 m`.
- Return period currently processed: `RP100`.

## Expected Outputs

All expected RP100 land-use outputs are present.

- `ANALYSE_FINAL_CORRIDOR/LANDUSE/BUILD_LANDUSE_EXPOSURE.R`
- `ANALYSE_FINAL_CORRIDOR/LANDUSE/BUILD_LANDUSE_MAPS.R`
- `ANALYSE_FINAL_CORRIDOR/LANDUSE/BUILD_LANDUSE_RESULTS.R`
- `outputs/tables/corridor_landuse_exposure_wide.csv`
- `outputs/gpkg/corridor_landuse_exposure.gpkg`
- `outputs/tables/corridor_landuse_protection_analysis_rp100.csv`
- `outputs/gpkg/corridor_landuse_protection_analysis_rp100.gpkg`
- `outputs/tables/table_landuse_exposure_correlation_comparison.csv`
- `outputs/tables/table_artificial_exposure_by_vulnerability_quintile.csv`
- `outputs/tables/table_high_exposure_overlap.csv`
- `outputs/tables/table_protection_by_artificial_exposure_tercile.csv`
- `outputs/maps/map_rp100_artificial_land_flood_share_context.png`
- `outputs/maps/plot_rp100_flooded_landcover_composition.png`
- `outputs/plots/plot_rp100_artificial_land_exposure_by_vulnerability_quintile.png`
- `outputs/plots/map_artificial_exposure_vulnerability_bivariate.png`

## Structural Checks

- Wide land-use table rows: `835`.
- Wide land-use table unique AGS: `835`.
- Land-use/protection analysis table rows: `835`.
- Land-use/protection analysis table unique AGS: `835`.
- Land-use GPKG rows: `835`.
- Land-use GPKG CRS: `EPSG:25832`.
- Missing AGS: `0`.
- Missing protection join in final land-use/protection table: `0`.
- Missing vulnerability index: `1` municipality.
- Missing artificial-land flood share: `1` municipality, because it has no DLR Artificial Land area.
- Bad artificial-land flood shares outside `[0, 1]`: `0`.
- Wide-table `.x`/`.y` suffix columns after cleanup: `0`.

## Row-Count Checks

- Class-level land-cover rows: `5,845`.
  - Expected: `835 municipalities x 7 DLR classes`.
- Group-level land-cover rows: `3,340`.
  - Expected: `835 municipalities x 4 thesis groups`.
- Maximum difference between grouped flooded-area sums and stored flooded land-cover totals: about `1.49e-08 m²`.

## Area Checks

The land-cover area ratio to municipality area ranges from `0.9999985` to `1.000001`.

This means the DLR class areas, after reprojection and aggregation to the corridor-cropped EFAS grid in `EPSG:25832`, reproduce municipal area essentially exactly at the precision expected from raster/polygon operations.

## RP100 Key Results

Land-cover composition of modeled RP100 flooded area:

- Artificial land: `499.4 km²`, `6.0%` of all flooded land-cover area.
- Open/seasonal land cover: `3,883.8 km²`, `46.5%`.
- Perennial vegetation: `3,289.5 km²`, `39.4%`.
- Water: `682.4 km²`, `8.2%`.

Artificial-land exposure:

- Municipalities with any DLR Artificial Land: `834`.
- Municipalities with positive RP100 artificial-land flooding: `747`.
- Municipalities with at least `10%` artificial-land flooding: `322`.
- Municipalities with at least `25%` artificial-land flooding: `165`.

## Main Diagnostic Findings

The land-use extension does not overturn the existing exposure-vulnerability finding.

- Vulnerability vs RP100 total-area flood share:
  - Pearson `-0.050`.
  - Spearman `-0.050`.
- Vulnerability vs RP100 artificial-land flood share:
  - Pearson `-0.088`.
  - Spearman `-0.063`.

Mean artificial-land flood share by vulnerability quintile:

- Q1 lowest vulnerability: `26.4%`.
- Q2: `17.2%`.
- Q3: `14.9%`.
- Q4: `11.9%`.
- Q5 highest vulnerability: `16.7%`.

The highest-vulnerability quintile is therefore not the corridor-wide maximum in artificial-land exposure. This supports the interpretation that the weak RP100 exposure-vulnerability relationship is not simply an artifact of using total municipal area instead of built/artificial land-cover exposure.

## Protection/Loss Cross-Check

Positive modeled loss share by artificial-land exposure tercile:

- Low artificial-land exposure: `24.5%`.
- Medium artificial-land exposure: `45.0%`.
- High artificial-land exposure: `31.3%`.

This is not monotonic. The provider protection/loss layer remains more socially differentiated by vulnerability than by artificial-land exposure alone.

## Interpretation Limits

- DLR class `Artificial Land` is a built/artificial land-cover proxy, not a cadastral building footprint.
- It includes roads, industrial surfaces, and other artificial surfaces, not only residential settlement.
- The analysis is still municipality-level and cannot identify which households or buildings are flooded.
- The DLR land-cover year does not perfectly match all socio-economic and model years.
- Only RP100 has currently been processed. Multi-RP land-use exposure could be added later if needed, but RP100 is enough for the central exposure-metric robustness check.
