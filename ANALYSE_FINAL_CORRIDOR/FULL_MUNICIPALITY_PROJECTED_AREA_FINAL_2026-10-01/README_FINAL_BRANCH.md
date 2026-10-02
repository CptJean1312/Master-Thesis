# Full Municipality Projected-Area Final Branch

## Status

This directory contains the authoritative post-audit municipality-level outputs created on 1 October 2026. It does not replace or delete any earlier output. The former analysis, the initial full-geometry audit, and the intermediate corrected branch remain available for comparison.

## Why This Branch Exists

The earlier exposure workflow selected the correct 835 RP500 corridor municipalities but retained basin-clipped geometries for 51 municipalities crossing the Elbe processing-frame boundary. The final branch restores complete BKG VG250 municipality polygons before calculating municipality areas and joining municipality-level social data.

Flooded area is calculated on the projected 100 m by 100 m `EPSG:25832` grid with exact polygon coverage. This avoids combining `terra::cellSize()` geodetic cell areas with planar `sf` municipality areas.

## Authoritative Outputs

- Exposure: `outputs/tables/full_municipality_exposure_all_RPs_projected_area.csv`
- Exposure QA: `outputs/tables/exposure_quality_summary_projected_area_with_numeric_tolerance.csv`
- Main corridor analysis: `outputs/tables/corridor_analysis_rp100.csv`
- Main spatial layer: `outputs/gpkg/corridor_wide_pca_rp100_analysis.gpkg`
- Land-cover refinement: `LANDUSE/outputs/tables/corridor_landuse_exposure_wide.csv`
- Land-cover grouped summary: `LANDUSE/outputs/tables/corridor_landuse_exposure_summary_by_group.csv`
- Exposure-curve metrics: `EXPOSURE_CURVES/outputs/tables/corridor_exposure_curve_metrics.csv`
- Modelled-loss analysis: `MODELLED_LOSS/outputs/tables/corridor_modelled_loss_analysis.csv`
- Modelled-loss correlations: `MODELLED_LOSS/outputs/tables/modelled_loss_correlations.csv`

## Confirmed Sample Sizes

- RP500 corridor municipalities: 835
- Municipalities with total-area exposure: 835
- Municipalities with Artificial Land exposure: 835
- Municipalities with matched INKAR vulnerability data: 834
- Positive modelled-loss municipalities in the corridor: 280
- No-event municipalities in the 5,000-year loss catalogue: 555

## Quality Checks

- All exposure shares pass the zero-to-one check with a numerical tolerance of `1e-8`.
- All flooded areas pass the municipality-area check with a numerical tolerance of 1 m2.
- Two small non-monotonic sequences remain in the raw exposure data: Plau am See and Fockendorf.
- The raw values are preserved. A cumulative maximum is applied only when calculating curve metrics.
- The 5,000-year modelled-loss probability and return-period formulas are reproduced exactly for all 301 positive-loss rows in the supplied extract.

## Interpretation Boundary

The JRC rasters measure modelled river-flood extent without physical flood defences. The RFM extract measures modelled non-zero loss occurrence and frequency. Neither dataset is an inventory of flood-protection infrastructure. The source field named `protection_return_period` is renamed in this branch as `finite_loss_occurrence_return_period` because it equals `5,000 / n_nonzero_years`.

## Figure Preservation

Any figures regenerated from this branch must be written to a new dated figure folder. Existing thesis figures are archived versions and must not be overwritten or deleted.
