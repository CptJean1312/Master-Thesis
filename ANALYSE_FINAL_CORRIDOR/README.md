# ANALYSE FINAL CORRIDOR

This folder contains the corridor-based wide PCA workflow for the revised thesis design.

## Core logic

- Sample: only municipalities intersecting the EFAS-derived `RP500` flood corridor
- Socio-economic method: same wide PCA logic as the original analysis
- Hazard focus for now: `RP100` exposure
- Protection is intentionally excluded at this stage

## Main script

- `ANALYSE_FINAL_CORRIDOR.R`

## Key outputs

- `outputs/gpkg/corridor_wide_pca_rp100_analysis.gpkg`
- `outputs/tables/corridor_analysis_rp100.csv`
- `outputs/tables/corridor_scree_table.csv`
- `outputs/tables/corridor_pca_top_loadings_top8_per_pc.csv`
- `outputs/maps/map_study_area_corridor.png`
- `outputs/maps/map_vulnerability_index_corridor.png`
- `outputs/maps/map_rp100_exposure_corridor.png`

## Notes

- The PCA is recalculated on the corridor sample (`n = 835`), not transferred from the old basin-wide run.
- `RP100` is joined from the new EFAS-based exposure pipeline in `analysev2`.
- This folder is intended as the new working base for the final corridor-based thesis analysis.
