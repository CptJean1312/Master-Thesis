# Exposure Curves Overview

Created: `2026-03-20 13:21:25`

## Purpose

This output prepares municipality-level exposure curves across return periods so that later protection-related interpretations can build on an already structured set of curve metrics.

## Input

- Exposure table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/municipality_flood_exposure_all_RPs.csv`
- Quality checks: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/tables/exposure_quality_checks.csv`
- Corridor geometry: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/analysev2/outputs_exposure_pipeline/corridor/municipalities_corridor.gpkg`

## Metric definitions

- `onset_rp_any`: first return period with any positive flood share
- `onset_rp_1pct`: first return period with flood share >= 1% of municipality area
- `onset_rp_5pct`: first return period with flood share >= 5% of municipality area
- `slope_logrp`: slope of monotone-adjusted flood share over log10(return period)
- `normalized_auc`: area under the normalized exposure curve, scaled to `[0,1]`; higher values indicate earlier realization of final exposure
- `realized_by_rp100`: share of final RP500 exposure already realized by RP100
- `late_growth_share`: share of final RP500 exposure added after RP100
- `max_jump_interval`: interval with the largest increase in flood share
- `max_jump_rel`: largest interval increase relative to final RP500 exposure

## Treatment of non-monotonic cases

- Observed non-monotonic municipalities repaired with a simple `cummax` monotone adjustment: `2`
- The adjusted curve is used for all derived curve metrics.

## Typology rules

- `early exposure`: at least 75% of final RP500 exposure already realized by RP100, or onset <= RP50 together with substantial early realization
- `delayed jump`: less than 50% of final RP500 exposure realized by RP100, with at least 40% of final exposure added after RP100 and the largest jump in RP100-RP200 or RP200-RP500
- `gradual increase`: all remaining municipalities

## Output files

- Metrics table: `outputs/tables/corridor_exposure_curve_metrics.csv`
- Long curve table: `outputs/tables/corridor_exposure_curves_long.csv`
- Type summary: `outputs/tables/corridor_exposure_curve_type_summary.csv`
- Onset summary: `outputs/tables/corridor_onset_1pct_summary.csv`
- Geometry output: `outputs/gpkg/corridor_exposure_curve_metrics.gpkg`
- Mean normalized curves plot: `outputs/plots/plot_mean_normalized_curves_by_type.png`
- Type-count plot: `outputs/plots/plot_exposure_curve_type_counts.png`
- Type map: `outputs/plots/map_exposure_curve_types.png`
