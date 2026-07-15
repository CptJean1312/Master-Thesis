# Thesis Draft Audit

Date: 2026-07-15

Active draft:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/THESIS_DRAFTS/THESIS DRAFT.md`

## Current Core Status

The thesis is now coherent at the level of data logic:

- study area: RP500 Elbe flood corridor
- municipalities in corridor: `835`
- vulnerability index available for: `834`
- missing INKAR municipality: `16076094` Berga-Wuenschendorf
- main exposure benchmark: `RP100` flood share
- final vulnerability index: 23 INKAR indicators, direction-coded, dimension-balanced
- PCA: diagnostic and robustness tool, not the main index
- protection/loss portfolio: `280` positive-loss municipalities with finite protection return periods; `555` no-event municipalities treated as modeled annual loss probability zero

## What Was Fixed in the Draft

- Research questions were sharpened from broad basin/protection wording to corridor-based exposure and modeled protection/loss outcomes.
- The old study-area figure placeholder was replaced with an actual map and output path.
- Figure 3.1 now uses the absolute path to the conceptual model.
- Chapter 4 now reflects the final vulnerability-index logic rather than the old 51-variable PCA workflow.
- Old `FINALISE` markers for exposure-curve use, no-event interpretation, and statistical strategy were replaced with final working decisions.
- Chapter 5 now has a full results-figure skeleton and each figure is introduced in the surrounding text.

## Important Remaining Writing Work

Chapter 5 is now structurally ready but still needs to be written as full prose. The figure order is set, but the text still needs:

- actual result interpretation for each figure;
- careful wording around the weak overall RP100 exposure-vulnerability correlation;
- stronger interpretation of the protection/loss portfolio pattern by vulnerability quintile;
- a clear mini-synthesis before the Discussion.

Chapter 6 is still missing as a real discussion chapter. It should probably be structured around:

- why the exposure-vulnerability relationship is not simply linear;
- what the bivariate map shows that the correlation hides;
- why protection/loss outcomes add a second justice-relevant layer;
- how to interpret no-event municipalities without overstating high protection;
- model limitations, especially smaller streams not represented in the provider simulation;
- scale limits of municipality-level analysis;
- implications for flood justice and German flood risk governance.

Chapter 7 is still missing as full prose.

## Figure Roadmap

### Chapter 2

Figure 2.1 Study area map:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/maps/map_study_area_corridor.png`

### Chapter 3

Figure 3.1 Conceptual model:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/THESIS_DRAFTS/FIGURES/FIGURE_3_1_CONCEPTUAL_MODEL.png`

### Chapter 5 Main Results

Figure 5.1 Study area map, optional repeat only if not used as Figure 2.1:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/maps/map_study_area_corridor.png`

Figure 5.2 RP100 exposure:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/maps/map_rp100_exposure_corridor.png`

Figure 5.3 Vulnerability index:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/maps/map_vulnerability_index_corridor.png`

Figure 5.4 Vulnerability dimensions:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/map_vulnerability_dimensions.png`

Figure 5.5 Vulnerability vs RP100 exposure:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/plot_vulnerability_vs_rp100_exposure.png`

Figure 5.6 RP100 exposure by vulnerability quintile:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/plot_rp100_exposure_by_vulnerability_quintile.png`

Figure 5.7 Bivariate exposure-vulnerability map:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/map_exposure_vulnerability_bivariate.png`

Figure 5.8 Exposure-curve typology:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/EXPOSURE_CURVES/outputs/plots/map_exposure_curve_types.png`

Figure 5.9 Exposure curves by vulnerability quintile:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/plot_exposure_curves_by_vulnerability_quintile.png`

Figure 5.10 Protection coverage:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/maps/map_corridor_protection_coverage.png`

Figure 5.11 Protection return periods:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/maps/map_corridor_protection_return_period.png`

Figure 5.12 Modeled annual loss probability:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/map_modeled_annual_loss_probability.png`

Figure 5.13 Provider loss status by vulnerability quintile:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/plot_protection_status_by_vulnerability_quintile.png`

Figure 5.14 Annual loss probability by vulnerability quintile:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/plot_annual_loss_probability_by_vulnerability_quintile.png`

Figure 5.15 Finite protection return periods by vulnerability quintile:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/results_figures/plot_protection_return_period_by_vulnerability_quintile.png`

## Candidate Appendix or Methods Diagnostics

These are useful but probably too detailed for the main narrative unless space allows:

- PCA scree plot: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/plots/scree_kaiser_corridor.png`
- cumulative variance plot: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/plots/cumulative_variance_corridor.png`
- PCA component maps: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/outputs/maps/map_PC1_corridor.png` through `map_PC4_corridor.png`
- no-event stream context map: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/maps/map_no_loss_stream_context.png`

## Most Important Next Step

Write Chapter 5 as full result prose using the figure order already inserted in the draft. This should come before polishing Chapter 1-4 again, because the final wording of the Discussion depends on what the Results chapter actually says.

