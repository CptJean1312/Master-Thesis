# Chapter 4 Skeleton: Data & Methods

Working chapter skeleton based on the current method protocol, the method audit, and the available literature.

Purpose:
Chapter 4 must explain how the conceptual framework of Chapter 3 is operationalised empirically. It should be technical enough to be reproducible, but it should not read like a complete processing diary. The chapter should focus on the active final workflow and mark unresolved elements explicitly.

Core rule:

- Chapter 2 explains the concepts.
- Chapter 3 explains the research design.
- Chapter 4 explains the data and methods used to operationalise that design.

---

# 4.1 Data Sources

## What this section must do

Introduce all datasets that are actually used in the active workflow and explain their role in the analysis.

## What should go in

- VG250 municipality boundaries as spatial unit of analysis
- Elbe basin polygon as hydrological context and processing mask
- EFAS / EU multi-return-period flood rasters as the main flood hazard source
- INKAR as the main socio-economic dataset
- GISD only if used later as validation or contextualisation

## What should stay out or be only briefly contextualised

- BfG `Überflutungstiefen-DE` as active hazard source
- PESeta protection standard as active protection layer
- BfG `ManMadeObject-DE` unless it becomes part of final protection analysis
- DEM / DOM / DGM unless actually used analytically

## Key details to include

- VG250 source: BKG
- INKAR source: BBSR, INKAR 07/2025
- INKAR territorial status: 31.12.2023
- EFAS source: exact official citation still needs final confirmation
- all spatial data harmonised to `EPSG:25832`

## Literature support

Must-have:

- de Moel et al. (2015) for the importance of spatial scale and assessment purpose in flood-risk assessment
- Sairam et al. (2021) for process-based flood risk assessment in Germany and spatial consistency
- Alfieri et al. / Rojas et al. for European flood hazard / return-period logic

Optional:

- Rojas et al. (2013) for pan-European return-period flood hazard and adaptation context
- Fekete (2009, 2010) for German vulnerability data and subnational vulnerability assessment

## Placeholder

`[INSERT FINAL EFAS DATASET CITATION AND ACCESS DETAILS]`

---

# 4.2 Spatial Preprocessing and Study Area Definition

## What this section must do

Explain how the spatial units and the empirical study population were prepared.

## What should go in

- all layers transformed or prepared in `EPSG:25832`
- VG250 municipality geometries filtered to `Geofaktor_GF = "mit Struktur Land"`
- Elbe basin polygon used as a processing mask
- EFAS rasters clipped to the Elbe basin
- RP500 flood extent used to define the empirical corridor
- municipalities retained if they have non-zero RP500 flooded area

## Key distinction

The Elbe basin is the broader hydrological and governance context. The empirical population is the RP500 flood corridor.

## Recommended paragraph structure

Paragraph 1:
Describe CRS harmonisation and why area consistency matters.

Paragraph 2:
Describe VG250 municipality preparation and the land-only geometry issue.

Paragraph 3:
Describe the transition from basin mask to RP500 corridor.

## Strong sentence to include

`The analytical sample consists of municipalities intersecting the RP500 flood extent, while the Elbe basin remains the broader contextual and processing frame.`

## Literature support

Must-have:

- de Moel et al. (2015): flood risk assessments vary by spatial scale and assessment use
- Poussard et al. (2021): inequalities in flood exposure depend strongly on scale
- Fekete (2009, 2010): vulnerability and exposure must be aligned to meaningful spatial units

---

# 4.3 Flood Hazard and Exposure Construction

## What this section must do

Explain how flood exposure is derived from the prepared flood rasters and how it becomes a municipality-level variable.

## What should go in

- EFAS flood rasters for `RP10`, `RP20`, `RP50`, `RP100`, `RP200`, `RP500`
- rasters contain positive flood-depth values in flooded pixels
- `NoData` represents non-flooded areas
- flood extent is defined by valid depth pixels
- no polygonisation required
- raster-based zonal extraction used for municipality-level flooded area
- flood share calculated as flooded area divided by municipality area

## Formula

```text
flood_share_rp = flood_area_rp_m2 / municipality_area_m2
```

## Return-period logic

- `RP500` defines the study corridor
- `RP100` serves as the current main exposure benchmark
- all return periods are retained for exposure curves and later protection-related interpretation

## Methodological choice

The current main exposure metric uses flood extent, not depth-weighted exposure.

Reason:

- stable across all return periods
- comparable across municipalities
- interpretable as municipality area share
- suitable for constructing exposure curves

Limitation:

- does not distinguish shallow from deep inundation
- depth-weighted metrics remain a possible sensitivity extension

## Literature support

Must-have:

- Rojas et al. (2013): return-period flood hazard in European flood-risk assessment
- Alfieri et al. / Dottori et al. style sources if exact EFAS/JRC flood hazard map source is confirmed
- Sairam et al. (2021): large-scale process-based flood-risk assessment for Germany
- de Moel et al. (2015): spatial scale and uncertainty in flood-risk assessments

Optional:

- Poussard et al. (2021): municipality exposure classification and scale sensitivity
- Qiang (2019): overlay of flood hazard and population/social characteristics

## Placeholder

`[INSERT FINAL WORDING ON EFAS DATA SOURCE AND WHETHER THE RASTERS INCLUDE/EXCLUDE PROTECTION ASSUMPTIONS]`

---

# 4.4 Exposure Quality Checks

## What this section must do

Show that the exposure dataset was checked systematically.

## What should go in

- share bounds check: `0 <= flood_share <= 1`
- area bounds check: `flood_area <= municipality_area`
- monotonicity check: `RP10 <= RP20 <= RP50 <= RP100 <= RP200 <= RP500`
- suspicious municipalities exported for inspection
- non-monotonic cases treated carefully in exposure curve metrics

## Why it matters

The exposure pipeline combines raster data and municipal geometries. Automated checks protect against invalid area denominators, raster alignment issues, and implausible return-period behaviour.

## Literature support

Must-have:

- de Moel et al. (2015) for uncertainty and validation problems across flood-risk assessments
- Rojas et al. (2013) or Alfieri et al. for uncertainty in large return periods

Optional:

- Sairam et al. (2021) for consistency in German large-scale flood-risk modelling

---

# 4.5 Socio-Economic Data Processing

## What this section must do

Explain how INKAR was transformed from raw long format into the municipality-level analysis table.

## What should go in

- INKAR 2025 / edition 07/2025
- municipality-level data
- long format imported with `data.table::fread()`
- filtered to relevant municipality records
- latest available observation per municipality and indicator retained
- pivoted to wide format
- AGS used as join key
- joined to corridor municipalities
- one missing municipality in latest-build raw INKAR table documented

## Key numbers

- full INKAR raw file approximately 6.7 GB
- corridor municipalities: `835`
- latest raw-INKAR records matched for `834`
- available original INKAR indicators in corridor: `176`
- missing municipality: `16076094`, Berga-Wünschendorf

## Wording caution

Do not make the AGS-prefix filtering sound like the final study-area definition. The final study area is defined spatially by RP500 corridor intersection. AGS filtering is a data-processing step.

## Literature support

Must-have:

- Fekete (2009, 2010): social vulnerability for river floods in Germany and subnational indicator construction
- Rufat et al. (2015): social vulnerability indicator selection is context-specific and multidimensional

Optional:

- GISD / Michalski et al. only if used as validation context

---

# 4.6 Indicator Screening and Candidate Vulnerability Variables

## What this section must do

Explain how the socio-economic indicator pool was reviewed before PCA.

## What should go in

- 176 available INKAR indicators reviewed
- grouped into substantive dimensions
- coverage rates calculated
- missingness assessed
- correlation matrix used diagnostically
- redundant, nested, low-coverage, and ambiguous variables flagged
- old 51-variable set retained as exploratory benchmark
- current 23-variable thesis candidate set documented but not yet final locked
- student variables treated as sensitivity block

## Why this section matters

It responds directly to supervisor concerns about counterintuitive definitions and redundant variables.

## Specific methodological points

- `a_ALGII_SGBII` and `a_Unterkunft_SGBII` should not both be treated as independent core deprivation indicators in the final main index.
- `a_BG1P` should not be interpreted as single-parent households.
- student indicators may reflect precariousness but also university-town structure.
- digital infrastructure should avoid using several bandwidth thresholds simultaneously.

## Literature support

Must-have:

- Rufat et al. (2015): indicator selection, scale, weighting and aggregation are central issues in social vulnerability modelling
- Tate (2013): uncertainty in vulnerability indices arises from methodological choices
- Schmidtlein et al. (2008): sensitivity of SoVI to variable selection and modelling choices
- Rufat et al. (2019): social vulnerability models are useful but validity should not be assumed

Optional:

- Fekete (2009): variable selection and factor analysis for German river-flood vulnerability

## Placeholder

`[INSERT FINAL PCA VARIABLE LIST ONCE DECIDED]`

---

# 4.7 PCA-Based Vulnerability Index Construction

## What this section must do

Explain the final construction of the socio-economic vulnerability index.

## Current status

This section cannot be fully final yet because the final PCA input set is still open.

## What should eventually go in

- final input variable list
- missing-value handling
- z-standardisation
- PCA with `prcomp()`
- retained component rule
- explained variance
- interpretation of leading PCs
- aggregation into composite index
- direction anchoring
- sensitivity analysis

## Current working version

The active corridor analysis script currently uses the broad `original_51` PCA input set.

The working preferred candidate for the thesis is `thesis_candidate_23`, but this has not yet replaced the active script.

## Index formula

The current broad working index uses a variance-weighted sum of retained PC scores:

```text
V_i = sum(PC_ik * w_k)
```

where:

- `V_i` is the vulnerability index for municipality `i`
- `PC_ik` is municipality `i`'s score on retained component `k`
- `w_k` is the explained-variance weight of component `k`

The resulting index is then standardised as a z-score.

## Wording caution

Do not present the current 51-variable index as final if the thesis will switch to `thesis_candidate_23`.

## Literature support

Must-have:

- Fekete (2009): PCA / factor analysis for river-flood social vulnerability in Germany
- Tate (2013): uncertainty in social vulnerability indices
- Schmidtlein et al. (2008): sensitivity analysis of SoVI
- Rufat et al. (2015): methodological cautions in social vulnerability measurement
- Rufat et al. (2019): validation of social vulnerability models

## Placeholder

`[FINALISE: main PCA input set, component retention rule, direction anchor, robustness specification]`

---

# 4.8 Multi-Return-Period Exposure Curves

## What this section must do

Explain how multi-RP flood shares are used to describe exposure dynamics and prepare protection-related interpretation.

## What should go in

- each municipality has a flood-share curve across `RP10`, `RP20`, `RP50`, `RP100`, `RP200`, `RP500`
- metrics derived from curves:
  - onset return period
  - slope over log return period
  - normalized AUC
  - realized share by RP100
  - late-growth share after RP100
  - maximum jump interval
- non-monotonic cases adjusted for metric calculation using monotone `cummax`
- first typology:
  - early exposure
  - gradual increase
  - delayed jump

## Conceptual role

Exposure curves do not directly prove protection. They describe how exposure accumulates across return periods and can later support inferred protection-related interpretation when additional simulation or damage-based data are available.

## Literature support

Must-have:

- Serra-Llobet et al. (2022): residual risk behind levees and protected areas
- Fu et al. (2023): residual flood risk and protective structures
- Tobin (1995): levee effect and false security
- Barendrecht et al. (2018): socio-hydrological model for the Elbe

Optional:

- Ferdous et al. (2020): structural protection, exposure accumulation, and levee effect

## Placeholder

`[INSERT FINAL PROTECTION-RELATED METRIC IF / WHEN NEW PROTECTION DATA ARRIVE]`

---

# 4.9 Statistical and Spatial Analysis

## What this section must do

Describe the quantitative comparison between flood exposure, vulnerability, and protection-related metrics.

## Current status

This part is not final because the final vulnerability index and protection variable are still pending.

## What should eventually go in

- descriptive statistics
- maps
- bivariate correlations
- regression model for RP100 exposure and vulnerability
- possible spatial diagnostics
- models for exposure-curve metrics
- protection-related comparison once data are available

## Recommended current placeholder

The first empirical comparison focuses on the relationship between the corridor-specific vulnerability index and `RP100` flood share. Multi-return-period exposure metrics are retained for later analysis of exposure dynamics and protection-related differentiation.

## Literature support

Must-have:

- Poussard et al. (2021): comparing exposure and deprivation across spatial scales
- Qiang (2019): disparity analysis of flood exposure
- Odersky and Löffler (2024): German differential flood exposure

Placeholder:

`[FINALISE MODEL SPECIFICATION AFTER INDEX DECISION]`

---

# 4.10 Software and Reproducibility

## What this section must do

Document the computational environment and reproducibility logic.

## What should go in

- R as main analysis environment
- `sf` for vector data
- `terra` for raster and zonal extraction
- `data.table` for large INKAR import
- `dplyr` / `tidyr` for data wrangling
- `ggplot2` / `geom_sf()` for mapping
- QGIS used for inspection and some base-layer preparation
- scripted outputs with organised folders

## Recommended wording

The workflow separates geometry-heavy preprocessing, raster-based exposure extraction, socio-economic data processing, and statistical analysis in a reproducible R-based workflow.

---

# Chapter 4 Writing Order

Recommended writing sequence:

1. Write `4.1 Data Sources`.
2. Write `4.2 Spatial Preprocessing and Study Area Definition`.
3. Write `4.3 Flood Hazard and Exposure Construction`.
4. Write `4.4 Exposure Quality Checks`.
5. Write `4.5 Socio-Economic Data Processing`.
6. Write `4.6 Indicator Screening and Candidate Vulnerability Variables`.
7. Leave `4.7 PCA-Based Vulnerability Index Construction` partly placeholder until final index choice.
8. Write `4.8 Multi-Return-Period Exposure Curves`.
9. Leave `4.9 Statistical and Spatial Analysis` partly placeholder until final models.
10. Write `4.10 Software and Reproducibility`.

This keeps the chapter moving while protecting the unresolved PCA and protection decisions.
