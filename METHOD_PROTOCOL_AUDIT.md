# Method Protocol Audit

Purpose:
This audit checks which parts of `METHOD_PROTOCOL_CURRENT.md` are suitable for the final Chapter 4 methods section, which parts need revision, and which parts belong only in the workflow history.

Status labels:

- `FINAL CORE`: should be used in Chapter 4
- `KEEP AS CONTEXT`: can be mentioned briefly, but not as a central method
- `REVISE`: method is relevant, but wording or scope must be corrected before Chapter 4
- `PLACEHOLDER`: method depends on a later decision or missing data
- `SUPERSEDED`: was part of earlier work, but should not be presented as active final method

---

# 1. High-Level Assessment

The current thesis workflow is coherent, but Chapter 4 must be disciplined about the distinction between:

- the final active exposure workflow
- exploratory or diagnostic socio-economic work
- optional validation datasets
- older flood/protection approaches that were abandoned
- planned protection inference that is not yet empirically implemented

The final methods chapter should therefore not read like a full diary of every explored path. It should explain the active final workflow and use the older steps only where they justify methodological choices.

---

# 2. Final Core Components

## 2.1 Coordinate Reference System and Spatial Harmonisation

Status: `FINAL CORE`

Use in Chapter 4:

- explain that all spatial layers were harmonised to `EPSG:25832`
- justify this because municipality-level flood shares require consistent area calculations
- keep it short

Do not overdo:

- no need for long QGIS project setup details in final prose

---

## 2.2 VG250 Municipality Boundaries and Land-Only Geometry

Status: `FINAL CORE`

Use in Chapter 4:

- explain that VG250 municipalities form the spatial unit of analysis
- explain the `Geofaktor_GF = "mit Struktur Land"` filter
- this is important because it affects area denominators and flood-share validity

This is a strong methods point and should definitely be included.

---

## 2.3 Elbe Basin as Processing Mask

Status: `REVISE`

Use in Chapter 4:

- describe the basin only as hydrological context and raster processing mask
- make clear that the empirical study population is not the whole basin

Required wording:

- `The Elbe basin is used as a processing and contextual frame, while the empirical sample is defined by the RP500 flood corridor.`

Avoid:

- language implying that all basin municipalities are included in the analysis

---

## 2.4 EFAS Multi-Return-Period Flood Rasters

Status: `FINAL CORE`

Use in Chapter 4:

- this is the current central hazard source
- explain return periods `RP10`, `RP20`, `RP50`, `RP100`, `RP200`, `RP500`
- explain positive depth values and `NoData`
- explain that valid depth pixels define flood extent

Important:

- add exact source citation once confirmed
- avoid presenting these rasters as official protected/unprotected flood maps

---

## 2.5 RP500 Corridor Definition

Status: `FINAL CORE`

Use in Chapter 4:

- this is one of the most important updates
- explain why RP500 is used to define the flood-relevant study corridor
- explain that municipalities with non-zero RP500 flooded area are retained

Recommended final wording:

- `The analytical sample consists of municipalities intersecting the RP500 flood extent.`

Avoid:

- saying that RP500 is the main hazard outcome
- the main benchmark is currently RP100

---

## 2.6 Municipality-Level Flood Exposure Across Return Periods

Status: `FINAL CORE`

Use in Chapter 4:

- explain flooded area and flood share per municipality
- explain all return periods are processed
- explain RP100 is current main benchmark
- explain multi-RP structure is retained for exposure curve / protection inference

Very important:

- include the formula:

```text
flood_share_rp = flooded_area_rp_m2 / municipality_area_m2
```

---

## 2.7 Exposure Quality Checks

Status: `FINAL CORE`

Use in Chapter 4:

- share bounds
- area bounds
- monotonicity
- suspicious municipality export

This strengthens the credibility of the exposure pipeline.

---

## 2.8 INKAR Socio-Economic Data Processing

Status: `FINAL CORE`

Use in Chapter 4:

- full INKAR long-format import
- municipality-level filtering
- latest available observation per municipality-indicator
- wide-format transformation
- AGS join to corridor municipalities

Need final care:

- the AGS-prefix filtering should be described carefully and not overemphasized if the final corridor join is the more important selection step

---

## 2.9 INKAR Indicator Inventory and Screening

Status: `FINAL CORE`

Use in Chapter 4:

- all 176 corridor indicators were reviewed
- indicators grouped into substantive domains
- correlation and coverage checks used to detect redundancy and problematic variables

Important:

- this is the bridge to the final index decision
- keep the final prose neutral until the final PCA set is locked

---

## 2.10 PCA-Based Vulnerability Index

Status: `REVISE / PLACEHOLDER`

Current situation:

- active corridor script still uses the original broad 51-variable PCA
- the 23-variable thesis candidate set is documented and likely more defensible
- final index choice is not locked

Use in Chapter 4 now:

- write as a placeholder section with current options
- do not claim the 23-variable set is already the final implemented index unless we update the script

Required next step:

- finalise PCA input set
- rerun final vulnerability pipeline
- update Chapter 4 accordingly

---

## 2.11 Multi-RP Exposure Curves

Status: `FINAL CORE / BRIDGE TO PROTECTION`

Use in Chapter 4:

- exposure curves are already computed
- metrics include onset, normalized AUC, realized-by-RP100, late-growth share, and jump intervals
- these metrics prepare the later protection inference

Important:

- do not overclaim that they already measure true protection
- phrase as exposure dynamics and protection-related interpretation

---

# 3. Components to Keep Only as Context or Optional

## 3.1 GISD

Status: `KEEP AS CONTEXT / OPTIONAL`

Use in Chapter 4:

- mention only if actually used as validation or robustness check
- otherwise keep out of the main methods chapter

Reason:

- the main vulnerability index is INKAR-based
- GISD was explicitly not meant to replace the INKAR index

---

## 3.2 BfG Überflutungstiefen-DE

Status: `SUPERSEDED`

Use in Chapter 4:

- not part of active final hazard workflow
- may be briefly mentioned only if explaining why the workflow changed

Reason:

- BfG vector-based flood zone logic was superseded by EFAS multi-return-period raster workflow
- Zone 3 protection logic was not robust enough

Avoid:

- detailed BfG geometry repair workflow in final methods chapter
- old zone 1 / zone 3 residual-risk framing

---

## 3.3 PESeta Protection Dataset

Status: `SUPERSEDED`

Use in Chapter 4:

- should not be part of the current final method

Reason:

- NUTS3 protection standard was too coarse
- direct binary masking of HQ100 exposure removed meaningful variation
- not suitable as direct municipality-level protected/unprotected layer

---

## 3.4 BfG ManMadeObject-DE

Status: `OPTIONAL / POSSIBLE FUTURE EXTENSION`

Use in Chapter 4:

- only include if later used for infrastructure distribution
- otherwise keep in data source inventory, not active methods

---

## 3.5 DEM / DGM / DOM Products

Status: `BACKGROUND / NOT ACTIVE`

Use in Chapter 4:

- only include if used in actual analysis
- current active workflow does not require DEM-derived slope, flow accumulation, or topographic covariates

Important correction:

- Copernicus DEM is a DSM, not a bare-earth DTM
- do not describe it as object-free ground surface

---

# 4. Critical Corrections Before Writing Chapter 4

## 4.1 Fix Section Numbering

Current protocol has:

- `114. Composite Vulnerability Index`

This should be:

- `14. Composite Vulnerability Index`

---

## 4.2 Fix Typos and Language

Known examples:

- `PYSICSAL` -> `PHYSICAL`
- `MUNCIPALITY` -> `MUNICIPALITY`
- `slood share` -> `flood share`
- `Wewhile` in thesis draft -> should be fixed in Chapter 2 / Chapter 3 transition

---

## 4.3 Clarify Current Analytical Population

The final methods chapter must consistently say:

- broader frame: Elbe basin / Elbe river basin
- empirical sample: municipalities intersecting RP500 flood extent

This distinction is central.

---

## 4.4 Clarify Current Exposure Outcome

The final methods chapter must consistently say:

- multi-RP exposure is calculated for all return periods
- RP100 is current main benchmark
- RP500 defines the corridor
- RP500 is not the main outcome

---

## 4.5 Clarify Current Protection Status

The final methods chapter must consistently say:

- protection is conceptually central
- direct protection data are not yet used in the active workflow
- exposure curves prepare later inferred protection analysis
- no final protection index exists yet

---

## 4.6 Clarify Vulnerability Index Status

The final methods chapter should currently treat the vulnerability index section as partly provisional.

Open decision:

- use original broad 51-variable PCA as final main index
- or switch to `thesis_candidate_23`
- or use one as main and the other as sensitivity

Recommended direction:

- use `thesis_candidate_23` as final main index if supervisors accept the curation logic
- keep `original_51` as exploratory / robustness benchmark
- keep `student_sensitivity_26` as sensitivity analysis

---

# 5. Recommended Chapter 4 Structure Based on This Audit

## 4.1 Data Sources

Include:

- VG250 municipality boundaries
- EFAS multi-return-period flood rasters
- Elbe basin polygon as processing mask
- INKAR
- optional GISD only if used

Keep out or mention only briefly:

- BfG flood zones
- PESeta
- ManMadeObject-DE
- DEMs unless used

---

## 4.2 Spatial Preprocessing and Study Area Definition

Include:

- CRS harmonisation
- land-only municipality geometries
- Elbe basin clipping mask
- RP500 corridor definition

---

## 4.3 Flood Exposure Construction

Include:

- valid depth pixel = flooded
- NoData = not flooded
- flooded area per return period
- flood share formula
- RP100 main benchmark
- all-RP exposure table
- quality checks

---

## 4.4 Socio-Economic Data Processing

Include:

- INKAR import
- latest-year selection
- wide transformation
- AGS join
- coverage and missingness

---

## 4.5 Indicator Screening and Vulnerability Index Construction

Include:

- full 176 indicator inventory
- correlation screening
- PCA candidate sets
- final PCA set once decided
- median imputation
- z-standardisation
- PCA
- index aggregation
- direction anchoring
- sensitivity checks

Placeholder:

- final PCA input set and final retained PC rule

---

## 4.6 Multi-Return-Period Exposure Curves and Protection-Related Interpretation

Include:

- onset metrics
- slope / AUC / late-growth share
- typology if used
- careful wording that these are exposure dynamics, not direct proof of protection

Placeholder:

- final protection index / inferred protection variable

---

## 4.7 Statistical and Spatial Analysis

Include later:

- descriptive maps and summary statistics
- bivariate comparisons
- regression or spatial models if used
- robustness checks

Placeholder:

- final model specification

---

## 4.8 Software and Reproducibility

Include:

- R
- `sf`
- `terra`
- `data.table`
- `dplyr`
- `ggplot2`
- `spdep` / spatial packages if used
- QGIS only for base inspection / preparation where relevant

---

# 6. Immediate Next Steps

1. Confirm exact EFAS data source citation and metadata.
2. Decide whether Chapter 4 should mention GISD at all.
3. Decide whether the active final PCA should be switched from `original_51` to `thesis_candidate_23`.
4. Write Chapter 4 skeleton from this audit.
5. Scrape literature for Chapter 4 method support.
6. Draft Chapter 4 with placeholders for the final PCA and protection index.
