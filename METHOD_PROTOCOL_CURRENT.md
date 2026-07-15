# Method Protocol Current

Working protocol for the current thesis workflow.

Project: **Elbe Flood Exposure, Vulnerability & Inferred Protection Analysis**

Status note:
This document records the current processing logic as a working protocol. It is not yet the final Chapter 4 text. Some steps were completed during earlier workflow stages and need to be audited before they are used in the final methods chapter.

Update 2026-07-15:
The vulnerability-index workflow has been updated after the INKAR/PCA review. The active specification is now documented in:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/VULNERABILITY_INDEX_FINAL_SPEC.md`

In short, the old broad `original_51` PCA is no longer the active main vulnerability index. The main index now uses the screened `thesis_candidate_23` indicator set and is constructed as a direction-coded, dimension-balanced composite across three dimensions: demographic and household structure, deprivation and labour-market position, and access and adaptive capacity. PCA is retained as a diagnostic and robustness tool. Older notes in Sections 13 and 14 that describe the `original_51` PCA as the active working index should be read as historical documentation of the previous implementation.

---

# A. GIS / Spatial Data Preparation

Base layers, study area definition, and hazard raster preparation.

---

## 1. Project Initialisation

Status: `DONE`

A dedicated GIS project folder structure was created to support transparency and reproducibility:

- `raw/`
- `processed/`
- `analysis/`
- `tables/`
- `maps/`

A unified project CRS was defined as:

- `ETRS89 / UTM Zone 32N`
- `EPSG:25832`

This CRS was fixed from the outset to ensure consistency in all area-based calculations and municipality-level overlay operations.

---

## 2. Import of Base Layers

Status: `DONE`

Core vector datasets were loaded and harmonised:

- Elbe main river polyline
- federal state boundaries
- Elbe basin polygon
- official administrative boundaries at VG250 municipality level

All relevant vector layers were stored in GeoPackage (`.gpkg`) format to standardise geometry handling and avoid format-specific limitations.

---

## 3. Elbe Basin Definition

Status: `DONE`

The Elbe River Basin District polygon was loaded and exported as a separate basin layer.

The basin polygon was reprojected to `EPSG:25832`.

In the current workflow, the Elbe basin serves primarily as a processing mask for clipping the EU flood rasters and reducing data volume.

---

## 4. Administrative Boundary Preparation

Status: `DONE`

VG250 municipality polygons were loaded from the official BKG dataset.

A critical issue was identified: some municipalities were duplicated by AGS due to the attribute `Geofaktor_GF`, distinguishing:

- `mit Struktur Land`
- `mit Struktur Gewässer`

Using the `Gewässer` geometries led to invalid municipality areas and inflated flood shares.

Therefore, the municipality layer was filtered to:

- `Geofaktor_GF = "mit Struktur Land"`

This land-only municipality geometry serves as the consistent spatial base for all subsequent exposure calculations.

---

## 5. EU Flood Hazard Raster Preparation

Status: `DONE`

Flood hazard maps are no longer derived from the earlier BfG vector-based flood zone workflow. Instead, the current workflow uses EU / EFAS flood hazard rasters for multiple return periods.

### 5.1 Input Data

Prepared raster layers:

- `RP10`
- `RP20`
- `RP50`
- `RP100`
- `RP200`
- `RP500`

These rasters contain:

- positive flood depth values in flooded cells
- `NoData` values outside flooded areas

### 5.2 Reprojection and Clipping

All EFAS flood rasters were preprocessed in R:

1. Reprojected from `EPSG:3035` (`ETRS89-extended / LAEA Europe`) to `EPSG:25832` (`ETRS89 / UTM Zone 32N`)
2. Clipped to the Elbe basin polygon

This preprocessing step:

- reduces data volume
- ensures spatial consistency with municipality geometries
- avoids repeated raster handling in QGIS

### 5.3 Conceptual Interpretation

The rasters are interpreted as flood-depth rasters, not binary flood masks.

Important implication:

- flooded areas are identified by the presence of valid flood-depth pixels
- non-flooded areas are represented by `NoData`
- no separate threshold-based binary conversion is required for the core extent workflow

---

## 6. Study Area Definition

Status: `UPDATED / DONE`

The study area is now defined as the Elbe flood corridor, operationalised as:

> all municipalities with non-zero flood exposure in the RP500 flood raster

This replaces the earlier logic of using the hydrological Elbe basin as the direct analytical population.

The new definition:

- focuses the analysis on municipalities that are actually flood-relevant
- avoids including basin municipalities with no meaningful flood exposure
- provides a hazard-based empirical study area

Implementation:

- the `RP500` raster is used to identify municipalities intersecting valid flood-depth cells
- municipalities with `rp500_flood_area_m2_exact > 0` are retained in the study population

Output:

- `municipalities_corridor.gpkg`

---

# B. R Workflow

Hazard exposure construction, socio-economic processing, vulnerability analysis.

---

## 1. Import of Large INKAR Dataset

Status: `DONE`

The full INKAR 2025 dataset, approximately 6.7 GB in long format, was imported using:

- `data.table::fread()`

Dataset structure:

- indicators
- spatial units
- time
- values

---

## 2. Filtering to Relevant Spatial Units

Status: `DONE / NEEDS FINAL WORDING`

Retained:

- `Bereich == "LRB"`
- `Raumbezug == "Gemeinden"`

Municipalities were filtered by AGS prefixes covering the Elbe-related federal states and edge municipalities:

- `01`, `02`, `03`, `09`, `12`, `13`, `14`, `15`, `16`

This step ensured compatibility with the municipality layer used in the spatial workflow.

Audit note:
The final methods chapter should clarify whether this AGS-prefix filtering is still needed after the corridor definition, or whether the final logic should be described primarily as a join between the corridor municipalities and the full municipality-level INKAR extract.

---

## 3. Selection of Most Recent Year per Indicator

Status: `DONE`

INKAR records were ordered by:

- municipality
- indicator
- descending year

Only the most recent available observation per municipality and indicator was retained.

This maximises data completeness while preserving temporal traceability.

---

## 4. Pivot to Wide Format

Status: `DONE`

Long-format INKAR data were transformed to wide format using:

- `dcast()`

Result:

- one row per municipality
- one column per indicator

A parallel wide table of reference years was created and merged back into the main data.

---

## 5. Indicator Curation and Harmonisation

Status: `UPDATED / DONE, NOT FINAL LOCKED`

All INKAR indicators available for corridor municipalities were systematically reviewed after restricting the study area to the RP500 corridor.

The resulting corridor-level indicator inventory contains:

- `176` variables
- generally high coverage

Indicators were grouped into substantive domains:

- basic population and structural counts
- demography and household structure
- economic structure and tourism
- education
- health, services and accessibility
- housing and land-use context
- income, deprivation and social transfers
- labour market and employment

A structured screening system was applied to classify variables as:

- `core_curated_pca`
- `context_or_control_only`
- `review_not_in_curated_set`
- `drop_*` categories, for example low coverage, redundancy, compositional overlap, nested age block

This step explicitly addressed:

- duplicated constructs
- compositional or complementary indicators
- nested demographic structures
- specialist-specific duplicates
- low-coverage variables

Variable names were harmonised into readable analysis-ready names.

`Kennziffer` was renamed to `AGS`.

Audit note:
The final PCA input set is not yet locked. The current documentation distinguishes the original 51-variable exploratory PCA, the full 176-variable diagnostic pool, and a 23-variable working thesis candidate set.

---

## 6. Integration of Socio-Economic Data with Municipalities

Status: `DONE`

The cleaned INKAR dataset was linked to the municipality layer through the `AGS` key.

Remaining join inconsistencies due to state filtering were corrected by rerunning INKAR extraction including Bavaria and Thuringia prefixes where needed.

---

## 7. Municipality Identifier Consolidation

Status: `DONE`

Where multiple AGS-like fields existed due to different joins, a robust unified key was created:

- `AGS_final = coalesce(Gemeindeschlüssel_AGS, AGS, Gemeindeschlüssel_AGS_2)`

The identifier was:

- padded to 8 digits
- retained as character

This key is used consistently throughout the analytical workflow.

---

## 8. Exposure Construction from EFAS Flood Rasters

Status: `NEW / DONE`

This is the central updated hazard workflow.

### 8.1 Core Principle

Flood exposure is derived directly from valid flood-depth pixels in the EFAS rasters.

Interpretation:

- valid depth pixel = flooded
- `NoData` = non-flooded

### 8.2 Corridor Derivation from RP500

A municipality-level flood corridor is constructed using `RP500`:

1. Municipal polygons are raster-aligned
2. A first raster-based zonal preselection identifies potentially intersecting municipalities
3. Exact flooded area is then extracted for preselected municipalities
4. Municipalities with non-zero `RP500` flooded area are retained

Output:

- `municipalities_corridor.gpkg`

### 8.3 Flooded Area Calculation per Return Period

For each return period raster:

- `RP10`
- `RP20`
- `RP50`
- `RP100`
- `RP200`
- `RP500`

the pipeline calculates:

- flooded area per municipality, e.g. `flood_area_rpXX_m2`
- flooded share per municipality, e.g. `flood_share_rpXX`

Method:

- municipality layer is converted to terra vector format
- cell area is calculated from raster resolution
- flooded area is extracted as the sum of valid flooded cell areas intersecting each municipality
- no raster polygonization is required

### 8.4 Final Exposure Dataset

All return periods are merged into one final municipality-level dataset containing:

- `AGS`
- `mun_name`
- `municipality_area_m2`
- `flood_area_rp10_m2`
- `flood_area_rp20_m2`
- `flood_area_rp50_m2`
- `flood_area_rp100_m2`
- `flood_area_rp200_m2`
- `flood_area_rp500_m2`
- `flood_share_rp10`
- `flood_share_rp20`
- `flood_share_rp50`
- `flood_share_rp100`
- `flood_share_rp200`
- `flood_share_rp500`

Outputs:

- `municipality_flood_exposure_all_RPs.csv`
- `municipalities_corridor_exposure_all_RPs.gpkg`

---

## 9. Exposure Quality Checks

Status: `NEW / DONE`

Automated quality checks are performed to validate the exposure pipeline.

### 9.1 Share Bounds Check

For each return period, flood shares must satisfy:

- `0 <= flood_share <= 1`

### 9.2 Area Bounds Check

For each return period:

- flooded area must not exceed municipality area

### 9.3 Monotonicity Check

For each municipality, the following logical sequence is checked:

- `RP10 <= RP20 <= RP50 <= RP100 <= RP200 <= RP500`

This checks that flood extent does not decrease with increasing return period.

### 9.4 Suspicious Municipalities

Municipalities violating any of the checks are flagged and exported separately for inspection.

Outputs:

- `exposure_quality_checks.csv`
- `exposure_quality_checks_suspicious_only.csv`
- `quality_check_summary.csv`

---

## 10. Methodological Note on Flood Extent vs. Flood Depth

Status: `NEW / CURRENT MAIN CHOICE`

At the current stage, the main exposure metric is based on flood extent, operationalised as the share of municipality area covered by valid flooded raster pixels.

This is a deliberate simplification.

Reason:

The primary analytical goal is to build a stable, comparable, reproducible municipality-level exposure curve across multiple return periods.

Area-based flood exposure is currently preferred because it:

- is robust across all return periods
- is directly interpretable
- forms the basis for later inference of effective protection levels

Limitation:

This approach does not yet use flood depth directly in the main metric. Thus, differences such as `0.1 m` versus `3 m` flood depth are not currently represented in the core exposure indicator.

Planned extension or sensitivity analysis:

- mean flood depth per municipality and return period
- maximum flood depth per municipality and return period
- depth-weighted exposure metrics

For the current phase, the priority is to complete the multi-return-period flood share dataset.

---

## 11. Integration of the German Index of Socioeconomic Deprivation

Status: `DONE / OPTIONAL ANALYTICAL EXTENSION`

The GISD dataset was prepared as an external deprivation reference.

GISD is not used to replace the INKAR-based vulnerability index. Instead, it may serve as:

- contextualisation tool
- external validation source
- robustness check for deprivation patterns

Audit note:
GISD should not appear as part of the main analysis unless it is actually used later. In Chapter 4 it should probably be framed as optional or omitted until needed.

---

## 12. Exploratory Data Analysis and Correlation Screening

Status: `UPDATED / DONE`

Before multivariate modelling, the full corridor-level socio-economic indicator pool was explored systematically.

### 12.1 Corridor-Wide Indicator Inventory

After restricting the sample to municipalities in the RP500 corridor:

- `176` INKAR indicators were available in the latest-build table
- coverage and thematic grouping were reviewed for all available corridor indicators
- one municipality, `16076094` Berga-Wünschendorf, remained missing in the current raw-INKAR latest-build table

### 12.2 Descriptive Statistics

Summary statistics were calculated for all retained variables:

- minimum
- median
- mean
- standard deviation

Missingness patterns were assessed and exported.

Coverage rates were used as one criterion for excluding weak indicators.

### 12.3 Correlation Analysis

A full 176-variable correlation matrix was produced as a diagnostic step.

The purpose was to detect:

- multicollinearity
- nested indicator structures
- compositional relationships
- duplicated constructs
- dominant structural dimensions

This full matrix was treated as a diagnostic tool, not as a direct model input.

### 12.4 Corridor-Specific Sample Logic

The socio-economic diagnostic stage is no longer based on all Elbe-basin municipalities, but only on municipalities with non-zero RP500 flood exposure.

This ensures that both vulnerability exploration and later PCA are estimated on the hazard-relevant population only.

---

## 13. Principal Component Analysis for Vulnerability Structure

Status: `UPDATED / DONE, FINAL CHOICE PENDING`

A wide exploratory PCA framework was implemented on the RP500 corridor municipalities only.

### 13.1 Available PCA Sets

Several PCA input sets were compared:

- `all_176`: all available corridor indicators after latest-build corridor filtering
- `original_51`: earlier broad exploratory PCA block from the original clean analysis
- `thesis_candidate_23`: current working thesis candidate set, reduced from the full corridor inventory
- `student_sensitivity_26`: thesis candidate set extended by student-related indicators as a sensitivity block

### 13.2 Role of the Original 51-Variable Set

The original 51-variable set was retained as an exploratory benchmark and historical reference.

It covered a broad vulnerability field including:

- poverty and welfare dependence
- labour market strain
- demographic composition
- household structure
- income groups
- healthcare access
- digital infrastructure
- accessibility
- density

However, this set was later identified as methodologically imperfect because it still included:

- duplicated constructs
- nested age blocks
- overlapping digital-access measures
- overlapping physician measures
- at least one definition-sensitive SGB-II pair

### 13.3 Current Working Thesis Candidate Set

The current working thesis candidate set contains `23` variables and represents the main thesis-oriented reduced PCA set at this stage.

Selection logic:

- retain broad but interpretable coverage of deprivation, labour-market strain, demographic sensitivity, household/social structure, accessibility, health access, and digital infrastructure
- remove direct duplicates or near-duplicates wherever possible
- remove variables with clearly problematic definitions for interpretation
- remain broad enough for PCA while avoiding obvious over-weighting of the same latent construct

### 13.4 Student Sensitivity Block

Student variables are treated separately as a sensitivity extension because:

- they may capture precarious socio-economic situations in some municipalities
- at municipality scale, they also strongly reflect university-centre structure

Therefore, they are not part of the current main thesis candidate set.

### 13.5 Data Handling

For the corridor-specific PCA workflow:

- all PCA variables were converted to numeric
- missing values were handled using median imputation
- indicators were standardised using z-scores
- PCA was performed using `prcomp()`

### 13.6 Corridor-Specific PCA Outputs

The current corridor PCA workflow produces:

- component scores
- scree plot
- cumulative variance plot
- top-loading tables
- missingness table for PCA inputs
- corridor-specific PCA variable list

### 13.7 PCA Comparison Across Candidate Sets

The following comparison was used to assess dimensional structure:

- `original_51`: 51 variables, 834 municipalities, PC1 variance 16.18%, cumulative PC1-PC4 43.40%, cumulative PC1-PC8 60.07%
- `all_176`: 176 variables, 834 municipalities, PC1 variance 15.14%, cumulative PC1-PC4 33.01%, cumulative PC1-PC8 44.14%
- `thesis_candidate_23`: 23 variables, 834 municipalities, PC1 variance 24.21%, cumulative PC1-PC4 54.81%, cumulative PC1-PC8 72.60%
- `student_sensitivity_26`: 26 variables, 834 municipalities, PC1 variance 22.14%, cumulative PC1-PC4 54.59%, cumulative PC1-PC8 72.14%

Interpretation:

- `all_176` is too broad to serve as the main thesis PCA input, but useful as a structural diagnostic
- `original_51` remains useful as an exploratory benchmark
- `thesis_candidate_23` is currently the most balanced working set for the thesis
- `student_sensitivity_26` is used to evaluate the impact of adding student-concentration variables

### 13.8 Active Vulnerability Input in the Corridor Analysis Script

The active corridor analysis script has been updated to use the screened `thesis_candidate_23` set.

The broad 51-variable PCA is retained only as historical exploratory context and is no longer the active main vulnerability implementation.

### 13.9 Interpretation

The first principal components were interpreted substantively using their dominant loadings.

The final main vulnerability index is not taken directly from these components. PCA is retained as a diagnostic and robustness tool because the corridor PCA strongly reflects the dominant spatial covariance structure, especially accessibility and settlement-periphery patterns.

---

## 14. Composite Vulnerability Index

Status: `UPDATED / DONE, SEMI-FINAL`

A corridor-specific composite vulnerability index was constructed from the screened 23-variable INKAR set.

### 14.1 Current Implementation

The active analysis script computes the main vulnerability index as a direction-coded, dimension-balanced composite.

The specification is documented in:

`/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/VULNERABILITY_INDEX_FINAL_SPEC.md`

### 14.2 Main Index

The current working index is based on:

- the screened `thesis_candidate_23` indicator set
- z-standardisation of all indicators
- reverse-coding of capacity indicators
- three equally weighted dimensions:
  - demographic and household structure
  - deprivation and labour-market position
  - access and adaptive capacity

This produces:

- `vuln_index_main`
- `vuln_index_main_z`

### 14.3 Sensitivity Version

PCA sensitivity indices are calculated using the same direction-coded 23-variable matrix:

- 8 components retained by the 70% cumulative-variance rule
- 6 components retained by the Kaiser criterion

This produces:

- `vuln_index_pca23_70pct`
- `vuln_index_pca23_70pct_z`
- `vuln_index_pca23_kaiser`
- `vuln_index_pca23_kaiser_z`

### 14.4 Direction Anchoring

Capacity indicators are reverse-coded before aggregation so that higher values consistently indicate higher social vulnerability or lower adaptive capacity.

### 14.5 Relation to the Thesis Candidate Set

The `thesis_candidate_23` set is now the active main input set.

### 14.6 Current Interpretation

The current vulnerability index should be understood as:

- a descriptive corridor-specific socio-economic vulnerability measure
- suitable for exposure-vulnerability and protection-related distributive comparisons
- conceptually balanced across three dimensions
- separate from protection/loss variables, which are analysed as flood-risk outcomes

---

## 15. Current Hazard Focus within the Corridor

Status: `NEW / DONE`

Although flood exposure was calculated for all return periods, the current analytical focus is temporarily placed on:

- `RP100` flood share within the corridor

Rationale:

- `RP100` is a widely used reference return period in flood-risk analysis
- it provides a consistent first comparison point for the vulnerability analysis
- it can later be complemented by multi-return-period change metrics once the protection-related simulation data are available

Current RP100 summary outputs:

- `flood_share_rp100`
- RP100 exposure maps
- RP100 summary statistics
- first descriptive association between vulnerability and RP100 exposure

This is an interim analytical step, not the final protection inference stage.

---

## 16. Final Corridor Analysis Layer

Status: `NEW / DONE`

A final corridor-level analytical layer was created by joining:

- corridor municipality geometries
- flood exposure results across all return periods
- PCA component scores
- current vulnerability indices

Outputs include:

- `corridor_analysis_rp100.csv`
- `corridor_wide_pca_rp100_analysis.gpkg`
- `corridor_wide_pca_rp100.rds`

The layer contains:

- `AGS`
- municipality name
- municipality area
- flood area per return period
- flood share per return period
- PCA scores
- main and sensitivity vulnerability indices

---

## 17. Cartographic Outputs for the Corridor

Status: `NEW / DONE`

A first corridor-specific mapping workflow was implemented.

Exported maps currently include:

- study area map showing municipalities intersecting the RP500 corridor
- vulnerability index map for corridor municipalities only
- RP100 exposure map
- PC1-PC4 interpretation maps

Maps were produced in R using:

- `ggplot2`
- `geom_sf()`
- `ggspatial`

Map elements include:

- north arrow
- scale bar
- Elbe overlay
- standardised captions and export settings

---

## 18. Future Protection Inference

Status: `UPDATED / PLANNED`

The current workflow does not yet estimate effective protection directly.

Instead, the present analytical setup is intended as the foundation for the next step:

- exposure curves across `RP10`-`RP500` are available
- vulnerability structure is estimated on the corridor sample
- the next planned step is to use additional simulation or damage-based data to infer effective protection behaviour per municipality

The conceptual logic remains:

- municipalities with low exposure at lower return periods
- but strong exposure increases at higher return periods
- may indicate thresholds at which protection ceases to be effective

Thus, the current corridor-wide exposure-vulnerability workflow is the preparatory stage for the later protection analysis.

---

# Methodological Rationale: GIS-R Interface

The updated workflow deliberately separates responsibilities.

GIS / spatial preprocessing:

- base layer preparation
- basin clipping
- CRS harmonisation
- municipal geometry preparation

R:

- raster-based exposure extraction
- municipality corridor definition
- automated multi-return-period exposure calculation
- quality checks
- socio-economic processing
- vulnerability analysis

This division:

- improves reproducibility
- reduces manual GIS repetition
- avoids unnecessary raster polygonization
- ensures computational feasibility
- allows cleaner scaling to multiple return periods

---

# Mini Changelog: What Fundamentally Changed

Compared to the earlier workflow, the following components were abandoned or superseded:

- BfG vector-based flood zone analysis
- zone 1 / zone 3 residual-risk framing
- PESeta-based direct hazard masking
- protection buffer proxy as central analytical measure

The current workflow now centres on:

- EFAS multi-return-period flood rasters
- RP500-based empirical flood corridor
- municipality-level flood share curves
- later inference of effective protection from simulated return-period differences
