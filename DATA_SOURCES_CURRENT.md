# Data Sources Current

Working inventory of datasets mentioned for the current thesis workflow.

Status note:
This document separates current core datasets from background, optional, or superseded datasets. It is a working basis for Chapter 4 and will need final source citations and access dates before submission.

---

# 1. Current Core Spatial Data

## VG250 Administrative Boundaries

Role in thesis:

- municipality boundaries
- unit of analysis
- area denominator for flood shares
- spatial join framework for INKAR and flood exposure

Source:

- Bundesamt für Kartographie und Geodäsie (BKG)
- Verwaltungsgebiete 1:250 000, Stand 31.12.
- URL: https://gdz.bkg.bund.de/index.php/default/verwaltungsgebiete-1-250-000-stand-31-12-vg250-31-12.html

Format / CRS note:

- georeferenced UTM 32
- GeoPackage used in the working pipeline

Processing note:

- municipality polygons were filtered to `Geofaktor_GF = "mit Struktur Land"`
- water-structure geometries were excluded to avoid invalid municipal area denominators

---

## Elbe River Basin District / Basin Polygon

Role in thesis:

- broader hydrological context
- processing mask for clipping EFAS flood rasters
- not the final empirical study population

Source:

- Bundesanstalt für Gewässerkunde (BfG)
- River basin district as ESRI file geodatabase, ETRS89, non-INSPIRE conformant
- URL: https://geoportal.bafg.de/inspire/download/AM/riverBasinDistrict/datasetfeed.xml

Current interpretation:

- the Elbe basin frames the thesis conceptually
- the final empirical study population is the RP500 flood corridor, not all basin municipalities

---

## Elbe / Watercourse Base Layers

Role in thesis:

- cartographic context
- Elbe overlay in maps
- hydrological orientation

Potential sources:

- water body for WFD / Wasserkörper-DE
- URL: https://geoportal.bafg.de/inspire/download/AM/waterBodyForWFD/datasetfeed.xml

Additional related BfG source:

- Catchment as GML / HYP Catchment
- URL: https://geoportal.bafg.de/inspire/download/HY/catchment/datasetfeed.xml

---

# 2. Current Core Flood Hazard Data

## EU / EFAS Multi-Return-Period Flood Hazard Rasters

Role in thesis:

- main hazard input for exposure construction
- `RP500` defines the Elbe flood corridor
- `RP100` serves as the current main exposure benchmark
- all return periods support multi-RP exposure curves

Return periods used:

- `RP10`
- `RP20`
- `RP50`
- `RP100`
- `RP200`
- `RP500`

Raster interpretation:

- positive flood-depth values indicate flooded cells
- `NoData` indicates non-flooded areas
- flood extent is derived from valid flood-depth pixels

Processing:

- reprojected from `EPSG:3035` to `EPSG:25832`
- clipped to the Elbe basin polygon

Methodological role:

- current central flood exposure source
- replaces earlier BfG vector-based flood-zone workflow for the main analysis

Source note:

- official JRC catalogue match: Baugh, Calum; Colonese, Juan; D'Angelo, Claudia; Dottori, Francesco; Neal, Jeffrey; Prudhomme, Christel; Salamon, Peter (2024): *River flood hazard maps for Europe and the Mediterranean Basin region*. European Commission, Joint Research Centre [Dataset]. DOI: `10.2905/1D128B6C-A4EE-4858-9E34-6210707F3C81`. PID: http://data.europa.eu/89h/1d128b6c-a4ee-4858-9e34-6210707f3c81
- official dataset page describes the maps as gridded flood water-depth maps for nine return periods from 1-in-10 to 1-in-500 years, created as part of the Copernicus Emergency Management Service
- related methodology references listed by the catalogue: Alfieri et al. (2014) and Dottori et al. (2022)
- locally downloaded files follow the catalogue naming pattern used in the workflow: `floodmap_EFAS_RP010_C.tif`, `floodmap_EFAS_RP020_C.tif`, `floodmap_EFAS_RP050_C.tif`, `floodmap_EFAS_RP100_C.tif`, `floodmap_EFAS_RP200_C.tif`, and `floodmap_EFAS_RP500_C.tif`

---

# 3. Current Core Land-Cover Data

## DLR Land Cover DE

Role in thesis:

- RP100 land-cover-sensitive exposure refinement
- composition of modeled flooded land by land-cover group
- flooded share of Artificial Land by municipality
- robustness test for the total municipal-area exposure result

Dataset:

- German Aerospace Center (DLR) (2020): *Land Cover DE - Sentinel-2 - Germany, 2015*
- DOI: https://doi.org/10.15489/1ccmlap3mn39
- official dataset page: https://geoservice.dlr.de/web/datasets/landcover
- official download directory: https://download.geoservice.dlr.de/LCC_DE/files/

Source characteristics:

- 10 m x 10 m categorical raster
- source CRS: `EPSG:3035` (ETRS89 / LAEA Europe)
- derived from Sentinel-2 imagery acquired between 2015 and 2017 and LUCAS 2015 reference data
- seven classes: Artificial Land, Open Soil, High Seasonal Vegetation, High Perennial Vegetation, Low Seasonal Vegetation, Low Perennial Vegetation, and Water
- license: CC BY-NC 4.0

Local external input:

- `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/LANDUSE/data_raw/Land_Cover_DE_2015.tif`
- metadata: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/LANDUSE/data_raw/Land_Cover_DE_2015_ISO_Metadata_8205e1c9-14c3-47a7-ac12-8c890225680a.xml`
- the large raster remains outside the Git repository

Processing and interpretation:

- cropped to the RP500 corridor processing extent
- each source class converted to a binary raster and reprojected to the 100 m EFAS/JRC template in `EPSG:25832` using area-weighted averaging
- class fractions multiplied by cell area and intersected with the RP100 flood mask
- DLR class 1, Artificial Land, used as a built/artificial land-cover proxy
- not interpreted as a cadastral building footprint, residential land alone, exposed population, or exposed asset value
- processed for RP100 in the current thesis workflow

Recommended dataset citation:

- German Aerospace Center (DLR). (2020). *Land Cover DE - Sentinel-2 - Germany, 2015*. https://doi.org/10.15489/1ccmlap3mn39

---

# 4. Current Core Socio-Economic Data

## INKAR

Role in thesis:

- main municipality-level socio-economic dataset
- basis for socio-economic vulnerability indicator selection
- basis for the socio-economic vulnerability index and PCA diagnostics

Description:

INKAR is the interactive online atlas of the Federal Institute for Research on Building, Urban Affairs and Spatial Development. It provides indicators on living conditions in Germany and Europe and supports city-rural comparisons and temporal analyses.

User-provided source description:

- current INKAR edition: `07/2025`
- approximately 600 indicators
- data on consistent territorial status as of `31.12.2023`

Citation suggestion:

- Laufende Raumbeobachtung des BBSR - INKAR, Ausgabe 07/2025. Hrsg.: Bundesinstitut für Bau-, Stadt- und Raumforschung (BBSR), Bonn.

URL:

- https://www.inkar.de

Current use:

- filtered to municipality level
- latest available year per municipality and indicator retained
- pivoted from long to wide format
- joined to corridor municipalities by `AGS`

---

# 5. Optional / Validation Socio-Economic Data

## German Index of Socioeconomic Deprivation (GISD)

Role in thesis:

- optional contextualisation or validation source
- not the main vulnerability index

Source:

- German Index of Socioeconomic Deprivation (GISD)
- URL: https://robert-koch-institut.github.io/German_Index_of_Socioeconomic_Deprivation_GISD/

Dataset citation:

- Michalski, Niels, Soliman, Lola Omar, Reis, Marvin, Tetzlaff, Fabian, Nowossadeck, Enno und Hoebel, Jens (2025): German Index of Socioeconomic Deprivation (GISD), Berlin: Zenodo. DOI: 10.5281/zenodo.14781119

Recommended publication citation:

- Michalski, N., Reis, M., Tetzlaff, F., Herber, M., Kroll, L. E., Hövener, C., Nowossadeck, E., & Hoebel, J. (2022). German Index of Socioeconomic Deprivation (GISD): Revision, Aktualisierung und Anwendungsbeispiele. Journal of Health Monitoring, 7(S5), 2-24. DOI: 10.25646/10640

Description:

GISD measures relative regional socioeconomic deprivation. It was developed at the Robert Koch Institute to make regional socioeconomic inequalities in health visible and to support explanations of regional differences in health opportunities, morbidity, and mortality.

Conceptual dimensions:

- education
- occupation
- income

Current status:

- prepared as optional analytical extension
- should not be presented as part of the final main workflow unless used in final robustness or validation analysis

---

# 6. Background / Potentially Relevant Physical Data

These datasets are relevant background sources or potential extensions, but they are not currently central to the active flood exposure workflow unless explicitly reintroduced.

## Copernicus DEM

Role:

- potential topographic background or hydrological modelling input
- not currently part of the active exposure calculation

Description:

The Copernicus DEM is a Digital Surface Model representing the surface of the Earth including buildings, infrastructure, and vegetation. It is provided in different instances, including `EEA-10`, `GLO-30`, and `GLO-90`.

Acquisition:

- TanDEM-X mission
- data acquired between 2011 and 2015
- older elevation models may have been used to fill gaps
- datasets made available for use in 2019 and maintained until 2026

Important conceptual note:

- Copernicus DEM is a DSM, not a pure DTM/DGM
- it includes buildings, infrastructure, and vegetation
- therefore, it should not be described as a bare-earth DEM unless another DTM product is used

---

## Digital Surface Model of Watercourses Elbe and Lower Havel

Role:

- potential detailed physical background for the Elbe / Lower Havel
- not currently part of active exposure calculation

Dataset:

- Digital surface model of the watercourses Elbe and Lower Havel (Germany), DGM-W Elbe project, DOM Elbe 2022

URL:

- https://zenodo.org/records/17458487

User note:

- `DOM1`

---

## Digital Terrain Model of Watercourses Elbe and Lower Havel

Role:

- potential detailed terrain reference for the Elbe / Lower Havel
- not currently part of active exposure calculation

Dataset:

- Digital terrain model of the watercourses Elbe and Lower Havel (Germany), DGM-W Elbe 2022 project

URL:

- https://zenodo.org/records/17378779

---

## Watershed Boundaries of GRDC Stations

Role:

- possible hydrological context or reference
- not currently part of active thesis workflow

URL:

- https://grdc.bafg.de/products/basin_layers/watershed_boundaries/

---

# 7. Superseded or Contextual Flood and Protection Data

These datasets were explored or may remain conceptually relevant, but they are not the current main exposure workflow.

## BfG Überflutungstiefen-DE

Role in earlier workflow:

- earlier BfG vector-based flood-depth / flood-zone analysis
- included zone and depth-class information

Source:

- BfG OpenData Dataset Überflutungstiefen-DE
- URL: https://geoportal.bafg.de/smartfinderClient/?lang=de#/datasets/iso/06ead55d-6348-45b0-bd80-dbd7d7912fd7

Dataset description:

The dataset contains water depths in flood areas for three flood scenarios:

- low probability
- medium probability
- high probability

Water depths are classified into five classes:

- `0-0.5 m`
- `0.5-1 m`
- `1-2 m`
- `2-4 m`
- `>4 m`

The data distinguish three area types:

- flood area
- informational flood area
- flood-protected area

These two classifications are combined in the attribute `T_class`.

Earlier issue:

- the protection-related `Zone 3` logic was not sufficiently consistent for the current analytical workflow
- the BfG vector workflow has been superseded by the EFAS multi-return-period raster workflow

State-level download feeds:

- Brandenburg: https://geoportal.bafg.de/download/opendata/ueberflutungstiefen/DEBB/datasetfeed.xml
- Hamburg: https://geoportal.bafg.de/download/opendata/ueberflutungstiefen/DEHH/datasetfeed.xml
- Mecklenburg-Vorpommern: https://geoportal.bafg.de/download/opendata/ueberflutungstiefen/DEMV/datasetfeed.xml
- Niedersachsen: https://geoportal.bafg.de/download/opendata/ueberflutungstiefen/DENI/datasetfeed.xml
- Sachsen: https://geoportal.bafg.de/download/opendata/ueberflutungstiefen/DESN/datasetfeed.xml
- Sachsen-Anhalt: https://geoportal.bafg.de/download/opendata/ueberflutungstiefen/DEST/datasetfeed.xml
- Schleswig-Holstein: https://geoportal.bafg.de/download/opendata/ueberflutungstiefen/DESH/datasetfeed.xml

---

## BfG ManMadeObject-DE

Role:

- potential source for actual flood protection infrastructure or hydraulic structures
- considered as alternative protection-related data source
- not currently central in the active exposure pipeline

URL:

- https://geoportal.bafg.de/inspire/download/HY/manmadeobject/datasetfeed.xml

Potential use:

- later contextual or protection-infrastructure analysis, if the thesis returns to infrastructure distribution explicitly

---

## BfG INSPIRE Natural Risk Zones

Role:

- contextual or alternative flood-risk data sources
- not current main workflow

Sources:

- Flood Hazard Area: https://geoportal.bafg.de/inspire/download/NZ/hazardArea_flood/datasetfeed.xml
- Flood Risk Zones: https://geoportal.bafg.de/inspire/download/NZ/riskZone_flood/datasetfeed.xml
- Flood Observed Events: https://geoportal.bafg.de/inspire/download/NZ/observedEvent_flood/datasetfeed.xml

---

## Flood Unit of Management

Role:

- governance / administrative flood-risk context
- not current core analysis layer

Source:

- Bewirtschaftungseinheiten für Hochwasserrisiken-DE
- URL: https://geoportal.bafg.de/inspire/download/AM/floodUnitOfManagement/datasetfeed.xml

---

## Floodplain Dataset

Role:

- contextual floodplain layer
- not current main exposure source

Source:

- Überschwemmungsgebiete-DE
- URL: https://geoportal.bafg.de/inspire/download/AM/floodplain/datasetfeed.xml
