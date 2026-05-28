# Chapter 4 Methods Literature Sweep

Purpose:
This file collects the most useful method-oriented literature for writing Chapter 4. It is based on a full-library sweep across the extracted text files, with emphasis on data, spatial scale, exposure construction, vulnerability index design, PCA, uncertainty, and residual-risk / protection interpretation.

Use:
This is not the final reference list. It is a writing aid for selecting citations in the methods chapter.

---

# 1. Spatial Scale and Flood-Risk Assessment

## de Moel et al. (2015) - Flood risk assessments at different spatial scales

Best use in Chapter 4:

- justify why scale matters in flood-risk assessment
- support the municipality-level scale choice
- support the statement that method, interpretation, and uncertainty depend on scale

Useful extracted hooks:

- spatial scale, method, and use of assessments are closely linked
- validation of flood-risk assessments is sparse at all scales
- uncertainty is relevant across scales
- return-period assumptions and spatial uniformity can affect flood-risk estimates

Where to use:

- `4.1 Data Sources`
- `4.2 Spatial Preprocessing and Study Area Definition`
- `4.4 Exposure Quality Checks`
- `4.9 Statistical and Spatial Analysis`

---

## Poussard et al. (2021) - Environmental Inequalities in Flood Exposure: A Matter of Scale

Best use in Chapter 4:

- support the scale-sensitive exposure inequality approach
- justify municipality-level comparison of exposure and socio-economic conditions
- support the use of exposure classifications or municipality-level exposure shares

Useful extracted hooks:

- flood exposure inequalities vary across spatial scales
- exposure and socio-economic status relationships can change when moving from broader regions to smaller units
- municipality-level exposure can be compared across socio-economic groups or deprivation classes

Where to use:

- `4.2 Spatial Preprocessing and Study Area Definition`
- `4.6 Indicator Screening and Candidate Vulnerability Variables`
- `4.9 Statistical and Spatial Analysis`

---

## Sairam et al. (2021) - Process-Based Flood Risk Assessment for Germany

Best use in Chapter 4:

- support the importance of spatial consistency in German flood-risk assessment
- connect the thesis workflow to German large-scale flood-risk modelling
- justify harmonised spatial processing and consistent exposure indicators

Where to use:

- `4.1 Data Sources`
- `4.2 Spatial Preprocessing and Study Area Definition`
- `4.3 Flood Hazard and Exposure Construction`

---

# 2. Flood Hazard, Return Periods, and Exposure Construction

## Rojas et al. (2013) - Climate change and river floods in the European Union

Best use in Chapter 4:

- support the use of multi-return-period flood hazard information
- support European flood-risk assessment framing
- support the interpretation that uncertainty increases especially for higher return periods

Useful extracted hooks:

- flood hazard was considered for different return periods
- return periods include values up to 500 years in European-scale assessment
- risk assessment combines hazard, exposure, and vulnerability
- large-scale approaches involve uncertainty in return levels, hydrology, and exposed assets

Where to use:

- `4.1 Data Sources`
- `4.3 Flood Hazard and Exposure Construction`
- `4.4 Exposure Quality Checks`

---

## Alfieri et al. / Dottori et al. style EFAS/JRC hazard-map sources

Best use in Chapter 4:

- final official source for EFAS / European flood hazard rasters
- define raster content, return periods, resolution, and whether protection assumptions are included

Status:

- exact dataset citation still needs to be confirmed from source metadata.

Where to use:

- `4.1 Data Sources`
- `4.3 Flood Hazard and Exposure Construction`

Placeholder:

`[INSERT EXACT EFAS/JRC FLOOD HAZARD MAP CITATION]`

---

## Modelling the socio-economic impact of river floods in Europe

Best use in Chapter 4:

- support the hazard-exposure-vulnerability framing
- support the use of several return periods
- support the link between flood extent/depth, exposed people/assets, and risk assessment

Useful extracted hooks:

- flood risk assessment highlights hazard, exposure, and vulnerability
- impacts were estimated for selected return periods
- uncertainty affects flood depth, extent, land use, and model chain stages

Where to use:

- `4.3 Flood Hazard and Exposure Construction`
- `4.4 Exposure Quality Checks`
- `4.9 Statistical and Spatial Analysis`

---

# 3. Social Vulnerability, Indicator Selection, and PCA

## Fekete (2009) - Validation of a social vulnerability index in context to river-floods in Germany

Best use in Chapter 4:

- support PCA / factor-analysis-based vulnerability construction in the German river-flood context
- support sensitivity and uncertainty awareness
- support the use of subnational socio-demographic indicators

Useful extracted hooks:

- principal component analysis used for data reduction
- factor analysis is sensitive to variable selection and exclusion
- validation and sensitivity analysis are important for vulnerability indices

Where to use:

- `4.5 Socio-Economic Data Processing`
- `4.6 Indicator Screening and Candidate Vulnerability Variables`
- `4.7 PCA-Based Vulnerability Index Construction`

---

## Fekete (2010) / Assessment of Social Vulnerability for River-Floods in Germany

Best use in Chapter 4:

- support vulnerability assessment with exposure, susceptibility, and capacities
- support the idea that social vulnerability indicators can be linked with exposure information
- support transparent standardisation, aggregation, and sensitivity logic

Useful extracted hooks:

- exposure can be calculated as the percentage of inundated area
- social vulnerability index construction involves standardisation, weighting, and aggregation
- exact weighting and aggregation carry uncertainty and should be interpreted cautiously
- constructing a social vulnerability index is a trade-off between evidence and comparability

Where to use:

- `4.6 Indicator Screening and Candidate Vulnerability Variables`
- `4.7 PCA-Based Vulnerability Index Construction`
- `4.8 Multi-Return-Period Exposure Curves`

---

## Rufat et al. (2015) - Social vulnerability to floods: Review of case studies and implications for measurement

Best use in Chapter 4:

- support the multidimensional and context-dependent character of vulnerability measurement
- justify transparent indicator selection
- justify caution around scale, weighting, and aggregation

Useful extracted hooks:

- social vulnerability modelling research addresses scale, temporal change, uncertainty, and validation
- documentation of analysis scale, weighting, and aggregation is often weak
- vulnerability index construction involves variable selection, scale, internal structure, weighting, and aggregation

Where to use:

- `4.6 Indicator Screening and Candidate Vulnerability Variables`
- `4.7 PCA-Based Vulnerability Index Construction`

---

## Schmidtlein et al. (2008) - A Sensitivity Analysis of the Social Vulnerability Index

Best use in Chapter 4:

- support sensitivity checks for vulnerability index construction
- justify comparing alternative PCA input sets
- justify caution about indicator choice and component retention

Where to use:

- `4.6 Indicator Screening and Candidate Vulnerability Variables`
- `4.7 PCA-Based Vulnerability Index Construction`

---

## Tate (2013) - Uncertainty Analysis for a Social Vulnerability Index

Best use in Chapter 4:

- support the statement that vulnerability index outcomes depend on methodological choices
- justify sensitivity analysis, robustness checks, and transparent index construction

Where to use:

- `4.7 PCA-Based Vulnerability Index Construction`
- `4.9 Statistical and Spatial Analysis`

---

## Tate (2019) - How Valid Are Social Vulnerability Models?

Best use in Chapter 4:

- justify treating vulnerability indices as descriptive analytical tools, not self-validating truth
- support the need for validation / robustness logic

Useful extracted hooks:

- social vulnerability models are increasingly used in hazard mitigation and recovery planning
- many social vulnerability indexes are inductive and use factor analysis / PCA
- empirical validation is limited
- caution is needed when interpreting social vulnerability indexes

Where to use:

- `4.6 Indicator Screening and Candidate Vulnerability Variables`
- `4.7 PCA-Based Vulnerability Index Construction`

---

# 4. Exposure Inequality and Statistical Comparison

## Qiang (2019) - Disparities of population exposed to flood hazards in the United States

Best use in Chapter 4:

- support comparing flood exposure with social disadvantage
- support the idea that flood exposure is assessed by spatial intersection / overlay
- support disparity-oriented statistical analysis

Useful extracted hooks:

- flood exposure is usually assessed by intersection of hazard zones with exposed populations or places
- economically disadvantaged populations can be more likely to reside in flood zones
- local and national scales can produce different findings

Where to use:

- `4.3 Flood Hazard and Exposure Construction`
- `4.9 Statistical and Spatial Analysis`

---

## Odersky and Löffler (2024) - Differential Exposure to Climate Change?

Best use in Chapter 4:

- support German local-scale exposure inequality analysis
- useful especially in results/discussion, but can motivate the statistical comparison

Where to use:

- `4.9 Statistical and Spatial Analysis`
- later Chapter 5 / Chapter 6

---

# 5. Residual Risk, Protection, and Exposure Curves

## Tobin (1995) - The Levee Love Affair

Best use in Chapter 4:

- justify why protection should not be treated as absolute safety
- support the interpretation of protection as modifying exposure and residual risk

Where to use:

- `4.8 Multi-Return-Period Exposure Curves`
- protection placeholder

---

## Serra-Llobet et al. (2022) - Managing residual flood risk behind levees

Best use in Chapter 4:

- support cautious treatment of protection and residual risk
- support not using binary protected/unprotected logic
- support later protection-related interpretation

Useful extracted hooks:

- areas behind levees are still subject to residual risk of breach or overtopping
- flood maps differ in whether they show areas behind levees as flood-prone
- managing residual risk requires acknowledging areas behind protection

Where to use:

- `4.8 Multi-Return-Period Exposure Curves`
- `4.9 Statistical and Spatial Analysis`
- final protection methods placeholder

---

## Fu et al. (2023) - Managing rising residual flood risk

Best use in Chapter 4:

- support residual-risk framing and caution around protection
- useful for later protection interpretation

Where to use:

- `4.8 Multi-Return-Period Exposure Curves`
- protection placeholder

---

## Ferdous et al. (2020) - Structural flood protection, population density, and flood mortality

Best use in Chapter 4:

- support the broader claim that structural protection can shape exposure and settlement dynamics
- useful for interpreting why exposure curves may reveal protection-related dynamics but do not directly prove protection

Useful extracted hooks:

- levees protect against frequent flooding but can contribute to severe losses
- structural protection can attract more people and assets into flood-prone areas
- the levee effect / residual risk / safety dilemma is central

Where to use:

- `4.8 Multi-Return-Period Exposure Curves`
- later Chapter 6 discussion

---

# 6. Data Sources That Need Non-Literature Citations

These are official dataset citations rather than academic literature.

## BKG VG250

Use for:

- municipality boundaries

Need:

- official citation / access date

---

## BBSR INKAR

Use for:

- socio-economic indicators

Current citation candidate:

- Laufende Raumbeobachtung des BBSR - INKAR, Ausgabe 07/2025. Hrsg.: Bundesinstitut für Bau-, Stadt- und Raumforschung (BBSR), Bonn.

Need:

- final bibliography formatting

---

## EFAS / EU Flood Hazard Rasters

Use for:

- multi-return-period flood hazard rasters

Need:

- exact dataset title
- provider
- version / publication year
- URL
- access date
- documentation on protection assumptions

This is the most important missing data-source citation.

---

## GISD

Use only if later included.

Dataset citation:

- Michalski et al. (2025), Zenodo, DOI: 10.5281/zenodo.14781119

Recommended publication:

- Michalski et al. (2022), Journal of Health Monitoring, DOI: 10.25646/10640

---

# 7. Strong Chapter 4 Citation Logic

Recommended minimal citation spine:

- `de Moel et al. (2015)` for scale and flood-risk assessment
- `Rojas et al. (2013)` plus final EFAS/JRC source for multi-RP flood hazard
- `Sairam et al. (2021)` for German large-scale flood-risk assessment
- `Fekete (2009, 2010)` for German river-flood vulnerability and PCA/factor-analysis logic
- `Rufat et al. (2015)` for vulnerability measurement caveats
- `Tate (2013)` and `Schmidtlein et al. (2008)` for uncertainty/sensitivity in vulnerability indices
- `Poussard et al. (2021)` for exposure inequality and scale
- `Serra-Llobet et al. (2022)` and `Tobin (1995)` for residual risk and protection caution

That set is enough to make Chapter 4 methodologically grounded without overloading the text.
