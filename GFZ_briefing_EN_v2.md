# INKAR Indicator

# 1. What does the Master's thesis do?

The thesis examines the social distribution of flood exposure and socio-economic vulnerability in the German part of the Elbe river basin, with a particular focus on questions of flood-protection justice.

In substantive terms, the project combines three components:

- a flood-hazard and exposure component based on BfG HQ100 flood data
- a socio-economic vulnerability component based on municipality-level INKAR indicators
- a protection-context component, currently informed by EU flood-protection data because the BfG Zone 3 information is not sufficiently consistent for the Elbe case

The central analytical question is whether municipalities with higher socio-economic vulnerability are also more exposed to flood hazards, and whether this exposure is distributed unevenly across different social contexts.

Methodologically, the project combines:

- GIS-based spatial preprocessing and area calculations
- municipality-level socio-economic data harmonisation
- principal component analysis for vulnerability structure
- index construction
- regression and spatial diagnostics

The current redesign keeps the overall project goal unchanged, but improves the socio-economic operationalisation by moving from a very broad exploratory PCA to a more selective and theory-led vulnerability design.

# 2. What is the INKAR dataset in this project?

For the socio-economic side of the thesis, the BBSR INKAR 2025 dataset was used. The raw file is a large CSV in long format.

- raw structure: long format
- raw size: `63,426,864` rows
- areas in the raw dataset: `EU`, `LRB`, `SDG`, `ZOM`
- number of spatial-reference categories: `45`

For this project, the dataset is therefore not Elbe-specific from the beginning. It only becomes Elbe-specific through the later spatial filtering.

# 3. INKAR: municipal indicators available in all 16 federal states

## Scope

This overview is based on the original INKAR 2025 raw CSV on municipality level after filtering to:

- `Bereich = LRB`
- `Raumbezug = Gemeinden`
- all 16 federal-state AGS prefixes (`01` to `16`)

Availability definition:

- `present in all 16 federal states`: an indicator has at least one municipality with a valid latest value in each of the 16 states
- `full coverage in all 16 federal states`: the latest value is available for every municipality in every state

## Headline result

- municipality-level indicators in total: `176`
- indicators present in all 16 federal states: `173`
- indicators with full municipality coverage in all 16 federal states: `32`
- indicators not present in all 16 federal states: `3`

The three indicators not present in all 16 states are all from public-finance reporting:

- `q_sach` - `Ausgaben für Sachinvestitionen`
- `q_schlüsselzuw` - `Schlüsselzuweisungen`
- `q_investZ` - `Zuweisungen für Investitionsfördermaßnahmen`

## Broad socio-economic dimensions

- `Labour market and employment`: `37` indicators, of which `7` have full coverage
- `Demography and household structure`: `35` indicators, of which `2` have full coverage
- `Health, services and accessibility`: `33` indicators, of which `6` have full coverage
- `Housing and land-use context`: `27` indicators, of which `8` have full coverage
- `Income, deprivation and social transfers`: `22` indicators, of which `0` have full coverage
- `Basic population and structural counts`: `9` indicators, of which `9` have full coverage
- `Education`: `7` indicators, of which `0` have full coverage
- `Economic structure and tourism`: `3` indicators, of which `0` have full coverage

Interpretation for the meeting:

- on the socio-economic side, the INKAR framework is in principle transferable across Germany
- the stronger limitation is not the existence of the indicators, but the completeness of municipality-level data and, even more importantly, the hazard and protection side of the overall flood-justice framework

# 4. Which INKAR subset was used for the thesis?

First, the municipal INKAR subset relevant for the project was defined by filtering:

- `Bereich == "LRB"`
- `Raumbezug == "Gemeinden"`

This municipality-level subset contains:

- `26,879,954` rows
- `10,978` municipalities
- `176` indicators
- reference years from `1995` to `2024`

Then the dataset was restricted to the AGS prefixes relevant for the Elbe basin study area.

The resulting Elbe subset contains:

- `16,286,960` rows
- `6,699` municipalities
- `176` indicators
- reference years from `1995` to `2024`

# 5. Variable naming logic in the project

Three naming layers should be kept apart clearly:

1. original raw INKAR code, for example `a_ALGII_SGBII`
2. cleaned internal project name, for example `share_alg2_sgb2`
3. later joined analysis name after the spatial merge, for example `exposure_share_alg2_sgb2`

For the meeting, the safest presentation order is:

- first the original INKAR code
- then the official INKAR indicator name
- then, if needed, the internal project label

# 6. INKAR variables used for the first wide PCA in Version 1

The first wide PCA in `Analyse_CLEAN.R` used a broad municipality-level variable pool. It was intentionally exploratory and much wider than the final redesign.

It included indicators from:

- social benefits and deprivation
- labour market disadvantage
- income and fiscal structure
- age structure and dependency ratios
- household structure
- education and student indicators
- health-care provision
- digital infrastructure
- distance-to-services indicators
- density and settlement-structure indicators

In total, the first wide PCA used `52` variables in one global PCA.

Why this was done at the beginning:

- to scan the broader socio-economic structure without narrowing the concept too early
- to identify redundancy and latent clusters in the indicator space
- to see whether a broad data-driven PCA would yield a meaningful first vulnerability pattern

Method in the first wide PCA:

- median imputation for missing values
- z-standardisation of all variables
- one single global PCA across the full variable set
- retention of the first `8` principal components for the main index
- variance-weighted aggregation of these components
- sign anchoring using `ALG II / SGB II` so that higher values indicate higher vulnerability

# 7. Kept variables in the redesigned main index

## `a_ALGII_SGBII`

- official INKAR indicator: `ALG II-Leistungen an SGBII`
- project label: `share_alg2_sgb2`
- human-readable name: Share of population receiving ALG II / SGB II benefits
- dimension: Socioeconomic deprivation
- literature-based justification: Poverty and welfare dependency consistently increase vulnerability because households with fewer financial resources have lower capacity for disaster preparedness, evacuation and recovery. Rufat et al. identify socioeconomic status as one of the most robust drivers of flood vulnerability. (Rufat et al., 2015)
- transformation: none

## `a_Unterkunft_SGBII`

- official INKAR indicator: `Leistungen für Unterkunft an SGBII`
- project label: `share_sgb2_with_housing_costs`
- human-readable name: Share of SGB II recipients receiving housing support
- dimension: Socioeconomic deprivation
- literature-based justification: Housing-related welfare dependence indicates structural economic disadvantage and limited financial buffers. Lower-income households face higher recovery burdens after flood events. (Fekete 2009; Cutter et al. 2003)
- transformation: none

## `a_aloLang`

- official INKAR indicator: `Langzeitarbeitslose`
- project label: `share_longterm_unemp`
- human-readable name: Share of long-term unemployed population
- dimension: Socioeconomic deprivation
- literature-based justification: Long-term unemployment represents persistent labour market exclusion and reduced financial resilience, which reduces the capacity to prepare for and recover from disasters. (Rufat et al., 2015)
- transformation: none

## `q_alo_u25_einw`

- official INKAR indicator: `Jüngere Arbeitslose`
- project label: `unemp_u25_per_1000`
- human-readable name: Youth unemployment per 1000 inhabitants
- dimension: Socioeconomic deprivation
- literature-based justification: Youth unemployment reflects structural economic fragility and limited labour market integration, which may signal broader regional vulnerability. (Cutter et al., 2003; Rufat et al., 2015)
- transformation: none

## `q_kaufkraft`

- official INKAR indicator: `Kaufkraft`
- project label: `purchasing_power`
- human-readable name: Purchasing power per capita
- dimension: Socioeconomic deprivation
- literature-based justification: Lower purchasing power indicates reduced economic resources and lower adaptive capacity of households in hazard situations. Economic resources strongly influence preparedness and recovery potential. (Rufat et al., 2015)
- transformation: inverted (low purchasing power = higher vulnerability)

## `a_hheink_niedrig`

- official INKAR indicator: `Haushalte mit niedrigem Einkommen`
- project label: `share_hh_income_low`
- human-readable name: Share of households with low income
- dimension: Socioeconomic deprivation
- literature-based justification: Low-income households consistently show higher disaster vulnerability due to limited access to insurance, savings and recovery resources. (Cutter et al., 2003; Rufat et al., 2015)
- transformation: none

## `a_bev65um`

- official INKAR indicator: `Einwohner 65 Jahre und älter`
- project label: `share_age_65plus`
- human-readable name: Share of population aged 65+
- dimension: Demographic sensitivity
- literature-based justification: Elderly populations show higher flood vulnerability due to reduced mobility, higher health risks and greater evacuation dependency. (Fekete 2009; Rufat et al., 2015)
- transformation: none

## `a_bev75um`

- official INKAR indicator: `Einwohner 75 Jahre und älter`
- project label: `share_age_75plus`
- human-readable name: Share of population aged 75+
- dimension: Demographic sensitivity
- literature-based justification: Advanced age further increases vulnerability due to higher care needs and limited self-reliance during hazard events. (Fekete 2009)
- transformation: none

## `q_abhg_alt`

- official INKAR indicator: `Abhängigenquote Alte`
- project label: `old_age_dependency`
- human-readable name: Old-age dependency ratio
- dimension: Demographic sensitivity
- literature-based justification: A higher dependency ratio indicates a larger share of elderly relative to working-age population, implying reduced social support capacity during disasters. (Cutter et al., 2003)
- transformation: none

## `a_BG1P`

- official INKAR indicator: `Einpersonen-Bedarfsgemeinschaften`
- current project label: `share_bg_single_parent`
- current human-readable project label: Share of single-parent households
- dimension: Household / social structure
- literature-based justification used in the project: Single-parent households often face higher economic stress and reduced internal support networks during disaster response and recovery. (Rufat et al., 2015)
- transformation: none
- critical note: the official INKAR indicator name and the current project rename do not match. This variable should not simply be presented as a single-parent-household indicator.

## `q_HH1`

- official INKAR indicator: `Einpersonenhaushalte`
- project label: `share_single_households`
- human-readable name: Share of single-person households
- dimension: Household / social structure
- literature-based justification: Single-person households may have limited social support and emergency assistance networks, increasing vulnerability in crisis situations. (Cutter et al., 2003)
- transformation: none

## `m_OEV20_DIST`

- official INKAR indicator: `Entfernung zur ÖV Haltestelle`
- project label: `dist_public_transport_m`
- human-readable name: Distance to nearest public transport stop
- dimension: Accessibility / adaptive capacity
- literature-based justification: Access to transport infrastructure influences evacuation capacity and access to services during disasters. Poor accessibility reduces adaptive capacity. (IPCC AR6 WGII; Rufat et al., 2015)
- transformation: none

## `m_Q07_HA_DIST`

- official INKAR indicator: `Entfernung zum Hausarzt`
- project label: `dist_gp_m`
- human-readable name: Distance to nearest general practitioner
- dimension: Accessibility / adaptive capacity
- literature-based justification: Access to healthcare services is an important component of disaster resilience, especially for vulnerable populations such as elderly or chronically ill residents. (Rufat et al., 2015)
- transformation: none

## `m_Q01_APO_DIST`

- official INKAR indicator: `Entfernung zur Apotheke`
- project label: `dist_pharmacy_m`
- human-readable name: Distance to nearest pharmacy
- dimension: Accessibility / adaptive capacity
- literature-based justification: Access to medication and health infrastructure affects recovery and coping capacity during and after flood events. (IPCC AR6 WGII)
- transformation: none

# 8. Optional variables for sensitivity / robustness checks

## `q_ärzte_bev`

- official INKAR indicator: `Ärzte`
- project label: `doctors_total_per_1000`
- human-readable name: Doctors per 1000 inhabitants
- dimension: Health / adaptive capacity
- literature-based justification: Healthcare capacity can influence disaster resilience and post-disaster recovery, especially for vulnerable populations. (Rufat et al., 2015)
- transformation: inverted (fewer doctors = higher vulnerability)

## `m_G02_SUP_DIST`

- official INKAR indicator: `Entfernung zum Supermarkt/Discounter`
- project label: `dist_supermarket_m`
- human-readable name: Distance to supermarket
- dimension: Accessibility / adaptive capacity
- literature-based justification: Access to food supply infrastructure influences everyday resilience and recovery capacity following disruptions. (IPCC AR6 WGII)
- transformation: none

## `a_bb_100Mbits`

- official INKAR indicator: `Bandbreitenverfügbarkeit mindestens 100 Mbit/s`
- project label: `share_bb_100mbit`
- human-readable name: Share of households with broadband access
- dimension: Digital access / adaptive capacity
- literature-based justification: No explicit literature note was entered in the Excel table. Conceptually, this variable is retained only as an optional robustness indicator rather than as a core index component.
- transformation: if used, invert so that lower broadband coverage indicates higher vulnerability

# 9. Excluded variables

## `q_bev_fl`

- official INKAR indicator: `Einwohnerdichte`
- project label: `pop_density_per_km2`
- dimension: Spatial structure
- literature-based justification: Population density reflects settlement patterns rather than social vulnerability directly. Effects may be positive or negative depending on infrastructure and services. (Cutter et al., 2003)
- reason for exclusion: treated as contextual variable rather than vulnerability indicator

## `q_bevsva_qkm`

- official INKAR indicator: `Einwohner-Arbeitsplatz-Dichte`
- project label: `employment_density_per_km2`
- dimension: Economic structure
- literature-based justification: Employment density reflects regional economic structure and urbanisation rather than household-level vulnerability. (Rufat et al., 2015)
- reason for exclusion: structural indicator

## `q_gewst_bev`

- official INKAR indicator: `Gewerbesteuer`
- project label: `trade_tax_per_capita`
- dimension: Fiscal capacity
- literature-based justification: Municipal fiscal indicators reflect local government finances rather than individual household vulnerability. (Cutter et al., 2003)
- reason for exclusion: not household-level vulnerability

## `d_steuereinnahme`

- official INKAR indicator: `Steuereinnahmen`
- project label: `tax_revenue_total`
- dimension: Fiscal capacity
- literature-based justification: Municipal revenue is strongly driven by economic structure and population size, not necessarily by social vulnerability.
- reason for exclusion: structural indicator

## `a_hheink_mittel`

- official INKAR indicator: `Haushalte mit mittlerem Einkommen`
- project label: `share_hh_income_medium`
- dimension: Income structure
- literature-based justification: Medium-income share is redundant when low-income indicators are already included in the index.
- reason for exclusion: redundancy

## `a_hheink_hoch`

- official INKAR indicator: `Haushalte mit hohem Einkommen`
- project label: `share_hh_income_high`
- dimension: Income structure
- literature-based justification: High-income share acts mainly as the inverse of low-income indicators and therefore does not add additional explanatory power.
- reason for exclusion: redundancy

## `q_stud`

- official INKAR indicator: `Studierende`
- project label: `students_total_per_1000`
- dimension: Demography / education
- literature-based justification: Student populations may be either highly mobile or socioeconomically privileged; the direction of vulnerability is unclear. (Rufat et al., 2015)
- reason for exclusion: ambiguous indicator

## `q_stud_1825`

- official INKAR indicator: `Studierende je 100 Einwohner 18 bis 25 Jahre`
- project label: `students_18_25_per_1000`
- dimension: Demography / education
- literature-based justification: Same conceptual ambiguity as general student population.
- reason for exclusion: ambiguous

## `q_stud_fh`

- official INKAR indicator: `Studierende an FH`
- project label: `students_fh_per_1000`
- dimension: Demography / education
- literature-based justification: Student populations do not consistently correlate with disaster vulnerability.
- reason for exclusion: ambiguous

## `i_wans`

- official INKAR indicator: `Gesamtwanderungssaldo`
- project label: `migration_balance`
- dimension: Demography
- literature-based justification: Migration flows capture regional dynamics rather than vulnerability; interpretation is context dependent.
- reason for exclusion: indirect indicator

## `i_saldo_nat`

- official INKAR indicator: `Natürlicher Saldo`
- project label: `natural_pop_change`
- dimension: Demography
- literature-based justification: Population growth or decline does not directly reflect vulnerability.
- reason for exclusion: indirect indicator

# 10. Why the first wide PCA was used at the beginning

The first wide PCA was used as an exploratory starting point.

Why this made sense at the beginning:

- it allowed a broad initial scan of the municipality-level socio-economic structure
- it helped identify redundancy and multicollinearity before narrowing the concept
- it gave a first data-driven picture of how deprivation, labour market variables, ageing, household structure, accessibility and infrastructure cluster together

In methodological terms, the first version did the following:

- one global PCA across `52` variables
- median imputation
- z-standardisation
- retention of the first `8` PCs for the main index
- variance-weighted aggregation of the retained PCs
- sign anchoring with `ALG II / SGB II`

# 11. What I learned from the first wide PCA

The first wide PCA was useful, but mainly as an exploratory step.

What it showed:

- the municipality-level indicator space clearly contains strong latent structure
- deprivation, ageing, accessibility and settlement structure are empirically interrelated
- several indicators are redundant or structurally overlapping

What became clear from that first version:

- a single global PCA mixes conceptually different dimensions too strongly
- fiscal and settlement-structure variables can enter the same statistical space as social vulnerability indicators even though they are not the same thing
- the resulting index is harder to defend conceptually in a meeting because the retained components are less transparent substantively

# 12. What the redesigned PCA now does

The redesigned workflow keeps PCA, but uses it in a more selective and theory-led way.

Main logic:

- the variable set is reduced before dimensional reduction
- variables are grouped into domains first
- PCA is applied within domains instead of across one very broad pool
- a separate `selected_pca` remains available as a sensitivity check, but not as the main conceptual backbone

The redesigned main domains are:

- deprivation
- demographic sensitivity
- household / social structure
- accessibility / adaptive capacity

Technical implementation:

- domain-specific PCA for each domain
- PC1 retained within each domain
- sign alignment so that higher values always mean higher vulnerability
- aggregation of domain scores into the redesigned main index
- additional reduced-set `selected_pca` only as robustness check

# 13. What the redesigned PCA says now

The redesigned version gives a cleaner and more defensible interpretation.

Domain PCA summary:

- deprivation domain: `6` variables, PC1 explains about `38.1%`
- age domain: `3` variables, PC1 explains about `90.5%`
- household domain: `2` variables, PC1 explains about `50.1%`
- access domain: `3` variables, PC1 explains about `71.3%`

Validation signals:

- redesigned main index correlation with `ALG II / SGB II`: about `0.54`
- redesign sensitivity PCA correlation with `ALG II / SGB II`: about `0.63`
- redesigned main index correlation with long-term unemployment: about `0.25`

Meeting interpretation:

- the redesign does not reject PCA
- it uses PCA more carefully and more transparently
- the redesigned index is easier to defend because it is more clearly tied to interpretable vulnerability dimensions

# 14. Can the INKAR-based vulnerability framework be scaled up to Germany?

On the socio-economic side, the short answer is: yes, in principle.

Why:

- the INKAR indicators used here are not specific to the Elbe basin
- the municipality-level INKAR framework exists at a Germany-wide scale
- `173` out of `176` municipality-level indicators are present in all 16 federal states
- the main socio-economic logic of deprivation, demographic sensitivity, household structure and accessibility can therefore be transferred beyond the Elbe case

What this means in practice:

- a Germany-wide vulnerability index based on INKAR is feasible
- the redesigned indicator selection is actually better suited for upscaling than the first wide PCA, because it is conceptually narrower and less dependent on highly local or structurally ambiguous variables

What still limits full project upscaling:

- the hazard side would need nationally consistent flood-exposure data
- the protection side would need nationally comparable protection information
- municipality-level completeness is still weaker than simple state-level indicator presence

So the careful formulation for the meeting is:

- the INKAR-based socio-economic vulnerability component is, in principle, scalable to Germany
- the full flood-justice framework is only scalable if the hazard and protection data can also be standardised nationally

# 15. Short meeting message

A short summary line for the meeting could be:

"The first version used a broad exploratory wide PCA to understand the overall structure of the socio-economic indicator space. The redesigned version keeps PCA, but applies it to a reduced and theory-led indicator set organised into domains. This makes the final vulnerability index more interpretable and more defensible, while still retaining a data-reduction step. On the socio-economic side, the framework is in principle transferable across Germany; the stronger limitation lies in the hazard and protection data, not in the INKAR variables themselves."

# 16. Is anything still missing?

Substantively, almost nothing. The one point that should remain visible in the final PDF is the `a_BG1P` caveat, because the official INKAR name and the current project rename do not match.
