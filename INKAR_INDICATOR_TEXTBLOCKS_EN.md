# INKAR Indicator Text Blocks

This document translates the final indicator list from `/Users/maxi_161/Desktop/INKAR INDEXE.xlsx` into concise English text blocks. The literature-based justifications below follow the final rationale documented in the Excel sheet.

## Core Indicators for the Main Flood Vulnerability Index

### `exposure_share_alg2_sgb2`

Human-readable name: Share of population receiving ALG II / SGB II benefits.

Dimension: Socioeconomic deprivation.

Literature-based justification: Poverty and welfare dependency consistently increase vulnerability because households with fewer financial resources have lower capacity for disaster preparedness, evacuation and recovery. Rufat et al. identify socioeconomic status as one of the most robust drivers of flood vulnerability. (Rufat et al., 2015)

Transformation: none.

### `exposure_share_sgb2_with_housing_costs`

Human-readable name: Share of SGB II recipients receiving housing support.

Dimension: Socioeconomic deprivation.

Literature-based justification: Housing-related welfare dependence indicates structural economic disadvantage and limited financial buffers. Lower-income households face higher recovery burdens after flood events. (Fekete 2009; Cutter et al. 2003)

Transformation: none.

### `exposure_share_longterm_unemp`

Human-readable name: Share of long-term unemployed population.

Dimension: Socioeconomic deprivation.

Literature-based justification: Long-term unemployment represents persistent labour market exclusion and reduced financial resilience, which reduces the capacity to prepare for and recover from disasters. (Rufat et al., 2015)

Transformation: none.

### `exposure_unemp_u25_per_1000`

Human-readable name: Youth unemployment per 1000 inhabitants.

Dimension: Socioeconomic deprivation.

Literature-based justification: Youth unemployment reflects structural economic fragility and limited labour market integration, which may signal broader regional vulnerability. (Cutter et al., 2003; Rufat et al., 2015)

Transformation: none.

### `exposure_purchasing_power`

Human-readable name: Purchasing power per capita.

Dimension: Socioeconomic deprivation.

Literature-based justification: Lower purchasing power indicates reduced economic resources and lower adaptive capacity of households in hazard situations. Economic resources strongly influence preparedness and recovery potential. (Rufat et al., 2015)

Transformation: inverted (low purchasing power = higher vulnerability).

### `exposure_share_hh_income_low`

Human-readable name: Share of households with low income.

Dimension: Socioeconomic deprivation.

Literature-based justification: Low-income households consistently show higher disaster vulnerability due to limited access to insurance, savings and recovery resources. (Cutter et al., 2003; Rufat et al., 2015)

Transformation: none.

### `exposure_share_age_65plus`

Human-readable name: Share of population aged 65+.

Dimension: Demographic sensitivity.

Literature-based justification: Elderly populations show higher flood vulnerability due to reduced mobility, higher health risks and greater evacuation dependency. (Fekete 2009; Rufat et al., 2015)

Transformation: none.

### `exposure_share_age_75plus`

Human-readable name: Share of population aged 75+.

Dimension: Demographic sensitivity.

Literature-based justification: Advanced age further increases vulnerability due to higher care needs and limited self-reliance during hazard events. (Fekete 2009)

Transformation: none.

### `exposure_old_age_dependency`

Human-readable name: Old-age dependency ratio.

Dimension: Demographic sensitivity.

Literature-based justification: A higher dependency ratio indicates a larger share of elderly relative to working-age population, implying reduced social support capacity during disasters. (Cutter et al., 2003)

Transformation: none.

### `exposure_share_bg_single_parent`

Human-readable name: Share of single-parent households.

Dimension: Household / social structure.

Literature-based justification: Single-parent households often face higher economic stress and reduced internal support networks during disaster response and recovery. (Rufat et al., 2015)

Transformation: none.

### `exposure_share_single_households`

Human-readable name: Share of single-person households.

Dimension: Household / social structure.

Literature-based justification: Single-person households may have limited social support and emergency assistance networks, increasing vulnerability in crisis situations. (Cutter et al., 2003)

Transformation: none.

### `exposure_dist_public_transport_m`

Human-readable name: Distance to nearest public transport stop.

Dimension: Accessibility / adaptive capacity.

Literature-based justification: Access to transport infrastructure influences evacuation capacity and access to services during disasters. Poor accessibility reduces adaptive capacity. (IPCC AR6 WGII; Rufat et al., 2015)

Transformation: none.

### `exposure_dist_gp_m`

Human-readable name: Distance to nearest general practitioner.

Dimension: Accessibility / adaptive capacity.

Literature-based justification: Access to healthcare services is an important component of disaster resilience, especially for vulnerable populations such as elderly or chronically ill residents. (Rufat et al., 2015)

Transformation: none.

### `exposure_dist_pharmacy_m`

Human-readable name: Distance to nearest pharmacy.

Dimension: Accessibility / adaptive capacity.

Literature-based justification: Access to medication and health infrastructure affects recovery and coping capacity during and after flood events. (IPCC AR6 WGII)

Transformation: none.

## Optional Variables for Sensitivity / Robustness Checks

### `exposure_doctors_total_per_1000`

Human-readable name: Doctors per 1000 inhabitants.

Dimension: Health / adaptive capacity.

Literature-based justification: Healthcare capacity can influence disaster resilience and post-disaster recovery, especially for vulnerable populations. (Rufat et al., 2015)

Transformation: inverted (fewer doctors = higher vulnerability).

### `exposure_dist_supermarket_m`

Human-readable name: Distance to supermarket.

Dimension: Accessibility / adaptive capacity.

Literature-based justification: Access to food supply infrastructure influences everyday resilience and recovery capacity following disruptions. (IPCC AR6 WGII)

Transformation: none.

### `exposure_share_bb_100mbit`

Human-readable name: Share of households with broadband access.

Dimension: Digital access / adaptive capacity.

Literature-based justification: No explicit literature note was entered in the Excel table. Conceptually, this variable is retained only as an optional robustness indicator rather than as a core index component.

Transformation: if used, invert so that lower broadband coverage indicates higher vulnerability.

## Excluded Variables

### `exposure_pop_density_per_km2`

Human-readable name: Population density.

Dimension: Spatial structure.

Literature-based justification: Population density reflects settlement patterns rather than social vulnerability directly. Effects may be positive or negative depending on infrastructure and services. (Cutter et al., 2003)

Reason for exclusion: treated as contextual variable rather than vulnerability indicator.

### `exposure_employment_density_per_km2`

Human-readable name: Employment density.

Dimension: Economic structure.

Literature-based justification: Employment density reflects regional economic structure and urbanisation rather than household-level vulnerability. (Rufat et al., 2015)

Reason for exclusion: structural indicator.

### `exposure_trade_tax_per_capita`

Human-readable name: Trade tax revenue per capita.

Dimension: Fiscal capacity.

Literature-based justification: Municipal fiscal indicators reflect local government finances rather than individual household vulnerability. (Cutter et al., 2003)

Reason for exclusion: not household-level vulnerability.

### `exposure_tax_revenue_total`

Human-readable name: Total municipal tax revenue.

Dimension: Fiscal capacity.

Literature-based justification: Municipal revenue is strongly driven by economic structure and population size, not necessarily by social vulnerability.

Reason for exclusion: structural indicator.

### `exposure_share_hh_income_medium`

Human-readable name: Share of medium-income households.

Dimension: Income structure.

Literature-based justification: Medium-income share is redundant when low-income indicators are already included in the index.

Reason for exclusion: redundancy.

### `exposure_share_hh_income_high`

Human-readable name: Share of high-income households.

Dimension: Income structure.

Literature-based justification: High-income share acts mainly as the inverse of low-income indicators and therefore does not add additional explanatory power.

Reason for exclusion: redundancy.

### `exposure_students_total_per_1000`

Human-readable name: Students per 1000 inhabitants.

Dimension: Demography.

Literature-based justification: Student populations may be either highly mobile or socioeconomically privileged; the direction of vulnerability is unclear. (Rufat et al., 2015)

Reason for exclusion: ambiguous indicator.

### `exposure_students_18_25_per_1000`

Human-readable name: Students aged 18-25.

Dimension: Demography.

Literature-based justification: Same conceptual ambiguity as general student population.

Reason for exclusion: ambiguous.

### `exposure_students_fh_per_1000`

Human-readable name: University students.

Dimension: Demography.

Literature-based justification: Student populations do not consistently correlate with disaster vulnerability.

Reason for exclusion: ambiguous.

### `exposure_migration_balance`

Human-readable name: Migration balance.

Dimension: Demography.

Literature-based justification: Migration flows capture regional dynamics rather than vulnerability; interpretation is context dependent.

Reason for exclusion: indirect indicator.

### `exposure_natural_pop_change`

Human-readable name: Natural population change.

Dimension: Demography.

Literature-based justification: Population growth or decline does not directly reflect vulnerability.

Reason for exclusion: indirect indicator.
