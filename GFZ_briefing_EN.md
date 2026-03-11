# INKAR-based Municipal Vulnerability Index

## Briefing note for the GFZ meeting

## 1. Purpose of this note

This briefing note summarises the socio-economic data basis of the thesis, the logic of the vulnerability index, the difference between the first wide PCA approach and the redesigned approach, and the implications for a possible Germany-wide transfer of the socio-economic part of the index.

A practical rule for the meeting is the following:

- when presenting variables, start with the original INKAR code and the official INKAR indicator name
- only afterwards refer to the internal project variable names
- make clear that the socio-economic index is conceptually more transferable across Germany than the hazard and protection components

## 2. The INKAR dataset in this project

The socio-economic part of the thesis is based on the BBSR INKAR 2025 dataset.

Key facts:

- raw format: large CSV in long format
- raw size: 63,426,864 rows
- top-level areas in the raw data: EU, LRB, SDG, ZOM
- number of spatial-reference categories in the raw data: 45

For the municipal analysis, the raw file was filtered to:

- `Bereich = LRB`
- `Raumbezug = Gemeinden`

This Germany-wide municipal subset contains:

- 26,879,954 rows
- 10,978 municipalities
- 176 municipality-level indicators
- reference years ranging from 1995 to 2024

For the Elbe application, the data were then restricted to the AGS prefixes relevant for the Elbe basin study area.

Important implication:

The INKAR part is not inherently Elbe-specific. It becomes Elbe-specific only through the later AGS-based spatial filtering and the join with flood-related layers.

## 3. Germany-wide availability of municipality indicators

After filtering to all 16 federal-state AGS prefixes and keeping the most recent available value for each municipality-indicator combination, the availability pattern is as follows:

- municipality indicators in total: 176
- indicators present in all 16 federal states: 173
- indicators with full municipality coverage in all 16 federal states: 32
- indicators not present in all 16 federal states: 3

The three indicators that are not present in all 16 states are all from public-finance reporting:

- `q_sach` - Ausgaben für Sachinvestitionen
- `q_schlüsselzuw` - Schlüsselzuweisungen
- `q_investZ` - Zuweisungen für Investitionsfördermaßnahmen

Interpretation for the meeting:

- the socio-economic side of the index is, in principle, transferable to Germany as a whole
- the real bottleneck is not basic indicator presence across states, but municipality-level completeness and the availability and comparability of hazard and protection data

## 4. Naming chain: original INKAR names versus project names

Three naming layers have to be distinguished clearly.

1. Original INKAR code in the raw dataset, for example `a_ALGII_SGBII`
2. Internal cleaned name after the first harmonisation step, for example `share_alg2_sgb2`
3. Final joined analysis name after the spatial merge, for example `exposure_share_alg2_sgb2`

For presentation purposes, the original INKAR code and the official INKAR indicator name should always come first.

## 5. Final redesigned index: retained core variables

The redesigned main index uses a reduced and theory-led set of core indicators. The official INKAR names are kept in German here in order to preserve the original source terminology.

| Original INKAR code | Official INKAR indicator name | Internal project name | Analytical role |
| --- | --- | --- | --- |
| `a_ALGII_SGBII` | `ALG II-Leistungen an SGBII` | `share_alg2_sgb2` | core deprivation |
| `a_Unterkunft_SGBII` | `Leistungen für Unterkunft an SGBII` | `share_sgb2_with_housing_costs` | core deprivation |
| `a_aloLang` | `Langzeitarbeitslose` | `share_longterm_unemp` | core deprivation |
| `q_alo_u25_einw` | `Jüngere Arbeitslose` | `unemp_u25_per_1000` | core deprivation |
| `q_kaufkraft` | `Kaufkraft` | `purchasing_power` | core deprivation |
| `a_hheink_niedrig` | `Haushalte mit niedrigem Einkommen` | `share_hh_income_low` | core deprivation |
| `a_bev65um` | `Einwohner 65 Jahre und älter` | `share_age_65plus` | core age sensitivity |
| `a_bev75um` | `Einwohner 75 Jahre und älter` | `share_age_75plus` | core age sensitivity |
| `q_abhg_alt` | `Abhängigenquote Alte` | `old_age_dependency` | core age sensitivity |
| `a_BG1P` | `Einpersonen-Bedarfsgemeinschaften` | `share_bg_single_parent` | current project label is likely misleading |
| `q_HH1` | `Einpersonenhaushalte` | `share_single_households` | core household structure |
| `m_OEV20_DIST` | `Entfernung zur ÖV Haltestelle` | `dist_public_transport_m` | core accessibility |
| `m_Q07_HA_DIST` | `Entfernung zum Hausarzt` | `dist_gp_m` | core accessibility |
| `m_Q01_APO_DIST` | `Entfernung zur Apotheke` | `dist_pharmacy_m` | core accessibility |

## 6. Optional sensitivity variables in the redesign

These variables were not kept in the core index but remain available for robustness checks.

| Original INKAR code | Official INKAR indicator name | Internal project name | Analytical role |
| --- | --- | --- | --- |
| `q_ärzte_bev` | `Ärzte` | `doctors_total_per_1000` | optional health-capacity proxy |
| `m_G02_SUP_DIST` | `Entfernung zum Supermarkt/Discounter` | `dist_supermarket_m` | optional accessibility proxy |
| `a_bb_100Mbits` | `Bandbreitenverfügbarkeit mindestens 100 Mbit/s` | `share_bb_100mbit` | optional digital-access proxy |

## 7. Variables not retained in the final redesigned index

These are the variables that should still appear in the document with their original names, because they were actively considered but not retained in the final redesigned index.

| Original INKAR code | Official INKAR indicator name | Internal project name | Why it was not retained |
| --- | --- | --- | --- |
| `q_bev_fl` | `Einwohnerdichte` | `pop_density_per_km2` | settlement structure rather than vulnerability in the narrow sense |
| `q_bevsva_qkm` | `Einwohner-Arbeitsplatz-Dichte` | `employment_density_per_km2` | structural indicator, not household-level vulnerability |
| `q_gewst_bev` | `Gewerbesteuer` | `trade_tax_per_capita` | municipal fiscal structure rather than social vulnerability |
| `d_steuereinnahme` | `Steuereinnahmen` | `tax_revenue_total` | municipal revenue structure rather than household vulnerability |
| `a_hheink_mittel` | `Haushalte mit mittlerem Einkommen` | `share_hh_income_medium` | redundant once low-income households are already included |
| `a_hheink_hoch` | `Haushalte mit hohem Einkommen` | `share_hh_income_high` | largely the inverse of low-income share |
| `q_stud` | `Studierende` | `students_total_per_1000` | conceptually ambiguous for flood vulnerability |
| `q_stud_1825` | `Studierende je 100 Einwohner 18 bis 25 Jahre` | `students_18_25_per_1000` | conceptually ambiguous for flood vulnerability |
| `q_stud_fh` | `Studierende an FH` | `students_fh_per_1000` | conceptually ambiguous for flood vulnerability |
| `i_wans` | `Gesamtwanderungssaldo` | `migration_balance` | demographic dynamics rather than a direct vulnerability marker |
| `i_saldo_nat` | `Natürlicher Saldo` | `natural_pop_change` | indirect demographic context, not a core vulnerability marker |

## 8. Why the first wide PCA was used in the original analysis

The first analysis in `Analyse_CLEAN.R` used a very broad exploratory PCA setup.

Core idea:

- start from a high-dimensional vulnerability space rather than from a tightly reduced indicator set
- allow the data to reveal latent structure across deprivation, labour market, age structure, household structure, health, digital access, accessibility and density-related variables
- build an initial vulnerability index before narrowing the concept more strongly

Method in the original wide PCA:

- 52 municipality-level variables entered one single PCA
- variables were median-imputed and z-standardised
- a global PCA was run on the full set
- the first 8 principal components were retained for the main index
- a variance-weighted composite score was then built from the retained PCs
- the sign of the final index was anchored to `ALG II / SGB II` so that higher values indicate higher vulnerability

Why this was a defensible first step:

- at an early stage, it allowed a broad exploratory scan of the socio-economic structure of the municipalities
- it gave an initial picture of which clusters of variables move together
- it was useful for identifying multicollinearity and latent structure before a stricter theory-led reduction

## 9. What the first wide PCA showed, and why it became problematic

The first wide PCA was useful for exploration, but less convincing as the final design.

What it showed:

- socio-economic disadvantage, labour-market exclusion, ageing, accessibility and settlement structure are all entangled in the municipality data
- a single global PCA can compress a very broad indicator space into a manageable structure
- the broad variable pool was useful to identify which indicators were likely to be redundant or overly structural

What became problematic:

- the components mixed conceptually different things such as deprivation, ageing, urbanity, accessibility and infrastructure
- several indicators were structurally redundant
- fiscal, density and settlement-structure variables risked dominating the index even though they are not the same thing as social vulnerability
- the resulting composite index was harder to explain in a supervisor meeting, because one global PCA makes the substantive meaning of each retained component less transparent

In short:

The first wide PCA was a good exploratory starting point, but it was not the strongest final argument for a thesis chapter that needs clear conceptual justification.

## 10. What the redesigned PCA framework does instead

The redesigned workflow in `Analyse_PROTOCOL_REDESIGN.R` keeps PCA, but in a much more structured way.

Main changes:

- the variable set was reduced to a theory-led core
- indicators were grouped into domains before dimensional reduction
- PCA is no longer treated as one single global black box
- interpretability was prioritised over maximal breadth

The redesigned core domains are:

- deprivation
- demographic sensitivity
- household / social structure
- accessibility / adaptive capacity

Technical logic of the redesign:

- a separate domain PCA is run within each domain
- the first component of each domain is used as the domain score
- domain scores are oriented so that higher values always mean higher vulnerability
- the four domain scores are combined into the main redesigned index
- a second `selected_pca` sensitivity index is still calculated on the reduced indicator set, but only as a robustness check

## 11. What the redesigned PCA says

The redesigned PCA gives a clearer and more defensible message.

Domain PCA summary:

- deprivation domain: 6 variables, PC1 explains about 38.1% of domain variance
- age domain: 3 variables, PC1 explains about 90.5%
- household domain: 2 variables, PC1 explains about 50.1%
- access domain: 3 variables, PC1 explains about 71.3%

Validation signals from the redesign:

- the redesigned main index correlates with `ALG II / SGB II` at about `0.54`
- the redesign sensitivity PCA correlates with `ALG II / SGB II` at about `0.63`
- the redesigned main index correlates with long-term unemployment at about `0.25`
- this means the redesigned index still tracks deprivation signals clearly, but does so in a more interpretable way than the first global PCA

Main substantive message for the meeting:

- the redesign does not reject PCA
- it uses PCA more selectively and more transparently
- the result is a vulnerability measure that is easier to defend conceptually, because it is less driven by indiscriminate breadth and more clearly tied to vulnerability dimensions

## 12. Critical caveat that should be stated explicitly

One project rename should be treated with caution:

- raw INKAR code: `a_BG1P`
- official INKAR name: `Einpersonen-Bedarfsgemeinschaften`
- current project name: `share_bg_single_parent`

These are not the same thing.

For the meeting, this variable should therefore not simply be presented as `single-parent households`.

If it comes up in discussion, the safe formulation is:

- the original INKAR indicator is `a_BG1P = Einpersonen-Bedarfsgemeinschaften`
- the current internal rename is likely misleading and should be corrected in the next clean-up step
- if a true single-parent indicator is needed, a different INKAR variable should be checked instead

## 13. What to emphasise in the meeting

A concise meeting line could be:

"The socio-economic part of the thesis is based on the nationwide INKAR 2025 municipality dataset. The first version used a broad exploratory wide PCA to understand the overall structure of the indicator space. The redesigned version keeps PCA, but reduces the indicator set and applies PCA within theory-led domains, which makes the final vulnerability measure more interpretable and easier to defend. On the socio-economic side, the framework is in principle transferable across Germany; the stronger limitation lies in the hazard and protection data, not in the INKAR variables themselves."

## 14. Is anything still missing?

Substantively, only three things would still matter before circulation as a PDF:

- make sure the `a_BG1P` naming caveat is left visible and not silently glossed over
- keep original INKAR codes and official names for all retained and excluded variables
- keep the old-versus-new PCA comparison short and strategic rather than overly technical

If those points are covered, the document is already a solid briefing paper for the professors and a good basis for the next PowerPoint step.
