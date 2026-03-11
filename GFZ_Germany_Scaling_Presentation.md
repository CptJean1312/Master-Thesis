# INKAR-based vulnerability index

- Germany scaling check for the GFZ meeting
- March 2026

# 1. What does the thesis do?

- Study area: German part of the Elbe river basin
- Core goal: analyse the social distribution of flood exposure and socio-economic vulnerability
- Justice angle: are more vulnerable municipalities also more exposed to flood hazards?
- Data logic:
  - flood exposure from BfG HQ100-based processing
  - socio-economic vulnerability from municipality-level INKAR indicators
  - protection context currently informed by EU flood-protection data
- Current task for the meeting:
  - assess whether the INKAR-based vulnerability component can be scaled up to Germany

# 2. Why INKAR matters for national scaling

- INKAR 2025 is not Elbe-specific from the start
- Raw dataset structure:
  - long-format CSV
  - 63,426,864 rows
  - 45 spatial-reference categories
  - top-level areas: EU, LRB, SDG, ZOM
- Germany-wide municipal subset used for the check:
  - `Bereich = LRB`
  - `Raumbezug = Gemeinden`
  - 10,978 municipalities
  - 176 municipality-level indicators

# 3. Germany-wide availability on AGS level

- Municipality-level indicators in total: `176`
- Present in all 16 federal states: `173`
- Fully complete on AGS level in all 16 federal states: `32`
- Not present in all 16 states: `3`
  - `q_sach`
  - `q_schluesselzuw`
  - `q_investZ`

These correspond to:

- capital-investment expenditure
- key allocations
- grants for investment support measures

Main message:

- the socio-economic part is broadly available nationwide
- but absolute municipality completeness is much stricter than state-level presence

# 4. Why strict full coverage alone is not enough

If I only keep variables with full AGS-level completeness, the pool becomes too narrow and too structural.

What remains is mostly:

- minijob structure
- broadband availability
- accessibility to central places
- density and structural counts
- land-use context

What largely drops out:

- social transfers
- low income
- long-term unemployment
- ageing and dependency
- single-person households
- everyday health-access indicators

Conclusion:

- `strict full coverage only` is clean in data terms
- but weak in substantive vulnerability terms

# 5. Proposed nationwide shortlist

Best current compromise: variables that

- are present in all 16 federal states
- have very high AGS-level availability
- still represent meaningful vulnerability dimensions

Core shortlist:

- Deprivation:
  - `a_ALGII_SGBII`
  - `a_Unterkunft_SGBII`
  - `a_aloLang`
  - `q_alo_u25_einw`
  - `q_kaufkraft`
  - `a_hheink_niedrig`
- Demographic sensitivity:
  - `a_bev65um`
  - `a_bev75um`
  - `q_abhg_alt`
- Household structure:
  - `q_HH1`
- Accessibility:
  - `m_OEV20_DIST`
  - `m_Q07_HA_DIST`
  - `m_Q01_APO_DIST`

# 6. AGS coverage of the proposed shortlist

![](/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_inkar_overview/map_germany_shortlist_coverage.png)

# 7. Coverage rates of the shortlist variables

![](/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/outputs_inkar_overview/plot_nationwide_shortlist_coverage_bar.png)

# 8. First wide PCA: why I started there

Initial version in `Analyse_CLEAN.R`:

- one global PCA across 52 variables
- intentionally exploratory
- broad indicator pool:
  - deprivation
  - labour market
  - age structure
  - household structure
  - education
  - health provision
  - digital access
  - service distances
  - density and structural indicators
- first 8 PCs retained for the main index
- variance-weighted aggregation
- sign anchored to `ALG II / SGB II`

Why that was useful:

- good first scan of the municipality-level socio-economic space
- good for identifying redundancy and latent structure

# 9. What changed in the redesign

Redesign in `Analyse_PROTOCOL_REDESIGN.R`:

- narrower and theory-led indicator set
- PCA kept, but applied more selectively
- indicators grouped into domains before reduction
- PC1 used within each domain
- domains combined into the redesigned main index
- additional reduced-set PCA kept only as sensitivity check

Redesigned main domains:

- deprivation
- demographic sensitivity
- household / social structure
- accessibility / adaptive capacity

Why this is better for the thesis:

- more interpretable
- less driven by structural overlap
- easier to defend conceptually in a supervisor meeting

# 10. Bottom line for the meeting

My current position:

- yes, an INKAR-based socio-economic vulnerability component can in principle be scaled up to Germany
- no, I would not base that on `strict full coverage only`
- the proposed shortlist is the strongest current compromise between:
  - nationwide availability
  - AGS-level coverage
  - substantive vulnerability relevance
- the stronger scaling bottleneck is not INKAR itself, but the national consistency of hazard and protection data

One caveat that should remain visible:

- `a_BG1P = single-person benefit units` (official INKAR category)
- this does not match the current project label `share_bg_single_parent`
- I would not use that variable in a national shortlist unless the naming and interpretation are cleaned up first
