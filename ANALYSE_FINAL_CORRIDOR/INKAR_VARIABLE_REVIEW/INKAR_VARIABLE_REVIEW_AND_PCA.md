# INKAR Corridor Variable Review and PCA Note

Created: `2026-03-19 14:55:55`

## 1. Why this note exists

This note does two things at once:
- it prepares a concrete answer to the supervisor email about potentially counterintuitive INKAR variables;
- it documents a cleaner decision path from the old exploratory wide PCA towards a curated, literature-guided vulnerability PCA for the EFAS corridor sample.

## 2. Main conclusion

- A wide PCA with many variables is still useful as an exploratory diagnostic tool.
- It should not automatically be the final vulnerability index if the variable set contains compositional counterparts, nested age blocks, duplicated access measures, or variables with counterintuitive definitions.
- The original 51-variable set was a broad exploratory screening set. It was not yet fully cleaned with respect to definition conflicts and duplicated constructs.
- For the final thesis, the stronger strategy is: keep the old wide PCA as an exploratory/sensitivity result, but build the main interpretable index from a curated set of variables with clear theoretical roles.

## 3. Direct answer to the email concern

- In the current corridor raw-INKAR build, `a_ALGII_SGBII` and `a_Unterkunft_SGBII` have `825` complete municipality pairs and a correlation of `-0.703`.
- That is a strong negative relationship, but it is not a perfect `-1.0` complement relationship in the current corridor data.
- So the safe interpretation is: they are conceptually overlapping and partly compositional, but they should not be treated as a clean mathematical complement pair without checking the exact official denominator definitions.
- For PCA interpretation, both variables should not remain together in the final curated index, because they would overweight one SGB-II-related dimension.
- My recommendation is to keep `a_ALGII_SGBII` as the more directly interpretable deprivation signal and to drop `a_Unterkunft_SGBII` from the curated main PCA.

## 4. Corridor full-INKAR layer

- Full original-code corridor layer:
  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/gpkg/corridor_full_inkar_original_codes.gpkg`
- One-row-per-municipality CSV:
  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_full_inkar_original_codes.csv`
- Corridor municipalities in geometry layer: `835`
- Municipalities with any INKAR row in the latest-build table: `834`
- Original INKAR indicators available in the corridor build: `176`
- One municipality (`16076094`, Berga-Wünschendorf) remains missing in the raw latest-build table.

## 5. What the old 51-variable PCA was actually doing

- The old 51-variable set was a manually assembled, broad exploratory screening set.
- It mixed poverty/welfare, detailed age composition, household structure, income composition, health-care supply, digital infrastructure, accessibility, and density in one block.
- That is useful to let the data show broad covariance structure.
- But it also means that some concepts appear several times, for example:
  - age structure through multiple nested age shares plus dependency ratios;
  - household income through low, medium, and high shares together;
  - health access through both total physicians and specialist-specific sub-indicators;
  - digital access through several broadband thresholds;
  - SGB-II deprivation through several related but definition-sensitive indicators.
- So the old 51-variable PCA is defensible as an exploratory starting point, but not ideal as the final, most interpretable vulnerability index.

## 6. Initial correlation diagnosis

- Original 51-variable correlation heatmap:
  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/old51_correlation_heatmap.png`
- All-176-variable correlation heatmap:
  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/all176_correlation_heatmap.png`
- Curated-variable correlation heatmap:
  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/curated_correlation_heatmap.png`

The correlation step should come first because it makes duplication visible before PCA loadings are interpreted.

Top correlation pairs from the old 51-variable set are saved here:
- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/old51_top_correlation_pairs.csv`

## 7. PCA comparison

Comparison across model versions:
- `original_51`:  51 variables, 834 municipalities, PC1 variance `16.18%`, cumulative PC1-PC4 `43.40%`, cumulative PC1-PC8 `60.07%`
- `all_176`: 176 variables, 834 municipalities, PC1 variance `15.14%`, cumulative PC1-PC4 `33.01%`, cumulative PC1-PC8 `44.14%`
- `thesis_candidate_23`:  23 variables, 834 municipalities, PC1 variance `24.21%`, cumulative PC1-PC4 `54.81%`, cumulative PC1-PC8 `72.60%`
- `student_sensitivity_26`:  26 variables, 834 municipalities, PC1 variance `22.14%`, cumulative PC1-PC4 `54.59%`, cumulative PC1-PC8 `72.14%`
- `curated_17`:  17 variables, 834 municipalities, PC1 variance `27.65%`, cumulative PC1-PC4 `63.44%`, cumulative PC1-PC8 `83.15%`

Interpretation:
- The all-176 PCA is useful as a broad structural diagnostic, not as the main final index.
- The original 51-variable PCA is a cleaner exploratory subset than all 176, but it still contains several duplicated constructs.
- The thesis-candidate set is the pragmatic middle ground: broader than the strict curated core, but already cleaned from the most obvious definition conflicts and duplicate blocks.
- The student-sensitivity set adds the student concentration indicators only as a robustness run, not as part of the main thesis PCA.
- The curated PCA sacrifices breadth for interpretability and remains the strictest option.

## 8. Thesis-candidate variable set from all available corridor indicators

This is the proposed working set for the thesis right now: broader than the strict 17-variable core, but already cleaned enough to remove the most problematic overlaps.
- `a_ALGII_SGBII` — ALG II-Leistungen an SGBII | coverage `98.80%` | Keep as core candidate in the curated vulnerability PCA.
- `q_newfBGu15_bev` — Kinderarmut | coverage `99.88%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_aloLang` — Langzeitarbeitslose | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_alo_u25_einw` — Jüngere Arbeitslose | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_alo_ü55_einw` — Ältere Arbeitslose | coverage `99.88%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_Minijobs` — Anteil Minijobs an den Beschäftigungsverhältnissen | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_svw` — Beschäftigtenquote | coverage `99.88%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_kaufkraft` — Kaufkraft | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_einkst_bev` — Einkommensteuer | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_hheink_niedrig` — Haushalte mit niedrigem Einkommen | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_bev65um` — Einwohner 65 Jahre und älter | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_bev_0006` — Einwohner unter 6 Jahre | coverage `99.88%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_HH1` — Einpersonenhaushalte | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_hh_kind` — Haushalte mit Kindern | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_ewfBG_allein` — Alleinerziehende erwerbsfähige Leistungsberechtigte | coverage `98.56%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `m_G02_SUP_DIST` — Entfernung zum Supermarkt/Discounter | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q01_APO_DIST` — Entfernung zur Apotheke | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q07_HA_DIST` — Entfernung zum Hausarzt | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_OEV20_DIST` — Entfernung zur ÖV Haltestelle | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_P01_PRIM_DIST` — Entfernung zur Grundschule | coverage `99.88%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_ärzte_bev` — Ärzte | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_bb_50Mbits` — Bandbreitenverfügbarkeit mindestens 50 Mbit/s | coverage `99.88%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_bb_4G` — 4G-Mobilfunkverfügbarkeit | coverage `99.88%` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.

Working logic of the thesis-candidate set:
- it keeps more nuance than the strict curated 17, especially within labour market strain and accessibility;
- it still removes the clearest duplicates and definition problems;
- it is broad enough for a meaningful PCA, but much easier to interpret than the old 51-variable block;
- it avoids keeping several variables that say almost the same thing at the same time.

Thesis-candidate correlation heatmap:
- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/thesis_candidate_correlation_heatmap.png`

Top correlation pairs from the thesis-candidate set:
- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_top_correlation_pairs.csv`

## 9. Student sensitivity block

Because student populations can be socio-economically precarious but are also strongly tied to urban university locations, student indicators are not included in the main thesis-candidate PCA.
- Instead, they are added only in a dedicated sensitivity run.
- Student sensitivity variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_variable_list.csv`
- Student sensitivity correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/student_sensitivity_correlation_heatmap.png`
- Student sensitivity top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_top_correlation_pairs.csv`
- Student sensitivity scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_scree.csv`

## 10. Proposed strict curated variable set for the main vulnerability PCA

These are the proposed core variables for the main interpretable PCA:
- `a_ALGII_SGBII` — ALG II-Leistungen an SGBII | coverage `98.80%` | Keep as core candidate in the curated vulnerability PCA.
- `a_aloLang` — Langzeitarbeitslose | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_alo_u25_einw` — Jüngere Arbeitslose | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_Minijobs` — Anteil Minijobs an den Beschäftigungsverhältnissen | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_kaufkraft` — Kaufkraft | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_einkst_bev` — Einkommensteuer | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_hheink_niedrig` — Haushalte mit niedrigem Einkommen | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_bev65um` — Einwohner 65 Jahre und älter | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_HH1` — Einpersonenhaushalte | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_hh_kind` — Haushalte mit Kindern | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_G02_SUP_DIST` — Entfernung zum Supermarkt/Discounter | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q01_APO_DIST` — Entfernung zur Apotheke | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q07_HA_DIST` — Entfernung zum Hausarzt | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `m_OEV20_DIST` — Entfernung zur ÖV Haltestelle | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_ärzte_bev` — Ärzte | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `q_abhg_alt` — Abhängigenquote Alte | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.
- `a_bb_100Mbits` — Bandbreitenverfügbarkeit mindestens 100 Mbit/s | coverage `99.88%` | Keep as core candidate in the curated vulnerability PCA.

Why this curated set is stronger:
- it keeps one representative for each main vulnerability mechanism instead of several near-duplicates;
- it avoids explicit definition conflicts like `a_BG1P` and the SGB-II housing-cost share problem;
- it reduces compositional over-weighting from nested age, income, physician, and broadband blocks;
- it stays close to the literature-based logic of deprivation, demographic sensitivity, household/social structure, and accessibility/adaptive capacity.

## 11. Selected variables that need explicit caution

- `a_ALGII_SGBII` — ALG II-Leistungen an SGBII | coverage `98.80%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_newfBGu15_bev` — Kinderarmut | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_aloLang` — Langzeitarbeitslose | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_alo_u25_einw` — Jüngere Arbeitslose | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_alo_ü55_einw` — Ältere Arbeitslose | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_Minijobs` — Anteil Minijobs an den Beschäftigungsverhältnissen | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_svw` — Beschäftigtenquote | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_kaufkraft` — Kaufkraft | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_einkst_bev` — Einkommensteuer | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_hheink_niedrig` — Haushalte mit niedrigem Einkommen | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_bev65um` — Einwohner 65 Jahre und älter | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_bev_0006` — Einwohner unter 6 Jahre | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_HH1` — Einpersonenhaushalte | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_hh_kind` — Haushalte mit Kindern | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_ewfBG_allein` — Alleinerziehende erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `m_G02_SUP_DIST` — Entfernung zum Supermarkt/Discounter | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q01_APO_DIST` — Entfernung zur Apotheke | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q07_HA_DIST` — Entfernung zum Hausarzt | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_OEV20_DIST` — Entfernung zur ÖV Haltestelle | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_P01_PRIM_DIST` — Entfernung zur Grundschule | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_ärzte_bev` — Ärzte | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_bb_50Mbits` — Bandbreitenverfügbarkeit mindestens 50 Mbit/s | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_bb_4G` — 4G-Mobilfunkverfügbarkeit | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_stud` — Studierende | coverage `99.88%` | status `student_sensitivity_only` | Use only in a student-sensitivity PCA. These variables mainly capture university-town structure rather than core municipality-level social vulnerability.
- `q_stud_1825` — Studierende je 100 Einwohner 18 bis 25 Jahre | coverage `99.88%` | status `student_sensitivity_only` | Use only in a student-sensitivity PCA. These variables mainly capture university-town structure rather than core municipality-level social vulnerability.
- `q_stud_fh` — Studierende an FH | coverage `99.88%` | status `student_sensitivity_only` | Use only in a student-sensitivity PCA. These variables mainly capture university-town structure rather than core municipality-level social vulnerability.
- `q_abhg_alt` — Abhängigenquote Alte | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_bb_100Mbits` — Bandbreitenverfügbarkeit mindestens 100 Mbit/s | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_Unterkunft_SGBII` — Leistungen für Unterkunft an SGBII | coverage `98.80%` | status `drop_overlap_with_algii` | Conceptually overlaps with a_ALGII_SGBII and is difficult to interpret as an independent vulnerability signal.
- `a_BG1P` — Einpersonen-Bedarfsgemeinschaften | coverage `98.44%` | status `drop_definition_mismatch` | Officially this is Einpersonen-Bedarfsgemeinschaften, not single-parent households; the old label was misleading.
- `a_BGKind` — Bedarfsgemeinschaften mit Kindern | coverage `98.44%` | status `review_sgbii_specific_household_measure` | Potentially informative, but highly specific to SGB II household composition and hard to generalize.
- `a_BG5um` — Große Bedarfsgemeinschaften | coverage `98.44%` | status `review_sgbii_specific_household_measure` | Potentially informative, but highly specific to SGB II household composition and hard to generalize.
- `a_bev75um` — Einwohner 75 Jahre und älter | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bb_1000Mbits` — Bandbreitenverfügbarkeit mindestens 1.000 Mbit/s | coverage `99.88%` | status `drop_duplicate_digital_measure` | Measures the same digital access construct at multiple thresholds; keep only one representative.
- `q_hausarzt_bev` — Hausärzte | coverage `99.88%` | status `drop_duplicate_health_measure` | Specialist-specific supply measure; keep total physician availability instead.
- `q_bev_fl` — Einwohnerdichte | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `q_bevsva_qkm` — Einwohner-Arbeitsplatz-Dichte | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.

## 12. Full corridor variable inventory

The complete review table is saved here:
- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/corridor_variable_review_table.csv`

Below, all corridor variables are listed by broad socio-economic dimension with a proposed handling note.

### Basic population and structural counts

- Indicators in corridor inventory: `9`
- Median coverage in this block: `99.88%`

- `TN23-kataster_qkm` — Bodenfläche gesamt qkm | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `alo` — Arbeitslose | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `bev_korr` — Bevölkerung (mit BBSR-Zensuskorrekturen) | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `ewf_1565_ges` — Erwerbsfähige Bevölkerung (15 bis unter 65 Jahre) | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `sva` — Sozialversicherungspflichtig Beschäftigte am Arbeitsort | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `svw` — Sozialversicherungspflichtig Beschäftigte am Wohnort | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `xbev` — Bevölkerung gesamt | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `xbevf` — Bevölkerung weiblich | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `xbevm` — Bevölkerung männlich | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.

### Demography and household structure

- Indicators in corridor inventory: `35`
- Median coverage in this block: `99.88%`

- `a_bev65um` — Einwohner 65 Jahre und älter | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_hh_kind` — Haushalte mit Kindern | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_HH1` — Einpersonenhaushalte | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_abhg_alt` — Abhängigenquote Alte | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_bev0003` — Einwohner unter 3 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev0306` — Einwohner von 3 bis unter 6 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev0618` — Einwohner von 6 bis unter 18 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev1825` — Einwohner von 18 bis unter 25 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev2530` — Einwohner von 25 bis unter 30 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev3050` — Einwohner von 30 bis unter 50 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev5065` — Einwohner von 50 bis unter 65 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev6575` — Einwohner von 65 bis unter 75 Jahren | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_bev75um` — Einwohner 75 Jahre und älter | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `i_saldo_nat` — Natürlicher Saldo | coverage `99.88%` | status `context_demographic_dynamics` | Demographic dynamics may be reported descriptively, but are not core vulnerability dimensions here.
- `i_wans` — Gesamtwanderungssaldo | coverage `99.88%` | status `context_demographic_dynamics` | Demographic dynamics may be reported descriptively, but are not core vulnerability dimensions here.
- `q_abhg_jung` — Abhängigenquote Junge | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `q_bev_fl` — Einwohnerdichte | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `q_bevsva_qkm` — Einwohner-Arbeitsplatz-Dichte | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.
- `a_bev1825_f` — Einwohner von 18 bis unter 25 Jahren, Frauen | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_bev2530_f` — Einwohner von 25 bis unter 30 Jahren, Frauen | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_bev6575_f` — Einwohner von 65 bis unter 75 Jahren, Frauen | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_bev65um_f` — Einwohner 65 Jahre und älter, Frauen | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_bev75um_f` — Einwohner 75 Jahre und älter, Frauen | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_bevMZ` — Bevölkerung in Mittelzentren | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_bevOZ` — Bevölkerung in Oberzentren | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_bev_0006` — Einwohner unter 6 Jahre | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_bevd150_` — Ländlichkeit | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_bevf` — Frauenanteil | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_bevf_2040` — Frauenanteil 20 bis unter 40 Jahre | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_geb_bev` — Geborene | coverage `99.88%` | status `context_demographic_dynamics` | Demographic dynamics may be reported descriptively, but are not core vulnerability dimensions here.
- `a_gest_bev` — Gestorbene | coverage `99.88%` | status `context_demographic_dynamics` | Demographic dynamics may be reported descriptively, but are not core vulnerability dimensions here.
- `e10_bev` — Bevölkerungsentwicklung (10 Jahre) | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `e5_bev` — Bevölkerungsentwicklung (5 Jahre) | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `m_bev_alter` — Durchschnittsalter der Bevölkerung | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_HH` — Haushaltsgröße | coverage `99.88%` | status `context_or_control_only` | Useful as descriptive context or control, but not as part of the main vulnerability index.

### Economic structure and tourism

- Indicators in corridor inventory: `3`
- Median coverage in this block: `99.88%`

- `m_übern` — Verweildauer in Beherbergungsbetrieben | coverage `44.67%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.
- `q_schlafg_bev` — Schlafgelegenheiten in Beherbergungsbetrieben | coverage `99.88%` | status `drop_tourism_context` | Not part of the social-vulnerability concept for the main index.
- `q_übern_bev` — Gästeübernachtungen in Beherbergungsbetrieben | coverage `99.88%` | status `drop_tourism_context` | Not part of the social-vulnerability concept for the main index.

### Education

- Indicators in corridor inventory: `7`
- Median coverage in this block: `3.95%`

- `q_stud` — Studierende | coverage `99.88%` | status `student_sensitivity_only` | Use only in a student-sensitivity PCA. These variables mainly capture university-town structure rather than core municipality-level social vulnerability.
- `q_stud_1825` — Studierende je 100 Einwohner 18 bis 25 Jahre | coverage `99.88%` | status `student_sensitivity_only` | Use only in a student-sensitivity PCA. These variables mainly capture university-town structure rather than core municipality-level social vulnerability.
- `q_stud_fh` — Studierende an FH | coverage `99.88%` | status `student_sensitivity_only` | Use only in a student-sensitivity PCA. These variables mainly capture university-town structure rather than core municipality-level social vulnerability.
- `a_stud_1` — Studierende im 1. Semester | coverage `3.95%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.
- `a_stud_a` — Ausländische Studierende | coverage `3.95%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.
- `a_stud_m` — Männliche Studierende | coverage `3.95%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.
- `a_stud_w` — Weibliche Studierende | coverage `3.95%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.

### Health, services and accessibility

- Indicators in corridor inventory: `33`
- Median coverage in this block: `99.88%`

- `a_bb_100Mbits` — Bandbreitenverfügbarkeit mindestens 100 Mbit/s | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_G02_SUP_DIST` — Entfernung zum Supermarkt/Discounter | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_OEV20_DIST` — Entfernung zur ÖV Haltestelle | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q01_APO_DIST` — Entfernung zur Apotheke | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `m_Q07_HA_DIST` — Entfernung zum Hausarzt | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_ärzte_bev` — Ärzte | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_bb_1000Mbits` — Bandbreitenverfügbarkeit mindestens 1.000 Mbit/s | coverage `99.88%` | status `drop_duplicate_digital_measure` | Measures the same digital access construct at multiple thresholds; keep only one representative.
- `a_bb_4G` — 4G-Mobilfunkverfügbarkeit | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_bb_50Mbits` — Bandbreitenverfügbarkeit mindestens 50 Mbit/s | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `m_P01_PRIM_DIST` — Entfernung zur Grundschule | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_allgemeinärzte_bev` — Allgemeinärzte | coverage `99.88%` | status `drop_duplicate_health_measure` | Specialist-specific supply measure; keep total physician availability instead.
- `q_hausarzt_bev` — Hausärzte | coverage `99.88%` | status `drop_duplicate_health_measure` | Specialist-specific supply measure; keep total physician availability instead.
- `q_internist_bev` — Internisten | coverage `99.88%` | status `drop_duplicate_health_measure` | Specialist-specific supply measure; keep total physician availability instead.
- `q_kinderarzt_kinder` — Kinderärzte | coverage `99.88%` | status `drop_duplicate_health_measure` | Specialist-specific supply measure; keep total physician availability instead.
- `a_G02_SUP_ANT` — Nahversorgung Supermarkt/Discounter | coverage `99.88%` | status `review_redundant_access_measure` | Access concept already captured by distance-based measures; use as optional robustness check only.
- `a_OEV20_ANT` — Nahversorgung ÖV Haltestelle | coverage `99.88%` | status `review_redundant_access_measure` | Access concept already captured by distance-based measures; use as optional robustness check only.
- `a_P01_PRIM_ANT` — Nahversorgung Grundschule | coverage `99.88%` | status `review_redundant_access_measure` | Access concept already captured by distance-based measures; use as optional robustness check only.
- `a_Q01_APO_ANT` — Nahversorgung Apotheke | coverage `99.88%` | status `review_redundant_access_measure` | Access concept already captured by distance-based measures; use as optional robustness check only.
- `a_Q07_HA_ANT` — Nahversorgung Hausarzt | coverage `99.88%` | status `review_redundant_access_measure` | Access concept already captured by distance-based measures; use as optional robustness check only.
- `a_ausp_svw` — Auspendler | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_einp_sva` — Einpendler | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_pend150` — Pendler mit Arbeitsweg 150 km und mehr | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_pend300` — Pendler mit Arbeitsweg 300 km und mehr | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_pend50` — Pendler mit Arbeitsweg 50 km und mehr | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `i_fz_mz` — Erreichbarkeit von Mittelzentren | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `i_fz_oz` — Erreichbarkeit von Oberzentren | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `m_ErrAir_bev` — Erreichbarkeit von Flughäfen | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `m_ErrBAB_bev` — Erreichbarkeit von Autobahnen | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `m_ErrIC_bev` — Erreichbarkeit von IC/EC/ICE-Bahnhöfen | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_pendlersaldo` — Pendlersaldo | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_unf_bev` — Straßenverkehrsunfälle | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_vunp_bev` — Verunglückte im Straßenverkehr | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_vunpt_bev` — Getötete im Straßenverkehr | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.

### Housing and land-use context

- Indicators in corridor inventory: `27`
- Median coverage in this block: `99.88%`

- `a_erh` — Erholungsfläche | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_fert_wg12` — Neue Ein- und Zweifamilienhäuser | coverage `99.04%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_fert_wo12` — Neubauwohnungen in Ein- und Zweifamilienhäusern | coverage `97.96%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_fert_wohn` — Fertiggestellte Wohnungen je Wohnung im Bestand | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_freifläche` — Freifläche | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_gen_wo12` — Baugenehmigungen für Wohnungen in Ein- und Zweifamilienhäusern | coverage `99.16%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_gen_wo3um` — Baugenehmigungen für Wohnungen in Mehrfamilienhäusern | coverage `99.16%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_landwirtschaft` — Landwirtschaftsfläche | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_naturnah` — Naturnähere Fläche | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_suv_fl` — Siedlungs- und Verkehrsfläche | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wald` — Waldfläche | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wasser` — Wasserfläche | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wg_wo12` — Ein- und Zweifamilienhäuser | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wg_wo3um` — Mehrfamilienhäuser | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wo_r12` — Ein- und Zweiraumwohnungen | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wo_r5um` — 5- und mehr Raum-Wohnungen | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wo_wg12` — Wohnungen in Ein- und Zweifamilienhäusern | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `a_wo_wg3um` — Wohnungen in Mehrfamilienhäusern | coverage `99.88%` | status `context_land_use` | Land-use and housing context variable, not a direct social-vulnerability indicator.
- `q_erhfl_bev` — Erholungsfläche je Einwohner | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_ew_suv_qkm` — Siedlungsdichte in km² | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_fert_wo12_bev` — Neubauwohnungen in Ein- und Zweifamilienhäusern je Einwohner | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_fert_wo3um_bev` — Neubauwohnungen in Mehrfamilienhäusern je Einwohner | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_fert_wo_bev` — Neubauwohnungen je Einwohner | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_freifläche_bev` — Freifläche je Einwohner | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_gen_wo_ew` — Baugenehmigungen für Wohnungenje Einwohner | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_naturnah_bev` — Naturnähere Fläche je Einwohner | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `v_suvfl_vorjahr` — Flächenneuinanspruchnahme | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.

### Income, deprivation and social transfers

- Indicators in corridor inventory: `25`
- Median coverage in this block: `98.80%`

- `a_ALGII_SGBII` — ALG II-Leistungen an SGBII | coverage `98.80%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_hheink_niedrig` — Haushalte mit niedrigem Einkommen | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_einkst_bev` — Einkommensteuer | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_kaufkraft` — Kaufkraft | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_BG1P` — Einpersonen-Bedarfsgemeinschaften | coverage `98.44%` | status `drop_definition_mismatch` | Officially this is Einpersonen-Bedarfsgemeinschaften, not single-parent households; the old label was misleading.
- `a_BG5um` — Große Bedarfsgemeinschaften | coverage `98.44%` | status `review_sgbii_specific_household_measure` | Potentially informative, but highly specific to SGB II household composition and hard to generalize.
- `a_BGKind` — Bedarfsgemeinschaften mit Kindern | coverage `98.44%` | status `review_sgbii_specific_household_measure` | Potentially informative, but highly specific to SGB II household composition and hard to generalize.
- `a_Unterkunft_SGBII` — Leistungen für Unterkunft an SGBII | coverage `98.80%` | status `drop_overlap_with_algii` | Conceptually overlaps with a_ALGII_SGBII and is difficult to interpret as an independent vulnerability signal.
- `a_hheink_hoch` — Haushalte mit hohem Einkommen | coverage `99.88%` | status `drop_complement_or_context` | Either compositional counterpart or broader fiscal context rather than a distinct vulnerability mechanism.
- `a_hheink_mittel` — Haushalte mit mittlerem Einkommen | coverage `99.88%` | status `drop_complement_or_context` | Either compositional counterpart or broader fiscal context rather than a distinct vulnerability mechanism.
- `d_steuereinnahme` — Steuereinnahmen | coverage `99.88%` | status `drop_complement_or_context` | Either compositional counterpart or broader fiscal context rather than a distinct vulnerability mechanism.
- `q_gewst_bev` — Gewerbesteuer | coverage `99.88%` | status `drop_complement_or_context` | Either compositional counterpart or broader fiscal context rather than a distinct vulnerability mechanism.
- `a_ewfBG` — Erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `review_sgbii_specific_household_measure` | Potentially informative, but highly specific to SGB II household composition and hard to generalize.
- `a_ewfBG_55um_` — Ältere erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_ewfBG_allein` — Alleinerziehende erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `a_ewfBG_f` — Erwerbsfähige Leistungsberechtigte (Frauen) | coverage `98.56%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_ewfBG_u25_` — Junge erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `q_PBG_bev` — Personen in Bedarfsgemeinschaften | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_einzelhandelskaufkraft` — Einzelhandelsrelevante Kaufkraft | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_gest_bev` — Steuerkraft | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_investZ` — Zuweisungen für Investitionsfördermaßnahmen | coverage `29.70%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.
- `q_newfBGu15_bev` — Kinderarmut | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_sach` — Ausgaben für Sachinvestitionen | coverage `30.06%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.
- `q_schlüsselzuw` — Schlüsselzuweisungen | coverage `29.34%` | status `drop_low_coverage` | Too sparse in the corridor sample for a stable main PCA.
- `q_umsst` — Umsatzsteuer | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.

### Labour market and employment

- Indicators in corridor inventory: `37`
- Median coverage in this block: `99.88%`

- `a_aloLang` — Langzeitarbeitslose | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_alo_u25_einw` — Jüngere Arbeitslose | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `a_Minijobs` — Anteil Minijobs an den Beschäftigungsverhältnissen | coverage `99.88%` | status `core_curated_pca` | Keep as core candidate in the curated vulnerability PCA.
- `q_alo_ü55_einw` — Ältere Arbeitslose | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `r_ewf_jungalt` — Verhältnis junge zu alte Erwerbsfähige | coverage `99.88%` | status `drop_nested_age_block` | Nested age-composition block; retaining many of these together overweights the same demographic structure.
- `a_aloLang_f` — Weibliche Langzeitarbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_aloLang_m` — Männliche Langzeitarbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_ausländer` — Ausländische Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_ausländer_f` — Ausländische weibliche Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_ausländer_m` — Ausländische männliche Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_experte` — Arbeitslose mit Anforderungsniveau Experte | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_alo_f` — Arbeitslose Frauen | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_fachkraft` — Arbeitslose mit Anforderungsniveau Fachkraft | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_alo_helfer` — Arbeitslose mit Anforderungsniveau Helfer | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_alo_m` — Arbeitslose Männer | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_oAusb` — Arbeitslose ohne Ausbildung | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_alo_spezialist` — Arbeitslose mit Anforderungsniveau Spezialist | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_alo_u25` — Anteil jüngere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_u25_f` — Anteil weibliche jüngere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_u25_m` — Anteil männliche jüngere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_ü55` — Anteil ältere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_ü55_f` — Anteil weibliche ältere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_alo_ü55_m` — Anteil männliche ältere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_gb_Frauen` — Anteil Minijobs (Frauen) | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_gb_Männer` — Anteil Minijobs (Männer) | coverage `99.88%` | status `drop_redundant_subgroup` | Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.
- `a_gb_knj` — Anteil Minijobs (ausschließlich) an geringfügig Beschäftigten | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_gb_knj_prBv` — Anteil Minijobs (ausschließlich) an den Beschäftigungsverhältnissen | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_gb_nj` — Anteil Minijobs (Nebenverdienst) | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `a_gb_nj_prBv` — Anteil Minijobs (Nebenverdienst) an den Beschäftigungsverhältnissen | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_alo_u25_einw_f` — Weibliche jüngere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_alo_u25_einw_m` — Männliche jüngere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_alo_ü55_einw_f` — Weibliche ältere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_alo_ü55_einw_m` — Männliche ältere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_svw` — Beschäftigtenquote | coverage `99.88%` | status `core_thesis_candidate` | Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.
- `q_svw_f` — Beschäftigtenquote Frauen | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_svw_m` — Beschäftigtenquote Männer | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.
- `q_svw_ü55` — Quote ältere Beschäftigte | coverage `99.88%` | status `review_not_in_curated_set` | Review individually; not obviously required for the core vulnerability index.

## 13. Recommended next step

- Keep the original 51-variable PCA as an exploratory comparison / appendix result.
- Use the all-176 PCA only as a diagnostic exercise, not as the main thesis index.
- Use the thesis-candidate set as the current working set for substantive thesis analyses.
- Use the student-sensitivity run only as a robustness check for how strongly university-center structure changes the PCA.
- Keep the strict curated 17-variable set as a robustness / interpretability check.
- If needed, we can next turn this directly into a keep/review/drop decision meeting note or a reply draft to the supervisors.
