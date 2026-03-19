# INKAR HQ500 OVERVIEW

Created: `2026-03-19 14:55:55`

## 1. Study area

The study area is defined as the HQ500 Elbe flood corridor.
- Municipalities were retained if they intersect the RP500 flood extent derived from the prepared EFAS flood rasters.
- This replaces the older basin-wide definition and focuses the analysis on municipalities that are actually flood-relevant under an extreme return period.
- Corridor municipalities in the current geometry layer: `835`.
- Municipalities with a latest available raw-INKAR record after corridor filtering: `834`.
- One municipality (`16076094`, Berga-Wünschendorf) remains missing in the current raw-INKAR latest-build table and should be documented explicitly.

Core corridor files:
- Geometry + all original INKAR codes: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/gpkg/corridor_full_inkar_original_codes.gpkg`
- Flat table + all original INKAR codes: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_full_inkar_original_codes.csv`

## 2. All available INKAR variables in the corridor

- Number of original INKAR indicators available in the corridor latest-build table: `176`.
- Inventory with official names, coverage, and broad dimensions: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_inkar_indicator_inventory_dimensioned.csv`
- Short summary by broad dimension: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_inkar_dimension_summary.csv`

All available corridor indicators grouped by broad dimension:

### Basic population and structural counts

- Indicators in this block: `9`
- Median coverage: `99.88%`

- `TN23-kataster_qkm` — Bodenfläche gesamt qkm | coverage `99.88%` | status `context_or_control_only`
- `alo` — Arbeitslose | coverage `99.88%` | status `context_or_control_only`
- `bev_korr` — Bevölkerung (mit BBSR-Zensuskorrekturen) | coverage `99.88%` | status `context_or_control_only`
- `ewf_1565_ges` — Erwerbsfähige Bevölkerung (15 bis unter 65 Jahre) | coverage `99.88%` | status `context_or_control_only`
- `sva` — Sozialversicherungspflichtig Beschäftigte am Arbeitsort | coverage `99.88%` | status `context_or_control_only`
- `svw` — Sozialversicherungspflichtig Beschäftigte am Wohnort | coverage `99.88%` | status `context_or_control_only`
- `xbev` — Bevölkerung gesamt | coverage `99.88%` | status `context_or_control_only`
- `xbevf` — Bevölkerung weiblich | coverage `99.88%` | status `context_or_control_only`
- `xbevm` — Bevölkerung männlich | coverage `99.88%` | status `context_or_control_only`

### Demography and household structure

- Indicators in this block: `35`
- Median coverage: `99.88%`

- `a_bev65um` — Einwohner 65 Jahre und älter | coverage `99.88%` | status `core_curated_pca`
- `a_hh_kind` — Haushalte mit Kindern | coverage `99.88%` | status `core_curated_pca`
- `q_HH1` — Einpersonenhaushalte | coverage `99.88%` | status `core_curated_pca`
- `q_abhg_alt` — Abhängigenquote Alte | coverage `99.88%` | status `core_curated_pca`
- `a_bev0003` — Einwohner unter 3 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev0306` — Einwohner von 3 bis unter 6 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev0618` — Einwohner von 6 bis unter 18 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev1825` — Einwohner von 18 bis unter 25 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev2530` — Einwohner von 25 bis unter 30 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev3050` — Einwohner von 30 bis unter 50 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev5065` — Einwohner von 50 bis unter 65 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev6575` — Einwohner von 65 bis unter 75 Jahren | coverage `99.88%` | status `drop_nested_age_block`
- `a_bev75um` — Einwohner 75 Jahre und älter | coverage `99.88%` | status `drop_nested_age_block`
- `i_saldo_nat` — Natürlicher Saldo | coverage `99.88%` | status `context_demographic_dynamics`
- `i_wans` — Gesamtwanderungssaldo | coverage `99.88%` | status `context_demographic_dynamics`
- `q_abhg_jung` — Abhängigenquote Junge | coverage `99.88%` | status `drop_nested_age_block`
- `q_bev_fl` — Einwohnerdichte | coverage `99.88%` | status `context_or_control_only`
- `q_bevsva_qkm` — Einwohner-Arbeitsplatz-Dichte | coverage `99.88%` | status `context_or_control_only`
- `a_bev1825_f` — Einwohner von 18 bis unter 25 Jahren, Frauen | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_bev2530_f` — Einwohner von 25 bis unter 30 Jahren, Frauen | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_bev6575_f` — Einwohner von 65 bis unter 75 Jahren, Frauen | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_bev65um_f` — Einwohner 65 Jahre und älter, Frauen | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_bev75um_f` — Einwohner 75 Jahre und älter, Frauen | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_bevMZ` — Bevölkerung in Mittelzentren | coverage `99.88%` | status `review_not_in_curated_set`
- `a_bevOZ` — Bevölkerung in Oberzentren | coverage `99.88%` | status `review_not_in_curated_set`
- `a_bev_0006` — Einwohner unter 6 Jahre | coverage `99.88%` | status `core_thesis_candidate`
- `a_bevd150_` — Ländlichkeit | coverage `99.88%` | status `review_not_in_curated_set`
- `a_bevf` — Frauenanteil | coverage `99.88%` | status `review_not_in_curated_set`
- `a_bevf_2040` — Frauenanteil 20 bis unter 40 Jahre | coverage `99.88%` | status `review_not_in_curated_set`
- `a_geb_bev` — Geborene | coverage `99.88%` | status `context_demographic_dynamics`
- `a_gest_bev` — Gestorbene | coverage `99.88%` | status `context_demographic_dynamics`
- `e10_bev` — Bevölkerungsentwicklung (10 Jahre) | coverage `99.88%` | status `review_not_in_curated_set`
- `e5_bev` — Bevölkerungsentwicklung (5 Jahre) | coverage `99.88%` | status `review_not_in_curated_set`
- `m_bev_alter` — Durchschnittsalter der Bevölkerung | coverage `99.88%` | status `review_not_in_curated_set`
- `q_HH` — Haushaltsgröße | coverage `99.88%` | status `context_or_control_only`

### Economic structure and tourism

- Indicators in this block: `3`
- Median coverage: `99.88%`

- `m_übern` — Verweildauer in Beherbergungsbetrieben | coverage `44.67%` | status `drop_low_coverage`
- `q_schlafg_bev` — Schlafgelegenheiten in Beherbergungsbetrieben | coverage `99.88%` | status `drop_tourism_context`
- `q_übern_bev` — Gästeübernachtungen in Beherbergungsbetrieben | coverage `99.88%` | status `drop_tourism_context`

### Education

- Indicators in this block: `7`
- Median coverage: `3.95%`

- `q_stud` — Studierende | coverage `99.88%` | status `student_sensitivity_only`
- `q_stud_1825` — Studierende je 100 Einwohner 18 bis 25 Jahre | coverage `99.88%` | status `student_sensitivity_only`
- `q_stud_fh` — Studierende an FH | coverage `99.88%` | status `student_sensitivity_only`
- `a_stud_1` — Studierende im 1. Semester | coverage `3.95%` | status `drop_low_coverage`
- `a_stud_a` — Ausländische Studierende | coverage `3.95%` | status `drop_low_coverage`
- `a_stud_m` — Männliche Studierende | coverage `3.95%` | status `drop_low_coverage`
- `a_stud_w` — Weibliche Studierende | coverage `3.95%` | status `drop_low_coverage`

### Health, services and accessibility

- Indicators in this block: `33`
- Median coverage: `99.88%`

- `a_bb_100Mbits` — Bandbreitenverfügbarkeit mindestens 100 Mbit/s | coverage `99.88%` | status `core_curated_pca`
- `m_G02_SUP_DIST` — Entfernung zum Supermarkt/Discounter | coverage `99.88%` | status `core_curated_pca`
- `m_OEV20_DIST` — Entfernung zur ÖV Haltestelle | coverage `99.88%` | status `core_curated_pca`
- `m_Q01_APO_DIST` — Entfernung zur Apotheke | coverage `99.88%` | status `core_curated_pca`
- `m_Q07_HA_DIST` — Entfernung zum Hausarzt | coverage `99.88%` | status `core_curated_pca`
- `q_ärzte_bev` — Ärzte | coverage `99.88%` | status `core_curated_pca`
- `a_bb_1000Mbits` — Bandbreitenverfügbarkeit mindestens 1.000 Mbit/s | coverage `99.88%` | status `drop_duplicate_digital_measure`
- `a_bb_4G` — 4G-Mobilfunkverfügbarkeit | coverage `99.88%` | status `core_thesis_candidate`
- `a_bb_50Mbits` — Bandbreitenverfügbarkeit mindestens 50 Mbit/s | coverage `99.88%` | status `core_thesis_candidate`
- `m_P01_PRIM_DIST` — Entfernung zur Grundschule | coverage `99.88%` | status `core_thesis_candidate`
- `q_allgemeinärzte_bev` — Allgemeinärzte | coverage `99.88%` | status `drop_duplicate_health_measure`
- `q_hausarzt_bev` — Hausärzte | coverage `99.88%` | status `drop_duplicate_health_measure`
- `q_internist_bev` — Internisten | coverage `99.88%` | status `drop_duplicate_health_measure`
- `q_kinderarzt_kinder` — Kinderärzte | coverage `99.88%` | status `drop_duplicate_health_measure`
- `a_G02_SUP_ANT` — Nahversorgung Supermarkt/Discounter | coverage `99.88%` | status `review_redundant_access_measure`
- `a_OEV20_ANT` — Nahversorgung ÖV Haltestelle | coverage `99.88%` | status `review_redundant_access_measure`
- `a_P01_PRIM_ANT` — Nahversorgung Grundschule | coverage `99.88%` | status `review_redundant_access_measure`
- `a_Q01_APO_ANT` — Nahversorgung Apotheke | coverage `99.88%` | status `review_redundant_access_measure`
- `a_Q07_HA_ANT` — Nahversorgung Hausarzt | coverage `99.88%` | status `review_redundant_access_measure`
- `a_ausp_svw` — Auspendler | coverage `99.88%` | status `review_not_in_curated_set`
- `a_einp_sva` — Einpendler | coverage `99.88%` | status `review_not_in_curated_set`
- `a_pend150` — Pendler mit Arbeitsweg 150 km und mehr | coverage `99.88%` | status `review_not_in_curated_set`
- `a_pend300` — Pendler mit Arbeitsweg 300 km und mehr | coverage `99.88%` | status `review_not_in_curated_set`
- `a_pend50` — Pendler mit Arbeitsweg 50 km und mehr | coverage `99.88%` | status `review_not_in_curated_set`
- `i_fz_mz` — Erreichbarkeit von Mittelzentren | coverage `99.88%` | status `review_not_in_curated_set`
- `i_fz_oz` — Erreichbarkeit von Oberzentren | coverage `99.88%` | status `review_not_in_curated_set`
- `m_ErrAir_bev` — Erreichbarkeit von Flughäfen | coverage `99.88%` | status `review_not_in_curated_set`
- `m_ErrBAB_bev` — Erreichbarkeit von Autobahnen | coverage `99.88%` | status `review_not_in_curated_set`
- `m_ErrIC_bev` — Erreichbarkeit von IC/EC/ICE-Bahnhöfen | coverage `99.88%` | status `review_not_in_curated_set`
- `q_pendlersaldo` — Pendlersaldo | coverage `99.88%` | status `review_not_in_curated_set`
- `q_unf_bev` — Straßenverkehrsunfälle | coverage `99.88%` | status `review_not_in_curated_set`
- `q_vunp_bev` — Verunglückte im Straßenverkehr | coverage `99.88%` | status `review_not_in_curated_set`
- `q_vunpt_bev` — Getötete im Straßenverkehr | coverage `99.88%` | status `review_not_in_curated_set`

### Housing and land-use context

- Indicators in this block: `27`
- Median coverage: `99.88%`

- `a_erh` — Erholungsfläche | coverage `99.88%` | status `review_not_in_curated_set`
- `a_fert_wg12` — Neue Ein- und Zweifamilienhäuser | coverage `99.04%` | status `context_land_use`
- `a_fert_wo12` — Neubauwohnungen in Ein- und Zweifamilienhäusern | coverage `97.96%` | status `context_land_use`
- `a_fert_wohn` — Fertiggestellte Wohnungen je Wohnung im Bestand | coverage `99.88%` | status `context_land_use`
- `a_freifläche` — Freifläche | coverage `99.88%` | status `context_land_use`
- `a_gen_wo12` — Baugenehmigungen für Wohnungen in Ein- und Zweifamilienhäusern | coverage `99.16%` | status `context_land_use`
- `a_gen_wo3um` — Baugenehmigungen für Wohnungen in Mehrfamilienhäusern | coverage `99.16%` | status `context_land_use`
- `a_landwirtschaft` — Landwirtschaftsfläche | coverage `99.88%` | status `context_land_use`
- `a_naturnah` — Naturnähere Fläche | coverage `99.88%` | status `context_land_use`
- `a_suv_fl` — Siedlungs- und Verkehrsfläche | coverage `99.88%` | status `context_land_use`
- `a_wald` — Waldfläche | coverage `99.88%` | status `context_land_use`
- `a_wasser` — Wasserfläche | coverage `99.88%` | status `context_land_use`
- `a_wg_wo12` — Ein- und Zweifamilienhäuser | coverage `99.88%` | status `context_land_use`
- `a_wg_wo3um` — Mehrfamilienhäuser | coverage `99.88%` | status `context_land_use`
- `a_wo_r12` — Ein- und Zweiraumwohnungen | coverage `99.88%` | status `context_land_use`
- `a_wo_r5um` — 5- und mehr Raum-Wohnungen | coverage `99.88%` | status `context_land_use`
- `a_wo_wg12` — Wohnungen in Ein- und Zweifamilienhäusern | coverage `99.88%` | status `context_land_use`
- `a_wo_wg3um` — Wohnungen in Mehrfamilienhäusern | coverage `99.88%` | status `context_land_use`
- `q_erhfl_bev` — Erholungsfläche je Einwohner | coverage `99.88%` | status `review_not_in_curated_set`
- `q_ew_suv_qkm` — Siedlungsdichte in km² | coverage `99.88%` | status `review_not_in_curated_set`
- `q_fert_wo12_bev` — Neubauwohnungen in Ein- und Zweifamilienhäusern je Einwohner | coverage `99.88%` | status `review_not_in_curated_set`
- `q_fert_wo3um_bev` — Neubauwohnungen in Mehrfamilienhäusern je Einwohner | coverage `99.88%` | status `review_not_in_curated_set`
- `q_fert_wo_bev` — Neubauwohnungen je Einwohner | coverage `99.88%` | status `review_not_in_curated_set`
- `q_freifläche_bev` — Freifläche je Einwohner | coverage `99.88%` | status `review_not_in_curated_set`
- `q_gen_wo_ew` — Baugenehmigungen für Wohnungenje Einwohner | coverage `99.88%` | status `review_not_in_curated_set`
- `q_naturnah_bev` — Naturnähere Fläche je Einwohner | coverage `99.88%` | status `review_not_in_curated_set`
- `v_suvfl_vorjahr` — Flächenneuinanspruchnahme | coverage `99.88%` | status `review_not_in_curated_set`

### Income, deprivation and social transfers

- Indicators in this block: `25`
- Median coverage: `98.80%`

- `a_ALGII_SGBII` — ALG II-Leistungen an SGBII | coverage `98.80%` | status `core_curated_pca`
- `a_hheink_niedrig` — Haushalte mit niedrigem Einkommen | coverage `99.88%` | status `core_curated_pca`
- `q_einkst_bev` — Einkommensteuer | coverage `99.88%` | status `core_curated_pca`
- `q_kaufkraft` — Kaufkraft | coverage `99.88%` | status `core_curated_pca`
- `a_BG1P` — Einpersonen-Bedarfsgemeinschaften | coverage `98.44%` | status `drop_definition_mismatch`
- `a_BG5um` — Große Bedarfsgemeinschaften | coverage `98.44%` | status `review_sgbii_specific_household_measure`
- `a_BGKind` — Bedarfsgemeinschaften mit Kindern | coverage `98.44%` | status `review_sgbii_specific_household_measure`
- `a_Unterkunft_SGBII` — Leistungen für Unterkunft an SGBII | coverage `98.80%` | status `drop_overlap_with_algii`
- `a_hheink_hoch` — Haushalte mit hohem Einkommen | coverage `99.88%` | status `drop_complement_or_context`
- `a_hheink_mittel` — Haushalte mit mittlerem Einkommen | coverage `99.88%` | status `drop_complement_or_context`
- `d_steuereinnahme` — Steuereinnahmen | coverage `99.88%` | status `drop_complement_or_context`
- `q_gewst_bev` — Gewerbesteuer | coverage `99.88%` | status `drop_complement_or_context`
- `a_ewfBG` — Erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `review_sgbii_specific_household_measure`
- `a_ewfBG_55um_` — Ältere erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `drop_redundant_subgroup`
- `a_ewfBG_allein` — Alleinerziehende erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `core_thesis_candidate`
- `a_ewfBG_f` — Erwerbsfähige Leistungsberechtigte (Frauen) | coverage `98.56%` | status `drop_redundant_subgroup`
- `a_ewfBG_u25_` — Junge erwerbsfähige Leistungsberechtigte | coverage `98.56%` | status `drop_redundant_subgroup`
- `q_PBG_bev` — Personen in Bedarfsgemeinschaften | coverage `99.88%` | status `review_not_in_curated_set`
- `q_einzelhandelskaufkraft` — Einzelhandelsrelevante Kaufkraft | coverage `99.88%` | status `review_not_in_curated_set`
- `q_gest_bev` — Steuerkraft | coverage `99.88%` | status `review_not_in_curated_set`
- `q_investZ` — Zuweisungen für Investitionsfördermaßnahmen | coverage `29.70%` | status `drop_low_coverage`
- `q_newfBGu15_bev` — Kinderarmut | coverage `99.88%` | status `core_thesis_candidate`
- `q_sach` — Ausgaben für Sachinvestitionen | coverage `30.06%` | status `drop_low_coverage`
- `q_schlüsselzuw` — Schlüsselzuweisungen | coverage `29.34%` | status `drop_low_coverage`
- `q_umsst` — Umsatzsteuer | coverage `99.88%` | status `review_not_in_curated_set`

### Labour market and employment

- Indicators in this block: `37`
- Median coverage: `99.88%`

- `a_aloLang` — Langzeitarbeitslose | coverage `99.88%` | status `core_curated_pca`
- `q_alo_u25_einw` — Jüngere Arbeitslose | coverage `99.88%` | status `core_curated_pca`
- `a_Minijobs` — Anteil Minijobs an den Beschäftigungsverhältnissen | coverage `99.88%` | status `core_curated_pca`
- `q_alo_ü55_einw` — Ältere Arbeitslose | coverage `99.88%` | status `core_thesis_candidate`
- `r_ewf_jungalt` — Verhältnis junge zu alte Erwerbsfähige | coverage `99.88%` | status `drop_nested_age_block`
- `a_aloLang_f` — Weibliche Langzeitarbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_aloLang_m` — Männliche Langzeitarbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_ausländer` — Ausländische Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_ausländer_f` — Ausländische weibliche Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_ausländer_m` — Ausländische männliche Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_experte` — Arbeitslose mit Anforderungsniveau Experte | coverage `99.88%` | status `review_not_in_curated_set`
- `a_alo_f` — Arbeitslose Frauen | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_fachkraft` — Arbeitslose mit Anforderungsniveau Fachkraft | coverage `99.88%` | status `review_not_in_curated_set`
- `a_alo_helfer` — Arbeitslose mit Anforderungsniveau Helfer | coverage `99.88%` | status `review_not_in_curated_set`
- `a_alo_m` — Arbeitslose Männer | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_oAusb` — Arbeitslose ohne Ausbildung | coverage `99.88%` | status `review_not_in_curated_set`
- `a_alo_spezialist` — Arbeitslose mit Anforderungsniveau Spezialist | coverage `99.88%` | status `review_not_in_curated_set`
- `a_alo_u25` — Anteil jüngere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_u25_f` — Anteil weibliche jüngere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_u25_m` — Anteil männliche jüngere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_ü55` — Anteil ältere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_ü55_f` — Anteil weibliche ältere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_alo_ü55_m` — Anteil männliche ältere Arbeitslose | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_gb_Frauen` — Anteil Minijobs (Frauen) | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_gb_Männer` — Anteil Minijobs (Männer) | coverage `99.88%` | status `drop_redundant_subgroup`
- `a_gb_knj` — Anteil Minijobs (ausschließlich) an geringfügig Beschäftigten | coverage `99.88%` | status `review_not_in_curated_set`
- `a_gb_knj_prBv` — Anteil Minijobs (ausschließlich) an den Beschäftigungsverhältnissen | coverage `99.88%` | status `review_not_in_curated_set`
- `a_gb_nj` — Anteil Minijobs (Nebenverdienst) | coverage `99.88%` | status `review_not_in_curated_set`
- `a_gb_nj_prBv` — Anteil Minijobs (Nebenverdienst) an den Beschäftigungsverhältnissen | coverage `99.88%` | status `review_not_in_curated_set`
- `q_alo_u25_einw_f` — Weibliche jüngere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set`
- `q_alo_u25_einw_m` — Männliche jüngere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set`
- `q_alo_ü55_einw_f` — Weibliche ältere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set`
- `q_alo_ü55_einw_m` — Männliche ältere Arbeitslose | coverage `99.88%` | status `review_not_in_curated_set`
- `q_svw` — Beschäftigtenquote | coverage `99.88%` | status `core_thesis_candidate`
- `q_svw_f` — Beschäftigtenquote Frauen | coverage `99.88%` | status `review_not_in_curated_set`
- `q_svw_m` — Beschäftigtenquote Männer | coverage `99.88%` | status `review_not_in_curated_set`
- `q_svw_ü55` — Quote ältere Beschäftigte | coverage `99.88%` | status `review_not_in_curated_set`

## 3. Correlation matrix of all available corridor variables

The full 176-variable correlation matrix is useful as a diagnostic step, not as a final model input on its own.
- Correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/all176_correlation_heatmap.png`
- Top absolute correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/all176_top_correlation_pairs.csv`

Why this step matters:
- it shows where several indicators measure the same underlying construct;
- it makes compositional counterparts visible;
- it helps identify definition-sensitive indicators before PCA loadings are interpreted.

## 4. Original 51-variable wide PCA set

The original 51-variable set came from the earlier INKAR preprocessing script and was designed as a broad exploratory screening block.
- Source logic: one manually assembled selection spanning poverty/welfare, unemployment, demography, dependency, household structure, income groups, health care, digital infrastructure, accessibility, and density.
- This was a reasonable exploratory starting point because it allowed the covariance structure of a broad vulnerability field to emerge from the data.
- At the same time, the block was not yet fully cleaned for duplicated constructs or definition conflicts.

Original 51-variable files:
- Variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/original_wide51_variable_list.csv`
- Correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/old51_correlation_heatmap.png`
- Top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/old51_top_correlation_pairs.csv`
- Scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/old51_scree.csv`

Main methodological issue with the old 51-variable set:
- several variables partly say the same thing;
- nested age blocks over-weight demographic composition;
- multiple broadband thresholds duplicate one digital-access construct;
- specialist-specific physician indicators overlap strongly with total physician supply;
- SGB-II-related indicators include at least one definition-sensitive pair and one naming mismatch (`a_BG1P`).

## 5. Working thesis candidate set from all available variables

This is the recommended working set for the thesis at the current stage: broader than the strict 17-variable core, but cleaned enough to avoid the clearest double-counting problems.

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

Selection logic of the thesis candidate set:
- keep a broad but interpretable coverage of deprivation, labour-market strain, demographic sensitivity, household/social structure, accessibility, health access, and digital infrastructure;
- remove direct duplicates or near-duplicates wherever possible;
- remove variables with clearly problematic definitions for interpretation;
- keep the set broad enough for PCA while avoiding obvious over-weighting of the same latent construct.

Thesis-candidate files:
- Variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_variable_list.csv`
- Correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/thesis_candidate_correlation_heatmap.png`
- Top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_top_correlation_pairs.csv`
- Scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_scree.csv`
- Top loadings: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_top_loadings.csv`

Student variables are intentionally excluded from this main set and are handled only in a separate sensitivity PCA.

## 6. Student sensitivity block

- Student sensitivity variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_variable_list.csv`
- Student sensitivity correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/student_sensitivity_correlation_heatmap.png`
- Student sensitivity top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_top_correlation_pairs.csv`
- Student sensitivity scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_scree.csv`

Rationale:
- student populations can be economically precarious;
- but at municipality scale, student indicators also capture university-center structure very strongly;
- therefore they are better treated as a sensitivity block than as part of the main thesis index.

## 7. PCA comparison across all sets

- `original_51`:  51 variables, 834 municipalities, PC1 variance `16.18%`, cumulative PC1-PC4 `43.40%`, cumulative PC1-PC8 `60.07%`
- `all_176`: 176 variables, 834 municipalities, PC1 variance `15.14%`, cumulative PC1-PC4 `33.01%`, cumulative PC1-PC8 `44.14%`
- `thesis_candidate_23`:  23 variables, 834 municipalities, PC1 variance `24.21%`, cumulative PC1-PC4 `54.81%`, cumulative PC1-PC8 `72.60%`
- `student_sensitivity_26`:  26 variables, 834 municipalities, PC1 variance `22.14%`, cumulative PC1-PC4 `54.59%`, cumulative PC1-PC8 `72.14%`
- `curated_17`:  17 variables, 834 municipalities, PC1 variance `27.65%`, cumulative PC1-PC4 `63.44%`, cumulative PC1-PC8 `83.15%`

Interpretation:
- `all 176` is too broad to be the main thesis PCA, but very useful as a structural diagnostic.
- `original 51` remains useful as an exploratory benchmark and historical reference.
- `thesis candidate 23` is currently the most balanced working set for the thesis.
- `student sensitivity 26` shows what changes once student concentration variables are allowed into the PCA.
- `curated 17` remains the stricter fallback / robustness set.

## 8. Immediate methodological takeaway

- Yes, the old wide PCA made sense as an exploratory first step.
- No, it should probably not remain the only final vulnerability index without additional variable curation.
- The safer thesis strategy is to document three layers clearly:
  - `all 176` for inventory and structural diagnosis,
  - `original 51` for exploratory comparison,
  - `thesis candidate 23` as the main working PCA set,
  - `student sensitivity 26` as an additional robustness run.

## 9. Related note

A more detailed review note with additional handling commentary is saved here:
- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/INKAR_VARIABLE_REVIEW_AND_PCA.md`
