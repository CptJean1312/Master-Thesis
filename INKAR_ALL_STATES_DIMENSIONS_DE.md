# INKAR: Verfuegbare Gemeindeindikatoren in allen 16 Bundeslaendern

## Abgrenzung

Grundlage ist der originale INKAR-2025-Rohdatensatz auf Gemeindeebene (`Bereich = LRB`, `Raumbezug = Gemeinden`) nach Filterung auf alle 16 AGS-Praefixe der Bundeslaender (`01` bis `16`).

Verfuegbarkeitsdefinition:
- `in allen 16 Bundeslaendern vorhanden`: Ein Indikator hat in jedem Bundesland mindestens eine Gemeinde mit gueltigem Letztwert.
- `vollstaendig in allen 16 Bundeslaendern`: Der Letztwert ist fuer jede Gemeinde in jedem Bundesland vorhanden.

## Kernergebnis

- Gemeindeindikatoren insgesamt: `176`
- In allen 16 Bundeslaendern vorhanden: `173`
- Vollstaendig in allen 16 Bundeslaendern: `32`
- Nicht in allen 16 Bundeslaendern vorhanden: `3` (alle drei aus dem Bereich Oeffentliche Finanzen)

## Grobe soziooekonomische Dimensionen

- `Labour market and employment`: 37 Indikatoren, davon 7 mit Vollabdeckung
- `Demography and household structure`: 35 Indikatoren, davon 2 mit Vollabdeckung
- `Health, services and accessibility`: 33 Indikatoren, davon 6 mit Vollabdeckung
- `Housing and land-use context`: 27 Indikatoren, davon 8 mit Vollabdeckung
- `Income, deprivation and social transfers`: 22 Indikatoren, davon 0 mit Vollabdeckung
- `Basic population and structural counts`: 9 Indikatoren, davon 9 mit Vollabdeckung
- `Education`: 7 Indikatoren, davon 0 mit Vollabdeckung
- `Economic structure and tourism`: 3 Indikatoren, davon 0 mit Vollabdeckung

## Nicht in allen Bundeslaendern verfuegbare Indikatoren

- `q_sach` - Ausgaben fuer Sachinvestitionen
- `q_schluesselzuw` - Schluesselzuweisungen
- `q_investZ` - Zuweisungen fuer Investitionsfoerdermassnahmen

## Indikatoren nach grober soziooekonomischer Dimension

### Labour market and employment

**Arbeitslosigkeit**

_Arbeitslosigkeit - Allgemein_

- `a_alo_f` - Arbeitslose Frauen
- `a_alo_m` - Arbeitslose Männer

_Arbeitslosigkeit - Altersgruppen_

- `a_alo_u25` - Anteil jüngere Arbeitslose
- `a_alo_u25_m` - Anteil männliche jüngere Arbeitslose
- `a_alo_ü55_m` - Anteil männliche ältere Arbeitslose
- `a_alo_u25_f` - Anteil weibliche jüngere Arbeitslose
- `a_alo_ü55_f` - Anteil weibliche ältere Arbeitslose
- `a_alo_ü55` - Anteil ältere Arbeitslose
- `q_alo_u25_einw` - Jüngere Arbeitslose
- `q_alo_u25_einw_m` - Männliche jüngere Arbeitslose
- `q_alo_ü55_einw_m` - Männliche ältere Arbeitslose
- `q_alo_u25_einw_f` - Weibliche jüngere Arbeitslose
- `q_alo_ü55_einw_f` - Weibliche ältere Arbeitslose
- `q_alo_ü55_einw` - Ältere Arbeitslose

_Arbeitslosigkeit - Struktur_

- `a_alo_experte` - Arbeitslose mit Anforderungsniveau Experte
- `a_alo_fachkraft` - Arbeitslose mit Anforderungsniveau Fachkraft
- `a_alo_helfer` - Arbeitslose mit Anforderungsniveau Helfer
- `a_alo_spezialist` - Arbeitslose mit Anforderungsniveau Spezialist
- `a_alo_oAusb` - Arbeitslose ohne Ausbildung
- `a_alo_ausländer` - Ausländische Arbeitslose
- `a_alo_ausländer_m` - Ausländische männliche Arbeitslose
- `a_alo_ausländer_f` - Ausländische weibliche Arbeitslose
- `a_aloLang` - Langzeitarbeitslose
- `a_aloLang_m` - Männliche Langzeitarbeitslose
- `a_aloLang_f` - Weibliche Langzeitarbeitslose

**Beschäftigung und Erwerbstätigkeit**

_Beschäftigung und Erwerbstätigkeit - Atypische Beschäftigung_

- `a_gb_Frauen` - Anteil Minijobs (Frauen) [vollstaendig]
- `a_gb_Männer` - Anteil Minijobs (Männer) [vollstaendig]
- `a_gb_nj` - Anteil Minijobs (Nebenverdienst) [vollstaendig]
- `a_gb_nj_prBv` - Anteil Minijobs (Nebenverdienst) an den Beschäftigungsverhältnissen [vollstaendig]
- `a_gb_knj_prBv` - Anteil Minijobs (ausschließlich) an den Beschäftigungsverhältnissen [vollstaendig]
- `a_gb_knj` - Anteil Minijobs (ausschließlich) an geringfügig Beschäftigten [vollstaendig]
- `a_Minijobs` - Anteil Minijobs an den Beschäftigungsverhältnissen [vollstaendig]

_Beschäftigung und Erwerbstätigkeit – Struktur_

- `q_svw` - Beschäftigtenquote
- `q_svw_f` - Beschäftigtenquote Frauen
- `q_svw_m` - Beschäftigtenquote Männer
- `q_svw_ü55` - Quote ältere Beschäftigte
- `r_ewf_jungalt` - Verhältnis junge zu alte Erwerbsfähige

### Demography and household structure

**Bevölkerung**

_Bevölkerung - Natürliche Bevölkerungsbewegungen_

- `a_geb_bev` - Geborene
- `a_gest_bev` - Gestorbene
- `i_saldo_nat` - Natürlicher Saldo

_Bevölkerung – Altersstruktur_

- `m_bev_alter` - Durchschnittsalter der Bevölkerung
- `a_bev65um` - Einwohner 65 Jahre und älter
- `a_bev65um_f` - Einwohner 65 Jahre und älter, Frauen
- `a_bev75um` - Einwohner 75 Jahre und älter
- `a_bev75um_f` - Einwohner 75 Jahre und älter, Frauen
- `a_bev0003` - Einwohner unter 3 Jahren
- `a_bev_0006` - Einwohner unter 6 Jahre
- `a_bev1825` - Einwohner von 18 bis unter 25 Jahren
- `a_bev1825_f` - Einwohner von 18 bis unter 25 Jahren, Frauen
- `a_bev2530` - Einwohner von 25 bis unter 30 Jahren
- `a_bev2530_f` - Einwohner von 25 bis unter 30 Jahren, Frauen
- `a_bev0306` - Einwohner von 3 bis unter 6 Jahren
- `a_bev3050` - Einwohner von 30 bis unter 50 Jahren
- `a_bev5065` - Einwohner von 50 bis unter 65 Jahren
- `a_bev0618` - Einwohner von 6 bis unter 18 Jahren
- `a_bev6575` - Einwohner von 65 bis unter 75 Jahren
- `a_bev6575_f` - Einwohner von 65 bis unter 75 Jahren, Frauen

_Bevölkerung – Bevölkerungsstruktur_

- `q_abhg_alt` - Abhängigenquote Alte
- `q_abhg_jung` - Abhängigenquote Junge
- `e10_bev` - Bevölkerungsentwicklung (10 Jahre)
- `e5_bev` - Bevölkerungsentwicklung (5 Jahre)
- `q_HH1` - Einpersonenhaushalte
- `a_bevf` - Frauenanteil
- `a_bevf_2040` - Frauenanteil 20 bis unter 40 Jahre
- `a_hh_kind` - Haushalte mit Kindern
- `q_HH` - Haushaltsgröße

_Bevölkerung – Wanderungen_

- `i_wans` - Gesamtwanderungssaldo

**Siedlungsstruktur**

- `a_bevMZ` - Bevölkerung in Mittelzentren
- `a_bevOZ` - Bevölkerung in Oberzentren
- `q_bevsva_qkm` - Einwohner-Arbeitsplatz-Dichte [vollstaendig]
- `q_bev_fl` - Einwohnerdichte [vollstaendig]
- `a_bevd150_` - Ländlichkeit

### Health, services and accessibility

**Erreichbarkeit**

- `a_bb_1000Mbits` - Bandbreitenverfügbarkeit mindestens 1.000 Mbit/s [vollstaendig]
- `a_bb_100Mbits` - Bandbreitenverfügbarkeit mindestens 100 Mbit/s [vollstaendig]
- `a_bb_50Mbits` - Bandbreitenverfügbarkeit mindestens 50 Mbit/s [vollstaendig]
- `m_Q07_HA_DIST` - Entfernung zum Hausarzt
- `m_G02_SUP_DIST` - Entfernung zum Supermarkt/Discounter
- `m_Q01_APO_DIST` - Entfernung zur Apotheke
- `m_P01_PRIM_DIST` - Entfernung zur Grundschule
- `m_OEV20_DIST` - Entfernung zur ÖV Haltestelle
- `m_ErrBAB_bev` - Erreichbarkeit von Autobahnen
- `m_ErrAir_bev` - Erreichbarkeit von Flughäfen
- `m_ErrIC_bev` - Erreichbarkeit von IC/EC/ICE-Bahnhöfen
- `i_fz_mz` - Erreichbarkeit von Mittelzentren [vollstaendig]
- `i_fz_oz` - Erreichbarkeit von Oberzentren [vollstaendig]
- `a_Q01_APO_ANT` - Nahversorgung Apotheke
- `a_P01_PRIM_ANT` - Nahversorgung Grundschule
- `a_Q07_HA_ANT` - Nahversorgung Hausarzt
- `a_G02_SUP_ANT` - Nahversorgung Supermarkt/Discounter
- `a_OEV20_ANT` - Nahversorgung ÖV Haltestelle

_Digitale Erreichbarkeit_

- `a_bb_4G` - 4G-Abdeckung (metadata label missing in overview sheet)

_Verkehr und Erreichbarkeit - Pendler_

- `a_ausp_svw` - Auspendler [vollstaendig]
- `a_einp_sva` - Einpendler
- `a_pend150` - Pendler mit Arbeitsweg 150 km und mehr
- `a_pend300` - Pendler mit Arbeitsweg 300 km und mehr
- `a_pend50` - Pendler mit Arbeitsweg 50 km und mehr
- `q_pendlersaldo` - Pendlersaldo

_Verkehr und Erreichbarkeit - Straßenverkehr_

- `q_vunpt_bev` - Getötete im Straßenverkehr
- `q_unf_bev` - Straßenverkehrsunfälle
- `q_vunp_bev` - Verunglückte im Straßenverkehr

**Medizinische Versorgung**

- `q_allgemeinärzte_bev` - Allgemeinärzte
- `q_hausarzt_bev` - Hausärzte
- `q_internist_bev` - Internisten
- `q_kinderarzt_kinder` - Kinderärzte
- `q_ärzte_bev` - Ärzte

### Housing and land-use context

**Bauen und Wohnen**

_Bauen und Wohnen - Baulandmarkt und Bautätigkeit_

- `a_gen_wo12` - Baugenehmigungen für Wohnungen in Ein- und Zweifamilienhäusern
- `a_gen_wo3um` - Baugenehmigungen für Wohnungen in Mehrfamilienhäusern
- `q_gen_wo_ew` - Baugenehmigungen für Wohnungenje Einwohner
- `a_fert_wohn` - Fertiggestellte Wohnungen je Wohnung im Bestand
- `a_fert_wo12` - Neubauwohnungen in Ein- und Zweifamilienhäusern
- `q_fert_wo12_bev` - Neubauwohnungen in Ein- und Zweifamilienhäusern je Einwohner
- `q_fert_wo3um_bev` - Neubauwohnungen in Mehrfamilienhäusern je Einwohner
- `q_fert_wo_bev` - Neubauwohnungen je Einwohner
- `a_fert_wg12` - Neue Ein- und Zweifamilienhäuser

_Bauen und Wohnen - Gebäude- und Wohnungsbestand_

- `a_wo_r5um` - 5- und mehr Raum-Wohnungen
- `a_wg_wo12` - Ein- und Zweifamilienhäuser
- `a_wo_r12` - Ein- und Zweiraumwohnungen
- `a_wg_wo3um` - Mehrfamilienhäuser
- `a_wo_wg12` - Wohnungen in Ein- und Zweifamilienhäusern
- `a_wo_wg3um` - Wohnungen in Mehrfamilienhäusern

**Flächennutzung**

- `a_erh` - Erholungsfläche [vollstaendig]
- `q_erhfl_bev` - Erholungsfläche je Einwohner
- `v_suvfl_vorjahr` - Flächenneuinanspruchnahme [vollstaendig]
- `a_freifläche` - Freifläche [vollstaendig]
- `q_freifläche_bev` - Freifläche je Einwohner
- `a_landwirtschaft` - Landwirtschaftsfläche [vollstaendig]
- `a_naturnah` - Naturnähere Fläche [vollstaendig]
- `q_naturnah_bev` - Naturnähere Fläche je Einwohner
- `a_suv_fl` - Siedlungs- und Verkehrsfläche [vollstaendig]
- `q_ew_suv_qkm` - Siedlungsdichte in km²
- `a_wald` - Waldfläche [vollstaendig]
- `a_wasser` - Wasserfläche [vollstaendig]

### Income, deprivation and social transfers

**Privateinkommen und private Schulden**

- `q_einzelhandelskaufkraft` - Einzelhandelsrelevante Kaufkraft
- `a_hheink_hoch` - Haushalte mit hohem Einkommen
- `a_hheink_mittel` - Haushalte mit mittlerem Einkommen
- `a_hheink_niedrig` - Haushalte mit niedrigem Einkommen
- `q_kaufkraft` - Kaufkraft

**Sozialleistungen**

_Sozialleistungen - Transferleistungen_

- `a_ALGII_SGBII` - ALG II-Leistungen an SGBII
- `a_Unterkunft_SGBII` - Leistungen für Unterkunft an SGBII

_Sozialleistungen – Bedarfsgemeinschaften_

- `a_ewfBG_allein` - Alleinerziehende erwerbsfähige Leistungsberechtigte
- `a_BGKind` - Bedarfsgemeinschaften mit Kindern
- `a_BG1P` - Einpersonen-Bedarfsgemeinschaften
- `a_ewfBG` - Erwerbsfähige Leistungsberechtigte
- `a_ewfBG_f` - Erwerbsfähige Leistungsberechtigte (Frauen)
- `a_BG5um` - Große Bedarfsgemeinschaften
- `a_ewfBG_u25_` - Junge erwerbsfähige Leistungsberechtigte
- `q_newfBGu15_bev` - Kinderarmut
- `q_PBG_bev` - Personen in Bedarfsgemeinschaften
- `a_ewfBG_55um_` - Ältere erwerbsfähige Leistungsberechtigte

**Öffentliche Finanzen**

- `q_einkst_bev` - Einkommensteuer
- `q_gewst_bev` - Gewerbesteuer
- `d_steuereinnahme` - Steuereinnahmen
- `q_gest_bev` - Steuerkraft
- `q_umsst` - Umsatzsteuer

### Basic population and structural counts

**Absolutzahlen**

- `alo` - Arbeitslose [vollstaendig]
- `bev_korr` - Bevölkerung (mit BBSR-Zensuskorrekturen) [vollstaendig]
- `xbev` - Bevölkerung gesamt [vollstaendig]
- `xbevm` - Bevölkerung männlich [vollstaendig]
- `xbevf` - Bevölkerung weiblich [vollstaendig]
- `TN23-kataster_qkm` - Bodenfläche gesamt qkm [vollstaendig]
- `ewf_1565_ges` - Erwerbsfähige Bevölkerung (15 bis unter 65 Jahre) [vollstaendig]
- `sva` - Sozialversicherungspflichtig Beschäftigte am Arbeitsort [vollstaendig]
- `svw` - Sozialversicherungspflichtig Beschäftigte am Wohnort [vollstaendig]

### Education

**Bildung**

_Bildung – Ausbildungsangebot_

- `a_stud_a` - Ausländische Studierende
- `a_stud_m` - Männliche Studierende
- `q_stud` - Studierende
- `q_stud_fh` - Studierende an FH
- `a_stud_1` - Studierende im 1. Semester
- `q_stud_1825` - Studierende je 100 Einwohner 18 bis 25 Jahre
- `a_stud_w` - Weibliche Studierende

### Economic structure and tourism

**Wirtschaft**

_Wirtschaft – Fremdenverkehr_

- `q_übern_bev` - Gästeübernachtungen in Beherbergungsbetrieben
- `q_schlafg_bev` - Schlafgelegenheiten in Beherbergungsbetrieben
- `m_übern` - Verweildauer in Beherbergungsbetrieben

