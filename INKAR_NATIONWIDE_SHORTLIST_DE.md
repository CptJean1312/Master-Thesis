# INKAR: bundesweite AGS-Shortlist fuer einen moeglichen Vulnerabilitaetsindex

## Ausgangspunkt

Nach Filterung auf Gemeindeebene (`Bereich = LRB`, `Raumbezug = Gemeinden`) und alle 16 Bundeslaender gilt:

- `176` Gemeindeindikatoren insgesamt
- `173` sind in allen 16 Bundeslaendern vorhanden
- aber nur `32` haben wirklich Vollabdeckung auf AGS-Ebene

Das ist methodisch wichtig:

- `bundesweit vorhanden` ist nicht dasselbe wie `fuer jede Gemeinde vollstaendig vorhanden`
- wenn man nur Variablen mit strikter Vollabdeckung zulassen wuerde, waere der Kandidatenpool fuer einen sozial sinnvollen Vulnerabilitaetsindex ziemlich duenn

## 1. Strikte bundesweite Vollabdeckung: Was ist wirklich ueberall da?

Diese Variablen sind fuer alle Gemeinden in allen Bundeslaendern verfuegbar.

### Potenziell sinnvoll fuer einen Index oder als Kern-Proxy

- `a_Minijobs` - Anteil Minijobs an den Beschaeftigungsverhaeltnissen
- `a_gb_Frauen` - Anteil Minijobs (Frauen)
- `a_gb_Maenner` / `a_gb_Männer` - Anteil Minijobs (Maenner)
- `a_gb_nj` - Anteil Minijobs (Nebenverdienst)
- `a_gb_nj_prBv` - Anteil Minijobs (Nebenverdienst) an den Beschaeftigungsverhaeltnissen
- `a_gb_knj_prBv` - Anteil Minijobs (ausschliesslich) an den Beschaeftigungsverhaeltnissen
- `a_gb_knj` - Anteil Minijobs (ausschliesslich) an geringfuegig Beschaeftigten
- `a_bb_100Mbits` - Bandbreitenverfuegbarkeit mindestens 100 Mbit/s
- `a_bb_50Mbits` - Bandbreitenverfuegbarkeit mindestens 50 Mbit/s
- `a_bb_1000Mbits` - Bandbreitenverfuegbarkeit mindestens 1.000 Mbit/s
- `i_fz_mz` - Erreichbarkeit von Mittelzentren
- `i_fz_oz` - Erreichbarkeit von Oberzentren

### Besser als Kontext- oder Kontrollvariablen als als Kernindex

- `q_bev_fl` - Einwohnerdichte
- `q_bevsva_qkm` - Einwohner-Arbeitsplatz-Dichte
- `xbev` / `bev_korr` - Bevoelkerung gesamt
- `xbevf`, `xbevm` - Geschlechterstruktur in Absolutzahlen
- `ewf_1565_ges` - Erwerbsfaehige Bevoelkerung
- `sva`, `svw`, `alo` - absolute Struktur- bzw. Bestandsgroessen

## 2. Warum das fuer einen Vulnerabilitaetsindex allein nicht reicht

Wenn man nur die Variablen mit echter Vollabdeckung nimmt, fehlen fast alle klassischen Vulnerabilitaetsindikatoren aus den Bereichen:

- materielle Deprivation
- Sozialleistungen
- Alterung und Abhaengigkeit
- Einkommensschwaeche
- Einpersonenhaushalte
- wohnortnahe Gesundheits- und Alltagsversorgung

Das heisst:

- eine `strict full coverage only`-Strategie ist datenlogisch sauber
- inhaltlich waere sie aber fuer einen guten Vulnerabilitaetsindex eher zu schwach und zu strukturorientiert

## 3. Pragmatic nationwide shortlist: fachlich sinnvoll, wenn kleine Missingness toleriert wird

Wenn geringe Fehlstellen akzeptiert und transparent behandelt werden, ist diese Shortlist deutlich sinnvoller fuer einen bundesweiten Vulnerabilitaetsindex:

### Deprivation

- `a_ALGII_SGBII` - ALG II-Leistungen an SGBII
- `a_Unterkunft_SGBII` - Leistungen fuer Unterkunft an SGBII
- `a_aloLang` - Langzeitarbeitslose
- `q_alo_u25_einw` - Juengere Arbeitslose
- `q_kaufkraft` - Kaufkraft
- `a_hheink_niedrig` - Haushalte mit niedrigem Einkommen

### Demographic sensitivity

- `a_bev65um` - Einwohner 65 Jahre und aelter
- `a_bev75um` - Einwohner 75 Jahre und aelter
- `q_abhg_alt` - Abhaengigenquote Alte

### Household / social structure

- `q_HH1` - Einpersonenhaushalte

### Accessibility / adaptive capacity

- `m_OEV20_DIST` - Entfernung zur OeV Haltestelle
- `m_Q07_HA_DIST` - Entfernung zum Hausarzt
- `m_Q01_APO_DIST` - Entfernung zur Apotheke

### Optional robustness variables

- `q_aerzte_bev` / `q_ärzte_bev` - Aerzte
- `m_G02_SUP_DIST` - Entfernung zum Supermarkt/Discounter
- `a_bb_100Mbits` - Bandbreitenverfuegbarkeit mindestens 100 Mbit/s

## 4. Eine Variable wuerde ich aktuell nicht ungeprueft in einen nationalen Index uebernehmen

- `a_BG1P`

Grund:

- offizieller INKAR-Name: `Einpersonen-Bedarfsgemeinschaften`
- aktueller Projektname: `share_bg_single_parent`

Diese Bedeutungen passen nicht zusammen. Solange das nicht sauber geklaert ist, wuerde ich diese Variable nicht in eine deutschlandweite Kernshortlist aufnehmen.

## 5. Meine klare Empfehlung

Wenn du im Meeting eine `erste deutschlandweite Vorschlagsliste` nennen willst, wuerde ich nicht die strikte Vollabdeckungsliste als Hauptvorschlag nehmen.

Ich wuerde stattdessen sagen:

- fuer einen bundesweiten Index gibt es eine kleine Zahl an Variablen mit kompletter Vollabdeckung
- diese sind aber fuer einen guten Vulnerabilitaetsindex inhaltlich zu schmal und oft zu strukturorientiert
- sinnvoller ist eine `pragmatic nationwide shortlist` mit fachlich starken Variablen und geringer, transparent behandelter Missingness

## 6. Kurzform fuer das Meeting

`If we require absolute completeness on AGS level across all of Germany, the candidate pool becomes too narrow and too structural. If we allow a small and transparent amount of missingness, we can build a much more defensible nationwide vulnerability index around deprivation, ageing, household structure and accessibility.`
