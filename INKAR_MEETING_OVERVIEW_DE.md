# INKAR: kurzer Datenueberblick fuer das Meeting

## 1. Was ist der INKAR-Datensatz in diesem Projekt?

Fuer die soziooekonomische Seite der Arbeit wurde der INKAR-2025-Datensatz des BBSR verwendet. Die Rohdatei liegt als grosse CSV im Long-Format vor:

- Pfad: `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/SOCIOECONOMIC.nosync/Downloads/inkar_2025/inkar_2025.csv`
- Dateistruktur: Long format
- Rohumfang: `63,426,864` Datenzeilen
- Bereiche im Rohdatensatz: `EU`, `LRB`, `SDG`, `ZOM`
- Anzahl Raumbezug-Kategorien: `45`

Fuer dieses Projekt ist der Datensatz also nicht von vornherein "Elbe-spezifisch". Elbe-spezifisch wird er erst durch die nachgelagerte raeumliche Auswahl.

## 2. Welche Teilmenge wurde fuer die Arbeit genutzt?

Zunaechst wurde aus dem Gesamtdatensatz die fuer das Projekt relevante kommunale Teilmenge gefiltert:

- `Bereich == "LRB"`
- `Raumbezug == "Gemeinden"`

Diese Teilmenge umfasst:

- `26,879,954` Zeilen
- `10,978` Gemeinden
- `176` Indikatoren
- Zeitbezug von `1995` bis `2024`

Danach wurde auf die fuer das Elbe-Einzugsgebiet relevanten Bundeslaender bzw. AGS-Praefixe eingeschraenkt:

- `01`, `02`, `03`, `09`, `11`, `12`, `13`, `14`, `15`, `16`

Die daraus entstandene Elbe-Teilmenge umfasst:

- `16,286,960` Zeilen
- `6,699` Gemeinden
- `176` Indikatoren
- Zeitbezug von `1995` bis `2024`

## 3. Wie wurde daraus der Arbeitsdatensatz gebaut?

Die weitere Verarbeitung lief in mehreren Schritten:

1. Pro `Gemeinde x Kuerzel` wurde immer die neueste verfuegbare Beobachtung behalten.
2. Danach wurde der Datensatz von Long nach Wide gepivotet.
3. Fuer jede Variable wurde zusaetzlich das jeweilige Bezugsjahr als `_year`-Spalte mitgenommen.
4. Anschliessend wurden nur die thematisch relevanten INKAR-Indikatoren behalten.
5. Erst danach wurden die Variablen im Projekt intern umbenannt.

Der daraus erzeugte bereinigte Arbeitsdatensatz ist:

- Datei: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/inkar_elbe_core_MAIN.csv`
- Umfang: `6,699` Zeilen
- Gesamtspalten: `118`
- Davon `58` Wertespalten
- Davon `58` Jahrspalten
- Plus `AGS` und `mun_name`

## 4. Wichtige Namenskette fuer die Praesentation

Fuer das Meeting ist wichtig, drei Ebenen auseinanderzuhalten:

1. Originaler INKAR-Kuerzelname im Rohdatensatz
2. Interner bereinigter Name nach dem Umbenennen
3. Spaeterer Analyse-Name mit `exposure_`-Praefix nach dem Join in den raeumlichen Exposure-Datensatz

Beispiel:

- Roh-INKAR: `a_ALGII_SGBII`
- bereinigter INKAR-Name: `share_alg2_sgb2`
- spaeter im Analyseobjekt: `exposure_share_alg2_sgb2`

Fuer das Word-Dokument und die Vorstellung vor den Profs sollte deshalb zuerst immer der originale INKAR-Name bzw. der offizielle Indikatorname genannt werden. Die umbenannten Analysevariablen sind nur Arbeitsnamen des Projekts.

## 5. Kernvariablen des Vulnerabilitaetsindex mit Originalnamen

| Rolle | Originales Kuerzel | Offizieller INKAR-Indikator | Interner Name | Spaeterer Analyse-Name |
| --- | --- | --- | --- | --- |
| core | `a_ALGII_SGBII` | `ALG II-Leistungen an SGBII` | `share_alg2_sgb2` | `exposure_share_alg2_sgb2` |
| core | `a_Unterkunft_SGBII` | `Leistungen fuer Unterkunft an SGBII` | `share_sgb2_with_housing_costs` | `exposure_share_sgb2_with_housing_costs` |
| core | `a_aloLang` | `Langzeitarbeitslose` | `share_longterm_unemp` | `exposure_share_longterm_unemp` |
| core | `q_alo_u25_einw` | `Juengere Arbeitslose` | `unemp_u25_per_1000` | `exposure_unemp_u25_per_1000` |
| core | `q_kaufkraft` | `Kaufkraft` | `purchasing_power` | `exposure_purchasing_power` |
| core | `a_hheink_niedrig` | `Haushalte mit niedrigem Einkommen` | `share_hh_income_low` | `exposure_share_hh_income_low` |
| core | `a_bev65um` | `Einwohner 65 Jahre und aelter` | `share_age_65plus` | `exposure_share_age_65plus` |
| core | `a_bev75um` | `Einwohner 75 Jahre und aelter` | `share_age_75plus` | `exposure_share_age_75plus` |
| core | `q_abhg_alt` | `Abhaengigenquote Alte` | `old_age_dependency` | `exposure_old_age_dependency` |
| core | `a_BG1P` | `Einpersonen-Bedarfsgemeinschaften` | `share_bg_single_parent` *(aktueller Projektname, fachlich fraglich)* | `exposure_share_bg_single_parent` |
| core | `q_HH1` | `Einpersonenhaushalte` | `share_single_households` | `exposure_share_single_households` |
| core | `m_OEV20_DIST` | `Entfernung zur OEV Haltestelle` | `dist_public_transport_m` | `exposure_dist_public_transport_m` |
| core | `m_Q07_HA_DIST` | `Entfernung zum Hausarzt` | `dist_gp_m` | `exposure_dist_gp_m` |
| core | `m_Q01_APO_DIST` | `Entfernung zur Apotheke` | `dist_pharmacy_m` | `exposure_dist_pharmacy_m` |

## 6. Optionale Robustheitsvariablen mit Originalnamen

| Rolle | Originales Kuerzel | Offizieller INKAR-Indikator | Interner Name | Spaeterer Analyse-Name |
| --- | --- | --- | --- | --- |
| optional | `q_ärzte_bev` | `Aerzte` | `doctors_total_per_1000` | `exposure_doctors_total_per_1000` |
| optional | `m_G02_SUP_DIST` | `Entfernung zum Supermarkt/Discounter` | `dist_supermarket_m` | `exposure_dist_supermarket_m` |
| optional | `a_bb_100Mbits` | `Bandbreitenverfuegbarkeit mindestens 100 Mbit/s` | `share_bb_100mbit` | `exposure_share_bb_100mbit` |

## 7. Kritischer Punkt fuer morgen: ein Rename ist sehr wahrscheinlich inhaltlich falsch

Beim Rueckcheck der Rohmetadaten ist ein wichtiger Widerspruch aufgefallen:

- `a_BG1P` ist im Roh-INKAR als `Einpersonen-Bedarfsgemeinschaften` bezeichnet.
- Im Projekt wurde diese Variable spaeter zu `share_bg_single_parent` umbenannt.

Das passt inhaltlich nicht zusammen.

`Einpersonen-Bedarfsgemeinschaften` sind nicht dasselbe wie `single-parent households`.

Das ist fuer das Meeting wichtig, weil du diese Variable in der aktuellen Form nicht einfach als "single-parent households" vorstellen solltest. Sauber waeren fuer morgen drei Optionen:

1. Die Variable im Meeting nur unter ihrem offiziellen INKAR-Namen vorstellen und die problematische interne Umbenennung offen benennen.
2. Die Variable vor der naechsten Analyse inhaltlich korrekt umbenennen.
3. Pruefen, ob stattdessen ein tatsaechlich passender Alleinerziehenden-Indikator verwendet werden soll, etwa `a_ewfBG_allein` oder ein anderer fachlich besser passender INKAR-Indikator.

## 8. Kann der Vulnerabilitaetsindex grundsaetzlich auf ganz Deutschland uebertragen werden?

Fuer den soziooekonomischen Teil lautet die kurze Antwort: im Prinzip ja.

Begruendung:

- Die verwendeten INKAR-Indikatoren sind nicht Elbe-spezifisch.
- Der Elbe-Bezug entsteht erst durch die raeumliche Filterung auf relevante AGS-Codes und spaeter durch den Join mit den Flood-Daten.
- Der Vulnerabilitaetsindex als soziooekonomischer Index kann deshalb prinzipiell auch fuer ganz Deutschland gebaut werden.

Die eigentliche Einschraenkung liegt nicht primaer im INKAR-Teil, sondern in der Hazard- und Protection-Seite:

- Welche Flood-Daten liegen bundesweit konsistent vor?
- Auf welcher raeumlichen Ebene liegen sie vor?
- Ist die Schutzinformation bundesweit vergleichbar?

## 9. Meeting-taugliche Kurzformulierung

Eine knappe Einleitung fuer morgen koennte so lauten:

"Der soziooekonomische Teil der Analyse basiert auf dem bundesweiten INKAR-2025-Datensatz des BBSR. Fuer die Arbeit wurde daraus zunaechst die kommunale LRB-Teilmenge gefiltert und anschliessend auf die fuer das Elbe-Einzugsgebiet relevanten Gemeinden eingeschraenkt. Die finalen Analysevariablen sind also projektinterne Umbenennungen; fuer die methodische Vorstellung sollten die originalen INKAR-Kuerzel und offiziellen Indikatornamen zuerst genannt werden."
