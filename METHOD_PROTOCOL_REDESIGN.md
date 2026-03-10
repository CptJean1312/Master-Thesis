# Method Protocol Redesign

## Ziel

Dieses Protokoll beschreibt die Logik hinter [Analyse_PROTOCOL_REDESIGN.R](/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/Analyse_PROTOCOL_REDESIGN.R).

Der zentrale Unterschied zum alten Stand ist:

- Vulnerabilitaet wird nicht mehr als ein grosser, schwer interpretierbarer Sammelindex aus sehr vielen gemischten INKAR-Variablen gebaut.
- Stattdessen wird sie als theoriebasierte Kombination aus vier Domänen operationalisiert:
  - soziooekonomische Deprivation
  - demografische Sensitivitaet
  - Haushalts- und Sozialstruktur
  - alltagsbezogene Erreichbarkeit und adaptive Kapazitaet
- GISD wird bewusst noch nicht verwendet. Es bleibt fuer einen spaeteren Validierungs- oder Robustheitsschritt reserviert.
- Die finale Variablenauswahl folgt jetzt dem zuletzt abgestimmten Excel-Stand: Kernindex, optionale Sensitivitaetsvariablen und ausgeschlossene Kontextindikatoren sind explizit getrennt.

## Konzeptuelle Basis

Die Neuausrichtung folgt einer sauberen Trennung zwischen Hazard, Exposure, Vulnerability und Protection Context.

- Hazard/Exposure:
  - Hauptoutcome bleibt `Zone 1 share` aus dem BfG-HQ100-Datensatz.
  - Das ist keine vollstaendige "Risk"-Messung, sondern eine Messung ungeschuetzter Flaechenexposition.
- Vulnerability:
  - soll die soziale Anfaelligkeit bzw. begrenzte Bewaeltigungs- und Anpassungskapazitaet abbilden
  - soll daher nicht mit reiner Urbanitaet, Gemeindegroesse oder fiskalischer Struktur vermischt werden
- Protection:
  - wird nicht mehr ueber den alten Buffer-Proxy in das neue Hauptmodell eingespeist
  - stattdessen wird der EU-Flood/NUTS3-Schutzstandard als Schutzkontext verwendet

Literaturbasis:

- IPCC AR6 WGII Chapter 1: Vulnerabilitaet umfasst Sensitivitaet/Suszeptibilitaet und mangelnde Kapazitaet zum coping und adapting.
  - https://www.ipcc.ch/report/ar6/wg2/chapter/chapter-1/
- Hinkel (2011): Vulnerabilitaetsindikatoren sollten problemorientiert und theoriegeleitet gebaut werden; ein unspezifischer "one-size-fits-all"-Index ist konzeptionell schwach.
  - https://doi.org/10.1016/j.gloenvcha.2010.08.002
- Rufat et al. (2015): In der Flood-Literatur sind demografische Merkmale, soziooekonomischer Status und Gesundheit die robustesten Treiber sozialer Vulnerabilitaet; coping capacity sollte explizit beruecksichtigt werden.
  - https://doi.org/10.1016/j.ijdrr.2015.09.013
- Fekete (2009): Fuer Deutschland zeigen sich bei Flussfluten besonders Aeltere und finanziell schwaechere Gruppen als relevante Risikogruppen.
  - https://doi.org/10.5194/nhess-9-393-2009
- Cutter et al. (2003): Soziale Vulnerabilitaet ist multidimensional; Faktoren koennen je nach Raum unterschiedlich wirken, deshalb ist Transparenz in Auswahl und Aggregation zentral.
  - https://doi.org/10.1111/1540-6237.8402002
- Tate (2012): Indexergebnisse sind stark sensitiv gegenueber Variablenauswahl und Gewichtung; deshalb sollte der Konstruktionsprozess transparent und mit Sensitivitaetschecks aufgebaut werden.
  - https://doi.org/10.1007/s11069-012-0152-2

## Ausgewaehlte INKAR-Indikatoren

### 1. Soziooekonomische Deprivation

Verwendete Variablen:

- `exposure_share_alg2_sgb2`
- `exposure_share_sgb2_with_housing_costs`
- `exposure_share_longterm_unemp`
- `exposure_unemp_u25_per_1000`
- `exposure_purchasing_power` (invertiert)
- `exposure_share_hh_income_low`

Begruendung:

- Diese Gruppe bildet materielle Ressourcenknappheit, Arbeitsmarktbenachteiligung und eingeschraenkte finanzielle Bewaeltigungskapazitaet ab.
- Genau diese Dimension ist in der Flood-Vulnerability-Literatur am stabilsten.
- Armut, Transferbezug und niedrige Kaufkraft sind plausible Marker fuer geringere Vorsorge-, Anpassungs- und Erholungskapazitaet.
- Langzeitarbeitslosigkeit und Jugend-Arbeitslosigkeit erfassen prekäre Erwerbslagen besser als eine einzige allgemeine Arbeitslosenquote.

Literaturbasis:

- Rufat et al. (2015)
- Fekete (2009)
- Cutter et al. (2003)

### 2. Demografische Sensitivitaet

Verwendete Variablen:

- `exposure_share_age_65plus`
- `exposure_share_age_75plus`
- `exposure_old_age_dependency`

Begruendung:

- Aeltere Bevoelkerung ist bei Flood-Ereignissen empirisch haeufig mit hoeherer Betroffenheit verbunden.
- Relevant sind Mobilitaetseinschraenkungen, hoehere gesundheitliche Belastbarkeitsschwellen und ein hoeherer Bedarf an externer Hilfe im Warn-, Evakuierungs- und Recovery-Prozess.
- Der Old-Age-Dependency-Wert ergaenzt die reinen Altersanteile um eine strukturelle Lastenrelation.

Literaturbasis:

- Rufat et al. (2015)
- Fekete (2009)
- Cutter et al. (2003)

### 3. Haushalts- und Sozialstruktur

Verwendete Variablen:

- `exposure_share_bg_single_parent`
- `exposure_share_single_households`

Begruendung:

- Diese Variablen dienen als Proxy fuer begrenzte innerhaeusliche Unterstuetzungsressourcen.
- Einpersonenhaushalte und Einelternhaushalte koennen in Warnung, Evakuierung, Alltagspuffer und Recovery weniger interne Redundanz haben.
- Diese Dimension ist nicht identisch mit Armut, sondern bildet soziale Organisation und Verfuegbarkeit privater Hilfe ab.

Literaturbasis:

- Cutter et al. (2003)
- Rufat et al. (2015)

### 4. Erreichbarkeit und adaptive Kapazitaet

Verwendete Variablen:

- `exposure_dist_public_transport_m`
- `exposure_dist_gp_m`
- `exposure_dist_pharmacy_m`

Begruendung:

- Diese Gruppe operationalisiert keine Armut, sondern alltagspraktische Anpassungs- und Bewaeltigungskapazitaet.
- Lange Wege zu OePNV, hausarztlicher Versorgung und Apotheke deuten auf geringere Alltagsresilienz und schlechtere Erreichbarkeit von Hilfs- und Gesundheitsinfrastruktur hin.
- Die reduzierte Fassung dieser Domäne vermeidet, dass der Kernindex zu stark laendliche Versorgungsstruktur oder allgemeine Infrastrukturmodernisierung abbildet.

Literaturbasis:

- IPCC AR6 WGII Chapter 1
- Rufat et al. (2015)
- Hinkel (2011)

## Optionale Robustheitsvariablen

Die folgenden Variablen werden bewusst nicht im Hauptindex verwendet, aber fuer Sensitivitaetschecks im Code vorgehalten:

- `exposure_doctors_total_per_1000` (invertiert)
- `exposure_dist_supermarket_m`
- `exposure_share_bb_100mbit` (invertiert)

Begruendung:

- Alle drei Variablen sind inhaltlich plausibel, aber konzeptionell indirekter als die Kernindikatoren.
- Sie koennen teilweise laendliche Lage, Versorgungsorganisation oder generelle Infrastrukturmodernisierung mitabbilden.
- Genau deshalb eignen sie sich besser fuer Robustheitsanalysen als fuer den definitorischen Kern des Hauptindex.

## Was bewusst nicht im neuen Hauptindex ist

Die folgenden Variablentypen wurden im Redesign bewusst nicht in den Hauptindex aufgenommen:

- `exposure_pop_density_per_km2`
- `exposure_employment_density_per_km2`
- `exposure_pop_total`
- `exposure_tax_revenue_total`
- `exposure_trade_tax_per_capita`
- `exposure_income_tax_per_capita`
- viele fein aufgeloeste Altersgruppen
- mehrere parallele Einkommens- oder Breitbandstufen
- mehrere hoch ueberlappende Arzt-Spezialisierungen
- Studierendenvariablen
- Migrations- und natuerliche Bevoelkerungsbilanz

Begruendung:

- Dichte und Bevoelkerung sind eher Siedlungs- bzw. Raumstruktur als soziale Vulnerabilitaet im engeren Sinn. Deshalb werden sie als Kontrollvariablen verwendet, nicht als Indexbestandteil.
- Fiskalvariablen auf Gemeindeebene sagen eher etwas ueber kommunale Struktur als ueber household-level vulnerability.
- Auch `income_tax_per_capita` ist eher ein allgemeiner Wohlstands- bzw. Strukturindikator und fuegt dem Kernset neben Kaufkraft und Niedrigeinkommensanteil keinen zwingenden Zusatznutzen hinzu.
- Viele Alters-, Einkommens- und Infrastrukturvariablen sind mathematisch oder inhaltlich redundant. Sie wuerden in einer globalen PCA vor allem Urbanitaet und Datenskalierung verstaerken.
- Studierendenvariablen sind fuer das Flood-Vulnerability-Konstrukt in deinem Setting zu ambivalent.
- Bevoelkerungsdynamik und Migration sind interessante Kontextindikatoren, aber keine sauberen Kernmarker fuer Vulnerabilitaet in dieser Arbeit.

## Warum Domain-PCA statt einer einzigen grossen PCA

Die neue Logik ist:

1. Erst theoriebasiert die relevanten Vulnerabilitaetsdimensionen bestimmen.
2. Dann innerhalb jeder Domäne mit PCA auf eine robuste Kernachse reduzieren.
3. Danach die Domänenscores transparent zu einem Gesamtindex kombinieren.

Begruendung:

- Eine einzige PCA ueber sehr viele Variablen mischt leicht Deprivation, Alterung, ländliche Lage, Urbanitaet und Infrastruktur in einer schwer interpretierbaren Hauptkomponente.
- Domain-PCA reduziert Redundanz, ohne das theoretische Signal zu verlieren.
- Die Domänenscores bleiben inhaltlich lesbar und kartierbar.
- Gleichzeitig bleibt ein datengetriebener Schritt erhalten, weil innerhalb jeder Domäne die gemeinsame Struktur empirisch verdichtet wird.

Methodische Basis:

- Cutter et al. (2003) zeigen die Nützlichkeit faktorenanalytischer Reduktion, betonen aber die Multidimensionalitaet.
- Tate (2012) zeigt, dass Indikatorauswahl und Gewichtung die Ergebnisse stark beeinflussen; deshalb ist ein transparenter, hierarchischer Aufbau robuster als eine Black-Box-Gesamt-PCA.
- Hinkel (2011) argumentiert gegen unspezifische Vulnerabilitaetsindikatoren ohne klar definierten Zweck.

## Warum es trotzdem noch eine PCA-Sensitivitaet im Code gibt

Zusatzlich zum Domain-Hauptindex wird eine `selected_pca`-Sensitivitaet gerechnet.

Logik:

- gleiche, bereits theoriegefilterte Indikatoren
- PCA ueber das reduzierte Set
- nur Komponenten mit `eigenvalue > 1` und inhaltlicher Naehe zum Deprivationsanker werden bevorzugt behalten

Das ist kein zweiter Hauptindex, sondern ein Robustheitstest:

- Wenn Hauptindex und Sensitivitaetsindex aehnliche raeumliche Muster und Modellresultate liefern, ist die Konstruktion weniger artefaktgetrieben.

## Schritt-fuer-Schritt: Was der neue Code macht und warum

### Schritt 1: Daten einlesen und harmonisieren

Was:

- `gemeinden_elbe_final_full.gpkg` laden
- AGS konsolidieren
- leere Strings zu `NA`
- `_year`-Spalten numerisch machen
- nur geeignete Spalten numerisch parsen

Warum:

- Das reduziert Join-Artefakte aus QGIS/R und stellt sicher, dass spaetere Modellierung nicht an gemischten Datentypen scheitert.

### Schritt 2: Hazardmetriken definieren

Was:

- `share_total_hazard`
- `share_zone1`
- `share_zone2`
- `share_zone3_bfg`
- Hauptoutcome `risk_zone1_share = share_zone1`

Warum:

- Das folgt deinem aktualisierten Methodenprotokoll.
- Zone 1 ist der sauberste und konsistenteste BfG-basierte Indikator fuer ungeschuetzte HQ100-Exposition.

### Schritt 3: EU-Schutzkontext integrieren

Was:

- `peseta4_protection_nuts3.shp` wird auf Gemeinden gejoint
- Join ueber `point-on-surface`
- Ergebnis: `eu_protection_years`, `eu_protection_class`

Warum:

- Der BfG-Zone-3-Layer ist laut deiner Dokumentation und Rueckmeldung inkonsistent.
- Der EU-Datensatz wird deshalb nicht als exakte lokale Schutzwirkung interpretiert, sondern als regionaler Schutzkontext.

### Schritt 4: Tiefensensitivitaet aus AGG1 bauen

Was:

- `AGG1_flood_by_zone_depth.gpkg` wird tabellarisch eingelesen
- fuer `depth_class 1-5` werden Tiefenmittelwerte hinterlegt
- daraus wird eine `zone1_depth_weighted_exposure` berechnet

Warum:

- Die reine Flaechenquote sagt nichts ueber Ereignisintensitaet.
- Die tiefengewichtete Sensitivitaet ergaenzt deshalb die Hauptmetrik um eine Intensitaetsdimension.
- Klassen `6` und `7` werden nicht in den Tiefenscore eingerechnet, sondern explizit als unknown-depth-Anteil dokumentiert.

### Schritt 5: INKAR-Indikatoren in analytische Domänen uebersetzen

Was:

- Rohvariablen werden auf theoretische Analysevariablen gemappt
- einige werden invertiert, damit immer gilt: hoeher = verletzlicher

Warum:

- Das schafft eine konsistente Richtung fuer die spaetere Aggregation.
- Gleichzeitig wird explizit dokumentiert, was jede Variable inhaltlich bedeuten soll.
- Dabei wird zwischen Kernvariablen fuer den Hauptindex und optionalen Robustheitsvariablen fuer Sensitivitaetslaeufe getrennt.

### Schritt 6: Year-Audit

Was:

- fuer jede verwendete INKAR-Variable wird das Jahrsspektrum aus den `_year`-Feldern exportiert

Warum:

- Deine INKAR-Daten sind heterogen aktualisiert.
- Das muss transparent sein und ist methodisch berichtspflichtig.

### Schritt 7: Domain-PCA rechnen

Was:

- je Domäne PCA auf standardisierten, median-imputierten Variablen
- falls nur eine Variable bleibt: einfacher z-Wert
- Vorzeichenfix, damit hohe Werte immer hohe Vulnerabilitaet bedeuten

Warum:

- So wird Redundanz reduziert, ohne die Theorie aufzugeben.

### Schritt 8: Hauptindex bauen

Was:

- Mittelwert der vier Domain-scores
- z-Standardisierung und 0-1-Reskalierung

Warum:

- Jede Domäne bekommt bewusst gleiches Gewicht.
- Das vermeidet, dass eine grosse, datenreiche Domäne automatisch den gesamten Index dominiert.

### Schritt 9: PCA-Sensitivitaet bauen

Was:

- zweite, datengetriebene Verdichtung ueber das reduzierte Kernset plus optionale Robustheitsvariablen
- nur als Sensitivitaetsindex

Warum:

- Reaktion auf die von Tate betonte Sensitivitaet gegenueber Designentscheidungen.

### Schritt 10: Deskriptive Justice-Tabellen

Was:

- Expositionsvergleich nach Vulnerabilitaetsquintilen
- Expositions- und Vulnerabilitaetsvergleich nach EU-Schutzklassen
- Diagnose `BfG Zone 3 vs EU protection`

Warum:

- Vor jeder Regressionslogik sollte sichtbar sein, ob ueberhaupt ein deskriptives Disparitaetsmuster vorliegt.

### Schritt 11: Two-part-Modelle

Was:

- Logit: `any Zone 1 exposure`
- lineares Modell: `amount of Zone 1 exposure | exposed`
- zusaetzlich ein Modell fuer `depth-weighted exposure`

Warum:

- Dein Outcome ist stark null-inflationiert.
- Ein einziges lineares Modell ueber alle Gemeinden ist dafuer nicht ideal.

### Schritt 12: Raeumliche Diagnostik

Was:

- Moran's I auf OLS-Residuen
- SAR/SEM als Sensitivitaetscheck
- kNN-Fallback bei Inseln/Subgraphen

Warum:

- Gemeinde-Flaechendaten im Flussraum sind raeumlich klar abhaengig.
- Das muss diagnostisch offengelegt werden.

## GISD-Status

GISD ist in diesem Redesign bewusst nicht Teil des laufenden Hauptworkflows.

Empfohlene spaetere Rolle:

- externer Plausibilitaetscheck des INKAR-Hauptindex
- alternative Operationalisierung fuer einen Sensitivitaetsanhang
- deskriptiver Vergleich der Exposition nach Deprivationsklassen

Nicht empfohlen fuer jetzt:

- GISD parallel schon im Hauptindex zu vermischen
- GISD gleichzeitig als Bestandteil und Validierungsziel zu verwenden

## Praktische Lesart fuer die Thesis

Wenn du diese Logik in die Arbeit uebernimmst, ist der wichtigste argumentative Punkt:

- Der Vulnerabilitaetsindex soll nicht "alles Soziale" abbilden.
- Er soll nur diejenigen sozialen Merkmale abbilden, die in der Flood-Vulnerability-Literatur konsistent mit erhoehter Sensitivitaet oder verringerter Anpassungs- und Bewaeltigungskapazitaet verbunden sind.
- Darum wurde die Indikatorzahl bewusst reduziert, in Domänen geordnet und transparent aggregiert.
