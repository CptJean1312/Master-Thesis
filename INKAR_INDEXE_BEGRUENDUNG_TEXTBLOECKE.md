# INKAR-Indikatoren fuer den Flood Vulnerability Index

## Kurzbewertung

Die aktuelle Auswahl ist insgesamt gut begruendbar. Der staerkste Punkt ist die jetzt klare Trennung zwischen:

- Kernindikatoren fuer den Hauptindex
- optionalen Robustheitsindikatoren
- bewusst ausgeschlossenen Strukturvariablen

Das ist methodisch sauberer als ein breiter Sammelindex, weil dadurch soziale Vulnerabilitaet enger an Sensitivitaet, materielle Ressourcen und alltagsbezogene Anpassungskapazitaet gebunden wird. Der einzige Punkt, den ich weiterhin vorsichtig formulieren wuerde: Erreichbarkeitsindikatoren koennen teilweise auch laendliche Lage oder Versorgungsstruktur abbilden. Genau deshalb ist es sinnvoll, `supermarket`, `doctors` und `broadband` nur als Sensitivitaetsvariablen zu fuehren und nicht im Kernindex zu verankern.

## Literaturbasierter Rahmen

Die Auswahl folgt einem Vulnerabilitaetsverstaendnis, in dem vor allem drei Aspekte relevant sind:

- materielle und soziale Benachteiligung
- demografische Sensitivitaet
- eingeschraenkte Bewaeltigungs- und Anpassungskapazitaet

Diese Logik ist anschlussfaehig an:

- IPCC AR6 WGII Chapter 1
- Cutter, Boruff and Shirley (2003)
- Rufat et al. (2015)
- Fekete (2009)
- Hinkel (2011)
- Tate (2012)

## 1. Kernindikatoren des Hauptindex

### `exposure_share_alg2_sgb2`

Diese Variable misst den Anteil der Bevoelkerung, der Leistungen nach ALG II bzw. SGB II bezieht. Sie ist ein direkter Indikator fuer materielle Benachteiligung und eingeschraenkte finanzielle Spielraeume. In der Literatur gilt soziooekonomischer Status als einer der robustesten Treiber sozialer Hochwasservulnerabilitaet, weil geringe Ressourcen Vorsorge, Evakuierung und Erholung erschweren. Die Variable wird ohne Transformation in den Index aufgenommen.

### `exposure_share_sgb2_with_housing_costs`

Diese Variable bildet den Anteil von SGB-II-Empfaengerinnen und -Empfaengern mit Wohnkostenunterstuetzung ab. Sie ergaenzt den allgemeinen Transferbezug um eine wohnungsbezogene Deprivationsdimension und weist auf geringe finanzielle Puffer im Haushalt hin. Gerade nach Flutereignissen wirken sich niedrige Ruecklagen und hohe Wohnkostenbelastung negativ auf die Recovery-Kapazitaet aus. Die Variable wird ohne Transformation verwendet.

### `exposure_share_longterm_unemp`

Langzeitarbeitslosigkeit steht fuer verfestigte Arbeitsmarktbenachteiligung. Im Unterschied zu kurzfristigen Schwankungen verweist sie auf stabile strukturelle Ausschlusslagen und damit auf eingeschraenkte materielle und soziale Resilienz. Die Variable ist deshalb ein plausibler Kernindikator fuer verringerte Bewaeltigungs- und Wiederaufbaukapazitaet. Sie wird ohne Transformation aufgenommen.

### `exposure_unemp_u25_per_1000`

Jugendarbeitslosigkeit pro 1000 Einwohner bildet strukturelle Fragilitaet des lokalen Arbeitsmarkts ab. Sie ist kein klassischer Katastrophenindikator, kann aber auf breitere soziale Problemlagen und geringe Integrationskapazitaet hinweisen. In Kombination mit den anderen Deprivationsindikatoren ist sie vertretbar, weil sie nicht Wohlstand, sondern prekare Zukunfts- und Erwerbslagen adressiert. Die Variable wird ohne Transformation genutzt.

### `exposure_purchasing_power`

Die Kaufkraft pro Kopf misst die verfuegbaren oekonomischen Ressourcen im lokalen Kontext. Hohe Kaufkraft spricht tendenziell fuer bessere Vorsorge-, Anpassungs- und Erholungsmoeglichkeiten; fuer den Vulnerabilitaetsindex muss die Richtung deshalb umgekehrt werden. Im Code wird die Variable invertiert, so dass niedrigere Kaufkraft hoeherer Vulnerabilitaet entspricht.

### `exposure_share_hh_income_low`

Diese Variable erfasst den Anteil einkommensschwacher Haushalte. Sie ist inhaltlich sehr nah an klassischer Deprivations- und Verwundbarkeitsforschung, weil niedrige Einkommen oft mit geringeren Versicherungsoptionen, schlechteren Reserven und erhoehten Erholungsbarrieren verbunden sind. Sie wird ohne Transformation in den Hauptindex aufgenommen.

### `exposure_share_age_65plus`

Der Anteil der Bevoelkerung ab 65 Jahren bildet eine zentrale demografische Sensitivitaetsdimension ab. Aeltere Menschen sind bei Flutereignissen haeufig staerker betroffen, unter anderem wegen eingeschraenkter Mobilitaet, gesundheitlicher Vorbelastungen und hoeherer Abhaengigkeit von Hilfe im Warn- und Evakuierungsprozess. Die Variable wird ohne Transformation genutzt.

### `exposure_share_age_75plus`

Diese Variable schaerft den Altersaspekt weiter, weil sehr hohe Altergruppen in der Regel noch staerker mit Pflegebedarf, eingeschraenkter Selbststaendigkeit und erhoehter gesundheitlicher Sensitivitaet verbunden sind. Sie ist deshalb als zusaetzlicher Hochaltrigkeitsindikator plausibel. Die Variable wird ohne Transformation verwendet.

### `exposure_old_age_dependency`

Die Old-Age-Dependency-Ratio beschreibt das Verhaeltnis aelterer zu erwerbsfaehigen Bevoelkerungsgruppen. Sie bildet nicht nur Alterung, sondern auch eine potenziell geringere Unterstuetzungs- und Tragekapazitaet auf Gemeindeebene ab. Damit ergaenzt sie die Altersanteile sinnvoll um eine strukturelle Perspektive. Die Variable wird ohne Transformation verwendet.

### `exposure_share_bg_single_parent`

Der Anteil von Einelternhaushalten ist ein sinnvoller Indikator fuer potenziell eingeschraenkte haushaltsinterne Unterstuetzungsressourcen. Einelternhaushalte sind haeufig staerker mit Zeitdruck, finanzieller Belastung und organisatorischen Engpaessen konfrontiert, was in Krisensituationen relevant werden kann. Die Variable wird ohne Transformation genutzt.

### `exposure_share_single_households`

Einpersonenhaushalte koennen ueber geringere private Unterstuetzungsnetzwerke und weniger interne Redundanz im Alltag verfuegen. Das ist in Warnung, Evakuierung und Recovery potenziell relevant. Die Variable ergaenzt deshalb die Einkommens- und Altersdimension um eine soziale Organisationsperspektive. Sie wird ohne Transformation aufgenommen.

### `exposure_dist_public_transport_m`

Die Distanz zur naechsten OePNV-Haltestelle dient als Proxy fuer alltagsbezogene Erreichbarkeit und Mobilitaetsoptionen. Schlechtere Erreichbarkeit kann Evakuierung, Zugang zu Dienstleistungen und die Bewaeltigung von Alltagsunterbrechungen erschweren. Gleichzeitig bildet die Variable nicht nur Vulnerabilitaet, sondern teilweise auch laendliche Lage ab; sie ist im Hauptindex dennoch vertretbar, solange dieser Punkt offen benannt wird. Die Variable wird ohne Transformation verwendet.

### `exposure_dist_gp_m`

Die Distanz zur allgemeinmedizinischen Versorgung ist ein plausibler Indikator fuer den Zugang zu zentraler Gesundheitsinfrastruktur. Gerade fuer aeltere oder gesundheitlich belastete Bevoelkerungsgruppen ist dies fuer Coping und Recovery relevant. Die Variable wird ohne Transformation in den Index aufgenommen.

### `exposure_dist_pharmacy_m`

Die Distanz zur Apotheke erfasst den Zugang zu Medikamenten und alltagsnaher Gesundheitsversorgung. Bei Flutereignissen und Nachsorgephasen kann dies fuer vulnerable Gruppen besonders bedeutsam sein. Die Variable passt deshalb gut in die Dimension adaptive Kapazitaet und wird ohne Transformation verwendet.

## 2. Optionale Variablen fuer Sensitivitaet und Robustheit

### `exposure_doctors_total_per_1000`

Die Arztzahl pro 1000 Einwohner kann als grober Gesundheits- und Versorgungskapazitaetsindikator gelesen werden. Konzeptionell ist die Variable plausibel, sie ist aber staerker system- und strukturbezogen als die Distanzindikatoren. Zudem kann sie regionale Versorgungsorganisation und Urbanitaet mitabbilden. Deshalb ist sie als Sensitivitaetsvariable besser aufgehoben als im Kernindex. Falls sie genutzt wird, sollte sie invertiert werden, so dass weniger Aertzte hoeherer Vulnerabilitaet entsprechen.

### `exposure_dist_supermarket_m`

Die Distanz zum Supermarkt bildet Grundversorgung und alltagspraktische Resilienz ab. Der Indikator ist nachvollziehbar, aber konzeptionell etwas indirekter als Gesundheit oder Mobilitaet. Er kann zudem stark mit laendlicher Lage zusammenhaengen. Als Robustheitscheck ist er sinnvoll; fuer den Hauptindex ist die vorsichtige Einordnung als optional angemessen. Die Variable benoetigt keine Transformation.

### `exposure_share_bb_100mbit`

Breitbandzugang kann Informationszugang, digitale Kommunikation und administrative Handlungsfaehigkeit beeinflussen. Gleichzeitig ist die theoretische Verbindung zur Flood Vulnerability weniger direkt als bei Armut, Alter oder Gesundheitszugang. Deshalb ist die Variable eher als moderne Kontext- oder Sensitivitaetsvariable geeignet. Falls sie verwendet wird, sollte sie invertiert werden, damit geringe digitale Abdeckung hoehere Vulnerabilitaet bedeutet.

## 3. Bewusst ausgeschlossene Variablen

### `exposure_pop_density_per_km2`

Bevoelkerungsdichte beschreibt vor allem Siedlungsstruktur. Hohe oder niedrige Dichte kann je nach Infrastruktur, Versorgung und sozialer Zusammensetzung sehr unterschiedliche Bedeutungen haben. Die Variable sollte deshalb eher als Kontrollvariable oder Kontextmerkmal behandelt werden, nicht als direkter Vulnerabilitaetsindikator.

### `exposure_employment_density_per_km2`

Beschaeftigungsdichte sagt primaer etwas ueber regionale Wirtschafts- und Raumstruktur aus. Der Bezug zu household-level vulnerability ist indirekt und inhaltlich unscharf. Als Hauptindexvariable ist sie deshalb nicht gut geeignet.

### `exposure_trade_tax_per_capita`

Gewerbesteuer pro Kopf bildet kommunale Ertragsstruktur ab, nicht direkt die Lebenslage der Bevoelkerung. Dadurch verschiebt sich die Analyse von sozialer Vulnerabilitaet hin zu kommunaler Finanz- und Wirtschaftsstruktur. Fuer den Hauptindex ist die Variable daher nicht passend.

### `exposure_tax_revenue_total`

Gesamte Steuerertraege sind stark von Groesse und Wirtschaftsstruktur der Gemeinde abhaengig. Die Variable ist kein sauberer Indikator fuer soziale Verwundbarkeit von Einwohnerinnen und Einwohnern und wurde deshalb ausgeschlossen.

### `exposure_share_hh_income_medium`

Der Anteil mittlerer Einkommen fuegt gegenueber dem Anteil niedriger Einkommen kaum eigenstaendige Information hinzu. Er ist vor allem als Komplementstruktur des Einkommensprofils lesbar und wurde daher aus Redundanzgruenden ausgeschlossen.

### `exposure_share_hh_income_high`

Der Anteil hoher Einkommen spiegelt weitgehend die inverse Seite von Einkommensschwaeche wider. Fuer einen fokussierten Vulnerabilitaetsindex liefert er daher nur begrenzten Zusatznutzen und wurde nicht aufgenommen.

### `exposure_students_total_per_1000`

Studierende sind fuer das Vulnerabilitaetskonstrukt schwer eindeutig zu interpretieren. Sie koennen mobil, ressourcenarm, ressourcenstark oder nur temporar vor Ort sein. Wegen dieser konzeptionellen Ambivalenz ist die Variable fuer den Kernindex ungeeignet.

### `exposure_students_18_25_per_1000`

Auch diese Variable bleibt in der Wirkungsrichtung unscharf. Sie bildet eher Bildungs- und Altersstruktur als konsistente Flood Vulnerability ab und wurde deshalb ausgeschlossen.

### `exposure_students_fh_per_1000`

Die Fachhochschul-Studierendenquote ist noch staerker kontextabhaengig und fuer das zentrale Konstrukt nicht robust genug begruendbar. Sie bleibt daher ausserhalb des Index.

### `exposure_migration_balance`

Migrationssalden zeigen raeumliche Dynamik, aber keine klare Richtung sozialer Vulnerabilitaet. Je nach Kontext kann Zu- oder Abwanderung sehr unterschiedlich zu interpretieren sein. Die Variable ist deshalb als Kernindikator ungeeignet.

### `exposure_natural_pop_change`

Natuerliche Bevoelkerungsentwicklung ist demografisch interessant, sagt aber nicht direkt etwas ueber materielle, soziale oder adaptive Vulnerabilitaet aus. Auch hier ist der Bezug zu Flood Vulnerability zu indirekt.

## Empfohlene methodische Formulierung

Fuer den Hauptindex sollten nur die Kernindikatoren verwendet werden. Die drei optionalen Variablen koennen in einem Sensitivitaetslauf ergaenzt werden, um zu pruefen, ob sich die raeumlichen Muster oder Modellresultate substanziell veraendern. Die ausgeschlossenen Variablen sollten, wenn ueberhaupt, als Kontext- oder Kontrollvariablen behandelt werden, nicht als Bestandteil des Vulnerabilitaetsindex.
