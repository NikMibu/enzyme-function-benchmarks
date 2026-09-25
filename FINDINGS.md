# Enzymfunktion vorhersagen: Gesamtfazit

Zwei abgeschlossene Experimente, 25. September 2026:

* [`enzyme-evidence/`](enzyme-evidence/): Homologie-Belege auf dem CARE-Benchmark (Task 1) auswerten
* [`enzyme-direct/`](enzyme-direct/): Funktion aus der Sequenz, mit berechnetem Kontext, mit Homologen und mit Protein-Embeddings

Details, Tabellen und Reproduktion stehen in den READMEs der beiden Ordner. Alle Zahlen hier
stammen aus deren `results/`-Dateien.

## Kernaussage

**Die Leistung kommt aus der Homologie-Suche, nicht aus dem Entscheidungsmodell.** Jev (TypeSafe,
`jev-1.13.0`) und Laya (Open-Weight) erreichen in keinem Versuch mehr als ein einfaches Verfahren
mit denselben Informationen. Der einzige signifikante Gewinn gegenüber reiner Homologie-Suche kommt
von einem Protein-Sprachmodell als Rückfall für Proteine ohne Treffer.

## enzyme-evidence: Homologe bewerten

CARE-Testsets (1.140 Proteine), exakte EC-Nummer, MMseqs2-Suche gegen CARE-train.

| Methode | Genauigkeit, gepoolt | Selektive Vorhersage (AURC, gepoolt, kleiner = besser) |
|---|---|---|
| Nächster Nachbar (MMseqs2) | 69,5 % | 0,173 |
| LightGBM über 19 Beleg-Merkmale | 69,8 % | **0,092** |
| Jev, zero-shot | 69,1 % | 0,221 |

* Keine Methode schlägt den Nächsten Nachbarn bei der Genauigkeit, obwohl etwa 8 Punkte
  Spielraum bestehen.
* Der einzige echte Gewinn ist LightGBMs Konfidenz: Auf dem schwersten Split (<30 % Identität)
  ist es bei 50 % beantworteten Proteinen zu 91,7 % richtig, gegen 80,6 % mit Identität ×
  Abdeckung. Jevs Konfidenz ist auf diesem Split schlechter als Identität × Abdeckung (AURC 0,269
  gegen 0,168).
* Die MMseqs2-Suche selbst übertrifft CAREs veröffentlichte Vergleichsverfahren (62,5 % gegen
  55,1 % für CLEAN auf dem <30 %-Split).

## enzyme-direct: Funktion ohne Annotation

Neuer Benchmark aus Swiss-Prot-Einträgen seit 2018, getrennt von den CARE-Testsets, eine
Proteinfamilie pro Eintrag: ec1 (6 EC-Klassen × 35, Zufall 16,7 %) und ec4 (11 exakte ECs × 10).

| Was das Verfahren bekommt | Jev | Beste Vergleichsmethode |
|---|---|---|
| nur Sequenz | 11,9 % (Laya 15,7–16,2 %) | ESM-2-Probe: **51,9 %** |
| + Eigenschaften & Zusammensetzung | 14,3 % | Logistische Regression: 29,0 % |
| + Motive | 27,1 % | Logistische Regression: 32,9 % |
| + Homologie-Belege | 87,1 % | Nächster Nachbar: 86,7 % |

EC-Klasse, ec1 (210 Proteine).

* **Rohsequenz:** Jev und Laya liegen auf Zufallsniveau und geben fast immer dieselbe Klasse aus.
* **Kontext:** Jev nutzt nur Hinweise, die eine Chemie beim Namen nennen (Motive), und bleibt auch
  dann unter der logistischen Regression.
* **Homologe:** Jev weicht nur bei 19 von 210 Proteinen vom Nächsten Nachbarn ab und verliert dabei
  bei der exakten EC 6 zu 1.
* **Protein-Sprachmodelle lesen Funktion aus der Sequenz:** ESM-2 (650M) mit linearer Probe
  erreicht 51,9 % ohne einen einzigen erkennbaren Homolog im Training.
* **Embedding-Rückfall:** Wenn MMseqs2 keinen Treffer findet, nimmt man den ProtT5-Nachbarn aus
  UniProts fertigen Einbettungen. Auf ec1 + ec4 (320 Proteine):

  | | MMseqs2-Nachbar | Mit Embedding-Rückfall |
  |---|---|---|
  | Exakte EC | 62,5 % | **64,4 %** (6 gewonnen, 0 verloren, p = 0,031) |
  | EC-Klasse | 80,3 % | **87,5 %** (p = 2×10⁻⁷) |

## Was daraus folgt

1. **Jev ist ein guter Leser, kein Denker.** Mit lesbar aufbereiteten Belegen wählt es die
   offensichtliche Antwort. Eigene Schlüsse jenseits der Hinweise zieht es nicht, und seine
   Wahrscheinlichkeiten sind schlecht kalibriert. Es ist für Textentscheidungen gebaut, nicht für
   Proteine.
2. **Das Werkzeug muss zur Datenart passen.** Für Verwandtschaft gibt es die Homologie-Suche, für
   Sequenzen ohne Verwandte Protein-Sprachmodelle. Ein Textmodell auf Aminosäureketten ist das
   falsche Werkzeug. Laya-Fine-Tuning wurde deshalb nicht mehr verfolgt.
3. **Mehrwert gibt es in zwei Nischen:** verlässliche Konfidenz (LightGBM) und Proteine ohne
   Homolog (Embedding-Rückfall). Beides sind kleine, gezielte Verbesserungen auf der
   Homologie-Suche.

## Methodik

**Stärken**

* Pass-Kriterien und Analysepläne vor jedem Lauf in Git festgeschrieben (enzyme-direct: Commits
  `1460eab`, `35e770b`, `987ac48`, `16e8082`, `ed4b12f`). Nachträgliche Analysen sind als solche
  markiert.
* Zu jeder Jev-Stufe eine Kontrolle mit genau derselben Information, gepaarte Tests auf denselben
  Proteinen.
* Leckage-Tests: kein Überlapp mit CARE-Testsets, keine Accessions, Namen oder Organismen in
  Prompts. Benchmarks als Dateien eingefroren.

**Grenzen**

* **Kleine Stichproben:** 210 bzw. 320 Proteine, 37 ohne Treffer. Der Embedding-Gewinn beruht auf
  6 Proteinen und sollte auf neuen Daten bestätigt werden.
* **Designänderungen vor den ersten Modellläufen:** CARE-train statt ganz CARE ausgeschlossen,
  35 statt 40 pro Klasse, Stichtag 2018 statt 2020. Das war nötig, weil sonst zu wenige EC-6-Familien
  übrig blieben, und geschah vor jedem Modelllauf.
* **Prompts:** Jev lief mit einem festen Prompt pro Stufe, bewusst ohne Tuning auf Testdaten.
  Bessere Prompts sind denkbar, würden die großen Abstände aber kaum schließen.
* **Schwache Kontrollen:** Die erste LogReg-Kontrolle hatte nur 4 Ligasen im Training. „Jev schlägt
  die Kontrolle nicht“ ist damit vorsichtig formuliert.
* **Embeddings:** ProtT5 wurde auf UniRef50 trainiert und kennt die Sequenzen womöglich, allerdings
  ohne Funktionsangaben.

## Offene Fragen

* Embedding-Rückfall auf einem größeren Benchmark mit vielen Proteinen ohne Treffer bestätigen.
* Kosinus-Ähnlichkeit als Konfidenz für Rückfall-Antworten (Median 0,72 ohne Treffer, 0,96 mit).
* LightGBMs Konfidenz aus enzyme-evidence auf dem neuen Benchmark prüfen.
