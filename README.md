# enzyme-function-benchmarks

Experimente zur Frage, wie gut sich die Funktion eines Enzyms (seine EC-Nummer) vorhersagen lässt:
mit Homologie-Suche, gelernten Integratoren, allgemeinen Entscheidungsmodellen (TypeSafe Jev,
Laya) und Protein-Sprachmodellen (ESM-2, ProtT5). Beide Experimente sind abgeschlossen.

**Gesamtfazit: [`FINDINGS.md`](FINDINGS.md)**

## Ergebnis in Kürze

* **Die Leistung kommt aus der Homologie-Suche.** Jev und Laya erreichen in keinem Versuch mehr als
  ein einfaches Verfahren mit denselben Informationen. Auf der reinen Sequenz liegen sie auf
  Zufallsniveau, mit Homologie-Belegen gleichauf mit dem Nächsten Nachbarn.
* **Protein-Sprachmodelle lesen Funktion aus der Sequenz:** Eine ESM-2-Probe erreicht 51,9 % bei
  der EC-Klasse, ohne einen einzigen erkennbaren Homolog im Training.
* **Einziger signifikanter Gewinn über die Homologie-Suche:** ProtT5-Embedding-Nachbarn als Rückfall
  für Proteine ohne MMseqs2-Treffer, 64,4 % statt 62,5 % bei der exakten EC (p = 0,031).
* **Verlässliche Konfidenz:** LightGBM über die Homologie-Belege halbiert grob die Fläche unter der
  Risiko-Abdeckungs-Kurve (0,092 gegen 0,173) und sagt damit besser, wann man einer Vorhersage
  trauen kann.

## Inhalt

| Ordner | Frage | Daten |
|---|---|---|
| [`enzyme-evidence/`](enzyme-evidence/) | Kann ein gelerntes Modell Homologie-Belege besser gewichten als der Nächste Nachbar? | CARE-Benchmark, Task 1 (1.140 Testproteine) |
| [`enzyme-direct/`](enzyme-direct/) | Können Jev und Laya die Funktion aus der Sequenz ablesen, mit Kontext, mit Homologen? Wie schneiden Protein-Sprachmodelle ab? | Eigener Benchmark aus Swiss-Prot seit 2018, getrennt von den CARE-Testsets (320 Proteine) |

Die beiden Ordner teilen keinen Code und keine Daten. Jeder hat seine eigene README mit Aufbau,
allen Tabellen, Einschränkungen und Reproduktion (`requirements.txt`, nummerierte Skripte,
`pytest -q`). Jev braucht `TYPESAFE_API_KEY`; alles andere läuft auf der CPU.

## Methodik

* Analysepläne und Pass-Kriterien sind vor jedem Lauf in Git festgeschrieben; nachträgliche
  Auswertungen sind als solche markiert.
* Jede Stufe hat eine Kontrolle mit genau derselben Information; Vergleiche nutzen gepaarte Tests.
* Leckage-Tests sichern die Trennung von Trainings- und Testdaten; Benchmarks sind als Dateien
  eingefroren.
