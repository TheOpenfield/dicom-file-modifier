# Umgebungswechsel Python 3.8 → 3.14 (Plan P0.2): Drift-Bericht

Stand 2026-09-30.

Verglichen wird derselbe Code in der alten und der neuen Umgebung. Es ist der Compat-Stand von P0.2, der unverändert auf beiden Umgebungen läuft. Werkzeug ist der Golden-Harness `tools/golden.py` mit 53 Szenarien:
- 43 synthetische Szenarien auf dem Demo-Fall,
- 10 Szenarien auf zwei lokalen echten Fällen, hier REAL1 und REAL2.

Der Bericht enthält nur Zahlen, keine Patientendaten.

## Ergebnis

- **Alle Kennzahlen sind bitgenau gleich.** Das betrifft:
  - Dosisindizes, Volumina, D-Werte, Abstände und Formmaße;
  - alle CT-Pixel beim Resampling in Ordnung 0, 1 und 3 sowie bei der Metadaten-Methode;
  - alle transformierten Kontur- und Positionskoordinaten.

  Außer bei den Hilfskonturen im nächsten Punkt findet der Migrationsvergleich in keiner JSON-, CSV- oder DICOM-Zahl eine Abweichung, auch nicht unterhalb der Toleranzen.
- **Eine fachliche Änderung betrifft die Hilfskonturen** im Isodosen-RTSTRUCT von `dose_indices`, also den unterdosierten Anteil des Ziels und den Anteil der Verschreibungsisodose außerhalb des Ziels.
  - GEOS 3.13 entfernt bei der Konturvereinfachung auch den Startpunkt eines Rings, wenn er innerhalb der Toleranz liegt. GEOS 3.11 hielt ihn immer fest.
  - Betroffen sind 337 Ringe mit je genau einem Punkt.
  - Der größte Hausdorff-Abstand beträgt 0,098 mm bei 0,1 mm Vereinfachungstoleranz.
  - Die Indizes werden auf den Masken berechnet und ändern sich nicht.
- **Sonst ändert sich nur die Darstellung:**
  - PNG-Abmessungen um bis zu 25 px (matplotlib 3.11);
  - `validation.html` wird etwa 0,4 MB kleiner (plotly.js 6.9);
  - in der Konsolenausgabe zwei Dateigrößen und ein Self-Test-Messwert.
- **Laufzeit:** Die CLI-Aufrufe aller 53 Szenarien laufen 20 % schneller (223 s → 178 s), die Transformation der echten Fälle 14 bis 36 %.
- **Speicher:** Der Speicherbedarf bleibt gleich. Das Resampling der echten CTs braucht in beiden Umgebungen 2,3 bis 2,4 GB physischen Speicher.
- **Determinismus:** Eine Wiederholung in jeder Umgebung ist bitgenau gleich.
- **Freigabe:** am 2026-09-30 erteilt. G3b ist die neue Referenz, danach gilt für jede Änderung an einem CLI-Modul exakte Parität.

## Umgebungen

| Paket | alt | neu |
|---|---|---|
| Python | 3.8.5 | 3.14.3 |
| numpy | 1.24.4 | 2.5.3 |
| scipy | 1.10.1 | 1.18.1 |
| pydicom | 2.4.4 | 3.0.2 |
| shapely (GEOS) | 2.0.7 (3.11.4) | 2.1.2 (3.13.1) |
| matplotlib | 3.7.5 | 3.11.2 |
| scikit-image | 0.21.0 | 0.26.0 |
| plotly | 6.6.0 | 6.9.0 |

- Die neue Umgebung ist in `uv.lock` festgeschrieben.
- Die alte steht in `tools/legacy/requirements-py38.txt`.
- Beide Läufe fanden auf demselben Windows-11-Rechner statt.

## Messstufen

| Lauf | Code | Umgebung | verglichen mit (Modus) | Ergebnis (53 Szenarien) |
|---|---|---|---|---|
| G0 | Original | alt | Referenz | – |
| G0b | Original | alt | G0 (exakt) | 53 identisch |
| G0c | Compat | alt | G0 (exakt, geplante Änderungen zugelassen) | 43 identisch, 10 geplante Änderungen |
| G1 | Compat | alt + `NPY_PROMOTION_STATE=weak` | G0c (exakt) | 53 identisch |
| G2 | Compat | alt + `PYDICOM_FUTURE=1` | G0c (exakt) | 53 identisch |
| G3 | Compat | neu | G0c (Migration) | Befunde unten |
| G0d | Compat | alt | G0c (exakt) | 52 identisch; 1 Unterschied durch die Harness-Korrektur beim VR-Enum (Abschnitt 3) |
| G3b | Compat | neu | G3 (exakt) | 53 identisch |
| | | | **G0d → G3b (Migration)** | **43 bestanden, 10 Warnungen (nur Darstellung), 0 Fehler** |

Die Läufe im Einzelnen:
- **G0b** wiederholt G0 und zeigt, dass der Harness deterministisch ist.
- **G1** stellt in der alten Umgebung die Typ-Promotion von numpy 2 ein (NEP 50).
- **G2** stellt das pydicom-3-Verhalten ein.
- Beide ändern nichts. Die Unterschiede in G3 stammen also aus den neuen Bibliotheksversionen selbst.
- **G0d und G3b** wiederholen G0c und G3 mit dem endgültigen Harness, dessen Korrekturen unten beschrieben sind. Sie liefern auch die Laufzeit- und Speicherwerte. G3b wird die neue Referenz.

Die 10 geplanten Änderungen des Compat-Commits (⚑ im Plan):
- `↔` wird in der Konsolenausgabe zu `<->` (7 Szenarien).
- Eine RLE-komprimierte CT-Eingabe wird als Explicit VR Little Endian gespeichert, statt beim Speichern abzustürzen (2).
- `Smallest/LargestImagePixelValue` entfallen nach dem Resampling (1 echter Fall).

## Befunde

### 1. Hilfskonturen: GEOS vereinfacht auch den Ring-Startpunkt

`dose.mask_to_contours` vereinfacht jeden Ring mit `LinearRing.simplify(0.1, preserve_topology=True)`. Liegt der Startpunkt des Rings innerhalb der Toleranz, entfernt ihn GEOS 3.13. GEOS 3.11 behielt ihn. Welcher Punkt Startpunkt wird, ergibt sich nur aus der Abtastreihenfolge von `find_contours` und hat keine geometrische Bedeutung.

```
Quadrat mit Startpunkt auf einer Kante, simplify(0.1):
GEOS 3.11.4: 5 Punkte, der Startpunkt bleibt
GEOS 3.13.1: 4 Punkte
```

| Kennzahl | synthetisch | echt |
|---|---|---|
| betroffene Szenarien | 9 | 1 |
| geänderte Ringe | 310 | 27 |
| entfernte Punkte je geändertem Ring | 1, immer der Startpunkt | 1, immer der Startpunkt |
| neuer Ring ist Punkt-Teilmenge des alten | alle | alle |
| größter Hausdorff-Abstand | 0,098 mm | 0,098 mm |
| größte Flächenänderung je Ring | 0,031 mm² | 0,031 mm² |

- **Betroffene ROIs:** nur `<Ziel>_minus_ISO100` und `ISO100_minus_<Ziel>`. Die Isodosen- und Schnittmengen-ROIs blieben in allen Szenarien gleich.
- **Nicht betroffen:** alle Läufe mit `--eclipse-compat`, denn sie vereinfachen nicht.
- **Flächenänderung:** 0,031 mm² entsprechen einem halben Feingitter-Pixel (0,25 × 0,25 mm).
  - Bei den kleinsten Ringen (0,125 mm²) sind das bis zu 25 % der Ringfläche.
  - Deshalb prüft der Harness bei reiner Punktentfernung den Hausdorff-Abstand und nicht die XOR-Fläche (Regel `contour` in `tools/golden_rules.json`).
  - Ringe, deren Punkte sich verschieben, müssen zusätzlich XOR ≤ 0,2 % einhalten.
- **Self-Test `rtstruct_writer`:** Das XOR-Volumen der Annulus-Konturen sinkt von 1,2650 auf 1,2644 cm³. Das Maskenvolumen beträgt 1,2650 cm³, die Toleranz des Tests 0,5 %.
- **Bewertung:** Die exportierten Hilfskonturen weichen um höchstens 0,1 mm ab, weniger als ein halbes Feingitter-Pixel. Ein überflüssiger Punkt fällt weg. Berechnete Werte, ROI-Namen, Ringzahl je Ebene und die Referenzen auf die CT-Schichten sind unverändert.

### 2. Darstellung und Konsolentext

| Ausgabe | alt | neu | Ursache |
|---|---|---|---|
| PNG-Abmessungen (8 Szenarien) | – | um bis zu 25 px verschieden | matplotlib 3.11, `bbox_inches="tight"` |
| `validation.html` (Konsole) | 10,1 MB / 6,8 MB | 9,7 MB / 6,4 MB | plotly.js 6.9 ist kleiner |
| `dose_overview.png` im Self-Test (Konsole) | 316 kB | 311 kB | PNG-Kodierung von matplotlib |

Der übrige Konsolentext aller CLIs ist identisch. Das gilt auch für die argparse-Oberfläche aller CLIs (Szenario `cli_surface`: Flags, Defaults, Choices, Typen).

### 3. Harness-Korrekturen (keine Änderung an den Ausgaben)

Der erste Vergleich G0c → G3 zeigte zunächst weitere Unterschiede. Die meisten entstanden, weil pydicom 2 und 3 dieselben Dateien verschieden einlesen. Einer entstand durch die zufällige Länge neuer UIDs. Inhaltlich unterschieden sich die Ausgaben in keinem dieser Fälle. Der Harness schreibt DICOM jetzt unabhängig von der pydicom-Version:

- **`FileMetaInformationGroupLength`** war in etwa 1 % der CT-Dateien 2 Byte kleiner.
  - pydicom 2 kürzte für zufällige UIDs einen SHA-512-Hash immer auf 64 Zeichen.
  - pydicom 3 hängt eine Zufallszahl aus `secrets.randbelow` an. Etwa jede zehnte UID ist daher kürzer, gut 1 % hat 62 Zeichen oder weniger. Erst dann ändert sich die auf gerade Länge aufgefüllte Elementlänge.
  - Die UIDs sind gültig. Die abgeleitete Gruppenlänge wird jetzt ausgeblendet.
- **Tags ohne Keyword** (private Tags) schreibt pydicom 3 als `(0009,0010)` statt `(0009, 0010)`. Der Harness schreibt Tags jetzt einheitlich.
- **Mehrdeutiger VR** (`US or SS`, z. B. `SmallestImagePixelValue`): pydicom 2.4 speichert den aufgelösten VR als Enum, `str()` ergibt unter Python 3.8 `VR.US`. Der Harness verwendet jetzt den Wert des VR.
- **Ganzzahlige Strings** (VR IS, Zählwerte in Text und CSV) werden im Migrationsmodus exakt verglichen. Vorher ließ die Standardregel ±1 in der letzten Stelle zu. Das war für `NumberOfContourPoints` zu großzügig.
- **Neu:** `compare --stats` hängt an den Bericht an:
  - alle numerischen Abweichungen, auch innerhalb der Toleranz;
  - Pixelunterschiede je Szenario;
  - Laufzeit und Speicher;
  - neue oder weggefallene stderr-Warnungen.
- **Neu:** Laufzeit und Peak-Speicher je CLI-Aufruf werden aufgezeichnet, aber nicht verglichen.
  - `python.exe` einer venv ist unter Windows nur ein Launcher. Er startet den Interpreter als eigenen Prozess.
  - Der Harness misst deshalb über ein Job-Objekt den ganzen Prozessbaum.

### 4. Warnungen

- **In den Läufen:** keine neuen Warnungen auf stderr.
- **In beiden Umgebungen gleich:** eine `UserWarning` zu `tight_layout` und eine `RuntimeWarning` („invalid value encountered in subtract“) in der Strukturvisualisierung. Beide bestanden schon vorher.
- **Voller Warnungslauf:** Alle 43 synthetischen Szenarien liefen zusätzlich mit `PYTHONWARNINGS=always::DeprecationWarning,always::FutureWarning,always::PendingDeprecationWarning`. Die Ausgaben sind identisch mit G3b.
  - Die einzige Deprecation-Quelle ist `skimage.measure.marching_cubes` in scikit-image 0.26 (`_marching_cubes_lewiner.py`, Zeilen 213 und 237): „Setting the shape on a NumPy array has been deprecated in NumPy 2.5“. Sie tritt in 8 Szenarien auf, dort, wo Oberflächen extrahiert werden (Sphärizität, Isodosen- und CT-Oberflächen).
  - Der eigene Code, pydicom 3, scipy, shapely, matplotlib und plotly lösen in diesen Szenarien keine Deprecation- oder Future-Warnung aus.

## Laufzeit und Speicher

Gemessen wurde je CLI-Aufruf auf demselben Rechner mit 24 logischen Kernen, ohne die Zeit des Harness für die Schnappschüsse. WS ist der Peak Working Set (physischer Speicher), Commit der Peak-Commit, jeweils des größten Prozesses im Prozessbaum. Aufgeführt sind die Szenarien ab 5 s Laufzeit oder 500 MB WS.

| Szenario | alt s | neu s | neu/alt | alt WS MB | neu WS MB | alt Commit MB | neu Commit MB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `cm_self_test` | 16,6 | 10,3 | 0,62 | 408 | 379 | 1902 | 1868 |
| `cm_res_verify` | 7,4 | 5,7 | 0,76 | 400 | 374 | 1892 | 1864 |
| `cm_meta_verify` | 6,0 | 4,2 | 0,69 | 400 | 318 | 1893 | 1808 |
| `viz_std` | 7,1 | 5,6 | 0,79 | 335 | 313 | 1852 | 1798 |
| REAL1 Strukturanalyse | 7,9 | 5,0 | 0,63 | 495 | 353 | 2013 | 1839 |
| REAL1 Visualisierung | 10,9 | 8,3 | 0,76 | 505 | 427 | 2022 | 1910 |
| REAL1 Dosisindizes | 3,0 | 2,8 | 0,93 | 502 | 514 | 2041 | 2017 |
| REAL1 Transformation metadata | 24,9 | 15,9 | 0,64 | 2212 | 1781 | 3708 | 3273 |
| REAL1 Transformation resample | 31,0 | 22,3 | 0,72 | 2358 | 2321 | 3854 | 3817 |
| REAL2 Visualisierung | 7,8 | 7,7 | 0,99 | 218 | 254 | 1733 | 1737 |
| REAL2 Transformation metadata | 11,5 | 8,2 | 0,72 | 1300 | 1197 | 2794 | 2687 |
| REAL2 Transformation resample | 17,1 | 14,8 | 0,86 | 2272 | 2394 | 3771 | 3890 |
| **alle 53 Szenarien** | **223** | **178** | **0,80** | | | | |

- **Commit-Grundlast von rund 1,5 GB:** In der neuen Umgebung stammt sie aus OpenBLAS. Beim Import sagen numpy und scipy (mit eigener OpenBLAS-Kopie) je Thread einen Puffer zu, je rund 30 MB pro logischem Kern.
  - `import scipy.ndimage` allein: 1518 MB Commit bei 59 MB WS.
  - Mit `OPENBLAS_NUM_THREADS=1`: 40 MB Commit.
  - Physischen Speicher belegen die Puffer kaum. Sie zählen aber gegen die Commit-Grenze aus RAM und Auslagerungsdatei.
- **Resampling:** Beide echten CTs brauchen in beiden Umgebungen 2,3 bis 2,4 GB WS (Abweichung ±5 %).
- **Metadaten-Methode:** Sie lädt heute ebenfalls alle Pixel und kommt so auf 1,2 bis 2,2 GB WS, obwohl sie die Pixel nicht verändert.

## Hinweise für die nächsten Schritte

- **Worker (P0.8):** Er setzt `OPENBLAS_NUM_THREADS`, wie im Plan vorgesehen. Die Messung liefert den Grund: Ohne Begrenzung sagt jeder Prozess bei 24 Kernen rund 1,5 GB zu. Die zeitkritischen Teile (`map_coordinates`, `cKDTree`) nutzen kein BLAS. Welcher Thread-Wert sinnvoll ist, wird im Worker-Schritt gemessen.
- **RAM-Prüfung vor dem Resampling:** Die Messwerte (2,3 bis 2,4 GB WS) passen zur Schätzung im Plan (2 bis 2,8 GB).
- **Metadaten-Methode:** Sie muss keine Pixel dekodieren. Hier können P0.5 und P0.6 Speicher sparen.
- **Abhängigkeits-Updates:** Die Deprecation in `marching_cubes` wird mit einer späteren numpy-Version zum Fehler, falls scikit-image sie bis dahin nicht behebt. Die Self-Tests von analyzer und dose_indices rufen `marching_cubes` auf und fangen das ab.

## Bewertung

- Die Migration ändert keine berechnete Kennzahl und keinen CT-Pixelwert.
- Die einzige inhaltliche Änderung betrifft die Geometrie exportierter Hilfskonturen. Sie ist kleiner als die bewusst gewählte Vereinfachungstoleranz und entfernt nur einen überflüssigen Punkt.
- Alle weiteren Unterschiede betreffen die Darstellung oder den Harness.
- **Freigegeben am 2026-09-30:**
  - G3b ist die Referenz.
  - Die P0.2-Einträge in `expected_changes` sind entfernt.
  - Ab P0.3 muss `compare G3b <neu> --mode exact` ohne FAIL bleiben.
- **Goldens werden nicht committet:**
  - Die synthetischen Schnappschüsse sind 47,5 MB groß.
  - Exakte Parität ist nur auf derselben CPU verlässlich.
  - Die Referenz bleibt lokal. Die CI (Phase 1) soll Basis- und PR-Commit auf demselben Runner rechnen und exakt vergleichen.

## Reproduktion

```bash
python tools/golden.py make-inputs
.venv/Scripts/python.exe    tools/golden.py run --label G0d --real data/<case-id> --real data/<case-id>
.venv314/Scripts/python.exe tools/golden.py run --label G3b --real data/<case-id> --real data/<case-id>
python tools/golden.py compare G0d G3b --mode migration --stats
```

Die Berichte liegen im Arbeitsordner `%USERPROFILE%\dfm-golden\reports`. Sie bleiben lokal, weil die echten Fälle Patientendaten enthalten.
