"""
gui - Desktop-App (PySide6): ``dfm gui [FALL]`` bzw. ``dicom-rt-toolkit``.

  app         Einstieg (QApplication, Hauptfenster)
  config      App-Einstellungen (Ergebnis-Stammordner, Job-Ordner)
  jobs        ``JobRunner`` (Worker-Prozess, JSON-Zeilen -> Signale), Hintergrund-Aufgaben
  widgets     Einstellungsformular aus ``FieldMeta``, Befundliste, Bildergalerie
  window      Hauptfenster: Fallzeile, Seitenleiste, Fortschritt, Protokoll
  structures_page  Seite Strukturanalyse

Regeln: Die GUI nutzt nur ``api``.  Berechnungen laufen im Worker-Prozess
(``dfm worker``); der GUI-Prozess laedt nie pyplot, ``visualizer``,
``dose_viz`` oder plotly (per Test erzwungen).
"""
