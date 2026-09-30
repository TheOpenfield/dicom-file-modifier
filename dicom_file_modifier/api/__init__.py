"""
api - stabile Schnittstelle zwischen dem Kern und seinen Oberflaechen (Plan P0.8).

Je Workflow ein Modul mit derselben Form:

  ``Settings``       Dataclass mit Feld-Metadaten; Defaults = CLI-Defaults
  ``inspect(sel)``   liest die Eingaben, rechnet nicht, druckt nichts
  ``preview(info, settings)``  Pruefung vor dem Start: effektive Werte, Issues
  ``command(settings, sel, out_dir)``  aequivalente ``dfm``-Befehle (argv-Listen)
  ``run(settings, sel, out_dir)``      fuehrt aus, schreibt nach ``out_dir``, liefert
                                       immer ein ``JobResult`` (wirft nicht)

  api.structures   Strukturanalyse (analyzer + visualizer)
  api.dose         Dosisindizes (dose_indices)
  api.transform    Transformation CT + RTSTRUCT oder nur CT (case_modifier, modifier)

Gemeinsame Bausteine: ``issues`` (Befunde, Ausnahmen als Klartext), ``fields``
(Feld-Metadaten, ``SettingsBase``, CLI-Abbildung), ``selection``
(``CaseSelection``), ``outputs`` (``OutputSpec``), ``results`` (``JobResult``,
Stufenzeiten, Manifest), ``sysinfo`` (freier Arbeitsspeicher).

Regeln (per Test erzwungen): Kein Kernmodul importiert ``api`` (ausser dem
Einstiegspunkt ``cli``); ``api`` importiert nie ``gui``; ``api`` laedt beim
Import und in ``inspect``/``preview`` weder pyplot noch ``visualizer`` oder
``dose_viz`` - die laufen erst in ``run``, also im Worker-Prozess.
"""

from __future__ import annotations

import importlib

WORKFLOWS = ("structures", "dose", "transform")


def workflow(name: str):
    """Workflow-Modul ``structures`` | ``dose`` | ``transform`` (erst bei Bedarf importiert)."""
    if name not in WORKFLOWS:
        raise ValueError(f"Unbekannter Workflow {name!r}; erlaubt: {', '.join(WORKFLOWS)}")
    return importlib.import_module(f"{__name__}.{name}")
