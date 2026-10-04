"""
Feld-Metadaten und Settings-Basisklasse (Plan: Kernbausteine, Settings).

Jede Einstellung eines Workflows ist ein Dataclass-Feld mit ``FieldMeta``:
Anzeige (Label, Hilfe, Einheit, Ebene Basis/Erweitert/Experte), Wertebereich
(Choices, Min/Max), Abhaengigkeiten (``enabled_if``) und die Abbildung auf die
CLI (Flag, invertierte Schalter, betroffene Befehle).  Die Defaults der
Dataclass sind die CLI-Defaults.  ``to_argv``/``from_namespace`` bilden beide
Richtungen ab; ein Paritaetstest vergleicht sie mit den argparse-Parsern.
Formulare, ``run.json`` und "Befehl kopieren" nutzen dieselben Metadaten.
"""

from __future__ import annotations

import dataclasses
import functools
import math
import subprocess
import types
import typing
from dataclasses import dataclass
from typing import Any, Callable, Optional

from .issues import Issue

META_KEY = "dfm"
FIELD_LEVELS = ("basic", "advanced", "expert")


class _Same:
    """Marker: CLI-Default = Default der Dataclass."""

    def __repr__(self) -> str:
        return "SAME"


SAME = _Same()


@dataclass(frozen=True)
class FieldMeta:
    label: str
    help: str = ""
    unit: str = ""
    level: str = "basic"             # basic | advanced | expert
    choices: tuple = ()
    min: Optional[float] = None      # inklusiv
    max: Optional[float] = None      # inklusiv
    cli_flag: Optional[str] = None   # z.B. "--grid"; None = keine eigene CLI-Option
    cli_invert: bool = False         # bool-Feld; das CLI-Flag (store_true) setzt es auf False
    cli_default: Any = SAME          # CLI-Default in CLI-Werten, falls != Dataclass-Default
    tools: tuple = ()                # nur fuer diese dfm-Befehle (leer = alle des Workflows)
    kind: str = ""                   # names (Liste <-> "A,B") | spec (Mini-Sprache) | path
    enabled_if: Optional[Callable] = None   # (settings, info) -> None | Grund, warum gesperrt

    @property
    def dest(self) -> Optional[str]:
        """argparse-``dest`` des Flags."""
        return self.cli_flag.lstrip("-").replace("-", "_") if self.cli_flag else None

    def applies_to(self, tool: Optional[str]) -> bool:
        return tool is None or not self.tools or tool in self.tools


def setting(default=dataclasses.MISSING, *, default_factory=dataclasses.MISSING, **meta):
    """``dataclasses.field`` mit ``FieldMeta`` (Schluesselwoerter wie ``FieldMeta``)."""
    if meta.get("level", "basic") not in FIELD_LEVELS:
        raise ValueError(f"level muss einer von {FIELD_LEVELS} sein")
    return dataclasses.field(default=default, default_factory=default_factory,
                             metadata={META_KEY: FieldMeta(**meta)})


def _is_union(hint) -> bool:
    return typing.get_origin(hint) in (typing.Union, types.UnionType)


@functools.cache
def _type_hints(cls) -> dict:
    """``typing.get_type_hints`` je Klasse einmal (wertet die Annotationstexte jedes Mal neu aus)."""
    return typing.get_type_hints(cls)


def _coerce(name: str, value, hint):
    """JSON-Wert -> Feldtyp (int -> float erlaubt, 2.0 -> 2); sonst ``ValueError``."""
    if _is_union(hint):
        args = typing.get_args(hint)
        if value is None and type(None) in args:
            return None
        inner = [a for a in args if a is not type(None)]
        return _coerce(name, value, inner[0])
    bad = ValueError(f"Einstellung {name!r}: ungueltiger Wert {value!r}")
    if hint is bool:
        if isinstance(value, bool):
            return value
        raise bad
    if hint is int:
        if isinstance(value, bool):
            raise bad
        if isinstance(value, int):
            return value
        if isinstance(value, float) and value.is_integer():
            return int(value)
        raise bad
    if hint is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise bad
        return float(value)
    if hint is str:
        if isinstance(value, str):
            return value
        raise bad
    if typing.get_origin(hint) is list:
        if not isinstance(value, (list, tuple)):
            raise bad
        (item,) = typing.get_args(hint) or (Any,)
        return [v if item is Any else _coerce(name, v, item) for v in value]
    return value


def _cli_text(value) -> str:
    if isinstance(value, float):
        return repr(value)
    return str(value)


class SettingsBase:
    """Gemeinsame Methoden der Settings-Dataclasses (Unterklassen sind ``@dataclass``)."""

    WORKFLOW: typing.ClassVar[str] = ""

    # -- Metadaten ------------------------------------------------------------
    @classmethod
    def field_meta(cls) -> dict:
        """``{feldname: FieldMeta}`` in Definitionsreihenfolge."""
        return {f.name: f.metadata[META_KEY] for f in dataclasses.fields(cls)
                if META_KEY in f.metadata}

    @classmethod
    def _defaults(cls) -> dict:
        out = {}
        for f in dataclasses.fields(cls):
            if f.default is not dataclasses.MISSING:
                out[f.name] = f.default
            elif f.default_factory is not dataclasses.MISSING:
                out[f.name] = f.default_factory()
        return out

    # -- Dict (run.json, Job-Datei) ---------------------------------------------
    def to_dict(self) -> dict:
        return {name: (list(v) if isinstance(v, tuple) else v)
                for name, v in ((f.name, getattr(self, f.name)) for f in dataclasses.fields(self))}

    @classmethod
    def unknown_keys(cls, d: dict) -> list:
        names = {f.name for f in dataclasses.fields(cls)}
        return sorted(k for k in d if k not in names)

    @classmethod
    def from_dict(cls, d: dict, *, strict: bool = True):
        """Umkehrung von ``to_dict``; fehlende Schluessel = Default.  ``strict``:
        unbekannte Schluessel sind ein Fehler (sonst ignoriert, z.B. alte ``run.json``)."""
        unknown = cls.unknown_keys(d)
        if unknown and strict:
            raise ValueError(f"Unbekannte Einstellungen fuer {cls.__name__}: {', '.join(unknown)}")
        hints = _type_hints(cls)
        kw = {f.name: _coerce(f.name, d[f.name], hints[f.name])
              for f in dataclasses.fields(cls) if f.name in d}
        return cls(**kw)

    def replace(self, **changes):
        return dataclasses.replace(self, **changes)

    def non_default(self) -> dict:
        """Felder, die vom Default abweichen (fuer "Auf Standard zuruecksetzen")."""
        defaults = self._defaults()
        return {k: v for k, v in self.to_dict().items() if defaults.get(k) != v}

    # -- CLI ------------------------------------------------------------------
    def _to_cli(self, meta: FieldMeta, value):
        if meta.cli_invert:
            return not value
        if meta.kind == "names":
            return ",".join(value) if value else None
        return value

    def _cli_default(self, name: str, meta: FieldMeta):
        if meta.cli_default is not SAME:
            return meta.cli_default
        return self._to_cli(meta, self._defaults()[name])

    def to_argv(self, tool: Optional[str] = None) -> list:
        """CLI-Optionen der Felder, die vom CLI-Default abweichen (``--flag=wert`` bei ``-``)."""
        out = []
        for name, meta in self.field_meta().items():
            if not meta.cli_flag or not meta.applies_to(tool):
                continue
            value = self._to_cli(meta, getattr(self, name))
            if value == self._cli_default(name, meta):
                continue
            if isinstance(value, bool):
                if value:
                    out.append(meta.cli_flag)
                continue
            if value is None:
                continue
            text = _cli_text(value)
            out += [f"{meta.cli_flag}={text}"] if text.startswith("-") else [meta.cli_flag, text]
        return out

    @classmethod
    def from_namespace(cls, ns, tool: Optional[str] = None):
        """Settings aus einem argparse-``Namespace`` (Felder ohne Flag bleiben Default)."""
        defaults = cls._defaults()
        kw = {}
        for name, meta in cls.field_meta().items():
            if not meta.cli_flag or not meta.applies_to(tool) or not hasattr(ns, meta.dest):
                continue
            v = getattr(ns, meta.dest)
            if meta.cli_invert:
                v = not v
            elif meta.kind == "names":
                v = [t.strip() for t in (v or "").split(",") if t.strip()] or None
            elif meta.cli_default is not SAME and v == meta.cli_default:
                v = defaults[name]
            kw[name] = v
        return cls(**kw)

    # -- Pruefung -------------------------------------------------------------
    def validate(self, info=None) -> list:
        """Wertebereiche aller Felder, dann ``_check`` der Unterklasse (felduebergreifend)."""
        issues = []
        for name, meta in self.field_meta().items():
            v = getattr(self, name)
            if v is None:
                continue
            if meta.choices and v not in meta.choices:
                allowed = ", ".join(_cli_text(c) for c in meta.choices)
                issues.append(Issue("error", "SETTINGS.CHOICE",
                                    f"{meta.label}: {v!r} ist nicht erlaubt (erlaubt: {allowed}).",
                                    field=name))
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                if not math.isfinite(v):
                    issues.append(Issue("error", "SETTINGS.RANGE", f"{meta.label}: keine Zahl.",
                                        field=name))
                elif meta.min is not None and v < meta.min:
                    issues.append(Issue("error", "SETTINGS.RANGE",
                                        f"{meta.label}: mindestens {meta.min:g}{_unit(meta)}.",
                                        field=name))
                elif meta.max is not None and v > meta.max:
                    issues.append(Issue("error", "SETTINGS.RANGE",
                                        f"{meta.label}: hoechstens {meta.max:g}{_unit(meta)}.",
                                        field=name))
        return issues + list(self._check(info))

    def _check(self, info) -> list:
        return []

    def disabled_fields(self, info=None) -> dict:
        """``{feld: Grund}`` fuer Felder, die in dieser Lage nicht wirken (Formular: gesperrt)."""
        out = {}
        for name, meta in self.field_meta().items():
            if meta.enabled_if is not None:
                reason = meta.enabled_if(self, info)
                if reason:
                    out[name] = reason
        return out


def _unit(meta: FieldMeta) -> str:
    return f" {meta.unit}" if meta.unit else ""


def command_string(argvs: list) -> str:
    """argv-Listen als kopierbarer Text (eine Zeile je Befehl, Windows-Quoting)."""
    return "\n".join(subprocess.list2cmdline([str(a) for a in argv]) for argv in argvs)
