"""Bausteine der Seiten: Einstellungsformular aus ``FieldMeta``, Befundliste, Bildergalerie."""

from __future__ import annotations

import os
import types
import typing
from pathlib import Path

from PySide6.QtCore import QSize, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QIcon, QPixmap
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QLineEdit,
                               QListView, QListWidget, QListWidgetItem, QSpinBox, QStyle, QWidget)

INT_MAX = 2**31 - 1
ZWSP = chr(0x200B)                      # Leerzeichen der Breite null: erlaubt Umbruch in Pfaden
STATUS_DE = {"ok": "Fertig", "ok_warnings": "Fertig, mit Hinweisen",
             "failed": "Fehlgeschlagen", "cancelled": "Abgebrochen"}


def breakable(path) -> str:
    """Pfad fuer QLabel mit Zeilenumbruch: nach jedem Trenner darf umbrochen werden."""
    return str(path).replace(os.sep, os.sep + ZWSP).replace("/", "/" + ZWSP)


def open_path(path) -> None:
    """Datei oder Ordner mit dem Standardprogramm von Windows oeffnen."""
    QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))


def _unwrap(hint) -> tuple:
    """``(Typ, optional)`` aus einem Typ-Hinweis wie ``Optional[float]``."""
    args = typing.get_args(hint)
    if typing.get_origin(hint) in (typing.Union, types.UnionType) and type(None) in args:
        rest = [a for a in args if a is not type(None)]
        return (rest[0] if len(rest) == 1 else str), True
    return hint, False


class SettingsForm(QWidget):
    """
    Formular fuer die Felder einer ``Settings``-Klasse: Label, Einheit, Hilfe
    (Tooltip), Bereich und Choices aus ``FieldMeta``.  ``settings()`` liefert
    die Werte als ``Settings`` (``ValueError`` bei ungueltiger Eingabe).
    """

    changed = Signal()

    def __init__(self, settings_cls, fields=None, parent=None):
        super().__init__(parent)
        self.cls = settings_cls
        metas = settings_cls.field_meta()
        hints = typing.get_type_hints(settings_cls)
        self._rows: dict = {}
        layout = QFormLayout(self)
        for name in fields or [n for n, m in metas.items() if m.level == "basic"]:
            meta = metas[name]
            typ, optional = _unwrap(hints[name])
            kind, w = self._widget(meta, typ, optional)
            w.setToolTip(meta.help)
            self._rows[name] = (kind, w, typ, optional, meta)
            label = meta.label + (f" [{meta.unit}]" if meta.unit else "")
            layout.addRow("" if kind == "bool" else label, w)
        self._base = settings_cls()
        self.set_settings(self._base)

    def _widget(self, meta, typ, optional) -> tuple:
        emit = lambda *_: self.changed.emit()   # noqa: E731
        if meta.choices:
            w = QComboBox()
            for c in meta.choices:
                w.addItem(str(c), c)
            w.currentIndexChanged.connect(emit)
            return "choice", w
        if typ is bool:
            w = QCheckBox(meta.label)
            w.toggled.connect(emit)
            return "bool", w
        if typ in (int, float) and not optional:
            w = QSpinBox() if typ is int else QDoubleSpinBox()
            lo = meta.min if meta.min is not None else -INT_MAX
            hi = meta.max if meta.max is not None else INT_MAX
            if typ is float:
                w.setDecimals(3)
                w.setRange(float(lo), float(hi))
            else:
                w.setRange(int(max(lo, -INT_MAX)), int(min(hi, INT_MAX)))
            w.valueChanged.connect(emit)
            return "number", w
        w = QLineEdit()
        if optional:
            w.setPlaceholderText("automatisch")
        w.textChanged.connect(emit)
        return ("names" if meta.kind == "names" else "text"), w

    def set_placeholder(self, name: str, text: str) -> None:
        w = self._rows[name][1]
        if isinstance(w, QLineEdit):
            w.setPlaceholderText(text)

    def set_settings(self, s) -> None:
        for name, (kind, w, _typ, _opt, _meta) in self._rows.items():
            v = getattr(s, name)
            w.blockSignals(True)
            if kind == "choice":
                w.setCurrentIndex(max(0, w.findData(v)))
            elif kind == "bool":
                w.setChecked(bool(v))
            elif kind == "number":
                w.setValue(v)
            elif kind == "names":
                w.setText(", ".join(v or []))
            else:
                w.setText("" if v is None else str(v))
            w.blockSignals(False)
        self._base = s
        self.changed.emit()

    def values(self) -> dict:
        out = {}
        for name, (kind, w, typ, optional, meta) in self._rows.items():
            if kind == "choice":
                out[name] = w.currentData()
            elif kind == "bool":
                out[name] = w.isChecked()
            elif kind == "number":
                out[name] = w.value()
            elif kind == "names":
                out[name] = [t.strip() for t in w.text().split(",") if t.strip()] or None
            else:
                text = w.text().strip()
                if not text:
                    out[name] = None if optional else ""
                elif typ in (int, float):
                    try:
                        out[name] = int(text) if typ is int else float(text.replace(",", "."))
                    except ValueError:
                        raise ValueError(f"{meta.label}: keine gültige Zahl ({text!r})") from None
                else:
                    out[name] = text
        return out

    def settings(self):
        return self.cls.from_dict({**self._base.to_dict(), **self.values()})


_LEVEL_ICONS = {"error": QStyle.StandardPixmap.SP_MessageBoxCritical,
                "warning": QStyle.StandardPixmap.SP_MessageBoxWarning,
                "info": QStyle.StandardPixmap.SP_MessageBoxInformation}


class IssueList(QListWidget):
    """Befunde (``Issue`` oder dict) mit Symbol je Stufe; Hinweis in der zweiten Zeile, Details im Tooltip."""

    def set_issues(self, issues) -> None:
        self.clear()
        for issue in issues:
            d = issue if isinstance(issue, dict) else issue.to_dict()
            text = d["message_de"] + (f"\n→ {d['hint_de']}" if d.get("hint_de") else "")
            item = QListWidgetItem(self.style().standardIcon(_LEVEL_ICONS[d["level"]]), text)
            item.setToolTip(d.get("detail") or d.get("code", ""))
            self.addItem(item)
        self.setVisible(bool(issues))


class Gallery(QListWidget):
    """Vorschaubilder; Doppelklick oeffnet das Bild in der Bildanzeige von Windows."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setViewMode(QListView.ViewMode.IconMode)
        self.setIconSize(QSize(240, 170))
        self.setResizeMode(QListView.ResizeMode.Adjust)
        self.setMovement(QListView.Movement.Static)
        self.setSpacing(8)
        self.setWordWrap(True)
        self.itemActivated.connect(lambda item: open_path(item.data(Qt.ItemDataRole.UserRole)))

    def set_images(self, paths) -> None:
        self.clear()
        for p in paths:
            pix = QPixmap(str(p))
            if pix.isNull():
                continue
            thumb = pix.scaled(self.iconSize(), Qt.AspectRatioMode.KeepAspectRatio,
                               Qt.TransformationMode.SmoothTransformation)
            item = QListWidgetItem(QIcon(thumb), Path(p).stem)
            item.setData(Qt.ItemDataRole.UserRole, str(p))
            item.setToolTip(f"{Path(p).name} - Doppelklick zum Öffnen")
            self.addItem(item)
