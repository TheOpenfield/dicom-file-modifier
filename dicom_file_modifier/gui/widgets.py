"""Bausteine der Seiten: Einstellungsformular aus ``FieldMeta``, Befundliste, Bildergalerie."""

from __future__ import annotations

import os
import types
import typing
from pathlib import Path

from PySide6.QtCore import QSize, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFontDatabase, QIcon, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox,
                               QFormLayout, QLineEdit, QListView, QListWidget, QListWidgetItem,
                               QPlainTextEdit, QSpinBox, QStyle, QTableWidget, QTableWidgetItem,
                               QToolButton, QVBoxLayout, QWidget)

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


def make_table(headers: list) -> QTableWidget:
    """Schreibgeschuetzte Tabelle ohne Zeilennummern."""
    t = QTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    t.verticalHeader().hide()
    t.horizontalHeader().setStretchLastSection(True)
    return t


def fill_table(table: QTableWidget, rows: list, fit: bool = False) -> None:
    """Zeilen setzen; ``fit``: Hoehe auf die Zeilen einpassen (hoechstens 8)."""
    table.setRowCount(len(rows))
    for i, row in enumerate(rows):
        for j, text in enumerate(row):
            table.setItem(i, j, QTableWidgetItem("" if text is None else str(text)))
    table.resizeColumnsToContents()
    if fit:
        n = min(max(len(rows), 1), 8)
        table.setFixedHeight(table.horizontalHeader().sizeHint().height() + n * table.verticalHeader().defaultSectionSize()
                             + table.horizontalScrollBar().sizeHint().height() + 2 * table.frameWidth())


def report_view() -> QPlainTextEdit:
    """Textbericht (statistics.txt, Indizes-Bericht) in Festbreitenschrift."""
    v = QPlainTextEdit()
    v.setReadOnly(True)
    v.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
    v.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
    return v


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

    def __init__(self, settings_cls, fields=None, advanced=None, parent=None):
        """``fields``: sichtbare Felder (Standard: Ebene basic); ``advanced``: aufklappbar unter "Erweitert"."""
        super().__init__(parent)
        self.cls = settings_cls
        metas = settings_cls.field_meta()
        hints = typing.get_type_hints(settings_cls)
        self._rows: dict = {}
        self._disabled: dict = {}
        self._errors: dict = {}
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        basic = QFormLayout()
        outer.addLayout(basic)
        sections = [(basic, fields or [n for n, m in metas.items() if m.level == "basic"])]
        self.advanced_button = self.advanced_box = None
        if advanced:
            self.advanced_button = QToolButton()
            self.advanced_button.setText("Erweitert")
            self.advanced_button.setCheckable(True)
            self.advanced_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
            self.advanced_button.setArrowType(Qt.ArrowType.RightArrow)
            self.advanced_button.toggled.connect(self._toggle_advanced)
            self.advanced_box = QWidget()
            self.advanced_box.hide()
            outer.addWidget(self.advanced_button)
            outer.addWidget(self.advanced_box)
            sections.append((QFormLayout(self.advanced_box), advanced))
        for layout, names in sections:
            for name in names:
                meta = metas[name]
                typ, optional = _unwrap(hints[name])
                kind, w = self._widget(meta, typ, optional)
                self._rows[name] = (kind, w, typ, optional, meta)
                self._update_row(name)
                label = meta.label + (f" [{meta.unit}]" if meta.unit else "")
                layout.addRow("" if kind == "bool" else label, w)
        self._base = settings_cls()
        self.set_settings(self._base)

    def _toggle_advanced(self, on: bool) -> None:
        self.advanced_box.setVisible(on)
        self.advanced_button.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)

    def _widget(self, meta, typ, optional) -> tuple:
        emit = lambda *_: self.changed.emit()   # noqa: E731
        if meta.choices:
            w = QComboBox()
            if optional:
                w.addItem("aus", None)
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
            w.setPlaceholderText("automatisch" if meta.kind == "names" or typ in (int, float) else "keine")
        w.textChanged.connect(emit)
        return ("names" if meta.kind == "names" else "text"), w

    def _update_row(self, name: str) -> None:
        """Sperre, Fehlerrahmen und Tooltip (Hilfe, Sperrgrund, Befund) eines Felds."""
        kind, w, _typ, _opt, meta = self._rows[name]
        reason, error = self._disabled.get(name), self._errors.get(name)
        w.setEnabled(reason is None)
        w.setStyleSheet(f"{type(w).__name__} {{ border: 1px solid #d32f2f; }}" if error else "")
        w.setToolTip("\n\n".join(t for t in (meta.help, reason and f"Gesperrt: {reason}", error) if t))

    def set_disabled(self, reasons: dict) -> None:
        """Felder sperren, die in dieser Lage nicht wirken (``Settings.disabled_fields``)."""
        self._disabled = {k: v for k, v in reasons.items() if k in self._rows}
        for name in self._rows:
            self._update_row(name)

    def mark_issues(self, issues) -> None:
        """Felder mit Fehler-Befund rot umranden (Befund im Tooltip)."""
        self._errors = {}
        for issue in issues:
            d = issue if isinstance(issue, dict) else issue.to_dict()
            if d.get("level") == "error" and d.get("field") in self._rows:
                self._errors.setdefault(d["field"], d["message_de"])
        for name in self._rows:
            self._update_row(name)

    def widget(self, name: str):
        """Eingabe-Widget eines Felds."""
        return self._rows[name][1]

    def set_placeholder(self, name: str, text: str) -> None:
        w = self._rows[name][1]
        if isinstance(w, QLineEdit):
            w.setPlaceholderText(text)

    def set_settings(self, s) -> None:
        for name, (kind, w, _typ, optional, _meta) in self._rows.items():
            v = getattr(s, name)
            w.blockSignals(True)
            if kind == "choice":
                w.setCurrentIndex(0 if v is None and optional else max(0, w.findData(v)))
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

    max_height = 160

    def set_issues(self, issues) -> None:
        self.clear()
        for issue in issues:
            d = issue if isinstance(issue, dict) else issue.to_dict()
            text = d["message_de"] + (f"\n→ {d['hint_de']}" if d.get("hint_de") else "")
            item = QListWidgetItem(self.style().standardIcon(_LEVEL_ICONS[d["level"]]), text)
            item.setToolTip(d.get("detail") or d.get("code", ""))
            self.addItem(item)
        self.setVisible(bool(issues))
        rows = sum(self.sizeHintForRow(i) for i in range(self.count()))
        self.setFixedHeight(min(rows + 2 * self.frameWidth() + 6, self.max_height))


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
