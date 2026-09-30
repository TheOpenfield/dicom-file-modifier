"""
Bausteine der Seiten: Einstellungsformular aus ``FieldMeta``, Statuszeile,
Befundliste, Bildansicht, Tabellen und Textbericht.
"""

from __future__ import annotations

import os
import types
import typing
from pathlib import Path

from PySide6.QtCore import QSize, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFontDatabase, QIcon, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox,
                               QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                               QListView, QListWidget, QListWidgetItem, QMessageBox,
                               QPlainTextEdit, QSizePolicy, QSpinBox, QStyle, QTableWidget,
                               QTableWidgetItem, QToolButton, QVBoxLayout, QWidget)

INT_MAX = 2**31 - 1
ZWSP = chr(0x200B)                      # Leerzeichen der Breite null: erlaubt Umbruch in Pfaden
STATUS_DE = {"ok": "Fertig", "ok_warnings": "Fertig, mit Hinweisen",
             "failed": "Fehlgeschlagen", "cancelled": "Abgebrochen"}
# Anzeige der Choice-Werte (Feldname -> Wert -> Text); fehlende Werte zeigen sich selbst
CHOICE_LABELS = {
    "eclipse_compat": {"high": "high – 1 CT-Pixel", "default": "default – 2 CT-Pixel"},
    "dose_interp": {"linear": "linear", "cubic": "cubic – B-Spline"},
    "volume_model": {"slab": "slab – volle Schichtdicke", "eclipse": "eclipse – Endschichten halb"},
    "piv_scope": {"component": "component – nur Isodosen-Anteile am Ziel",
                  "global": "global – ganze Isodose"},
    "iso_contours": {"mask": "mask – Maskenkanten", "field": "field – Isolinie wie Eclipse"},
    "transfer_syntax": {"explicit": "Explicit VR Little Endian", "implicit": "Implicit VR Little Endian"},
    "method": {"resample": "resample – neu abtasten, axiale Schichten",
               "metadata": "metadata – nur Lage ändern, Pixel unverändert"},
    "order": {0: "0 – nächster Nachbar (HU exakt)", 1: "1 – linear", 3: "3 – kubisch"},
}
SPEC_EXAMPLES = {"isodose": "z. B. 100,80,50,12Gy", "eclipse_values": "keine – z. B. TV=1.20,PIV=1.45",
                 "center": "volume, marker:NAME oder x,y,z"}
PATH_FILTERS = {"append_csv": "CSV-Dateien (*.csv)"}
IMAGE_TITLES = {"volumes": "Volumina", "shape_metrics": "Formmetriken",
                "distances": "Abstände Ziel – Risikoorgan", "centroids_3d": "Schwerpunkte 3D",
                "proximity_matrix": "Nähe-Matrix", "nearest_critical_oar": "Nächstes kritisches Organ",
                "sphericity_vs_elongation": "Sphärizität und Elongation",
                "gtv_ptv_margin": "GTV–PTV-Saum", "dose_overview": "Dosisübersicht",
                "transform_overview": "Vorher/Nachher", "displacement": "Verschiebung je ROI"}


def breakable(path) -> str:
    """Pfad fuer QLabel mit Zeilenumbruch: nach jedem Trenner darf umbrochen werden."""
    return str(path).replace(os.sep, os.sep + ZWSP).replace("/", "/" + ZWSP)


def open_path(path) -> None:
    """Datei oder Ordner mit dem Standardprogramm von Windows oeffnen."""
    QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))


def std_icon(widget: QWidget, name: str) -> QIcon:
    return widget.style().standardIcon(getattr(QStyle.StandardPixmap, name))


# ---------------------------------------------------------------------------
# Tabellen und Textbericht
# ---------------------------------------------------------------------------

class _Item(QTableWidgetItem):
    """Zelle, die nach ``sort_value`` sortiert (Zahlen numerisch)."""

    def __init__(self, text: str, sort_value=None):
        super().__init__(text)
        self.sort_value = sort_value

    def __lt__(self, other) -> bool:
        a, b = self.sort_value, getattr(other, "sort_value", None)
        if a is not None and b is not None:
            return a < b
        return self.text().lower() < other.text().lower()


def make_table(headers: list, sortable: bool = False) -> QTableWidget:
    """Schreibgeschuetzte Tabelle ohne Zeilennummern."""
    t = QTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    t.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    t.verticalHeader().hide()
    t.horizontalHeader().setStretchLastSection(True)
    t.setProperty("sortable", sortable)
    t.horizontalHeader().setSortIndicator(-1, Qt.SortOrder.AscendingOrder)   # erst auf Klick sortieren
    return t


def fill_table(table: QTableWidget, rows: list, fit: bool = False) -> None:
    """
    Zeilen setzen.  Eine Zelle ist Text oder ``(text, zahl)``: Zahlen stehen
    rechtsbuendig und sortieren numerisch.  ``fit``: Hoehe auf die Zeilen
    einpassen (hoechstens 8).
    """
    sortable = bool(table.property("sortable"))
    table.setSortingEnabled(False)
    table.setRowCount(len(rows))
    for i, row in enumerate(rows):
        for j, cell in enumerate(row):
            text, value = cell if isinstance(cell, tuple) else ("" if cell is None else str(cell), None)
            item = _Item(text, value)
            if value is not None:
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            table.setItem(i, j, item)
    table.resizeColumnsToContents()
    table.setSortingEnabled(sortable)
    if fit:
        n = min(max(len(rows), 1), 8)
        table.setFixedHeight(table.horizontalHeader().sizeHint().height()
                             + n * table.verticalHeader().defaultSectionSize()
                             + table.horizontalScrollBar().sizeHint().height() + 2 * table.frameWidth())


def report_view() -> QPlainTextEdit:
    """Textbericht (statistics.txt, Indizes-Bericht) in Festbreitenschrift."""
    v = QPlainTextEdit()
    v.setReadOnly(True)
    v.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
    v.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
    return v


# ---------------------------------------------------------------------------
# Einstellungsformular
# ---------------------------------------------------------------------------

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
    Gesperrte Felder zeigen den Grund, Felder mit Fehler ein rotes Label (und
    ein Symbol im Textfeld); keine Stylesheets auf Eingabefeldern, damit der
    native Windows-Stil erhalten bleibt.
    """

    changed = Signal()

    def __init__(self, settings_cls, fields=None, advanced=None, parent=None):
        """``fields``: sichtbare Felder (Standard: Ebene basic); ``advanced``: aufklappbar unter "Erweitert"."""
        super().__init__(parent)
        self.cls = settings_cls
        metas = settings_cls.field_meta()
        hints = typing.get_type_hints(settings_cls)
        self._rows: dict = {}
        self._labels: dict = {}
        self._error_actions: dict = {}
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
                kind, w = self._widget(name, meta, typ, optional)
                self._rows[name] = (kind, w, typ, optional, meta)
                label = QLabel("" if kind == "bool" else meta.label + (f" [{meta.unit}]" if meta.unit else ""))
                self._labels[name] = label
                layout.addRow(label, self._with_browse(name, w) if meta.kind == "path" else w)
                self._update_row(name)
        self._base = settings_cls()
        self.set_settings(self._base)

    def _toggle_advanced(self, on: bool) -> None:
        self.advanced_box.setVisible(on)
        self.advanced_button.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)

    def _widget(self, name, meta, typ, optional) -> tuple:
        emit = lambda *_: self.changed.emit()   # noqa: E731
        if meta.choices:
            w = QComboBox()
            if optional:
                w.addItem("aus", None)
            labels = CHOICE_LABELS.get(name, {})
            for c in meta.choices:
                w.addItem(labels.get(c, str(c)), c)
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
                w.setDecimals(2)
                w.setRange(float(lo), float(hi))
            else:
                w.setRange(int(max(lo, -INT_MAX)), int(min(hi, INT_MAX)))
            w.valueChanged.connect(emit)
            return "number", w
        w = QLineEdit()
        if name in SPEC_EXAMPLES:
            w.setPlaceholderText(SPEC_EXAMPLES[name])
        elif optional:
            w.setPlaceholderText("automatisch" if meta.kind == "names" or typ in (int, float) else "keine")
        w.textChanged.connect(emit)
        return ("names" if meta.kind == "names" else "text"), w

    def _with_browse(self, name: str, edit: QLineEdit) -> QWidget:
        """Textfeld mit "…"-Knopf fuer eine Datei (``kind="path"``)."""
        box = QWidget()
        row = QHBoxLayout(box)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(edit, 1)
        button = QToolButton()
        button.setText("…")
        button.setToolTip("Datei wählen")

        def browse() -> None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Datei wählen", edit.text(), PATH_FILTERS.get(name, "Alle Dateien (*)"),
                options=QFileDialog.Option.DontConfirmOverwrite)
            if path:
                edit.setText(str(Path(path)))

        button.clicked.connect(browse)
        row.addWidget(button)
        return box

    def _update_row(self, name: str) -> None:
        """Sperre, Fehlermarkierung und Tooltip (Hilfe, Sperrgrund, Befund) eines Felds."""
        kind, w, _typ, _opt, meta = self._rows[name]
        reason, error = self._disabled.get(name), self._errors.get(name)
        w.setEnabled(reason is None)
        label = self._labels.get(name)
        if label is not None:
            label.setStyleSheet("color: #c62828;" if error else "")
        if isinstance(w, QLineEdit):
            action = self._error_actions.get(name)
            if error and action is None:
                self._error_actions[name] = w.addAction(std_icon(self, "SP_MessageBoxCritical"),
                                                        QLineEdit.ActionPosition.TrailingPosition)
            elif not error and action is not None:
                w.removeAction(self._error_actions.pop(name))
        w.setToolTip("\n\n".join(t for t in (meta.help, reason and f"Gesperrt: {reason}", error) if t))

    def set_disabled(self, reasons: dict) -> None:
        """Felder sperren, die in dieser Lage nicht wirken (``Settings.disabled_fields``)."""
        self._disabled = {k: v for k, v in reasons.items() if k in self._rows}
        for name in self._rows:
            self._update_row(name)

    def mark_issues(self, issues) -> None:
        """Felder mit Fehler-Befund markieren (Befund im Tooltip)."""
        self._errors = {}
        for issue in issues:
            d = issue if isinstance(issue, dict) else issue.to_dict()
            if d.get("level") == "error" and d.get("field") in self._rows:
                self._errors.setdefault(d["field"], d["message_de"])
        for name in self._rows:
            self._update_row(name)

    def error_fields(self) -> set:
        return set(self._errors)

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

    def label_of(self, name: str) -> str:
        return self.cls.field_meta()[name].label


# ---------------------------------------------------------------------------
# Statuszeile, Befunde, Bilder
# ---------------------------------------------------------------------------

_LEVEL_ICONS = {"error": "SP_MessageBoxCritical", "warning": "SP_MessageBoxWarning",
                "info": "SP_MessageBoxInformation", "ok": "SP_DialogApplyButton"}


class StateLine(QWidget):
    """Symbol und fetter Text: Zustand der Pruefung ("Bereit", "Nicht startbar", ...)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        self.icon = QLabel()
        self.text = QLabel()
        font = self.text.font()
        font.setBold(True)
        self.text.setFont(font)
        row.addWidget(self.icon)
        row.addWidget(self.text, 1)

    def set_state(self, level, text: str) -> None:
        self.text.setText(text)
        self.icon.setPixmap(std_icon(self, _LEVEL_ICONS[level]).pixmap(16, 16) if level else QPixmap())


class IssueList(QListWidget):
    """Befunde (``Issue`` oder dict) mit Symbol je Stufe; Hinweis in der zweiten Zeile, Doppelklick: Details."""

    max_height = 160

    def __init__(self, parent=None):
        super().__init__(parent)
        self.itemDoubleClicked.connect(self._details)

    def set_issues(self, issues) -> None:
        self.clear()
        for issue in issues:
            d = issue if isinstance(issue, dict) else issue.to_dict()
            text = d["message_de"] + (f"\n→ {d['hint_de']}" if d.get("hint_de") else "")
            item = QListWidgetItem(std_icon(self, _LEVEL_ICONS[d["level"]]), text)
            item.setData(Qt.ItemDataRole.UserRole, d)
            item.setToolTip("Doppelklick: Details")
            self.addItem(item)
        self.setVisible(bool(issues))
        rows = sum(self.sizeHintForRow(i) for i in range(self.count()))
        self.setFixedHeight(min(rows + 2 * self.frameWidth() + 6, self.max_height))

    def _details(self, item: QListWidgetItem) -> None:
        d = item.data(Qt.ItemDataRole.UserRole)
        box = QMessageBox(self)
        box.setIcon({"error": QMessageBox.Icon.Critical, "warning": QMessageBox.Icon.Warning}
                    .get(d["level"], QMessageBox.Icon.Information))
        box.setWindowTitle("Befund")
        box.setText(d["message_de"])
        box.setInformativeText("\n".join(t for t in (d.get("hint_de"), f"Code: {d.get('code', '')}") if t))
        if d.get("detail"):
            box.setDetailedText(d["detail"])
        box.exec()


class _ClickLabel(QLabel):
    double_clicked = Signal()

    def mouseDoubleClickEvent(self, event) -> None:
        self.double_clicked.emit()


class ImageViewer(QWidget):
    """
    Grosse Vorschau des gewaehlten Bilds, darunter die Vorschaubilder (bei mehr
    als einem Bild).  Doppelklick oeffnet das Bild in der Bildanzeige von Windows.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.view = _ClickLabel()
        self.view.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.view.setMinimumSize(200, 150)
        self.view.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.view.setToolTip("Doppelklick: in der Bildanzeige öffnen")
        self.view.double_clicked.connect(lambda: self._path and open_path(self._path))
        self.strip = QListWidget()
        self.strip.setViewMode(QListView.ViewMode.IconMode)
        self.strip.setFlow(QListView.Flow.LeftToRight)
        self.strip.setWrapping(False)
        self.strip.setMovement(QListView.Movement.Static)
        self.strip.setIconSize(QSize(120, 85))
        self.strip.setSpacing(4)
        self.strip.setFixedHeight(140)
        self.strip.currentItemChanged.connect(self._select)
        self.strip.itemActivated.connect(lambda item: open_path(item.data(Qt.ItemDataRole.UserRole)))
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.view, 1)
        lay.addWidget(self.strip)
        self._pix = self._path = None

    def count(self) -> int:
        return self.strip.count()

    def clear(self) -> None:
        self.strip.clear()
        self.view.clear()
        self._pix = self._path = None

    def set_images(self, paths) -> None:
        self.clear()
        for p in paths:
            pix = QPixmap(str(p))
            if pix.isNull():
                continue
            thumb = pix.scaled(self.strip.iconSize(), Qt.AspectRatioMode.KeepAspectRatio,
                               Qt.TransformationMode.SmoothTransformation)
            item = QListWidgetItem(QIcon(thumb), IMAGE_TITLES.get(Path(p).stem, Path(p).stem))
            item.setData(Qt.ItemDataRole.UserRole, str(p))
            item.setToolTip(f"{Path(p).name} – Doppelklick: in der Bildanzeige öffnen")
            self.strip.addItem(item)
        self.strip.setVisible(self.strip.count() > 1)
        if self.strip.count():
            self.strip.setCurrentRow(0)

    def _select(self, item, _previous=None) -> None:
        if item is None:
            return
        self._path = item.data(Qt.ItemDataRole.UserRole)
        self._pix = QPixmap(self._path)
        self._rescale()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._rescale()

    def _rescale(self) -> None:
        if self._pix is not None and not self._pix.isNull():
            self.view.setPixmap(self._pix.scaled(self.view.size(), Qt.AspectRatioMode.KeepAspectRatio,
                                                 Qt.TransformationMode.SmoothTransformation))
