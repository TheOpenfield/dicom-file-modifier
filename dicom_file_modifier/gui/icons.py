"""
App-Icon, Symbole der Oberflaeche und die Farben der Marke.  Die Dateien liegen
in ``assets/``: ``app.png``/``app.ico`` (Fenster, Taskleiste, EXE) und
SVG-Symbole fuer Seitenleiste und Werkzeugleiste (scharf in jeder Skalierung).
"""

from __future__ import annotations

from pathlib import Path

from PySide6 import QtSvg  # noqa: F401  # Qt6Svg ins Bundle: PyInstaller sammelt nur importierte Module
from PySide6.QtCore import QRect, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QIcon, QPainter, QPixmap

ASSETS = Path(__file__).parent / "assets"
# Farben des App-Icons: Kachel, Rotationsbogen, Zielvolumen, Isodosen, CT-Schicht
NAVY, TEAL, CORAL, ORANGE, YELLOW, SLICE = "#19355E", "#30C5BB", "#E35C35", "#F0962F", "#FAD641", "#D6D9E1"
PAGE_ICONS = {"structures": "page_structures", "dose": "page_dose", "transform": "page_transform"}


def app_icon() -> QIcon:
    return QIcon(str(ASSETS / "app.png"))


def icon(name: str) -> QIcon:
    """Symbol aus ``assets/<name>.svg`` (open, results, log, info, page_*)."""
    return QIcon(str(ASSETS / f"{name}.svg"))


def page_icon(workflow: str) -> QIcon:
    return icon(PAGE_ICONS[workflow])


def placeholder_pixmap(symbol: QIcon, text: str, font: QFont, color: QColor,
                       dpr: float = 1.0, size: int = 96) -> QPixmap:
    """Blasses Symbol ueber einem Hinweistext, z.B. "Noch kein Ergebnis" in der Bildansicht."""
    fm = QFontMetrics(font)
    w = max(size, fm.horizontalAdvance(text)) + 16
    h = size + 10 + fm.height()
    pix = QPixmap(round(w * dpr), round(h * dpr))
    pix.setDevicePixelRatio(dpr)
    pix.fill(Qt.GlobalColor.transparent)
    p = QPainter(pix)
    p.setOpacity(0.3)
    symbol.paint(p, QRect((w - size) // 2, 0, size, size))
    p.setOpacity(1.0)
    p.setFont(font)
    p.setPen(color)
    p.drawText(QRect(0, size + 10, w, fm.height()), Qt.AlignmentFlag.AlignCenter, text)
    p.end()
    return pix
