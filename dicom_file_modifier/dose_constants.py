"""
dose_constants.py - Gemeinsame Konstanten der Dosisindex-Module.

Eigenes Modul, damit ``rtstruct_writer`` und ``dose_viz`` die Konstanten ohne
Import von ``dose_indices`` nutzen koennen (sonst entstuende ein Import-Zyklus).
Das Modul importiert nur die Paketversion.
"""

from . import __version__

TOOL_NAME = "dose_indices"
TOOL_VERSION = __version__          # eine Nummer fuer Paket, Werkzeuge und Bundle

# Farbvorschlaege (RGB) fuer Isodosen-ROIs nach Prozent-Level
LEVEL_COLORS = {
    100: (255, 0, 255), 95: (255, 128, 255), 90: (255, 105, 180), 80: (255, 165, 0),
    70: (255, 215, 0), 60: (0, 200, 100), 50: (0, 255, 255), 30: (0, 128, 255), 20: (0, 0, 255),
}
HELPER_COLORS = {"inter": (0, 255, 0), "under": (255, 255, 0), "spill": (255, 0, 0)}
