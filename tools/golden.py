#!/usr/bin/env python
"""
golden.py - Golden-Output-Harness fuer die CLI-Paritaet (Plan P0.1).

Ruft die CLIs des Pakets als Subprozess mit dem aktuellen Interpreter auf,
normalisiert die Ausgaben (UIDs, Zeitstempel, Pfade) zu Schnappschuessen und
vergleicht zwei Laeufe exakt oder mit Migrationstoleranzen.

  python tools/golden.py make-inputs [--work DIR]
  python tools/golden.py run --label G0 [--work DIR] [--only MUSTER] [--real CASE ...]
  python tools/golden.py compare G0 G3 [--work DIR] [--mode exact|migration] [--expected] [--stats]

Arbeitsordner (Default ``%USERPROFILE%\\dfm-golden``) liegt bewusst ausserhalb
des Repos: ``inputs/`` (synthetische Demo-Faelle + Manifest), ``runs/<label>/``
(Schnappschuesse, stdout/stderr, Pixel-Arrays).  Echte Faelle (``--real``)
werden nur lokal ausgefuehrt; ihr Ordnername wird in Schnappschuessen und
Berichten durch ``REAL<n>`` ersetzt, Textdiffs zeigen dort keinen Inhalt.
Je CLI-Aufruf werden Laufzeit und (unter Windows) Peak-Speicher mitgeschrieben;
``compare --stats`` haengt Driftstatistik, Laufzeit/Speicher und
stderr-Warnungen an den Bericht an.

Laeuft auf Python >= 3.8 (auch in der alten Umgebung, siehe Plan P0.2).
Konsolenausgabe ist ASCII.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import io
import json
import math
import os
import platform
import re
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

REPO = Path(__file__).resolve().parents[1]
RULES_PATH = Path(__file__).resolve().with_name("golden_rules.json")
DEFAULT_WORK = Path.home() / "dfm-golden"

# ---------------------------------------------------------------------------
# Eingaben und Szenarien
# ---------------------------------------------------------------------------

VARIANTS = {
    "std": [],
    "flat": ["--layout", "flat"],
    "beams": ["--rd-set", "plan+2beams"],
    "2plans": ["--rd-set", "2plans"],
    "nodvh": ["--no-dvh"],
    "evr": ["--explicit-vr"],
    "rle": ["--rle"],
    "ecl1": ["--eclipse-ref", "ptv1"],
    "ecl2": ["--eclipse-ref", "both"],
}

T6 = ["--tx", "5", "--ty", "-3", "--tz", "2", "--rx", "2", "--ry", "-3", "--rz", "10"]


def _sc(name, module, args=None, case=None, steps=None, tags=(), timeout=900, script=None):
    return {"name": name, "module": module, "args": args or [], "case": case,
            "steps": steps, "tags": list(tags), "timeout": timeout, "script": script}


SCENARIOS = [
    # argparse-Oberflaeche aller CLIs (Flags, Defaults, Choices, Typen)
    _sc("cli_surface", None, ["{out}"], script="cli_surface.py"),
    # Self-Tests (ohne Eingaben)
    _sc("ana_self_test", "analyzer", ["--self-test"]),
    _sc("dose_self_test", "dose_indices", ["--self-test"]),
    _sc("rsw_self_test", "rtstruct_writer", ["--self-test"]),
    # analyzer / visualizer
    _sc("ana_std", "analyzer", ["{rs}", "--output", "{out}"], case="std"),
    _sc("ana_list", "analyzer", ["{rs}", "--list"], case="std"),
    _sc("ana_filter", "analyzer", ["{rs}", "--targets", "PTV,GTV", "--oars", "Hirnstamm,Chiasma",
                                   "--output", "{out}"], case="std"),
    _sc("viz_std", "visualizer", ["{rs}", "--output", "{out}"], case="std"),
    # modifier (CT-only)
    _sc("mod_res_o0", "modifier", ["{ct}"] + T6 + ["--order", "0", "--no-viz", "--output", "{out}"], case="std"),
    _sc("mod_res_o1", "modifier", ["{ct}"] + T6 + ["--order", "1", "--no-viz", "--output", "{out}"], case="std"),
    _sc("mod_res_o3", "modifier", ["{ct}"] + T6 + ["--order", "3", "--no-viz", "--output", "{out}"], case="std"),
    _sc("mod_meta", "modifier", ["{ct}"] + T6 + ["--method", "metadata", "--no-viz", "--output", "{out}"],
        case="std"),
    _sc("mod_viz", "modifier", ["{ct}", "--tx", "5", "--rz", "10", "--output", "{out}"], case="std"),
    _sc("mod_flat", "modifier", ["{case}", "--tx", "3", "--no-viz", "--output", "{out}"], case="flat",
        tags=["known_bug"]),
    _sc("mod_rle", "modifier", ["{ct}", "--tx", "3", "--no-viz", "--output", "{out}"], case="rle"),
    # case_modifier
    _sc("cm_list_markers", "case_modifier", ["{case}", "--list-markers"], case="std"),
    _sc("cm_dry", "case_modifier", ["{case}"] + T6 + ["--center", "marker:HS1", "--dry-run", "--output", "{out}"],
        case="std"),
    _sc("cm_meta_verify", "case_modifier", ["{case}"] + T6 + ["--method", "metadata", "--center", "marker:HS1",
                                                               "--verify", "--no-viz", "--output", "{out}"],
        case="std"),
    _sc("cm_res_verify", "case_modifier", ["{case}"] + T6 + ["--order", "1", "--center", "volume",
                                                              "--non-interactive", "--verify", "--output", "{out}"],
        case="std"),
    _sc("cm_res_newfor", "case_modifier", ["{case}", "--tx", "-4", "--rz", "-7", "--center", "0,0,0",
                                           "--new-frame-of-reference", "--label", "_NF", "--no-viz",
                                           "--output", "{out}"], case="std"),
    _sc("cm_self_test", "case_modifier", ["{case}", "--self-test"], case="std"),
    _sc("cm_flat", "case_modifier", ["{case}", "--tx", "3", "--center", "volume", "--no-viz",
                                     "--output", "{out}"], case="flat"),
    _sc("cm_rle", "case_modifier", ["{case}", "--tx", "3", "--center", "volume", "--no-viz",
                                    "--output", "{out}"], case="rle"),
    # dose_indices
    _sc("dose_default", "dose_indices", ["{case}", "--output", "{out}"], case="std"),
    _sc("dose_list", "dose_indices", ["{case}", "--list"], case="std"),
    _sc("dose_ecl_high", "dose_indices", ["{case}", "--eclipse-compat", "high", "--label", "_ECL",
                                          "--output", "{out}"], case="std"),
    _sc("dose_ecl_default", "dose_indices", ["{case}", "--eclipse-compat", "default", "--label", "_ECD",
                                             "--no-viz", "--output", "{out}"], case="std"),
    _sc("dose_no_rs", "dose_indices", ["{case}", "--no-rs", "--no-viz", "--output", "{out}"], case="std"),
    _sc("dose_levels", "dose_indices", ["{case}", "--isodose", "100,80,50,12Gy", "--no-viz",
                                        "--output", "{out}"], case="std"),
    _sc("dose_variant", "dose_indices", ["{case}", "--grid", "0.5", "--dose-interp", "cubic",
                                         "--volume-model", "eclipse", "--piv-scope", "global", "--no-viz",
                                         "--output", "{out}"], case="std"),
    _sc("dose_two_targets_fail", "dose_indices", ["{case}", "--target", "PTV_1,PTV_2", "--no-viz",
                                                  "--output", "{out}"], case="std"),
    _sc("dose_two_targets", "dose_indices", ["{case}", "--target", "PTV_1,PTV_2", "--rx", "20", "--no-viz",
                                             "--output", "{out}"], case="std"),
    _sc("dose_append_csv", "dose_indices", case="std", steps=[
        ["{case}", "--no-viz", "--no-rs", "--append-csv", "{tmp}/sammel.csv", "--output", "{out}"],
        ["{case}", "--no-viz", "--no-rs", "--label", "_B", "--append-csv", "{tmp}/sammel.csv",
         "--output", "{out}"]]),
    _sc("dose_ecl_values", "dose_indices", ["{case}", "--eclipse-values", "PIV=5.65,V10Gy=14.4", "--no-viz",
                                            "--output", "{out}"], case="std"),
    _sc("dose_rx_gy_crash", "dose_indices", ["{case}", "--isodose", "20Gy", "--no-viz", "--output", "{out}"],
        case="std"),
    _sc("dose_rx_pct", "dose_indices", ["{case}", "--rx-pct-of-max", "80", "--no-viz", "--output", "{out}"],
        case="std"),
    _sc("dose_flat", "dose_indices", ["{case}", "--no-viz", "--output", "{out}"], case="flat"),
    _sc("dose_beams", "dose_indices", ["{case}", "--no-viz", "--output", "{out}"], case="beams"),
    _sc("dose_2plans", "dose_indices", ["{case}", "--no-viz", "--output", "{out}"], case="2plans"),
    _sc("dose_nodvh", "dose_indices", ["{case}", "--no-viz", "--output", "{out}"], case="nodvh"),
    _sc("dose_evr", "dose_indices", ["{case}", "--no-viz", "--output", "{out}"], case="evr"),
    _sc("dose_ecl1", "dose_indices", ["{case}", "--no-viz", "--output", "{out}"], case="ecl1"),
    _sc("dose_ecl2", "dose_indices", ["{case}", "--no-viz", "--output", "{out}"], case="ecl2"),
]

# Szenarien je echtem Fall (nur lokal; Voraussetzungen werden geprueft)
REAL_SCENARIOS = [
    _sc("ana", "analyzer", ["{rs}", "--output", "{out}"], tags=["needs_rs"]),
    _sc("viz", "visualizer", ["{rs}", "--output", "{out}"], tags=["needs_rs"]),
    _sc("dose_default", "dose_indices", ["{case}", "--output", "{out}"], tags=["needs_rd"]),
    _sc("dose_ecl_high", "dose_indices", ["{case}", "--eclipse-compat", "high", "--label", "_ECL",
                                          "--output", "{out}"], tags=["needs_rd", "needs_ct"]),
    _sc("cm_meta", "case_modifier", ["{case}"] + T6 + ["--method", "metadata", "--center", "volume",
                                                        "--non-interactive", "--verify", "--no-viz",
                                                        "--output", "{out}"], tags=["needs_rs", "needs_ct"]),
    _sc("cm_res", "case_modifier", ["{case}"] + T6 + ["--order", "1", "--center", "volume",
                                                       "--non-interactive", "--verify", "--no-viz",
                                                       "--output", "{out}"], tags=["needs_rs", "needs_ct"],
        timeout=1800),
]

# ---------------------------------------------------------------------------
# Hilfen
# ---------------------------------------------------------------------------


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _versions() -> dict:
    out = {"python": sys.version.split()[0], "platform": platform.platform(),
           "executable_prefix": Path(sys.prefix).name}
    try:
        from importlib.metadata import version as _v   # py >= 3.8
    except ImportError:                                  # pragma: no cover
        return out
    for pkg in ("numpy", "scipy", "pydicom", "shapely", "matplotlib", "scikit-image", "plotly"):
        try:
            out[pkg] = _v(pkg)
        except Exception:
            out[pkg] = None
    for var in ("NPY_PROMOTION_STATE", "PYDICOM_FUTURE"):
        if os.environ.get(var):
            out["env_" + var] = os.environ[var]
    return out


def _child_env() -> dict:
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = "0"
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONPATH"] = str(REPO) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env.setdefault("MPLBACKEND", "Agg")
    return env


def _ascii(s: str) -> str:
    return s.encode("ascii", "backslashreplace").decode("ascii")


# ---------------------------------------------------------------------------
# Normalisierung
# ---------------------------------------------------------------------------

# UID auch in Dateinamen wie "RS.<uid>.dcm" (Punkt davor nur, wenn keine Ziffer davor steht)
_UID_RE = re.compile(r"(?<![0-9])(?<![0-9]\.)[0-2](?:\.[0-9]+){5,}(?![0-9]|\.[0-9])")
_ISO_TS_RE = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}(?::\d{2})?(?:\.\d+)?")
_STANDARD_UID_PREFIX = "1.2.840.10008."


def _base_replacements() -> List[tuple]:
    """Pfade, die in jedem Lauf anders sein koennen (Repo, Interpreter, Temp, Home)."""
    import tempfile
    out = [(str(REPO), "<REPO>"), (sys.prefix, "<PREFIX>"), (sys.base_prefix, "<BASE_PREFIX>"),
           (tempfile.gettempdir(), "<SYSTMP>"), (str(Path.home()), "<HOME>")]
    return [(p, ph) for p, ph in out if p]


# Zufallssuffix von tempfile.mkdtemp/mkstemp (8 Zeichen) unter <SYSTMP>
_TMPNAME_RE = re.compile(r"(<SYSTMP>[/\\]+)([^/\\\s'\"]*?)[a-z0-9_]{8}(?=[/\\\s'\"]|$)")


class Normalizer:
    """UIDs -> ``UID<n>`` (erste Fundstelle), Pfade -> Platzhalter, Zeiten maskiert."""

    def __init__(self, replacements: List[tuple]):
        # laengste Pfade zuerst; Slash-Varianten und repr()-Form mit doppelten
        # Backslashes; case-insensitiv (Windows)
        pairs = []
        for raw, ph in list(replacements) + _base_replacements():
            back = raw.replace("/", "\\")
            for variant in {raw, raw.replace("\\", "/"), back, back.replace("\\", "\\\\")}:
                if variant:
                    pairs.append((variant, ph))
        pairs.sort(key=lambda t: -len(t[0]))
        self._path_res = [(re.compile(re.escape(v), re.IGNORECASE), ph) for v, ph in pairs]
        self.uids: Dict[str, str] = {}

    def uid(self, u: str) -> str:
        u = str(u).strip()
        if not u or u.startswith(_STANDARD_UID_PREFIX):
            return u
        if u not in self.uids:
            self.uids[u] = f"UID{len(self.uids) + 1}"
        return self.uids[u]

    def text(self, s: str) -> str:
        for rx, ph in self._path_res:
            s = rx.sub(ph.replace("\\", "\\\\"), s)
        s = _TMPNAME_RE.sub(lambda m: m.group(1) + m.group(2) + "*", s)
        s = _ISO_TS_RE.sub("<TS>", s)
        s = _UID_RE.sub(lambda m: self.uid(m.group(0)), s)
        return s.replace("\\", "/")


def _norm_json(obj, norm: Normalizer, key: str = ""):
    if isinstance(obj, dict):
        return {norm.text(str(k)): _norm_json(v, norm, str(k)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_norm_json(v, norm, key) for v in obj]
    if isinstance(obj, str):
        return "<TS>" if key in ("timestamp", "run_timestamp") else norm.text(obj)
    return obj


# ---------------------------------------------------------------------------
# Schnappschuss je Datei
# ---------------------------------------------------------------------------

def _png_dims(p: Path):
    try:
        with p.open("rb") as fh:
            head = fh.read(24)
        if head[:8] == b"\x89PNG\r\n\x1a\n":
            return list(struct.unpack(">II", head[16:24]))
    except OSError:
        pass
    return None


def _dicom_value(elem, norm: Normalizer):
    # pydicom 2.4 legt aufgeloeste mehrdeutige VRs als Enum ab; str() ergibt dort
    # unter Python < 3.11 'VR.US' statt 'US'
    vr = str(getattr(elem.VR, "value", elem.VR))
    val = elem.value
    if val is None or val == "":
        return None
    if vr == "UI":
        if isinstance(val, (list, tuple)) or type(val).__name__ == "MultiValue":
            return [norm.uid(v) for v in val]
        return norm.uid(val)
    if vr in ("DA", "TM", "DT"):
        return "<" + vr + ">"
    if vr in ("OB", "OW", "OF", "OD", "OL", "OV", "UN") or isinstance(val, (bytes, bytearray)):
        return {"bytes_sha256": _sha256_bytes(bytes(val)), "len": len(val)}
    if vr in ("DS", "IS"):
        seq = val if (isinstance(val, (list, tuple)) or type(val).__name__ == "MultiValue") else [val]
        out = [str(getattr(v, "original_string", None) or v).strip() for v in seq]
        return out if len(out) != 1 else out[0]
    seq = val if (isinstance(val, (list, tuple)) or type(val).__name__ == "MultiValue") else [val]
    # Mehrdeutiger VR ("US or SS"): pydicom 3 loest ihn beim Lesen auf, pydicom 2 nicht
    ambiguous_int = " or " in vr and all(isinstance(v, int) and not isinstance(v, bool) for v in seq)
    if vr in ("FL", "FD", "SL", "SS", "UL", "US", "SV", "UV") or ambiguous_int:
        out = [float(v) if vr in ("FL", "FD") else int(v) for v in seq]
        return out if len(out) != 1 else out[0]
    if vr == "AT":
        out = [_tag_key(str(v)) for v in seq]
        return out if len(out) != 1 else out[0]
    if isinstance(val, (list, tuple)) or type(val).__name__ == "MultiValue":
        return [norm.text(str(v)) for v in val]
    return norm.text(str(val))


# Abgeleitete Laenge: haengt von den UID-Laengen ab, und pydicom 3 erzeugt
# zufaellige UIDs mit 62-64 statt immer 64 Zeichen.  Die Elemente selbst
# werden einzeln verglichen.
_META_IGNORE = ("FileMetaInformationGroupLength",)


_TAG_KEY_RE = re.compile(r"^\(([0-9A-Fa-f]{4}),\s*([0-9A-Fa-f]{4})\)$")


def _tag_key(key: str) -> str:
    """Tag ohne Keyword einheitlich '(gggg,eeee)' (pydicom 2: '(0009, 0010)', 3: '(0009,0010)')."""
    m = _TAG_KEY_RE.match(key)
    return f"({m.group(1).lower()},{m.group(2).lower()})" if m else key


def _canon_tag_keys(obj):
    """Tag-Schluessel und AT-Werte alter Schnappschuesse einheitlich schreiben."""
    if isinstance(obj, dict):
        return {_tag_key(k): _canon_tag_keys(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_canon_tag_keys(v) for v in obj]
    return _tag_key(obj) if isinstance(obj, str) else obj


def _dicom_elements(ds, norm: Normalizer) -> dict:
    out = {}
    for elem in ds:
        key = elem.keyword or _tag_key(str(elem.tag))
        if key == "PixelData":
            continue
        if key == "ImplementationVersionName":
            out[key] = "<IVN>"
            continue
        if key in _META_IGNORE:
            out[key] = "<GL>"
            continue
        if str(elem.VR) == "SQ":
            out[key] = [_dicom_elements(item, norm) for item in elem.value]
        else:
            out[key] = _dicom_value(elem, norm)
    return out


def _dicom_snapshot(p: Path, norm: Normalizer, arrays: Optional[dict], arr_key: str) -> dict:
    import numpy as np
    import pydicom

    ds = pydicom.dcmread(str(p), force=True)
    snap = {"kind": "dicom", "size": p.stat().st_size,
            "file_meta": _dicom_elements(ds.file_meta, norm) if getattr(ds, "file_meta", None) else {},
            "dataset": _dicom_elements(ds, norm)}
    if "PixelData" in ds:
        try:
            arr = ds.pixel_array
            snap["pixels"] = {"sha256": _sha256_bytes(np.ascontiguousarray(arr).tobytes()),
                              "dtype": str(arr.dtype), "shape": list(arr.shape),
                              "min": float(arr.min()), "max": float(arr.max()),
                              "mean": float(arr.mean())}
            if arrays is not None and str(ds.get("Modality", "")).upper() == "CT":
                arrays[arr_key] = arr
        except Exception as e:  # noqa: BLE001 - Dekodierfehler sind Teil des Befunds
            raw = ds.PixelData if isinstance(ds.PixelData, (bytes, bytearray)) else b""
            snap["pixels"] = {"decode_error": f"{type(e).__name__}: {e}"[:300],
                              "raw_sha256": _sha256_bytes(bytes(raw))}
    return snap


def _file_snapshot(p: Path, norm: Normalizer, arrays: Optional[dict], arr_key: str) -> dict:
    suf = p.suffix.lower()
    size = p.stat().st_size
    try:
        if suf == ".dcm":
            return _dicom_snapshot(p, norm, arrays, arr_key)
        if suf == ".json":
            return {"kind": "json", "size": size,
                    "data": _norm_json(json.loads(p.read_text(encoding="utf-8")), norm)}
        if suf == ".csv":
            rows = list(csv.reader(io.StringIO(p.read_text(encoding="utf-8-sig"))))
            header = rows[0] if rows else []
            norm_rows = []
            for r in rows:
                norm_rows.append(["<TS>" if (i < len(header) and header[i] == "run_timestamp" and r is not rows[0])
                                  else norm.text(c) for i, c in enumerate(r)])
            return {"kind": "csv", "size": size, "rows": norm_rows}
        if suf in (".txt", ".log", ".md"):
            return {"kind": "text", "size": size,
                    "lines": norm.text(p.read_text(encoding="utf-8", errors="replace")).splitlines()}
        if suf == ".png":
            return {"kind": "png", "size": size, "dims": _png_dims(p)}
        if suf in (".html", ".htm"):
            txt = p.read_text(encoding="utf-8", errors="replace")
            # Inline-Skripte (z.B. eingebettetes plotly.js) enthalten URL-Strings im Code;
            # geprueft werden nur Tag-Attribute, die tatsaechlich nachladen wuerden.
            tags_only = re.sub(r"(?is)(<script\b[^>]*>).*?(</script>)", r"\1\2", txt)
            ext = re.search(r"""(?i)<(?:script|link|img|iframe)\b[^>]*\b(?:src|href)\s*=\s*["']https?://""",
                            tags_only)
            return {"kind": "html", "size": size, "external_src": bool(ext)}
    except Exception as e:  # noqa: BLE001
        return {"kind": "error", "size": size, "error": f"{type(e).__name__}: {e}"[:300]}
    return {"kind": "other", "size": size, "sha256": _sha256_file(p)}


# ---------------------------------------------------------------------------
# make-inputs
# ---------------------------------------------------------------------------

def cmd_make_inputs(args) -> int:
    inputs = Path(args.work) / "inputs"
    if inputs.exists():
        if not args.force:
            print(f"Eingaben existieren bereits: {inputs}  (--force zum Neuerzeugen)")
            return 2
        shutil.rmtree(inputs)
    inputs.mkdir(parents=True)
    env = _child_env()
    manifest = {"created": time.strftime("%Y-%m-%dT%H:%M:%S"), "versions": _versions(), "variants": {}}
    for name, vargs in VARIANTS.items():
        t0 = time.perf_counter()
        p = subprocess.run([sys.executable, "-m", "dicom_file_modifier.demo", str(inputs / name)] + vargs,
                           cwd=str(REPO), env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if p.returncode != 0:
            print(f"FEHLER beim Erzeugen von {name}:\n{p.stderr.decode('utf-8', 'replace')}")
            return 1
        files = {str(f.relative_to(inputs / name)).replace("\\", "/"): _sha256_file(f)
                 for f in sorted((inputs / name).rglob("*")) if f.is_file()}
        manifest["variants"][name] = {"args": vargs, "files": files}
        print(f"  {name:8s} {len(files):4d} Dateien  {time.perf_counter() - t0:5.1f} s")
    (inputs / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Manifest: {inputs / 'manifest.json'}")
    return 0


def _verify_inputs(inputs: Path) -> Optional[str]:
    mf = inputs / "manifest.json"
    if not mf.is_file():
        return f"Kein Manifest in {inputs} (zuerst make-inputs)."
    manifest = json.loads(mf.read_text(encoding="utf-8"))
    for name, v in manifest["variants"].items():
        for rel, sha in v["files"].items():
            f = inputs / name / rel
            if not f.is_file() or _sha256_file(f) != sha:
                return f"Eingabe veraendert oder fehlend: {name}/{rel}"
    return None


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------

def _case_paths(case_dir: Path) -> dict:
    rs = sorted(case_dir.glob("RS*.dcm"))
    rd = sorted(case_dir.glob("RD*.dcm"))
    ct = case_dir / "CT"
    return {"case": case_dir, "rs": rs[0] if len(rs) == 1 else None, "has_rd": bool(rd),
            "ct": ct if ct.is_dir() else case_dir}


def _expand(arg: str, ctx: dict) -> str:
    out = arg
    for key in ("case", "rs", "ct", "out", "tmp"):
        if "{" + key + "}" in out:
            out = out.replace("{" + key + "}", str(ctx[key]))
    return out


class _JobMeter:
    """
    Peak-Speicher eines Kindprozesses samt Enkeln ueber ein Windows-Job-Objekt.
    Noetig, weil ``python.exe`` einer venv nur ein Launcher ist, der den
    Interpreter als eigenen Prozess startet.  Das Kind startet angehalten
    (CREATE_SUSPENDED), wird dem Job zugeordnet und erst dann fortgesetzt, damit
    kein Enkel am Job vorbei entsteht.  Gemessen werden der Peak Working Set
    (physischer Speicher; ein Thread oeffnet dazu jeden Prozess des Jobs,
    solange er laeuft) und der Peak-Commit.  Der Commit enthaelt die Puffer,
    die OpenBLAS beim Import je Thread zusagt: je rund 30 MB pro logischem Kern
    fuer numpy und fuer scipy (eigene Kopie), bei 24 Kernen gut 1,5 GB schon
    bei kleinen Laeufen; mit OPENBLAS_NUM_THREADS=1 nur rund 40 MB.
    Ausserhalb von Windows oder bei einem Fehler bleibt die Messung leer; der
    Lauf selbst ist davon nie betroffen.
    """

    CREATE_SUSPENDED = 0x00000004
    KILL_ON_JOB_CLOSE = 0x00002000
    BASIC_PROCESS_ID_LIST = 3
    EXTENDED_LIMIT_INFORMATION = 9
    PROCESS_QUERY_LIMITED_INFORMATION = 0x1000

    def __init__(self):
        self.job = None
        self._handles: Dict[int, int] = {}
        self._thread = None
        if os.name != "nt":
            return
        try:
            import ctypes
            from ctypes import wintypes

            class _Basic(ctypes.Structure):
                _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64), ("PerJobUserTimeLimit", ctypes.c_int64),
                            ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                            ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                            ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD),
                            ("SchedulingClass", wintypes.DWORD)]

            class _Extended(ctypes.Structure):
                _fields_ = [("BasicLimitInformation", _Basic), ("IoInfo", ctypes.c_uint64 * 6),
                            ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                            ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]

            k32 = ctypes.WinDLL("kernel32", use_last_error=True)
            k32.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
            k32.CreateJobObjectW.restype = wintypes.HANDLE
            for fn in (k32.SetInformationJobObject, k32.QueryInformationJobObject):
                fn.restype = wintypes.BOOL
            k32.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p,
                                                    wintypes.DWORD]
            k32.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p,
                                                      wintypes.DWORD, ctypes.c_void_p]
            class _PidList(ctypes.Structure):
                _fields_ = [("NumberOfAssignedProcesses", wintypes.DWORD),
                            ("NumberOfProcessIdsInList", wintypes.DWORD), ("ProcessIdList", ctypes.c_size_t * 64)]

            class _Counters(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

            k32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
            k32.AssignProcessToJobObject.restype = wintypes.BOOL
            k32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
            k32.OpenProcess.restype = wintypes.HANDLE
            k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
            k32.K32GetProcessMemoryInfo.restype = wintypes.BOOL
            k32.CloseHandle.argtypes = [wintypes.HANDLE]
            ntdll = ctypes.WinDLL("ntdll")
            ntdll.NtResumeProcess.argtypes = [wintypes.HANDLE]
            ntdll.NtResumeProcess.restype = ctypes.c_long
            self._ct, self._k32, self._ntdll = ctypes, k32, ntdll
            self._Extended, self._PidList, self._Counters = _Extended, _PidList, _Counters
            job = k32.CreateJobObjectW(None, None)
            if job:
                info = _Extended()
                info.BasicLimitInformation.LimitFlags = self.KILL_ON_JOB_CLOSE   # Baum stirbt mit dem Harness
                k32.SetInformationJobObject(job, self.EXTENDED_LIMIT_INFORMATION, ctypes.byref(info),
                                            ctypes.sizeof(info))
                self.job = job
        except Exception:  # noqa: BLE001 - Messung ist optional
            self.job = None

    def popen(self, args: list, **kw):
        if not self.job:
            return subprocess.Popen(args, **kw)
        proc = subprocess.Popen(args, creationflags=self.CREATE_SUSPENDED, **kw)
        handle = int(proc._handle)
        try:
            assigned = bool(self._k32.AssignProcessToJobObject(self.job, handle))
        except Exception:  # noqa: BLE001
            assigned = False
        if self._ntdll.NtResumeProcess(handle) != 0:
            # sollte nie passieren; ohne Messung neu starten
            proc.kill()
            proc.communicate()
            self.close()
            return subprocess.Popen(args, **kw)
        if not assigned:
            self.close()
            return proc
        import threading
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._watch, daemon=True)
        self._thread.start()
        return proc

    def _watch(self) -> None:
        """Jeden Prozess des Jobs oeffnen, solange er laeuft (Handle haelt die Zaehler)."""
        ct, k32 = self._ct, self._k32
        lst = self._PidList()
        while True:
            if k32.QueryInformationJobObject(self.job, self.BASIC_PROCESS_ID_LIST, ct.byref(lst),
                                             ct.sizeof(lst), None):
                for pid in lst.ProcessIdList[:lst.NumberOfProcessIdsInList]:
                    if pid not in self._handles:
                        h = k32.OpenProcess(self.PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
                        if h:
                            self._handles[pid] = h
            if self._stop.wait(0.01):
                return

    def peak_mb(self) -> tuple:
        """(Peak Working Set, Peak-Commit) des groessten Prozesses im Baum, in MB."""
        if not self.job:
            return None, None
        if self._thread is not None:
            self._stop.set()
            self._thread.join()
            self._thread = None
        ws = None
        for h in self._handles.values():
            c = self._Counters()
            c.cb = self._ct.sizeof(c)
            if self._k32.K32GetProcessMemoryInfo(h, self._ct.byref(c), c.cb):
                ws = max(ws or 0.0, c.PeakWorkingSetSize / 2 ** 20)
        info = self._Extended()
        ok = self._k32.QueryInformationJobObject(self.job, self.EXTENDED_LIMIT_INFORMATION,
                                                 self._ct.byref(info), self._ct.sizeof(info), None)
        commit = info.PeakProcessMemoryUsed / 2 ** 20 if ok else None
        return (None if ws is None else round(ws, 1)), (None if commit is None else round(commit, 1))

    def close(self) -> None:
        if self._thread is not None:
            self._stop.set()
            self._thread.join()
            self._thread = None
        for h in self._handles.values():
            self._k32.CloseHandle(h)
        self._handles.clear()
        if self.job:
            self._k32.CloseHandle(self.job)
            self.job = None


def _run_one(sc: dict, sdir: Path, ctx_paths: dict, replacements: list, save_arrays: bool) -> dict:
    import numpy as np

    if sdir.exists():
        shutil.rmtree(sdir)
    out, tmp, cwd = sdir / "out", sdir / "tmp", sdir / "cwd"
    for d in (out, tmp, cwd):
        d.mkdir(parents=True)
    ctx = dict(ctx_paths, out=out, tmp=tmp)
    env = _child_env()
    norm = Normalizer(replacements + [(str(sdir), "<RUN>")])
    steps = sc["steps"] or [sc["args"]]
    step_results = []
    for i, argv_t in enumerate(steps):
        argv = [_expand(a, ctx) for a in argv_t]
        head = ([sys.executable, str(Path(__file__).resolve().with_name(sc["script"]))] if sc.get("script")
                else [sys.executable, "-m", "dicom_file_modifier." + sc["module"]])
        meter = _JobMeter()
        t0 = time.perf_counter()
        try:
            proc = meter.popen(head + argv, cwd=str(cwd), env=env, stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            try:
                so, se = proc.communicate(timeout=sc["timeout"])
                code = proc.returncode
            except subprocess.TimeoutExpired:
                proc.kill()
                so, se = proc.communicate()
                code = "timeout"
            dt = time.perf_counter() - t0
            peak = meter.peak_mb()
        finally:
            meter.close()
        (sdir / f"stdout_{i}.txt").write_bytes(so)
        (sdir / f"stderr_{i}.txt").write_bytes(se)
        step_results.append({
            "argv": [norm.text(a) for a in argv], "exit_code": code, "duration_s": round(dt, 2),
            "peak_ws_mb": peak[0], "peak_commit_mb": peak[1],
            "stdout": norm.text(so.decode("utf-8", "replace")).splitlines(),
            "stderr_tail": norm.text(se.decode("utf-8", "replace")).splitlines()[-25:],
        })
    arrays = {} if save_arrays else None
    files = {}
    for base_name, base in (("out", out), ("tmp", tmp), ("cwd", cwd)):
        for f in sorted(base.rglob("*")):
            if not f.is_file():
                continue
            rel = base_name + "/" + str(f.relative_to(base)).replace("\\", "/")
            key = norm.text(rel)
            files[key] = _file_snapshot(f, norm, arrays, key)
    if arrays:
        np.savez_compressed(str(sdir / "arrays.npz"),
                            **{re.sub(r"[^A-Za-z0-9_]", "_", k): v for k, v in arrays.items()})
    snap = {"scenario": sc["name"], "module": sc["module"], "tags": sc["tags"], "steps": step_results,
            "files": files, "array_keys": sorted(arrays) if arrays else []}
    (sdir / "snapshot.json").write_text(json.dumps(snap, indent=1, ensure_ascii=True), encoding="utf-8")
    return snap


def cmd_run(args) -> int:
    work = Path(args.work)
    inputs = work / "inputs"
    err = _verify_inputs(inputs)
    if err:
        print(err)
        return 2
    run_dir = work / "runs" / args.label
    if run_dir.exists() and not args.only and not args.force:
        print(f"Lauf {args.label} existiert bereits: {run_dir}  (--force zum Ueberschreiben)")
        return 2
    run_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for sc in SCENARIOS:
        if sc["case"] is None:
            ctx = {"case": "", "rs": "", "ct": ""}
        else:
            cp = _case_paths(inputs / sc["case"])
            ctx = {"case": cp["case"], "rs": cp["rs"] or "", "ct": cp["ct"]}
        jobs.append((sc, ctx, [(str(inputs), "<INPUTS>")]))
    for n, case in enumerate(args.real or [], start=1):
        cdir = Path(case).resolve()
        cp = _case_paths(cdir)
        label = f"REAL{n}"
        for rsc in REAL_SCENARIOS:
            need = set(rsc["tags"])
            if ("needs_rs" in need and cp["rs"] is None) or ("needs_rd" in need and not cp["has_rd"]) \
                    or ("needs_ct" in need and cp["ct"] == cdir):
                continue
            sc = dict(rsc, name=f"real{n}_{rsc['name']}", tags=rsc["tags"] + ["real"])
            jobs.append((sc, {"case": cdir, "rs": cp["rs"] or "", "ct": cp["ct"]},
                         [(str(cdir), f"<{label}>"), (cdir.name, label)]))
    if args.only:
        pats = args.only.split(",")
        jobs = [j for j in jobs if any(fnmatch.fnmatch(j[0]["name"], p) for p in pats)]
    meta_path = run_dir / "run.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.is_file() else {}
    meta.update({"label": args.label, "versions": _versions(), "repo_head": _git_head(),
                 "inputs_manifest_sha256": _sha256_file(inputs / "manifest.json")})
    meta.setdefault("scenarios", {})
    t_all = time.perf_counter()
    for sc, ctx, repl in jobs:
        t0 = time.perf_counter()
        snap = _run_one(sc, run_dir / sc["name"], ctx, repl, save_arrays=not args.no_arrays)
        codes = [s["exit_code"] for s in snap["steps"]]
        dur, peak_ws, peak_commit = _step_usage(snap)
        meta["scenarios"][sc["name"]] = {"exit_codes": codes, "tags": sc["tags"], "duration_s": round(dur, 1),
                                         "peak_ws_mb": peak_ws, "peak_commit_mb": peak_commit}
        print(f"  {sc['name']:24s} exit {','.join(str(c) for c in codes):8s} "
              f"{time.perf_counter() - t0:6.1f} s  {len(snap['files']):4d} Dateien"
              + (f"  {peak_ws:7.0f} MB" if peak_ws is not None else ""))
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Lauf {args.label}: {len(jobs)} Szenarien in {time.perf_counter() - t_all:.0f} s -> {run_dir}")
    return 0


def _git_head() -> Optional[str]:
    try:
        p = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=str(REPO),
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        head = p.stdout.decode().strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=str(REPO),
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.strip()
        return head + ("+dirty" if dirty else "")
    except Exception:
        return None


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------

_NUM_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
_INT_RE = re.compile(r"^\s*[-+]?\d+\s*$")


def _decimals(s: str) -> int:
    s = s.strip().lower()
    mant, _, exp = s.partition("e")
    dec = len(mant.split(".")[1]) if "." in mant else 0
    return dec - (int(exp) if exp.lstrip("+-").isdigit() else 0)


def _num(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


class Rules:
    def __init__(self, data: dict, mode: str):
        self.mode = mode
        mig = data.get("migration", {})
        self.float_rules = mig.get("float_rules", [])
        self.pixel_rules = mig.get("pixel_rules", [])
        self.contour = {k: v for k, v in mig.get("contour", {}).items() if not k.startswith("_")}
        self.expected = {k: v for k, v in data.get("expected_changes", {}).items() if not k.startswith("_")}
        # --stats: jede numerische Abweichung (auch innerhalb der Toleranz) als
        # (file, path, a, b, regel, toleranz, ok)
        self.recorder: Optional[list] = None

    def tol_ok(self, a: float, b: float, file: str, path: str, a_str: str = "", b_str: str = "") -> bool:
        if a == b:
            return True
        ok, label, tol = self._judge(a, b, file, path, a_str, b_str)
        if self.recorder is not None:
            self.recorder.append((file, path, float(a), float(b), label, tol, ok))
        return ok

    def _judge(self, a: float, b: float, file: str, path: str, a_str: str, b_str: str) -> tuple:
        """(ok, Regel, Toleranz) fuer zwei verschiedene Zahlen."""
        if self.mode == "exact":
            return False, "exact", 0.0
        key = path.rsplit(".", 1)[-1] if path else ""
        label = "Default (letzte Stelle)"
        for r in self.float_rules:
            if "file" in r and not fnmatch.fnmatch(file, r["file"]):
                continue
            if "key" in r and not fnmatch.fnmatch(key, r["key"]):
                continue
            if "path" in r and not fnmatch.fnmatch(path, r["path"]):
                continue
            if r.get("last_digit"):
                label = _rule_label(r)
                break
            tol = max(float(r.get("abs", 0.0)), float(r.get("rel", 0.0)) * max(abs(a), abs(b)))
            return abs(a - b) <= tol * (1 + 1e-9), _rule_label(r), tol
        # Default: 1 Einheit der letzten gedruckten Stelle; Ganzzahlen (IS,
        # Zaehlwerte in Text/CSV) exakt
        sa = a_str or repr(float(a))
        sb = b_str or repr(float(b))
        if _INT_RE.match(sa) and _INT_RE.match(sb):
            return False, "Ganzzahl (exakt)", 0.0
        tol = 10.0 ** (-max(_decimals(sa), _decimals(sb)))
        return abs(a - b) <= tol * (1 + 1e-9), label, tol

    def pixel_rule(self, scenario: str) -> dict:
        for r in self.pixel_rules:
            if fnmatch.fnmatch(scenario, r.get("scenario", "*")):
                return r
        return {"max_abs": 0, "max_frac": 0.0}


def _rule_label(r: dict) -> str:
    sel = " ".join(f"{k}={r[k]}" for k in ("file", "key", "path") if k in r)
    if r.get("last_digit"):
        return sel + " (letzte Stelle)"
    return sel + " (" + ", ".join(f"{k} {float(r[k]):g}" for k in ("abs", "rel") if k in r) + ")"


class Diff:
    def __init__(self, redact: bool, count: bool = False):
        self.redact = redact
        self.fails: List[str] = []
        self.warns: List[str] = []
        # Pixelstatistik ueber alle CT-Dateien des Szenarios (Migrationsmodus)
        self.pix_total = 0
        self.pix_diff = 0
        self.pix_max = 0
        self.pix_files = 0
        # --stats: numerische Werte in Dateien (auch gleiche) und alle Abweichungen
        self.count = count
        self.n_num = 0
        self.drift: List[tuple] = []

    def fail(self, msg: str):
        self.fails.append(msg)

    def warn(self, msg: str):
        self.warns.append(msg)


def _contour_len_differs(a: dict, b: dict) -> bool:
    ca, cb = a.get("ContourData"), b.get("ContourData")
    return isinstance(ca, list) and isinstance(cb, list) and len(ca) != len(cb)


def _cmp_ring(a: dict, b: dict, rules: Rules, file: str, path: str, diff: Diff) -> None:
    """
    Migrationsmodus: ein Konturring mit anderer Punktzahl (z.B. andere
    Vereinfachung durch eine neue GEOS-Version) wird geometrisch verglichen.
    Gleiche Ebene; Hausdorff <= simplify_mm + extra_mm; ist ein Ring keine
    Vertex-Teilmenge des anderen, zusaetzlich XOR-Flaeche <= xor_frac.
    """
    p = f"{path}.ContourData" if path else "ContourData"
    sa, sb = a["ContourData"], b["ContourData"]
    va, vb = [_num(v) for v in sa], [_num(v) for v in sb]
    if None in va or None in vb or len(va) % 3 or len(vb) % 3 or len(va) < 9 or len(vb) < 9:
        diff.fail(f"{file}: {p} Laenge {len(sa)} vs {len(sb)}")
        return
    for item, vals, lab in ((a, va, "A"), (b, vb, "B")):
        ncp = _num(item.get("NumberOfContourPoints"))
        if ncp is not None and int(ncp) != len(vals) // 3:
            diff.fail(f"{file}: {path}.NumberOfContourPoints passt nicht zu ContourData ({lab})")
    za, zb = {round(v, 4) for v in va[2::3]}, {round(v, 4) for v in vb[2::3]}
    if za != zb or len(za) != 1:
        diff.fail(f"{file}: {p} {len(va) // 3} vs {len(vb) // 3} Punkte, Ebene verschieden oder nicht planar")
        return
    from shapely.geometry import LinearRing, Polygon

    xy_a = [(va[i], va[i + 1]) for i in range(0, len(va), 3)]
    xy_b = [(vb[i], vb[i + 1]) for i in range(0, len(vb), 3)]
    cr = rules.contour
    tol = float(cr.get("simplify_mm", 0.1)) + float(cr.get("extra_mm", 0.01))
    hd = LinearRing(xy_a).hausdorff_distance(LinearRing(xy_b))
    ok = hd <= tol * (1 + 1e-9)
    # Vertex-Teilmenge textuell (gleiche DS-Strings = gleicher Punkt)
    ta = {(sa[i], sa[i + 1]) for i in range(0, len(sa), 3)}
    tb = {(sb[i], sb[i + 1]) for i in range(0, len(sb), 3)}
    subset = tb <= ta or ta <= tb
    note = ""
    if not subset:
        pa, pb = Polygon(xy_a).buffer(0), Polygon(xy_b).buffer(0)
        frac = pa.symmetric_difference(pb).area / max(pa.area, 1e-12)
        xmax = float(cr.get("xor_frac", 0.002))
        ok = ok and frac <= xmax
        note = f", XOR {frac:.2e} der Flaeche (erlaubt {xmax:g})"
    kind = "Vertex-Teilmenge" if subset else "Punkte verschoben"
    if rules.recorder is not None:
        rules.recorder.append((file, p, 0.0, hd, f"ContourData mit anderer Punktzahl ({kind}): Hausdorff mm",
                               tol, ok))
    if not ok:
        diff.fail(f"{file}: {p} {len(va) // 3} vs {len(vb) // 3} Punkte ({kind}): Hausdorff {hd:.4f} mm "
                  f"(erlaubt {tol:g}){note}")


def _cmp_value(a, b, rules: Rules, file: str, path: str, diff: Diff):
    if isinstance(a, dict) and isinstance(b, dict):
        skip = ()
        if rules.mode == "migration" and _contour_len_differs(a, b):
            _cmp_ring(a, b, rules, file, path, diff)
            skip = ("ContourData", "NumberOfContourPoints")
        for k in sorted(set(a) | set(b)):
            if k in skip:
                continue
            p = f"{path}.{k}" if path else k
            if k not in a or k not in b:
                diff.fail(f"{file}: {p} nur in {'B' if k not in a else 'A'}")
                continue
            _cmp_value(a[k], b[k], rules, file, p, diff)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diff.fail(f"{file}: {path} Laenge {len(a)} vs {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            _cmp_value(x, y, rules, file, f"{path}[{i}]" if not path.endswith("]") else f"{path}[{i}]", diff)
        return
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None:
        if a != b:
            diff.fail(f"{file}: {path} {_show(a, diff)} vs {_show(b, diff)}")
        return
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        diff.n_num += 1
        if not rules.tol_ok(float(a), float(b), file, _strip_idx(path)):
            diff.fail(f"{file}: {path} {a!r} vs {b!r} (d={float(b) - float(a):+.3g})")
        return
    if isinstance(a, str) and isinstance(b, str):
        if a == b:
            if diff.count and _num(a) is not None:
                diff.n_num += 1
            return
        na, nb = _num(a), _num(b)
        if na is not None and nb is not None:
            diff.n_num += 1
            if rules.tol_ok(na, nb, file, _strip_idx(path), a, b):
                return
        diff.fail(f"{file}: {path} {_show(a, diff)} vs {_show(b, diff)}")
        return
    if a != b:
        diff.fail(f"{file}: {path} {_show(a, diff)} vs {_show(b, diff)}")


def _strip_idx(path: str) -> str:
    return re.sub(r"\[\d+\]", "", path)


def _show(v, diff: Diff) -> str:
    if diff.redact and isinstance(v, str) and _num(v) is None:
        return "<text>"
    s = repr(v)
    return _ascii(s if len(s) <= 80 else s[:77] + "...")


def _cmp_lines(a: list, b: list, rules: Rules, where: str, diff: Diff, numeric_warn: bool):
    """Zeilenvergleich; im Migrationsmodus Zahlen mit 1 Einheit der letzten Stelle."""
    if len(a) != len(b):
        diff.fail(f"{where}: {len(a)} vs {len(b)} Zeilen")
    for i, (x, y) in enumerate(zip(a, b)):
        if x == y:
            continue
        tx, ty = _NUM_RE.split(x), _NUM_RE.split(y)
        nx, ny = _NUM_RE.findall(x), _NUM_RE.findall(y)
        if tx == ty and len(nx) == len(ny):
            bad = [(u, v) for u, v in zip(nx, ny) if not rules.tol_ok(float(u), float(v), where, "", u, v)]
            if not bad:
                continue
            msg = f"{where} Z.{i + 1}: Zahlen " + ", ".join(f"{u}->{v}" for u, v in bad[:4])
            (diff.warn if numeric_warn else diff.fail)(msg)
            continue
        content = "" if diff.redact else f": {_ascii(x[:70])!r} vs {_ascii(y[:70])!r}"
        diff.fail(f"{where} Z.{i + 1} Text weicht ab{content}")


def _cmp_pixels(fa: dict, fb: dict, name: str, file: str, rules: Rules, arrays: tuple, diff: Diff):
    pa, pb = fa.get("pixels"), fb.get("pixels")
    if pa is None and pb is None:
        return
    if (pa is None) != (pb is None):
        diff.fail(f"{file}: PixelData nur in {'A' if pb is None else 'B'}")
        return
    if "decode_error" in pa or "decode_error" in pb:
        if pa != pb:
            diff.fail(f"{file}: Pixel-Dekodierung {pa.get('decode_error', 'ok')} vs {pb.get('decode_error', 'ok')}")
        return
    if pa["sha256"] == pb["sha256"]:
        n = 1
        for s in pa.get("shape", []):
            n *= int(s)
        diff.pix_total += n          # identische Dateien zaehlen im Nenner mit
        return
    if pa["shape"] != pb["shape"] or pa["dtype"] != pb["dtype"]:
        diff.fail(f"{file}: Pixel {pa['dtype']}{pa['shape']} vs {pb['dtype']}{pb['shape']}")
        return
    if rules.mode == "exact":
        diff.fail(f"{file}: Pixel-Hash verschieden")
        return
    arr_a, arr_b = arrays
    key = re.sub(r"[^A-Za-z0-9_]", "_", file)
    if arr_a is None or arr_b is None or key not in arr_a.files or key not in arr_b.files:
        diff.warn(f"{file}: Pixel-Hash verschieden, keine Arrays fuer Detailvergleich "
                  f"(mean {pa['mean']:.4f} vs {pb['mean']:.4f})")
        return
    import numpy as np
    d = np.abs(arr_a[key].astype(np.int64) - arr_b[key].astype(np.int64))
    diff.pix_total += int(d.size)
    diff.pix_diff += int(np.count_nonzero(d))
    diff.pix_max = max(diff.pix_max, int(d.max()))
    diff.pix_files += 1


def _judge_pixels(name: str, rules: Rules, diff: Diff):
    """Pixeltoleranz ueber alle Dateien des Szenarios (Anteil bezogen aufs Volumen)."""
    if not diff.pix_files:
        return
    rule = rules.pixel_rule(name)
    frac = diff.pix_diff / max(diff.pix_total, 1)
    msg = (f"Pixel in {diff.pix_files} Dateien: |d|max {diff.pix_max}, Anteil {frac:.2e} "
           f"(erlaubt {rule.get('max_abs')}, {rule.get('max_frac')})")
    if diff.pix_max > rule.get("max_abs", 0) or frac > rule.get("max_frac", 0.0):
        diff.fail(msg)
    else:
        diff.warn(msg + " -> innerhalb Toleranz")


def _compare_snap(a: dict, b: dict, rules: Rules, arrays: tuple, stats: bool = False) -> Diff:
    diff = Diff(redact="real" in a.get("tags", []), count=stats)
    rules.recorder = diff.drift if stats else None
    try:
        _compare_snap_into(a, b, rules, arrays, diff)
    finally:
        rules.recorder = None
    return diff


def _compare_snap_into(a: dict, b: dict, rules: Rules, arrays: tuple, diff: Diff) -> None:
    name = a["scenario"]
    if len(a["steps"]) != len(b["steps"]):
        diff.fail(f"Schritte {len(a['steps'])} vs {len(b['steps'])}")
    for i, (sa, sb) in enumerate(zip(a["steps"], b["steps"])):
        if sa["exit_code"] != sb["exit_code"]:
            diff.fail(f"Schritt {i}: Exit {sa['exit_code']} vs {sb['exit_code']}")
        # Zahlen im Konsolentext: im Migrationsmodus nur Warnung, exakt ein Fehler
        _cmp_lines(sa["stdout"], sb["stdout"], rules, f"stdout[{i}]", diff,
                   numeric_warn=rules.mode != "exact")
    fa, fb = a["files"], b["files"]
    for f in sorted(set(fa) | set(fb)):
        if f not in fa or f not in fb:
            diff.fail(f"Datei nur in {'B' if f not in fa else 'A'}: {f}")
            continue
        x, y = fa[f], fb[f]
        if x.get("kind") != y.get("kind"):
            diff.fail(f"{f}: Art {x.get('kind')} vs {y.get('kind')}")
            continue
        kind = x["kind"]
        if kind == "json":
            _cmp_value(x["data"], y["data"], rules, f, "", diff)
        elif kind == "csv":
            if len(x["rows"]) != len(y["rows"]):
                diff.fail(f"{f}: {len(x['rows'])} vs {len(y['rows'])} Zeilen")
            header = x["rows"][0] if x["rows"] else []
            for ri, (ra, rb) in enumerate(zip(x["rows"], y["rows"])):
                if len(ra) != len(rb):
                    diff.fail(f"{f}: Zeile {ri} {len(ra)} vs {len(rb)} Spalten")
                    continue
                for ci, (ca, cb) in enumerate(zip(ra, rb)):
                    col = header[ci] if ci < len(header) else str(ci)
                    _cmp_value(ca, cb, rules, f, col, diff)
        elif kind == "text":
            _cmp_lines(x["lines"], y["lines"], rules, f, diff, numeric_warn=False)
        elif kind == "dicom":
            meta_a = {k: v for k, v in x["file_meta"].items() if k not in _META_IGNORE}
            meta_b = {k: v for k, v in y["file_meta"].items() if k not in _META_IGNORE}
            _cmp_value(meta_a, meta_b, rules, f, "meta", diff)
            _cmp_value(_canon_tag_keys(x["dataset"]), _canon_tag_keys(y["dataset"]), rules, f, "", diff)
            _cmp_pixels(x, y, name, f, rules, arrays, diff)
        elif kind == "html":
            if x.get("external_src") or y.get("external_src"):
                diff.fail(f"{f}: externe src-URL (offline-Pflicht) A={x.get('external_src')} B={y.get('external_src')}")
        elif kind == "png":
            if x.get("dims") != y.get("dims"):
                diff.warn(f"{f}: PNG-Groesse {x.get('dims')} vs {y.get('dims')}")
        elif kind in ("other", "error"):
            if x != y:
                diff.fail(f"{f}: {kind} weicht ab")
    _judge_pixels(name, rules, diff)


def _load_arrays(p: Path):
    if not p.is_file():
        return None
    import numpy as np
    return np.load(str(p))


# Python-Warnungszeile "pfad/datei.py:123: DeprecationWarning: text"
_WARN_RE = re.compile(r"^(?P<src>\S.*?\.py):\d+: (?P<cat>[A-Za-z]*(?:Warning|Error)): (?P<msg>.*)$")


def _stderr_warnings(snap: dict, redact: bool) -> set:
    """Warnungen aus den stderr-Enden; bei echten Faellen ohne Meldungstext."""
    out = set()
    for st in snap.get("steps", []):
        for line in st.get("stderr_tail", []):
            m = _WARN_RE.match(line.strip())
            if m:
                src = re.sub(r"^.*?/site-packages/", "", m.group("src")).replace("<REPO>/", "")
                src = src.replace("numpy/core/", "numpy/_core/")   # numpy 2 hat core -> _core umbenannt
                out.add(f"{src}: {m.group('cat')}" + ("" if redact else ": " + m.group("msg")[:110]))
    return out


def _step_usage(snap: dict) -> tuple:
    """(Sekunden der CLI-Aufrufe, Peak Working Set MB, Peak-Commit MB) ueber alle Schritte."""
    steps = snap.get("steps", [])
    ws = [s["peak_ws_mb"] for s in steps if s.get("peak_ws_mb") is not None]
    cm = [s["peak_commit_mb"] for s in steps if s.get("peak_commit_mb") is not None]
    return (sum(float(s.get("duration_s") or 0.0) for s in steps),
            max(ws) if ws else None, max(cm) if cm else None)


def _fmt(x: float) -> str:
    return "-" if math.isinf(x) else ("0" if x == 0 else f"{x:.3g}")


def _md(s: str) -> str:
    return _ascii(s).replace("|", "\\|")


def _mb(x) -> str:
    return "-" if x is None else f"{x:.0f}"


def _stats_report(results: list) -> tuple:
    """
    Drift-Statistik aus [(szenario, Diff, warnungen_a, warnungen_b, nutzung_a,
    nutzung_b)] -> (Markdown-Zeilen, Konsolenzeilen).
    """
    agg: Dict[tuple, dict] = {}
    per_scen = []
    for name, d, *_ in results:
        n_file = n_stdout = 0
        worst = 0.0
        for file, path, a, b, label, tol, ok in d.drift:
            src = "stdout" if file.startswith("stdout[") else "Datei"
            g = agg.setdefault((src, label), {"n": 0, "over": 0, "abs": 0.0, "rel": -1.0, "ratio": -1.0,
                                              "at": "", "keys": {}, "scen": set()})
            dabs = abs(b - a)
            ratio = dabs / tol if tol > 0 else math.inf
            g["n"] += 1
            g["over"] += 0 if ok else 1
            g["abs"] = max(g["abs"], dabs)
            if min(abs(a), abs(b)) > 0:        # relativ zu 0 ist nichtssagend
                g["rel"] = max(g["rel"], dabs / max(abs(a), abs(b)))
            g["scen"].add(name)
            if ratio > g["ratio"]:
                g["ratio"], g["at"] = ratio, name
            if not d.redact and path:
                k = path.rsplit(".", 1)[-1]
                g["keys"][k] = g["keys"].get(k, 0) + 1
            if src == "stdout":
                n_stdout += 1
            else:
                n_file += 1
            worst = max(worst, ratio)
        if d.drift or d.pix_files:
            per_scen.append((name, d, n_file, n_stdout, worst))

    md = ["## Drift-Statistik", "",
          "Alle numerischen Abweichungen A -> B, auch innerhalb der Toleranz. Quote = |d| / Toleranz "
          "(<= 1 besteht). Schluessel nur aus synthetischen Szenarien.", "",
          "| Quelle | Regel | Werte | ueber Tol. | Szenarien | max abs(d) | max rel(d) | max Quote | bei "
          "| haeufigste Schluessel |",
          "|---|---|---:|---:|---:|---:|---:|---:|---|---|"]
    con = ["Drift-Statistik (Quelle, Regel, Werte, max |d|, max Quote):"]
    for (src, label), g in sorted(agg.items(), key=lambda kv: (kv[0][0], -kv[1]["ratio"], kv[0][1])):
        keys = ", ".join(f"{k} ({c})" for k, c in sorted(g["keys"].items(), key=lambda kc: -kc[1])[:4])
        rel = _fmt(g["rel"]) if g["rel"] >= 0 else "-"
        md.append(f"| {src} | {_md(label)} | {g['n']} | {g['over']} | {len(g['scen'])} | {_fmt(g['abs'])} "
                  f"| {rel} | {_fmt(g['ratio'])} | {g['at']} | {_md(keys)} |")
        con.append(f"  {src:6s} {_ascii(label)[:52]:52s} {g['n']:8d}  {_fmt(g['abs']):>9s}  {_fmt(g['ratio']):>7s}")
    if not agg:
        md.append("| - | keine numerischen Abweichungen | 0 | 0 | 0 | - | - | - | - | - |")
        con.append("  keine numerischen Abweichungen")

    md += ["", "### Je Szenario mit Abweichungen", "",
           "| Szenario | Zahlen in Dateien | davon geaendert | stdout-Zahlen geaendert | max Quote "
           "| CT-Dateien mit Pixeldiff | Pixel geaendert | Anteil | max abs(d) Pixel |",
           "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for name, d, n_file, n_stdout, worst in per_scen:
        frac = d.pix_diff / d.pix_total if d.pix_total else 0.0
        md.append(f"| {name} | {d.n_num} | {n_file} | {n_stdout} | {_fmt(worst) if d.drift else '0'} "
                  f"| {d.pix_files} | {d.pix_diff} | {_fmt(frac)} | {d.pix_max} |")

    md += ["", "### Laufzeit und Peak-Speicher der CLI-Aufrufe (Szenarien ab 5 s oder 500 MB Working Set)", "",
           "WS = Peak Working Set, Commit = Peak-Commit, jeweils des groessten Prozesses im Prozessbaum.", "",
           "| Szenario | A s | B s | B/A | A WS MB | B WS MB | A Commit MB | B Commit MB |",
           "|---|---:|---:|---:|---:|---:|---:|---:|"]
    t_a = t_b = 0.0
    for name, _, _, _, ua, ub in results:
        t_a += ua[0]
        t_b += ub[0]
        if max(ua[0], ub[0]) >= 5 or max(ua[1] or 0, ub[1] or 0) >= 500:
            t_ratio = f"{ub[0] / ua[0]:.2f}" if ua[0] else "-"
            md.append(f"| {name} | {ua[0]:.1f} | {ub[0]:.1f} | {t_ratio} | {_mb(ua[1])} | {_mb(ub[1])} "
                      f"| {_mb(ua[2])} | {_mb(ub[2])} |")
    md.append(f"| Summe ({len(results)} Szenarien) | {t_a:.0f} | {t_b:.0f} | "
              f"{(t_b / t_a if t_a else 0):.2f} | | | | |")
    con.append(f"  Laufzeit der CLI-Aufrufe: A {t_a:.0f} s, B {t_b:.0f} s")

    only_a: Dict[str, int] = {}
    only_b: Dict[str, int] = {}
    both: Dict[str, int] = {}
    for _, _, wa, wb, _, _ in results:
        for w in wa | wb:
            tgt = both if (w in wa and w in wb) else (only_a if w in wa else only_b)
            tgt[w] = tgt.get(w, 0) + 1
    md += ["", "### Warnungen auf stderr (Anzahl Szenarien)", ""]
    if not (only_a or only_b or both):
        md.append("keine")
    for title, dct in (("nur in B", only_b), ("nur in A", only_a), ("in beiden", both)):
        for w, c in sorted(dct.items(), key=lambda wc: (-wc[1], wc[0])):
            md.append(f"- {title}: {c}x `{_ascii(w)}`")
            if title != "in beiden":
                con.append(f"  stderr {title}: {c}x {_ascii(w)[:100]}")
    return md, con


def cmd_compare(args) -> int:
    work = Path(args.work)
    ra, rb = work / "runs" / args.run_a, work / "runs" / args.run_b
    for r in (ra, rb):
        if not r.is_dir():
            print(f"Lauf nicht gefunden: {r}")
            return 2
    rules = Rules(json.loads(RULES_PATH.read_text(encoding="utf-8")), args.mode)
    names = sorted({p.name for p in ra.iterdir() if (p / "snapshot.json").is_file()}
                   | {p.name for p in rb.iterdir() if (p / "snapshot.json").is_file()})
    if args.only:
        pats = args.only.split(",")
        names = [n for n in names if any(fnmatch.fnmatch(n, p) for p in pats)]
    lines = [f"# Golden-Vergleich {args.run_a} vs {args.run_b} (Modus {args.mode})", ""]
    va = json.loads((ra / "run.json").read_text(encoding="utf-8")).get("versions", {}) if (ra / "run.json").is_file() else {}
    vb = json.loads((rb / "run.json").read_text(encoding="utf-8")).get("versions", {}) if (rb / "run.json").is_file() else {}
    for k in sorted(set(va) | set(vb)):
        if va.get(k) != vb.get(k):
            lines.append(f"- {k}: {va.get(k)} -> {vb.get(k)}")
    lines.append("")
    n_fail = n_warn = n_exp = 0
    stats_results = []
    for n in names:
        pa, pb = ra / n / "snapshot.json", rb / n / "snapshot.json"
        if not pa.is_file() or not pb.is_file():
            status, detail = "FAIL", [f"Szenario nur in {'A' if pa.is_file() else 'B'}"]
        else:
            a = json.loads(pa.read_text(encoding="utf-8"))
            b = json.loads(pb.read_text(encoding="utf-8"))
            arrays = (_load_arrays(ra / n / "arrays.npz"), _load_arrays(rb / n / "arrays.npz"))
            d = _compare_snap(a, b, rules, arrays, stats=args.stats)
            if args.stats:
                stats_results.append((n, d, _stderr_warnings(a, d.redact), _stderr_warnings(b, d.redact),
                                      _step_usage(a), _step_usage(b)))
            status = "FAIL" if d.fails else ("WARN" if d.warns else "PASS")
            fails = d.fails
            spec = next((s for pat, s in rules.expected.items() if fnmatch.fnmatch(n, pat)), None)
            head = []
            if status == "FAIL" and args.expected and spec is not None:
                # Eintrag: Text (ganzes Szenario) oder {"reason", "only": [Muster der Abweichungen]}
                reason = spec if isinstance(spec, str) else spec.get("reason", "")
                only = [] if isinstance(spec, str) else list(spec.get("only", []))
                unexpected = [f for f in fails if only and not any(fnmatch.fnmatch(f, p) for p in only)]
                if unexpected:
                    head = [f"erwartete Aenderung ({reason}), dazu {len(unexpected)} UNERWARTETE Abweichungen:"]
                    fails = unexpected
                else:
                    status = "EXPECTED"
                    head = [f"erwartete Aenderung: {reason}"]
            detail = head + fails[:args.max_details] + [f"(Warnung) {w}" for w in d.warns[:args.max_details]]
            if len(fails) > args.max_details:
                detail.append(f"... {len(fails) - args.max_details} weitere Abweichungen")
        n_fail += status == "FAIL"
        n_warn += status == "WARN"
        n_exp += status == "EXPECTED"
        lines.append(f"## {n}: {status}")
        lines.extend(f"- {_ascii(x)}" for x in detail)
        lines.append("")
        print(f"  {status:8s} {n}")
    summary = f"{len(names)} Szenarien: {len(names) - n_fail - n_warn - n_exp} PASS, {n_warn} WARN, " \
              f"{n_exp} EXPECTED, {n_fail} FAIL"
    lines.insert(1, summary)
    if args.stats:
        md, con = _stats_report(stats_results)
        lines += md
        print("\n".join(con))
    report = Path(args.report) if args.report else work / "reports" / f"{args.run_a}_vs_{args.run_b}_{args.mode}.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary)
    print(f"Bericht: {report}")
    return 1 if n_fail else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Golden-Output-Harness (CLI-Paritaet, Plan P0.1)")
    p.add_argument("--work", default=str(DEFAULT_WORK), help=f"Arbeitsordner (Default {DEFAULT_WORK})")
    sub = p.add_subparsers(dest="cmd")
    mi = sub.add_parser("make-inputs", help="Synthetische Eingaben erzeugen (einmalig je Vergleichsreihe)")
    mi.add_argument("--force", action="store_true", help="Vorhandene Eingaben ersetzen")
    rn = sub.add_parser("run", help="Alle Szenarien ausfuehren und Schnappschuesse schreiben")
    rn.add_argument("--label", required=True, help="Name des Laufs, z.B. G0")
    rn.add_argument("--only", help="Komma-getrennte Namensmuster (fnmatch)")
    rn.add_argument("--real", action="append", help="Echter Fallordner (nur lokal; mehrfach moeglich)")
    rn.add_argument("--no-arrays", action="store_true", help="Keine Pixel-Arrays speichern")
    rn.add_argument("--force", action="store_true", help="Vorhandenen Lauf ueberschreiben")
    cp = sub.add_parser("compare", help="Zwei Laeufe vergleichen")
    cp.add_argument("run_a")
    cp.add_argument("run_b")
    cp.add_argument("--mode", choices=("exact", "migration"), default="exact")
    cp.add_argument("--expected", action="store_true",
                    help="Abweichungen in Szenarien aus 'expected_changes' als EXPECTED werten")
    cp.add_argument("--only", help="Komma-getrennte Namensmuster (fnmatch)")
    cp.add_argument("--report", help="Pfad des Markdown-Berichts")
    cp.add_argument("--max-details", type=int, default=12)
    cp.add_argument("--stats", action="store_true",
                    help="Drift-Statistik anhaengen: alle numerischen Abweichungen (auch innerhalb der "
                         "Toleranz), Pixel je Szenario, stderr-Warnungen")
    args = p.parse_args(argv)
    if args.cmd == "make-inputs":
        return cmd_make_inputs(args)
    if args.cmd == "run":
        return cmd_run(args)
    if args.cmd == "compare":
        return cmd_compare(args)
    p.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
