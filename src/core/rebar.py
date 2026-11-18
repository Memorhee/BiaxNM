"""Standard reinforcing bar sizes: nominal area [mm^2] by designation.

Standards
    "KS"      KS D 3504 deformed bars D4 ... D57 (Table 4).  The nominal diameters are the standard's; the nominal
              area is pi/4 * d^2 of that diameter (D10 = 71.33, D25 = 506.7 mm^2), which agrees with the standard's
              own areas to 0.03 %.
    "ASTM"    ASTM A615 inch-pound designations #3 ... #18 (#3 = 0.11 in^2 ... #18 = 4.00 in^2, 1 in = 25.4 mm).
    "ASTM_SI" ASTM A615M metric designations #10 ... #57 (the same bars as #3 ... #18), with the areas the
              standard states in mm^2 (Table 1; e.g. #16 = 199 mm^2).
    "BS"      Nominal sizes 6 ... 50 mm (also written T16, H16, B16; area = pi/4 * d^2).  6 and 50 mm are
              non-preferred sizes in the standard's table.
    "custom"  no table: give ``area`` or ``dia`` for every bar.

Checked against the tables of KS D 3504 (Table 4), ASTM A615/A615M-12 (Table 1) and the BS area table
(cross-sectional area and mass; areas 28.3 ... 1963.5 mm^2).
"""
from __future__ import annotations

import math
import re
from typing import Dict, Tuple

_KS_DIA = {"D4": 4.23, "D5": 5.29, "D6": 6.35, "D7": 7.00, "D8": 7.94, "D10": 9.53, "D13": 12.7, "D16": 15.9,
           "D19": 19.1, "D22": 22.2, "D25": 25.4, "D29": 28.6, "D32": 31.8, "D35": 34.9, "D38": 38.1, "D41": 41.3,
           "D43": 43.0, "D51": 50.8, "D57": 57.3}

_ASTM_IN2 = {"#3": 0.11, "#4": 0.20, "#5": 0.31, "#6": 0.44, "#7": 0.60, "#8": 0.79, "#9": 1.00, "#10": 1.27,
             "#11": 1.56, "#14": 2.25, "#18": 4.00}
# A615M states the metric areas itself (they are not an exact conversion, e.g. #16 = 199)
_ASTM_SI_MM2 = {"#10": 71, "#13": 129, "#16": 199, "#19": 284, "#22": 387, "#25": 510, "#29": 645,
                "#32": 819, "#36": 1006, "#43": 1452, "#57": 2581}

_BS_DIA = (6, 8, 10, 12, 16, 20, 25, 32, 40, 50)


def _circle(d: float) -> float:
    return math.pi * d * d / 4.0


def _table() -> Dict[str, Dict[str, float]]:
    ks = {name: _circle(d) for name, d in _KS_DIA.items()}
    astm = {name: a * 25.4 ** 2 for name, a in _ASTM_IN2.items()}
    astm_si = {si: float(a) for si, a in _ASTM_SI_MM2.items()}
    bs = {str(d): _circle(d) for d in _BS_DIA}
    return {"KS": ks, "ASTM": astm, "ASTM_SI": astm_si, "BS": bs}


TABLES = _table()
STANDARDS = ("KS", "ASTM", "ASTM_SI", "BS", "custom")


def _normalise(standard: str, size) -> str:
    s = str(size).strip().upper().replace(" ", "")
    if standard == "BS":
        m = re.fullmatch(r"[THRB]?(\d+)", s)           # 16, T16, H16, R16, B16
        return m.group(1) if m else s
    if standard in ("ASTM", "ASTM_SI"):
        return s if s.startswith("#") else "#" + s     # "8" -> "#8"
    return s


def sizes(standard: str):
    """Designations available in a standard, in table order."""
    return list(TABLES[str(standard).upper()])


def rebar_area(standard: str, size) -> float:
    """Nominal area [mm^2] of the bar ``size`` of ``standard`` ("KS", "ASTM", "ASTM_SI", "BS")."""
    std = str(standard).upper()
    if std not in TABLES:
        raise ValueError(f"unknown rebar standard {standard!r} (choose from: {', '.join(STANDARDS)})")
    key = _normalise(std, size)
    if key not in TABLES[std]:
        raise ValueError(f"{std} has no bar {size!r} (available: {', '.join(TABLES[std])})")
    return TABLES[std][key]


def circle_area(dia: float) -> float:
    """Area of a bar of diameter ``dia`` [mm]: pi/4 d^2."""
    if not (isinstance(dia, (int, float)) and math.isfinite(dia) and dia > 0):
        raise ValueError(f"bar diameter must be a positive finite number, got {dia!r}")
    return _circle(float(dia))
