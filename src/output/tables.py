"""Plain-text tables of a PM curve (kN, kN*m)."""
from __future__ import annotations

import csv

import numpy as np

TABLE_HEADER = ["c_mm", "theta_deg", "Pn_kN", "Mnx_kNm", "Mny_kNm", "Mn_kNm", "eps_t", "phi",
                "phiPn_kN", "phiMn_kNm"]


def curve_table(curve) -> np.ndarray:
    """Rows of the curve in the units of ``TABLE_HEADER``."""
    return np.column_stack([curve.c, curve.theta_deg, curve.Pn / 1e3, curve.Mx / 1e6, curve.My / 1e6,
                            curve.Mn / 1e6, curve.eps_t, curve.phi, curve.phiPn / 1e3, curve.phiMn / 1e6])


def save_csv(curve, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(TABLE_HEADER)
        for row in curve_table(curve):
            w.writerow([f"{v:.6g}" for v in row])


def _write(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([v if isinstance(v, str) else f"{v:.6g}" for v in r])


def _ring_rows(name, ring):
    return [(name, i, float(x), float(y)) for i, (x, y) in enumerate(ring)]


def save_geometry(section, prefix) -> list:
    """Write the section input as CSV files ``<prefix>_<part>.csv`` (mm, mm^2); return the paths written.

    ``section`` (concrete outline, holes as ring/vertex/x/y), ``bars`` and ``tendons`` (x, y, area),
    ``steel`` (structural-steel polygons and holes, same layout as ``section``). Empty parts are skipped.
    """
    out = []
    rows = _ring_rows("outer", section.outer)
    for i, h in enumerate(section.holes):
        rows += _ring_rows(f"hole{i + 1}", h)
    parts = [("section", ["ring", "vertex", "x_mm", "y_mm"], rows)]
    parts.append(("bars", ["id", "x_mm", "y_mm", "area_mm2", "dia_mm"],
                  [(i + 1, b.x, b.y, b.area, (4 * b.area / 3.141592653589793) ** 0.5)
                   for i, b in enumerate(section.bars)]))
    parts.append(("tendons", ["id", "x_mm", "y_mm", "area_mm2", "fse_MPa"],
                  [(i + 1, t.x, t.y, t.area, t.fse) for i, t in enumerate(section.tendons)]))
    srows = []
    for i, s in enumerate(section.steel):
        srows += _ring_rows(f"steel{i + 1}", s.polygon)
        for j, h in enumerate(s.holes):
            srows += _ring_rows(f"steel{i + 1}_hole{j + 1}", h)
    parts.append(("steel", ["ring", "vertex", "x_mm", "y_mm"], srows))
    for name, header, r in parts:
        if r:
            path = f"{prefix}_{name}.csv"
            _write(path, header, r)
            out.append(path)
    return out
