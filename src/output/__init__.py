"""Plot and CSV output (kN, kN*m)."""
from __future__ import annotations

import os

from .plot import plot_curve
from .state import plot_state
from .tables import TABLE_HEADER, curve_table, save_csv, save_geometry


def save_outputs(obj, outdir, basename="pm"):
    """Write ``<basename>.png`` and ``<basename>.csv`` into ``outdir``; return the two paths."""
    curve = getattr(obj, "curve", obj)
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{basename}.png")
    csv_path = os.path.join(outdir, f"{basename}.csv")
    plot_curve(obj, png)
    save_csv(curve, csv_path)
    section = getattr(curve, "section", None)
    if section is not None:                            # section input: outline, bars, tendons, steel
        save_geometry(section, os.path.join(outdir, basename))
    if hasattr(obj, "capacity"):                       # a load check: also draw the state at its capacity point
        plot_state(obj, os.path.join(outdir, f"{basename}_state.png"), title=basename)
    return png, csv_path


__all__ = ["plot_curve", "save_csv", "save_geometry", "save_outputs", "curve_table", "TABLE_HEADER"]
