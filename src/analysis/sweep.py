"""Neutral-axis sweep: solve theta for each depth c and collect the strain-state rows."""
from __future__ import annotations

import math
from typing import List, NamedTuple, Optional, Tuple

import numpy as np

from ..core.integration import Integrator
from ..core.solver import AngleSolver, Point, elastic_theta


class Row(NamedTuple):
    """One point of the curve before the design code is applied."""
    c: float
    theta: float
    P: float
    Mx: float
    My: float
    ext: dict            # largest tensile strain per reinforcement kind
    yld: dict            # yield strain of the most tensioned element per kind


def sweep_neutral_axis(solver: AngleSolver, section, lam_u: float, cs: np.ndarray,
                       full_scan: bool = False) -> Tuple[List[Point], List[float]]:
    """Solve every depth in ``cs`` (descending), starting at the one nearest 0.5 D and moving both ways.

    Returns the solved points ordered by descending c and the depths that were left out.
    """
    n_c = len(cs)
    D = section.diag
    i0 = int(np.argmin(np.abs(cs - 0.5 * D)))
    th_start = elastic_theta(section, lam_u)

    solved = {}
    skipped: List[float] = []

    def run(indices):
        # After three failures in a row (typically the P0 plateau, where lambda_n
        # cannot reach lambda_u) the expensive 0..360 deg scan is skipped until a
        # depth is solved again.
        th, fails = th_start, 0
        for i in indices:
            pt = solver.solve(cs[i], th, force_global=full_scan, allow_global=full_scan or fails < 3)
            if pt is None:
                skipped.append(float(cs[i]))
                fails += 1
                continue
            solved[i] = pt
            th, fails = pt.theta, 0

    run(range(i0, n_c))
    run(range(i0 - 1, -1, -1))
    return [solved[i] for i in sorted(solved)], skipped


def rows_from_points(integ: Integrator, points: List[Point]) -> Tuple[List[Row], List[str]]:
    """Rows for the solved points plus the (de-duplicated) table-range warnings."""
    rows: List[Row] = []
    warnings: List[str] = []
    seen = set()
    for p in points:
        ext = integ.tension_extremes(p.theta, p.a, p.k)
        rows.append(Row(p.c, p.theta, p.P, p.Mx, p.My, ext, integ.governing_yield(p.theta, p.a, p.k)))
        for msg in integ.rupture_flags(p.theta, p.a, p.k):
            key = msg.split(":")[0] + (" below" if "below" in msg else " above")
            if key not in seen:
                seen.add(key)
                warnings.append(f"{msg} (first at c = {p.c:.1f} mm)")
    return rows, warnings
