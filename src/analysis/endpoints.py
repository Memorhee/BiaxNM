"""End points of the curve: pure compression (top) and pure tension (bottom)."""
from __future__ import annotations

import math
from typing import Tuple

from ..core.integration import Integrator
from .sweep import Row


def top_row(integ: Integrator, ref: Tuple[float, float]) -> Row:
    """Uniform compression at the pivot's uniform strain."""
    eps_u = integ.code.pivot.uniform_strain()
    P0, Qx, Qy = integ.uniform()
    return Row(math.inf, math.nan, P0, Qx - ref[1] * P0, Qy - ref[0] * P0,
               integ.tension_extremes(0.0, eps_u, 0.0), integ.governing_yield(0.0, eps_u, 0.0))


def bottom_row(integ: Integrator, ref: Tuple[float, float]) -> Row:
    """Pure tension: every reinforcement element at its tensile limit (uniform strain)."""
    e_b = integ.tension_limit_strain()
    Pb, Qxb, Qyb = integ.resultants(0.0, e_b, 0.0)
    return Row(0.0, math.nan, Pb, Qxb - ref[1] * Pb, Qyb - ref[0] * Pb,
               integ.tension_extremes(0.0, e_b, 0.0), integ.governing_yield(0.0, e_b, 0.0))
