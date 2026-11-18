"""Load check: demand/capacity ratio along a constant-eccentricity ray."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .angles import lambda_u_deg
from .diagram import PMCurve, pm_curve


@dataclass
class CheckResult:
    curve: PMCurve
    Pu: float
    Mux: float
    Muy: float
    Mu: float
    dcr: float
    capacity: Optional[Tuple[float, float]]   # (phiMn, phiPn) where the ray meets the curve, N*mm / N

    @property
    def ok(self) -> bool:
        return self.dcr <= 1.0


def _ray_hit(mp: np.ndarray, d: Tuple[float, float]) -> Optional[float]:
    """Smallest t > 0 with t*d on the polyline mp (columns: M, P)."""
    best = None
    for A, B in zip(mp[:-1], mp[1:]):
        e = B - A
        cross_de = d[0] * e[1] - d[1] * e[0]
        if abs(cross_de) < 1e-12 * (math.hypot(*d) * max(math.hypot(*e), 1e-30)):
            continue
        t = (A[0] * e[1] - A[1] * e[0]) / cross_de
        s = (A[0] * d[1] - A[1] * d[0]) / cross_de
        if -1e-9 <= s <= 1 + 1e-9 and t > 0.0 and (best is None or t < best):
            best = t
    return best


def check_load(curve: PMCurve, Pu: float, Mux: float, Muy: float, tol_deg: float = 0.01) -> CheckResult:
    """DCR for a load whose moment direction equals the curve's lambda_u."""
    Mu = math.hypot(Mux, Muy)
    if Mu > 0:
        diff = (lambda_u_deg(Mux, Muy) - curve.lambda_u_deg + 180.0) % 360.0 - 180.0
        if abs(diff) > tol_deg:
            raise ValueError(f"load direction {lambda_u_deg(Mux, Muy):.4f} deg differs from the curve's "
                             f"lambda_u = {curve.lambda_u_deg:.4f} deg; use check_section() instead")
    if Mu == 0 and Pu == 0:
        return CheckResult(curve, Pu, Mux, Muy, 0.0, 0.0, None)

    mp = np.column_stack([curve.phiMn, curve.phiPn])
    if mp[-1, 0] > 1e-9:                      # close the polygon down the P axis ...
        mp = np.vstack([mp, [0.0, mp[-1, 1]]])
    if mp[0, 0] > 1e-9:                       # ... and up it (moments about a non-plastic reference)
        mp = np.vstack([[0.0, mp[0, 1]], mp])
    t = _ray_hit(mp, (Mu, Pu))
    if t is None:
        # the ray never leaves the capacity region only when the load lies outside it
        # (e.g. tension on a section without any tension capacity)
        return CheckResult(curve, Pu, Mux, Muy, Mu, math.inf, None)
    return CheckResult(curve, Pu, Mux, Muy, Mu, 1.0 / t, (t * Mu, t * Pu))


def check_section(section, code, Pu: float, Mux: float, Muy: float, **kwargs) -> CheckResult:
    """Build the PM curve in the direction of (Mux, Muy) and check the load against it."""
    curve = pm_curve(section, code, lambda_u_deg(Mux, Muy), **kwargs)
    return check_load(curve, Pu, Mux, Muy)
