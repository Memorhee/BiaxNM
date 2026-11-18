"""Strain/stress state at the capacity point where a load's ray meets the design curve."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from scipy.optimize import brentq

from ..core.integration import Integrator
from ..core.solver import AngleSolver


@dataclass
class Component:
    name: str
    F: float            # N (compression +)
    x: float            # point of action [mm]
    y: float


@dataclass
class StateInfo:
    """Nominal state at the capacity point of a load (compression positive, N, mm)."""
    integ: Integrator
    theta: float                      # compression direction angle [rad]
    c: float
    a: float
    k: float
    ref: Tuple[float, float]
    lam: float                        # moment direction [rad]
    components: List[Component]
    P: float                          # nominal resultants about ``ref``
    Mx: float
    My: float
    phi: float
    P_design: float
    Mx_design: float
    My_design: float
    capped: bool                      # the design axial force is limited by the axial cap
    exact: bool                       # False when the state is the nearest solved curve point

    @property
    def M(self) -> float:
        return math.hypot(self.Mx, self.My)


def _ray_segment(Md, Pd, Mu, Pu) -> Optional[Tuple[int, float]]:
    """Index and parameter of the nearest crossing of the ray t*(Mu, Pu) with the polyline (Md, Pd)."""
    best = None
    for i in range(len(Md) - 1):
        A = np.array([Md[i], Pd[i]])
        e = np.array([Md[i + 1], Pd[i + 1]]) - A
        den = Mu * e[1] - Pu * e[0]
        if abs(den) < 1e-12 * math.hypot(Mu, Pu) * max(math.hypot(*e), 1e-30):
            continue
        t = (A[0] * e[1] - A[1] * e[0]) / den
        s = (A[0] * Pu - A[1] * Mu) / den
        if -1e-9 <= s <= 1 + 1e-9 and t > 0 and (best is None or t < best[0]):
            best = (t, i, min(max(s, 0.0), 1.0))
    return None if best is None else (best[1], best[2])


def solve_state(result) -> Optional[StateInfo]:
    """State of the section at the capacity point of a :class:`~biaxnm.core.check.CheckResult`.

    The point where the load's ray meets the design curve lies between two solved curve points; the neutral
    axis depth there is found by a one-dimensional root search (phi is interpolated in c).  Returns ``None``
    when the load has no capacity point (DCR infinite, or no curve crossing).
    """
    curve = result.curve
    if result.capacity is None or not math.isfinite(result.dcr) or result.Mu <= 0:
        return None
    section, code = curve.section, curve.code
    if section is None or code is None:
        return None
    hit = _ray_segment(curve.phiMn, curve.phiPn, result.Mu, result.Pu)
    if hit is None:
        return None
    i, s = hit
    n = len(curve.c)
    lam = math.radians(curve.lambda_u_deg)
    integ = Integrator(section, code)
    ref = curve.ref_point
    scale = max(abs(curve.P0), 1.0)
    solver = AngleSolver(integ, lam, ref, m_min=1e-9 * scale * section.diag)
    cap = curve.phiPn_max
    th = np.radians(curve.theta_deg)

    exact = True
    if i == 0 or i + 1 == n - 1:                  # next to a pure-compression / pure-tension end point
        j = i + 1 if i == 0 else i
        c, theta, exact = float(curve.c[j]), float(th[j]), False
        phi = float(curve.phi[j])
    else:
        c_hi, c_lo = float(curve.c[i]), float(curve.c[i + 1])
        phi_hi, phi_lo = float(curve.phi[i]), float(curve.phi[i + 1])

        def at(c_try):
            w = (c_try - c_lo) / (c_hi - c_lo)
            pt = solver.solve(c_try, float(th[i + 1] + w * (th[i] - th[i + 1])))
            if pt is None:
                raise ValueError("no solution")
            phi_try = phi_lo + w * (phi_hi - phi_lo)
            Mn = pt.Mx * math.cos(lam) + pt.My * math.sin(lam)
            Pd = phi_try * pt.P if cap is None else min(phi_try * pt.P, cap)
            return pt, phi_try, phi_try * Mn * result.Pu - Pd * result.Mu

        try:
            c = brentq(lambda v: at(v)[2], c_lo, c_hi, xtol=1e-9 * c_hi, rtol=1e-12, maxiter=100)
            theta = at(c)[0].theta
        except ValueError:                         # no bracket or no solution: use the segment fraction
            c = c_hi + s * (c_lo - c_hi)
            theta = float(th[i] + s * (th[i + 1] - th[i]))
            exact = False
        w = (c - c_lo) / (c_hi - c_lo)
        phi = phi_lo + w * (phi_hi - phi_lo)

    a, k, _ = integ.plane(theta, c)
    comps, P, Qx, Qy = [], 0.0, 0.0, 0.0
    for name, F, qx, qy in integ.components(theta, a, k):
        P += F; Qx += qx; Qy += qy
        if abs(F) > 1e-9 * max(abs(curve.P0), 1.0):
            comps.append(Component(name, F, qy / F, qx / F))
    Mx, My = Qx - ref[1] * P, Qy - ref[0] * P
    Pd = phi * P if cap is None else min(phi * P, cap)
    return StateInfo(integ, theta, c, a, k, ref, lam, comps, P, Mx, My, phi, Pd, phi * Mx, phi * My,
                     cap is not None and phi * P > cap, exact)
