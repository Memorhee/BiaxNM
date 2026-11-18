"""Apply the design code (phi, axial cap) to the nominal rows and build the arrays of a curve."""
from __future__ import annotations

import math
from typing import Dict, List

import numpy as np

from ..codes.base import DesignCode, TensionState
from .sweep import Row


def apply_design(rows: List[Row], lam_u: float, code: DesignCode, P0: float, Ag: float = 0.0) -> Dict[str, object]:
    """Return the array fields of a :class:`~biaxnm.diagram.PMCurve` for ``rows`` (top ... bottom)."""
    c_arr = np.array([r.c for r in rows])
    th_arr = np.degrees(np.array([r.theta for r in rows]))
    P_unc = np.array([r.P for r in rows])
    Mx = np.array([r.Mx for r in rows])
    My = np.array([r.My for r in rows])
    # end points: take the moment component along lambda_u (never negative)
    for j in (0, len(rows) - 1):
        m = Mx[j] * math.cos(lam_u) + My[j] * math.sin(lam_u)
        m = max(m, 0.0)
        Mx[j], My[j] = m * math.cos(lam_u), m * math.sin(lam_u)
    Mn = np.hypot(Mx, My)

    eps_t = np.array([max([v for v in r.ext.values() if v is not None], default=math.nan) for r in rows])
    eps_t_display = eps_t.copy()
    eps_t_display[-1] = math.nan       # the tension end point has no meaningful strain
    states = [TensionState(r.ext, r.yld) for r in rows]
    phi = np.array([code.phi_for(s) for s in states], dtype=float)
    if np.all(np.isfinite(phi)):
        phi = np.asarray(code.refine_phi(phi, P_unc, states, Ag), dtype=float)
    if not np.all(np.isfinite(phi)):
        raise ValueError("the design code returned a non-finite strength reduction factor phi")
    if not (np.all(np.isfinite(P_unc)) and np.all(np.isfinite(Mx)) and np.all(np.isfinite(My))):
        raise ValueError("non-finite resultants: check the material laws")

    Pn_max, phiPn_max = code.axial_cap(P0)
    Pn = np.minimum(P_unc, Pn_max) if Pn_max is not None else P_unc.copy()
    phiPn = phi * P_unc
    if phiPn_max is not None:
        phiPn = np.minimum(phiPn, phiPn_max)
    phiMn = phi * Mn
    return dict(c=c_arr, theta_deg=th_arr, Pn=Pn, Pn_uncapped=P_unc, Mx=Mx, My=My, Mn=Mn,
                eps_t=eps_t_display, phi=phi, phiPn=phiPn, phiMn=phiMn, Pn_max=Pn_max, phiPn_max=phiPn_max)
