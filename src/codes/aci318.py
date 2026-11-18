"""ACI 318-25, SI units: MPa, mm.

Checked against the code text: Table 21.2.2 (phi), 21.2.2.1/21.2.2.2 (eps_ty), 21.2.2.3 (phi cap), Table 22.2.2.4.3
(beta1, here with the SI-rounded limits 28 and 55 MPa), Table 22.4.2.1 and 22.4.2.2/22.4.2.3 (maximum axial strength, P0).
"""
from __future__ import annotations

import math
from typing import Dict, Optional

import numpy as np

from ..core.laws import elastic_plastic, whitney_step
from .base import DesignCode, TopFiberPivot


def beta1(fc: float) -> float:
    if fc <= 28.0:
        return 0.85
    if fc >= 55.0:
        return 0.65
    return 0.85 - 0.05 * (fc - 28.0) / 7.0


class ACI318(DesignCode):
    """ACI 318-25 rules for compression-controlled / flexural sections.

    fc, fy        concrete strength and rebar yield strength (MPa)
    Fy            structural steel yield strength (needed when the section has steel)
    fpy, Ep       tendon yield strength and modulus (needed when it has tendons)
    transverse    ``"tied"`` or ``"spiral"``
    phi_strain    ``"all"``: eps_t is taken at the most tensioned of rebar, steel
                  and tendon; ``"rebar"``: rebar only (falls back to ``"all"`` if
                  the section has no rebar)
    gamma_c/gamma_s  optional material factors (default 1, as in ACI)
    """

    def __init__(self, fc: float, fy: float, Es: float = 200000.0, Fy: Optional[float] = None,
                 Es_steel: float = 200000.0, fpy: Optional[float] = None, Ep: float = 195000.0,
                 transverse: str = "tied", phi_strain: str = "all", eps_cu: float = 0.003,
                 gamma_c: float = 1.0, gamma_s: float = 1.0):
        if transverse not in ("tied", "spiral"):
            raise ValueError("transverse must be 'tied' or 'spiral'")
        if phi_strain not in ("all", "rebar"):
            raise ValueError("phi_strain must be 'all' or 'rebar'")
        for name, v in (("fc", fc), ("fy", fy), ("Es", Es), ("Es_steel", Es_steel), ("Ep", Ep),
                        ("eps_cu", eps_cu), ("gamma_c", gamma_c), ("gamma_s", gamma_s)):
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
                raise ValueError(f"ACI318: {name} must be a positive finite number, got {v!r}")
        for name, v in (("Fy", Fy), ("fpy", fpy)):
            if v is not None and not (math.isfinite(v) and v > 0):
                raise ValueError(f"ACI318: {name} must be a positive finite number, got {v!r}")
        self.fc, self.fy, self.Es = fc, fy, Es
        self.Fy, self.Es_steel, self.fpy, self.Ep = Fy, Es_steel, fpy, Ep
        self.transverse, self.phi_strain = transverse, phi_strain
        self.gamma_c, self.gamma_s = gamma_c, gamma_s
        self.pivot = TopFiberPivot(eps_cu)
        self.beta1 = beta1(fc)

    # ------------------------------------------------------------------ laws
    def concrete_law(self):
        return whitney_step(self.fc, self.beta1, self.pivot.eps_cu, 0.85).scaled(1.0 / self.gamma_c)

    def rebar_law(self):
        return elastic_plastic(self.Es, self.fy).scaled(1.0 / self.gamma_s)

    def steel_law(self):
        if self.Fy is None:
            raise ValueError("ACI318(..., Fy=...) is required for sections with structural steel")
        return elastic_plastic(self.Es_steel, self.Fy).scaled(1.0 / self.gamma_s)

    def tendon_law(self):
        if self.fpy is None:
            raise ValueError("ACI318(..., fpy=...) is required for sections with tendons")
        return elastic_plastic(self.Ep, self.fpy).scaled(1.0 / self.gamma_s)

    # ------------------------------------------------------------- strength
    @property
    def phi_compression(self) -> float:
        return 0.65 if self.transverse == "tied" else 0.75

    def _yield_strain(self, kind: str, eps_y: Optional[Dict[str, Optional[float]]] = None) -> float:
        """Yield strain used for phi: the material's own law when it states one, else the code's value."""
        own = (eps_y or {}).get(kind)
        if own is not None:
            return own
        if kind == "rebar":
            return self.fy / self.Es
        if kind == "steel":
            # without Fy the yield strain is unknown: be conservative (phi stays at the compression value)
            return math.inf if self.Fy is None else self.Fy / self.Es_steel
        return 0.002  # prestressing steel (ACI 318-25 21.2.2.2)

    def _governing(self, ext, eps_y=None):
        """(eps_t, eps_ty) of the most tensioned reinforcement, or ``None`` when the section has none."""
        cands = [(e, self._yield_strain(k, eps_y)) for k, e in ext.items() if e is not None]
        if self.phi_strain == "rebar" and ext.get("rebar") is not None:
            cands = [(ext["rebar"], self._yield_strain("rebar", eps_y))]
        return max(cands, key=lambda t: t[0]) if cands else None

    def phi(self, ext: Dict[str, Optional[float]], eps_y: Optional[Dict[str, Optional[float]]] = None) -> float:
        gov = self._governing(ext, eps_y)
        phi_c = self.phi_compression
        if gov is None:
            return phi_c
        eps_t, eps_ty = gov
        if eps_t <= eps_ty:
            return phi_c
        if eps_t >= eps_ty + 0.003:
            return 0.90
        return phi_c + (0.90 - phi_c) * (eps_t - eps_ty) / 0.003

    def refine_phi(self, phi, P, states, Ag):
        """ACI 318-25 21.2.2.3: between 0.1 f'c Ag and Pn,bal, phi must not exceed the line from 0.90 at
        0.1 f'c Ag to phi_cc at Pn,bal.  Pn,bal is the nominal axial force where the governing tensile strain
        equals its yield strain (interpolated along the curve).

        The code text applies this to nonprestressed members only.  It is applied to every section here, also
        to those with tendons, which can only lower phi (conservative)."""
        phi = np.array(phi, dtype=float)
        P = np.asarray(P, dtype=float)
        d = np.full(len(states), np.nan)
        for i, s in enumerate(states):
            gov = self._governing(s.strain, s.yield_strain)
            if gov is not None and math.isfinite(gov[1]):
                d[i] = gov[0] - gov[1]
        p_bal = None
        for i in range(len(d) - 1):                       # curve runs from compression to tension
            if d[i] <= 0.0 < d[i + 1]:
                t = -d[i] / (d[i + 1] - d[i])
                p_bal = P[i] + t * (P[i + 1] - P[i])
                break
        p_lo = 0.1 * self.fc * Ag
        if p_bal is None or not p_bal > p_lo:
            return phi
        inside = (P >= p_lo) & (P <= p_bal)
        line = 0.90 + (self.phi_compression - 0.90) * (P - p_lo) / (p_bal - p_lo)
        phi[inside] = np.minimum(phi[inside], line[inside])
        return phi

    def axial_cap(self, P0):
        factor = 0.80 if self.transverse == "tied" else 0.85
        return (factor * P0, self.phi_compression * factor * P0)
