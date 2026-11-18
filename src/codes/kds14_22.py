"""KDS 14 20 20 : 2022 (concrete structures, flexure and compression; strength design), SI units: MPa, mm.

This code also stands for KDS 24 14 20 (bridges, strength design): its flexure and compression rules refer to
KDS 14 20 20 (4.4.3) and its strength reduction factors (4.4.2(2)(2)) equal those of KDS 14 20 10 (4.2.3(2)).

Strength reduction factor phi (KDS 14 20 10 4.2.3(2)): tension-controlled 0.85, compression-controlled 0.70
(spirals) / 0.65 (other), linear in between.  The limits of the net tensile strain eps_t are those of
KDS 14 20 20 4.1.2(3) and (4); eps_t is that of the outermost tendon when the section has tendons and of the
outermost rebar otherwise (a section never consists of structural steel alone).

Concrete (4.1.1):
    eps_cu = 0.0033 - (fck - 40) / 100000, at most 0.0033 (fck above 90 MPa needs test-based values: rejected)
    ``stress_block = "rectangular"``  4.1.1(8): stress eta * 0.85 fck over the depth beta1 * c (Table 4.1-2)
    ``stress_block = "parabola"``     4.1.1(7): parabola-rectangle curve (eq. 4.1-1 ... 4.1-5)
Table values are interpolated linearly between the listed strengths (40 MPa and below: constant).
"""
from __future__ import annotations

import math
from typing import Dict, Optional

import numpy as np

from ..core.laws import FunctionLaw, elastic_plastic
from .base import DesignCode, TopFiberPivot

_FCK_NODES = (40.0, 50.0, 60.0, 70.0, 80.0, 90.0)
_ETA = (1.00, 0.97, 0.95, 0.91, 0.87, 0.84)       # Table 4.1-2
_BETA1 = (0.80, 0.80, 0.76, 0.74, 0.72, 0.70)     # Table 4.1-2
_PS_EPS_TY = 0.002                                # prestressed reinforcement: compression-controlled limit
_EPS_TENSION_CONTROLLED = 0.005                   # 4.1.2(4)
_PHI_TENSION = 0.85                               # KDS 14 20 10 4.2.3(2)
FCK_MAX = 90.0


def eps_cu(fck: float) -> float:
    """Ultimate compressive strain of the concrete, eq. (4.1-5)."""
    return min(0.0033, 0.0033 - (fck - 40.0) / 100000.0)


def eps_co(fck: float) -> float:
    """Strain at the peak stress of the parabola-rectangle curve, eq. (4.1-4)."""
    return max(0.002, 0.002 + (fck - 40.0) / 100000.0)


def exponent_n(fck: float) -> float:
    """Exponent of the rising branch, eq. (4.1-3)."""
    return min(2.0, 1.2 + 1.5 * ((100.0 - fck) / 60.0) ** 4)


def block_factors(fck: float):
    """``(eta, beta1)`` of the equivalent rectangular block (Table 4.1-2)."""
    return (float(np.interp(fck, _FCK_NODES, _ETA)), float(np.interp(fck, _FCK_NODES, _BETA1)))


def parabola_rectangle(fck: float, peak: float) -> FunctionLaw:
    """Parabola-rectangle law (4.1-1, 4.1-2) with the peak stress ``peak`` (0.85 fck, or 0.85 phi_c fck)."""
    eco, n = eps_co(fck), exponent_n(fck)

    def law(e):
        x = np.clip(e / eco, 0.0, 1.0)
        return np.where(e > 0.0, peak * (1.0 - (1.0 - x) ** n), 0.0)
    # for n != 2 the rising branch is not a polynomial of degree 2: refine the integration
    return FunctionLaw(law, breakpoints=[0.0, eco], n_sub=1 if n == 2.0 else 12, eps_range=(0.0, eco))


def check_fck(fck: float, who: str) -> None:
    if fck > FCK_MAX:
        raise ValueError(f"{who}: fck above {FCK_MAX:g} MPa needs test-based strain limits (KDS 14 20 20 4.1.1(3)); "
                         "not supported")


class _KDSPivot(TopFiberPivot):
    """Top-fibre pivot whose pure-compression state is just beyond eps_cu.

    The rectangular block carries eta * 0.85 fck, but the pure axial strength is P0 = 0.85 fck (Ag - Ast) + fy Ast
    (eq. 4.1-16, 4.1-17).  The concrete law returns 0.85 fck for strains beyond eps_cu, which only the pure-compression
    state reaches (every other state has eps_cu at the extreme fibre).
    """
    UNIFORM = 1.0 + 1.0e-4

    def uniform_strain(self):
        return self.eps_cu * self.UNIFORM


class KDS14_22(DesignCode):
    """KDS 14 20 20 : 2022 strength design for sections under flexure and axial force.

    fck, fy       concrete strength and rebar yield strength (MPa); fck up to 90 MPa
    Fy            structural steel yield strength (needed when the section has steel)
    fpy, Ep       tendon yield strength and modulus (needed when it has tendons)
    transverse    ``"tied"`` (phi_c 0.65) or ``"spiral"`` (phi_c 0.70)
    stress_block  ``"rectangular"`` (4.1.1(8), default) or ``"parabola"`` (4.1.1(7))
    phi_max       optional upper limit of phi (e.g. 0.75 for sections near the end of a pretensioned member,
                  KDS 14 20 10 4.2.3(2), without entering transfer and development lengths)
    """

    def __init__(self, fck: float, fy: float, Es: float = 200000.0, Fy: Optional[float] = None,
                 Es_steel: float = 200000.0, fpy: Optional[float] = None, Ep: float = 195000.0,
                 transverse: str = "tied", stress_block: str = "rectangular", phi_max: Optional[float] = None):
        if transverse not in ("tied", "spiral"):
            raise ValueError("transverse must be 'tied' or 'spiral'")
        if stress_block not in ("rectangular", "parabola"):
            raise ValueError("stress_block must be 'rectangular' or 'parabola'")
        for name, v in (("fck", fck), ("fy", fy), ("Es", Es), ("Es_steel", Es_steel), ("Ep", Ep)):
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
                raise ValueError(f"KDS14_22: {name} must be a positive finite number, got {v!r}")
        for name, v in (("Fy", Fy), ("fpy", fpy)):
            if v is not None and not (math.isfinite(v) and v > 0):
                raise ValueError(f"KDS14_22: {name} must be a positive finite number, got {v!r}")
        check_fck(fck, "KDS14_22")
        if phi_max is not None and not (isinstance(phi_max, (int, float)) and 0.0 < phi_max <= _PHI_TENSION):
            raise ValueError(f"KDS14_22: phi_max must be in (0, {_PHI_TENSION}], got {phi_max!r}")
        self.fck, self.fy, self.Es = fck, fy, Es
        self.Fy, self.Es_steel, self.fpy, self.Ep = Fy, Es_steel, fpy, Ep
        self.transverse = transverse
        self.stress_block, self.phi_max = stress_block, phi_max
        self.eta, self.beta1 = block_factors(fck)
        self.pivot = _KDSPivot(eps_cu(fck))

    # ------------------------------------------------------------------ laws
    def concrete_law(self):
        fck, ecu = self.fck, self.pivot.eps_cu
        top = 0.85 * fck
        if self.stress_block == "parabola":
            return parabola_rectangle(fck, top)
        thr = (1.0 - self.beta1) * ecu
        level = self.eta * top
        plateau = ecu * (1.0 + 0.5 * (_KDSPivot.UNIFORM - 1.0))       # only the pure-compression state is beyond

        def block(e):
            return np.where(e > plateau, top, np.where(e >= thr, level, 0.0))
        return FunctionLaw(block, breakpoints=[thr, plateau])

    def rebar_law(self):
        return elastic_plastic(self.Es, self.fy)

    def steel_law(self):
        if self.Fy is None:
            raise ValueError("KDS14_22(..., Fy=...) is required for sections with structural steel")
        return elastic_plastic(self.Es_steel, self.Fy)

    def tendon_law(self):
        if self.fpy is None:
            raise ValueError("KDS14_22(..., fpy=...) is required for sections with tendons")
        return elastic_plastic(self.Ep, self.fpy)

    # ------------------------------------------------------------- strength
    @property
    def phi_compression(self) -> float:
        return 0.65 if self.transverse == "tied" else 0.70

    def _yield_strain(self, kind: str, eps_y: Optional[Dict[str, Optional[float]]] = None) -> float:
        own = (eps_y or {}).get(kind)
        if own is not None:
            return own
        if kind == "rebar":
            return self.fy / self.Es
        if kind == "steel":
            return math.inf if self.Fy is None else self.Fy / self.Es_steel
        return _PS_EPS_TY

    def _tension_controlled_limit(self, kind: str, eps_ty: float) -> float:
        """4.1.2(4): 0.005, or 2.5 times the yield strain when the yield strength exceeds 400 MPa."""
        if kind == "tendon" or not math.isfinite(eps_ty):
            return _EPS_TENSION_CONTROLLED
        fy = eps_ty * (self.Es if kind == "rebar" else self.Es_steel)
        return _EPS_TENSION_CONTROLLED if fy <= 400.0 else 2.5 * eps_ty

    def phi(self, ext: Dict[str, Optional[float]], eps_y: Optional[Dict[str, Optional[float]]] = None) -> float:
        # eps_t: outermost tendon if there are tendons, else outermost rebar (structural steel only as a last resort)
        kind = next((k for k in ("tendon", "rebar", "steel") if ext.get(k) is not None), None)
        phi_c = self.phi_compression
        if kind is None:
            phi = phi_c
        else:
            eps_t = ext[kind]
            eps_ty = self._yield_strain(kind, eps_y)
            eps_tc = self._tension_controlled_limit(kind, eps_ty)
            if eps_t <= eps_ty or not math.isfinite(eps_ty):
                phi = phi_c
            elif eps_t >= eps_tc:
                phi = _PHI_TENSION
            else:
                phi = phi_c + (_PHI_TENSION - phi_c) * (eps_t - eps_ty) / (eps_tc - eps_ty)
        return phi if self.phi_max is None else min(phi, self.phi_max)

    def axial_cap(self, P0):
        """Eq. (4.1-16), (4.1-17): phi Pn,max = 0.85 (spirals) or 0.80 (ties) of phi P0."""
        factor = 0.80 if self.transverse == "tied" else 0.85
        return (factor * P0, self.phi_compression * factor * P0)
