"""KDS design by material factors (no strength reduction factor), SI units: MPa, mm.

``KDS14_22App``  KDS 14 20 20 : 2022, appendix "flexure and compression with material factors"
``KDS24_25``     KDS 24 14 21, concrete bridges, limit state design

Both analyse the section with design material strengths and use the result directly as the design strength:
    concrete   0.85 fck replaced by 0.85 phi_c fck in the parabola-rectangle curve (KDS 14 20 20 4.1.1(7),
               eq. 4.1-1 ... 4.1-5, KDS 24 14 21 3.1.2.5(2)); the appendix may use the equivalent rectangular
               block instead (``stress_block = "rectangular"``, 4.1.1(8), the default there)
    rebar      elastic - perfectly plastic
               KDS14_22App  stress = phi_s min(Es eps, fy)       (appendix eq. 3-2, 3-3: the factor multiplies the
                            whole curve, the yield strain stays fy / Es)
               KDS24_25     stress = min(Es eps, phi_s fy)       (fyd = phi_s fy, 3.2-1: modulus Es unchanged)
    tendon, steel   the same with fpy, Ep and Fy, Es_steel
    factors    ultimate load combinations phi_c = 0.65, phi_s = 0.90 (KDS 24 14 21 Table 1.4-1, appendix 2.2(1)); the
               other combinations of KDS 24 14 21 (extreme, service, fatigue) use 1.0
Pure compression is the state with the uniform strain eps_co, which gives P0 = phi_c 0.85 fck (Ag - Ast) +
phi_s fy Ast (appendix eq. 3-2) or, when the rebar has not yielded at eps_co, phi_s eps_co Es Ast (eq. 3-3).  There is
no 0.80 / 0.85 cap; the minimum eccentricity of the load is not applied here.

Strain plane (parabola): eps_cu at the extreme compression fibre while the neutral axis is inside the section
(c <= h).  For a larger c the plane turns about the point C at the depth h (1 - eps_co / eps_cu) below the extreme
fibre, where the strain stays eps_co (KDS 24 14 21 4.1.1.2(2), Figure 4.1-1); it is uniform eps_co for c -> infinity.
Strain plane (rectangular block): eps_cu stays at the extreme fibre for every c, as in KDS 14 20 20 4.1.1(8).
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

from ..core.laws import FunctionLaw, elastic_plastic
from .base import DesignCode, Geom, Pivot, TopFiberPivot
from .kds14_22 import block_factors, check_fck, eps_co, eps_cu, parabola_rectangle

# material factors (phi_c, phi_s) of the load combinations of KDS 24 14 21 Table 1.4-1
COMBINATIONS = {"ultimate": (0.65, 0.90), "extreme": (1.0, 1.0), "service": (1.0, 1.0), "fatigue": (1.0, 1.0)}
FY_MAX_BRIDGE = 600.0          # KDS 24 14 21 3.2.2(2): rebar strength above this is not valid


class CPointPivot(Pivot):
    """Top-fibre strain eps_cu while c <= h, rotation about the point C (strain eps_co) for c > h."""

    def __init__(self, eps_cu: float, eps_co: float):
        self.eps_cu, self.eps_co = eps_cu, eps_co

    def curvature(self, c: float, geom: Geom) -> float:
        h = geom.d_max - geom.d_min
        if c <= h:
            return self.eps_cu / c
        y_c = h * (1.0 - self.eps_co / self.eps_cu)       # depth of the point C below the extreme fibre
        return self.eps_co / (c - y_c)

    def uniform_strain(self) -> float:
        return self.eps_co


class _BlockPivot(TopFiberPivot):
    """Top-fibre strain eps_cu for every c; the pure-compression state is the uniform strain eps_co."""

    def __init__(self, eps_cu: float, eps_co: float):
        super().__init__(eps_cu)
        self.eps_co = eps_co

    def uniform_strain(self) -> float:
        return self.eps_co


class _KDSFactored(DesignCode):
    """Common implementation of the material-factor design (see the module docstring).

    fck, fy       concrete strength and rebar yield strength (MPa); fck up to 90 MPa
    Fy            structural steel yield strength (needed when the section has steel)
    fpy, Ep       tendon yield strength and modulus (needed when it has tendons)
    phi_c, phi_s  material factors of the concrete and of rebar, steel and tendons
    stress_block  ``"parabola"`` or ``"rectangular"`` (4.1.1(8))
    """

    scale_modulus = False                  # True: the factor phi_s multiplies the whole elastic branch too

    def __init__(self, fck: float, fy: float, Es: float = 200000.0, Fy: Optional[float] = None,
                 Es_steel: float = 200000.0, fpy: Optional[float] = None, Ep: float = 195000.0,
                 phi_c: float = 0.65, phi_s: float = 0.90, stress_block: str = "parabola"):
        who = type(self).__name__
        if stress_block not in ("parabola", "rectangular"):
            raise ValueError(f"{who}: stress_block must be 'parabola' or 'rectangular'")
        for name, v in (("fck", fck), ("fy", fy), ("Es", Es), ("Es_steel", Es_steel), ("Ep", Ep),
                        ("phi_c", phi_c), ("phi_s", phi_s)):
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
                raise ValueError(f"{who}: {name} must be a positive finite number, got {v!r}")
        for name, v in (("Fy", Fy), ("fpy", fpy)):
            if v is not None and not (math.isfinite(v) and v > 0):
                raise ValueError(f"{who}: {name} must be a positive finite number, got {v!r}")
        for name, v in (("phi_c", phi_c), ("phi_s", phi_s)):
            if v > 1.0:
                raise ValueError(f"{who}: {name} must not exceed 1, got {v!r}")
        check_fck(fck, who)
        self.fck, self.fy, self.Es = fck, fy, Es
        self.Fy, self.Es_steel, self.fpy, self.Ep = Fy, Es_steel, fpy, Ep
        self.phi_c, self.phi_s = phi_c, phi_s
        self.stress_block = stress_block
        self.eta, self.beta1 = block_factors(fck)
        if stress_block == "rectangular":
            self.pivot = _BlockPivot(eps_cu(fck), eps_co(fck))
        else:
            self.pivot = CPointPivot(eps_cu(fck), eps_co(fck))

    def concrete_law(self):
        top = 0.85 * self.phi_c * self.fck
        if self.stress_block == "parabola":
            return parabola_rectangle(self.fck, top)
        ecu, eco = self.pivot.eps_cu, self.pivot.eps_co
        thr, level = (1.0 - self.beta1) * ecu, self.eta * top
        # Pure compression (uniform strain eco) carries the full 0.85 phi_c fck (eq. 3-2): a vanishing window of
        # strain around eco gives that state the full stress; in bending it is a sliver of zero thickness.
        w = 1.0e-9

        def block(e):
            return np.where(np.abs(e - eco) <= w, top, np.where(e >= thr, level, 0.0))
        return FunctionLaw(block, breakpoints=[thr, eco - w, eco + w])

    def _steel_curve(self, E: float, f: float):
        if self.scale_modulus:
            return elastic_plastic(self.phi_s * E, self.phi_s * f)
        return elastic_plastic(E, self.phi_s * f)

    def rebar_law(self):
        return self._steel_curve(self.Es, self.fy)

    def steel_law(self):
        if self.Fy is None:
            raise ValueError(f"{type(self).__name__}(..., Fy=...) is required for sections with structural steel")
        return self._steel_curve(self.Es_steel, self.Fy)

    def tendon_law(self):
        if self.fpy is None:
            raise ValueError(f"{type(self).__name__}(..., fpy=...) is required for sections with tendons")
        return self._steel_curve(self.Ep, self.fpy)


class KDS14_22App(_KDSFactored):
    """KDS 14 20 20 : 2022, appendix: flexure and compression with material factors.

    stress_block  ``"rectangular"`` (default, 4.1.1(8)) or ``"parabola"`` (4.1.1(7))
    """
    scale_modulus = True

    def __init__(self, fck: float, fy: float, Es: float = 200000.0, Fy: Optional[float] = None,
                 Es_steel: float = 200000.0, fpy: Optional[float] = None, Ep: float = 195000.0,
                 phi_c: float = 0.65, phi_s: float = 0.90, stress_block: str = "rectangular"):
        super().__init__(fck, fy, Es, Fy, Es_steel, fpy, Ep, phi_c, phi_s, stress_block)

    def axial_limit(self, section):
        """Design axial strength P_do of the appendix 3.2(1), eq. (3-2) / (3-3) (no prestress):

            P_do = phi_c 0.85 fck (Ag - Ast) + sum A_i s_i,   s_i = phi_s fy_i   if phi_s fy_i / E_i <= eps_co
                                                              s_i = phi_s eps_co E_i   otherwise

        Ast is the area of rebar and structural steel.  ``phi_c 0.85 fck`` times the section's compressive strength
        is also the cap of phi Pn (no further factor).  Sections with tendons (3.2(1)③) or with elements that
        carry their own stress-strain law return ``None``: the pure-compression point of the analysis is used.
        """
        from ..section.ringgeom import signed_area
        if section.tendons or any(b.law is not None for b in section.bars) \
                or any(s.law is not None for s in section.steel):
            return None
        eco = self.pivot.eps_co

        def stress(E, f):
            return self.phi_s * f if self.phi_s * f / E <= eco else self.phi_s * eco * E
        A_el = F = 0.0
        for b in section.bars:
            A_el += b.area
            F += b.area * stress(self.Es, self.fy)
        for s in section.steel:
            if self.Fy is None:
                raise ValueError("KDS14_22App(..., Fy=...) is required for sections with structural steel")
            a = abs(signed_area(np.asarray(s.polygon, float))) \
                - sum(abs(signed_area(np.asarray(h, float))) for h in s.holes)
            A_el += a
            F += a * stress(self.Es_steel, self.Fy)
        return self.phi_c * 0.85 * self.fck * (section.area - A_el) + F

    def axial_cap(self, P0):
        """phi Pn must not exceed P_do (the factors are in the material strengths already)."""
        return (P0, P0)


class KDS24_25(_KDSFactored):
    """KDS 24 14 21: concrete bridges, limit state design (parabola-rectangle curve only).

    combination   ``"ultimate"`` (default, phi_c 0.65 and phi_s 0.90), ``"extreme"``, ``"service"`` or
                  ``"fatigue"`` (both factors 1.0); phi_c and phi_s given explicitly override the combination,
                  e.g. for the increased factors allowed under good quality control (1.4.2.3(3))
    fy            at most 600 MPa (3.2.2(2))
    """

    def __init__(self, fck: float, fy: float, Es: float = 200000.0, Fy: Optional[float] = None,
                 Es_steel: float = 200000.0, fpy: Optional[float] = None, Ep: float = 195000.0,
                 combination: str = "ultimate", phi_c: Optional[float] = None, phi_s: Optional[float] = None):
        if combination not in COMBINATIONS:
            raise ValueError(f"KDS24_25: combination must be one of {', '.join(COMBINATIONS)}, got {combination!r}")
        if isinstance(fy, (int, float)) and math.isfinite(fy) and fy > FY_MAX_BRIDGE:
            raise ValueError(f"KDS24_25: fy above {FY_MAX_BRIDGE:g} MPa is not valid (KDS 24 14 21 3.2.2(2)), "
                             f"got {fy!r}")
        c, s = COMBINATIONS[combination]
        super().__init__(fck, fy, Es, Fy, Es_steel, fpy, Ep, c if phi_c is None else phi_c,
                         s if phi_s is None else phi_s, "parabola")
        self.combination = combination
