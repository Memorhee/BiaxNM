"""Eurocode 2 (first generation) with National Annex presets, SI units: MPa, mm.

``EC2_1_1``  EN 1992-1-1:2004 (buildings), configuration ``type = "ec2-04-1-1"``
``EC2_2``    EN 1992-2:2005 (concrete bridges), configuration ``type = "ec2-05-2"``

The two parts share the section analysis (EN 1992-2 refers to 3.1.7, 3.2.7 and 6.1 of EN 1992-1-1) and differ in
the nationally determined parameters, checked against the text of EN 1992-2:2005:
    alpha_cc  3.1.6(101)P  recommended 0.85 (EN 1992-1-1 3.1.6(1)P: 1.0), range 0.80 ... 1.00
    fck       3.1.2(102)P  recommended C30/37 ... C70/85 (EN 1992-1-1: C12/15 ... C90/105 of Table 3.1)
    ductility 3.2.4(101)P  recommended reinforcement classes B and C (EN 1992-1-1: A, B, C)
gamma_c, gamma_s (2.4.2.4) and the fyk limit (3.2.2(3)P) are those of EN 1992-1-1 and its National Annex.

Design by material factors, like :mod:`kds_factored`: the section is analysed with design strengths and the result
is the design resistance directly (there is no strength reduction factor phi, ``phi`` stays 1.0).

Concrete (3.1.6, 3.1.7)
    fcd = alpha_cc fck / gamma_c                                  (3.15)
    ``stress_block = "parabola"``     parabola-rectangle, eq. (3.17):  sigma = fcd [1 - (1 - eps/eps_c2)^n]
    ``stress_block = "rectangular"``  stress eta fcd over the depth lambda x (3.1.7(3), eq. 3.19-3.22)
    Table 3.1 (fck 12 ... 90 MPa), above C50:
        eps_c2  = 2.0 + 0.085 (fck - 50)^0.53 per mille          (2.0 per mille up to C50)
        eps_cu2 = 2.6 + 35 ((90 - fck) / 100)^4 per mille         (3.5 per mille up to C50)
        n       = 1.4 + 23.4 ((90 - fck) / 100)^4                 (2.0 up to C50)
        lambda  = 0.8 - (fck - 50) / 400,   eta = 1.0 - (fck - 50) / 200
Reinforcement and tendons (3.2.7, 3.3.6)
    fyd = fyk / gamma_s, Es = 200 000 MPa; ``steel_branch = "horizontal"`` (no strain limit) or ``"inclined"``
    (top branch rising to k fyd at eps_ud, k = ft/fy of the ductility class A/B/C of Annex C, eps_ud = 0.9 eps_uk).
    Structural steel: Fy / gamma_a (gamma_a = 1.0, the EN 1994 value), elastic - perfectly plastic.
    Tendons: ``fpy`` is read as fp0.1k and divided by gamma_s.
Strain plane (6.1, Figure 6.1)
    pivot B: eps_cu2 at the extreme compression fibre while the neutral axis is inside the section;
    pivot C: for a larger depth the plane turns about the point C (strain eps_c2 at the depth h (1 - eps_c2/eps_cu2),
             3/7 h up to C50); uniform eps_c2 for c -> infinity;
    pivot A: with ``steel_branch = "inclined"`` the tension of the most tensioned bar is limited to eps_ud.
    The rectangular block keeps eps_cu3 (= eps_cu2) at the extreme fibre for every c.  Its pure-compression state
    (uniform eps_c2) carries the full fcd, as in kds_factored.
National Annex
    ``NA_PRESETS_1_1`` and ``NA_PRESETS_2`` hold the NDPs this analysis needs: gamma_c, gamma_s, alpha_cc for
    persistent / transient and accidental situations, the upper limit of fyk, eps_ud, the fck range and the
    reinforcement classes.  Part 1-1 presets: "EN" (recommended values), "UK", "DE", "FR", "ES"; part 2 presets:
    "EN", "UK", "FR", "ES", "CZ".  The country values are taken from the Scia summary of the National Annexes
    (see the presets), not from the National Annexes themselves; for any other country use ``na = "custom"`` and
    give the values of that National Annex.
Not applied: minimum eccentricity e0 (6.1(4)), slenderness effects (5.8), confinement (3.1.9), the fyk range 400-600
lower limit (3.2.2(3)), the 0.9 eps_uk value of eps_ud for tendons (the pivot A limit uses the strain plane, which
excludes the prestrain).
Not checked against the original text of EN 1992-1-1: Figure 3.8 end point of the inclined design branch
(k fyk / gamma_s at eps_ud), the Annex C ductility class values below, the accidental-situation factors.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Dict, Optional

import numpy as np

from ..core.laws import FunctionLaw, TabularLaw, elastic_plastic
from .base import DesignCode, Geom, Pivot

ES = 200000.0                           # 3.2.7(4)
DUCTILITY = {"A": (1.05, 0.025), "B": (1.08, 0.050), "C": (1.15, 0.075)}   # Annex C: (k = ft/fy, eps_uk)
EPS_UD_RATIO = 0.9                      # 3.2.7(2) NOTE: recommended eps_ud = 0.9 eps_uk


@dataclass(frozen=True)
class NationalAnnex:
    """The nationally determined parameters used here (EN 1992-1-1 2.4.2.4 Table 2.1N, 3.1.6(1)P)."""
    name: str
    gamma_c: float = 1.5                # persistent and transient
    gamma_s: float = 1.15
    alpha_cc: float = 1.0               # 3.1.6(1)P: between 0.8 and 1.0
    gamma_c_acc: float = 1.2            # accidental
    gamma_s_acc: float = 1.0
    fy_max: float = 600.0               # 3.2.2(3)P: upper limit of fyk
    eps_ud: Optional[float] = None      # 3.2.7(2): None = 0.9 eps_uk
    fck_min: float = 12.0               # strength classes: Table 3.1 (C12/15 ... C90/105), part 2: 3.1.2(102)P
    fck_max: float = 90.0
    ductility_classes: tuple = ("A", "B", "C")     # part 2 3.2.4(101)P: recommended B and C


# Source of the country values: Scia, "Theoretical Background - National Annexes to EN 1992" (2016), which lists the
# NDPs of each National Annex (DIN EN 1992-1-1/NA:2011-01, NF EN 1992-1-1/NA:2007-03, UNE-EN 1992-1-1/AN:2015-10,
# BS EN 1992-2/NA:2007-12, ...).  Not the National Annexes themselves; check against the current edition before using
# a result for a design.
NA_PRESETS_1_1: Dict[str, NationalAnnex] = {
    "EN": NationalAnnex("EN"),                           # recommended values of EN 1992-1-1 (accidental: own memory)
    "UK": NationalAnnex("UK", alpha_cc=0.85),            # alpha_cc 0.85 for flexure and axial load, 1.0 otherwise
    "DE": NationalAnnex("DE", alpha_cc=0.85, gamma_c_acc=1.3, fy_max=500.0, eps_ud=0.025),
    "FR": NationalAnnex("FR", fy_max=500.0),             # all other NDPs: EN values
    "ES": NationalAnnex("ES", gamma_c_acc=1.3, fy_max=500.0),
}

# EN 1992-2: alpha_cc, fck range and the reinforcement classes are part 2 NDPs; gamma_c, gamma_s, the fyk limit and
# eps_ud come from EN 1992-1-1 and the same country's National Annex to it.  Reinforcement classes: the National
# Annexes of these countries use the EN 1992-2 recommendation (B and C).
_B_C = ("B", "C")
NA_PRESETS_2: Dict[str, NationalAnnex] = {
    "EN": NationalAnnex("EN", alpha_cc=0.85, fck_min=30.0, fck_max=70.0, ductility_classes=_B_C),
    "UK": NationalAnnex("UK", alpha_cc=0.85, fck_min=25.0, fck_max=70.0, ductility_classes=_B_C),
    "FR": NationalAnnex("FR", alpha_cc=1.0, fy_max=500.0, fck_min=20.0, fck_max=90.0, ductility_classes=_B_C),
    "ES": NationalAnnex("ES", alpha_cc=1.0, gamma_c_acc=1.3, fy_max=500.0, fck_min=25.0, fck_max=90.0,
                        ductility_classes=_B_C),
    "CZ": NationalAnnex("CZ", alpha_cc=0.90, fck_min=20.0, fck_max=90.0, ductility_classes=_B_C),
}


def eps_c2(fck: float) -> float:
    return 0.002 if fck <= 50.0 else (2.0 + 0.085 * (fck - 50.0) ** 0.53) / 1000.0


def eps_cu2(fck: float) -> float:
    return 0.0035 if fck <= 50.0 else (2.6 + 35.0 * ((90.0 - fck) / 100.0) ** 4) / 1000.0


def exponent_n(fck: float) -> float:
    return 2.0 if fck <= 50.0 else 1.4 + 23.4 * ((90.0 - fck) / 100.0) ** 4


def block_factors(fck: float):
    """``(lambda, eta)`` of the rectangular block, eq. (3.19) - (3.22)."""
    if fck <= 50.0:
        return 0.8, 1.0
    return 0.8 - (fck - 50.0) / 400.0, 1.0 - (fck - 50.0) / 200.0


class EC2Pivot(Pivot):
    """Pivot B (eps_cu at the top fibre), pivot C (point C, optional) and pivot A (tension limit, optional)."""

    def __init__(self, eps_cu: float, eps_co: float, point_c: bool = True, eps_ud: Optional[float] = None):
        self.eps_cu, self.eps_co, self.point_c, self.eps_ud = eps_cu, eps_co, point_c, eps_ud

    def curvature(self, c: float, geom: Geom) -> float:
        h = geom.d_max - geom.d_min
        if self.point_c and c > h:
            k = self.eps_co / (c - h * (1.0 - self.eps_co / self.eps_cu))
        else:
            k = self.eps_cu / c
        if self.eps_ud is not None:
            span = (geom.d_max - geom.d_rein) - c           # lever of the most tensioned bar about the neutral axis
            if span > 0.0:
                k = min(k, self.eps_ud / span)
        return k

    def uniform_strain(self) -> float:
        return self.eps_co


def _positive(who: str, **values) -> None:
    for name, v in values.items():
        if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
            raise ValueError(f"{who}: {name} must be a positive finite number, got {v!r}")


class EC2_1_1(DesignCode):
    """EN 1992-1-1:2004 design resistance of a section under axial force and biaxial bending.

    ``EC2_2`` (EN 1992-2:2005) differs only in the class attributes ``NA_PRESETS`` and ``STANDARD``.

    fck, fy       concrete characteristic cylinder strength (12 ... 90 MPa, National Annex) and rebar yield
                  strength fyk (<= 600 MPa, National Annex)
    Fy            structural steel yield strength (needed when the section has steel)
    fpy, Ep       tendon fp0.1k and modulus (needed when it has tendons)
    na            ``"EN"`` (recommended values, default), ``"UK"``, ``"DE"``, ``"FR"``, ``"ES"`` or ``"custom"``
                  (gamma_c, gamma_s, alpha_cc required); see ``NA_PRESETS`` of the class for the presets
    situation     ``"persistent"`` (persistent and transient, default) or ``"accidental"``
    gamma_c, gamma_s, alpha_cc   override the values of the National Annex (gamma_c, gamma_s apply as given, so an
                  accidental situation with explicit values uses them unchanged)
    stress_block  ``"parabola"`` (3.1.7(1), default) or ``"rectangular"`` (3.1.7(3))
    steel_branch  ``"horizontal"`` (default, no strain limit) or ``"inclined"`` (strain limit eps_ud, pivot A)
    ductility     ``"A"``, ``"B"`` (default) or ``"C"``: k = ft/fy and eps_uk of the inclined branch
    eps_ud        strain limit of the inclined branch; default the National Annex value (DE 0.025) or 0.9 eps_uk
    gamma_a       partial factor of structural steel (default 1.0)
    """

    STANDARD = "EN 1992-1-1:2004"
    NA_PRESETS = NA_PRESETS_1_1

    def __init__(self, fck: float, fy: float, Es: float = ES, Fy: Optional[float] = None,
                 Es_steel: float = 200000.0, fpy: Optional[float] = None, Ep: float = 195000.0,
                 na: str = "EN", situation: str = "persistent", gamma_c: Optional[float] = None,
                 gamma_s: Optional[float] = None, alpha_cc: Optional[float] = None, stress_block: str = "parabola",
                 steel_branch: str = "horizontal", ductility: str = "B", eps_ud: Optional[float] = None,
                 gamma_a: float = 1.0):
        who = type(self).__name__
        if stress_block not in ("parabola", "rectangular"):
            raise ValueError(f"{who}: stress_block must be 'parabola' or 'rectangular'")
        if steel_branch not in ("horizontal", "inclined"):
            raise ValueError(f"{who}: steel_branch must be 'horizontal' or 'inclined'")
        if situation not in ("persistent", "accidental"):
            raise ValueError(f"{who}: situation must be 'persistent' or 'accidental'")
        if ductility not in DUCTILITY:
            raise ValueError(f"{who}: ductility must be one of A, B, C, got {ductility!r}")
        _positive(who, fck=fck, fy=fy, Es=Es, Es_steel=Es_steel, Ep=Ep, gamma_a=gamma_a)
        for name, v in (("Fy", Fy), ("fpy", fpy), ("eps_ud", eps_ud)):
            if v is not None:
                _positive(who, **{name: v})
        presets = self.NA_PRESETS
        key = str(na).upper() if str(na).lower() != "custom" else "custom"
        if key == "custom":
            if gamma_c is None or gamma_s is None or alpha_cc is None:
                raise ValueError(f"{who}: na = 'custom' needs gamma_c, gamma_s and alpha_cc")
            ann = replace(presets["EN"], name="custom", gamma_c=gamma_c, gamma_s=gamma_s, alpha_cc=alpha_cc,
                          gamma_c_acc=gamma_c, gamma_s_acc=gamma_s, fy_max=600.0, eps_ud=None)
        elif key in presets:
            ann = presets[key]
        else:
            raise ValueError(f"{who}: unknown National Annex {na!r}; presets {', '.join(presets)} or 'custom' "
                             "with gamma_c, gamma_s, alpha_cc")
        if not ann.fck_min <= fck <= ann.fck_max:
            raise ValueError(f"{who}: fck must be between {ann.fck_min:g} and {ann.fck_max:g} MPa for the National "
                             f"Annex {ann.name} ({self.STANDARD}), got {fck!r}")
        if ductility not in ann.ductility_classes:
            raise ValueError(f"{who}: reinforcement class {ductility!r} is not allowed for the National Annex "
                             f"{ann.name} ({self.STANDARD}); allowed: {', '.join(ann.ductility_classes)}")
        if fy > ann.fy_max:
            raise ValueError(f"{who}: fy above {ann.fy_max:g} MPa is not valid for the National Annex {ann.name} "
                             f"(3.2.2(3)P), got {fy!r}")
        acc = situation == "accidental"
        self.gamma_c = gamma_c if gamma_c is not None else (ann.gamma_c_acc if acc else ann.gamma_c)
        self.gamma_s = gamma_s if gamma_s is not None else (ann.gamma_s_acc if acc else ann.gamma_s)
        self.alpha_cc = alpha_cc if alpha_cc is not None else ann.alpha_cc
        _positive(who, gamma_c=self.gamma_c, gamma_s=self.gamma_s, alpha_cc=self.alpha_cc)
        if self.gamma_c < 1.0 or self.gamma_s < 1.0:
            raise ValueError(f"{who}: gamma_c and gamma_s must be at least 1.0")
        if not 0.8 <= self.alpha_cc <= 1.0:
            raise ValueError(f"{who}: alpha_cc must be between 0.8 and 1.0 (3.1.6(1)P), got {self.alpha_cc!r}")

        self.fck, self.fy, self.Es = fck, fy, Es
        self.Fy, self.Es_steel, self.fpy, self.Ep = Fy, Es_steel, fpy, Ep
        self.na, self.situation = ann.name, situation
        self.stress_block, self.steel_branch, self.ductility, self.gamma_a = stress_block, steel_branch, ductility, gamma_a
        self.fcd = self.alpha_cc * fck / self.gamma_c
        self.eps_c2, self.eps_cu2, self.n = eps_c2(fck), eps_cu2(fck), exponent_n(fck)
        self.lam, self.eta = block_factors(fck)
        self.k_steel, eps_uk = DUCTILITY[ductility]
        self.eps_ud = None
        if steel_branch == "inclined":
            self.eps_ud = eps_ud if eps_ud is not None else (
                ann.eps_ud if ann.eps_ud is not None else EPS_UD_RATIO * eps_uk)
            if self.eps_ud <= self.gamma_s_yield(fy, Es):
                raise ValueError(f"{who}: eps_ud must exceed the design yield strain "
                                 f"{self.gamma_s_yield(fy, Es):.5f}, got {self.eps_ud!r}")
        self.pivot = EC2Pivot(self.eps_cu2, self.eps_c2, point_c=stress_block == "parabola", eps_ud=self.eps_ud)

    def gamma_s_yield(self, f: float, E: float) -> float:
        return f / self.gamma_s / E

    # ------------------------------------------------------------------ laws
    def concrete_law(self):
        fcd, eco, n = self.fcd, self.eps_c2, self.n
        if self.stress_block == "parabola":
            def law(e):
                x = np.clip(e / eco, 0.0, 1.0)
                return np.where(e > 0.0, fcd * (1.0 - (1.0 - x) ** n), 0.0)
            return FunctionLaw(law, breakpoints=[0.0, eco], n_sub=1 if n == 2.0 else 12, eps_range=(0.0, eco))
        thr, level = (1.0 - self.lam) * self.eps_cu2, self.eta * fcd
        w = 1.0e-9           # the uniform state eps_c2 carries the full fcd; in bending it is a sliver of zero thickness

        def block(e):
            return np.where(np.abs(e - eco) <= w, fcd, np.where(e >= thr, level, 0.0))
        return FunctionLaw(block, breakpoints=[thr, eco - w, eco + w])

    def _reinforcement(self, E: float, f: float):
        fd = f / self.gamma_s
        if self.steel_branch == "horizontal":
            return elastic_plastic(E, fd)
        ey, eu, top = fd / E, self.eps_ud, self.k_steel * fd
        return TabularLaw([-1.0, -eu, -ey, ey, eu, 1.0], [-top, -top, -fd, fd, top, top], beyond="hold",
                          yield_strain=ey)

    def rebar_law(self):
        return self._reinforcement(self.Es, self.fy)

    def steel_law(self):
        if self.Fy is None:
            raise ValueError(f"{type(self).__name__}(..., Fy=...) is required for sections with structural steel")
        return elastic_plastic(self.Es_steel, self.Fy / self.gamma_a)

    def tendon_law(self):
        if self.fpy is None:
            raise ValueError(f"{type(self).__name__}(..., fpy=...) is required for sections with tendons")
        return elastic_plastic(self.Ep, self.fpy / self.gamma_s)


class EC2_2(EC2_1_1):
    """EN 1992-2:2005 (concrete bridges): the part 2 values of alpha_cc (recommended 0.85), the fck range
    (C30/37 ... C70/85) and the reinforcement classes (B, C).  Arguments as for :class:`EC2_1_1`."""
    STANDARD = "EN 1992-2:2005"
    NA_PRESETS = NA_PRESETS_2
