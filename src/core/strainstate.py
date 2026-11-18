"""Strain-state diagnostics for the design code: tensile extremes, yield strains, rupture, tension limit."""
from __future__ import annotations

import math
from typing import Dict, List, Optional

import numpy as np

from .materials import ResolvedMaterials


class StrainReport:
    """Reads a strain plane (theta, a, k) against the resolved materials; never integrates stresses."""

    def __init__(self, materials: ResolvedMaterials):
        self.materials = materials

    def __getattr__(self, name):
        if name == "materials":
            raise AttributeError(name)
        return getattr(self.materials, name)

    def _strain_sets(self, theta, a, k):
        """Yield (kind, law, law_strains, plane_strains) for the concrete and every
        reinforcement group.  ``law_strains`` are the strains the law is evaluated at
        (tendons: minus the effective prestrain); ``plane_strains`` are the strain
        plane values eps(d), which exclude the prestrain (used for eps_t)."""
        c, s = math.cos(theta), math.sin(theta)
        out = []
        dc = self.conc_pts[:, 0] * c + self.conc_pts[:, 1] * s
        e = a + k * dc
        out.append(("concrete", self.conc, e, e))
        if len(self.bar_x):
            d = self.bar_x * c + self.bar_y * s
            for law, idx in self.bar_groups:
                e = a + k * d[idx]
                out.append(("rebar", law, e, e))
        for _, law, poly in self.steel:
            d = poly[:, 0] * c + poly[:, 1] * s
            e = a + k * d
            out.append(("steel", law, e, e))
        if len(self.t_x):
            d = self.t_x * c + self.t_y * s
            for law, idx in self.t_groups:
                e = a + k * d[idx]
                out.append(("tendon", law, e - self.t_pe[idx], e))
        return out

    def tension_extremes(self, theta, a, k) -> Dict[str, Optional[float]]:
        """Largest tensile strain per reinforcement kind (negative if all compressed).

        Prestrain is excluded: this is the net tensile strain eps_t of ACI 318."""
        ext: Dict[str, Optional[float]] = {"rebar": None, "steel": None, "tendon": None}
        for kind, _, _, eps in self._strain_sets(theta, a, k):
            if kind == "concrete":
                continue
            t = float(-eps.min())
            ext[kind] = t if ext[kind] is None else max(ext[kind], t)
        return ext

    def governing_yield(self, theta, a, k) -> Dict[str, Optional[float]]:
        """Yield strain (from the material's own law) of the most tensioned element of each kind.

        ``None`` when that law does not state a yield strain; the design code then uses its default.
        """
        c, s = math.cos(theta), math.sin(theta)
        out: Dict[str, Optional[float]] = {"rebar": None, "steel": None, "tendon": None}
        if len(self.bar_x):
            e = a + k * (self.bar_x * c + self.bar_y * s)
            out["rebar"] = self.bar_laws[int(np.argmin(e))].yield_strain()
        if len(self.t_x):
            e = a + k * (self.t_x * c + self.t_y * s)
            out["tendon"] = self.t_laws[int(np.argmin(e))].yield_strain()
        if self.steel:
            lows = [float((a + k * (poly[:, 0] * c + poly[:, 1] * s)).min()) for _, _, poly in self.steel]
            out["steel"] = self.steel[int(np.argmin(lows))][1].yield_strain()
        return out

    def tension_limit_strain(self) -> float:
        """Uniform strain (negative) of the pure-tension state.

        Every reinforcement element is taken to its law's tensile rupture limit; the least
        extreme limit governs so that no element is pushed past its table.  Laws that never
        rupture put the section at a large tensile strain where they hold their end stress.
        """
        limits = []
        pairs = [(law, 0.0) for law in self.bar_laws] + [(law, 0.0) for _, law, _ in self.steel]
        pairs += [(law, float(pe)) for law, pe in zip(self.t_laws, self.t_pe)]
        for law, pe in pairs:
            lo = law.rupture_limits()[0]
            if lo is not None:
                limits.append(lo + pe)           # tendon: law strain = plane strain - prestrain
        return max(limits) if limits else -0.5

    def rupture_flags(self, theta, a, k) -> List[str]:
        flags = []
        for kind, law, eps, _ in self._strain_sets(theta, a, k):
            lo, hi = law.rupture_limits()
            if lo is not None and eps.min() < lo:
                flags.append(f"{kind}: strain {eps.min():.5f} below the law's table range ({lo:.5f})")
            if hi is not None and eps.max() > hi:
                flags.append(f"{kind}: strain {eps.max():.5f} above the law's table range ({hi:.5f})")
        return flags
