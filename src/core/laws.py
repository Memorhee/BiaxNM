"""Stress-strain laws.

Sign convention: compression is positive for both strain and stress.
Every law exposes:

    stress(eps)         vectorised stress for an array of strains
    breakpoints()       strains where the slope changes or the stress jumps
                        (the integrator splits its strips at these strains)
    scaled(sf, ef)      copy with stress multiplied by ``sf`` and strain by ``ef``
    rupture_limits()    (low, high) strains beyond which the law drops to zero,
                        ``None`` for a side that never ruptures
    initial_modulus()   slope of the law at zero strain
"""
from __future__ import annotations

import csv
from typing import Callable, Iterable, Sequence

import numpy as np


class StressStrainLaw:
    def stress(self, eps):  # pragma: no cover - interface
        raise NotImplementedError

    def breakpoints(self) -> np.ndarray:
        return np.empty(0)

    def scaled(self, stress_factor: float = 1.0, strain_factor: float = 1.0):
        raise NotImplementedError

    def rupture_limits(self):
        return (None, None)

    def initial_modulus(self) -> float:
        d = 1e-7
        return float((self.stress(np.array([d]))[0] - self.stress(np.array([-d]))[0]) / (2 * d))

    def yield_strain(self):
        """Tensile yield strain of the material, or ``None`` when the law does not say."""
        return None

    def prestrain_for(self, fse: float) -> float:
        """Smallest tensile strain e >= 0 at which the law's tensile stress reaches ``fse`` (MPa).

        This is the strain that a tendon stressed to ``fse`` carries before the section is
        loaded.  It is found from the tension branch itself, so one-sided laws and laws with
        a flat spot at zero strain are handled; it equals fse/E for an elastic branch.
        """
        if not np.isfinite(fse) or fse < 0:
            raise ValueError(f"effective prestress must be finite and >= 0, got {fse}")
        if fse == 0:
            return 0.0

        def tension(e):
            return -float(self.stress(np.array([-e]))[0])

        hi = 1e-6
        while tension(hi) < fse * (1 - 1e-9):
            hi *= 2.0
            if hi > 1.0:
                raise ValueError(f"effective prestress {fse} MPa is above the tensile strength of the tendon law")
        grid = np.linspace(0.0, hi, 4001)
        vals = np.array([tension(e) for e in grid])
        first = int(np.argmax(vals >= fse * (1 - 1e-12)))
        if first == 0:
            return 0.0
        lo, up = float(grid[first - 1]), float(grid[first])
        for _ in range(100):
            mid = 0.5 * (lo + up)
            if tension(mid) >= fse:
                up = mid
            else:
                lo = mid
        return 0.5 * (lo + up)


def _as_pair(beyond):
    if beyond is None:
        return (None, None)
    if isinstance(beyond, str):
        return (beyond, beyond)
    low, high = beyond
    return (low, high)


class TabularLaw(StressStrainLaw):
    """Piecewise-linear law from (strain, stress) points.

    ``beyond`` controls the behaviour outside the table: ``"hold"`` keeps the end
    stress, ``"zero"`` drops to zero (rupture / crushing).  A single string applies
    to both sides, a pair ``(low, high)`` sets the tension-side (low strain) and
    compression-side (high strain) behaviour.  ``None`` means "unspecified"; the
    design code then fills in its own default via :meth:`with_default_beyond`.
    """

    def __init__(self, strains: Sequence[float], stresses: Sequence[float], beyond=None, yield_strain=None):
        e = np.asarray(strains, dtype=float)
        s = np.asarray(stresses, dtype=float)
        if e.ndim != 1 or e.shape != s.shape or e.size < 2:
            raise ValueError("strains and stresses must be 1-D arrays of equal length >= 2")
        if not (np.all(np.isfinite(e)) and np.all(np.isfinite(s))):
            raise ValueError("strains/stresses contain NaN or inf")
        if np.any(np.diff(e) <= 0):
            raise ValueError("strains must be strictly increasing")
        lo, hi = _as_pair(beyond)
        for b in (lo, hi):
            if b not in (None, "hold", "zero"):
                raise ValueError("beyond must be 'hold' or 'zero'")
        self.e, self.s = e, s
        self.beyond = (lo, hi)
        self._eps_y = None if yield_strain is None else float(yield_strain)

    def yield_strain(self):
        return self._eps_y

    @classmethod
    def from_csv(cls, path, beyond=None, delimiter=",", skip_header=None):
        rows = []
        with open(path, newline="", encoding="utf-8-sig") as f:
            for row in csv.reader(f, delimiter=delimiter):
                if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                    continue
                try:
                    rows.append((float(row[0]), float(row[1])))
                except ValueError:
                    if skip_header is False:
                        raise
                    continue  # header line
        if len(rows) < 2:
            raise ValueError(f"{path}: need at least two numeric (strain, stress) rows")
        a = np.array(rows)
        return cls(a[:, 0], a[:, 1], beyond=beyond)

    def with_default_beyond(self, low: str, high: str) -> "TabularLaw":
        lo = self.beyond[0] or low
        hi = self.beyond[1] or high
        return TabularLaw(self.e, self.s, beyond=(lo, hi), yield_strain=self._eps_y)

    def _resolved(self):
        return (self.beyond[0] or "hold", self.beyond[1] or "hold")

    def stress(self, eps):
        eps = np.asarray(eps, dtype=float)
        out = np.interp(eps, self.e, self.s)
        lo, hi = self._resolved()
        if lo == "zero":
            out = np.where(eps < self.e[0], 0.0, out)
        if hi == "zero":
            out = np.where(eps > self.e[-1], 0.0, out)
        return out

    def breakpoints(self):
        return self.e.copy()

    def scaled(self, stress_factor=1.0, strain_factor=1.0):
        eps_y = None if self._eps_y is None else self._eps_y * strain_factor
        return TabularLaw(self.e * strain_factor, self.s * stress_factor, beyond=self.beyond, yield_strain=eps_y)

    def rupture_limits(self):
        lo, hi = self._resolved()
        return (self.e[0] if lo == "zero" else None, self.e[-1] if hi == "zero" else None)


class FunctionLaw(StressStrainLaw):
    """Law defined by a vectorised Python function ``func(eps_array) -> stress_array``.

    ``breakpoints`` lists strains where the function is non-smooth or
    discontinuous.  For smooth non-polynomial curves, ``n_sub > 1`` adds an even
    grid of ``n_sub`` sub-intervals over ``eps_range`` to refine the integration.
    The function is evaluated right-continuously at breakpoints.
    """

    def __init__(self, func: Callable, breakpoints: Iterable[float] = (), n_sub: int = 1,
                 eps_range=(-0.01, 0.01), rupture=(None, None)):
        self.func = func
        self._bp = np.asarray(list(breakpoints), dtype=float)
        self.n_sub = int(n_sub)
        self.eps_range = tuple(eps_range)
        self._rupture = tuple(rupture)
        self._sf, self._ef = 1.0, 1.0

    def stress(self, eps):
        eps = np.asarray(eps, dtype=float)
        return self._sf * np.asarray(self.func(eps / self._ef), dtype=float)

    def breakpoints(self):
        pts = self._bp * self._ef
        if self.n_sub > 1:
            pts = np.concatenate([pts, np.linspace(self.eps_range[0] * self._ef,
                                                   self.eps_range[1] * self._ef, self.n_sub + 1)])
        return np.unique(pts)

    def scaled(self, stress_factor=1.0, strain_factor=1.0):
        new = FunctionLaw(self.func, self._bp, self.n_sub, self.eps_range, self._rupture)
        new._sf = self._sf * stress_factor
        new._ef = self._ef * strain_factor
        return new

    def rupture_limits(self):
        lo, hi = self._rupture
        return (None if lo is None else lo * self._ef, None if hi is None else hi * self._ef)


def elastic_plastic(E: float, fy: float, eps_big: float = 1.0) -> TabularLaw:
    """Symmetric elastic-perfectly-plastic law (holds the yield stress forever)."""
    ey = fy / E
    return TabularLaw([-eps_big, -ey, ey, eps_big], [-fy, -fy, fy, fy], beyond="hold", yield_strain=ey)


def whitney_step(fc: float, beta1: float, eps_cu: float, factor: float = 0.85) -> FunctionLaw:
    """ACI equivalent rectangular block expressed as a step in strain.

    With the top-fibre pivot the strain at depth ``a = beta1*c`` equals
    ``(1 - beta1) * eps_cu``, so the block is exactly a step law in strain.
    """
    thr = (1.0 - beta1) * eps_cu
    level = factor * fc
    return FunctionLaw(lambda e: np.where(e >= thr, level, 0.0), breakpoints=[thr])
