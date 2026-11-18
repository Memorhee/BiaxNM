"""Find the neutral-axis angle theta that satisfies lambda_n(theta) = lambda_u.

For a fixed neutral-axis depth c the strain plane is a function of theta only.
f(theta) = wrap(lambda_n(theta) - lambda_u) is solved with a bracketed Brent
iteration seeded from a starting guess (the elastic estimate for the first c,
the previous solution afterwards = warm start).  A coarse scan over 0..360 deg
is the fallback when no bracket is found.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Tuple

import numpy as np
from scipy.optimize import brentq

from .integration import Integrator


def wrap(a: float) -> float:
    """Wrap an angle to (-pi, pi]."""
    return (a + math.pi) % (2 * math.pi) - math.pi


@dataclass
class Point:
    theta: float
    c: float
    P: float
    Mx: float
    My: float
    a: float
    k: float

    @property
    def M(self) -> float:
        return math.hypot(self.Mx, self.My)


def elastic_theta(section, lambda_u: float) -> float:
    """Neutral-axis direction of the elastic section for a moment angle lambda_u."""
    I = np.array([[section.Ix, section.Ixy], [section.Ixy, section.Iy]])
    kap = np.linalg.solve(I, np.array([math.cos(lambda_u), math.sin(lambda_u)]))
    return math.atan2(kap[0], kap[1])


class AngleSolver:
    def __init__(self, integ: Integrator, lambda_u: float, ref: Tuple[float, float], m_min: float,
                 tol_deg: float = 1e-6, check_deg: float = 0.01, scan_deg: float = 5.0):
        self.integ, self.lam_u, self.ref = integ, lambda_u, ref
        self.m_min = m_min
        self.xtol = math.radians(tol_deg)
        self.check = math.radians(check_deg)
        self.scan = math.radians(scan_deg)
        self.n_global = 0
        self._cache: Dict[Tuple[float, float], Point] = {}

    # ------------------------------------------------------------------ f(theta)
    def point(self, theta: float, c: float) -> Point:
        key = (theta, c)
        p = self._cache.get(key)
        if p is None:
            a, k, _ = self.integ.plane(theta, c)
            P, Qx, Qy = self.integ.resultants(theta, a, k)
            p = Point(theta, c, P, Qx - self.ref[1] * P, Qy - self.ref[0] * P, a, k)
            if len(self._cache) > 4000:
                self._cache.clear()
            self._cache[key] = p
        return p

    def f(self, theta: float, c: float) -> float:
        p = self.point(theta, c)
        if p.M <= self.m_min:
            return math.nan
        return wrap(math.atan2(p.My, p.Mx) - self.lam_u)

    # ------------------------------------------------------------------ brackets
    @staticmethod
    def _genuine(f1: float, f2: float) -> bool:
        return (not math.isnan(f1)) and (not math.isnan(f2)) and f1 * f2 <= 0 and abs(f1 - f2) < math.pi

    def _bracket(self, c: float, th0: float):
        f0 = self.f(th0, c)
        if not math.isnan(f0) and f0 == 0.0:
            return (th0, th0)
        prev = {1: (th0, f0), -1: (th0, f0)}
        w = math.radians(2.0)
        while w <= math.pi + 1e-9:
            for sgn in (1, -1):
                th = th0 + sgn * w
                fv = self.f(th, c)
                pth, pf = prev[sgn]
                if self._genuine(pf, fv):
                    return (min(pth, th), max(pth, th))
                prev[sgn] = (th, fv)
            w *= 2.0
        return None

    def _refine(self, c: float, bracket) -> Optional[float]:
        lo, hi = bracket
        if lo == hi:
            return lo
        try:
            return brentq(lambda t: self.f(t, c), lo, hi, xtol=self.xtol, rtol=1e-12, maxiter=200)
        except (ValueError, RuntimeError):
            return None

    def _accept(self, c: float, th: Optional[float]) -> bool:
        if th is None:
            return False
        fv = self.f(th, c)
        return (not math.isnan(fv)) and abs(fv) <= self.check

    def _global(self, c: float, th0: float) -> Optional[float]:
        self.n_global += 1
        n = max(8, int(round(2 * math.pi / self.scan)))
        grid = th0 + np.linspace(0.0, 2 * math.pi, n + 1)
        vals = [self.f(t, c) for t in grid]
        best, best_dist = None, math.inf
        for i in range(n):
            if self._genuine(vals[i], vals[i + 1]):
                root = self._refine(c, (grid[i], grid[i + 1]))
                if self._accept(c, root):
                    dist = abs(wrap(root - th0))
                    if dist < best_dist:
                        best, best_dist = root, dist
        if best is None:
            best = self._global_close_pairs(c, th0, grid[:n], vals[:n])
        return best

    def _global_close_pairs(self, c: float, th0: float, grid, vals) -> Optional[float]:
        """Two roots closer than the scan step give no sign change on the coarse grid.

        Such a pair sits where |f| has a shallow minimum, so every local minimum of |f| is
        searched again on a grid five times finer than the scan step.
        """
        n = len(grid)
        a = np.array([abs(v) if not math.isnan(v) else math.inf for v in vals])
        best, best_dist = None, math.inf
        for i in range(n):
            if not (a[i] < 0.35 and a[i] <= a[(i - 1) % n] and a[i] <= a[(i + 1) % n]):
                continue
            sub = grid[i] + np.linspace(-self.scan, self.scan, 41)
            sv = [self.f(t, c) for t in sub]
            for j in range(40):
                if self._genuine(sv[j], sv[j + 1]):
                    root = self._refine(c, (sub[j], sub[j + 1]))
                    if self._accept(c, root):
                        dist = abs(wrap(root - th0))
                        if dist < best_dist:
                            best, best_dist = root, dist
        return best

    # ------------------------------------------------------------------ public
    def solve(self, c: float, th0: float, force_global: bool = False,
              allow_global: bool = True) -> Optional[Point]:
        th = None
        if not force_global:
            br = self._bracket(c, th0)
            if br is not None:
                root = self._refine(c, br)
                if self._accept(c, root):
                    th = root
        if th is None and allow_global:
            th = self._global(c, th0)
        if th is None:
            return None
        return self.point(th, c)
