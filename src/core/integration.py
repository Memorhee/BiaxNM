"""Exact strip integration of stress resultants over polygonal regions.

For a compression direction n = (cos t, sin t) every region is rotated into
(d, u) coordinates (d along n).  Between consecutive breakpoints in d (polygon
vertices and strains where a stress-strain law kinks) the section width b(d)
is linear and the first moment of the width is quadratic, so 3-point
Gauss-Legendre quadrature is exact for linear and parabolic stress laws.
Holes are rings with clockwise orientation and subtract automatically.

Strain plane:  eps(d) = a + k*d   (compression positive)
Resultants:    P = sum F,  Qx = sum F*y,  Qy = sum F*x  (about the origin)
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..codes.base import DesignCode, Geom
from .materials import ResolvedMaterials, edge_arrays, group_by_law  # noqa: F401
from ..section import PolygonSection, clip_halfplane
from .strainstate import StrainReport

_XG, _WG = np.polynomial.legendre.leggauss(3)

def _integrate_region(edges, cs, sigma, bps, a, k):
    """Return (N, Fd, Fu) = integrals of sigma, sigma*d, sigma*(first moment in u)."""
    X1, Y1, X2, Y2 = edges
    c, s = cs
    d1 = X1 * c + Y1 * s
    d2 = X2 * c + Y2 * s
    u1 = -X1 * s + Y1 * c
    u2 = -X2 * s + Y2 * c
    dd = d2 - d1
    ok = np.abs(dd) > 1e-10
    d1, d2, u1, u2, dd = d1[ok], d2[ok], u1[ok], u2[ok], dd[ok]
    sgn = -np.sign(dd)
    lo, hi = np.minimum(d1, d2), np.maximum(d1, d2)

    pts = np.concatenate([d1, d2])
    if abs(k) > 1e-15 and len(bps):
        bd = (bps - a) / k
        bd = bd[(bd > pts.min()) & (bd < pts.max())]
        pts = np.concatenate([pts, bd])
    pts = np.sort(pts)
    pts = pts[np.concatenate(([True], np.diff(pts) > 1e-9))]
    left, right = pts[:-1], pts[1:]
    mid, h = 0.5 * (left + right), right - left
    t = mid[:, None] + 0.5 * h[:, None] * _XG[None, :]
    w = 0.5 * h[:, None] * _WG[None, :]

    span = ((lo[:, None] <= mid[None, :]) & (mid[None, :] <= hi[:, None]))[:, :, None]
    slope = (u2 - u1) / dd
    u = u1[:, None, None] + slope[:, None, None] * (t[None, :, :] - d1[:, None, None])
    sg = sgn[:, None, None]
    b = np.sum(np.where(span, sg * u, 0.0), axis=0)
    q = np.sum(np.where(span, 0.5 * sg * u * u, 0.0), axis=0)
    sig = sigma(a + k * t)
    return float(np.sum(sig * b * w)), float(np.sum(sig * b * t * w)), float(np.sum(sig * q * w))


# kept under their former names for callers inside the package and the tests
_edge_arrays = edge_arrays
_group = group_by_law


class Integrator:
    """Evaluates stress resultants of a section for a strain plane.

    Material binding lives in :class:`ResolvedMaterials` and strain diagnostics in
    :class:`StrainReport`; their attributes and methods stay reachable from here.
    """

    def __init__(self, section: PolygonSection, code: DesignCode):
        self.section, self.code = section, code
        self.materials = ResolvedMaterials(section, code)
        self.report = StrainReport(self.materials)
        self.n_evals = 0

    def __getattr__(self, name):
        if name in ("materials", "report"):
            raise AttributeError(name)
        return getattr(self.materials, name)

    # diagnostics (see strainstate.py)
    def tension_extremes(self, theta, a, k):
        return self.report.tension_extremes(theta, a, k)

    def governing_yield(self, theta, a, k):
        return self.report.governing_yield(theta, a, k)

    def tension_limit_strain(self) -> float:
        return self.report.tension_limit_strain()

    def rupture_flags(self, theta, a, k):
        return self.report.rupture_flags(theta, a, k)

    def _strain_sets(self, theta, a, k):
        return self.report._strain_sets(theta, a, k)

    # ---------------------------------------------------------------- geometry
    def geometry(self, theta: float) -> Geom:
        c, s = math.cos(theta), math.sin(theta)
        dc = self.conc_pts[:, 0] * c + self.conc_pts[:, 1] * s
        d_max, d_min = float(dc.max()), float(dc.min())
        if self.rein_pts is None:
            d_rein = d_min
        else:
            d_rein = float((self.rein_pts[:, 0] * c + self.rein_pts[:, 1] * s).min())
        return Geom(d_max, d_min, d_rein)

    def plane(self, theta: float, c: float) -> Tuple[float, float, Geom]:
        g = self.geometry(theta)
        k = self.code.pivot.curvature(c, g)
        return -k * (g.d_max - c), k, g

    # --------------------------------------------------------------- resultants
    def resultants(self, theta: float, a: float, k: float) -> Tuple[float, float, float]:
        self.n_evals += 1
        c, s = math.cos(theta), math.sin(theta)
        conc = self.conc
        N, Fd, Fu = _integrate_region(self.conc_edges, (c, s), conc.stress, self.conc_bps, a, k)

        for edges, law, _ in self.steel:
            def sig(e, law=law):
                return law.stress(e) - conc.stress(e)
            bps = np.union1d(law.breakpoints(), self.conc_bps)
            n2, fd2, fu2 = _integrate_region(edges, (c, s), sig, bps, a, k)
            N += n2; Fd += fd2; Fu += fu2

        Qx = s * Fd + c * Fu
        Qy = c * Fd - s * Fu
        P = N

        if len(self.bar_x):
            d = self.bar_x * c + self.bar_y * s
            eps = a + k * d
            sig_c = conc.stress(eps)
            f = np.empty_like(eps)
            for law, idx in self.bar_groups:
                f[idx] = law.stress(eps[idx])
            F = self.bar_A * (f - sig_c)
            P += F.sum(); Qx += (F * self.bar_y).sum(); Qy += (F * self.bar_x).sum()

        if len(self.t_x):
            d = self.t_x * c + self.t_y * s
            eps = a + k * d
            sig_c = conc.stress(eps)
            f = np.empty_like(eps)
            for law, idx in self.t_groups:
                f[idx] = law.stress(eps[idx] - self.t_pe[idx])
            F = self.t_A * f - self.t_Ad * sig_c
            P += F.sum(); Qx += (F * self.t_y).sum(); Qy += (F * self.t_x).sum()
        return P, Qx, Qy

    def components(self, theta: float, a: float, k: float):
        """Stress resultants split by material: list of ``(name, F, Qx, Qy)`` about the origin.

        ``F`` is the force, ``Qx = sum F*y`` and ``Qy = sum F*x`` (so the force acts at ``x = Qy/F``,
        ``y = Qx/F``).  Steel, rebar and tendon forces are net of the concrete they displace, exactly as in
        :meth:`resultants`; the components add up to its result.

        So that every force has a meaningful point of action, rebar and tendons are listed as a compression
        group "(C)" and a tension group "(T)" (by the sign of each bar's force), and a steel region that
        crosses the neutral axis is cut along it into a compressed and a tensioned part.
        """
        c, s = math.cos(theta), math.sin(theta)
        n = np.array([c, s])
        conc = self.conc
        out = []
        N, Fd, Fu = _integrate_region(self.conc_edges, (c, s), conc.stress, self.conc_bps, a, k)
        out.append(("concrete", N, s * Fd + c * Fu, c * Fd - s * Fu))
        d_na = -a / k if k != 0 else None
        for i, (edges, law, _) in enumerate(self.steel):
            label = f"steel {i + 1}" if len(self.steel) > 1 else "steel"

            def sig(e, law=law):
                return law.stress(e) - conc.stress(e)
            bps = np.union1d(law.breakpoints(), self.conc_bps)
            region = self.section.steel[i]
            parts = []
            if d_na is not None:
                rings = [region.polygon] + list(region.holes)
                for tag, nn, d0 in (("C", n, d_na), ("T", -n, -d_na)):
                    cl = [clip_halfplane(r, nn, d0) for r in rings]
                    if cl[0] is not None:
                        parts.append((tag, edge_arrays([r for r in cl if r is not None])))
            if len(parts) != 2:
                parts = [("", edges)]
            for tag, e_arr in parts:
                n2, fd2, fu2 = _integrate_region(e_arr, (c, s), sig, bps, a, k)
                out.append((f"{label} ({tag})" if tag else label, n2, s * fd2 + c * fu2, c * fd2 - s * fu2))
        if len(self.bar_x):
            eps = a + k * (self.bar_x * c + self.bar_y * s)
            f = np.empty_like(eps)
            for law, idx in self.bar_groups:
                f[idx] = law.stress(eps[idx])
            out += self._by_sign("rebar", self.bar_A * (f - conc.stress(eps)), self.bar_x, self.bar_y)
        if len(self.t_x):
            eps = a + k * (self.t_x * c + self.t_y * s)
            f = np.empty_like(eps)
            for law, idx in self.t_groups:
                f[idx] = law.stress(eps[idx] - self.t_pe[idx])
            out += self._by_sign("tendon", self.t_A * f - self.t_Ad * conc.stress(eps), self.t_x, self.t_y)
        return out

    @staticmethod
    def _by_sign(name, F, x, y):
        res = []
        for tag, mask in (("C", F > 0), ("T", F < 0)):
            if mask.any():
                res.append((f"{name} ({tag})", float(F[mask].sum()), float((F * y)[mask].sum()),
                            float((F * x)[mask].sum())))
        return res

    def uniform(self) -> Tuple[float, float, float]:
        """Pure compression state: (P0, Qx, Qy) at the pivot's uniform strain."""
        return self.resultants(0.0, self.code.pivot.uniform_strain(), 0.0)

    def plastic_centroid(self) -> Tuple[float, float]:
        P0, Qx, Qy = self.uniform()
        if P0 <= 0:
            return self.section.centroid
        return (Qy / P0, Qx / P0)
