"""The section object: cleaned geometry plus area, centroid and inertia of the concrete."""
from __future__ import annotations

import numpy as np

from .elements import Bar, SteelRegion, Tendon
from .ringgeom import _clean_ring, _orient, ring_properties
from .validation import validate_section


class PolygonSection:
    """Reinforced / composite section built from polygons.

    outer     concrete outline (list of vertices)
    holes     list of hole outlines
    bars      list of :class:`Bar`
    steel     list of :class:`SteelRegion` (encased structural steel)
    tendons   list of :class:`Tendon`
    """

    def __init__(self, outer, holes=(), bars=(), steel=(), tendons=()):
        self.outer = _orient(_clean_ring(outer, "outer"), True)
        self.holes = [_orient(_clean_ring(h, f"hole {i}"), False) for i, h in enumerate(holes)]
        self.bars = list(bars)
        self.steel = []
        for i, s in enumerate(steel):
            poly = _orient(_clean_ring(s.polygon, f"steel {i}"), True)
            sholes = [_orient(_clean_ring(h, f"steel {i} hole {j}"), False) for j, h in enumerate(s.holes)]
            self.steel.append(SteelRegion(poly, sholes, s.law))
        self.tendons = list(tendons)

        allp = np.vstack([self.outer])
        self.bbox = (allp[:, 0].min(), allp[:, 1].min(), allp[:, 0].max(), allp[:, 1].max())
        self.diag = float(np.hypot(self.bbox[2] - self.bbox[0], self.bbox[3] - self.bbox[1]))
        self._tol = 1e-9 * max(self.diag, 1.0) ** 2
        validate_section(self)
        self._properties()

    def _properties(self):
        # work about the middle of the bounding box so that sections far from the origin keep their precision
        ref = np.array([0.5 * (self.bbox[0] + self.bbox[2]), 0.5 * (self.bbox[1] + self.bbox[3])])
        A = Sx = Sy = Ixx = Iyy = Ixy = 0.0
        for ring in [self.outer] + self.holes:
            a, sx, sy, ixx, iyy, ixy = ring_properties(ring - ref)
            A += a; Sx += sx; Sy += sy; Ixx += ixx; Iyy += iyy; Ixy += ixy
        self.area = A
        cx, cy = Sy / A, Sx / A                       # centroid relative to ref
        self.centroid = (float(ref[0] + cx), float(ref[1] + cy))
        self.Ix = Ixx - A * cy * cy
        self.Iy = Iyy - A * cx * cx
        self.Ixy = Ixy - A * cx * cy

    @property
    def has_tendons(self):
        return bool(self.tendons)

    @property
    def has_steel(self):
        return bool(self.steel)
