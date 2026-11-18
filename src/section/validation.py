"""Geometric validation of a :class:`PolygonSection` (rings, steel, bars, tendons)."""
from __future__ import annotations

import math

import numpy as np

from .ringgeom import _inside, _on_boundary, _rings_touch, _self_intersects


def _steel_pair_ok(a, b, tol) -> bool:
    """Two steel regions may be disjoint, or one may sit entirely inside a void of the other."""
    if _rings_touch(a.polygon, b.polygon, tol):
        return False
    if any(_rings_touch(a.polygon, h, tol) for h in b.holes) or any(_rings_touch(b.polygon, h, tol) for h in a.holes):
        return False
    if _inside(a.polygon[:1], b.polygon)[0]:
        return any(_inside(a.polygon, h).all() for h in b.holes)
    if _inside(b.polygon[:1], a.polygon)[0]:
        return any(_inside(b.polygon, h).all() for h in a.holes)
    return True


def validate_section(sec) -> None:
    """Raise ``ValueError`` when the (already cleaned and oriented) section is inconsistent."""
    tol = sec._tol
    rings = [("outer", sec.outer)] + [(f"hole {i}", h) for i, h in enumerate(sec.holes)]
    for name, r in rings:
        if _self_intersects(r, tol):
            raise ValueError(f"{name} is self-intersecting")
    for i, h in enumerate(sec.holes):
        if not _inside(h, sec.outer).all() or _rings_touch(h, sec.outer, tol):
            raise ValueError(f"hole {i} is not strictly inside the outer outline")
        for j in range(i):
            g = sec.holes[j]
            if _rings_touch(h, g, tol) or _inside(h, g).any() or _inside(g, h).any():
                raise ValueError(f"holes {j} and {i} overlap")

    for i, s in enumerate(sec.steel):
        if _self_intersects(s.polygon, tol):
            raise ValueError(f"steel {i} is self-intersecting")
        if not _inside(s.polygon, sec.outer).all() or _rings_touch(s.polygon, sec.outer, tol):
            raise ValueError(f"steel {i} is not strictly inside the concrete outline")
        for k, h in enumerate(sec.holes):
            if _rings_touch(s.polygon, h, tol) or _inside(s.polygon, h).any():
                raise ValueError(f"steel {i} touches or lies in concrete hole {k}")
            if any(_rings_touch(sh, h, tol) for sh in s.holes):
                raise ValueError(f"steel {i}: the boundary of a steel void meets concrete hole {k}")
            if _inside(h[:1], s.polygon)[0] and not any(_inside(h, sh).all() for sh in s.holes):
                raise ValueError(f"concrete hole {k} lies inside steel {i} (outside its voids)")
        for j, sh in enumerate(s.holes):
            if _self_intersects(sh, tol) or not _inside(sh, s.polygon).all() \
                    or _rings_touch(sh, s.polygon, tol):
                raise ValueError(f"steel {i} hole {j} is not strictly inside the steel outline")
            for m in range(j):
                g = s.holes[m]
                if _rings_touch(sh, g, tol) or _inside(sh, g).any() or _inside(g, sh).any():
                    raise ValueError(f"steel {i} holes {m} and {j} overlap")
        for j in range(i):
            if not _steel_pair_ok(s, sec.steel[j], tol):
                raise ValueError(f"steel regions {j} and {i} overlap")

    btol = 1e-9 * max(sec.diag, 1.0)          # a point this close to an outline counts as on it
    for kind, items in (("bar", sec.bars), ("tendon", sec.tendons)):
        for i, b in enumerate(items):
            if not (math.isfinite(b.area) and b.area > 0):
                raise ValueError(f"{kind} {i}: area must be positive and finite")
            if not (math.isfinite(b.x) and math.isfinite(b.y)):
                raise ValueError(f"{kind} {i}: coordinates must be finite")
            pt = np.array([[b.x, b.y]])
            # the outline itself belongs to the section: a point on it is accepted on every edge
            if not (_inside(pt, sec.outer)[0] or _on_boundary(pt, sec.outer, btol)[0]):
                raise ValueError(f"{kind} {i} at ({b.x}, {b.y}) is outside the concrete outline")
            for k, h in enumerate(sec.holes):
                if _inside(pt, h)[0] and not _on_boundary(pt, h, btol)[0]:
                    raise ValueError(f"{kind} {i} at ({b.x}, {b.y}) lies inside hole {k}")
            for k, s in enumerate(sec.steel):
                if _inside(pt, s.polygon)[0] and not _on_boundary(pt, s.polygon, btol)[0] \
                        and not any(_inside(pt, sh)[0] or _on_boundary(pt, sh, btol)[0] for sh in s.holes):
                    raise ValueError(f"{kind} {i} at ({b.x}, {b.y}) overlaps steel region {k}")
    for i, t in enumerate(sec.tendons):
        if not (math.isfinite(t.fse) and t.fse >= 0):
            raise ValueError(f"tendon {i}: effective prestress fse must be finite and >= 0, got {t.fse}")
        if t.eps_pe is not None and not math.isfinite(t.eps_pe):
            raise ValueError(f"tendon {i}: eps_pe must be finite")
        if t.duct_area is not None and not (math.isfinite(t.duct_area) and t.duct_area >= 0):
            raise ValueError(f"tendon {i}: duct_area must be finite and >= 0, got {t.duct_area}")
        for j, b in enumerate(sec.bars):
            if math.hypot(t.x - b.x, t.y - b.y) <= btol:
                raise ValueError(f"tendon {i} coincides with bar {j}")

