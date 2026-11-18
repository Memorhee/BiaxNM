"""Polygon helpers: circle(), rounded() (corner fillets) and h_section()."""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from ..section import dedupe_ring


def circle(cx: float, cy: float, r: float, n: int = 64, equal_area: bool = False):
    """Regular n-gon approximating a circle (counter-clockwise).

    Vertices lie on the circle of radius ``r``; with ``equal_area=True`` the
    radius is inflated so the polygon area equals pi*r^2.
    """
    if n < 3 or r <= 0:
        raise ValueError("circle needs n >= 3 and r > 0")
    if equal_area:
        a = 2 * math.pi / n
        r = r * math.sqrt(a / math.sin(a))
    t = 2 * math.pi * np.arange(n) / n
    return [(cx + r * math.cos(a), cy + r * math.sin(a)) for a in t]


def rounded(points: Sequence[Sequence[float]], r, n_arc: int = 12):
    """Round the corners of a polygon given by its *unrounded* vertices.

    ``r`` is a scalar or a per-vertex list (``0`` leaves a vertex sharp).
    ``n_arc`` is the number of arc segments per 90 degrees of turn.  Convex and
    concave corners are both handled.  Raises ValueError when the tangent
    lengths of two neighbouring fillets exceed the edge between them.
    """
    pts = np.asarray(points, dtype=float)
    ext = max(float(np.ptp(pts[:, 0])), float(np.ptp(pts[:, 1])), 1e-300) if pts.ndim == 2 and len(pts) else 1.0
    if len(pts) > 1 and np.max(np.abs(pts[0] - pts[-1])) <= 1e-10 * ext:      # closed ring: drop the repeat
        pts = pts[:-1]
    n = len(pts)
    if n < 3:
        raise ValueError("a polygon needs at least 3 vertices")
    radii = np.full(n, float(r)) if np.isscalar(r) else np.asarray(r, dtype=float)
    if radii.shape != (n,):
        raise ValueError("r must be a scalar or one value per vertex")
    if np.any(radii < 0):
        raise ValueError("fillet radius must be >= 0")

    tangent = np.zeros(n)
    geo = {}
    for i in range(n):
        p0, p1, p2 = pts[i - 1], pts[i], pts[(i + 1) % n]
        v1, v2 = p0 - p1, p2 - p1
        l1, l2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if l1 < 1e-12 or l2 < 1e-12:
            raise ValueError(f"duplicate vertices around index {i}")
        u1, u2 = v1 / l1, v2 / l2
        phi = math.acos(max(-1.0, min(1.0, float(u1 @ u2))))  # interior angle between edges
        if radii[i] == 0 or abs(math.pi - phi) < 1e-9:
            continue
        if phi < 1e-9:
            raise ValueError(f"degenerate spike at vertex {i}")
        t = radii[i] / math.tan(phi / 2)
        tangent[i] = t
        bis = u1 + u2
        bis /= np.linalg.norm(bis)
        center = p1 + bis * radii[i] / math.sin(phi / 2)
        geo[i] = (p1 + u1 * t, p1 + u2 * t, center, radii[i])

    for i in range(n):
        j = (i + 1) % n
        edge = np.linalg.norm(pts[j] - pts[i])
        if tangent[i] + tangent[j] > edge + 1e-9:
            raise ValueError(
                f"fillet radii at vertices {i} and {j} need {tangent[i] + tangent[j]:.3f} "
                f"but the edge between them is only {edge:.3f} long")

    out = []
    for i in range(n):
        if i not in geo:
            out.append(tuple(pts[i]))
            continue
        t1, t2, c, rad = geo[i]
        a1 = math.atan2(*(t1 - c)[::-1])
        a2 = math.atan2(*(t2 - c)[::-1])
        delta = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi
        m = max(1, math.ceil(abs(delta) / (math.pi / 2) * n_arc))
        for k in range(m + 1):
            a = a1 + delta * k / m
            out.append((c[0] + rad * math.cos(a), c[1] + rad * math.sin(a)))
    # two fillets that exactly use up the edge between them end and start at the same point
    return [tuple(p) for p in dedupe_ring(np.array(out))]


def h_section(d: float, bf: float, tw: float, tf: float, cx: float = 0.0, cy: float = 0.0,
              rotation_deg: float = 0.0, r: float = 0.0, n_arc: int = 8):
    """Rolled H-section polygon (counter-clockwise).

    Web runs along y (strong axis = x) before rotation.  ``r`` adds fillets at the
    four web/flange junctions.  The section is rotated counter-clockwise by
    ``rotation_deg`` about its centre ``(cx, cy)``.
    """
    if not (0 < tw < bf and 0 < 2 * tf < d):
        raise ValueError("invalid H-section proportions")
    hd, hb, ht = d / 2, bf / 2, tw / 2
    pts = [(-hb, -hd), (hb, -hd), (hb, -hd + tf), (ht, -hd + tf), (ht, hd - tf), (hb, hd - tf),
           (hb, hd), (-hb, hd), (-hb, hd - tf), (-ht, hd - tf), (-ht, -hd + tf), (-hb, -hd + tf)]
    if r > 0:
        radii = [0, 0, 0, r, r, 0, 0, 0, 0, r, r, 0]
        pts = rounded(pts, radii, n_arc=n_arc)
    a = math.radians(rotation_deg)
    ca, sa = math.cos(a), math.sin(a)
    return [(cx + x * ca - y * sa, cy + x * sa + y * ca) for x, y in pts]


# ----------------------------------------------------------------------------
# typical sections used by the configuration files
# ----------------------------------------------------------------------------
def _rotate(pts, rotation_deg: float, cx: float, cy: float):
    if not rotation_deg:
        return [tuple(map(float, p)) for p in pts]
    a = math.radians(rotation_deg)
    ca, sa = math.cos(a), math.sin(a)
    return [(cx + (x - cx) * ca - (y - cy) * sa, cy + (x - cx) * sa + (y - cy) * ca) for x, y in pts]


def rectangle(x: float, y: float, width: float, height: float, r: float = 0.0, n_arc: int = 12):
    """Rectangle with its lower-left corner at (x, y); ``r`` rounds all four corners."""
    if not (width > 0 and height > 0):
        raise ValueError("rectangle needs width > 0 and height > 0")
    pts = [(x, y), (x + width, y), (x + width, y + height), (x, y + height)]
    return rounded(pts, r, n_arc) if r > 0 else pts


def track(cx: float, cy: float, length: float, width: float, rotation_deg: float = 0.0, n_arc: int = 12):
    """Track (stadium) section: two semicircles of diameter ``width`` joined by straight sides.

    ``length`` is the overall length along x before rotation (``length > width``); the section is rotated
    counter-clockwise by ``rotation_deg`` about its centre (cx, cy).
    """
    if not (width > 0 and length > width):
        raise ValueError("track needs width > 0 and length > width")
    hl, hw = length / 2, width / 2
    pts = rounded([(cx - hl, cy - hw), (cx + hl, cy - hw), (cx + hl, cy + hw), (cx - hl, cy + hw)], hw, n_arc)
    return _rotate(pts, rotation_deg, cx, cy)


def t_section(x: float, y: float, bf: float, h: float, tf: float, tw: float, rotation_deg: float = 0.0):
    """T-section, flange on top: bounding box lower-left corner (x, y), flange width ``bf``, total depth ``h``,
    flange thickness ``tf`` and web thickness ``tw``; rotated about the box centre."""
    if not (0 < tw < bf and 0 < tf < h):
        raise ValueError("invalid T-section proportions (need 0 < tw < bf and 0 < tf < h)")
    m = x + bf / 2
    pts = [(m - tw / 2, y), (m + tw / 2, y), (m + tw / 2, y + h - tf), (x + bf, y + h - tf), (x + bf, y + h),
           (x, y + h), (x, y + h - tf), (m - tw / 2, y + h - tf)]
    return _rotate(pts, rotation_deg, x + bf / 2, y + h / 2)


def l_section(x: float, y: float, width: float, height: float, tx: float, ty: float, r: float = 0.0,
              rotation_deg: float = 0.0, n_arc: int = 12):
    """L-section: bounding box lower-left corner (x, y), overall ``width`` x ``height``.

    The vertical leg (on the left) is ``tx`` wide and the horizontal leg (at the bottom) is ``ty`` high.
    ``r`` rounds the re-entrant corner; the section is rotated about the box centre.
    """
    if not (0 < tx < width and 0 < ty < height):
        raise ValueError("invalid L-section proportions (need 0 < tx < width and 0 < ty < height)")
    pts = [(x, y), (x + width, y), (x + width, y + ty), (x + tx, y + ty), (x + tx, y + height), (x, y + height)]
    if r > 0:
        pts = rounded(pts, [0, 0, 0, r, 0, 0], n_arc)
    return _rotate(pts, rotation_deg, x + width / 2, y + height / 2)


def hollow_rect(cx: float, cy: float, width: float, height: float, t: float, rotation_deg: float = 0.0):
    """Hollow rectangular (box) steel section of uniform wall thickness ``t``, centred at (cx, cy).

    Returns ``(outer, hole)``: the outer outline and the void.
    """
    if not (t > 0 and 2 * t < min(width, height)):
        raise ValueError("invalid hollow rectangle (need 0 < t and 2 t < min(width, height))")
    def box(w, h):
        return [(cx - w / 2, cy - h / 2), (cx + w / 2, cy - h / 2), (cx + w / 2, cy + h / 2), (cx - w / 2, cy + h / 2)]
    return (_rotate(box(width, height), rotation_deg, cx, cy),
            _rotate(box(width - 2 * t, height - 2 * t), rotation_deg, cx, cy))