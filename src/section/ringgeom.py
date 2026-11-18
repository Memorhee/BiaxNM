"""Pure polygon geometry on (n, 2) vertex arrays: cleaning, orientation, intersection, containment, properties."""
from __future__ import annotations

import numpy as np


def dedupe_ring(p: np.ndarray) -> np.ndarray:
    """Drop consecutive duplicate vertices and a repeated closing vertex.

    Two vertices are duplicates when they differ by less than 1e-10 times the
    extent of the ring (not of the coordinates, so a ring far from the origin
    keeps its short edges).
    """
    p = np.asarray(p, dtype=float)
    if len(p) < 2:
        return p
    tol = 1e-10 * max(float(np.ptp(p[:, 0])), float(np.ptp(p[:, 1])), 1e-300)
    keep = [0]
    for i in range(1, len(p)):
        if np.max(np.abs(p[i] - p[keep[-1]])) > tol:
            keep.append(i)
    while len(keep) > 1 and np.max(np.abs(p[keep[-1]] - p[keep[0]])) <= tol:
        keep.pop()
    return p[keep]


def _clean_ring(ring, name="ring") -> np.ndarray:
    p = np.asarray(ring, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"{name}: expected a list of (x, y) vertices")
    if not np.all(np.isfinite(p)):
        raise ValueError(f"{name}: NaN/inf coordinates")
    p = dedupe_ring(p)
    if len(p) < 3:
        raise ValueError(f"{name}: need at least 3 vertices")
    return p


def signed_area(p: np.ndarray) -> float:
    x, y = p[:, 0], p[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _orient(p: np.ndarray, ccw: bool) -> np.ndarray:
    a = signed_area(p)
    if abs(a) < 1e-12:
        raise ValueError("polygon has zero area")
    return p if (a > 0) == ccw else p[::-1].copy()


_REL_TOL = 1e-12      # orientation tests are relative to the lengths involved, not to the section size


def _orient_sign(p, q, r):
    """Sign of (q-p) x (r-p); 0 when r lies on the line pq to within a relative 1e-12."""
    dq0, dq1 = q[..., 0] - p[..., 0], q[..., 1] - p[..., 1]
    dr0, dr1 = r[..., 0] - p[..., 0], r[..., 1] - p[..., 1]
    cr = dq0 * dr1 - dq1 * dr0
    lq = np.hypot(dq0, dq1)
    scale = lq * (np.hypot(dr0, dr1) + lq)
    return np.where(np.abs(cr) <= _REL_TOL * scale, 0, np.sign(cr))


def _seg_hits(A1, A2, B1, B2, tol=None):
    """Boolean (nA, nB) matrix: segment A_i touches or crosses segment B_j (``tol`` is unused)."""
    a1, a2 = A1[:, None, :], A2[:, None, :]
    b1, b2 = B1[None, :, :], B2[None, :, :]
    s1 = _orient_sign(b1, b2, a1)
    s2 = _orient_sign(b1, b2, a2)
    s3 = _orient_sign(a1, a2, b1)
    s4 = _orient_sign(a1, a2, b2)
    general = (s1 * s2 <= 0) & (s3 * s4 <= 0)
    colin = (s1 == 0) & (s2 == 0) & (s3 == 0) & (s4 == 0)
    ov = ((np.minimum(a1[..., 0], a2[..., 0]) <= np.maximum(b1[..., 0], b2[..., 0]))
          & (np.minimum(b1[..., 0], b2[..., 0]) <= np.maximum(a1[..., 0], a2[..., 0]))
          & (np.minimum(a1[..., 1], a2[..., 1]) <= np.maximum(b1[..., 1], b2[..., 1]))
          & (np.minimum(b1[..., 1], b2[..., 1]) <= np.maximum(a1[..., 1], a2[..., 1])))
    return np.where(colin, ov, general)


def _rings_touch(P: np.ndarray, Q: np.ndarray, tol) -> bool:
    return bool(_seg_hits(P, np.roll(P, -1, axis=0), Q, np.roll(Q, -1, axis=0), tol).any())


def _self_intersects(P: np.ndarray, tol) -> bool:
    n = len(P)
    hits = _seg_hits(P, np.roll(P, -1, axis=0), P, np.roll(P, -1, axis=0), tol)
    i, j = np.indices((n, n))
    adjacent = (i == j) | ((i - j) % n == 1) | ((j - i) % n == 1)
    return bool((hits & ~adjacent).any())


def _inside(points: np.ndarray, ring: np.ndarray) -> np.ndarray:
    """Even-odd point-in-polygon test for an (m, 2) array."""
    pts = np.atleast_2d(points)
    x, y = pts[:, 0][:, None], pts[:, 1][:, None]
    x1, y1 = ring[:, 0][None, :], ring[:, 1][None, :]
    x2, y2 = np.roll(ring[:, 0], -1)[None, :], np.roll(ring[:, 1], -1)[None, :]
    straddle = (y1 > y) != (y2 > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        xint = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
    return (np.sum(straddle & (x < xint), axis=1) % 2) == 1


def _on_boundary(points: np.ndarray, ring: np.ndarray, tol: float) -> np.ndarray:
    """True for each point of an (m, 2) array lying within ``tol`` of the ring's outline."""
    pts = np.atleast_2d(points)[:, None, :]
    a = ring[None, :, :]
    ab = np.roll(ring, -1, axis=0)[None, :, :] - a
    ap = pts - a
    L2 = (ab * ab).sum(-1)
    t = np.clip((ap * ab).sum(-1) / np.where(L2 > 0, L2, 1.0), 0.0, 1.0)
    d = np.linalg.norm(ap - t[..., None] * ab, axis=-1)
    return (d <= tol).any(axis=1)


def clip_halfplane(ring, n, d0):
    """Sutherland-Hodgman clip of a closed ring to the half-plane ``{p : p . n >= d0}``.

    Returns the clipped ring (orientation preserved) or ``None`` when fewer than three vertices remain.
    """
    pts = [np.asarray(p, dtype=float) for p in ring]
    n = np.asarray(n, dtype=float)
    out = []
    for i, p in enumerate(pts):
        q = pts[(i + 1) % len(pts)]
        dp, dq = p @ n - d0, q @ n - d0
        if dp >= 0:
            out.append(p)
        if (dp >= 0) != (dq >= 0):
            out.append(p + (q - p) * (dp / (dp - dq)))
    return np.array(out) if len(out) >= 3 else None


def ring_properties(p: np.ndarray):
    """Signed area, first moments and second moments about the origin."""
    x, y = p[:, 0], p[:, 1]
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y2 - x2 * y
    A = 0.5 * cr.sum()
    Sx = ((y + y2) * cr).sum() / 6.0            # integral of y dA
    Sy = ((x + x2) * cr).sum() / 6.0            # integral of x dA
    Ixx = ((y * y + y * y2 + y2 * y2) * cr).sum() / 12.0
    Iyy = ((x * x + x * x2 + x2 * x2) * cr).sum() / 12.0
    Ixy = ((x * y2 + 2 * x * y + 2 * x2 * y2 + x2 * y) * cr).sum() / 24.0
    return A, Sx, Sy, Ixx, Iyy, Ixy
