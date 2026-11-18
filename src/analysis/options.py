"""Validation of the arguments of :func:`biaxnm.diagram.pm_curve`."""
from __future__ import annotations

import math

import numpy as np

from ..core.integration import Integrator


def _resolve_reference(integ: Integrator, reference):
    if isinstance(reference, str):
        if reference == "plastic":
            return integ.plastic_centroid()
        if reference == "geometric":
            return integ.section.centroid
        raise ValueError("reference must be 'plastic', 'geometric' or (x, y)")
    try:
        x, y = reference
        x, y = float(x), float(y)
    except (TypeError, ValueError):
        raise ValueError("reference must be 'plastic', 'geometric' or (x, y)") from None
    if not (math.isfinite(x) and math.isfinite(y)):
        raise ValueError("reference point must be finite")
    return (x, y)


def _check_args(lambda_u_deg, n_c, c_min_ratio, c_max_ratio, tol_deg, scan_deg):
    def finite(v):
        return isinstance(v, (int, float, np.integer, np.floating)) and math.isfinite(v)

    if not finite(lambda_u_deg):
        raise ValueError(f"lambda_u_deg must be a finite number, got {lambda_u_deg!r}")
    if not (isinstance(n_c, (int, np.integer)) and n_c >= 2):
        raise ValueError(f"n_c must be an integer >= 2, got {n_c!r}")
    if not (finite(c_min_ratio) and finite(c_max_ratio) and 0 < c_min_ratio < c_max_ratio):
        raise ValueError(f"need 0 < c_min_ratio < c_max_ratio, got {c_min_ratio!r}, {c_max_ratio!r}")
    if not (finite(tol_deg) and tol_deg > 0):
        raise ValueError(f"tol_deg must be a positive finite number, got {tol_deg!r}")
    if not (finite(scan_deg) and scan_deg > 0):
        raise ValueError(f"scan_deg must be a positive finite number, got {scan_deg!r}")
