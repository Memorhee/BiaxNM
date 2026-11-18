"""PM interaction curve along a fixed moment direction lambda_u."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from ..analysis import apply_design, bottom_row, rows_from_points, sweep_neutral_axis, top_row
from ..analysis.options import _check_args, _resolve_reference  # noqa: F401  (re-exported)
from ..codes.base import DesignCode
from .integration import Integrator
from ..output.tables import TABLE_HEADER, curve_table
from ..section import PolygonSection
from .solver import AngleSolver


@dataclass
class PMCurve:
    """Points are ordered from pure compression (top) to pure tension (bottom).

    Forces in N, moments in N*mm.  ``Mn`` is the resultant moment, which lies in
    the direction ``lambda_u_deg`` (Mx = Mn cos(lambda), My = Mn sin(lambda)).
    """
    lambda_u_deg: float
    c: np.ndarray
    theta_deg: np.ndarray
    Pn: np.ndarray          # nominal axial force, capped at Pn_max
    Pn_uncapped: np.ndarray
    Mx: np.ndarray
    My: np.ndarray
    Mn: np.ndarray
    eps_t: np.ndarray       # largest tensile strain among reinforcement (nan if none)
    phi: np.ndarray
    phiPn: np.ndarray
    phiMn: np.ndarray
    P0: float
    Pn_max: Optional[float]
    phiPn_max: Optional[float]
    ref_point: Tuple[float, float]
    n_evals: int = 0
    n_global_scans: int = 0
    skipped_c: List[float] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    section: Optional[PolygonSection] = None
    code: Optional[DesignCode] = None

    def table(self) -> np.ndarray:
        """Rows in the units of ``TABLE_HEADER`` (kept here for compatibility; see ``output.tables``)."""
        return curve_table(self)

    TABLE_HEADER = TABLE_HEADER


def pm_curve(section: PolygonSection, code: DesignCode, lambda_u_deg: float, n_c: int = 80,
             c_min_ratio: float = 0.02, c_max_ratio: float = 10.0, reference="geometric",
             tol_deg: float = 1e-6, scan_deg: float = 5.0, full_scan: bool = False) -> PMCurve:
    """Compute the nominal and design PM curve for the moment direction ``lambda_u_deg``.

    c is swept geometrically between ``c_min_ratio`` and ``c_max_ratio`` times the
    section's bounding-box diagonal, starting at 0.5 D and moving down and up.
    ``full_scan=True`` skips the warm start and always scans 0..360 deg (slow;
    for verification).
    """
    _check_args(lambda_u_deg, n_c, c_min_ratio, c_max_ratio, tol_deg, scan_deg)
    lam_u = math.radians(lambda_u_deg)
    integ = Integrator(section, code)
    P0, _, _ = integ.uniform()
    if not math.isfinite(P0):
        raise ValueError("the pure-compression force P0 is not finite: check the material laws and parameters")
    ref = _resolve_reference(integ, reference)
    D = section.diag

    scale = max(abs(P0), 1.0)
    solver = AngleSolver(integ, lam_u, ref, m_min=1e-9 * scale * D, tol_deg=tol_deg, scan_deg=scan_deg)
    cs = np.geomspace(c_max_ratio * D, c_min_ratio * D, n_c)       # descending

    points, skipped = sweep_neutral_axis(solver, section, lam_u, cs, full_scan)
    mid_rows, warnings = rows_from_points(integ, points)
    top = top_row(integ, ref)
    rows = [top] + mid_rows + [bottom_row(integ, ref)]

    if skipped:
        warnings.append(f"{len(skipped)} neutral-axis depths were left out: no theta gives lambda_n = lambda_u "
                        f"there (typically the P0 plateau, where the moment direction is nearly fixed)")
    if abs(P0) < 1e-12:
        warnings.append("P0 is zero")

    p_limit = code.axial_limit(section)            # closed-formula axial strength of the code (e.g. P_do), if any
    fields = apply_design(rows, lam_u, code, top.P if p_limit is None else p_limit, section.area)
    return PMCurve(lambda_u_deg=lambda_u_deg, P0=top.P, ref_point=ref, n_evals=integ.n_evals,
                   n_global_scans=solver.n_global, skipped_c=skipped, warnings=warnings,
                   section=section, code=code, **fields)
