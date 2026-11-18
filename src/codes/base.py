"""Design-code interface.

A :class:`DesignCode` supplies everything that differs between design codes:
the material laws (including material factors), the compression pivot that
fixes the strain plane for a given neutral-axis depth, the strength reduction
factor and the axial-load cap.  The integration and solver code never contain
code-specific constants.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from ..core.laws import StressStrainLaw


@dataclass(frozen=True)
class Geom:
    """Extent of the section measured along the compression direction ``d``."""
    d_max: float    # most compressed concrete fibre
    d_min: float    # most tensioned concrete fibre
    d_rein: float   # most tensioned reinforcement / steel / tendon (smallest d)


class Pivot:
    """Defines the strain plane eps(d) = k * (d - d_na) with d_na = d_max - c."""

    def curvature(self, c: float, geom: Geom) -> float:  # pragma: no cover - interface
        raise NotImplementedError

    def uniform_strain(self) -> float:  # pragma: no cover - interface
        """Strain of the pure-compression state (used for P0 and the plastic centroid)."""
        raise NotImplementedError


class TopFiberPivot(Pivot):
    """Strain at the extreme compression fibre is fixed at ``eps_cu`` (ACI 318)."""

    def __init__(self, eps_cu: float = 0.003):
        self.eps_cu = eps_cu

    def curvature(self, c, geom):
        return self.eps_cu / c

    def uniform_strain(self):
        return self.eps_cu


@dataclass(frozen=True)
class TensionState:
    """Tensile state seen by the strength reduction factor.

    ``strain`` maps ``"rebar" | "steel" | "tendon"`` to the largest tensile strain of that material
    (negative when compressed, ``None`` when absent); ``yield_strain`` maps the same kinds to the yield
    strain stated by the law of the most tensioned element (``None`` when the law states none).
    """
    strain: Dict[str, Optional[float]]
    yield_strain: Dict[str, Optional[float]]


class DesignCode:
    """A design code is four parts: material laws (``*_law``), a compression :class:`Pivot`, the strength
    reduction factor (``phi``) and the axial cap (``axial_cap``).  Integration and solving never see
    code-specific constants."""
    pivot: Pivot

    def concrete_law(self) -> StressStrainLaw:  # pragma: no cover - interface
        raise NotImplementedError

    def rebar_law(self) -> StressStrainLaw:
        raise ValueError(f"{type(self).__name__} has no rebar law")

    def steel_law(self) -> StressStrainLaw:
        raise ValueError(f"{type(self).__name__} has no structural steel law")

    def tendon_law(self) -> StressStrainLaw:
        raise ValueError(f"{type(self).__name__} has no tendon law")

    def phi(self, ext: Dict[str, Optional[float]], eps_y: Optional[Dict[str, Optional[float]]] = None) -> float:
        """Strength reduction factor.

        ``ext`` maps ``"rebar" | "steel" | "tendon"`` to the largest *tensile*
        strain of that material (negative when it is in compression) or ``None``
        when the section has none of it.  Tendon strains exclude the effective
        prestrain.  ``eps_y`` maps the same kinds to the yield strain stated by the
        law of the most tensioned element (``None`` where the law does not state one).
        """
        return 1.0

    def phi_for(self, state: TensionState) -> float:
        """Entry point used by the analysis; forwards to ``phi`` so existing overrides keep working."""
        return self.phi(state.strain, state.yield_strain)

    def refine_phi(self, phi, P, states, Ag: float):
        """Curve-level adjustment of the strength reduction factors (default: none).

        ``phi`` and ``P`` are arrays over the curve points (top to bottom; ``P`` is the nominal axial force),
        ``states`` the matching :class:`TensionState` list and ``Ag`` the gross concrete area.  Return the
        adjusted ``phi`` array.
        """
        return phi

    def axial_limit(self, section) -> Optional[float]:
        """Axial strength given by a closed formula of the code for ``section`` (e.g. P_do), used in place of the
        pure-compression point P0 when the cap is evaluated; ``None`` to use P0 of the analysis."""
        return None

    def axial_cap(self, P0: float) -> Tuple[Optional[float], Optional[float]]:
        """Return ``(Pn_max, phiPn_max)`` or ``(None, None)`` for no cap."""
        return (None, None)
