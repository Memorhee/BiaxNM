"""Research code: user-supplied laws, pivot, phi and axial cap."""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from ..core.laws import StressStrainLaw, TabularLaw
from .base import DesignCode, Pivot


def _default(law: Optional[StressStrainLaw], low: str, high: str):
    if isinstance(law, TabularLaw):
        return law.with_default_beyond(low, high)
    return law


class CustomCode(DesignCode):
    """Everything is supplied by the caller.

    concrete / rebar / steel / tendon   stress-strain laws (compression positive)
    pivot                               :class:`~biaxnm.codes.base.Pivot`
    phi                                 constant or ``f(ext) -> float``
    axial_cap                           ``None`` or ``f(P0) -> (Pn_max, phiPn_max)``

    Tabular laws without an explicit ``beyond`` default to: concrete holds on the
    tension side and ruptures (stress 0) beyond the compression end of the table;
    reinforcement, steel and tendons hold their end stress.
    """

    def __init__(self, concrete: StressStrainLaw, pivot: Pivot, rebar: Optional[StressStrainLaw] = None,
                 steel: Optional[StressStrainLaw] = None, tendon: Optional[StressStrainLaw] = None,
                 phi=1.0, axial_cap: Optional[Callable] = None):
        if isinstance(concrete, TabularLaw) and concrete.e[0] > 0:
            # the table starts at a compressive strain: below it the concrete carries no stress
            concrete = TabularLaw(np.concatenate([[-1.0, 0.0], concrete.e]),
                                  np.concatenate([[0.0, 0.0], concrete.s]), beyond=concrete.beyond)
        self._concrete = _default(concrete, "hold", "zero")
        self._rebar = _default(rebar, "hold", "hold")
        self._steel = _default(steel, "hold", "hold")
        self._tendon = _default(tendon, "hold", "hold")
        self.pivot = pivot
        self._phi = phi
        self._cap = axial_cap

    def concrete_law(self):
        return self._concrete

    def rebar_law(self):
        if self._rebar is None:
            raise ValueError("CustomCode: no rebar law was given")
        return self._rebar

    def steel_law(self):
        if self._steel is None:
            raise ValueError("CustomCode: no structural steel law was given")
        return self._steel

    def tendon_law(self):
        if self._tendon is None:
            raise ValueError("CustomCode: no tendon law was given")
        return self._tendon

    def phi(self, ext, eps_y=None):
        return self._phi(ext) if callable(self._phi) else float(self._phi)

    def axial_cap(self, P0):
        return self._cap(P0) if self._cap is not None else (None, None)
