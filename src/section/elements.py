"""Input records: reinforcement, tendons and encased structural steel."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

from ..core.laws import StressStrainLaw


@dataclass
class Bar:
    x: float
    y: float
    area: float
    law: Optional[StressStrainLaw] = None


@dataclass
class Tendon:
    x: float
    y: float
    area: float
    fse: float = 0.0                 # effective prestress (tension), MPa
    law: Optional[StressStrainLaw] = None
    eps_pe: Optional[float] = None   # overrides fse/Ep when given (tension strain)
    duct_area: Optional[float] = None  # duct area Apd [mm^2] that displaces concrete; None -> same as area (Apd = Apt)


@dataclass
class SteelRegion:
    polygon: Sequence[Sequence[float]]
    holes: Sequence[Sequence[Sequence[float]]] = ()
    law: Optional[StressStrainLaw] = None
