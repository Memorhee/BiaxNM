"""Moment-direction angle helper."""
from __future__ import annotations

import math
import warnings


def lambda_u_deg(Mux: float, Muy: float) -> float:
    """Angle of the moment vector, atan2(Muy, Mux), in degrees."""
    if Mux == 0 and Muy == 0:
        warnings.warn("Mux = Muy = 0: lambda_u is taken as 0 degrees")
        return 0.0
    return math.degrees(math.atan2(Muy, Mux))
