"""biaxnm: biaxial P-M interaction diagrams for arbitrary polygonal sections.

The implementation lives in :mod:`biaxnm.core`; this module only re-exports the public names.
"""
from .core.angles import lambda_u_deg
from .core.check import CheckResult, check_load, check_section
from .codes import ACI318, EC2_1_1, EC2_2, KDS14_22, KDS14_22App, KDS24_25, CustomCode, DesignCode, Pivot, TopFiberPivot
from .core.config import Model, build_model, load_model, run_model
from .core.diagram import PMCurve, pm_curve
from .core.laws import FunctionLaw, StressStrainLaw, TabularLaw, elastic_plastic
from .core.rebar import circle_area, rebar_area
from .output import plot_curve, save_csv, save_outputs
from .section import Bar, PolygonSection, SteelRegion, Tendon
from .core.shapes import circle, h_section, hollow_rect, l_section, rectangle, rounded, t_section, track

__all__ = [
    "ACI318", "EC2_1_1", "EC2_2", "KDS14_22", "KDS14_22App", "KDS24_25", "CustomCode", "DesignCode", "Pivot", "TopFiberPivot",
    "Bar", "PolygonSection", "SteelRegion", "Tendon", "circle", "h_section", "rounded",
    "rebar_area", "circle_area", "rectangle", "track", "t_section", "l_section", "hollow_rect",
    "FunctionLaw", "StressStrainLaw", "TabularLaw", "elastic_plastic",
    "PMCurve", "pm_curve", "CheckResult", "check_load", "check_section", "lambda_u_deg",
    "Model", "build_model", "load_model", "run_model",
    "plot_curve", "save_csv", "save_outputs",
]
