"""Model input from a configuration file (TOML or JSON).

Section and material units are mm, mm^2 and MPa; loads are given in kN and kN*m (compression positive).
Relative paths inside the file (``csv`` laws, ``[output] dir``) are taken relative to the file.
Results go to ``out/<config file name>/`` (under the working folder) unless ``--out`` or ``[output] dir`` says otherwise.
See ``examples/rect_aci318-25.toml`` (every option is listed there)::

    [code]        type = "aci318-25" | "kds14-22" | "kds14-22-app" | "kds24-25" | "ec2-04-1-1" | "ec2-05-2" | "custom" plus that code's parameters
    [section]     outer, holes, bars, steel, tendons
    [analysis]    lambda_u_deg, n_c, reference, ...      (arguments of pm_curve)
    [loads.NAME]  Pu [kN], Mux [kN*m], Muy [kN*m]        (optional: DCR for each load)
    [output]      dir, basename                          (optional)
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from .check import CheckResult, check_section
from ..codes import ACI318, EC2_1_1, EC2_2, KDS14_22, KDS14_22App, KDS24_25, CustomCode, DesignCode, TopFiberPivot
from .diagram import PMCurve, pm_curve
from .laws import StressStrainLaw, TabularLaw, elastic_plastic
from ..output import save_outputs
from ..section import Bar, PolygonSection, SteelRegion, Tendon
from .rebar import STANDARDS, circle_area, rebar_area
from .shapes import circle, h_section, hollow_rect, l_section, rectangle, t_section, track


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def _take(d: Any, allowed, where: str, required=()) -> Dict[str, Any]:
    """Check that ``d`` is a table with only ``allowed`` keys (and all ``required`` ones)."""
    if not isinstance(d, dict):
        raise ValueError(f"{where}: expected a table, got {type(d).__name__}")
    extra = sorted(set(d) - set(allowed))
    if extra:
        raise ValueError(f"{where}: unknown key {extra[0]!r} (allowed: {', '.join(sorted(allowed))})")
    for k in required:
        if k not in d:
            raise ValueError(f"{where}: missing required key {k!r}")
    return d


def _num(v, where: str) -> float:
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ValueError(f"{where}: expected a finite number, got {v!r}")
    return float(v)


def _point(v, where: str):
    if not (isinstance(v, (list, tuple)) and len(v) == 2):
        raise ValueError(f"{where}: expected [x, y], got {v!r}")
    return (_num(v[0], where), _num(v[1], where))


# ----------------------------------------------------------------------------
# polygons
# ----------------------------------------------------------------------------
# shape name -> (builder, {key: type}, required keys, defaults)
_CONCRETE_SHAPES = {
    "rect": (rectangle, dict(x="num", y="num", width="num", height="num", r="num", n_arc="int"),
             ("width", "height"), dict(x=0.0, y=0.0)),
    "circle": (circle, dict(cx="num", cy="num", r="num", n="int", equal_area="bool"), ("r",), dict(cx=0.0, cy=0.0)),
    "track": (track, dict(cx="num", cy="num", length="num", width="num", rotation_deg="num", n_arc="int"),
              ("length", "width"), dict(cx=0.0, cy=0.0)),
    "t_section": (t_section, dict(x="num", y="num", bf="num", h="num", tf="num", tw="num", rotation_deg="num"),
                  ("bf", "h", "tf", "tw"), dict(x=0.0, y=0.0)),
    "l_section": (l_section, dict(x="num", y="num", width="num", height="num", tx="num", ty="num", r="num",
                                  rotation_deg="num", n_arc="int"), ("width", "height", "tx", "ty"),
                  dict(x=0.0, y=0.0)),
}
_STEEL_SHAPES = {
    "h_section": (h_section, dict(d="num", bf="num", tw="num", tf="num", cx="num", cy="num", rotation_deg="num",
                                  r="num", n_arc="int"), ("d", "bf", "tw", "tf"), {}),
    "hollow_rect": (hollow_rect, dict(cx="num", cy="num", width="num", height="num", t="num", rotation_deg="num"),
                    ("width", "height", "t"), dict(cx=0.0, cy=0.0)),
}


def _shape(spec, where: str, table: Dict[str, Any]):
    """Build a typical section from ``{shape = "...", ...}``; ``table`` says which shapes are allowed here."""
    _take(spec, set(spec) if isinstance(spec, dict) else (), where, required=("shape",))
    name = spec["shape"]
    if name not in table:
        raise ValueError(f"{where}: unknown shape {name!r} (choose from: {', '.join(table)})")
    builder, keys, required, defaults = table[name]
    s = {k: v for k, v in spec.items() if k != "shape"}
    extra = sorted(set(s) - set(keys))
    if extra:
        raise ValueError(f"{where}: shape {name!r} does not use {extra[0]!r} (keys: {', '.join(keys)})")
    for k in required:
        if k not in s:
            raise ValueError(f"{where}: shape {name!r} needs {k!r}")
    kw = dict(defaults)
    for k, v in s.items():
        if keys[k] == "num":
            kw[k] = _num(v, f"{where}.{k}")
        elif keys[k] == "int":
            kw[k] = int(_num(v, f"{where}.{k}"))
        else:
            kw[k] = bool(v)
    return builder(**kw)


def _polygon(spec, where: str):
    """Concrete outline or hole: a list of [x, y] vertices, or a table with ``shape`` =
    rect | circle | track | t_section | l_section."""
    if isinstance(spec, list):
        return [_point(p, f"{where}[{i}]") for i, p in enumerate(spec)]
    return _shape(spec, where, _CONCRETE_SHAPES)


def _steel_polygon(spec, where: str):
    """Steel outline and the voids that come with the shape: ``(outline, [holes])``.

    ``spec`` is a list of vertices or a table with ``shape`` = h_section | hollow_rect.
    """
    if isinstance(spec, list):
        return [_point(p, f"{where}[{i}]") for i, p in enumerate(spec)], []
    out = _shape(spec, where, _STEEL_SHAPES)
    if spec["shape"] == "hollow_rect":
        return out[0], [out[1]]
    return out, []

# ----------------------------------------------------------------------------
# laws
# ----------------------------------------------------------------------------
def _beyond(v, where: str):
    if v is None:
        return None
    if isinstance(v, str):
        return v
    if isinstance(v, list) and len(v) == 2 and all(isinstance(x, str) for x in v):
        return tuple(v)
    raise ValueError(f"{where}: beyond must be a string or a pair of strings, got {v!r}")


def _law(spec, where: str, base_dir: str) -> StressStrainLaw:
    """``{strain=[..], stress=[..]}``, ``{csv="file"}`` or ``{type="elastic_plastic", E=, fy=}``."""
    _take(spec, {"type", "strain", "stress", "csv", "beyond", "delimiter", "E", "fy", "yield_strain"}, where)
    kind = spec.get("type") or ("csv" if "csv" in spec else "table" if "strain" in spec else None)
    if kind == "elastic_plastic":
        _take(spec, {"type", "E", "fy"}, where, required=("E", "fy"))
        return elastic_plastic(_num(spec["E"], where), _num(spec["fy"], where))
    beyond = _beyond(spec.get("beyond"), where)
    if kind == "csv":
        _take(spec, {"type", "csv", "beyond", "delimiter"}, where, required=("csv",))
        path = spec["csv"] if os.path.isabs(spec["csv"]) else os.path.join(base_dir, spec["csv"])
        return TabularLaw.from_csv(path, beyond=beyond, delimiter=spec.get("delimiter", ","))
    if kind in ("table", None) and "strain" in spec:
        _take(spec, {"type", "strain", "stress", "beyond", "yield_strain"}, where, required=("strain", "stress"))
        kw = {}
        if "yield_strain" in spec:
            kw["yield_strain"] = _num(spec["yield_strain"], where)
        return TabularLaw([_num(v, where) for v in spec["strain"]], [_num(v, where) for v in spec["stress"]],
                          beyond=beyond, **kw)
    raise ValueError(f"{where}: give strain/stress lists, a csv file or type = \"elastic_plastic\"")


# ----------------------------------------------------------------------------
# design code
# ----------------------------------------------------------------------------
_ACI_KEYS = {"fc", "fy", "Es", "Fy", "Es_steel", "fpy", "Ep", "transverse", "phi_strain", "eps_cu",
             "gamma_c", "gamma_s"}


_KDS_COMMON = {"fck", "fy", "Es", "Fy", "Es_steel", "fpy", "Ep"}
_KDS_KEYS = {"kds14-22": _KDS_COMMON | {"transverse", "stress_block", "phi_max"},
             "kds14-22-app": _KDS_COMMON | {"phi_c", "phi_s", "stress_block"},
             "kds24-25": _KDS_COMMON | {"phi_c", "phi_s", "combination"}}
_KDS_CLASSES = {"kds14-22": KDS14_22, "kds14-22-app": KDS14_22App, "kds24-25": KDS24_25}
_KDS_TEXT_KEYS = ("transverse", "stress_block", "combination")


_EC2_KEYS = _KDS_COMMON | {"na", "situation", "gamma_c", "gamma_s", "alpha_cc", "stress_block", "steel_branch",
                           "ductility", "eps_ud", "gamma_a"}
_EC2_TEXT_KEYS = ("na", "situation", "stress_block", "steel_branch", "ductility")
_EC2_CLASSES = {"ec2-04-1-1": EC2_1_1, "ec2-05-2": EC2_2}


def build_code(spec: Dict[str, Any], base_dir: str = ".") -> DesignCode:
    kind = str(_take(spec, _ACI_KEYS | set().union(*_KDS_KEYS.values()) | _EC2_KEYS | {"type", "concrete", "rebar", "steel", "tendon",
                                                                            "phi"}, "code",
                     required=("type",))["type"]).lower()
    if kind in _KDS_CLASSES:
        _take(spec, _KDS_KEYS[kind] | {"type"}, "code", required=("fck", "fy"))
        kw = {k: (v if k in _KDS_TEXT_KEYS else _num(v, f"code.{k}")) for k, v in spec.items() if k != "type"}
        return _KDS_CLASSES[kind](**kw)
    if kind in _EC2_CLASSES:
        _take(spec, _EC2_KEYS | {"type"}, "code", required=("fck", "fy"))
        return _EC2_CLASSES[kind](**{k: (v if k in _EC2_TEXT_KEYS else _num(v, f"code.{k}"))
                                     for k, v in spec.items() if k != "type"})
    if kind == "aci318-25":
        _take(spec, _ACI_KEYS | {"type"}, "code")
        kw = {k: v for k, v in spec.items() if k != "type"}
        for k in kw:
            if k not in ("transverse", "phi_strain"):
                kw[k] = _num(kw[k], f"code.{k}")
        return ACI318(**kw)
    if kind == "custom":
        _take(spec, {"type", "concrete", "rebar", "steel", "tendon", "eps_cu", "phi"}, "code", required=("concrete",))
        laws = {k: _law(spec[k], f"code.{k}", base_dir) for k in ("concrete", "rebar", "steel", "tendon") if k in spec}
        return CustomCode(pivot=TopFiberPivot(_num(spec.get("eps_cu", 0.003), "code.eps_cu")),
                          phi=_num(spec.get("phi", 1.0), "code.phi"), **laws)
    raise ValueError(f"code.type: unknown design code {kind!r} (aci318-25, kds14-22, kds14-22-app, kds24-25, ec2-04-1-1, ec2-05-2, custom)")


# ----------------------------------------------------------------------------
# section
# ----------------------------------------------------------------------------
def _expand_points(entry, where: str):
    """One bar/tendon entry -> list of (x, y, extra-dict)."""
    if isinstance(entry, list):
        if len(entry) != 3:
            raise ValueError(f"{where}: expected [x, y, area], got {entry!r}")
        return [(_num(entry[0], where), _num(entry[1], where), {"area": _num(entry[2], where)})]
    if not isinstance(entry, dict):
        raise ValueError(f"{where}: expected [x, y, area] or a table")
    gen_line, gen_circle = "line" in entry, "circle" in entry
    if gen_line and gen_circle:
        raise ValueError(f"{where}: give either line or circle, not both")
    if gen_line:
        e = dict(entry)
        a, b = (_point(p, where) for p in _line_ends(e.pop("line"), where))
        n = int(e.pop("n", 0))
        if n < 1:
            raise ValueError(f"{where}: line needs n >= 1")
        ts = [0.5] if n == 1 else [i / (n - 1) for i in range(n)]
        return [(a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), e) for t in ts]
    if gen_circle:
        e = dict(entry)
        cx, cy, r = e.pop("circle")
        n = int(e.pop("n", 0))
        start = _num(e.pop("start_deg", 0.0), where)
        if n < 1:
            raise ValueError(f"{where}: circle needs n >= 1")
        return [(_num(cx, where) + _num(r, where) * math.cos(math.radians(start + 360.0 * i / n)),
                 _num(cy, where) + _num(r, where) * math.sin(math.radians(start + 360.0 * i / n)), dict(e))
                for i in range(n)]
    e = dict(entry)
    if "x" not in e or "y" not in e:
        raise ValueError(f"{where}: needs x and y (or line / circle)")
    x, y = _num(e.pop("x"), where), _num(e.pop("y"), where)
    return [(x, y, e)]


def _line_ends(v, where):
    if not (isinstance(v, list) and len(v) == 2):
        raise ValueError(f"{where}: line must be [[x1, y1], [x2, y2]]")
    return v


def _bar_area(e: Dict[str, Any], where: str, standard: str) -> float:
    """Area of one bar/tendon entry from exactly one of ``area`` [mm^2], ``dia`` [mm] or ``size`` (table)."""
    given = [k for k in ("area", "dia", "size") if k in e]
    if len(given) != 1:
        raise ValueError(f"{where}: give exactly one of area, dia, size (got {', '.join(given) or 'none'})")
    k = given[0]
    if k == "area":
        return _num(e["area"], where + ".area")
    if k == "dia":
        return circle_area(_num(e["dia"], where + ".dia"))
    if standard == "custom":
        raise ValueError(f"{where}: 'size' needs section.rebar_standard = \"KS\" | \"ASTM\" | \"ASTM_SI\" | \"BS\"; "
                         f"with \"custom\" give area or dia")
    try:
        return rebar_area(standard, e["size"])
    except ValueError as err:
        raise ValueError(f"{where}: {err}") from None


def _bars(specs, base_dir: str, standard: str = "custom") -> List[Bar]:
    out = []
    for i, entry in enumerate(specs):
        where = f"section.bars[{i}]"
        for x, y, e in _expand_points(entry, where):
            _take(e, {"area", "dia", "size", "law"}, where)
            law = _law(e["law"], where + ".law", base_dir) if "law" in e else None
            out.append(Bar(x, y, _bar_area(e, where, standard), law))
    return out


def _tendons(specs, base_dir: str) -> List[Tendon]:
    out = []
    for i, entry in enumerate(specs):
        where = f"section.tendons[{i}]"
        for x, y, e in _expand_points(entry, where):
            _take(e, {"area", "dia", "fse", "eps_pe", "duct_area", "law"}, where)
            law = _law(e["law"], where + ".law", base_dir) if "law" in e else None
            kw = {k: _num(e[k], where) for k in ("fse", "eps_pe", "duct_area") if k in e}
            out.append(Tendon(x, y, _bar_area(e, where, "custom"), law=law, **kw))
    return out

def build_section(spec: Dict[str, Any], base_dir: str = ".") -> PolygonSection:
    _take(spec, {"outer", "holes", "bars", "steel", "tendons", "rebar_standard"}, "section", required=("outer",))
    standard = str(spec.get("rebar_standard", "custom"))
    standard = "custom" if standard.lower() == "custom" else standard.upper()
    if standard not in STANDARDS:
        raise ValueError(f"section.rebar_standard: unknown standard {spec['rebar_standard']!r} (choose from: {', '.join(STANDARDS)})")
    outer = _polygon(spec["outer"], "section.outer")
    holes = [_polygon(h, f"section.holes[{i}]") for i, h in enumerate(spec.get("holes", []))]
    steel = []
    for i, s in enumerate(spec.get("steel", [])):
        where = f"section.steel[{i}]"
        _take(s, {"polygon", "holes", "law"}, where, required=("polygon",))
        law = _law(s["law"], where + ".law", base_dir) if "law" in s else None
        poly, shape_holes = _steel_polygon(s["polygon"], where + ".polygon")
        extra_holes = [_polygon(h, f"{where}.holes[{j}]") for j, h in enumerate(s.get("holes", []))]
        steel.append(SteelRegion(poly, shape_holes + extra_holes, law))
    return PolygonSection(outer, holes=holes, bars=_bars(spec.get("bars", []), base_dir, standard), steel=steel,
                          tendons=_tendons(spec.get("tendons", []), base_dir))


# ----------------------------------------------------------------------------
# whole model
# ----------------------------------------------------------------------------
_ANALYSIS_KEYS = {"lambda_u_deg", "n_c", "c_min_ratio", "c_max_ratio", "reference", "tol_deg", "scan_deg",
                   "full_scan"}


KN = 1.0e3          # kN -> N
KNM = 1.0e6         # kN*m -> N*mm


@dataclass
class Load:
    """A factored load in N and N*mm (the config file gives kN and kN*m)."""
    name: str
    Pu: float
    Mux: float
    Muy: float


@dataclass
class Model:
    section: PolygonSection
    code: DesignCode
    analysis: Dict[str, Any] = field(default_factory=dict)     # keyword arguments of pm_curve
    lambda_u_deg: Optional[float] = None
    loads: List[Load] = field(default_factory=list)
    out_dir: str = "out"
    basename: str = "pm"


def build_model(data: Dict[str, Any], base_dir: str = ".", name: str = "pm") -> Model:
    _take(data, {"code", "section", "analysis", "loads", "output"}, "config", required=("code", "section"))
    code = build_code(data["code"], base_dir)
    section = build_section(data["section"], base_dir)

    a = dict(_take(data.get("analysis", {}), _ANALYSIS_KEYS, "analysis"))
    lam = a.pop("lambda_u_deg", None)
    if lam is not None:
        lam = _num(lam, "analysis.lambda_u_deg")
    if "reference" in a and not isinstance(a["reference"], str):
        a["reference"] = _point(a["reference"], "analysis.reference")

    # loads: {name: {Pu [kN], Mux [kN*m], Muy [kN*m]}}; stored in N and N*mm
    raw_loads = data.get("loads", {})
    if not isinstance(raw_loads, dict):
        raise ValueError("loads: expected a table of named loads, e.g. [loads.L1] (a list is not accepted)")
    loads = []
    for load_name, ld in raw_loads.items():
        where = f"loads.{load_name}"
        _take(ld, {"Pu", "Mux", "Muy"}, where, required=("Pu", "Mux", "Muy"))
        loads.append(Load(str(load_name),_num(ld["Pu"], where) * KN, _num(ld["Mux"], where) * KNM,
                          _num(ld["Muy"], where) * KNM))
    if not loads and lam is None:
        raise ValueError("config: give analysis.lambda_u_deg or at least one load, e.g. [loads.L1]")

    out = _take(data.get("output", {}), {"dir", "basename"}, "output")
    if "dir" in out:
        out_dir = out["dir"]
        if not os.path.isabs(out_dir):
            out_dir = os.path.join(base_dir, out_dir)
    else:
        out_dir = os.path.abspath(os.path.join("out", name))     # out/<config file name>/ under the working folder
    return Model(section, code, a, lam, loads, out_dir, str(out.get("basename", name)))


def load_model(path) -> Model:
    """Read a ``.toml`` or ``.json`` file and build the section, design code and analysis settings."""
    path = os.fspath(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".toml":
        import tomllib
        with open(path, encoding="utf-8-sig") as f:        # tolerate the BOM some Windows editors add
            data = tomllib.loads(f.read())
    elif ext == ".json":
        with open(path, encoding="utf-8-sig") as f:
            data = json.load(f)
    else:
        raise ValueError(f"unsupported config file type {ext!r} (use .toml or .json)")
    return build_model(data, os.path.dirname(os.path.abspath(path)), os.path.splitext(os.path.basename(path))[0])


def run_model(model: Model, out_dir: Optional[str] = None, save: bool = True):
    """Compute the curve (no loads) or the DCR of every load; write PNG and CSV when ``save`` is true.

    Returns a list of :class:`PMCurve` or :class:`CheckResult`.
    """
    out_dir = out_dir or model.out_dir
    results: list = []
    if model.loads:
        for ld in model.loads:
            res = check_section(model.section, model.code, ld.Pu, ld.Mux, ld.Muy, **model.analysis)
            results.append(res)
            if save:
                save_outputs(res, out_dir, f"{model.basename}_{ld.name}")
    else:
        curve = pm_curve(model.section, model.code, model.lambda_u_deg, **model.analysis)
        results.append(curve)
        if save:
            save_outputs(curve, out_dir, model.basename)
    return results
