"""Section drawing and PM curve plot (matplotlib is imported lazily and never configured globally)."""
from __future__ import annotations

import math

import numpy as np

# material colours
CONCRETE = "#cfcfcf"
STEEL = "#4f81bd"
REBAR = "#c0392b"
TENDON = "#2e8b57"


def _path(rings):
    """One compound path: counter-clockwise outline plus clockwise holes fill as a ring."""
    from matplotlib.path import Path

    verts, codes = [], []
    for r in rings:
        r = np.asarray(r, dtype=float)
        verts.extend(r.tolist() + [r[0].tolist()])
        codes.extend([Path.MOVETO] + [Path.LINETO] * (len(r) - 1) + [Path.CLOSEPOLY])
    return Path(verts, codes)


def draw_section(ax, section, ref_point=None, extra_handles=()):
    """Fill the section by material: concrete (holes left open), steel, rebar, tendons."""
    from matplotlib.patches import Circle, PathPatch

    ax.add_patch(PathPatch(_path([section.outer] + list(section.holes)), facecolor=CONCRETE,
                           edgecolor="#555555", lw=1.0, zorder=1))
    for s in section.steel:
        ax.add_patch(PathPatch(_path([s.polygon] + list(s.holes)), facecolor=STEEL, edgecolor="#1f3b63",
                               lw=0.8, zorder=2))
    r_min = 0.004 * section.diag                         # keep very small bars visible
    for items, colour in ((section.bars, REBAR), (section.tendons, TENDON)):
        for b in items:
            r = max(math.sqrt(b.area / math.pi), r_min)
            ax.add_patch(Circle((b.x, b.y), r, facecolor=colour, edgecolor="#222222", lw=0.5, zorder=3))
    if ref_point is not None:
        ax.plot([ref_point[0]], [ref_point[1]], "+", color="k", ms=10, mew=1.5, zorder=4)

    handles = [_swatch(CONCRETE, "concrete")]
    if section.steel:
        handles.append(_swatch(STEEL, "steel"))
    if section.bars:
        handles.append(_swatch(REBAR, f"rebar ({len(section.bars)})"))
    if section.tendons:
        handles.append(_swatch(TENDON, f"tendon ({len(section.tendons)})"))
    handles += list(extra_handles)
    from matplotlib.transforms import ScaledTranslation
    below = ax.transAxes + ScaledTranslation(0, -40 / 72, ax.figure.dpi_scale_trans)   # 40 pt under the axes box
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), bbox_transform=below, borderaxespad=0.0, ncol=len(handles), fontsize=8,
              frameon=False)
    x0, y0, x1, y1 = section.bbox
    pad = 0.06 * max(x1 - x0, y1 - y0)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title("section  (+ = moment reference)")
    ax.grid(alpha=0.25)


def _swatch(colour, label):
    from matplotlib.patches import Patch

    return Patch(facecolor=colour, edgecolor="#333333", label=label)


def plot_curve(obj, path=None, title=None):
    from matplotlib.figure import Figure

    from ..core.check import CheckResult

    curve = obj.curve if isinstance(obj, CheckResult) else obj
    section = getattr(curve, "section", None)
    if section is not None:
        fig = Figure(figsize=(13, 6))
        ax_s, ax = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1.0, 1.15]})
        draw_section(ax_s, section, getattr(curve, "ref_point", None))
    else:
        fig = Figure(figsize=(7, 6))
        ax = fig.subplots()
    ax.plot(curve.Mn / 1e6, curve.Pn / 1e3, "--", color="tab:gray", label="Pn - Mn (nominal)")
    ax.plot(curve.phiMn / 1e6, curve.phiPn / 1e3, "-", color="tab:blue", lw=2, label="phiPn - phiMn (design)")
    if isinstance(obj, CheckResult):
        ax.plot([obj.Mu / 1e6], [obj.Pu / 1e3], "o", color="tab:red", label=f"load  DCR = {obj.dcr:.3f}")
        if obj.capacity is not None:
            ax.plot([0, obj.capacity[0] / 1e6], [0, obj.capacity[1] / 1e3], ":", color="tab:red")
            ax.plot([obj.capacity[0] / 1e6], [obj.capacity[1] / 1e3], "x", color="tab:red")
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("M along lambda_u  [kN*m]")
    ax.set_ylabel("P  [kN]  (compression +)")
    ax.set_title(title or f"PM interaction, lambda_u = {curve.lambda_u_deg:.2f} deg")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight")
    return fig
