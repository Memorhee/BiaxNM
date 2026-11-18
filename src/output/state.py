"""Figure of the strain/stress state at the capacity point of a load.

Top row: section with the neutral axis and the compressed zone; strain diagram; stress diagram (concrete on the top
axis, reinforcement and steel on the bottom axis) with the resultants C and T and the lever arm z.
Bottom row: the resultants as force x arm sums, and a card with the load, the strengths and the DCR.
"""
from __future__ import annotations

import math

import numpy as np

from ..section import clip_halfplane
from .plot import CONCRETE, REBAR, STEEL, TENDON, _path, draw_section

COMP = "#1f5fa8"            # compression (blue)
TENS = "#c0392b"            # tension (red)
NA = "#d35400"              # neutral axis, compressed zone (orange)
INK = "#1f2a44"             # dark text, table header
PAPER = "#f4f6fa"           # card and table stripe
LINE = "#d5dbe8"


def _compressed_start(integ, a, k, d_lo, d_hi):
    """Smallest d from which the concrete carries compressive stress up to the extreme fibre."""
    d = np.linspace(d_hi, d_lo, 4001)
    sig = integ.conc.stress(a + k * d)
    pos = sig > 0
    if pos.all():
        return d_lo
    if not pos[0]:
        return d_hi
    return float(d[np.argmin(pos)])               # first fibre (going down) without compressive stress


def _fit_axis(ax, vals, pad=0.12):
    """Symmetric x limits around zero (so zero is always in the middle: compression right, tension left); every
    axis has its own scale."""
    m = max(abs(v) for v in vals) or 1.0
    ax.set_xlim(-(1.0 + pad) * m, (1.0 + pad) * m)


def _resultants(st, cth, sth, d_max):
    """Total compression C and tension T of the components and their depths below the extreme compression fibre."""
    out = {}
    for key, sel in (("C", lambda f: f > 0), ("T", lambda f: f < 0)):
        cs = [c for c in st.components if sel(c.F)]
        F = sum(c.F for c in cs)
        if cs and F != 0.0:
            d = sum(c.F * (c.x * cth + c.y * sth) for c in cs) / F
            out[key] = (F, d_max - d)
    return out


def _style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#8a93a6")
    ax.tick_params(labelsize=9, colors="#333b4d")
    ax.grid(alpha=0.25)


def _draw_table(ax, st, xr, yr, caption):
    """Resultants as force x arm: header, one row per component (C blue, T red), bold sum row."""
    ax.text(0.0, 1.0, caption, transform=ax.transAxes, fontsize=9, color="#5b6478", va="bottom", ha="left")
    head = ["component", "F\n[kN]", "x\n[mm]", "y\n[mm]", "arm x − xr\n[mm]", "arm y − yr\n[mm]",
            "My = F·ax\n[kN·m]", "Mx = F·ay\n[kN·m]"]
    rows, colours = [], []
    sP = sMx = sMy = 0.0
    for c in st.components:
        mx, my = c.F * (c.y - yr) / 1e6, c.F * (c.x - xr) / 1e6
        sP += c.F / 1e3
        sMx += mx
        sMy += my
        rows.append([c.name, f"{c.F / 1e3:,.1f}", f"{c.x:,.1f}", f"{c.y:,.1f}", f"{c.x - xr:,.1f}",
                     f"{c.y - yr:,.1f}", f"{my:,.1f}", f"{mx:,.1f}"])
        colours.append(COMP if c.F > 0 else TENS)
    rows.append(["sum  (Pn, Mn)", f"{sP:,.1f}", "", "", "", "", f"{sMy:,.1f}", f"{sMx:,.1f}"])
    tbl = ax.table(cellText=rows, colLabels=head, loc="upper left", cellLoc="right", colLoc="right",
                   bbox=[0.0, 0.0, 1.0, 0.92], colWidths=[0.16, 0.11, 0.1, 0.1, 0.12, 0.12, 0.145, 0.145])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    n_rows = len(rows)
    for (r, col), cell in tbl.get_celld().items():
        cell.set_edgecolor(LINE)                             # closed edges: a cell is only filled with its full outline
        cell.set_linewidth(0.8)
        cell.PAD = 0.04
        if col == 0:
            cell.set_text_props(ha="left")
            cell._loc = "left"
        if r == 0:
            cell.set_facecolor(INK)
            cell.set_edgecolor(INK)
            cell.set_text_props(color="white", fontweight="bold")
        elif r == n_rows:                                    # sum row
            cell.set_facecolor("#dfe6f3")
            cell.set_text_props(fontweight="bold", color=INK)
        else:
            cell.set_facecolor(PAPER if r % 2 == 0 else "white")
            cell.set_text_props(color=colours[r - 1] if col == 0 else "#1c2333",
                                fontweight="bold" if col == 0 else "normal")
    return math.hypot(sMx, sMy)


def _draw_card(ax, result, st, Mn_nom, factored):
    """Result card: DCR badge, load, nominal and design strength."""
    from matplotlib.patches import FancyBboxPatch

    ax.add_patch(FancyBboxPatch((0.0, 0.0), 1.0, 1.0, boxstyle="round,pad=0,rounding_size=0.035", transform=ax.transAxes,
                                facecolor=PAPER, edgecolor=LINE, lw=1.2))
    ok = result.dcr <= 1.0
    badge = "#2e8b57" if ok else "#c0392b"
    ax.add_patch(FancyBboxPatch((0.05, 0.70), 0.90, 0.24, boxstyle="round,pad=0,rounding_size=0.06", transform=ax.transAxes,
                                facecolor=badge, edgecolor="none"))
    ax.text(0.5, 0.865, f"DCR  {result.dcr:.3f}", transform=ax.transAxes, ha="center", va="center", color="white",
            fontsize=21, fontweight="bold")
    ax.text(0.5, 0.745, "OK" if ok else "NG  (over capacity)", transform=ax.transAxes, ha="center", va="center",
            color="white", fontsize=10, fontweight="bold")
    kn, knm = 1e3, 1e6
    items = [
        ("Load Pu", f"{result.Pu / kn:,.0f} kN"),
        ("Load Mux / Muy", f"{result.Mux / knm:,.0f} / {result.Muy / knm:,.0f} kN·m"),
        ("Load Mu", f"{result.Mu / knm:,.0f} kN·m"),
        None,
        ("Pn  (nominal)", f"{sum(c.F for c in st.components) / kn:,.0f} kN"),
        ("Mn  (nominal)", f"{Mn_nom:,.0f} kN·m"),
        ("phi", f"{st.phi:.3f}" + ("  material factors" if factored else "")),
        ("phi Pn", f"{st.P_design / kn:,.0f} kN" + ("  (axial cap)" if st.capped else "")),
        ("phi Mn", f"{st.phi * Mn_nom:,.0f} kN·m"),
    ]
    y = 0.62
    for it in items:
        if it is None:
            ax.plot([0.07, 0.93], [y + 0.012, y + 0.012], transform=ax.transAxes, color=LINE, lw=1.2)
            y -= 0.03
            continue
        ax.text(0.07, y, it[0], transform=ax.transAxes, fontsize=10, color="#5b6478", va="center")
        ax.text(0.93, y, it[1], transform=ax.transAxes, fontsize=10.5, color=INK, fontweight="bold", va="center",
                ha="right")
        y -= 0.075
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def plot_state(result, path=None, title=None):
    """Draw the state of a :class:`~biaxnm.core.check.CheckResult`; ``None`` if the load has no capacity point."""
    from matplotlib.figure import Figure
    from matplotlib.patches import Circle, PathPatch
    from matplotlib.transforms import blended_transform_factory

    from ..analysis.state import solve_state

    st = solve_state(result)
    if st is None:
        return None
    curve = result.curve
    section, integ = curve.section, st.integ
    cth, sth = math.cos(st.theta), math.sin(st.theta)
    n = np.array([cth, sth])
    g = integ.geometry(st.theta)
    d_na = g.d_max - st.c
    d_s = _compressed_start(integ, st.a, st.k, g.d_min, g.d_max)

    fig = Figure(figsize=(15, 10.8), facecolor="white")
    gs = fig.add_gridspec(2, 3, height_ratios=[3.4, 1.7], width_ratios=[1.15, 0.85, 0.95], hspace=0.34, wspace=0.30,
                          top=0.855, bottom=0.04, left=0.06, right=0.975)
    ax_s = fig.add_subplot(gs[0, 0])
    ax_e = fig.add_subplot(gs[0, 1])
    ax_f = fig.add_subplot(gs[0, 2], sharey=ax_e)
    ax_t = fig.add_subplot(gs[1, 0:2])
    ax_k = fig.add_subplot(gs[1, 2])
    ax_t.axis("off")

    # ---- section, neutral axis, compressed zone
    from matplotlib.patches import Patch
    draw_section(ax_s, section, st.ref, extra_handles=[Patch(facecolor="white", edgecolor=NA, hatch="////",
                                                              label="compressed zone")])
    rings = [section.outer] + list(section.holes)
    clipped = [clip_halfplane(r, n, d_s) for r in rings]
    if clipped[0] is not None:
        ax_s.add_patch(PathPatch(_path([r for r in clipped if r is not None]), facecolor="none",
                                 edgecolor=NA, hatch="////", lw=0.0, zorder=2.5))
    x0, y0, x1, y1 = section.bbox
    tvec = np.array([-sth, cth])
    base = n * d_na + tvec * (0.5 * np.array([x0 + x1, y0 + y1]) @ tvec)
    p1, p2 = base - tvec * section.diag, base + tvec * section.diag
    ax_s.plot([p1[0], p2[0]], [p1[1], p2[1]], "-", color=NA, lw=2.0, zorder=5)
    ctr = 0.5 * np.array([x0 + x1, y0 + y1])
    tip = ctr + n * 0.32 * section.diag
    ax_s.annotate("", xy=tip, xytext=ctr, arrowprops=dict(arrowstyle="-|>", color=NA, lw=1.8), zorder=6)
    ax_s.text(tip[0], tip[1], "  compression", color=NA, fontsize=9, fontweight="bold", va="center")
    ax_s.set_title("")                                   # draw_section's own title; the common titles come below
    titles = [(ax_s, f"θ = {math.degrees(st.theta):.1f}°,   c = {st.c:,.1f} mm")]     # placed together below
    ax_s.set_xlim(x0 - 0.12 * (x1 - x0), x1 + 0.12 * (x1 - x0))
    ax_s.set_ylim(y0 - 0.12 * (y1 - y0), y1 + 0.12 * (y1 - y0))
    ax_s.tick_params(labelsize=9)

    # ---- strain diagram versus depth from the extreme compression fibre
    d_lo = min(g.d_rein, g.d_min)
    d_grid = np.linspace(g.d_max, d_lo, 400)
    depth = g.d_max - d_grid
    eps = st.a + st.k * d_grid
    e3 = eps * 1e3
    ax_e.fill_betweenx(depth, 0, e3, where=e3 >= 0, color=COMP, alpha=0.18, lw=0)
    ax_e.fill_betweenx(depth, 0, e3, where=e3 <= 0, color=TENS, alpha=0.18, lw=0)
    ax_e.plot(e3, depth, color=INK, lw=2.2)
    ax_e.axvline(0, color="k", lw=0.8)
    ax_e.annotate(f"{e3[0]:.2f}", (e3[0], 0), textcoords="offset points", xytext=(6, -14), color=COMP, fontsize=10,
                  fontweight="bold")
    ax_e.annotate(f"{e3[-1]:.2f}", (e3[-1], depth[-1]), textcoords="offset points", xytext=(6, 6),
                  color=TENS if e3[-1] < 0 else COMP, fontsize=10, fontweight="bold")
    titles.append((ax_e, "strain [‰]   (compression +)"))
    ax_e.set_ylabel("depth below the top compression fibre [mm]", fontsize=10)
    _fit_axis(ax_e, list(e3))

    # ---- stress diagram: concrete on the top axis, reinforcement and steel on the bottom axis (their stresses differ
    # by a factor of 10 to 50, so a shared axis would flatten the concrete profile)
    sig_c = integ.conc.stress(eps)
    ax_c = ax_f.twiny()
    ax_c.fill_betweenx(depth, 0, sig_c, color="#9db4cf", alpha=0.65, edgecolor=COMP, lw=1.6, label="concrete")
    ax_f.set_zorder(ax_c.get_zorder() + 1)
    ax_f.patch.set_visible(False)
    ax_f.axvline(0, color="k", lw=0.8)
    s_vals = [0.0]
    ring_marks = []                                      # (x, y, radius, sign) for the section panel
    for items_x, items_y, areas, colour, label in ((integ.bar_x, integ.bar_y, integ.bar_A, REBAR, "rebar"),
                                                   (integ.t_x, integ.t_y, integ.t_A, TENDON, "tendon")):
        if len(items_x):
            d = items_x * cth + items_y * sth
            e = st.a + st.k * d
            if label == "tendon":
                e = e - integ.t_pe
            law_sig = np.empty_like(e)
            groups = integ.bar_groups if label == "rebar" else integ.t_groups
            for law, idx in groups:
                law_sig[idx] = law.stress(e[idx])
            edge = [COMP if s > 0 else TENS for s in law_sig]
            ax_f.scatter(law_sig, g.d_max - d, s=46, color=colour, edgecolor=edge, linewidths=1.8, zorder=3,
                         label=label)
            s_vals += list(law_sig)
            hi, lo = int(np.argmax(law_sig)), int(np.argmin(law_sig))
            for q in {hi, lo}:
                ax_f.annotate(f"{law_sig[q]:.0f}", (law_sig[q], g.d_max - d[q]), textcoords="offset points",
                              xytext=(7, 4), fontsize=9, fontweight="bold", color=edge[q])
            for xx, yy, aa, ee in zip(items_x, items_y, areas, edge):
                ring_marks.append((xx, yy, math.sqrt(aa / math.pi), ee))
    for _, law, poly in integ.steel:
        d = poly[:, 0] * cth + poly[:, 1] * sth
        ds = np.linspace(d.max(), d.min(), 100)
        s_steel = law.stress(st.a + st.k * ds)
        ax_f.fill_betweenx(g.d_max - ds, 0, s_steel, where=s_steel >= 0, color=COMP, alpha=0.15, lw=0)
        ax_f.fill_betweenx(g.d_max - ds, 0, s_steel, where=s_steel <= 0, color=TENS, alpha=0.15, lw=0)
        ax_f.plot(s_steel, g.d_max - ds, color=STEEL, lw=2.4, label="steel")
        s_vals += list(s_steel)
    top = float(sig_c[0])
    ax_c.annotate(f"{top:.1f} MPa", (top, 0), textcoords="offset points", xytext=(6, -14), fontsize=9.5,
                  fontweight="bold", color=COMP)
    _fit_axis(ax_c, list(sig_c))
    _fit_axis(ax_f, s_vals)
    ax_c.set_xlabel("concrete [MPa]", fontsize=10, color=COMP)
    ax_c.tick_params(axis="x", colors=COMP, labelsize=9)
    ax_f.set_xlabel("rebar / steel / tendon [MPa]", fontsize=10)
    titles.append((ax_f, "stress   (compression +)"))
    handles = ax_c.get_legend_handles_labels()[0] + ax_f.get_legend_handles_labels()[0]
    labels = ax_c.get_legend_handles_labels()[1] + ax_f.get_legend_handles_labels()[1]
    ax_f.legend(handles, labels, fontsize=9, loc="lower right", frameon=True, framealpha=0.9, edgecolor=LINE)
    plt_depth_max = float(depth[-1])
    for ax in (ax_e, ax_f):
        ax.axhline(st.c, color=NA, lw=1.4, ls="--")
        _style_axes(ax)
        ax.set_ylim(plt_depth_max * 1.04, -0.04 * plt_depth_max)
    if st.c <= plt_depth_max:
        ax_e.text(0.02, st.c, "  N.A.", color=NA, fontsize=9, fontweight="bold", va="bottom",
                  transform=ax_e.get_yaxis_transform())
    ax_f.tick_params(labelleft=False)

    # ring around every bar in the section panel: blue = compressed, red = in tension (as in the stress panel)
    for xx, yy, rr, ee in ring_marks:
        ax_s.add_patch(Circle((xx, yy), 1.55 * rr, facecolor="none", edgecolor=ee, lw=1.8, zorder=3.5))

    # resultants C and T with the lever arm z, drawn on the stress diagram (x in axes fraction, y = depth)
    res = _resultants(st, cth, sth, g.d_max)
    tr = blended_transform_factory(ax_f.transAxes, ax_f.transData)
    label_box = dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.85)
    for key, colour, x_end, ha in (("C", COMP, 0.93, "right"), ("T", TENS, 0.07, "left")):
        if key in res:
            F, dd = res[key]
            ax_f.annotate("", xy=(x_end, dd), xytext=(0.5, dd), xycoords=tr, textcoords=tr,
                          arrowprops=dict(arrowstyle="-|>", color=colour, lw=2.6, mutation_scale=16), zorder=6)
            x_txt = 0.5 + 0.5 * (x_end - 0.5)                # label above the middle of the arrow, clear of the head
            ax_f.annotate(f"{key} = {abs(F) / 1e3:,.0f} kN", xy=(x_txt, dd), xycoords=tr, xytext=(0, 7),
                          textcoords="offset points", color=colour, fontsize=10, fontweight="bold", ha="center",
                          va="bottom", bbox=label_box, zorder=7)
    if "C" in res and "T" in res:
        (_, dc), (_, dt) = res["C"], res["T"]
        ax_f.annotate("", xy=(0.5, dt), xytext=(0.5, dc), xycoords=tr, textcoords=tr,
                      arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=1.3, ls="--", mutation_scale=11), zorder=6)
        ax_f.text(0.52, 0.5 * (dc + dt), f"z = {abs(dt - dc):,.0f} mm", transform=tr, color=INK, fontsize=10,
                  fontweight="bold", va="center", ha="left", bbox=label_box, zorder=7)

    # ---- table and result card
    xr, yr = st.ref
    factored = hasattr(curve.code, "phi_c") and hasattr(curve.code, "phi_s")      # material-factor design
    kind = "Design strength (material factors applied)" if factored else "Nominal strength"
    caption = (f"{kind} at the capacity point of the load  ·  ε_cu at the extreme compression fibre  ·  "
               f"moment reference ({xr + 0.0:.1f}, {yr + 0.0:.1f}) mm"
               + ("" if st.exact else "  ·  approximate: nearest curve point"))
    Mn_nom = _draw_table(ax_t, st, xr, yr, caption)
    _draw_card(ax_k, result, st, Mn_nom, factored)

    # panel titles on one common baseline: above the highest panel and the tick labels of the top (concrete) axis
    ax_s.set_anchor("N")                                 # a wide section stays under its title, not centred
    ax_s.apply_aspect()                                  # the equal-aspect section panel may be shorter than the others
    y_top = max(ax.get_position().y1 for ax, _ in titles)
    for ax, text in titles:
        p = ax.get_position()
        fig.text(0.5 * (p.x0 + p.x1), y_top + 0.055, text, ha="center", va="bottom", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.suptitle(title or "Strain and stress state at the capacity point", fontsize=19, fontweight="bold", color=INK,
                 y=0.985)
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight")
    return fig
