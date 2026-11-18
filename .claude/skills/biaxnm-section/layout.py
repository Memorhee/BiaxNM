"""Rebar layout candidates for a rectangular or circular column (helper of the biaxnm-section skill).

    python .claude/skills/biaxnm-section/layout.py rect   --width 400 --height 600 --cover 40 --tie 10 --spacing 150 --ratio 1.7
    python .claude/skills/biaxnm-section/layout.py circle --diameter 600 --cover 40 --tie 10 --ratio 2.0

Units mm.  Perimeter bars only.  Rectangle: lower-left corner at (0, 0), as the default of ``shape = "rect"``.
Circle: centre at (0, 0).  ``--spacing`` is the largest allowed centre-to-centre spacing; ``--ratio`` is the target
steel ratio [%] of the gross area (the diameter is varied to approach it, it is never hit exactly).
Prints candidates and TOML ``bars`` lines.  Bar areas come from ``src/core/rebar.py`` (same tables as the program).
"""
import argparse
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("rebar", os.path.join(HERE, "..", "..", "..", "src", "core", "rebar.py"))
rebar = importlib.util.module_from_spec(_spec)
sys.modules["rebar"] = rebar
_spec.loader.exec_module(rebar)

DEFAULT_SIZES = {"KS": "D13 D16 D19 D22 D25 D29 D32 D35", "ASTM": "#5 #6 #7 #8 #9 #10 #11",
                 "ASTM_SI": "#16 #19 #22 #25 #29 #32 #36", "BS": "12 16 20 25 32 40"}


def clear_flag(clear, db):
    # reference only: the usual minimum clear spacing max(1.5 db, 40 mm); not checked against the code text
    return " (clear < max(1.5db, 40): check)" if clear < max(1.5 * db, 40.0) else ""


def rect_candidates(a):
    out = []
    for size in a.sizes:
        area = rebar.rebar_area(a.standard, size)
        db = math.sqrt(4 * area / math.pi)
        c = a.cover + a.tie + db / 2
        w, h = a.width - 2 * c, a.height - 2 * c
        if w <= 0 or h <= 0:
            continue
        if a.spacing:
            grids = [(max(1, math.ceil(w / a.spacing - 1e-9)), max(1, math.ceil(h / a.spacing - 1e-9)))]
        else:
            grids = [(nx, ny) for nx in range(1, 25) for ny in range(1, 25)]
        best = None
        for nx, ny in grids:
            n = 2 * (nx + ny)
            ratio = 100 * n * area / (a.width * a.height)
            sx, sy = w / nx, h / ny
            key = (abs(ratio - a.ratio) if a.ratio else 0.0, abs(sx - sy))
            if best is None or key < best[0]:
                best = (key, nx, ny, n, ratio, sx, sy)
        _, nx, ny, n, ratio, sx, sy = best
        toml = [f'  {{ line = [[{c:g}, {c:g}], [{a.width - c:g}, {c:g}]], n = {nx + 1}, size = "{size}" }},   # bottom',
                f'  {{ line = [[{c:g}, {a.height - c:g}], [{a.width - c:g}, {a.height - c:g}]], n = {nx + 1}, size = "{size}" }},   # top']
        if ny > 1:
            toml += [f'  {{ line = [[{c:g}, {c + sy:g}], [{c:g}, {a.height - c - sy:g}]], n = {ny - 1}, size = "{size}" }},   # left',
                     f'  {{ line = [[{a.width - c:g}, {c + sy:g}], [{a.width - c:g}, {a.height - c - sy:g}]], n = {ny - 1}, size = "{size}" }},   # right']
        out.append(dict(size=size, db=db, c=c, n=n, ratio=ratio, sp=(sx, sy), clear=min(sx, sy) - db, toml=toml,
                        info=f"{nx} x {ny} segments, spacing x {sx:.1f} / y {sy:.1f} mm"))
    return out


def circle_candidates(a):
    out = []
    for size in a.sizes:
        area = rebar.rebar_area(a.standard, size)
        db = math.sqrt(4 * area / math.pi)
        c = a.cover + a.tie + db / 2
        r = a.diameter / 2 - c
        if r <= 0:
            continue
        if a.spacing:
            ns = [max(4, math.ceil(2 * math.pi * r / a.spacing - 1e-9))]
        else:
            ns = range(4, 81)
        gross = math.pi * a.diameter ** 2 / 4
        n = min(ns, key=lambda k: abs(100 * k * area / gross - a.ratio) if a.ratio else 0.0)
        sp = 2 * math.pi * r / n
        out.append(dict(size=size, db=db, c=c, n=n, ratio=100 * n * area / gross, sp=(sp, sp),
                        clear=2 * r * math.sin(math.pi / n) - db,
                        toml=[f'  {{ circle = [0, 0, {r:g}], n = {n}, size = "{size}" }},'],
                        info=f"{n} bars on r = {r:.1f} mm, arc spacing {sp:.1f} mm"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("shape", choices=["rect", "circle"])
    ap.add_argument("--width", type=float)
    ap.add_argument("--height", type=float)
    ap.add_argument("--diameter", type=float)
    ap.add_argument("--cover", type=float, required=True, help="clear cover to the tie [mm]")
    ap.add_argument("--tie", type=float, required=True, help="tie / spiral bar diameter [mm]")
    ap.add_argument("--spacing", type=float, help="largest centre-to-centre spacing [mm]")
    ap.add_argument("--ratio", type=float, help="target steel ratio [%%]")
    ap.add_argument("--standard", default="KS", choices=["KS", "ASTM", "ASTM_SI", "BS"])
    ap.add_argument("--sizes", help="comma separated designations (default: common sizes of the standard)")
    a = ap.parse_args()
    if not (a.spacing or a.ratio):
        ap.error("give --spacing and/or --ratio")
    a.sizes = a.sizes.split(",") if a.sizes else DEFAULT_SIZES[a.standard].split()
    if a.shape == "rect" and not (a.width and a.height):
        ap.error("rect needs --width and --height")
    if a.shape == "circle" and not a.diameter:
        ap.error("circle needs --diameter")
    cands = rect_candidates(a) if a.shape == "rect" else circle_candidates(a)
    if not cands:
        print("no candidate: the cover and tie leave no room for bars", file=sys.stderr)
        return 2
    if a.ratio:
        cands.sort(key=lambda d: abs(d["ratio"] - a.ratio))
    print(f"standard = {a.standard}, cover {a.cover:g} + tie {a.tie:g}"
          + (f", spacing <= {a.spacing:g}" if a.spacing else "") + (f", target ratio {a.ratio:g} %" if a.ratio else ""))
    for i, d in enumerate(cands, 1):
        warn = " (outside 1..8 %: check)" if not 1.0 <= d["ratio"] <= 8.0 else ""
        print(f"\n[{i}] {d['size']} (db {d['db']:.1f} mm), {d['n']} bars, ratio {d['ratio']:.2f} %{warn}\n"
              f"    bar centre at {d['c']:.1f} mm from the surface, {d['info']}, clear {d['clear']:.1f} mm"
              f"{clear_flag(d['clear'], d['db'])}")
        print("    rebar_standard = \"%s\"\n    bars = [" % a.standard)
        print("\n".join("    " + t for t in d["toml"]))
        print("    ]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
