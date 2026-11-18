"""Write the 8 comparison cases as BiaxNM TOML files (type = "ec2-05-2") (verification/cases/biaxnm/input/caseN_*.toml).

    python/python.exe verification/scripts/make_toml.py

Geometry and rebar rules are in verification/input_summary.md.  Bars are listed point by point: every closed bar path is
split into segments, each segment contributes its first point and its interior points (the end point belongs to the
next segment), so no corner bar is doubled.
"""
import math
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "verification", "cases", "biaxnm", "input")
R2 = math.sqrt(2.0)

HEADER = """# comparison case {no}: {title}
# Inputs are summarized in verification/input_summary.md
# Material: fck 30, fy 500, gamma_c 1.5, gamma_s 1.15, alpha_cc 0.85 (EN 1992-2, na = "UK")
[code]
type = "ec2-05-2"
fck = 30
fy = 500
na = "UK"

"""
TAIL = """
[analysis]
n_c = 200
reference = "geometric"

[loads.LC1]
Pu = 0
Mux = 100
Muy = 100
"""


def path_points(vertices, divs):
    """Closed polygon path: divs[i] equal segments on the edge vertices[i] -> vertices[i+1]."""
    pts = []
    for i, d in enumerate(divs):
        (x0, y0), (x1, y1) = vertices[i], vertices[(i + 1) % len(vertices)]
        for k in range(d):
            t = k / d
            pts.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))
    return pts


def rect_vertices(hx, hy):
    return [(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy)]


def track_points(r, straight, div_arc, div_line):
    h = straight / 2.0
    pts = [(r, -h + straight * k / div_line) for k in range(div_line)]
    pts += [(r * math.cos(math.pi * k / div_arc), h + r * math.sin(math.pi * k / div_arc)) for k in range(div_arc)]
    pts += [(-r, h - straight * k / div_line) for k in range(div_line)]
    pts += [(r * math.cos(math.pi + math.pi * k / div_arc), -h + r * math.sin(math.pi + math.pi * k / div_arc))
            for k in range(div_arc)]
    return pts


def hex_vertices(d):
    """Hexagon of the cases 7, 8 (flats +-250, apex +-600) with every edge moved inward by d."""
    yf = 250.0 - d
    xa = 600.0 - d * R2
    xt = xa - yf
    return [(xt, -yf), (xa, 0.0), (xt, yf), (-xt, yf), (-xa, 0.0), (-xt, -yf)]


def fmt_points(pts, dia):
    rows = ",\n".join(f"  {{ x = {x:.4f}, y = {y:.4f}, dia = {dia} }}" for x, y in pts)
    return "bars = [\n" + rows + ",\n]\n"


def fmt_polygon(vs):
    return "[" + ", ".join(f"[{x:.4f}, {y:.4f}]" for x, y in vs) + "]"


CASES = []


def add(no, name, title, section_text):
    CASES.append((no, name, title, section_text))


# 1 solid circle
add(1, "solid_circle", "solid circle D500, 16 D16, dc 50",
    '[section]\nouter = { shape = "circle", r = 250, n = 96 }\n'
    'bars = [ { circle = [0, 0, 200], n = 16, dia = 16 } ]\n')

# 2 hollow circle
add(2, "hollow_circle", "hollow circle D500 T100, 16+16 D12, dc 30",
    '[section]\nouter = { shape = "circle", r = 250, n = 96 }\n'
    'holes = [ { shape = "circle", r = 150, n = 96 } ]\n'
    'bars = [\n  { circle = [0, 0, 220], n = 16, dia = 12 },\n  { circle = [0, 0, 180], n = 16, dia = 12 },\n]\n')

# 3 solid rectangle
add(3, "solid_rect", "solid rectangle 800 x 300, 12 D16, dc 50",
    '[section]\nouter = { shape = "rect", x = -400, y = -150, width = 800, height = 300 }\n'
    + fmt_points(path_points(rect_vertices(350, 100), [3, 3, 3, 3]), 16))

# 4 hollow rectangle
add(4, "hollow_rect", "hollow rectangle 800 x 300, opening 600 x 100, 12+12 D12, dc 30",
    '[section]\nouter = { shape = "rect", x = -400, y = -150, width = 800, height = 300 }\n'
    'holes = [ { shape = "rect", x = -300, y = -50, width = 600, height = 100 } ]\n'
    + fmt_points(path_points(rect_vertices(370, 120), [3, 3, 3, 3])
                 + path_points(rect_vertices(330, 80), [3, 3, 3, 3]), 12))

# 5 solid track (width 500, straight 500 -> total length 1000)
add(5, "solid_track", "solid track B500 H500, 22 D20, dc 50",
    '[section]\nouter = { shape = "track", length = 1000, width = 500, rotation_deg = 90, n_arc = 48 }\n'
    + fmt_points(track_points(200, 500, 6, 5), 20))

# 6 hollow track
add(6, "hollow_track", "hollow track B500 H500 T100, 22+22 D12, dc 30",
    '[section]\nouter = { shape = "track", length = 1000, width = 500, rotation_deg = 90, n_arc = 48 }\n'
    'holes = [ { shape = "track", length = 800, width = 300, rotation_deg = 90, n_arc = 48 } ]\n'
    + fmt_points(track_points(220, 500, 6, 5) + track_points(180, 500, 6, 5), 12))

# 7 solid hexagon
add(7, "solid_hexagon", "solid hexagon B1 700 B2 250 H2 250, 26 D16, dc 50",
    '[section]\nouter = ' + fmt_polygon(hex_vertices(0)) + '\n'
    + fmt_points(path_points(hex_vertices(50), [3, 3, 7, 3, 3, 7]), 16))

# 8 hollow hexagon
add(8, "hollow_hexagon", "hollow hexagon T100, 26+26 D12, dc 30",
    '[section]\nouter = ' + fmt_polygon(hex_vertices(0)) + '\nholes = [ ' + fmt_polygon(hex_vertices(100)) + ' ]\n'
    + fmt_points(path_points(hex_vertices(30), [3, 3, 7, 3, 3, 7]) + path_points(hex_vertices(70), [3, 3, 7, 3, 3, 7]), 12))


def main():
    os.makedirs(OUT, exist_ok=True)
    for no, name, title, section in CASES:
        path = os.path.join(OUT, f"case{no}_{name}.toml")
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(HEADER.format(no=no, title=title) + section + TAIL)
        print(path)


if __name__ == "__main__":
    main()
