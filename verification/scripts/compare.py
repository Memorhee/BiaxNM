"""Compare BiaxNM (Eurocode 2) with concreteproperties 0.8.0.

    python/python.exe verification/scripts/compare.py [rc [case numbers]] [phc] [src] [plot]      (no argument: everything)

rc    8 reinforced concrete sections (ec2-05-2, na "UK"), load direction Mx = My (45 deg), N from pure tension to pure compression.
phc   hollow circular PHC pile with 24 prestressing bars (ec2-04-1-1, na "EN", fck 80).  PrestressedSection only bends about x
      and must be symmetric about y, so for a load direction phi the section is rotated by -phi and the positive uniaxial
      capacity is taken.  Also an independent 1 mm grid integration with the EC2 strain planes (pivot B / pivot C).
src   SRC column (encased H section + rebar, the shipped example), three load directions (0-360 deg), biaxial bending.
      ConcreteSection with meshed steel and lumped bars.
plot  cases/compare/rc.png (2 x 4 panels) and cases/compare/phc_src.png (1 x 2 panels) from the stored csv files.

Geometry (concrete polygons, bars, tendons, steel) is read from the BiaxNM result csv files, so both programs analyse the same
section.  Writes cases/biaxnm/input, cases/biaxnm/result, cases/concreteproperties/*.csv and cases/compare/*.png.
concreteproperties convention: compression positive, N and M in N, N.mm.
Moment reference: rc - the geometric centroid (origin); phc / src - the centroid of the concrete outline with holes removed
(BiaxNM reference = "geometric"), given to concreteproperties as moment_centroid.
"""
import csv
import glob
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402
from shapely import affinity  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402
from shapely.ops import unary_union  # noqa: E402
from sectionproperties.pre import CompoundGeometry, Geometry  # noqa: E402

import concreteproperties.stress_strain_profile as ssp  # noqa: E402
from concreteproperties import (Concrete, ConcreteSection, PrestressedSection, Steel, SteelBar,  # noqa: E402
                                SteelStrand)
from concreteproperties.pre import add_bar  # noqa: E402


def _fast_get_stress(self, strain):
    """np.interp version of StressStrainProfile.get_stress (the original builds a scipy interp1d on every call, ~55% of run time).
    Same result: piecewise linear with linear extrapolation beyond both ends."""
    x, y = self.strains, self.stresses
    v = np.interp(strain, x, y)
    if strain < x[0] and x[1] > x[0]:
        v = y[0] + (strain - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
    elif strain > x[-1] and x[-1] > x[-2]:
        v = y[-1] + (strain - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return v


ssp.StressStrainProfile.get_stress = _fast_get_stress

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CASES = os.path.join(ROOT, "verification", "cases")
INPUT = os.path.join(CASES, "biaxnm", "input")
RESULT = os.path.join(CASES, "biaxnm", "result")
CPOUT = os.path.join(CASES, "concreteproperties")
IMG = os.path.join(CASES, "compare")
plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

GC, GS, EP, ES = 1.5, 1.15, 200000.0, 200000.0        # EN: gamma_c, gamma_s, Ep, Es


# ------------------------------------------------------------------ shared helpers
def make_concrete(fcd, e2, eu, n, E=30000.0, n_points=12):
    return Concrete(name="C", density=2.4e-6, stress_strain_profile=ssp.ConcreteLinear(elastic_modulus=E),
                    ultimate_stress_strain_profile=ssp.EurocodeParabolicUltimate(
                        compressive_strength=fcd, compressive_strain=e2, ultimate_strain=eu, n=n, n_points=n_points),
                    flexural_tensile_strength=0.0, colour="lightgrey")


def ec2_concrete(fck):
    """fcd (alpha_cc 1.0), eps_c2, eps_cu2, n of EN 1992-1-1 Table 3.1."""
    fcd = fck / GC
    if fck <= 50:
        return fcd, 0.002, 0.0035, 2.0
    return (fcd, (2.0 + 0.085 * (fck - 50.0) ** 0.53) / 1000.0, (2.6 + 35.0 * ((90.0 - fck) / 100.0) ** 4) / 1000.0,
            1.4 + 23.4 * ((90.0 - fck) / 100.0) ** 4)


def run_biaxnm(toml, out):
    os.makedirs(out, exist_ok=True)
    subprocess.run([sys.executable, os.path.join(ROOT, "src", "main.py"), toml, "--out", out], check=True, capture_output=True)


def read_rings(path):
    rings = {}
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rings.setdefault(r["ring"], []).append((float(r["x_mm"]), float(r["y_mm"])))
    return rings


def read_curve(path, pure_compression=False):
    """BiaxNM curve csv -> sorted (Pn kN, Mn kN.m); ``pure_compression`` keeps the c = inf row (P0, 0)."""
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["c_mm"] == "inf" and not pure_compression:
                continue
            rows.append((float(r["Pn_kN"]), float(r["Mn_kNm"])))
    rows.sort()
    return rows


def interp(pairs, n):
    for (n0, m0), (n1, m1) in zip(pairs, pairs[1:]):
        if n0 <= n <= n1 and n1 > n0:
            return m0 + (m1 - m0) * (n - n0) / (n1 - n0)
    return None


def stats(rows, label):
    d = np.abs([r[-1] for r in rows])
    print(f"{label}: {len(d)} points  mean |diff| {d.mean():.3f}%  median {np.median(d):.3f}%  max {d.max():.3f}%")


def cp_capacity(cs, n, phi):
    """Capacity along direction phi (radians, atan2(My, Mx)) at axial force n (N): scan the neutral-axis angle for a sign change of
    the moment direction, refine with brentq and reject wrap-around jumps.  Returns |M| (kN.m) or None."""
    cache = {}

    def f(th):
        if th not in cache:
            r = cs.ultimate_bending_capacity(theta=float(th), n=n)
            cache[th] = (r.m_x, r.m_y)
        return cache[th]

    def h(th):
        mx, my = f(th)
        d = math.atan2(my, mx) - phi
        return (d + math.pi) % (2 * math.pi) - math.pi

    ths = np.linspace(-np.pi, np.pi, 49)
    vals = []
    for th in ths:
        try:
            vals.append(h(th))
        except Exception:
            vals.append(None)
    roots = []
    for i in range(len(ths) - 1):
        a, b = vals[i], vals[i + 1]
        if a is None or b is None or a * b > 0:
            continue
        th0 = brentq(h, ths[i], ths[i + 1], xtol=1e-9)
        if abs(h(th0)) > 1e-4:          # a wrap-around jump, not a root
            continue
        mx, my = f(th0)
        roots.append(math.hypot(mx, my) / 1e6)
    return max(roots) if roots else None


# ================================================================== rc: 8 reinforced concrete sections
NAMES = {1: "Solid Circle", 2: "Hollow Circle", 3: "Solid Rectangle", 4: "Hollow Rectangle",
         5: "Solid Track", 6: "Hollow Track", 7: "Solid Hexagon", 8: "Hollow Hexagon"}
FCD, FYD, EC = 17.0, 500.0 / 1.15, 32837.0          # ec2-05-2, na "UK": 0.85 x 30 / 1.5, 500 / 1.15
EPS_C2, EPS_CU2 = 0.002, 0.0035


def rc_csv(no):
    return os.path.join(CPOUT, f"case{no}_{NAMES[no].replace(' ', '_')}.csv")


def rc_geometry(no):
    sec = glob.glob(os.path.join(RESULT, f"case{no}", f"case{no}_*_LC1_section.csv"))[0]
    rings = read_rings(sec)
    holes = [v for k, v in rings.items() if k.startswith("hole")]
    bars = []
    with open(sec.replace("_section", "_bars"), newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            bars.append((float(r["x_mm"]), float(r["y_mm"]), float(r["area_mm2"])))
    return Polygon(rings["outer"], holes), bars


def rc_section(no):
    outer, bars = rc_geometry(no)
    concrete = make_concrete(FCD, EPS_C2, EPS_CU2, 2.0, E=EC, n_points=10)
    steel = SteelBar(name="B500", density=7.85e-6, colour="grey",
                     stress_strain_profile=ssp.SteelElasticPlastic(yield_strength=FYD, elastic_modulus=ES, fracture_strain=1.0))
    geom = Geometry(geom=outer, material=concrete)
    for x, y, a in bars:
        geom = add_bar(geometry=geom, area=a, material=steel, x=x, y=y, n=8)
    return ConcreteSection(geom, moment_centroid=(0.0, 0.0))


def at_45(cs, n, guess=None):
    """Capacity (Mx, My) at axial force n (N) with Mx = My (45 deg load direction). Returns (theta, mx, my) or None.

    ``guess``: neutral-axis angle of a neighbouring N; the root is searched near it first (warm start), then over the full circle.
    """
    cache = {}

    def f(th):
        if th not in cache:
            r = cs.ultimate_bending_capacity(theta=float(th), n=n)
            cache[th] = (r.m_x, r.m_y)
        return cache[th]

    def h(th):
        mx, my = f(th)
        return np.arctan2(my, mx) - np.pi / 4

    def search(ths):
        prev = None
        for th in ths:
            try:
                cur = (th, h(th))
            except Exception:
                prev = None
                continue
            if prev and abs(prev[1]) < np.pi / 2 and abs(cur[1]) < np.pi / 2 and prev[1] * cur[1] <= 0:
                th0 = brentq(h, prev[0], th, xtol=1e-6)
                return th0, f(th0)[0], f(th0)[1]
            prev = cur
        return None

    if guess is not None:
        for half, k in ((0.3, 5), (0.8, 9)):
            r = search(np.linspace(guess - half, guess + half, k))
            if r:
                return r
    return search(np.linspace(-np.pi, np.pi, 73))


def rc_case(no):
    cs = rc_section(no)
    gp = cs.gross_properties
    # kappa0 = uniform strain eps_cu2 over the whole section (concreteproperties' pure compression point)
    p0 = cs.moment_interaction_diagram(theta=0, limits=[("kappa0", 0.0), ("d_n", 1e-6)], progress_bar=False).results[0].n
    pt = -gp.reinf_lumped_area * FYD     # pure tension: all bars at fyd, concrete cracked
    print(f"\n=== case {no} {NAMES[no]}  Ac {gp.concrete_area:.0f}  As {gp.reinf_lumped_area:.1f}  "
          f"P0(kappa0) {p0 / 1e3:.1f} kN  tension {pt / 1e3:.1f} kN")
    ns = np.concatenate([np.linspace(pt * 0.98, 0.0, 9), np.linspace(0.0, p0 * 0.995, 18)[1:]])
    pts, found, guess = [], {}, None
    for n in sorted([n for n in ns if n >= 0]) + sorted([n for n in ns if n < 0], reverse=True):
        if n < 0 and (n == max(m for m in ns if m < 0)):
            guess = found.get(0.0, guess)
        r = at_45(cs, n, guess)
        if r is None:
            continue
        th, mx, my = r
        guess = th
        found[float(n)] = th
        pts.append((n / 1e3, np.hypot(mx, my) / 1e6, mx / 1e6, my / 1e6))
        print(f"  N {n / 1e3:9.1f} kN  Mx=My {mx / 1e6:8.2f}  Mn {np.hypot(mx, my) / 1e6:8.2f} kN.m", flush=True)
    pts.sort()
    os.makedirs(CPOUT, exist_ok=True)
    with open(rc_csv(no), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N_kN", "Mn_kNm", "Mnx_kNm", "Mny_kNm"])
        w.writerows(pts)
    bx = read_curve(glob.glob(os.path.join(RESULT, f"case{no}", f"case{no}_*_LC1.csv"))[0])
    errs, mref = [], max(m for _, m in bx)
    for n, m, _, _ in pts:
        mb = interp(bx, n)
        if mb is not None:
            errs.append(abs(mb - m))
    m0_cp = at_45(cs, 0.0)
    m0_cp = np.hypot(m0_cp[1], m0_cp[2]) / 1e6 if m0_cp else float("nan")
    return dict(no=no, name=NAMES[no], mean=100 * np.mean(errs) / mref, mx=100 * np.max(errs) / mref, p0_cp=p0 / 1e3,
                p0_b=max(n for n, _ in bx), pt_cp=pt / 1e3, pt_b=min(n for n, _ in bx), m0_cp=m0_cp, m0_b=interp(bx, 0.0))


def run_rc(nos):
    results = [rc_case(n) for n in nos]
    print("\n| # | 단면 | 곡선 오차 평균 / 최대 | P0 CP / BiaxNM (kN) | 순인장 CP / BiaxNM (kN) | Md(N=0) CP / BiaxNM (kN·m) |")
    print("|---|---|---|---|---|---|")
    for r in results:
        print(f"| {r['no']} | {r['name']} | {r['mean']:.2f}% / {r['mx']:.2f}% | {r['p0_cp']:.0f} / {r['p0_b']:.0f} "
              f"({100 * (r['p0_b'] / r['p0_cp'] - 1):+.1f}%) | {r['pt_cp']:.0f} / {r['pt_b']:.0f} "
              f"({100 * (r['pt_b'] / r['pt_cp'] - 1):+.1f}%) | {r['m0_cp']:.1f} / {r['m0_b']:.1f} "
              f"({100 * (r['m0_b'] / r['m0_cp'] - 1):+.1f}%) |")


# ================================================================== phc: PHC pile (prestressing bars)
PHC_FCK, PHC_FPY, PHC_FSE = 80.0, 1275.0, 700.0
PHC_DIRS = [0.0, 7.5, 11.0, 22.5, 33.3, 45.0]
PHC_NS = [-1500, -1000, -500, 0, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]      # kN
PHC_TOML = """# 중공 원형 PHC 말뚝: 프리스트레스 강봉 24개, Eurocode 2 (EN 1992-1-1:2004), 비교용.
[code]
type = "ec2-04-1-1"
na = "EN"
fck = 80
fy = 400
fpy = 1275
Ep = 200000

[section]
outer = { shape = "circle", cx = 0, cy = 0, r = 300, n = 96 }
holes = [ { shape = "circle", cx = 0, cy = 0, r = 170, n = 96 } ]
tendons = [ { circle = [0, 0, 235], n = 24, area = 71, fse = 700 } ]

"""


def phc_section(rot_deg, d):
    rings = read_rings(os.path.join(d, "phc_pile_P1_section.csv"))
    holes = [v for k, v in rings.items() if k.startswith("hole")]
    poly = affinity.rotate(Polygon(rings["outer"], holes), rot_deg, origin=(0, 0))
    a = math.radians(rot_deg)
    ten = []
    with open(os.path.join(d, "phc_pile_P1_tendons.csv"), newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            x, y = float(r["x_mm"]), float(r["y_mm"])
            ten.append((x * math.cos(a) - y * math.sin(a), x * math.sin(a) + y * math.cos(a), float(r["area_mm2"])))
    fpd = PHC_FPY / GS
    strand = SteelStrand(name="PC", density=7.85e-6, colour="grey", prestress_stress=PHC_FSE,
                         stress_strain_profile=ssp.StrandHardening(yield_strength=fpd, elastic_modulus=EP, fracture_strain=1.0,
                                                                   breaking_strength=fpd))
    g = Geometry(geom=poly, material=make_concrete(*ec2_concrete(PHC_FCK)))
    for x, y, a_ in ten:
        g = add_bar(geometry=g, area=a_, material=strand, x=x, y=y, n=8)
    return PrestressedSection(g, moment_centroid=(0.0, 0.0))


def phc_job(args):
    phi, d = args
    cs = phc_section(-phi, d)
    out = []
    for n in PHC_NS:
        r = cs.ultimate_bending_capacity(n=n * 1e3)
        out.append((n, math.hypot(r.m_x, r.m_y) / 1e6))
    return phi, out


def phc_pivot_check():
    """Independent 1 mm grid integration of the PHC section with the EC2 strain planes (pivot B for c <= h, pivot C above)."""
    R, r, h = 300.0, 170.0, 600.0
    g = 1.0
    xs = np.arange(-R + g / 2, R, g)
    X, Y = np.meshgrid(xs, xs)
    rho = np.hypot(X, Y)
    y = Y[(rho <= R) & (rho >= r)]
    ty = 235 * np.sin(np.arange(24) * 2 * np.pi / 24)
    fcd, e2, eu, nx = ec2_concrete(PHC_FCK)
    fpd, At = PHC_FPY / GS, 71.0

    def conc(e):
        x = np.clip(e / e2, 0, 1)
        return np.where(e > 0, fcd * (1 - (1 - x) ** nx), 0.0)

    def forces(c):
        if c <= h:
            eps = lambda yy: eu * (c - (R - yy)) / c
        else:
            dp = (1 - e2 / eu) * h
            eps = lambda yy: e2 * (c - (R - yy)) / (c - dp)
        sc = conc(eps(y))
        sp = np.clip(EP * (eps(ty) - PHC_FSE / EP), -fpd, fpd)
        Ft = At * (sp - conc(eps(ty)))
        return ((sc * g * g).sum() + Ft.sum()) / 1e3, ((sc * g * g * y).sum() + (Ft * ty).sum()) / 1e6

    out = {}
    for nt in (6000, 7000, 8000, 9000):
        f = lambda c: forces(c)[0] - nt
        cs = np.geomspace(60, 20000, 400)
        v = [f(c) for c in cs]
        c0 = next(brentq(f, cs[i], cs[i + 1]) for i in range(len(cs) - 1) if v[i] * v[i + 1] < 0)
        out[nt] = float(forces(c0)[1])
    return out


def run_phc():
    d = os.path.join(RESULT, "phc")
    toml = os.path.join(INPUT, "phc_pile.toml")
    s = PHC_TOML
    for k, p in enumerate(PHC_DIRS, 1):
        s += f"[loads.P{k}]\nPu = 3000\nMux = {100 * math.cos(math.radians(p)):.6f}\nMuy = {100 * math.sin(math.radians(p)):.6f}\n\n"
    open(toml, "w", encoding="utf-8").write(s)
    run_biaxnm(toml, d)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=min(len(PHC_DIRS), os.cpu_count() or 2)) as ex:
        res = list(ex.map(phc_job, [(p, d) for p in PHC_DIRS]))
    print(f"[phc] concreteproperties {time.time() - t0:.0f}s")
    rows = []
    for k, (phi, out) in enumerate(res, 1):
        bx = read_curve(os.path.join(d, f"phc_pile_P{k}.csv"))
        for n, m in out:
            mb = interp(bx, n)
            if mb is not None:
                rows.append((phi, n, mb, m, 100 * (mb / m - 1)))
    os.makedirs(CPOUT, exist_ok=True)
    with open(os.path.join(CPOUT, "phc.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["dir_deg", "N_kN", "BiaxNM_kNm", "CP_kNm", "diff_pct"])
        w.writerows(rows)
    p0 = max(n for n, _ in read_curve(os.path.join(d, "phc_pile_P1.csv"), pure_compression=True))
    stats(rows, "[phc] all")
    stats([r for r in rows if r[1] <= 0.5 * p0], f"[phc] N <= 0.5 P0 ({0.5 * p0:.0f} kN)")
    chk = phc_pivot_check()
    print("[phc] pivot-C grid check (kN.m):", {n: round(m, 2) for n, m in chk.items()})
    with open(os.path.join(CPOUT, "phc_pivot.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N_kN", "Mn_kNm"])
        w.writerows(chk.items())


# ================================================================== src: SRC columns (encased steel, biaxial bending)
SRC_HEAD = """# {title}, Eurocode 2 (EN 1992-1-1:2004), 비교용.
[code]
type = "ec2-04-1-1"
na = "EN"
fck = {fck}
fy = {fy}
Fy = {Fy}

"""
SRC = {
    "src": dict(title="600 x 600 SRC 기둥 (예제): 매입 H형강 + 철근 8개", fck=30, fy=400, Fy=355, dirs=[32.9, 104.7, 205.3], body="""
[section]
rebar_standard = "KS"
outer = { shape = "rect", width = 600, height = 600 }
bars = [
  { line = [[60, 60],  [540, 60]],  n = 3, size = "D25" },
  { line = [[60, 540], [540, 540]], n = 3, size = "D25" },
  { line = [[60, 300], [540, 300]], n = 2, size = "D25" },
]
[[section.steel]]
polygon = { shape = "h_section", d = 300, bf = 200, tw = 10, tf = 16, cx = 300, cy = 300, r = 13 }
"""),
}


def src_section(name, d):
    spec = SRC[name]
    rebar = SteelBar(name="R", density=7.85e-6, colour="grey",
                     stress_strain_profile=ssp.SteelElasticPlastic(yield_strength=spec["fy"] / GS, elastic_modulus=ES,
                                                                   fracture_strain=1.0))
    steel = Steel(name="S", density=7.85e-6, colour="blue",
                  stress_strain_profile=ssp.SteelElasticPlastic(yield_strength=spec["Fy"], elastic_modulus=ES, fracture_strain=1.0))
    stem = os.path.join(d, f"{name}_D1")
    rings = read_rings(stem + "_section.csv")
    holes = [v for k, v in rings.items() if k.startswith("hole")]
    outer = Polygon(rings["outer"], holes)
    srings = read_rings(stem + "_steel.csv")
    spolys = []
    for k, v in srings.items():
        if "hole" not in k:
            spolys.append(Polygon(v, [w for kk, w in srings.items() if kk.startswith(k) and "hole" in kk]))
    cp = outer.difference(unary_union(spolys))
    parts = list(cp.geoms) if hasattr(cp, "geoms") else [cp]        # a box interior is a separate concrete piece
    conc = make_concrete(*ec2_concrete(spec["fck"]))
    g = CompoundGeometry([Geometry(geom=p, material=conc) for p in parts] + [Geometry(geom=p, material=steel) for p in spolys])
    with open(stem + "_bars.csv", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            g = add_bar(geometry=g, area=float(r["area_mm2"]), material=rebar, x=float(r["x_mm"]), y=float(r["y_mm"]), n=8)
    ref = outer.centroid        # BiaxNM reference "geometric": centroid of the concrete outline, holes removed
    return ConcreteSection(g, moment_centroid=(ref.x, ref.y))


def src_job(args):
    name, k, d, ns, deg = args
    cs = src_section(name, d)
    return name, k, deg, [(n, cp_capacity(cs, n * 1e3, math.radians(deg))) for n in ns]


def run_src():
    jobs, curves = [], {}
    for name, spec in SRC.items():
        toml = os.path.join(INPUT, f"{name}.toml")
        s = SRC_HEAD.format(**spec) + spec["body"] + "\n"
        for k, deg in enumerate(spec["dirs"], 1):
            s += f"[loads.D{k}]\nPu = 1000\nMux = {100 * math.cos(math.radians(deg)):.6f}\nMuy = {100 * math.sin(math.radians(deg)):.6f}\n\n"
        open(toml, "w", encoding="utf-8").write(s)
        d = os.path.join(RESULT, name)
        run_biaxnm(toml, d)
        for k, deg in enumerate(spec["dirs"], 1):
            bx = read_curve(os.path.join(d, f"{name}_D{k}.csv"))
            curves[(name, k)] = bx
            pt, p0 = bx[0][0], bx[-1][0]
            ns = sorted(set(round(x) for x in [0.0] + list(np.linspace(0.85 * pt, 0.8 * p0, 8))))
            jobs.append((name, k, d, ns, deg))
    t0 = time.time()
    with ProcessPoolExecutor() as ex:
        res = list(ex.map(src_job, jobs))
    print(f"[src] concreteproperties {time.time() - t0:.0f}s")
    rows = []
    for name, k, deg, out in res:
        for n, m in out:
            mb = interp(curves[(name, k)], n)
            if m is None or mb is None:
                print(f"  {name} dir {deg} N {n}: no value (BiaxNM {mb}, CP {m})")
                continue
            rows.append((name, deg, n, mb, m, 100 * (mb / m - 1)))
    os.makedirs(CPOUT, exist_ok=True)
    with open(os.path.join(CPOUT, "src.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["case", "dir_deg", "N_kN", "BiaxNM_kNm", "CP_kNm", "diff_pct"])
        w.writerows(rows)
    stats(rows, "[src] all")
    for name in SRC:
        stats([r for r in rows if r[0] == name], f"[src] {name}")


# ================================================================== plot: 2 x 4 panels from the stored csv files
def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def style(ax, title):
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Mn (kN·m)", fontsize=8)
    ax.set_ylabel("Pn (kN)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)


def grid(path, draw, n_used, handles, rows, cols, size):
    fig, axs = plt.subplots(rows, cols, figsize=size)
    axs = np.atleast_1d(axs).ravel()
    for i in range(n_used):
        draw(i, axs[i])
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=9, frameon=False)
    os.makedirs(IMG, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_all():
    from matplotlib.lines import Line2D
    h_b = Line2D([], [], color="#1f77b4", lw=1.5, label="BiaxNM")
    h_c = Line2D([], [], color="r", marker=".", ls="", ms=7, label="concreteproperties")
    h_g = Line2D([], [], color="g", marker="^", ls="", ms=6, label="격자 적분 (pivot C)")

    def draw_rc(i, ax):
        no = i + 1
        bx = read_curve(glob.glob(os.path.join(RESULT, f"case{no}", f"case{no}_*_LC1.csv"))[0], pure_compression=True)
        cp = read_csv(rc_csv(no))
        ax.plot([m for _, m in bx], [n for n, _ in bx], "-", color="#1f77b4", lw=1.4)
        ax.plot([float(r["Mn_kNm"]) for r in cp], [float(r["N_kN"]) for r in cp], "r.", ms=4)
        style(ax, f"{no}. {NAMES[no]}")

    grid(os.path.join(IMG, "rc.png"), draw_rc, 8, [h_b, h_c], 2, 4, (16, 7.4))

    phc_rows = read_csv(os.path.join(CPOUT, "phc.csv"))
    piv = read_csv(os.path.join(CPOUT, "phc_pivot.csv"))
    src_rows = read_csv(os.path.join(CPOUT, "src.csv"))
    colors = ["#1f77b4", "#2ca02c", "#9467bd"]

    def draw_ps(i, ax):
        if i == 0:
            bx = read_curve(os.path.join(RESULT, "phc", "phc_pile_P1.csv"), pure_compression=True)
            ax.plot([m for _, m in bx], [n for n, _ in bx], "-", color="#1f77b4", lw=1.4)
            ax.plot([float(r["CP_kNm"]) for r in phc_rows], [float(r["N_kN"]) for r in phc_rows], "r.", ms=4)
            ax.plot([float(r["Mn_kNm"]) for r in piv], [float(r["N_kN"]) for r in piv], "g^", ms=5)
            style(ax, "PHC 말뚝")
            return
        for k, deg in enumerate(SRC["src"]["dirs"], 1):
            bx = read_curve(os.path.join(RESULT, "src", f"src_D{k}.csv"), pure_compression=True)
            ax.plot([m for _, m in bx], [n for n, _ in bx], "-", color=colors[k - 1], lw=1.3, label=f"BiaxNM {deg:g}°")
            pts = [r for r in src_rows if float(r["dir_deg"]) == deg]
            ax.plot([float(r["CP_kNm"]) for r in pts], [float(r["N_kN"]) for r in pts], "r.", ms=4)
        ax.legend(fontsize=7, loc="upper right")
        style(ax, "SRC 기둥")

    grid(os.path.join(IMG, "phc_src.png"), draw_ps, 2, [h_b, h_c, h_g], 1, 2, (9, 4.6))


if __name__ == "__main__":
    args = sys.argv[1:]
    nos = [int(a) for a in args if a.isdigit()]
    stages = [a for a in args if a in ("rc", "phc", "src", "plot")]
    if not stages:
        stages = ["rc", "plot"] if nos else ["rc", "phc", "src", "plot"]
    if "rc" in stages:
        run_rc(nos or list(NAMES))
    if "phc" in stages:
        run_phc()
    if "src" in stages:
        run_src()
    if "plot" in stages:
        plot_all()
