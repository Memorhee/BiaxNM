"""Command line entry point:  python/python.exe src/main.py examples/rect_aci318-25.toml [--out DIR] [--no-save]

Runs one config file (.toml or .json): prints the DCR of every load (or the curve summary when the file
has no loads) and writes PNG and CSV files.
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import warnings

if importlib.util.find_spec("biaxnm") is None:      # local Python: register ./src under the name "biaxnm"
    _src = os.path.dirname(os.path.abspath(__file__))
    _spec = importlib.util.spec_from_file_location("biaxnm", os.path.join(_src, "__init__.py"),
                                                   submodule_search_locations=[_src])
    _mod = importlib.util.module_from_spec(_spec)
    sys.modules["biaxnm"] = _mod
    _spec.loader.exec_module(_mod)

from biaxnm.core.check import CheckResult
from biaxnm.core.config import load_model, run_model


def model_summary(model) -> str:
    """Parsed model as text (kN, kN*m, mm), to compare a config file with what was intended."""
    s = model.section
    x0, y0, x1, y1 = s.bbox
    bar_area = sum(b.area for b in s.bars)
    lines = [f"code: {type(model.code).__name__}",
             f"concrete area: {s.area:,.0f} mm^2  bbox: x {x0:g}..{x1:g}, y {y0:g}..{y1:g} mm",
             f"bars: {len(s.bars)}, total {bar_area:,.1f} mm^2, rho = {100 * bar_area / s.area:.2f} %",
             f"steel regions: {len(s.steel)} / tendons: {len(s.tendons)}"]
    if model.loads:
        lines += [f"load {ld.name}: Pu = {ld.Pu / 1e3:g} kN, Mux = {ld.Mux / 1e6:g} kN*m, Muy = {ld.Muy / 1e6:g} kN*m"
                  for ld in model.loads]
    else:
        lines.append(f"curve only: lambda_u = {model.lambda_u_deg:g} deg")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="biaxnm", description="PM interaction diagram / DCR from a config file")
    ap.add_argument("config", help="model file (.toml or .json)")
    ap.add_argument("--out", help="output directory (default: [output] dir of the file)")
    ap.add_argument("--no-save", action="store_true", help="print only; write no PNG/CSV")
    ap.add_argument("--summary", action="store_true", help="print the parsed model and exit (no analysis)")
    args = ap.parse_args(argv)
    warnings.simplefilter("ignore")       # the command line reports errors only
    try:
        model = load_model(args.config)
        if args.summary:
            print(model_summary(model))
            return 0
        results = run_model(model, args.out, save=not args.no_save)
    except (ValueError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    for ld, res in zip(model.loads or [None], results):
        c = res.curve if isinstance(res, CheckResult) else res
        head = f"{ld.name}: DCR = {res.dcr:.3f} ({'OK' if res.ok else 'NG'})  " if isinstance(res, CheckResult) else ""
        print(f"{head}lambda_u = {c.lambda_u_deg:.2f} deg  P0 = {c.P0 / 1e3:.0f} kN  points = {len(c.c)}")
    if not args.no_save:
        print(f"saved to {args.out or model.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
