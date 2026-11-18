"""Section preview before the PM run:  <PY> .claude/skills/biaxnm/preview.py out/NAME/NAME.toml

Draws the parsed section (concrete, steel, rebar, tendons) to NAME_section.png next to
the config file and prints the path. No analysis is run. Run from the repository root.
"""
from __future__ import annotations

import importlib.util
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if importlib.util.find_spec("biaxnm") is None:      # register ./src under the name "biaxnm" (same as src/main.py)
    _src = os.path.join(ROOT, "src")
    _spec = importlib.util.spec_from_file_location("biaxnm", os.path.join(_src, "__init__.py"),
                                                   submodule_search_locations=[_src])
    _mod = importlib.util.module_from_spec(_spec)
    sys.modules["biaxnm"] = _mod
    _spec.loader.exec_module(_mod)

from biaxnm.core.config import load_model  # noqa: E402
from biaxnm.output.plot import draw_section  # noqa: E402


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: preview.py CONFIG.toml", file=sys.stderr)
        return 2
    cfg = argv[0]
    try:
        model = load_model(cfg)
    except (ValueError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = model.section
    fig, ax = plt.subplots(figsize=(6, 6))
    draw_section(ax, s)
    ax.set_title(f"{os.path.basename(cfg)}  section")
    out = os.path.splitext(cfg)[0] + "_section.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
