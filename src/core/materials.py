"""Bind a design code's material laws to the elements of a section.

``ResolvedMaterials`` holds the per-element laws (code defaults filled in), the elements' coordinates as
arrays, the laws grouped by identity (so a law is evaluated once per group) and the tendon prestrains.
It does no integration.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from ..codes.base import DesignCode
from ..section import PolygonSection


def edge_arrays(rings: List[np.ndarray]):
    X1 = np.concatenate([r[:, 0] for r in rings])
    Y1 = np.concatenate([r[:, 1] for r in rings])
    X2 = np.concatenate([np.roll(r[:, 0], -1) for r in rings])
    Y2 = np.concatenate([np.roll(r[:, 1], -1) for r in rings])
    return X1, Y1, X2, Y2


def group_by_law(laws) -> List[Tuple[object, np.ndarray]]:
    groups: Dict[int, Tuple[object, list]] = {}
    for i, law in enumerate(laws):
        groups.setdefault(id(law), (law, []))[1].append(i)
    return [(law, np.array(idx)) for law, idx in groups.values()]


class ResolvedMaterials:
    def __init__(self, section: PolygonSection, code: DesignCode):
        self.section, self.code = section, code
        self.conc = code.concrete_law()
        self.conc_bps = self.conc.breakpoints()
        self.conc_edges = edge_arrays([section.outer] + section.holes)
        self.conc_pts = section.outer

        # structural steel regions
        self.steel = []
        default_steel = code.steel_law() if any(s.law is None for s in section.steel) else None
        for s in section.steel:
            law = s.law or default_steel
            self.steel.append((edge_arrays([s.polygon] + list(s.holes)), law, s.polygon))

        # bars
        self.bar_x = np.array([b.x for b in section.bars], dtype=float)
        self.bar_y = np.array([b.y for b in section.bars], dtype=float)
        self.bar_A = np.array([b.area for b in section.bars], dtype=float)
        default_rebar = code.rebar_law() if any(b.law is None for b in section.bars) else None
        bar_laws = [b.law or default_rebar for b in section.bars]
        self.bar_laws = bar_laws
        self.bar_groups = group_by_law(bar_laws)

        # tendons
        self.t_x = np.array([t.x for t in section.tendons], dtype=float)
        self.t_y = np.array([t.y for t in section.tendons], dtype=float)
        self.t_A = np.array([t.area for t in section.tendons], dtype=float)
        # area of concrete displaced by each tendon (the duct); defaults to the tendon area
        self.t_Ad = np.array([t.area if t.duct_area is None else t.duct_area for t in section.tendons], dtype=float)
        default_tendon = code.tendon_law() if any(t.law is None for t in section.tendons) else None
        t_laws = [t.law or default_tendon for t in section.tendons]
        self.t_laws = t_laws
        self.t_groups = group_by_law(t_laws)
        self.t_pe = np.zeros(len(section.tendons))
        for law, idx in self.t_groups:
            for i in idx:
                t = section.tendons[i]
                # the strain at which the tendon's own law gives the effective prestress
                self.t_pe[i] = t.eps_pe if t.eps_pe is not None else law.prestrain_for(t.fse)

        # extreme reinforcement positions for the pivot geometry
        pts = [np.column_stack([self.bar_x, self.bar_y])] if len(self.bar_x) else []
        pts += [np.column_stack([self.t_x, self.t_y])] if len(self.t_x) else []
        pts += [np.asarray(s.polygon) for s in section.steel]
        self.rein_pts = np.vstack(pts) if pts else None
