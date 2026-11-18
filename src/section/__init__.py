"""Cross-section description and geometric validation.

Units: mm, mm^2, MPa.  Coordinates are global (x, y); compression is positive.
"""
from .elements import Bar, SteelRegion, Tendon
from .polygon import PolygonSection
from .ringgeom import (_clean_ring, _inside, _on_boundary, _orient, _orient_sign, _rings_touch, _seg_hits,
                       _self_intersects, clip_halfplane, dedupe_ring, ring_properties, signed_area)
from .validation import _steel_pair_ok, validate_section

__all__ = ["Bar", "Tendon", "SteelRegion", "PolygonSection", "dedupe_ring", "signed_area", "ring_properties",
           "validate_section"]
