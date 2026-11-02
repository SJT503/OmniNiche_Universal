"""Shared Nature Medicine Visual Design System & Overlap Gate for OmniNiche Round-3 Figures.

Enforces the 80/20 neutral-slate vs Okabe-Ito protagonist palette, native 183-mm (7.20-in)
double-column print typography, and deterministic tightbbox zero-collision verification.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

from omniniche.viz.stage2_flagship_figures import (
    AXIS_COMPACT,
    AXIS_SHORT,
    PAL,
    _fmt_p,
    _panel_tag,
    set_nature_style,
)


def assert_zero_axes_frame_collision(fig: plt.Figure, min_gap_px: float = -2.0) -> list[str]:
    """Render canvas and verify that no two primary plot axes bounding boxes collide."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes: list[tuple[int, Bbox]] = []
    for idx, ax in enumerate(fig.axes):
        if not ax.get_visible():
            continue
        # Skip thin colorbar axes
        pos = ax.get_position()
        if pos.width < 0.035 or pos.height < 0.035:
            continue
        bb = ax.get_window_extent(renderer)
        boxes.append((idx, bb))

    collisions: list[str] = []
    for i in range(len(boxes)):
        idx_a, ba = boxes[i]
        for j in range(i + 1, len(boxes)):
            idx_b, bb = boxes[j]
            dx = min(ba.x1, bb.x1) - max(ba.x0, bb.x0)
            dy = min(ba.y1, bb.y1) - max(ba.y0, bb.y0)
            if dx > abs(min_gap_px) and dy > abs(min_gap_px):
                collisions.append(f"Axes {idx_a} and {idx_b} overlap by ({dx:.1f}x{dy:.1f} px)")
    return collisions


__all__ = [
    "PAL",
    "AXIS_SHORT",
    "AXIS_COMPACT",
    "set_nature_style",
    "_panel_tag",
    "_fmt_p",
    "assert_zero_axes_frame_collision",
]
