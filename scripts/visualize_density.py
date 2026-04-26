#!/usr/bin/env python3
"""
Visualize a Curtis table report as a *density* chart: a uniform bidegree
grid where each entry is a small colored dot and labels are dropped.

Usage:
    python3 visualize_density.py report.json [output.png]

Expects the JSON format emitted by `cargo run --bin table -- --json`.
Writes a PNG (or any format matplotlib supports, inferred from extension).

Layout: x-axis = stem (t-s), y-axis = filtration (n).  Each (stem, row)
cell is a uniform 1×1 box (no horizontal or vertical stretching) — the
goal is for *density and connections* to be the visual story, not
class identity.  Within a cell, entries are arranged in a sunflower /
golden-angle spiral so dots cluster naturally without overlapping.

For the labeled / class-identity view, see `visualize_classes.py`.
"""

from __future__ import annotations
import json
import math
import sys
from pathlib import Path
from collections import defaultdict


ROLE_COLOR = {
    "cycle":  "#ff5f05",   # orange — survives, H*(Λ) candidate
    "source": "#13294b",   # navy — has outgoing differential
    "target": "#707372",   # gray — hit by a differential
}


def main():
    argv = sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__.strip())
        sys.exit(0 if argv else 2)

    in_path = Path(argv[0])
    out_path = Path(argv[1]) if len(argv) >= 2 else in_path.with_suffix(".png")

    report = json.loads(in_path.read_text())

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("matplotlib not installed.  `pip install matplotlib`", file=sys.stderr)
        sys.exit(1)

    max_stem = report["max_stem"]
    # Drop the artifact (λ_0-tail bookkeeping) entries / diffs — they're
    # only there for resume round-tripping, not for visualization.
    entries = [e for e in report["entries"] if not e.get("artifact", False)]
    diffs   = [d for d in report["differentials"] if not d.get("artifact", False)]

    # Bucket entries by bidegree so we can spread them within a cell.
    by_cell = defaultdict(list)   # (stem, row) -> [entry, ...]
    for e in entries:
        by_cell[(e["stem"], e["row"])].append(e)

    max_row = max((r for _, r in by_cell.keys()), default=1)

    # Sunflower / phyllotaxis layout within each cell.  Cells are uniform
    # 1×1 boxes; cell (stem, row) has center (stem + 0.5, row - 0.5).
    # Radius grows as √i (constant area density), angle increments by the
    # golden angle so dots distribute without obvious lattice artefacts.
    GOLDEN = math.pi * (3 - math.sqrt(5))   # ≈ 137.5°

    pos = {}  # (stem, row, tuple(seq)) -> (x, y)
    for (stem, row), bucket in by_cell.items():
        # Deterministic in-cell ordering (lex on sequence).
        bucket.sort(key=lambda e: e["seq"])
        n = len(bucket)
        cx = stem + 0.5
        cy = row - 0.5
        if n == 1:
            pos[(stem, row, tuple(bucket[0]["seq"]))] = (cx, cy)
            continue
        # Spread saturates as the cell fills:
        #   n=2  → ~0.18 (just off-center)
        #   n=4  → ~0.25
        #   n=16 → ~0.40
        #   n→∞  → 0.50  (cell boundary)
        spread = 0.50 * (1.0 - 1.0 / (1.0 + 0.35 * n))
        for i, e in enumerate(bucket):
            # +0.5 keeps the i=0 dot off the exact center when siblings orbit.
            r = spread * math.sqrt((i + 0.5) / n)
            theta = i * GOLDEN
            x = cx + r * math.cos(theta)
            y = cy + r * math.sin(theta)
            pos[(stem, row, tuple(e["seq"]))] = (x, y)

    # Figure: uniform cell size, so width and height are simple linear
    # functions of stem count and max row.  Bigger inches-per-cell here
    # than in `visualize_classes.py` doesn't help — cells stay 1×1, so
    # we just want enough resolution to distinguish individual dots.
    INCHES_PER_CELL = 0.35
    fig_w = max(8, INCHES_PER_CELL * (max_stem + 1))
    fig_h = max(5, INCHES_PER_CELL * max_row)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Plot dots — small enough to make density legible at high stems
    # where cells get crowded, opaque enough to read individually in
    # sparse cells.
    for e in entries:
        x, y = pos[(e["stem"], e["row"], tuple(e["seq"]))]
        color = ROLE_COLOR.get(e["role"], "#444")
        ax.plot(x, y, "o", color=color, markersize=3.2, alpha=0.92, zorder=3)

    # Differentials: straight, headless gray lines.  Same as the labeled
    # viz — the eye should follow source ↔ target without arrowhead noise.
    for d in diffs:
        src_stem = sum(d["src"])
        tgt_stem = sum(d["tgt"])
        sk = (src_stem, d["src_row"], tuple(d["src"]))
        tk = (tgt_stem, d["tgt_row"], tuple(d["tgt"]))
        if sk not in pos or tk not in pos:
            continue
        x1, y1 = pos[sk]
        x2, y2 = pos[tk]
        ax.plot(
            [x1, x2], [y1, y2],
            color="#9A9A9A",
            alpha=0.5,
            linewidth=0.55,
            solid_capstyle="round",
            zorder=2,
        )

    # Axes + grid.  Major ticks label the cell centers; minor ticks at
    # integer boundaries draw the cell edges.
    ax.set_xlabel("stem  (t − s)")
    ax.set_ylabel("filtration  n")
    ax.set_title(f"Curtis table density through stem {max_stem}")

    ax.set_xticks([k + 0.5 for k in range(max_stem + 1)])
    ax.set_xticklabels([str(k) for k in range(max_stem + 1)], fontsize=7)
    ax.set_yticks([r - 0.5 for r in range(1, max_row + 1)])
    ax.set_yticklabels([str(r) for r in range(1, max_row + 1)], fontsize=7)
    ax.tick_params(which="major", length=0)

    ax.set_xticks(range(0, max_stem + 2), minor=True)
    ax.set_yticks(range(0, max_row + 2), minor=True)
    ax.tick_params(which="minor", length=0)

    ax.set_xlim(0, max_stem + 1)
    # CS chart convention: filtration increases downward (row 1 on top).
    ax.set_ylim(max_row, 0)
    ax.grid(which="minor", alpha=0.3, linewidth=0.5)
    ax.grid(which="major", alpha=0)
    ax.set_axisbelow(True)

    # Legend.
    handles = [
        mpatches.Patch(color=ROLE_COLOR["cycle"],  label="cycle (survivor)"),
        mpatches.Patch(color=ROLE_COLOR["source"], label="source of differential"),
        mpatches.Patch(color=ROLE_COLOR["target"], label="target of differential"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
