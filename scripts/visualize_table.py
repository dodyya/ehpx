#!/usr/bin/env python3
"""
Visualize a Curtis table report as a bidegree chart.

Usage:
    python3 visualize_table.py report.json [output.png]

Expects the JSON format emitted by `cargo run --bin table -- --json`.
Writes a PNG (or any format matplotlib supports, inferred from extension).

Layout: x-axis = stem (t-s), y-axis = filtration (n).  Each entry is a
dot; multiple entries in the same bidegree stack vertically inside the
cell, ordered by admissible-sequence length with the shortest monomials
at the visual bottom of the cell (matches the Adams spectral sequence
convention).  Differentials are straight, headless gray lines.
"""

from __future__ import annotations
import json
import sys
from pathlib import Path
from collections import defaultdict


ROLE_COLOR = {
    "cycle":  "#ff5f05",   # orange — survives, H*(Λ) candidate
    "source": "#13294b",   # navy — has outgoing differential
    "target": "#707372",   # gray — hit by a differential
}


def format_seq(seq):
    if not seq:
        return "1"
    return "λ(" + ",".join(str(x) for x in seq) + ")"


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
        matplotlib.use("Agg")  # non-interactive
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("matplotlib not installed.  `pip install matplotlib`", file=sys.stderr)
        sys.exit(1)

    max_stem = report["max_stem"]
    # Full-state JSONs (so the `--from` resume flag can round-trip) include
    # λ_0-tail bookkeeping entries flagged with `artifact: true`.  They're
    # needed for resume, not for visualization — skip them here so the chart
    # stays clean.  Older (pre-artifact-field) JSONs just have no such flag
    # and pass through unchanged.
    entries = [e for e in report["entries"] if not e.get("artifact", False)]
    diffs = [d for d in report["differentials"] if not d.get("artifact", False)]

    # Bucket entries by bidegree so we can space them out within a cell.
    by_cell = defaultdict(list)   # (stem, row) -> [entry, ...]
    for e in entries:
        by_cell[(e["stem"], e["row"])].append(e)

    # Map each entry to a plot (x, y) position so arrows find the right dots.
    # Points sit INSIDE cells (not at grid intersections): stem k cell spans
    # x ∈ [k, k+1] with center x = k + 0.5; filtration row n cell spans
    # y ∈ [n-1, n] with center y = n - 0.5.
    #
    # Arrangement within a cell: linear vertical stack, ordered by monomial
    # length so the shortest sequences sit at the visual bottom of the cell
    # (matches the Adams SS convention).  Recall the y-axis is inverted
    # (`set_ylim(max_row, 0)`), so "visual bottom" corresponds to a *larger*
    # data y inside the cell's [row-1, row] interval.

    pos = {}  # (stem, row, tuple(seq)) -> (x, y)
    CELL_MARGIN = 0.15  # padding from the cell edges
    for (stem, row), bucket in by_cell.items():
        # Primary key: sequence length ascending (shortest = index 0 = bottom).
        # Secondary: lex on the sequence itself, for deterministic ordering
        # among same-length monomials.
        bucket.sort(key=lambda e: (len(e["seq"]), e["seq"]))
        n = len(bucket)
        cx = stem + 0.5
        cy = row - 0.5
        if n == 1:
            pos[(stem, row, tuple(bucket[0]["seq"]))] = (cx, cy)
            continue
        # Map index 0 → visual bottom = data y = row - margin.
        # Map index n-1 → visual top = data y = row - 1 + margin.
        y_bot = row - CELL_MARGIN          # larger data y, visually lower
        y_top = (row - 1) + CELL_MARGIN    # smaller data y, visually higher
        for i, e in enumerate(bucket):
            t = i / (n - 1)  # 0 (shortest) .. 1 (longest)
            y = y_bot + t * (y_top - y_bot)
            pos[(stem, row, tuple(e["seq"]))] = (cx, y)

    # Figure: scale width with max_stem, height with max row observed.
    max_row = max((r for _, r in (k[:2] for k in ((e["stem"], e["row"]) for e in entries))), default=1)
    fig_w = max(10, 1.0 * (max_stem + 1))
    fig_h = max(6, 0.55 * max_row)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Plot dots + labels.  Labels sit just to the right of each dot, vertically
    # centered against it — with vertical stacking, this keeps each label on
    # its own row (instead of climbing diagonally and colliding with neighbours).
    for e in entries:
        x, y = pos[(e["stem"], e["row"], tuple(e["seq"]))]
        color = ROLE_COLOR.get(e["role"], "#444")
        ax.plot(x, y, "o", color=color, markersize=5, zorder=3)
        label = format_seq(e["seq"])
        ax.annotate(
            label,
            (x, y),
            xytext=(5, 0),
            textcoords="offset points",
            fontsize=5.5,
            color="#333",
            zorder=4,
            alpha=0.9,
            va="center",
        )

    # Differentials: straight, headless lines from source → target.  The
    # density at high stems makes arrowheads + curves illegible; a thin
    # gray segment is enough to pair source with target by eye.
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
            alpha=0.55,
            linewidth=0.7,
            solid_capstyle="round",
            zorder=2,
        )

    # Axes + grid.  Labels sit at cell centers; the actual grid lines run
    # along the cell boundaries (minor ticks) so each (stem, row) bidegree
    # reads as a distinct box.
    ax.set_xlabel("stem  (t − s)")
    ax.set_ylabel("filtration  n")
    ax.set_title(f"Curtis table through stem {max_stem}")

    # Major ticks: labels at cell centers.
    ax.set_xticks([k + 0.5 for k in range(max_stem + 1)])
    ax.set_xticklabels([str(k) for k in range(max_stem + 1)])
    ax.set_yticks([r - 0.5 for r in range(1, max_row + 1)])
    ax.set_yticklabels([str(r) for r in range(1, max_row + 1)])
    ax.tick_params(which="major", length=0)  # hide major tick marks

    # Minor ticks: cell-boundary gridlines.
    ax.set_xticks(range(0, max_stem + 2), minor=True)
    ax.set_yticks(range(0, max_row + 2), minor=True)
    ax.tick_params(which="minor", length=0)

    ax.set_xlim(0, max_stem + 1)
    # CS chart convention: filtration increases downward (row 1 on top).
    ax.set_ylim(max_row, 0)
    ax.grid(which="minor", alpha=0.35, linewidth=0.6)
    ax.grid(which="major", alpha=0)
    ax.set_axisbelow(True)

    # Legend.
    handles = [
        mpatches.Patch(color=ROLE_COLOR["cycle"], label="cycle (survivor)"),
        mpatches.Patch(color=ROLE_COLOR["source"], label="source of differential"),
        mpatches.Patch(color=ROLE_COLOR["target"], label="target of differential"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
