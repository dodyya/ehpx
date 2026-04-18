#!/usr/bin/env python3
"""
Visualize a Curtis table report as a bidegree chart.

Usage:
    python3 visualize_table.py report.json [output.png]

Expects the JSON format emitted by `cargo run --bin table -- --json`.
Writes a PNG (or any format matplotlib supports, inferred from extension).

Layout: x-axis = stem (t-s), y-axis = filtration (n).  Each entry is a
dot; multiple entries in the same bidegree stack with a small x-offset
and are labelled with their admissible sequence.  Differentials are
vertical arrows within a column (they preserve stem and lower filtration).
"""

from __future__ import annotations
import json
import sys
from pathlib import Path
from collections import defaultdict


ROLE_COLOR = {
    "cycle":  "#2e7d32",   # green — survives, H*(Λ) candidate
    "source": "#c62828",   # red — has outgoing differential
    "target": "#f9a825",   # amber — hit by a differential
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
    entries = report["entries"]
    diffs = report["differentials"]

    # Bucket entries by bidegree so we can space them out within a cell.
    by_cell = defaultdict(list)   # (stem, row) -> [entry, ...]
    for e in entries:
        by_cell[(e["stem"], e["row"])].append(e)

    # Map each entry to a plot (x, y) position so arrows find the right dots.
    # Points sit INSIDE cells (not at grid intersections): stem k cell spans
    # x ∈ [k, k+1] with center x = k + 0.5; filtration row n cell spans
    # y ∈ [n-1, n] with center y = n - 0.5.  Entries in the same bidegree
    # cluster organically around the cell center via sunflower / golden-
    # angle phyllotaxis — the radius grows as √i so area density stays
    # uniform, and the total spread saturates well before the cell edge.
    import math

    GOLDEN = math.pi * (3 - math.sqrt(5))  # ≈ 137.5° in radians

    pos = {}  # (stem, row, tuple(seq)) -> (x, y)
    for (stem, row), bucket in by_cell.items():
        bucket.sort(key=lambda e: e["seq"])
        n = len(bucket)
        cx = stem + 0.5
        cy = row - 0.5
        if n == 1:
            for e in bucket:
                pos[(stem, row, tuple(e["seq"]))] = (cx, cy)
            continue
        # Spread saturates as the cell fills:
        #   n=2  → ~0.18    (just off-center)
        #   n=4  → ~0.25
        #   n=16 → ~0.40
        #   n→∞  → 0.50     (cell boundary)
        spread = 0.50 * (1.0 - 1.0 / (1.0 + 0.35 * n))
        for i, e in enumerate(bucket):
            # Offset by 0.5 so no point lands exactly at the center (which
            # would look out of place when other siblings are orbiting it).
            r = spread * math.sqrt((i + 0.5) / n)
            theta = i * GOLDEN
            x = cx + r * math.cos(theta)
            y = cy + r * math.sin(theta)
            pos[(stem, row, tuple(e["seq"]))] = (x, y)

    # Figure: scale width with max_stem, height with max row observed.
    max_row = max((r for _, r in (k[:2] for k in ((e["stem"], e["row"]) for e in entries))), default=1)
    fig_w = max(10, 1.0 * (max_stem + 1))
    fig_h = max(6, 0.55 * max_row)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Plot dots + labels.
    for e in entries:
        x, y = pos[(e["stem"], e["row"], tuple(e["seq"]))]
        color = ROLE_COLOR.get(e["role"], "#444")
        ax.plot(x, y, "o", color=color, markersize=7, zorder=3)
        label = format_seq(e["seq"])
        ax.annotate(
            label,
            (x, y),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=6,
            color="#333",
            zorder=4,
            alpha=0.85,
        )

    # Differentials: curved arrow from source → target.  In the lambda
    # algebra, d lowers total degree by 1 (source stem = sum(src_seq),
    # target stem = sum(tgt_seq) = src_stem - 1) and also lowers filtration.
    # So arrows are generally diagonal (left + down), not vertical.
    for i, d in enumerate(diffs):
        src_stem = sum(d["src"])
        tgt_stem = sum(d["tgt"])
        sk = (src_stem, d["src_row"], tuple(d["src"]))
        tk = (tgt_stem, d["tgt_row"], tuple(d["tgt"]))
        if sk not in pos or tk not in pos:
            continue
        x1, y1 = pos[sk]
        x2, y2 = pos[tk]
        rad = 0.18 if (i % 2 == 0) else -0.18
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="->,head_width=0.35,head_length=0.6",
                color="#c62828",
                alpha=0.7,
                lw=1.3,
                shrinkA=6,
                shrinkB=6,
                connectionstyle=f"arc3,rad={rad}",
            ),
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
