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
    # x ∈ [k, k+1] with center x = k + 0.5.
    #
    # Vertical layout: rows are *not* uniform height.  Each row's height is
    # set by the densest cell on that row — a row whose busiest cell has 18
    # entries is given enough y-space that those 18 dots + labels don't pile
    # on top of each other.  Sparse rows stay compact.  This is much more
    # space-efficient than expanding every row globally to fit the worst case.
    #
    # Within a cell, entries stack linearly along y, sorted by monomial
    # length so the shortest sequences sit at the visual bottom of the cell
    # (matches the Adams SS convention).  Recall the y-axis is inverted
    # (`set_ylim(total_h, 0)`), so "visual bottom" = *larger* data y.

    CELL_MARGIN = 0.15      # padding from cell edges (data-y units)
    PER_ENTRY = 0.28        # extra row height per entry beyond the first
    MIN_ROW_HEIGHT = 1.0    # baseline row height for sparse rows

    # How tall does each row need to be?  Driven by the densest cell on it.
    row_max_count = defaultdict(int)
    for (s, r), bucket in by_cell.items():
        if len(bucket) > row_max_count[r]:
            row_max_count[r] = len(bucket)
    max_row = max((r for _, r in by_cell.keys()), default=1)

    def row_height(r):
        n = row_max_count.get(r, 0)
        needed = 2 * CELL_MARGIN + max(0, n - 1) * PER_ENTRY
        return max(MIN_ROW_HEIGHT, needed)

    # Lay rows out top-to-bottom with cumulative offsets.  `row_top[r]` is
    # the smaller data-y boundary of row r's cell band; `row_bot[r]` is the
    # larger.  Row 1 starts at y=0.
    row_top: dict[int, float] = {}
    row_bot: dict[int, float] = {}
    y_cursor = 0.0
    for r in range(1, max_row + 1):
        row_top[r] = y_cursor
        y_cursor += row_height(r)
        row_bot[r] = y_cursor
    total_h = y_cursor

    # Dots sit near the LEFT edge of each cell (not the center) so the
    # right-of-dot label has the full remaining cell width to extend into
    # before bumping the next column's content.
    DOT_X_OFFSET = 0.10  # data-x units from left edge of the cell

    pos = {}  # (stem, row, tuple(seq)) -> (x, y)
    for (stem, row), bucket in by_cell.items():
        # Primary key: sequence length ascending (shortest = index 0 = bottom).
        # Secondary: lex on the sequence itself, for deterministic ordering
        # among same-length monomials.
        bucket.sort(key=lambda e: (len(e["seq"]), e["seq"]))
        n = len(bucket)
        cx = stem + DOT_X_OFFSET
        cy_top = row_top[row]
        cy_bot = row_bot[row]
        if n == 1:
            pos[(stem, row, tuple(bucket[0]["seq"]))] = (cx, (cy_top + cy_bot) / 2)
            continue
        # Map index 0 → visual bottom (larger data y); index n-1 → visual top.
        y_bot = cy_bot - CELL_MARGIN
        y_top = cy_top + CELL_MARGIN
        for i, e in enumerate(bucket):
            t = i / (n - 1)  # 0 (shortest) .. 1 (longest)
            y = y_bot + t * (y_top - y_bot)
            pos[(stem, row, tuple(e["seq"]))] = (cx, y)

    # Figure: width scales with stems; height scales with the *total* y
    # span (sum of variable row heights) so dense rows render at a usable
    # vertical resolution without bloating sparse ones.
    fig_w = max(10, 1.0 * (max_stem + 1))
    fig_h = max(6, 0.55 * total_h)
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

    # Major ticks: row labels at the center of each variable-height row.
    ax.set_xticks([k + 0.5 for k in range(max_stem + 1)])
    ax.set_xticklabels([str(k) for k in range(max_stem + 1)])
    ax.set_yticks([(row_top[r] + row_bot[r]) / 2 for r in range(1, max_row + 1)])
    ax.set_yticklabels([str(r) for r in range(1, max_row + 1)])
    ax.tick_params(which="major", length=0)  # hide major tick marks

    # Minor ticks: cell-boundary gridlines.  Y boundaries follow the
    # variable row layout; X boundaries stay at integer stem values.
    ax.set_xticks(range(0, max_stem + 2), minor=True)
    y_boundaries = [row_top[1]] + [row_bot[r] for r in range(1, max_row + 1)]
    ax.set_yticks(y_boundaries, minor=True)
    ax.tick_params(which="minor", length=0)

    ax.set_xlim(0, max_stem + 1)
    # CS chart convention: filtration increases downward (row 1 on top).
    ax.set_ylim(total_h, 0)
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
    fig.savefig(out_path, dpi=120)
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
