#!/usr/bin/env python3
"""
Plot the empirical complexity of the Curtis table: total entries per stem.

Usage:
    python3 plot_complexity.py state.json [output.png]

The naive analysis of the population step is exponential (each new stem
multiplies survivors by λ_{n-1}, and there are many of them per row).
The actual growth — read off the table itself — is much tamer; this
script shows by how much.

If `/tmp/classic_curtis.txt` (Bill's classic Curtis table) is present,
its entries-per-stem are overlaid for comparison.  Source:
    https://williamb.info/lambda/classic-curtis-table.txt
"""
from __future__ import annotations
import json
import re
import sys
from collections import Counter
from pathlib import Path


def count_per_stem_ours(state_path: Path) -> dict[int, int]:
    s = json.loads(state_path.read_text())
    entries = [e for e in s["entries"] if not e.get("artifact", False)]
    return Counter(e["stem"] for e in entries)


def count_per_stem_classic(txt_path: Path) -> dict[int, int]:
    """Parse Bill's table.  Each line `((stem filt) #(...) n #(...) n)` lists
    one or two admissible words; we count every #(...) token by total degree."""
    counts: Counter[int] = Counter()
    for line in txt_path.read_text().splitlines():
        if not line.lstrip().startswith("(("):
            continue
        for m in re.findall(r"#\(([^)]+)\)", line):
            counts[sum(int(d) for d in m.split())] += 1
    return counts


def main():
    argv = sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__.strip())
        sys.exit(0 if argv else 2)

    in_path = Path(argv[0])
    out_path = Path(argv[1]) if len(argv) >= 2 else Path("curtis_complexity.png")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ours = count_per_stem_ours(in_path)
    stems_ours = sorted(ours)
    counts_ours = [ours[k] for k in stems_ours]

    classic_path = Path("/tmp/classic_curtis.txt")
    classic = count_per_stem_classic(classic_path) if classic_path.exists() else {}
    if classic:
        stems_cl = sorted(classic)
        counts_cl = [classic[k] for k in stems_cl]

    fig, (ax_lin, ax_log) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # ── Linear scale ─────────────────────────────────────────────────────
    ax_lin.plot(stems_ours, counts_ours, "o-", color="#13294b",
                markersize=4, linewidth=1.4, label="ehpx (this repo)")
    if classic:
        ax_lin.plot(stems_cl, counts_cl, "s--", color="#ff5f05",
                    markersize=3, linewidth=1.0, alpha=0.8,
                    label="classic (Bill's table)")
    ax_lin.set_ylabel("entries at stem k")
    ax_lin.set_title(f"Curtis table complexity through stem {max(stems_ours)}")
    ax_lin.grid(alpha=0.3)
    ax_lin.legend(loc="upper left", fontsize=9)

    # ── Log scale, with reference growth rates ───────────────────────────
    # Plot the same data on a log y-axis and overlay reference curves
    # anchored at stem 20 so the reader can eyeball whether growth is
    # roughly linear, polynomial (which would be a straight line in
    # loglog), or exponential (straight in semilog-y).
    ax_log.plot(stems_ours, counts_ours, "o-", color="#13294b",
                markersize=4, linewidth=1.4, label="ehpx (this repo)")
    if classic:
        ax_log.plot(stems_cl, counts_cl, "s--", color="#ff5f05",
                    markersize=3, linewidth=1.0, alpha=0.8,
                    label="classic (Bill's table)")

    # Reference: anchor at stem 20.
    anchor_k = 20
    if anchor_k in ours:
        anchor = ours[anchor_k]
        ks = list(range(1, max(stems_ours) + 1))
        # Linear:  c * k
        lin = [anchor * (k / anchor_k) for k in ks]
        # Quadratic: c * k^2
        quad = [anchor * (k / anchor_k) ** 2 for k in ks]
        # Cubic: c * k^3
        cub = [anchor * (k / anchor_k) ** 3 for k in ks]
        # Exponential 2^k (would saturate y axis, so scaled to anchor)
        expo = [anchor * 2 ** (k - anchor_k) for k in ks]
        ax_log.plot(ks, lin,  ":", color="#888", linewidth=1.0, label="∝ k (linear)")
        ax_log.plot(ks, quad, ":", color="#555", linewidth=1.0, label="∝ k²")
        ax_log.plot(ks, cub,  ":", color="#222", linewidth=1.0, label="∝ k³")
        ax_log.plot(ks, expo, ":", color="#a00", linewidth=1.0, label="∝ 2ᵏ (exponential)")
        ax_log.set_ylim(0.5, max(counts_ours) * 4)

    ax_log.set_yscale("log")
    ax_log.set_xlabel("stem  k = t − s")
    ax_log.set_ylabel("entries at stem k (log)")
    ax_log.grid(alpha=0.3, which="both")
    ax_log.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    print(f"wrote {out_path}", file=sys.stderr)

    # Also dump a plain-text table for the curious.
    print(f"\nstem  ehpx{'  classic' if classic else ''}", file=sys.stderr)
    for k in stems_ours:
        line = f"{k:>4}  {ours[k]:>4}"
        if classic and k in classic:
            line += f"  {classic[k]:>4}"
        print(line, file=sys.stderr)


if __name__ == "__main__":
    main()
