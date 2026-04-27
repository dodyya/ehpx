#!/usr/bin/env python3
"""
Parse the `table` binary's per-degree timing log into JSON, then plot
elapsed time vs degree on linear and log axes.

Usage:
    python3 plot_timing.py timings.json [output.png]
    python3 plot_timing.py --ingest run.stderr.log timings.json

The first form reads an already-parsed `timings.json`
    {"timings": [{"degree": 0, "seconds": 0.00}, ...],
     "command": "table 40 --json ..."}
and renders a plot.

The second form scrapes lines like
    "  degree 12: 0.00s  → /path/state_12.json (12007 bytes)"
out of the table binary's stderr (also matches the older "  stem 12: ..."
emit so historical logs still parse) and writes a freshly-built JSON.
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path


# Match either "degree NN: T.TTs" (current) or "stem NN: T.TTs" (legacy).
TIMING_RE = re.compile(
    r"^\s*(?:degree|stem)\s+(\d+):\s+(\d+(?:\.\d+)?)s",
    re.IGNORECASE,
)


def ingest(stderr_path: Path) -> dict:
    text = stderr_path.read_text()
    timings = []
    for line in text.splitlines():
        m = TIMING_RE.match(line)
        if m:
            timings.append({"degree": int(m.group(1)),
                            "seconds": float(m.group(2))})
    timings.sort(key=lambda t: t["degree"])
    return {"timings": timings, "source": str(stderr_path)}


def plot(timings_path: Path, out_path: Path):
    data = json.loads(timings_path.read_text())
    timings = data["timings"]
    if not timings:
        print(f"no timings in {timings_path}", file=sys.stderr)
        sys.exit(1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    degrees = [t["degree"] for t in timings]
    seconds = [t["seconds"] for t in timings]

    fig, (ax_lin, ax_log) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # ── Linear ───────────────────────────────────────────────────────────
    ax_lin.bar(degrees, seconds, color="#13294b", alpha=0.85)
    ax_lin.set_ylabel("step time (seconds)")
    max_deg = max(degrees)
    ax_lin.set_title(f"Per-degree compute time, table 0…{max_deg}")
    ax_lin.grid(axis="y", alpha=0.3)

    # Annotate the largest bar with its exact time so the eye doesn't
    # have to interpolate against the grid.
    for d, s in zip(degrees, seconds):
        if s >= 0.5 * max(seconds):
            ax_lin.annotate(f"{s:.1f}s", xy=(d, s),
                            xytext=(0, 4), textcoords="offset points",
                            ha="center", fontsize=8, color="#333")

    # ── Log scale + reference growth rates ───────────────────────────────
    # The table binary prints timings rounded to 2 decimals (0.01s
    # resolution), so anything quicker shows up as 0.00s — those points
    # would render as a flat plateau at the bottom of the log axis
    # carrying no real information.  Drop them, and instead draw a
    # horizontal guide at the resolution floor so the reader can see
    # what's been clipped.
    TIMER_FLOOR = 0.01
    sig_degrees = [d for d, s in zip(degrees, seconds) if s >= TIMER_FLOOR]
    sig_seconds = [s for s in seconds if s >= TIMER_FLOOR]
    ax_log.plot(sig_degrees, sig_seconds,
                "o-", color="#13294b", markersize=4, linewidth=1.4,
                label="measured")
    ax_log.axhline(TIMER_FLOOR, color="#888", linestyle=":", linewidth=0.8,
                   label=f"timer resolution ({TIMER_FLOOR}s)")

    if sig_seconds:
        anchor_k = sig_degrees[-1]
        anchor_t = sig_seconds[-1]
        ks = [k for k in sig_degrees if k >= 1]
        # Polynomial references: c · k^p, normalized to pass through anchor.
        for p, color, label in [
            (3, "#888", "∝ k³"),
            (5, "#555", "∝ k⁵"),
            (8, "#222", "∝ k⁸"),
        ]:
            curve = [anchor_t * (k / anchor_k) ** p for k in ks]
            ax_log.plot(ks, curve, ":", color=color, linewidth=1.0, label=label)
        # Exponential reference: c · 2^k normalized to anchor.
        expo = [anchor_t * 2 ** (k - anchor_k) for k in ks]
        ax_log.plot(ks, expo, ":", color="#a00", linewidth=1.0, label="∝ 2ᵏ (exponential)")
        # Pin the y-range so the timer-floor guide and references are visible
        # but we don't drag the axis down to ridiculous extremes.
        ax_log.set_ylim(TIMER_FLOOR / 3, max(sig_seconds) * 5)

    ax_log.set_yscale("log")
    ax_log.set_xlabel("degree  k")
    ax_log.set_ylabel("step time (s, log)")
    ax_log.grid(alpha=0.3, which="both")
    ax_log.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    print(f"wrote {out_path}", file=sys.stderr)


def main():
    argv = sys.argv[1:]
    if argv and argv[0] == "--ingest":
        if len(argv) < 3:
            print("usage: plot_timing.py --ingest run.stderr.log timings.json", file=sys.stderr)
            sys.exit(2)
        result = ingest(Path(argv[1]))
        Path(argv[2]).write_text(json.dumps(result, indent=2))
        print(f"wrote {argv[2]} ({len(result['timings'])} entries)", file=sys.stderr)
        return

    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__.strip())
        sys.exit(0 if argv else 2)

    in_path = Path(argv[0])
    out_path = Path(argv[1]) if len(argv) >= 2 else in_path.with_suffix(".png")
    plot(in_path, out_path)


if __name__ == "__main__":
    main()
