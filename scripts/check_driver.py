#!/usr/bin/env python3
"""
Wire-format shim that wraps lambda.py for use by check_diff.

argv[1] = absolute path to lambda.py

stdin:  one element per line — space-separated monomials.
        Monomial: comma-separated non-negative ints, '1' for the unit.
        Zero element: '0' or empty line.
stdout: the differential of each input element, same format.
"""
import sys
import os
import importlib.util

lambda_path = sys.argv[1]

# Load lambda.py by path (can't use `import lambda` — keyword clash).
spec = importlib.util.spec_from_file_location("_lambda", lambda_path)
mod  = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
diff_monomial = mod.diff_monomial   # diff_monomial(tuple) -> set of tuples


# ── wire format ───────────────────────────────────────────────────────────────

def parse_line(line: str) -> set:
    line = line.strip()
    if not line or line == '0':
        return set()
    result: set = set()
    for token in line.split():
        seq = () if token == '1' else tuple(int(x) for x in token.split(','))
        result ^= {seq}
    return result


def format_elem(elem: set) -> str:
    if not elem:
        return '0'
    parts = []
    for seq in sorted(elem):
        parts.append('1' if not seq else ','.join(map(str, seq)))
    return ' '.join(parts)


def diff_element(seqs: set) -> set:
    """d(Σ monomials) = Σ d(monomial)  (F₂)."""
    result: set = set()
    for seq in seqs:
        result ^= diff_monomial(seq)
    return result


# ── main ──────────────────────────────────────────────────────────────────────

for line in sys.stdin:
    elem = parse_line(line)
    print(format_elem(diff_element(elem)), flush=True)
