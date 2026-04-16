#!/usr/bin/env python3
"""
Lambda algebra differential — reference implementation over F₂.

Wire format (stdin → stdout):
  One element per line, written as space-separated monomials.
  Monomial: comma-separated non-negative ints, e.g.  "2,1"  for λ₂λ₁.
  Unit monomial (empty sequence): "1".
  Zero element:  "0"  (or empty line).

Output: the differential of each input element, in the same canonical format
        (monomials sorted ascending by lex order on their index sequences).
"""
import sys


# ── F₂ binomial ───────────────────────────────────────────────────────────────

def binom_mod2(n: int, k: int) -> bool:
    """C(n, k) mod 2 via Lucas' theorem: 1 iff (k & n) == k."""
    if k < 0 or k > n:
        return False
    return (n & k) == k


# ── Lambda algebra ────────────────────────────────────────────────────────────
#
# An *Element* is a frozenset of admissible monomial tuples.
# Unit = ()  (empty tuple).   F₂ add = symmetric_difference.

def _xor(fs: frozenset, m: tuple) -> frozenset:
    """Toggle monomial m inside a frozenset (F₂)."""
    return fs.symmetric_difference({m})


def adem_element(n: int, b: int) -> frozenset:
    """
    Adem relation:  λ_n · λ_b  for  b > 2n  (inadmissible junction).
    Returns a frozenset of 2-tuples  (l, r)  such that
        λ_n · λ_b = Σ λ_l · λ_r  (mod 2).
    Formula: i = b − 2n − 1;
      for j = 0 … ⌊(i−1)/2⌋:  if C(i−j−1, j) is odd, add (n+i−j,  2n+j+1).
    """
    i = b - 2 * n - 1       # ≥ 0 because b > 2n
    if i == 0:
        return frozenset()
    result: frozenset = frozenset()
    for j in range((i - 1) // 2 + 1):
        if binom_mod2(i - j - 1, j):
            result = _xor(result, (n + i - j, 2 * n + j + 1))
    return result


def mul_mono(s1: tuple, s2: tuple) -> frozenset:
    """
    Multiply two admissible monomials.  Applies Adem relations at the junction
    recursively until the product is admissible.  Returns a frozenset of
    admissible monomial tuples.
    """
    if not s1 or not s2:
        return frozenset({s1 + s2})
    if s2[0] <= 2 * s1[-1]:            # admissible junction → just concatenate
        return frozenset({s1 + s2})
    # Inadmissible junction: reduce via Adem
    z, v = s1[-1], s2[0]
    prefix, suffix = s1[:-1], s2[1:]
    middle = adem_element(z, v)         # frozenset of 2-tuples
    result: frozenset = frozenset()
    for mid in middle:                  # mid is an admissible 2-monomial
        pm = mul_mono(prefix, mid)      # prefix * mid  (may recurse)
        for m in pm:
            for t in mul_mono(m, suffix):
                result = _xor(result, t)
    return result


def mul_elem(a: frozenset, b: frozenset) -> frozenset:
    """Multiply two Elements."""
    result: frozenset = frozenset()
    for s1 in a:
        for s2 in b:
            for t in mul_mono(s1, s2):
                result = _xor(result, t)
    return result


# ── Differential ─────────────────────────────────────────────────────────────

def diff_gen(i: int) -> frozenset:
    """
    d(λ_i) = Σ_{j=1}^{⌊i/2⌋}  C(i−j, j) · λ_{i−j} · λ_{j−1}
    """
    result: frozenset = frozenset()
    for j in range(1, i // 2 + 1):
        if binom_mod2(i - j, j):
            for t in mul_mono((i - j,), (j - 1,)):
                result = _xor(result, t)
    return result


def diff_mono(seq: tuple) -> frozenset:
    """
    d(λ_{a₁} ··· λ_{aₙ})  via iterative Leibniz:
        Σ_k  λ_{a₁}···λ_{a_{k−1}}  ·  d(λ_{aₖ})  ·  λ_{a_{k+1}}···λ_{aₙ}
    """
    if not seq:
        return frozenset()
    if len(seq) == 1:
        return diff_gen(seq[0])
    result: frozenset = frozenset()
    for k in range(len(seq)):
        dg = diff_gen(seq[k])
        if not dg:
            continue
        prefix_elem = frozenset({seq[:k]})
        suffix_elem = frozenset({seq[k + 1:]})
        term = mul_elem(mul_elem(prefix_elem, dg), suffix_elem)
        result = result.symmetric_difference(term)
    return result


def diff_elem(elem: frozenset) -> frozenset:
    """d(Σ monomials) = Σ d(monomial)."""
    result: frozenset = frozenset()
    for seq in elem:
        result = result.symmetric_difference(diff_mono(seq))
    return result


# ── Wire format ───────────────────────────────────────────────────────────────

def parse_line(line: str) -> frozenset:
    line = line.strip()
    if not line or line == '0':
        return frozenset()
    result: frozenset = frozenset()
    for token in line.split():
        seq: tuple = () if token == '1' else tuple(int(x) for x in token.split(','))
        result = _xor(result, seq)
    return result


def format_elem(elem: frozenset) -> str:
    if not elem:
        return '0'
    parts = []
    for seq in sorted(elem):            # lex order on tuples matches Rust SmallVec lex
        parts.append('1' if not seq else ','.join(map(str, seq)))
    return ' '.join(parts)


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    for line in sys.stdin:
        elem = parse_line(line)
        print(format_elem(diff_elem(elem)), flush=True)


if __name__ == '__main__':
    main()
