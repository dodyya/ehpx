# ehpx

Lambda algebra calculator plus a Curtis-algorithm driver for the EHP spectral
sequence. Part of *IML Spring 2026: Computing with the EHP sequence*.

## What it does

Arithmetic in the lambda algebra Λ over **F₂** — the free associative graded
algebra on generators λᵢ (written `[i]`) modulo the admissibility condition
and the Adem relations — together with the differential that makes it a DGA.

- **Admissible monomials** — sequences `[s₁, …, sₙ]` with `sₖ₊₁ ≤ 2·sₖ`
- **Multiplication** — reduces non-admissible products via Adem
- **Differential** — `d([i]) = Σ C(i−j, j) [i−j, j−1]` (mod 2), Leibniz-extended
- **Curtis algorithm** — fills the Curtis table by stem/filtration, computes
  differentials, and prints the surviving cycles in H\*(Λ)

## Build

```sh
cargo build --release
```

(Install Rust via `rustup` from https://rustup.rs if you don't have it.)

## Binaries

### `ehpx` — interactive REPL

```sh
cargo run --bin ehpx
```

```
> [3]
λ_(3)
> [2] * [2]
λ_(3, 1)
> d([4])
λ_(2, 1) + λ_(3, 0)
> [2] * [3] + [4, 1]
λ_(3, 1, 1) + λ_(4, 1)
> quit
```

Syntax: `[a,b,c]` admissible monomial · `+` F₂ addition · `*` multiplication
(Adem applied automatically) · `d(…)` differential · parens for grouping ·
`quit` / `exit` to leave.

### `table` — Curtis table

```sh
cargo run --release --bin table [max_stem]
```

Runs the Curtis algorithm through the given stem (default 12), prints the
dot chart (cycles / sources / targets per bidegree), the list of
differentials, a per-stem detail view, and the final survivor list.

### `check_diff` — correctness test bench

```sh
cargo run --release --bin check_diff [max_deg]
```

Exercises the Rust differential against the reference implementation in
`lambda.py` (found by walking up from CWD). Runs every admissible monomial
of degree ≤ `max_deg`, every 2-element F₂ sum, and capped samples of 3-
and 4-element sums (~2200 tests at `max_deg=8`). Requires `python3`;
exits nonzero on any mismatch.

## Layout

- `src/lib.rs` — algebra: `Admissible`, `Monomial`, `Element`, mul, diff, Adem
- `src/curtis.rs` — Curtis algorithm & table rendering
- `src/repl.rs` — parser + interpreter for the REPL
- `src/bin/table.rs` — `table` binary entry point
- `src/bin/check_diff.rs` — correctness test bench
- `scripts/check_driver.py` — wire-format shim around `lambda.py`
- `lambda.py` — Python reference implementation of the differential
