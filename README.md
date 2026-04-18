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
λ_(2, 2)
> d([4])
λ_(2, 1) + λ_(3, 0)
> [1,2,4,8,16,32]*[3]*[9]*[27]+d([55,22,11,12]+[33])
λ_(17, 15) + λ_(25, 7) + λ_(29, 3) + λ_(31, 1) + λ_(31, 23, 22, 11, 12) + λ_(47, 11, 18, 11, 12) + λ_(47, 13, 16, 11, 12) + λ_(47, 14, 15, 11, 12) + λ_(55, 11, 10, 11, 12) + λ_(55, 13, 8, 11, 12) + λ_(55, 14, 7, 11, 12) + λ_(55, 15, 6, 11, 12) + λ_(55, 19, 6, 7, 12) + λ_(55, 19, 8, 6, 11) + λ_(55, 21, 6, 6, 11) + λ_(55, 21, 8, 6, 9) + λ_(55, 21, 8, 7, 8) + λ_(55, 21, 8, 8, 7) + λ_(55, 21, 10, 6, 7) + λ_(55, 21, 10, 8, 5) + λ_(55, 21, 10, 9, 4) + λ_(55, 21, 10, 10, 3) + λ_(55, 22, 7, 6, 9) + λ_(55, 22, 7, 7, 8) + λ_(55, 22, 7, 8, 7) + λ_(55, 22, 11, 6, 5) + λ_(55, 22, 11, 7, 4) + λ_(55, 22, 11, 10, 1) + λ_(55, 22, 11, 11, 0) + λ_(1, 2, 4, 8, 16, 32, 9, 11, 19) + λ_(1, 2, 4, 8, 16, 32, 13, 11, 15)
> quit
```

Syntax: `[a,b,c]` admissible monomial · `+` F₂ addition · `*` multiplication
(Adem applied automatically) · `d(…)` differential · parens for grouping ·
`quit` / `exit` to leave.

### `table` — Curtis table

```sh
cargo run --release --bin table [MAX_STEM] [options]
```

Runs the Curtis algorithm through the given stem (default 12) and prints
the list of differentials, a per-stem detail view, and the final survivor
list.

Output style auto-detects: ANSI color for an interactive terminal, plain
text when piped or writing to a file.  Override with flags:

- `--plain` — no ANSI; Unicode (λ, →).  Good for `.txt` files in modern editors.
- `--ascii` — no ANSI; strict ASCII (`l(...)`, `->`).  Maximum portability.
- `--color` — force ANSI color even when output is redirected.
- `--json` — emit a machine-readable report (consumed by the visualizer).
- `-o FILE` — write to `FILE` instead of stdout.

Examples:

```sh
cargo run --release --bin table 12 --ascii -o report.txt
cargo run --release --bin table 12 --json  -o report.json
python3 scripts/visualize_table.py report.json curtis.png
```

### `visualize_table.py` — bidegree chart

```sh
python3 scripts/visualize_table.py report.json [output.png]
```

Reads the JSON emitted by `table --json` and renders the Curtis table as
a bidegree chart: stem on the x-axis, filtration on the y-axis; entries
are colored by role (cycle / source / target); differentials are drawn
as red arrows from source to target.  Requires `matplotlib`.

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
- `src/curtis.rs` — Curtis algorithm, table rendering, interactive debugger
- `src/repl.rs` — parser + interpreter for the REPL
- `src/bin/table.rs` — `table` binary (flags + file output)
- `src/bin/check_diff.rs` — correctness test bench
- `scripts/check_driver.py` — wire-format shim around `lambda.py`
- `scripts/visualize_table.py` — matplotlib bidegree chart from JSON
- `lambda.py` — Python reference implementation of the differential

## References

Douglas C. Ravenel, *Complex Cobordism and Stable Homotopy Groups of Spheres*, 2nd ed., AMS Chelsea Publishing, 2004. ISBN 978-0-8218-2967-7.

Martin C. Tangora, *Computing the Homology of the Lambda Algebra*, Memoirs of the American Mathematical Society, vol. 58, no. 337, AMS, 1985. ISBN 978-0-8218-2338-5.

Keita Allen, "Computing the Homology of the C-Motivic Lambda Algebra," University of Chicago Mathematics REU, 2022. http://math.uchicago.edu/~may/REU2022/REUPapers/Allen.pdf
