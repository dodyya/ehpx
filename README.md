<img width="3500" height="1100" alt="curtis" src="https://github.com/user-attachments/assets/1e49f4cd-cfbf-4ed5-ae7d-be7c60dd7286" />

# ehpx

Lambda algebra calculator plus a Curtis-algorithm implementation for the EHP spectral
sequence. Part of *IML Spring 2026: Computing with the EHP sequence*.

## What?

We implement arithmetic and differentiation for the Lambda algebra (per Bousfield-Curtis-Kan-Quillen-Rector-Schlesinger). This is a differential filtered graded algebra over 𝔽₂, whose homology gives the E² page of the Adams spectral sequence. Essentially, it's an algebraic object whose structure tells us something about homotopy groups of spheres, which is a (the!!) central open problem in algebraic topology. 

Elements of the algebra are "admissible" finite sequences of natural numbers, ones where every element has to be no more than double the previous. Multiplication is concatenation of these sequences, except when the "seam" would make the result inadmissible, in which case the product decomposes into a sum of admissible sequences. Addition happens modulo 2, so every element can be thought of as a set of monomials. Differentiation is essentially left-multiplying by -1. 

Monomials admit a grading by taking the sum of all elements in a sequence. Accordingly, differentiation decreases this value by 1. This is how the rows of the Curtis table are arranged. Additionally, we can filter all monomials by their leading term, which is how the columns work. The "point" of the Curtis table calculation is to calculate the homology of this algebra. We can avoid computing and testing all admissible sequences (entirely intractable) and instead restrict to elements that could be representatives of homology classes, by looking at the Algebraic EHP spectral sequence. It is this spectral sequence that the differentials come from; in computing the Curtis table up to some n, we are computing the part of the sseq that converges to the homology of the n-th filtration of Lambda. 

## TL;DR
This is kind of just using dynamic programming, enabled by some mathematical considerations (AEHP) to efficiently compute a property (homology) of an esoteric object (Lambda-algebra) that tells us something about another object (stable homotopy groups of spheres), that is somewhat close to an object of interest (homotopy groups of spheres) for algebraic topologists (nerds).

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
list. 25 is about as far as I got on my machine.

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

### `check_diff` — correctness test bench; not for public use

## Layout

- `src/lib.rs` — algebra: `Admissible`, `Monomial`, `Element`, mul, diff, Adem
- `src/curtis.rs` — Curtis algorithm, table rendering, interactive debugger
- `src/repl.rs` — parser + interpreter for the REPL
- `src/bin/table.rs` — `table` binary (flags + file output)
- `src/bin/check_diff.rs` — correctness test bench
- `scripts/check_driver.py` — wire-format shim around `lambda.py` for ground truth
- `scripts/visualize_table.py` — matplotlib bidegree chart from JSON
- `lambda.py` — Python reference implementation of the differential

## References

Allen, Keita. "Computing the Homology of the C-Motivic Lambda Algebra." University of Chicago Mathematics REU, 2022. http://math.uchicago.edu/~may/REU2022/REUPapers/Allen.pdf.

Bousfield, A., E. Curtis, D. Kan, D. Quillen, D. L. Rector, and J. W. Schlesinger. "The mod-*p* lower central series and the Adams spectral sequence." *Topology* 5, no. 4 (1966): 331–342.

Ravenel, Douglas C. *Complex Cobordism and Stable Homotopy Groups of Spheres*. 2nd ed. AMS Chelsea Publishing, 2004. ISBN 978-0-8218-2967-7.

Tangora, Martin C. *Computing the Homology of the Lambda Algebra*. Memoirs of the American Mathematical Society, vol. 58, no. 337. AMS, 1985. ISBN 978-0-8218-2338-5.

