# ehpx

Lambda algebra calculator. Part of **IML Spring 2026: Computing with the EHP sequence**.

## What it does

Computes in the Lambda algebra over **F₂**: the free associative graded algebra on generators λ_i (written `[i]`) with the admissibility condition and Adem relations, together with the differential that makes it a DGA.

- **Admissible monomials** — sequences `[s₁, s₂, ..., sₙ]` with each sₖ₊₁ ≤ 2·sₖ
- **Multiplication** — reduces to admissible form via Adem relations
- **Differential** — `d([i]) = Σ C(i−j, j) [i−j, j−1]` (mod 2), extended by the Leibniz rule
- **REPL** — interactive expression evaluator
- **`gen` binary** — enumerates all admissible monomials by degree, writes to a file

## Getting Rust

If you don't have Rust installed:

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Follow the prompts (the defaults are fine). Then open a new terminal, or run:

```sh
source "$HOME/.cargo/env"
```

Verify it worked:

```sh
rustc --version   # should print something like: rustc 1.85.0 (...)
cargo --version
```

> **Windows**: download and run the installer from https://rustup.rs — it handles everything.
> **Linux**: same `curl` command above works; you may also need `build-essential` / `gcc` if prompted.

## Building and running

Clone the repo and enter the directory, then:

```sh
cargo build
```

### REPL

```sh
cargo run --bin ehpx
```

At the `>` prompt you can type expressions:

```
> [3]
λ_(3)
> [2] * [2]
λ_(3, 1)
> d([4])
λ_(2, 1) + λ_(3)
> d([2] * [3])
λ_(1, 2) + λ_(2, 1)
> [2] * [3] + [4, 1]
λ_(3, 1, 1) + λ_(4, 1)
> quit
```

**Syntax:**
- `[a, b, c]` — admissible monomial λ_a λ_b λ_c
- `+` — addition (mod 2, so `x + x = 0`)
- `*` — multiplication (applies Adem relations automatically)
- `d(expr)` — differential
- Parentheses for grouping
- `quit` or `exit` to leave

### Enumerate admissible monomials

```sh
cargo run --bin gen
```

Writes `admissible_monomials.txt` to the working directory, listing every admissible monomial by degree up to degree 20. To change the cutoff, edit `N_MAX` in `src/bin/gen.rs`.
