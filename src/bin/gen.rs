use std::fs::File;
use std::io::{BufWriter, Write};

const N_MAX: usize = 20;

fn main() {
    let path = "admissible_monomials.txt";
    let file = File::create(path).expect("cannot create output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "Admissible monomials up to degree {N_MAX}").unwrap();
    writeln!(w, "Condition: each term ≤ 2 × previous term, all terms ≥ 1").unwrap();
    writeln!(w, "Format: λ_(s1, s2, ..., sk)").unwrap();

    // degree 0: the unit (empty sequence)
    writeln!(w, "\ndeg 0  [1 monomial]").unwrap();
    writeln!(w, "  1").unwrap();

    for d in 1..=N_MAX {
        let mut all: Vec<Vec<usize>> = Vec::new();
        enumerate(d, usize::MAX, &mut Vec::new(), &mut all);
        // sort: by length ascending, then lex
        all.sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        writeln!(w, "\ndeg {d}  [{} monomials]", all.len()).unwrap();
        for seq in &all {
            let inner: Vec<String> = seq.iter().map(|x| x.to_string()).collect();
            writeln!(w, "  λ_({})", inner.join(", ")).unwrap();
        }
    }

    eprintln!("written to {path}");
}

/// Recursively enumerate all admissible sequences summing to `remaining`,
/// where the next appended element must be ≤ `max_next` (= 2 × last, or ∞ for first).
fn enumerate(remaining: usize, max_next: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if remaining == 0 {
        out.push(cur.clone());
        return;
    }
    let cap = remaining.min(max_next);
    for next in 1..=cap {
        cur.push(next);
        enumerate(remaining - next, 2 * next, cur, out);
        cur.pop();
    }
}
