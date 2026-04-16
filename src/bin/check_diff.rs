//! Differential correctness checker — compares every Rust answer against
//! the lambda.py reference implementation.
//!
//! Usage:
//!   cargo run --bin check_diff [max_deg]
//!
//! `max_deg` (default 8) bounds the total degree of monomials used in test
//! cases.  The binary finds `lambda.py` by walking up the directory tree from
//! the current working directory (so it works from the worktree or the repo
//! root).  `python3` must be available on PATH.
//!
//! Test suite:
//!   • Every admissible monomial of total degree 0 … max_deg       (singletons)
//!   • A handful of λ₀-containing monomials                         (singletons)
//!   • All 2-element F₂ sums from monomials of degree ≤ combo2_deg
//!   • All 3-element F₂ sums (capped at 500)  — degree ≤ combo3_deg
//!   • All 4-element F₂ sums (capped at 200)  — degree ≤ combo4_deg
//!
//! Exit code 0 = all pass, 1 = at least one mismatch, 2 = setup error.

use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};

use ehpx::{Admissible, Element, Monomial, Seq};

// Thin wire-format shim around lambda.py — embedded so the binary is
// self-contained apart from needing lambda.py on disk.
const DRIVER: &str = include_str!("../../scripts/check_driver.py");

// ── helpers ───────────────────────────────────────────────────────────────────

fn make_mono(seq: &[usize]) -> Monomial {
    Monomial {
        seq: Admissible(Seq::from_slice(seq)),
        deg: seq.iter().sum(),
    }
}

fn singleton(seq: &[usize]) -> Element {
    Element::singleton(make_mono(seq))
}

fn format_mono(m: &Monomial) -> String {
    if m.seq.0.is_empty() {
        "1".to_string()
    } else {
        m.seq.0.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
    }
}

/// Canonical wire-format for an Element: space-separated monomials in
/// ascending lex order (Element invariant), or "0" for the zero element.
fn format_elem(e: &Element) -> String {
    if e.0.is_empty() {
        return "0".to_string();
    }
    e.0.iter().map(format_mono).collect::<Vec<_>>().join(" ")
}

/// Walk up the directory tree from `start` looking for `lambda.py`.
fn find_lambda_py(start: &std::path::Path) -> Option<PathBuf> {
    let mut dir = start.to_path_buf();
    loop {
        let candidate = dir.join("lambda.py");
        if candidate.exists() {
            return Some(candidate);
        }
        if !dir.pop() {
            return None;
        }
    }
}

// ── monomial enumeration ──────────────────────────────────────────────────────

/// All admissible monomials with positive generators (≥ 1) of total degree
/// 0 … max_deg.  Includes the unit (empty sequence).
fn enumerate_admissible(max_deg: usize) -> Vec<Vec<usize>> {
    let mut out = vec![vec![]]; // unit
    enum_rec(max_deg, usize::MAX, &mut vec![], &mut out);
    out
}

fn enum_rec(remaining: usize, max_next: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if remaining == 0 {
        return;
    }
    let cap = remaining.min(max_next);
    for next in 1..=cap {
        cur.push(next);
        out.push(cur.clone());
        enum_rec(remaining - next, 2 * next, cur, out);
        cur.pop();
    }
}

// ── test-case assembly ────────────────────────────────────────────────────────

fn build_test_cases(max_deg: usize) -> Vec<Element> {
    let all_monos = enumerate_admissible(max_deg);

    let deg = |s: &Vec<usize>| -> usize { s.iter().sum() };

    let combo2_deg = max_deg.min(6);
    let combo3_deg = max_deg.min(5);
    let combo4_deg = max_deg.min(4);

    let combo2: Vec<&Vec<usize>> = all_monos.iter().filter(|s| deg(s) <= combo2_deg).collect();
    let combo3: Vec<&Vec<usize>> = all_monos.iter().filter(|s| deg(s) <= combo3_deg).collect();
    let combo4: Vec<&Vec<usize>> = all_monos.iter().filter(|s| deg(s) <= combo4_deg).collect();

    let mut cases: Vec<Element> = Vec::new();

    // ── singletons: all positive-generator monomials ──────────────────────────
    for seq in &all_monos {
        cases.push(singleton(seq));
    }

    // ── singletons: λ₀-containing monomials ─────────────────────────────────
    for extra in &[
        vec![0usize],
        vec![0, 0],
        vec![1, 0],
        vec![2, 0],
        vec![3, 0],
        vec![2, 1, 0],
        vec![3, 1, 0],
        vec![4, 0],
        vec![4, 2, 0],
    ] {
        cases.push(singleton(extra));
    }

    // ── 2-element sums ────────────────────────────────────────────────────────
    let n2 = combo2.len();
    for i in 0..n2 {
        for j in (i + 1)..n2 {
            let mut e = singleton(combo2[i]);
            e.add_mono(make_mono(combo2[j]));
            cases.push(e);
        }
    }

    // ── 3-element sums (capped) ───────────────────────────────────────────────
    let n3 = combo3.len();
    let mut cnt3 = 0usize;
    'l3: for i in 0..n3 {
        for j in (i + 1)..n3 {
            for k in (j + 1)..n3 {
                let mut e = singleton(combo3[i]);
                e.add_mono(make_mono(combo3[j]));
                e.add_mono(make_mono(combo3[k]));
                cases.push(e);
                cnt3 += 1;
                if cnt3 >= 500 {
                    break 'l3;
                }
            }
        }
    }

    // ── 4-element sums (capped) ───────────────────────────────────────────────
    let n4 = combo4.len();
    let mut cnt4 = 0usize;
    'l4: for i in 0..n4 {
        for j in (i + 1)..n4 {
            for k in (j + 1)..n4 {
                for l in (k + 1)..n4 {
                    let mut e = singleton(combo4[i]);
                    e.add_mono(make_mono(combo4[j]));
                    e.add_mono(make_mono(combo4[k]));
                    e.add_mono(make_mono(combo4[l]));
                    cases.push(e);
                    cnt4 += 1;
                    if cnt4 >= 200 {
                        break 'l4;
                    }
                }
            }
        }
    }

    cases
}

// ── main ──────────────────────────────────────────────────────────────────────

fn main() {
    let max_deg: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    // Locate lambda.py by searching up from CWD.
    let cwd = std::env::current_dir().expect("cannot read CWD");
    let lambda_path = find_lambda_py(&cwd).unwrap_or_else(|| {
        eprintln!("error: lambda.py not found in {cwd:?} or any parent directory.");
        eprintln!("Run from within the ehpx repository tree.");
        std::process::exit(2);
    });
    eprintln!("check_diff: using reference {}", lambda_path.display());

    // Write the embedded driver shim to a temp file.
    let tmp_driver = std::env::temp_dir().join("ehpx_check_driver.py");
    std::fs::write(&tmp_driver, DRIVER).unwrap_or_else(|e| {
        eprintln!("error: could not write driver to {}: {e}", tmp_driver.display());
        std::process::exit(2);
    });

    // Build test cases and compute Rust answers.
    let cases = build_test_cases(max_deg);
    let rust_pairs: Vec<(String, String)> = cases
        .iter()
        .map(|elem| {
            let input = format_elem(elem);
            let d = elem.clone().diff();
            (input, format_elem(&d))
        })
        .collect();

    eprintln!(
        "check_diff: {} test cases (max_deg={max_deg})",
        rust_pairs.len()
    );

    // Spawn Python driver with lambda.py path as argument.
    let mut child = Command::new("python3")
        .arg(&tmp_driver)
        .arg(&lambda_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .unwrap_or_else(|e| {
            eprintln!("error: failed to spawn python3: {e}");
            std::process::exit(2);
        });

    // Feed all inputs.
    {
        let mut stdin = child.stdin.take().unwrap();
        for (input, _) in &rust_pairs {
            writeln!(stdin, "{input}").unwrap();
        }
    } // drop → EOF

    // Collect Python outputs.
    let reader = BufReader::new(child.stdout.take().unwrap());
    let python_answers: Vec<String> = reader
        .lines()
        .map(|l| l.expect("error reading python stdout"))
        .collect();

    child.wait().unwrap();

    // Compare.
    let mut passed = 0usize;
    let mut failed = 0usize;

    for (i, (input, rust_out)) in rust_pairs.iter().enumerate() {
        let python_out = python_answers.get(i).map(String::as_str).unwrap_or("<missing>");
        if rust_out == python_out {
            passed += 1;
        } else {
            if failed < 20 {
                eprintln!("FAIL  d({input})");
                eprintln!("      rust  : {rust_out}");
                eprintln!("      python: {python_out}");
            } else if failed == 20 {
                eprintln!("… (further failures suppressed)");
            }
            failed += 1;
        }
    }

    println!(
        "{} passed, {} failed  (out of {} tests, max_deg={max_deg})",
        passed,
        failed,
        rust_pairs.len()
    );

    if failed > 0 {
        std::process::exit(1);
    }
}
