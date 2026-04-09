// ── Curtis algorithm ─────────────────────────────────────────────────────────
//
// Implements the Curtis algorithm for computing H*(Λ) via the EHP spectral
// sequence.  Two phases per stem k:
//
//   Phase 1 — populate column k of the Curtis table from previously computed
//             stems (all rows n ≥ 3; row 2 is deferred).
//   Phase 2 — compute differentials out of column k by "completing the cocycle."
//             Then fill row 2 using the differentials just found.
//
// Visualization prints the table to the terminal.

use std::collections::{BTreeMap, HashMap};
use crate::{Admissible, Element, Monomial};

// ── Table types ──────────────────────────────────────────────────────────────

/// A single entry in the Curtis table, living in row `n`, column `stem`.
#[derive(Clone, Debug)]
pub struct Entry {
    /// The element, stored as its admissible sequence.
    pub seq: Vec<usize>,
    /// Row (filtration) n.
    pub row: usize,
    /// Column (stem) t-s = k.
    pub stem: usize,
}

/// A recorded differential: source entry  →  target entry.
#[derive(Clone, Debug)]
pub struct Differential {
    pub source: Vec<usize>,
    pub target: Vec<usize>,
    pub source_row: usize,
    pub target_row: usize,
    pub stem: usize,
}

/// The full Curtis table.
pub struct CurtisTable {
    /// entries[stem][row] = list of admissible sequences
    pub entries: BTreeMap<usize, BTreeMap<usize, Vec<Vec<usize>>>>,
    /// differentials[stem] = list of differentials originating in that stem
    pub differentials: BTreeMap<usize, Vec<Differential>>,
    /// Quick lookup: sequence → true if it is the *target* of some differential
    target_set: HashMap<Vec<usize>, (usize, usize)>,  // seq → (stem, source_row)
    /// Quick lookup: sequence → true if it is the *source* of some differential
    source_set: HashMap<Vec<usize>, Vec<usize>>,  // source_seq → target_seq
    /// Maximum stem computed so far.
    pub max_stem: usize,
}

impl CurtisTable {
    pub fn new() -> Self {
        Self {
            entries: BTreeMap::new(),
            differentials: BTreeMap::new(),
            target_set: HashMap::new(),
            source_set: HashMap::new(),
            max_stem: 0,
        }
    }

    // ── helpers ──────────────────────────────────────────────────────────────

    fn insert_entry(&mut self, stem: usize, row: usize, seq: Vec<usize>) {
        self.entries
            .entry(stem)
            .or_default()
            .entry(row)
            .or_default()
            .push(seq);
    }

    fn record_differential(&mut self, stem: usize, src_row: usize, src: Vec<usize>, tgt_row: usize, tgt: Vec<usize>) {
        self.source_set.insert(src.clone(), tgt.clone());
        self.target_set.insert(tgt.clone(), (stem, src_row));
        self.differentials
            .entry(stem)
            .or_default()
            .push(Differential {
                source: src,
                target: tgt,
                source_row: src_row,
                target_row: tgt_row,
                stem,
            });
    }

    /// Is `seq` already the target of some differential?
    fn is_target(&self, seq: &[usize]) -> bool {
        self.target_set.contains_key(seq)
    }

    /// Is `seq` already the source of some differential?
    fn is_source(&self, seq: &[usize]) -> bool {
        self.source_set.contains_key(seq)
    }

    /// Does `seq` appear anywhere in the table?
    fn is_in_table(&self, seq: &[usize]) -> bool {
        let deg: usize = seq.iter().sum();
        if let Some(rows) = self.entries.get(&deg) {
            for entries in rows.values() {
                if entries.iter().any(|s| s == seq) {
                    return true;
                }
            }
        }
        false
    }

    /// Find which row `seq` lives in (if any).
    fn row_of(&self, seq: &[usize]) -> Option<usize> {
        let deg: usize = seq.iter().sum();
        if let Some(rows) = self.entries.get(&deg) {
            for (&row, entries) in rows {
                if entries.iter().any(|s| s == seq) {
                    return Some(row);
                }
            }
        }
        None
    }

    /// Read survivors from column `stem`, restricting to rows ≤ `max_row`.
    /// A survivor is an entry that is neither a source nor a target of a differential.
    fn survivors(&self, stem: usize, max_row: usize) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        if let Some(rows) = self.entries.get(&stem) {
            for (&row, entries) in rows {
                if row > max_row {
                    continue;
                }
                for seq in entries {
                    if !self.is_source(seq) && !self.is_target(seq) {
                        result.push(seq.clone());
                    }
                }
            }
        }
        result
    }

    // ── enumeration: all admissible monomials of degree `deg` with first
    //    generator ≤ `max_first` ──────────────────────────────────────────────

    fn enumerate_admissible(deg: usize, max_first: usize) -> Vec<Vec<usize>> {
        let mut out = Vec::new();
        Self::enum_rec(deg, max_first, &mut Vec::new(), &mut out);
        out
    }

    fn enum_rec(remaining: usize, max_next: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if remaining == 0 {
            out.push(cur.clone());
            return;
        }
        let cap = remaining.min(max_next);
        for next in 1..=cap {
            cur.push(next);
            Self::enum_rec(remaining - next, 2 * next, cur, out);
            cur.pop();
        }
    }

    // ── Phase 1: populate column k ──────────────────────────────────────────

    fn populate_column(&mut self, k: usize) {
        // Row 1: stem k portion of Λ(1).
        // Λ(1) is spanned by admissible monomials whose first generator is ≤ 0,
        // i.e. just the unit in degree 0.  For k=0 the unit lives here.
        if k == 0 {
            // The unit (empty sequence) lives in row 1, stem 0.
            self.insert_entry(0, 1, vec![]);
        }

        // Rows n ≥ 3: H_{k-n+1}(Λ(2n-1)), multiply survivors by λ_{n-1}.
        // n ranges so that k-n+1 ≥ 0  ⟹  n ≤ k+1.
        // Also n ≥ 3.
        let n_max = k + 1;
        for n in 3..=n_max {
            let prev_stem = (k + 1) - n;
            let cutoff = 2 * n - 1;     // rows ≤ 2n-1
            let survs = self.survivors(prev_stem, cutoff);
            for s in survs {
                // Multiply on the left by λ_{n-1}: prepend (n-1) to the sequence,
                // then reduce to admissible form via the algebra multiplication.
                let prefix = Monomial {
                    seq: Admissible(vec![n - 1]),
                    deg: n - 1,
                };
                let body = seq_to_element(&s);
                let product = Element::from(prefix) * body;
                // Each monomial in the product is an entry in row n.
                for mono in &product.0 {
                    self.insert_entry(k, n, mono.seq.0.clone());
                }
            }
        }
    }

    // ── Phase 2: compute differentials ──────────────────────────────────────

    fn compute_differentials(&mut self, k: usize) {
        // Collect the rows present in column k, sorted ascending (top to bottom = 3,4,5,…).
        let rows: Vec<usize> = if let Some(col) = self.entries.get(&k) {
            col.keys().filter(|&&r| r >= 3).copied().collect()
        } else {
            vec![]
        };

        for row in rows {
            // Snapshot the entries for this row (we may modify the table as we go).
            let entries: Vec<Vec<usize>> = self.entries
                .get(&k)
                .and_then(|col| col.get(&row))
                .cloned()
                .unwrap_or_default();

            for seq in entries {
                if self.is_source(&seq) || self.is_target(&seq) {
                    continue; // already paired
                }
                self.complete_cocycle(k, row, seq);
            }
        }
    }

    /// The "complete the cocycle" subroutine.
    ///
    /// Given `x` living in row `row` of column `stem`, compute d(x).
    /// - If d(x) = 0, x is a cycle candidate; return.
    /// - If the leading term y of d(x) is in the table and free, record x → y.
    /// - Otherwise, find a tail of y that's already a differential target,
    ///   splice in the source, and recurse.
    fn complete_cocycle(&mut self, stem: usize, row: usize, seq: Vec<usize>) {
        // Guard against infinite recursion with a bounded iteration.
        let mut current = seq;
        for _guard in 0..200 {
            let elem = seq_to_element(&current);
            let boundary = elem.diff();

            if boundary.0.is_empty() {
                return; // d(x) = 0, permanent cycle candidate
            }

            // Leading term = lexicographically first monomial.
            let leading = lex_leading(&boundary);

            // Check: is the leading term in the table and free?
            if self.is_in_table(&leading) && !self.is_target(&leading) && !self.is_source(&leading) {
                let tgt_row = self.row_of(&leading).unwrap();
                self.record_differential(stem, row, current, tgt_row, leading);
                return;
            }

            // Otherwise: "complete the cocycle."
            // Find the longest tail of `leading` that is already the target of
            // some differential.
            let mut found = false;
            for ell in 1..=leading.len() {
                let tail = &leading[ell..];
                if tail.is_empty() {
                    break;
                }
                if let Some((_stem, _src_row)) = self.target_set.get(tail) {
                    // `tail` is hit by `z` under some differential.
                    let z = self.find_source_of_target(tail).unwrap();
                    // x ← x + prefix · z
                    let prefix = &leading[..ell];
                    let prefix_elem = seq_to_element(prefix);
                    let z_elem = seq_to_element(&z);
                    let patch = prefix_elem * z_elem;

                    let old_elem = seq_to_element(&current);
                    let new_elem = old_elem + patch;

                    if new_elem.0.is_empty() {
                        // x became zero — it was a boundary all along
                        return;
                    }
                    current = lex_leading_element(&new_elem);
                    found = true;
                    break;
                }
            }

            if !found {
                // No tail found — leading term is not in table, and no sub-tail
                // is a differential target.  This can happen when the boundary
                // lives entirely outside the truncation.  x is a cycle.
                return;
            }
        }
    }

    /// Find the source sequence that maps to `target` under a differential.
    fn find_source_of_target(&self, target: &[usize]) -> Option<Vec<usize>> {
        // Reverse lookup through source_set.
        for (src, tgt) in &self.source_set {
            if tgt.as_slice() == target {
                return Some(src.clone());
            }
        }
        None
    }

    // ── Phase 1b: fill row 2 after differentials are known ──────────────────

    fn fill_row_2(&mut self, k: usize) {
        if k == 0 {
            return;
        }
        let prev_stem = k - 1;
        let cutoff = 3; // rows ≤ 3
        let survs = self.survivors(prev_stem, cutoff);
        for s in survs {
            let prefix = Monomial {
                seq: Admissible(vec![1]),
                deg: 1,
            };
            let body = seq_to_element(&s);
            let product = Element::from(prefix) * body;
            for mono in &product.0 {
                self.insert_entry(k, 2, mono.seq.0.clone());
            }
        }
    }

    // ── Main driver ─────────────────────────────────────────────────────────

    /// Run the Curtis algorithm up to and including stem `max_stem`.
    pub fn compute(max_stem: usize) -> Self {
        let mut table = Self::new();

        for k in 0..=max_stem {
            // Phase 1: populate (rows ≥ 3 and row 1)
            table.populate_column(k);

            // Phase 2: differentials
            table.compute_differentials(k);

            // Phase 1b: row 2 (needs differentials from this column)
            table.fill_row_2(k);

            table.max_stem = k;
        }

        table
    }

    // ── Visualization ───────────────────────────────────────────────────────

    /// Render the Curtis table to a string for terminal display.
    pub fn display(&self, max_stem: usize) -> String {
        let mut out = String::new();

        // ANSI colour codes
        const RED: &str = "\x1b[31m";
        const GRN: &str = "\x1b[32m";
        const YEL: &str = "\x1b[33m";
        const CYN: &str = "\x1b[36m";
        const DIM: &str = "\x1b[2m";
        const BLD: &str = "\x1b[1m";
        const RST: &str = "\x1b[0m";

        // ── Title ────────────────────────────────────────────────────────────
        out.push_str(&format!(
            "\n{BLD}{CYN}Curtis table  ·  EHP spectral sequence  ·  stems 0..{}{RST}\n",
            max_stem
        ));
        out.push_str(&format!(
            "{DIM}  {GRN}● cycle{RST}  {DIM}{RED}─▸ source{RST}  {DIM}{YEL}▸─ target{RST}\n\n"
        ));

        // ── Column-by-column view ────────────────────────────────────────────
        for k in 0..=max_stem {
            let rows = match self.entries.get(&k) {
                Some(col) => col,
                None => continue,
            };

            out.push_str(&format!("{BLD}stem {k}{RST}"));

            // Differentials originating here
            let _diffs: Vec<&Differential> = self.differentials
                .get(&k)
                .map(|v| v.iter().collect())
                .unwrap_or_default();

            // Column entries, sorted by row descending (high filtration first)
            let mut row_keys: Vec<usize> = rows.keys().copied().collect();
            row_keys.sort_unstable();
            row_keys.reverse();

            let mut first = true;
            for &row in &row_keys {
                let entries = &rows[&row];
                for seq in entries {
                    let sep = if first { "  " } else { ", " };
                    first = false;

                    let s = format_seq(seq);
                    if self.is_source(seq) {
                        // Find the target
                        let tgt = self.source_set.get(seq).unwrap();
                        let tgt_row = self.row_of(tgt).unwrap_or(0);
                        out.push_str(&format!(
                            "{sep}{RED}{s}{RST}{DIM}(n={row}){RST} {RED}─▸{RST} {YEL}{}{RST}{DIM}(n={tgt_row}){RST}",
                            format_seq(tgt),
                        ));
                    } else if self.is_target(seq) {
                        // already printed as part of the source above
                    } else {
                        out.push_str(&format!("{sep}{GRN}{s}{RST}{DIM}(n={row}){RST}"));
                    }
                }
            }
            out.push('\n');
        }

        // ── Survivors summary ────────────────────────────────────────────────
        out.push_str(&format!("\n{BLD}{CYN}Survivors  ·  H*(Λ) cycle candidates{RST}\n"));

        for k in 0..=max_stem {
            let survs = self.survivors(k, usize::MAX);
            if survs.is_empty() {
                continue;
            }
            let list: Vec<String> = survs.iter().map(|s| {
                let r = self.row_of(s).unwrap_or(0);
                format!("{GRN}{}{RST}{DIM}(n={r}){RST}", format_seq(s))
            }).collect();
            out.push_str(&format!("  {BLD}k={k}{RST}  {}\n", list.join("  ")));
        }

        // ── Grid view (compact) ─────────────────────────────────────────────
        let mut max_row = 0usize;
        for col in self.entries.values() {
            for &r in col.keys() {
                max_row = max_row.max(r);
            }
        }
        if max_row == 0 { max_row = 1; }

        out.push_str(&format!("\n{BLD}{CYN}Grid{RST}  {DIM}(row n × stem k){RST}\n"));

        // Determine column widths
        let label_w = 5;
        let mut col_widths: Vec<usize> = Vec::new();
        for k in 0..=max_stem {
            let mut w = format!("{}", k).len().max(1);
            for row in 1..=max_row {
                let cell_w = self.cell_text(k, row).chars().count();
                w = w.max(cell_w);
            }
            col_widths.push(w + 2); // padding
        }

        // Header
        out.push_str(&format!("{:>label_w$} ", "n\\k"));
        for (k, &cw) in col_widths.iter().enumerate() {
            out.push_str(&format!("{BLD}{:^cw$}{RST}", k));
        }
        out.push('\n');

        // Separator
        out.push_str(&format!("{:─>label_w$}─", ""));
        for &cw in &col_widths {
            out.push_str(&"─".repeat(cw));
        }
        out.push('\n');

        // Rows
        for row in (1..=max_row).rev() {
            out.push_str(&format!("{DIM}{:>label_w$}{RST} ", row));
            for (k, &cw) in col_widths.iter().enumerate() {
                let cell = self.cell_text(k, row);
                if cell == "·" {
                    out.push_str(&format!("{DIM}{:^cw$}{RST}", "·"));
                } else {
                    out.push_str(&format!("{:^cw$}", self.cell_coloured(k, row)));
                }
            }
            out.push('\n');
        }
        out.push('\n');
        out
    }

    /// Plain-text cell content (no ANSI).
    fn cell_text(&self, stem: usize, row: usize) -> String {
        let entries = self.entries
            .get(&stem)
            .and_then(|col| col.get(&row))
            .cloned()
            .unwrap_or_default();
        if entries.is_empty() {
            return "·".to_string();
        }
        let mut parts: Vec<String> = entries.iter()
            .filter(|seq| !self.is_target(seq))   // targets shown inline
            .map(|seq| format_seq(seq))
            .collect();
        parts.sort();
        if parts.is_empty() {
            // all entries are targets — show them dimmed
            let mut t: Vec<String> = entries.iter().map(|seq| format_seq(seq)).collect();
            t.sort();
            t.join(" ")
        } else {
            parts.join(" ")
        }
    }

    /// ANSI-coloured cell content for the grid.
    fn cell_coloured(&self, stem: usize, row: usize) -> String {
        const RED: &str = "\x1b[31m";
        const GRN: &str = "\x1b[32m";
        const YEL: &str = "\x1b[33m";
        const DIM: &str = "\x1b[2m";
        const RST: &str = "\x1b[0m";

        let entries = self.entries
            .get(&stem)
            .and_then(|col| col.get(&row))
            .cloned()
            .unwrap_or_default();
        if entries.is_empty() {
            return format!("{DIM}·{RST}");
        }
        let mut parts: Vec<String> = Vec::new();
        let mut sorted: Vec<Vec<usize>> = entries;
        sorted.sort();
        for seq in &sorted {
            let s = format_seq(seq);
            if self.is_source(seq) {
                parts.push(format!("{RED}{s}→{RST}"));
            } else if self.is_target(seq) {
                parts.push(format!("{YEL}→{s}{RST}"));
            } else {
                parts.push(format!("{GRN}{s}{RST}"));
            }
        }
        parts.join(" ")
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn seq_to_element(seq: &[usize]) -> Element {
    if seq.is_empty() {
        // The unit: single monomial with empty sequence, degree 0.
        let mono = Monomial { seq: Admissible(vec![]), deg: 0 };
        let mut s = std::collections::HashSet::new();
        s.insert(mono);
        Element(s)
    } else {
        let mono = Monomial {
            seq: Admissible(seq.to_vec()),
            deg: seq.iter().sum(),
        };
        let mut s = std::collections::HashSet::new();
        s.insert(mono);
        Element(s)
    }
}

impl From<Monomial> for Element {
    fn from(m: Monomial) -> Self {
        let mut s = std::collections::HashSet::new();
        s.insert(m);
        Element(s)
    }
}

/// Lexicographically leading monomial sequence from an Element.
fn lex_leading(elem: &Element) -> Vec<usize> {
    elem.0
        .iter()
        .map(|m| m.seq.0.clone())
        .min()  // lex order on Vec<usize> — smallest = leading
        .unwrap_or_default()
}

/// Return the lex-leading monomial's sequence from an Element.
fn lex_leading_element(elem: &Element) -> Vec<usize> {
    lex_leading(elem)
}

fn format_seq(seq: &[usize]) -> String {
    if seq.is_empty() {
        "1".to_string()
    } else {
        let inner: Vec<String> = seq.iter().map(|x| x.to_string()).collect();
        format!("λ({})", inner.join(","))
    }
}

// ── Binary entry point ──────────────────────────────────────────────────────

/// Run as `cargo run --bin curtis [max_stem]`
pub fn run_curtis(max_stem: usize) {
    eprintln!("Computing Curtis table through stem {}...", max_stem);
    let table = CurtisTable::compute(max_stem);
    println!("{}", table.display(max_stem));
}
