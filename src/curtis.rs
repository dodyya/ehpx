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
use crate::{Admissible, Element, Monomial, Seq};
use smallvec::smallvec;

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

    fn populate_column(&mut self, k: usize, max_stem: usize) {
        // Row 1: stem k portion of Λ(1).
        // Λ(1) is spanned by admissible monomials whose first generator is ≤ 0,
        // i.e. the unit and λ_0^p for all p ≥ 1.
        if k == 0 {
            self.insert_entry(0, 1, vec![]);           // unit
            // Include λ_0^p tails — enough for genuine pairings plus a buffer.
            // Each "real" differential typically targets entries with ≤1 trailing
            // zero; the rest pair off mechanically via [k,0^p]→[k-1,0^{p+1}].
            let tail_cap = (max_stem / 4).max(2);
            for p in 1..=tail_cap {
                self.insert_entry(0, 1, vec![0; p]);   // λ_0^p
            }
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
                    seq: Admissible(smallvec![n - 1]),
                    deg: n - 1,
                };
                let body = seq_to_element(&s);
                let product = Element::from(prefix) * body;
                // Each monomial in the product is an entry in row n.
                for mono in &product.0 {
                    self.insert_entry(k, n, mono.seq.0.to_vec());
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
    ///   splice in the source, and recurse with the full adjusted element.
    fn complete_cocycle(&mut self, stem: usize, row: usize, seq: Vec<usize>) {
        let original = seq.clone();
        // Track the full (multi-term) element, not just the leading monomial.
        let mut current_elem = seq_to_element(&seq);

        for _guard in 0..200 {
            let boundary = current_elem.clone().diff();

            if boundary.0.is_empty() {
                return; // d(x) = 0, permanent cycle candidate
            }

            // Leading term = highest-filtration monomial (largest first generator).
            let leading = filtration_leading(&boundary);

            // Check: is the leading term in the table and free?
            if self.is_in_table(&leading) && !self.is_target(&leading) && !self.is_source(&leading) {
                let tgt_row = self.row_of(&leading).unwrap();
                self.record_differential(stem, row, original, tgt_row, leading);
                return;
            }

            // Otherwise: "complete the cocycle."
            // Find a tail of `leading` that is already the target of
            // some differential, and adjust the element to cancel that term.
            let mut found = false;
            for ell in 1..=leading.len() {
                let tail = &leading[ell..];
                if tail.is_empty() {
                    break;
                }
                if let Some((_stem, _src_row)) = self.target_set.get(tail) {
                    // `tail` is hit by `z` under some differential.
                    let z = self.find_source_of_target(tail).unwrap();
                    // x ← x + prefix · z  (cancels the leading·tail component)
                    let prefix = &leading[..ell];
                    let prefix_elem = seq_to_element(prefix);
                    let z_elem = seq_to_element(&z);
                    let patch = prefix_elem * z_elem;

                    current_elem = current_elem + patch;

                    if current_elem.0.is_empty() {
                        // x became zero — it was a boundary all along
                        return;
                    }
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
                seq: Admissible(smallvec![1]),
                deg: 1,
            };
            let body = seq_to_element(&s);
            let product = Element::from(prefix) * body;
            for mono in &product.0 {
                self.insert_entry(k, 2, mono.seq.0.to_vec());
            }
        }
    }

    // ── Main driver ─────────────────────────────────────────────────────────

    /// Run the Curtis algorithm up to and including stem `max_stem`.
    pub fn compute(max_stem: usize) -> Self {
        let mut table = Self::new();

        for k in 0..=max_stem {
            // Phase 1: populate (rows ≥ 3 and row 1)
            table.populate_column(k, max_stem);

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
        // ANSI
        const RED: &str = "\x1b[31m";
        const GRN: &str = "\x1b[32m";
        const YEL: &str = "\x1b[33m";
        const CYN: &str = "\x1b[36m";
        const DIM: &str = "\x1b[2m";
        const BLD: &str = "\x1b[1m";
        const RST: &str = "\x1b[0m";

        let mut out = String::new();
        let mut max_row = 0usize;
        for col in self.entries.values() {
            for &r in col.keys() {
                max_row = max_row.max(r);
            }
        }
        if max_row == 0 { max_row = 1; }

        // ── 1. Dot chart  (n vs k, one char per generator) ──────────────────
        //
        //  Each cell gets a count of generators:
        //    .  = empty
        //    o  = 1 cycle (green)
        //    2  = 2 cycles, etc.
        //    x  = 1 source (red), X = mixed source+cycle
        //    *  = 1 target (yellow)
        //  This keeps columns to a fixed 4-char width so alignment is trivial.

        out.push_str(&format!(
            "\n{BLD}{CYN} Curtis table — EHP spectral sequence — stems 0..{}{RST}\n",
            max_stem,
        ));
        out.push_str(&format!(
            " {DIM}{GRN}o{RST}{DIM}=cycle  {RED}x{RST}{DIM}=source  {YEL}*{RST}{DIM}=target  number=count{RST}\n\n",
        ));

        let cw = 4usize; // fixed column width

        // Header row
        out.push_str(&format!("{:>4} |", "n\\k"));
        for k in 0..=max_stem {
            out.push_str(&format!("{:^cw$}", k));
        }
        out.push('\n');

        // Separator
        out.push_str("-----+");
        for _ in 0..=max_stem {
            out.push_str("----");
        }
        out.push('\n');

        // Body rows, high n first
        for row in (1..=max_row).rev() {
            out.push_str(&format!("{:>4} |", row));
            for k in 0..=max_stem {
                let glyph = self.cell_glyph(k, row);
                // glyph is (visible_char, ansi_colour)
                let pad_l = (cw - 1) / 2;
                let pad_r = cw - 1 - pad_l;
                out.push_str(&format!(
                    "{}{}{}{}{RST}{}",
                    " ".repeat(pad_l),
                    glyph.1, // colour
                    glyph.0, // char
                    if glyph.1.is_empty() { "" } else { RST },
                    " ".repeat(pad_r),
                ));
            }
            out.push('\n');
        }
        out.push('\n');

        // ── 2. Differentials ────────────────────────────────────────────────
        out.push_str(&format!("{BLD}{CYN} Differentials{RST}\n"));
        let mut any = false;
        for k in 0..=max_stem {
            if let Some(diffs) = self.differentials.get(&k) {
                for d in diffs {
                    if is_zero_tail_artifact(&d.source) || is_zero_tail_artifact(&d.target) {
                        continue;
                    }
                    any = true;
                    out.push_str(&format!(
                        "  {DIM}k={}{RST}  {RED}{}{RST} {DIM}(n={})  -->  {RST}{YEL}{}{RST} {DIM}(n={}){RST}\n",
                        d.stem,
                        format_seq(&d.source), d.source_row,
                        format_seq(&d.target), d.target_row,
                    ));
                }
            }
        }
        if !any { out.push_str(&format!("  {DIM}(none){RST}\n")); }
        out.push('\n');

        // ── 3. Per-stem detail ──────────────────────────────────────────────
        out.push_str(&format!("{BLD}{CYN} Detail by stem{RST}\n"));
        for k in 0..=max_stem {
            let col = match self.entries.get(&k) {
                Some(c) => c,
                None => continue,
            };
            // Collect all entries with status
            let mut lines: Vec<String> = Vec::new();
            let mut row_keys: Vec<usize> = col.keys().copied().collect();
            row_keys.sort_unstable();
            for &row in &row_keys {
                for seq in &col[&row] {
                    if is_zero_tail_artifact(seq) { continue; }
                    let s = format_seq(seq);
                    if self.is_source(seq) {
                        let tgt = self.source_set.get(seq).unwrap();
                        if is_zero_tail_artifact(tgt) { continue; }
                        let tr = self.row_of(tgt).unwrap_or(0);
                        lines.push(format!(
                            "    {RED}{s}{RST} {DIM}n={row}  -->  {RST}{YEL}{}{RST} {DIM}n={tr}{RST}",
                            format_seq(tgt),
                        ));
                    } else if self.is_target(seq) {
                        // printed with its source
                    } else {
                        lines.push(format!(
                            "    {GRN}{s}{RST} {DIM}n={row}{RST}",
                        ));
                    }
                }
            }
            if lines.is_empty() { continue; }
            out.push_str(&format!("  {BLD}k={k}{RST}\n"));
            for l in &lines {
                out.push_str(l);
                out.push('\n');
            }
        }
        out.push('\n');

        // ── 4. Survivors ────────────────────────────────────────────────────
        out.push_str(&format!("{BLD}{CYN} Survivors — H*(Lambda) cycle candidates{RST}\n"));
        for k in 0..=max_stem {
            let survs = self.survivors(k, usize::MAX);
            let survs: Vec<_> = survs.into_iter()
                .filter(|s| !is_zero_tail_artifact(s))
                .collect();
            if survs.is_empty() { continue; }
            let mut items: Vec<String> = survs.iter().map(|s| {
                let r = self.row_of(s).unwrap_or(0);
                format!("{GRN}{}{RST}{DIM}(n={r}){RST}", format_seq(s))
            }).collect();
            items.sort();
            out.push_str(&format!("  {BLD}k={k:<3}{RST} {}\n", items.join("  ")));
        }
        out.push('\n');

        out
    }

    /// Produce a (char, ansi_colour) pair for one cell in the dot chart.
    fn cell_glyph(&self, stem: usize, row: usize) -> (char, &'static str) {
        const RED: &str = "\x1b[31m";
        const GRN: &str = "\x1b[32m";
        const YEL: &str = "\x1b[33m";
        const DIM: &str = "\x1b[2m";

        let entries = match self.entries.get(&stem).and_then(|c| c.get(&row)) {
            Some(v) => v,
            None => return ('.', DIM),
        };
        if entries.is_empty() { return ('.', DIM); }

        let mut n_cycle = 0u32;
        let mut n_src = 0u32;
        let mut n_tgt = 0u32;
        for seq in entries {
            if is_zero_tail_artifact(seq) { continue; }
            if self.is_source(seq) { n_src += 1; }
            else if self.is_target(seq) { n_tgt += 1; }
            else { n_cycle += 1; }
        }
        let total = n_cycle + n_src + n_tgt;
        if total == 0 { return ('.', DIM); }

        if n_src > 0 && n_cycle == 0 && n_tgt == 0 {
            if n_src == 1 { ('x', RED) } else { (char_digit(n_src), RED) }
        } else if n_tgt > 0 && n_cycle == 0 && n_src == 0 {
            if n_tgt == 1 { ('*', YEL) } else { (char_digit(n_tgt), YEL) }
        } else if n_cycle > 0 && n_src == 0 && n_tgt == 0 {
            if n_cycle == 1 { ('o', GRN) } else { (char_digit(n_cycle), GRN) }
        } else {
            // mixed
            (char_digit(total), GRN)
        }
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn seq_to_element(seq: &[usize]) -> Element {
    let mono = Monomial {
        seq: Admissible(Seq::from_slice(seq)),
        deg: seq.iter().sum(),
    };
    Element::singleton(mono)
}

impl From<Monomial> for Element {
    fn from(m: Monomial) -> Self {
        Element::singleton(m)
    }
}

/// Leading monomial sequence from an Element — highest filtration first.
/// Element is maintained sorted ascending, so `.last()` is the max in O(1).
fn filtration_leading(elem: &Element) -> Vec<usize> {
    elem.0.last().map(|m| m.seq.0.to_vec()).unwrap_or_default()
}

fn char_digit(n: u32) -> char {
    if n <= 9 { char::from_digit(n, 10).unwrap() } else { '#' }
}

/// True if `seq` is a λ_0-tail truncation artifact that should be hidden
/// from the main display.  Keeps: the unit, single [0], entries with ≤1
/// trailing zero, and any entry whose non-zero prefix is "interesting."
fn is_zero_tail_artifact(seq: &[usize]) -> bool {
    if seq.is_empty() { return false; }           // unit
    // Pure λ_0^k (all zeros): artifact for k ≥ 2
    if seq.iter().all(|&x| x == 0) { return seq.len() >= 2; }
    // Count trailing zeros
    let trailing = seq.iter().rev().take_while(|&&x| x == 0).count();
    trailing >= 2
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
