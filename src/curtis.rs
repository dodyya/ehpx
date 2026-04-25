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

/// What `propose_cocycle` decided.  The serial committer applies it (or
/// re-proposes against fresher state if a sibling commit invalidated it).
#[derive(Clone, Debug)]
enum Outcome {
    /// Record `src → tgt` (target lives at `tgt_row`).  Stable as long as
    /// `tgt` hasn't been claimed by some other entry that committed first.
    Diff { src: Vec<usize>, tgt: Vec<usize>, tgt_row: usize },
    /// Artifact (zero-tail) virtual differential — `tgt_row` is synthetic.
    Artifact { src: Vec<usize>, tgt: Vec<usize> },
    /// No differential to record.  `last_leading=None` means the boundary
    /// went to zero (initial d(seq)=0 or elem zeroed via patches) — stable
    /// under sibling commits.  `last_leading=Some(L)` means propose gave up
    /// because no tail of L was a target; if a sibling commit later added
    /// such a tail to `target_set`, we must re-propose.
    NoDiff { src: Vec<usize>, last_leading: Option<Vec<usize>> },
}

/// Read-only view of the table, sharable across threads via `&Snapshot`.
/// Holds borrows of just the maps `propose_cocycle` reads — `&CurtisTable`
/// would also work, but `Snapshot` documents the read set and lets the
/// borrow checker confirm we never call any mutating method during scope.
struct Snapshot<'a> {
    entries: &'a BTreeMap<usize, BTreeMap<usize, Vec<Vec<usize>>>>,
    target_set: &'a HashMap<Vec<usize>, (usize, usize)>,
    target_to_source: &'a HashMap<Vec<usize>, Vec<usize>>,
    source_set: &'a HashMap<Vec<usize>, Vec<usize>>,
}

impl<'a> Snapshot<'a> {
    fn from_table(t: &'a CurtisTable) -> Self {
        Self {
            entries: &t.entries,
            target_set: &t.target_set,
            target_to_source: &t.target_to_source,
            source_set: &t.source_set,
        }
    }
    fn is_target(&self, seq: &[usize]) -> bool { self.target_set.contains_key(seq) }
    fn is_source(&self, seq: &[usize]) -> bool { self.source_set.contains_key(seq) }
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
    fn find_source_of_target(&self, target: &[usize]) -> Option<Vec<usize>> {
        self.target_to_source.get(target).cloned()
    }
}

/// Pure version of `complete_cocycle`: produces an `Outcome` against `snap`
/// without touching the table.  Same algorithm and termination invariants —
/// see comments on `complete_cocycle`.  Safe to call from worker threads.
fn propose_cocycle(snap: &Snapshot<'_>, stem: usize, row: usize, seq: Vec<usize>) -> Outcome {
    let original = seq.clone();
    let mut current_elem = seq_to_element(&seq);
    let mut current_boundary = current_elem.diff_ref();

    let mut prev_leading: Option<Vec<usize>> = None;
    const SAFETY_CAP: usize = 100_000;
    for _guard in 0..SAFETY_CAP {
        if current_boundary.0.is_empty() {
            return Outcome::NoDiff { src: original, last_leading: None };
        }
        let leading = filtration_leading(&current_boundary);
        if let Some(prev) = &prev_leading {
            debug_assert!(
                leading.as_slice() < prev.as_slice(),
                "propose_cocycle regression: leading {:?} did not strictly decrease below {:?} (stem {}, row {}, original {:?})",
                leading, prev, stem, row, original
            );
        }
        prev_leading = Some(leading.clone());

        if is_zero_tail_artifact(&leading) {
            if !snap.is_target(&leading) && !snap.is_source(&leading) {
                return Outcome::Artifact { src: original, tgt: leading };
            }
        } else if snap.is_in_table(&leading) && !snap.is_target(&leading) && !snap.is_source(&leading) {
            let tgt_row = snap.row_of(&leading).unwrap();
            return Outcome::Diff { src: original, tgt: leading, tgt_row };
        }

        let mut found = false;
        for ell in 0..leading.len() {
            let tail = &leading[ell..];
            if snap.target_set.contains_key(tail) {
                let z = snap.find_source_of_target(tail).unwrap();
                let prefix = &leading[..ell];
                let prefix_elem = seq_to_element(prefix);
                let z_elem = seq_to_element(&z);
                let patch = prefix_elem * z_elem;
                let patch_boundary = patch.diff_ref();
                current_elem = current_elem + patch;
                current_boundary = current_boundary + patch_boundary;
                if current_elem.0.is_empty() {
                    return Outcome::NoDiff { src: original, last_leading: None };
                }
                found = true;
                break;
            }
        }
        if !found {
            return Outcome::NoDiff { src: original, last_leading: Some(leading) };
        }
    }
    // Hit safety cap.  Treat as a (suspicious) cycle outcome; the debug
    // assert above would have already fired in a debug build.
    Outcome::NoDiff { src: original, last_leading: prev_leading }
}

/// The full Curtis table.
#[derive(Debug)]
pub struct CurtisTable {
    /// entries[stem][row] = list of admissible sequences
    pub entries: BTreeMap<usize, BTreeMap<usize, Vec<Vec<usize>>>>,
    /// differentials[stem] = list of differentials originating in that stem
    pub differentials: BTreeMap<usize, Vec<Differential>>,
    /// Quick lookup: sequence → true if it is the *target* of some differential
    target_set: HashMap<Vec<usize>, (usize, usize)>,  // seq → (stem, source_row)
    /// Quick lookup: sequence → true if it is the *source* of some differential
    source_set: HashMap<Vec<usize>, Vec<usize>>,  // source_seq → target_seq
    /// Reverse index: target_seq → source_seq.  Needed by `complete_cocycle`
    /// to patch via `prefix · z` where z is the source hitting the matched
    /// tail.  Without this, the lookup is a linear scan over `source_set`
    /// — hot path, and `source_set` grows quadratically with stems.
    target_to_source: HashMap<Vec<usize>, Vec<usize>>,
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
            target_to_source: HashMap::new(),
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
        self.target_to_source.insert(tgt.clone(), src.clone());
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

    /// Read survivors from column `stem`, restricting to info from rows ≤ `max_row`.
    /// A survivor is an entry that is neither a source nor a target of a differential, as far as we're aware.
    fn survivors(&self, stem: usize, max_row: usize) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        if let Some(rows) = self.entries.get(&stem) {
            for (&row, entries) in rows {
                if row > max_row {
                    continue;
                }
                for seq in entries {
                    if !self.is_source(seq) && self.target_set.get(seq).is_none_or(|(_stem,source_row)|{*source_row>max_row}) {// Needs to allow being a target, as long as the source of the differential is > max_row.


                        result.push(seq.clone());
                    }
                }
            }
        }
        result
    }



    // ── Phase 1: populate column k ──────────────────────────────────────────

    fn populate_column(&mut self, k: usize, _max_stem: usize) {
        // Row 1: stem k portion of Λ(1).
        // Λ(1) is spanned by admissible monomials whose first generator is ≤ 0,
        // i.e. the unit and λ_0^p for all p ≥ 1.
        //
        // Optimization: λ_0^p for p ≥ 2 are `is_zero_tail_artifact` —
        // pure bookkeeping that never shows up in the published output
        // and (empirically verified against the pre-optimization oracle)
        // doesn't affect non-artifact cocycle completion either.  Keeping
        // them would blow up later stems' work via the survivor cascade
        // (`λ_{n-1} · λ_0^p → λ_{n-1,0^p}`, each another artifact…).  So
        // we only insert the unit and λ_0 here, and filter artifacts out
        // of every product below.  `_max_stem` is kept in the signature
        // for ABI stability — it no longer drives a `tail_cap`.
        if k == 0 {
            self.insert_entry(0, 1, vec![]);     // unit
            self.insert_entry(0, 1, vec![0]);    // λ_0
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
                // Each non-artifact monomial in the product is an entry in row n.
                for mono in &product.0 {
                    if is_zero_tail_artifact(&mono.seq.0) { continue; }
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

        // Spawn-per-task threshold.  Below this, parallel-propose's
        // overhead (thread spawn, snapshot capture, commit re-checks) is
        // bigger than the saved work.  Profiled: low stems ≤ 30 spend
        // microseconds total here, so going parallel is pure overhead.
        const PARALLEL_THRESHOLD: usize = 4;

        // Cap how many proposers we run concurrently.  18 is the largest
        // row-size we observe through stem 40; more cores than that gives
        // diminishing returns and oversubscribes the OS scheduler when
        // we're also doing work in the main thread.
        let n_workers: usize = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8);

        for row in rows {
            // Snapshot the entries for this row (we may modify the table as we go).
            let entries: Vec<Vec<usize>> = self.entries
                .get(&k)
                .and_then(|col| col.get(&row))
                .cloned()
                .unwrap_or_default();
            let unpaired: Vec<Vec<usize>> = entries
                .into_iter()
                .filter(|seq| !self.is_source(seq) && !self.is_target(seq))
                .collect();

            if unpaired.len() < PARALLEL_THRESHOLD {
                for seq in unpaired {
                    self.complete_cocycle(k, row, seq);
                }
                continue;
            }

            // ── Parallel propose ────────────────────────────────────────
            // Borrow the read-set immutably; share via std::thread::scope
            // so we don't have to clone HashMaps.  Workers pull tasks from
            // a shared atomic index; each writes its outcome to a Vec slot
            // pre-sized to `unpaired.len()`.  No locks on the hot path.
            //
            // Note: `snap` and `next` MUST be declared outside the scope's
            // closure body — the borrow checker requires their lifetime to
            // outlive `'scope`, which the closure's locals don't satisfy.
            let n_tasks = unpaired.len();
            let n_threads = n_workers.min(n_tasks).max(1);
            let snap = Snapshot::from_table(self);
            let next = std::sync::atomic::AtomicUsize::new(0);
            let outcomes: Vec<Outcome> = {
                let snap_ref = &snap;
                let tasks_ref = &unpaired;
                let next_ref = &next;
                std::thread::scope(|s| {
                    let handles: Vec<_> = (0..n_threads)
                        .map(|_| {
                            s.spawn(move || {
                                use std::sync::atomic::Ordering::Relaxed;
                                let mut local: Vec<(usize, Outcome)> = Vec::new();
                                loop {
                                    let i = next_ref.fetch_add(1, Relaxed);
                                    if i >= tasks_ref.len() { break; }
                                    let seq = tasks_ref[i].clone();
                                    let out = propose_cocycle(snap_ref, k, row, seq);
                                    local.push((i, out));
                                }
                                local
                            })
                        })
                        .collect();

                    let mut slots: Vec<Option<Outcome>> = (0..n_tasks).map(|_| None).collect();
                    for h in handles {
                        for (i, out) in h.join().unwrap() {
                            slots[i] = Some(out);
                        }
                    }
                    slots.into_iter().map(|o| o.expect("every task should produce an outcome")).collect()
                })
            };
            drop(snap);  // explicit: end the immutable borrow of self before mutating below.

            // ── Serial commit, in row order ─────────────────────────────
            // Each commit may invalidate later proposals: a Diff outcome
            // claims its target, so a later sibling that proposed the same
            // target needs a re-propose; a NoDiff with `last_leading=Some(L)`
            // is invalidated only if some tail of L is newly in `target_set`.
            for outcome in outcomes {
                self.commit_outcome(k, row, outcome);
            }
        }
    }

    /// Commit a single proposal back to the live table, redoing it serially
    /// against current state if a sibling commit invalidated it.
    fn commit_outcome(&mut self, stem: usize, row: usize, outcome: Outcome) {
        match outcome {
            Outcome::Diff { src, tgt, tgt_row } => {
                if self.is_source(&src) || self.is_target(&src) {
                    // src already paired by some earlier commit (shouldn't
                    // normally happen — src was unpaired at row start —
                    // but cheap to guard).
                    return;
                }
                if self.is_target(&tgt) || self.is_source(&tgt) {
                    // Target was claimed by an earlier sibling commit.
                    // Re-run propose against current state.
                    self.complete_cocycle(stem, row, src);
                    return;
                }
                self.record_differential(stem, row, src, tgt_row, tgt);
            }
            Outcome::Artifact { src, tgt } => {
                if self.is_source(&src) || self.is_target(&src) { return; }
                if self.is_target(&tgt) || self.is_source(&tgt) {
                    self.complete_cocycle(stem, row, src);
                    return;
                }
                // Synthetic tgt_row=0 — see complete_cocycle for rationale.
                self.record_differential(stem, row, src, 0, tgt);
            }
            Outcome::NoDiff { src, last_leading } => {
                if self.is_source(&src) || self.is_target(&src) { return; }
                if let Some(l) = last_leading {
                    // Was the "no tail found" decision invalidated by a
                    // sibling commit that added a new tail-target?  Cheap
                    // O(L) check against current `target_set`.  If yes,
                    // re-propose against current state.
                    let now_has = (0..l.len()).any(|ell| self.target_set.contains_key(&l[ell..]));
                    if now_has {
                        self.complete_cocycle(stem, row, src);
                    }
                }
                // last_leading=None: stable (boundary went to zero) — no-op.
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
    ///
    /// The loop body lives in `propose_cocycle` (pure, snapshot-based).
    /// `complete_cocycle` just runs propose against current state and
    /// commits the outcome — no conflict detection needed since it's
    /// the single-threaded path used by `compute_partial` and by the
    /// commit-time fallback in `compute_differentials`.
    fn complete_cocycle(&mut self, stem: usize, row: usize, seq: Vec<usize>) {
        let outcome = {
            let snap = Snapshot::from_table(self);
            propose_cocycle(&snap, stem, row, seq)
        };
        match outcome {
            Outcome::Diff { src, tgt, tgt_row } => {
                self.record_differential(stem, row, src, tgt_row, tgt);
            }
            Outcome::Artifact { src, tgt } => {
                // Synthetic tgt_row=0: artifact entries aren't stored in
                // `entries`, and emit_json / the Python visualizer filter
                // artifact differentials out via `is_zero_tail_artifact`.
                self.record_differential(stem, row, src, 0, tgt);
            }
            Outcome::NoDiff { .. } => {
                // x was a cycle (d(x)=0 or boundary-all-along) or its
                // leading had no tail-target in the table.  Either way,
                // no diff to record — the entry stays as a survivor.
            }
        }
    }

    /// Find the source sequence that maps to `target` under a differential.
    /// O(1) via the `target_to_source` reverse index.
    fn find_source_of_target(&self, target: &[usize]) -> Option<Vec<usize>> {
        self.target_to_source.get(target).cloned()
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
                if is_zero_tail_artifact(&mono.seq.0) { continue; }
                self.insert_entry(k, 2, mono.seq.0.to_vec());
            }
        }
    }

    // ── Main driver ─────────────────────────────────────────────────────────

    /// Compute a single stem `k`.  Assumes stems `0..k` have already been
    /// processed (or `k == 0` on a fresh table).  Afterwards `max_stem == k`.
    /// `compute` and `extend_to` are thin loops over this.
    pub fn step(&mut self, k: usize) {
        // Phase 1: populate (rows ≥ 3 and row 1)
        self.populate_column(k, 0);
        // Phase 2: differentials
        self.compute_differentials(k);
        // Phase 1b: row 2 (needs differentials from this column)
        self.fill_row_2(k);
        self.max_stem = k;
    }

    /// Run the Curtis algorithm up to and including stem `max_stem`.
    pub fn compute(max_stem: usize) -> Self {
        let mut table = Self::new();
        for k in 0..=max_stem {
            table.step(k);
        }
        table
    }

    /// Run the Curtis algorithm but pause right before calling `complete_cocycle`
    /// on the entry matching `(stop_stem, stop_seq)` — whatever row it's in.
    /// Returns `(table, row_it_lives_in)` with all earlier state exactly as
    /// the normal driver would have produced it.  If `stop_seq` is never hit
    /// (or was already paired), falls through to full computation and returns
    /// `None` for the row.
    pub fn compute_partial(
        max_stem: usize,
        stop_stem: usize,
        stop_seq: &[usize],
    ) -> (Self, Option<usize>) {
        let mut table = Self::new();

        for k in 0..=max_stem {
            table.populate_column(k, max_stem);

            // Mirror compute_differentials, but watch for the stop entry.
            let rows: Vec<usize> = if let Some(col) = table.entries.get(&k) {
                col.keys().filter(|&&r| r >= 3).copied().collect()
            } else {
                vec![]
            };

            for row in rows {
                let entries: Vec<Vec<usize>> = table
                    .entries
                    .get(&k)
                    .and_then(|col| col.get(&row))
                    .cloned()
                    .unwrap_or_default();

                for seq in entries {
                    if k == stop_stem && seq.as_slice() == stop_seq {
                        table.max_stem = k;
                        return (table, Some(row));
                    }
                    if table.is_source(&seq) || table.is_target(&seq) {
                        continue;
                    }
                    table.complete_cocycle(k, row, seq);
                }
            }

            table.fill_row_2(k);
            table.max_stem = k;
        }

        (table, None)
    }

    // ── Visualization ───────────────────────────────────────────────────────

    /// Render the Curtis table to a string.
    pub fn display(&self, max_stem: usize, style: RenderStyle) -> String {
        let c = Codes::from(style);
        let lam = c.lam;
        let fmt_seq = |seq: &[usize]| format_seq_with(seq, lam);
        let section = |title: &str| -> String {
            match style {
                RenderStyle::Ansi => format!("{}{} {}{}\n", c.bld, c.cyn, title, c.rst),
                _ => {
                    let bar = "═".repeat(title.chars().count() + 2);
                    let bar = match style {
                        RenderStyle::Ascii => "=".repeat(title.chars().count() + 2),
                        _ => bar,
                    };
                    format!("{bar}\n {title}\n{bar}\n")
                }
            }
        };

        let mut out = String::new();

        // ── 1. Differentials ────────────────────────────────────────────────
        out.push_str(&section("Differentials"));
        let mut any = false;
        for k in 0..=max_stem {
            if let Some(diffs) = self.differentials.get(&k) {
                for d in diffs {
                    if is_zero_tail_artifact(&d.source) || is_zero_tail_artifact(&d.target) {
                        continue;
                    }
                    any = true;
                    out.push_str(&format!(
                        "  {dim}k={k}{rst}  {red}{src}{rst} {dim}(n={sr})  {arrow}  {rst}{yel}{tgt}{rst} {dim}(n={tr}){rst}\n",
                        dim = c.dim, rst = c.rst, red = c.red, yel = c.yel, arrow = c.arrow,
                        k = d.stem, src = fmt_seq(&d.source), sr = d.source_row,
                        tgt = fmt_seq(&d.target), tr = d.target_row,
                    ));
                }
            }
        }
        if !any { out.push_str(&format!("  {}(none){}\n", c.dim, c.rst)); }
        out.push('\n');

        // ── 2. Per-stem detail ──────────────────────────────────────────────
        out.push_str(&section("Detail by stem"));
        for k in 0..=max_stem {
            let col = match self.entries.get(&k) {
                Some(c) => c,
                None => continue,
            };
            // Gather per-line rendering, with raw monomial text for width alignment.
            let mut raw: Vec<(String, String)> = Vec::new(); // (lhs_raw, suffix_colored)
            let mut row_keys: Vec<usize> = col.keys().copied().collect();
            row_keys.sort_unstable();
            for &row in &row_keys {
                for seq in &col[&row] {
                    if is_zero_tail_artifact(seq) { continue; }
                    let s = fmt_seq(seq);
                    if self.is_source(seq) {
                        let tgt = self.source_set.get(seq).unwrap();
                        if is_zero_tail_artifact(tgt) { continue; }
                        let tr = self.row_of(tgt).unwrap_or(0);
                        let lhs_raw = s.clone();
                        let lhs = format!("{}{}{}", c.red, s, c.rst);
                        let suffix = format!(
                            " {}n={}  {}  {}{}{} {}n={}{}",
                            c.dim, row, c.arrow, c.yel, fmt_seq(tgt), c.rst,
                            c.dim, tr, c.rst,
                        );
                        raw.push((lhs_raw, format!("{}{}", lhs, suffix)));
                    } else if self.is_target(seq) {
                        // printed with its source
                    } else {
                        let lhs_raw = s.clone();
                        let lhs = format!("{}{}{}", c.grn, s, c.rst);
                        let suffix = format!(" {}n={}{}", c.dim, row, c.rst);
                        raw.push((lhs_raw, format!("{}{}", lhs, suffix)));
                    }
                }
            }
            if raw.is_empty() { continue; }
            out.push_str(&format!("  {}k={}{}\n", c.bld, k, c.rst));
            // Column-align the left monomial for nicer reading.
            let w = raw.iter().map(|(r, _)| r.chars().count()).max().unwrap_or(0);
            for (lhs_raw, line) in &raw {
                let pad = w.saturating_sub(lhs_raw.chars().count());
                // We can't just pad `line` (it has ANSI codes); reconstruct with pad.
                // Easier: split around the raw monomial position isn't clean, so
                // insert padding after the colored-lhs in `line`.  Since `line`
                // begins with color+seq+reset, append spaces right after the
                // first reset sequence.  For Plain/Ascii where c.rst is empty,
                // we splice after the seq text.
                let padded_line = if c.rst.is_empty() {
                    // find end of lhs_raw in line
                    match line.find(lhs_raw) {
                        Some(pos) => {
                            let end = pos + lhs_raw.len();
                            let mut s = String::new();
                            s.push_str(&line[..end]);
                            s.push_str(&" ".repeat(pad));
                            s.push_str(&line[end..]);
                            s
                        }
                        None => line.clone(),
                    }
                } else {
                    match line.find(c.rst) {
                        Some(pos) => {
                            let end = pos + c.rst.len();
                            let mut s = String::new();
                            s.push_str(&line[..end]);
                            s.push_str(&" ".repeat(pad));
                            s.push_str(&line[end..]);
                            s
                        }
                        None => line.clone(),
                    }
                };
                out.push_str("    ");
                out.push_str(&padded_line);
                out.push('\n');
            }
        }
        out.push('\n');

        // ── 3. Survivors ────────────────────────────────────────────────────
        out.push_str(&section("Survivors — H*(Lambda) cycle candidates"));
        for k in 0..=max_stem {
            let survs = self.survivors(k, usize::MAX);
            let survs: Vec<_> = survs.into_iter()
                .filter(|s| !is_zero_tail_artifact(s))
                .collect();
            if survs.is_empty() { continue; }
            // Present each survivor on its own line so wide stems don't overflow.
            out.push_str(&format!("  {}k={}{}\n", c.bld, k, c.rst));
            // Sort for deterministic output.
            let mut sorted: Vec<(String, usize)> = survs.iter()
                .map(|s| (fmt_seq(s), self.row_of(s).unwrap_or(0)))
                .collect();
            sorted.sort();
            let w = sorted.iter().map(|(s, _)| s.chars().count()).max().unwrap_or(0);
            for (s, r) in &sorted {
                let pad = w.saturating_sub(s.chars().count());
                out.push_str(&format!(
                    "    {}{}{}{}{} {}n={}{}\n",
                    c.grn, s, c.rst, " ".repeat(pad), "", c.dim, r, c.rst
                ));
            }
        }
        out.push('\n');

        out
    }

    /// Emit a machine-readable JSON report of the table.
    ///
    /// Only **non-artifact** entries and differentials are written —
    /// λ_0-tail bookkeeping (`is_zero_tail_artifact`) is stripped.  This
    /// is both what the visualizer wants and what `from_json` needs to
    /// resume: the algorithm's non-artifact output is insensitive to
    /// which λ_0-tails happen to be in the table at any given stem
    /// (empirically verified against full-state snapshots), and
    /// `from_json` reconstructs the handful of λ_0^p entries the
    /// resumption actually needs from `max_stem` alone.
    ///
    /// For a stem-24 run, full state would weigh ~80 MB (99.9% artifact);
    /// the filtered form is ~100 KB.
    ///
    /// Schema:
    ///
    /// ```json
    /// {
    ///   "max_stem": 12,
    ///   "entries": [
    ///     {"stem": 11, "row": 10, "seq": [9,1,1], "role": "source"},
    ///     {"stem": 11, "row": 3,  "seq": [2,2,3,3], "role": "target"},
    ///     {"stem": 7,  "row": 8,  "seq": [7], "role": "cycle"}
    ///   ],
    ///   "differentials": [
    ///     {"stem": 11, "src_row": 10, "src": [9,1,1],
    ///                  "tgt_row": 3,  "tgt": [2,2,3,3]}
    ///   ],
    ///   "survivors": [
    ///     {"stem": 7, "row": 8, "seq": [7]}
    ///   ]
    /// }
    /// ```
    pub fn emit_json(&self, max_stem: usize) -> String {
        let mut out = String::new();
        out.push_str("{\n");
        out.push_str(&format!("  \"max_stem\": {},\n", max_stem));

        // entries — non-artifact only.
        out.push_str("  \"entries\": [\n");
        let mut first = true;
        for k in 0..=max_stem {
            let col = match self.entries.get(&k) {
                Some(c) => c,
                None => continue,
            };
            let mut row_keys: Vec<usize> = col.keys().copied().collect();
            row_keys.sort_unstable();
            for &row in &row_keys {
                for seq in &col[&row] {
                    if is_zero_tail_artifact(seq) { continue; }
                    let role = if self.is_source(seq) {
                        // Skip sources whose target is an artifact — the
                        // pairing isn't meaningful for consumers, and on
                        // resume we rebuild via the cocycle loop anyway.
                        if self.source_set.get(seq).is_some_and(|t| is_zero_tail_artifact(t)) {
                            continue;
                        }
                        "source"
                    } else if self.is_target(seq) {
                        "target"
                    } else {
                        "cycle"
                    };
                    if !first { out.push_str(",\n"); }
                    first = false;
                    out.push_str(&format!(
                        "    {{\"stem\": {}, \"row\": {}, \"seq\": {}, \"role\": \"{}\"}}",
                        k, row, json_seq(seq), role
                    ));
                }
            }
        }
        out.push_str("\n  ],\n");

        // differentials — non-artifact only.
        out.push_str("  \"differentials\": [\n");
        let mut first = true;
        for k in 0..=max_stem {
            if let Some(diffs) = self.differentials.get(&k) {
                for d in diffs {
                    if is_zero_tail_artifact(&d.source) || is_zero_tail_artifact(&d.target) {
                        continue;
                    }
                    if !first { out.push_str(",\n"); }
                    first = false;
                    out.push_str(&format!(
                        "    {{\"stem\": {}, \"src_row\": {}, \"src\": {}, \"tgt_row\": {}, \"tgt\": {}}}",
                        d.stem, d.source_row, json_seq(&d.source),
                        d.target_row, json_seq(&d.target)
                    ));
                }
            }
        }
        out.push_str("\n  ],\n");

        // survivors — non-artifact only; this field is for downstream reading
        // and is recomputed on rehydration anyway.
        out.push_str("  \"survivors\": [\n");
        let mut first = true;
        for k in 0..=max_stem {
            let survs = self.survivors(k, usize::MAX);
            for seq in survs {
                if is_zero_tail_artifact(&seq) { continue; }
                let row = self.row_of(&seq).unwrap_or(0);
                if !first { out.push_str(",\n"); }
                first = false;
                out.push_str(&format!(
                    "    {{\"stem\": {}, \"row\": {}, \"seq\": {}}}",
                    k, row, json_seq(&seq)
                ));
            }
        }
        out.push_str("\n  ]\n");

        out.push_str("}\n");
        out
    }

    /// Rebuild a `CurtisTable` from a JSON report previously emitted by
    /// `emit_json`.  The JSON holds non-artifact state only; since the
    /// algorithm no longer produces or consults artifacts internally
    /// (see `populate_column`), there is nothing to reconstruct beyond
    /// replaying the stored entries and differentials.
    pub fn from_json(src: &str) -> Result<Self, String> {
        let root = json::Parser::new(src).parse_value()?;
        let max_stem = root.get("max_stem")?.as_num()? as usize;

        let mut table = Self::new();
        table.max_stem = max_stem;

        for e in root.get("entries")?.as_arr()? {
            let stem = e.get("stem")?.as_num()? as usize;
            let row = e.get("row")?.as_num()? as usize;
            let seq = json_seq_to_vec(e.get("seq")?)?;
            table.insert_entry(stem, row, seq);
        }

        for d in root.get("differentials")?.as_arr()? {
            let stem = d.get("stem")?.as_num()? as usize;
            let src_row = d.get("src_row")?.as_num()? as usize;
            let tgt_row = d.get("tgt_row")?.as_num()? as usize;
            let src = json_seq_to_vec(d.get("src")?)?;
            let tgt = json_seq_to_vec(d.get("tgt")?)?;
            table.record_differential(stem, src_row, src, tgt_row, tgt);
        }

        Ok(table)
    }

    /// Extend an already-built table to cover stems up to `new_max_stem`.
    /// No-op if `new_max_stem <= self.max_stem`.
    pub fn extend_to(&mut self, new_max_stem: usize) {
        if new_max_stem <= self.max_stem {
            return;
        }
        for k in (self.max_stem + 1)..=new_max_stem {
            self.step(k);
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
    format_seq_with(seq, "λ")
}

fn format_seq_with(seq: &[usize], lam: &str) -> String {
    if seq.is_empty() {
        "1".to_string()
    } else {
        let inner: Vec<String> = seq.iter().map(|x| x.to_string()).collect();
        format!("{}({})", lam, inner.join(","))
    }
}

fn json_seq(seq: &[usize]) -> String {
    let inner: Vec<String> = seq.iter().map(|x| x.to_string()).collect();
    format!("[{}]", inner.join(","))
}

fn json_seq_to_vec(v: &json::Value) -> Result<Vec<usize>, String> {
    v.as_arr()?
        .iter()
        .map(|x| x.as_num().map(|n| n as usize))
        .collect()
}

// ── Minimal JSON parser ─────────────────────────────────────────────────────
//
// Only exists to round-trip `emit_json` output — no generic-grade JSON
// compliance (no exponents, unicode escapes, etc.).  Adding a serde
// dependency for a few hundred lines of ASCII feels wasteful.

#[allow(dead_code)]
mod json {
    #[derive(Debug)]
    pub enum Value {
        Null,
        Bool(bool),
        Num(i64),
        Str(String),
        Arr(Vec<Value>),
        Obj(Vec<(String, Value)>),
    }

    impl Value {
        pub fn as_arr(&self) -> Result<&[Value], String> {
            match self {
                Value::Arr(a) => Ok(a),
                _ => Err("expected array".to_string()),
            }
        }
        pub fn as_num(&self) -> Result<i64, String> {
            match self {
                Value::Num(n) => Ok(*n),
                _ => Err("expected number".to_string()),
            }
        }
        #[allow(dead_code)]
        pub fn as_str(&self) -> Result<&str, String> {
            match self {
                Value::Str(s) => Ok(s),
                _ => Err("expected string".to_string()),
            }
        }
        pub fn get(&self, key: &str) -> Result<&Value, String> {
            match self {
                Value::Obj(pairs) => pairs
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, v)| v)
                    .ok_or_else(|| format!("missing key {:?}", key)),
                _ => Err(format!("expected object (looking up {:?})", key)),
            }
        }
    }

    pub struct Parser<'a> {
        src: &'a [u8],
        pos: usize,
    }

    impl<'a> Parser<'a> {
        pub fn new(src: &'a str) -> Self {
            Self { src: src.as_bytes(), pos: 0 }
        }

        fn peek(&self) -> Option<u8> {
            self.src.get(self.pos).copied()
        }

        fn skip_ws(&mut self) {
            while let Some(c) = self.peek() {
                if matches!(c, b' ' | b'\t' | b'\n' | b'\r') {
                    self.pos += 1;
                } else {
                    break;
                }
            }
        }

        pub fn parse_value(&mut self) -> Result<Value, String> {
            self.skip_ws();
            match self.peek() {
                Some(b'{') => self.parse_obj(),
                Some(b'[') => self.parse_arr(),
                Some(b'"') => self.parse_str().map(Value::Str),
                Some(b't') | Some(b'f') => self.parse_bool(),
                Some(b'n') => {
                    self.expect_lit(b"null")?;
                    Ok(Value::Null)
                }
                Some(c) if c == b'-' || c.is_ascii_digit() => self.parse_num(),
                Some(c) => Err(format!("unexpected {:?} at byte {}", c as char, self.pos)),
                None => Err("unexpected eof".to_string()),
            }
        }

        fn parse_obj(&mut self) -> Result<Value, String> {
            self.pos += 1; // consume '{'
            let mut pairs = Vec::new();
            self.skip_ws();
            if self.peek() == Some(b'}') {
                self.pos += 1;
                return Ok(Value::Obj(pairs));
            }
            loop {
                self.skip_ws();
                let k = self.parse_str()?;
                self.skip_ws();
                if self.peek() != Some(b':') {
                    return Err(format!("expected ':' at byte {}", self.pos));
                }
                self.pos += 1;
                let v = self.parse_value()?;
                pairs.push((k, v));
                self.skip_ws();
                match self.peek() {
                    Some(b',') => {
                        self.pos += 1;
                    }
                    Some(b'}') => {
                        self.pos += 1;
                        break;
                    }
                    Some(c) => {
                        return Err(format!(
                            "expected ',' or '}}' got {:?} at byte {}",
                            c as char, self.pos
                        ))
                    }
                    None => return Err("unexpected eof".to_string()),
                }
            }
            Ok(Value::Obj(pairs))
        }

        fn parse_arr(&mut self) -> Result<Value, String> {
            self.pos += 1; // consume '['
            let mut items = Vec::new();
            self.skip_ws();
            if self.peek() == Some(b']') {
                self.pos += 1;
                return Ok(Value::Arr(items));
            }
            loop {
                let v = self.parse_value()?;
                items.push(v);
                self.skip_ws();
                match self.peek() {
                    Some(b',') => {
                        self.pos += 1;
                    }
                    Some(b']') => {
                        self.pos += 1;
                        break;
                    }
                    Some(c) => {
                        return Err(format!(
                            "expected ',' or ']' got {:?} at byte {}",
                            c as char, self.pos
                        ))
                    }
                    None => return Err("unexpected eof".to_string()),
                }
            }
            Ok(Value::Arr(items))
        }

        fn parse_str(&mut self) -> Result<String, String> {
            self.skip_ws();
            if self.peek() != Some(b'"') {
                return Err(format!("expected '\"' at byte {}", self.pos));
            }
            self.pos += 1;
            let start = self.pos;
            while let Some(c) = self.peek() {
                if c == b'"' {
                    let s = std::str::from_utf8(&self.src[start..self.pos])
                        .map_err(|e| e.to_string())?
                        .to_string();
                    self.pos += 1;
                    return Ok(s);
                }
                if c == b'\\' {
                    // skip the escape byte too
                    self.pos += 2;
                } else {
                    self.pos += 1;
                }
            }
            Err("unterminated string".to_string())
        }

        fn parse_num(&mut self) -> Result<Value, String> {
            let start = self.pos;
            if self.peek() == Some(b'-') {
                self.pos += 1;
            }
            while let Some(c) = self.peek() {
                if c.is_ascii_digit() {
                    self.pos += 1;
                } else {
                    break;
                }
            }
            let slice = std::str::from_utf8(&self.src[start..self.pos])
                .map_err(|e| e.to_string())?;
            slice
                .parse::<i64>()
                .map(Value::Num)
                .map_err(|e| e.to_string())
        }

        fn parse_bool(&mut self) -> Result<Value, String> {
            if self.src[self.pos..].starts_with(b"true") {
                self.pos += 4;
                Ok(Value::Bool(true))
            } else if self.src[self.pos..].starts_with(b"false") {
                self.pos += 5;
                Ok(Value::Bool(false))
            } else {
                Err(format!("expected bool at byte {}", self.pos))
            }
        }

        fn expect_lit(&mut self, lit: &[u8]) -> Result<(), String> {
            if self.src[self.pos..].starts_with(lit) {
                self.pos += lit.len();
                Ok(())
            } else {
                Err(format!("literal mismatch at byte {}", self.pos))
            }
        }
    }
}

// ── Render styles ──────────────────────────────────────────────────────────

/// Rendering style for `display`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RenderStyle {
    /// Full ANSI colors + Unicode (good for terminals).
    Ansi,
    /// No colors, Unicode (λ, →).  Good for .txt files viewed in modern editors.
    Plain,
    /// No colors, strict ASCII ("l" for lambda, "->").  Good for anywhere.
    Ascii,
}

/// Resolved render codes — empty strings for everything in plain/ascii mode.
struct Codes {
    red: &'static str,
    grn: &'static str,
    yel: &'static str,
    cyn: &'static str,
    dim: &'static str,
    bld: &'static str,
    rst: &'static str,
    lam: &'static str,
    arrow: &'static str,
}

impl From<RenderStyle> for Codes {
    fn from(s: RenderStyle) -> Self {
        match s {
            RenderStyle::Ansi => Codes {
                red: "\x1b[31m", grn: "\x1b[32m", yel: "\x1b[33m", cyn: "\x1b[36m",
                dim: "\x1b[2m",  bld: "\x1b[1m",  rst: "\x1b[0m",
                lam: "λ", arrow: "-->",
            },
            RenderStyle::Plain => Codes {
                red: "", grn: "", yel: "", cyn: "", dim: "", bld: "", rst: "",
                lam: "λ", arrow: "→",
            },
            RenderStyle::Ascii => Codes {
                red: "", grn: "", yel: "", cyn: "", dim: "", bld: "", rst: "",
                lam: "l", arrow: "->",
            },
        }
    }
}

// ── Binary entry point ──────────────────────────────────────────────────────

/// Legacy entry point: compute and print with ANSI colors to stdout.
/// The `table` binary handles flag parsing / file output / styles directly.
pub fn run_curtis(max_stem: usize) {
    eprintln!("Computing Curtis table through stem {}...", max_stem);
    let table = CurtisTable::compute(max_stem);
    println!("{}", table.display(max_stem, RenderStyle::Ansi));
}

// ── Interactive debug harness ───────────────────────────────────────────────
//
// A stepper for `complete_cocycle`.  Builds the table up to a chosen entry,
// then drops into a REPL that lets you inspect the running element, query
// the state of any sequence, single-step the cocycle loop, or run it all the
// way out with a per-iteration trace.
//
// Kept inside curtis.rs because it reaches directly into the private
// internals (target_set, source_set, is_in_table, row_of, …) — and because
// it's meant to be surgical: a single reason to change, removed when the
// investigation ends.

/// Run the debug REPL for one cocycle-completion site.
///
/// Builds the Curtis table up to `max_stem`, stopping just before
/// `complete_cocycle` would be called on the entry matching
/// `(stop_stem, stop_seq)`.  The row is discovered automatically.
pub fn run_debug(max_stem: usize, stop_stem: usize, stop_seq: Vec<usize>) {
    eprintln!(
        "Building Curtis table through stem {} (stopping at stem={}, seq={})...",
        max_stem,
        stop_stem,
        format_seq(&stop_seq),
    );
    let (table, stop_row_opt) = CurtisTable::compute_partial(max_stem, stop_stem, &stop_seq);

    let stop_row = match stop_row_opt {
        Some(r) => {
            eprintln!("found {} at (stem={}, row={})", format_seq(&stop_seq), stop_stem, r);
            r
        }
        None => {
            eprintln!(
                "WARNING: {} was never reached at stem={} — either it's not an entry there, \
                 or it got pre-paired before its turn came up.  The debug session will still \
                 run against the fully-computed final table.",
                format_seq(&stop_seq),
                stop_stem,
            );
            0
        }
    };

    let mut session = DebugSession {
        current_elem: seq_to_element(&stop_seq),
        original_seq: stop_seq,
        stop_stem,
        stop_row,
        iter_count: 0,
        table,
    };

    session.banner();
    session.repl();
}

struct DebugSession {
    table: CurtisTable,
    stop_stem: usize,
    stop_row: usize,
    original_seq: Vec<usize>,
    current_elem: Element,
    iter_count: usize,
}

impl DebugSession {
    fn banner(&self) {
        println!();
        println!("=== Curtis debug harness ===");
        println!(
            "Stopped just before complete_cocycle(stem={}, row={}, seq={})",
            self.stop_stem,
            self.stop_row,
            format_seq(&self.original_seq),
        );
        println!();
        println!("Commands:");
        println!("  show                 current running element x");
        println!("  d | diff             show d(x)");
        println!("  lead                 filtration_leading(d(x)) + its in-table / source / target status");
        println!("  tails                enumerate tails of the leading term with target-match status");
        println!("  step                 execute one iteration of the cocycle loop (verbose)");
        println!("  auto [N]             run to completion, full trace (default cap N=200)");
        println!("  check <seq>          inspect any sequence, e.g. `check 2,2,3,3`");
        println!("  where <seq>          which row does seq live in?");
        println!("  diffs                all recorded differentials so far");
        println!("  row <k> <n>          list entries in column k, row n");
        println!("  reset                reset x back to the original seq");
        println!("  quit | exit          leave");
        println!();
    }

    fn repl(&mut self) {
        use std::io::{self, BufRead, Write};
        let stdin = io::stdin();
        let mut stdout = io::stdout();
        loop {
            print!("debug> ");
            stdout.flush().ok();
            let mut line = String::new();
            if stdin.lock().read_line(&mut line).unwrap_or(0) == 0 {
                break;
            }
            let line = line.trim().to_string();
            if line.is_empty() {
                continue;
            }
            let (cmd, rest) = match line.split_once(char::is_whitespace) {
                Some((c, r)) => (c, r.trim()),
                None => (line.as_str(), ""),
            };
            match cmd {
                "quit" | "exit" => break,
                "show" => self.cmd_show(),
                "d" | "diff" => self.cmd_d(),
                "lead" => self.cmd_lead(),
                "tails" => self.cmd_tails(),
                "step" => self.cmd_step(),
                "auto" => {
                    let n = rest.parse::<usize>().unwrap_or(200);
                    self.cmd_auto(n);
                }
                "check" => self.cmd_check(rest),
                "where" => self.cmd_where(rest),
                "diffs" => self.cmd_diffs(),
                "row" => self.cmd_row(rest),
                "reset" => self.cmd_reset(),
                "help" | "?" => self.banner(),
                _ => println!("unknown command: {:?}  (try `help`)", cmd),
            }
        }
    }

    // ── parsing ──────────────────────────────────────────────────────────────

    fn parse_seq(s: &str) -> Option<Vec<usize>> {
        let t = s.trim();
        if t.is_empty() {
            return Some(vec![]);
        }
        t.split(|c: char| c == ',' || c.is_whitespace() || c == '[' || c == ']')
            .map(|x| x.trim())
            .filter(|x| !x.is_empty())
            .map(|x| x.parse::<usize>().ok())
            .collect()
    }

    // ── commands ─────────────────────────────────────────────────────────────

    fn cmd_show(&self) {
        println!("x (iter {}):", self.iter_count);
        if self.current_elem.0.is_empty() {
            println!("  0");
        } else {
            for m in &self.current_elem.0 {
                println!("  {}", format_seq(&m.seq.0));
            }
        }
    }

    fn cmd_d(&self) {
        let boundary = self.current_elem.clone().diff();
        println!("d(x):");
        if boundary.0.is_empty() {
            println!("  0");
        } else {
            for m in &boundary.0 {
                println!("  {}", format_seq(&m.seq.0));
            }
        }
    }

    fn cmd_lead(&self) {
        let boundary = self.current_elem.clone().diff();
        if boundary.0.is_empty() {
            println!("d(x) = 0; no leading term.  x is a cycle candidate.");
            return;
        }
        let leading = filtration_leading(&boundary);
        let s = format_seq(&leading);
        let in_table = self.table.is_in_table(&leading);
        let is_source = self.table.is_source(&leading);
        let is_target = self.table.is_target(&leading);
        let row = self.table.row_of(&leading);
        println!("leading = {}", s);
        println!("  in_table:  {}", in_table);
        if let Some(r) = row {
            println!("  row:       {}", r);
        }
        println!("  is_source: {}", is_source);
        println!("  is_target: {}", is_target);
        if is_source {
            let tgt = self.table.source_set.get(&leading).unwrap();
            println!("  → differential lands on: {}", format_seq(tgt));
        }
        if is_target {
            let (st, sr) = self.table.target_set.get(&leading).unwrap();
            if let Some(src) = self.table.find_source_of_target(&leading) {
                println!(
                    "  ← hit by: {} (stem={}, row={})",
                    format_seq(&src),
                    st,
                    sr
                );
            }
        }
        if in_table && !is_source && !is_target {
            println!("  ⇒ RECORD: {} → {}", format_seq(&self.original_seq), s);
        } else {
            println!("  ⇒ need cocycle completion (see `tails`)");
        }
    }

    fn cmd_tails(&self) {
        let boundary = self.current_elem.clone().diff();
        if boundary.0.is_empty() {
            println!("d(x) = 0; no leading term.");
            return;
        }
        let leading = filtration_leading(&boundary);
        println!("leading = {}", format_seq(&leading));
        if leading.is_empty() {
            println!("  (leading is empty)");
            return;
        }
        let mut first_match: Option<usize> = None;
        for ell in 0..leading.len() {
            let prefix = &leading[..ell];
            let tail = &leading[ell..];
            let is_tgt = self.table.target_set.contains_key(tail);
            let marker = if is_tgt { "✓" } else { " " };
            let note = if is_tgt {
                let (st, sr) = self.table.target_set.get(tail).unwrap();
                let src = self
                    .table
                    .find_source_of_target(tail)
                    .map(|s| format_seq(&s))
                    .unwrap_or_else(|| "?".into());
                format!("   ← source {} (stem={}, row={})", src, st, sr)
            } else {
                String::new()
            };
            println!(
                "  {} ell={} prefix={:<14} tail={:<14}{}",
                marker,
                ell,
                format_seq(prefix),
                format_seq(tail),
                note
            );
            if is_tgt && first_match.is_none() {
                first_match = Some(ell);
            }
        }
        match first_match {
            Some(ell) => {
                let prefix = &leading[..ell];
                let tail = &leading[ell..];
                let z = self.table.find_source_of_target(tail).unwrap();
                let prefix_elem = seq_to_element(prefix);
                let z_elem = seq_to_element(&z);
                let patch = prefix_elem * z_elem;
                println!();
                println!("→ algorithm picks ell={} (first match, ascending)", ell);
                println!("  prefix = {}", format_seq(prefix));
                println!("  z      = {}  (source whose target is tail {})", format_seq(&z), format_seq(tail));
                println!("  patch  = prefix · z =");
                if patch.0.is_empty() {
                    println!("    0");
                } else {
                    for m in &patch.0 {
                        println!("    {}", format_seq(&m.seq.0));
                    }
                }
            }
            None => {
                println!();
                println!("(no tail is a differential target — cocycle loop would exit without recording)");
            }
        }
    }

    fn cmd_step(&mut self) {
        let boundary = self.current_elem.clone().diff();
        if boundary.0.is_empty() {
            println!("d(x) = 0; x is a cycle candidate.  No step.");
            return;
        }
        let leading = filtration_leading(&boundary);
        println!("── iteration {} ─────────────────────────", self.iter_count);
        println!(
            "d(x) has {} term(s); leading = {}",
            boundary.0.len(),
            format_seq(&leading)
        );

        let in_table = self.table.is_in_table(&leading);
        let is_source = self.table.is_source(&leading);
        let is_target = self.table.is_target(&leading);

        if in_table && !is_source && !is_target {
            let tgt_row = self.table.row_of(&leading).unwrap();
            println!(
                "leading free — RECORD: {} (row {}) → {} (row {})",
                format_seq(&self.original_seq),
                self.stop_row,
                format_seq(&leading),
                tgt_row
            );
            self.table.record_differential(
                self.stop_stem,
                self.stop_row,
                self.original_seq.clone(),
                tgt_row,
                leading,
            );
            return;
        }

        println!(
            "leading NOT free: in_table={} is_source={} is_target={}",
            in_table, is_source, is_target
        );
        for ell in 0..leading.len() {
            let tail = &leading[ell..];
            if self.table.target_set.contains_key(tail) {
                let prefix = &leading[..ell];
                let z = self.table.find_source_of_target(tail).unwrap();
                let prefix_elem = seq_to_element(prefix);
                let z_elem = seq_to_element(&z);
                let patch = prefix_elem * z_elem;
                println!(
                    "match at ell={}: tail={} z={}",
                    ell,
                    format_seq(tail),
                    format_seq(&z)
                );
                println!("patching x ← x + ({}) · ({})", format_seq(prefix), format_seq(&z));
                self.current_elem = self.current_elem.clone() + patch;
                self.iter_count += 1;
                if self.current_elem.0.is_empty() {
                    println!("x became zero — boundary.  No differential recorded.");
                } else {
                    let new_d = self.current_elem.clone().diff();
                    let new_lead = if new_d.0.is_empty() {
                        "0".into()
                    } else {
                        format_seq(&filtration_leading(&new_d))
                    };
                    println!("new d(x) has {} term(s); new leading = {}", new_d.0.len(), new_lead);
                }
                return;
            }
        }
        println!("no tail-match — cocycle loop terminates, no differential recorded.");
    }

    fn cmd_auto(&mut self, max_iters: usize) {
        let mut steps = 0;
        while steps < max_iters {
            let boundary = self.current_elem.clone().diff();
            if boundary.0.is_empty() {
                println!("[iter {}] d(x) = 0; x is a cycle candidate.", self.iter_count);
                return;
            }
            let leading = filtration_leading(&boundary);
            let in_table = self.table.is_in_table(&leading);
            let is_source = self.table.is_source(&leading);
            let is_target = self.table.is_target(&leading);

            if in_table && !is_source && !is_target {
                let tgt_row = self.table.row_of(&leading).unwrap();
                println!(
                    "[iter {}] RECORD lead={} free → row {}",
                    self.iter_count,
                    format_seq(&leading),
                    tgt_row
                );
                self.table.record_differential(
                    self.stop_stem,
                    self.stop_row,
                    self.original_seq.clone(),
                    tgt_row,
                    leading,
                );
                return;
            }

            let mut matched = None;
            for ell in 0..leading.len() {
                let tail = &leading[ell..];
                if self.table.target_set.contains_key(tail) {
                    matched = Some(ell);
                    break;
                }
            }
            match matched {
                Some(ell) => {
                    let prefix = &leading[..ell];
                    let tail = &leading[ell..];
                    let z = self.table.find_source_of_target(tail).unwrap();
                    let prefix_elem = seq_to_element(prefix);
                    let z_elem = seq_to_element(&z);
                    let patch = prefix_elem * z_elem;
                    let lead_s = format_seq(&leading);
                    let tail_s = format_seq(tail);
                    let z_s = format_seq(&z);
                    let prefix_s = format_seq(prefix);
                    self.current_elem = self.current_elem.clone() + patch;
                    self.iter_count += 1;
                    steps += 1;
                    let new_lead = if self.current_elem.0.is_empty() {
                        "0".into()
                    } else {
                        format_seq(&filtration_leading(&self.current_elem.clone().diff()))
                    };
                    println!(
                        "[iter {}] lead={} not free (in_tbl={} src={} tgt={}); ell={} tail={} z={} patch={}·z  →  new_lead={}",
                        self.iter_count,
                        lead_s,
                        in_table,
                        is_source,
                        is_target,
                        ell,
                        tail_s,
                        z_s,
                        prefix_s,
                        new_lead
                    );
                    if self.current_elem.0.is_empty() {
                        println!("  x became zero — boundary.");
                        return;
                    }
                }
                None => {
                    println!(
                        "[iter {}] lead={} not free and no tail-match — halt, no differential.",
                        self.iter_count,
                        format_seq(&leading)
                    );
                    return;
                }
            }
        }
        println!("hit iteration cap {} without convergence", max_iters);
    }

    fn cmd_check(&self, s: &str) {
        let seq = match Self::parse_seq(s) {
            Some(v) => v,
            None => {
                println!("couldn't parse seq from {:?}", s);
                return;
            }
        };
        println!("seq = {}  (deg {})", format_seq(&seq), seq.iter().sum::<usize>());
        println!("  in_table:  {}", self.table.is_in_table(&seq));
        if let Some(r) = self.table.row_of(&seq) {
            println!("  row:       {}", r);
        }
        println!("  is_source: {}", self.table.is_source(&seq));
        println!("  is_target: {}", self.table.is_target(&seq));
        if let Some(tgt) = self.table.source_set.get(&seq) {
            println!("  → differential to: {}", format_seq(tgt));
        }
        if let Some((st, sr)) = self.table.target_set.get(&seq) {
            if let Some(src) = self.table.find_source_of_target(&seq) {
                println!(
                    "  ← hit by: {} (stem={}, row={})",
                    format_seq(&src),
                    st,
                    sr
                );
            }
        }
        // Also show d(seq) as a convenience.
        let d = seq_to_element(&seq).diff();
        if d.0.is_empty() {
            println!("  d(seq):    0");
        } else {
            println!("  d(seq):");
            for m in &d.0 {
                println!("    {}", format_seq(&m.seq.0));
            }
        }
    }

    fn cmd_where(&self, s: &str) {
        let seq = match Self::parse_seq(s) {
            Some(v) => v,
            None => {
                println!("couldn't parse seq from {:?}", s);
                return;
            }
        };
        let deg: usize = seq.iter().sum();
        match self.table.row_of(&seq) {
            Some(r) => println!("{} is at stem={}, row={}", format_seq(&seq), deg, r),
            None => println!("{} is NOT in the table (would be in stem={})", format_seq(&seq), deg),
        }
    }

    fn cmd_diffs(&self) {
        let mut any = false;
        for (k, diffs) in &self.table.differentials {
            for d in diffs {
                if is_zero_tail_artifact(&d.source) || is_zero_tail_artifact(&d.target) {
                    continue;
                }
                any = true;
                println!(
                    "  k={}  {} (row {}) → {} (row {})",
                    k,
                    format_seq(&d.source),
                    d.source_row,
                    format_seq(&d.target),
                    d.target_row
                );
            }
        }
        if !any {
            println!("  (none)");
        }
    }

    fn cmd_row(&self, s: &str) {
        let parts: Vec<&str> = s.split_whitespace().collect();
        if parts.len() != 2 {
            println!("usage: row <k> <n>");
            return;
        }
        let k: usize = match parts[0].parse() {
            Ok(v) => v,
            Err(_) => {
                println!("bad k: {:?}", parts[0]);
                return;
            }
        };
        let n: usize = match parts[1].parse() {
            Ok(v) => v,
            Err(_) => {
                println!("bad n: {:?}", parts[1]);
                return;
            }
        };
        match self.table.entries.get(&k).and_then(|col| col.get(&n)) {
            Some(entries) => {
                for seq in entries {
                    let flags = {
                        let mut f = String::new();
                        if self.table.is_source(seq) {
                            f.push_str(" [source]");
                        }
                        if self.table.is_target(seq) {
                            f.push_str(" [target]");
                        }
                        f
                    };
                    println!("  {}{}", format_seq(seq), flags);
                }
            }
            None => println!("  (empty)"),
        }
    }

    fn cmd_reset(&mut self) {
        self.current_elem = seq_to_element(&self.original_seq);
        self.iter_count = 0;
        println!("x reset to {}", format_seq(&self.original_seq));
    }
}
