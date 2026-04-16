#![allow(dead_code)]

pub mod curtis;

use std::fmt;
use std::ops::{Add, Mul};
use smallvec::{smallvec, SmallVec};

/// Inline capacity for admissible sequences.  Most monomials encountered
/// through moderate stems have length ≤ 8; longer ones spill to heap.
pub type Seq = SmallVec<[usize; 8]>;

#[derive(PartialEq, Eq, Hash, Clone, Debug, PartialOrd, Ord)]
pub struct Admissible(pub Seq);

/// Ordered by sequence (lex), then degree.  Lex on Vec<usize> has the property
/// that the largest first generator gives the largest element — which is
/// exactly "highest filtration" in the EHP spectral sequence.
#[derive(PartialEq, Eq, Hash, Clone, Debug, PartialOrd, Ord)]
pub struct Monomial {
    pub seq: Admissible,
    pub deg: usize,
}

/// An element of Λ over F₂, stored as a sorted Vec of distinct monomials
/// (ascending by `seq`).  Invariant:
///   - No duplicates (F₂ coefficients mean two copies cancel).
///   - Strictly increasing order by monomial sequence.
///
/// This layout makes the common hot paths very fast:
///   - `filtration_leading` → O(1) via `.last()`.
///   - Addition → linear merge with cancellation.
///   - Iteration is cache-friendly; no hashing.
#[derive(PartialEq, Eq, Clone)]
pub struct Element(pub Vec<Monomial>);

impl Admissible {
    pub fn new(seq: &[usize]) -> Option<Self> {
        for (i, j) in seq.iter().zip(seq.iter().skip(1)) {
            if *j > 2 * i {
                return None;
            }
        }
        Some(Admissible(Seq::from_slice(seq)))
    }
}

impl Monomial {
    pub fn new(seq: &[usize]) -> Option<Self> {
        Some(Self {
            seq: Admissible::new(seq)?,
            deg: seq.iter().sum(),
        })
    }
}

impl Element {
    pub fn new(seq: &[usize]) -> Option<Self> {
        Some(Self(vec![Monomial::new(seq)?]))
    }

    pub fn zero() -> Self {
        Self(Vec::new())
    }

    pub fn singleton(mono: Monomial) -> Self {
        Self(vec![mono])
    }

    /// Toggle `mono` in/out of the sum (F₂ semantics) while preserving sort.
    pub fn add_mono(&mut self, mono: Monomial) {
        match self.0.binary_search(&mono) {
            Ok(idx) => { self.0.remove(idx); }       // cancels (mod 2)
            Err(idx) => { self.0.insert(idx, mono); }
        }
    }

    pub fn diff(self) -> Self {
        let mut result = Element::zero();
        for mono in self.0 {
            result = result + diff_mono(&mono.seq.0);
        }
        result
    }
}

// d(λ_i) = Σ_{j=1}^{⌊i/2⌋} C(i-j, j) · λ_{i-j} · λ_{j-1}
// All terms are admissible: j-1 ≤ 2(i-j) follows from j ≤ i/2.
pub fn diff_generator(i: usize) -> Element {
    let mut result = Element::zero();
    for j in 1..=(i / 2) {
        if binom_mod2(i - j, j) {
            let term = Element::new(&[i - j, j - 1]).unwrap();
            result = result + term;
        }
    }
    result
}

// Iterative Leibniz:
//     d(a_1 … a_n) = Σ_i  a_1 … a_{i-1} · d(a_i) · a_{i+1} … a_n
// This avoids building intermediate Elements via recursive splitting and
// keeps prefix/suffix as singletons, hitting the fast path in Element::mul.
pub fn diff_mono(seq: &[usize]) -> Element {
    if seq.is_empty() {
        return Element::zero();
    }
    if seq.len() == 1 {
        return diff_generator(seq[0]);
    }
    let mut result = Element::zero();
    for i in 0..seq.len() {
        let d_gen = diff_generator(seq[i]);
        if d_gen.0.is_empty() {
            continue;
        }
        let prefix = &seq[..i];
        let suffix = &seq[i + 1..];
        let prefix_elem = Element::singleton(Monomial {
            seq: Admissible(Seq::from_slice(prefix)),
            deg: prefix.iter().sum(),
        });
        let suffix_elem = Element::singleton(Monomial {
            seq: Admissible(Seq::from_slice(suffix)),
            deg: suffix.iter().sum(),
        });
        result = result + prefix_elem * d_gen * suffix_elem;
    }
    result
}

impl Add for Element {
    type Output = Element;

    /// Merge two sorted monomial lists with F₂ cancellation — O(n + m).
    fn add(self, rhs: Self) -> Self::Output {
        let Element(a) = self;
        let Element(b) = rhs;
        let mut out = Vec::with_capacity(a.len() + b.len());
        let mut ai = a.into_iter();
        let mut bi = b.into_iter();
        let mut a_cur = ai.next();
        let mut b_cur = bi.next();
        loop {
            match (a_cur.take(), b_cur.take()) {
                (None, None) => break,
                (Some(x), None) => {
                    out.push(x);
                    out.extend(ai);
                    break;
                }
                (None, Some(y)) => {
                    out.push(y);
                    out.extend(bi);
                    break;
                }
                (Some(x), Some(y)) => {
                    use std::cmp::Ordering::*;
                    match x.cmp(&y) {
                        Less => {
                            out.push(x);
                            a_cur = ai.next();
                            b_cur = Some(y);
                        }
                        Greater => {
                            out.push(y);
                            a_cur = Some(x);
                            b_cur = bi.next();
                        }
                        Equal => {
                            // cancel (mod 2)
                            a_cur = ai.next();
                            b_cur = bi.next();
                        }
                    }
                }
            }
        }
        Element(out)
    }
}

impl Monomial {
    /// Reference-taking multiplication — avoids redundant clones when the
    /// same operand appears in many pairwise products (the hot Element::Mul loop).
    pub fn mul_ref(&self, rhs: &Monomial) -> Element {
        let admissible_junction = match (self.seq.0.last(), rhs.seq.0.first()) {
            (Some(&l), Some(&r)) => r <= 2 * l,
            _ => true,
        };
        if admissible_junction {
            let mut seq: Seq = SmallVec::with_capacity(self.seq.0.len() + rhs.seq.0.len());
            seq.extend_from_slice(&self.seq.0);
            seq.extend_from_slice(&rhs.seq.0);
            Element::singleton(Monomial {
                deg: self.deg + rhs.deg,
                seq: Admissible(seq),
            })
        } else {
            let z = *self.seq.0.last().unwrap();
            let v = *rhs.seq.0.first().unwrap();
            let left_prefix = Monomial {
                seq: Admissible(Seq::from_slice(&self.seq.0[..self.seq.0.len() - 1])),
                deg: self.deg - z,
            };
            let right_suffix = Monomial {
                seq: Admissible(Seq::from_slice(&rhs.seq.0[1..])),
                deg: rhs.deg - v,
            };
            let middle = adem(z, v);
            Element::singleton(left_prefix) * middle * Element::singleton(right_suffix)
        }
    }
}

impl Mul for Monomial {
    type Output = Element;
    fn mul(self, rhs: Self) -> Self::Output { self.mul_ref(&rhs) }
}

pub fn binom_mod2(n: usize, k: usize) -> bool {
    n & k == k
}

pub fn adem(n: usize, b: usize) -> Element {
    let i = b - 2 * n - 1;
    let mut result = Element::zero();
    if i == 0 {
        return result;
    }
    for j in 0..=(i - 1) / 2 {
        if binom_mod2(i - j - 1, j) {
            let left = Monomial { seq: Admissible(smallvec![n + i - j]), deg: n + i - j };
            let right = Monomial { seq: Admissible(smallvec![2 * n + j + 1]), deg: 2 * n + j + 1 };
            result = result + left * right;
        }
    }
    result
}

impl Mul for Element {
    type Output = Element;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.0.is_empty() || rhs.0.is_empty() {
            return Element::zero();
        }
        // Fast path: both sides are singletons (very common via recursive
        // Monomial multiplication) — skip bucket sort overhead.
        if self.0.len() == 1 && rhs.0.len() == 1 {
            return self.0[0].mul_ref(&rhs.0[0]);
        }
        // General path: collect products, sort, cancel F₂ runs in place.
        let mut bucket: Vec<Monomial> = Vec::with_capacity(self.0.len() * rhs.0.len());
        for l in &self.0 {
            for r in &rhs.0 {
                bucket.extend(l.mul_ref(r).0);
            }
        }
        bucket.sort();
        // In-place compaction: keep odd-length runs, drop even.
        let n = bucket.len();
        let mut write = 0;
        let mut i = 0;
        while i < n {
            let mut j = i + 1;
            while j < n && bucket[j] == bucket[i] { j += 1; }
            if (j - i) & 1 == 1 {
                if write != i { bucket.swap(write, i); }
                write += 1;
            }
            i = j;
        }
        bucket.truncate(write);
        Element(bucket)
    }
}

// ── Display ──────────────────────────────────────────────────────────────────

impl fmt::Display for Monomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.seq.0.is_empty() {
            return write!(f, "1");
        }
        let inner: Vec<String> = self.seq.0.iter().map(|n| n.to_string()).collect();
        write!(f, "λ_({})", inner.join(", "))
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            return write!(f, "0");
        }
        // Element is already sorted by seq; reorder by degree for nicer printing.
        let mut terms: Vec<&Monomial> = self.0.iter().collect();
        terms.sort_by(|a, b| a.deg.cmp(&b.deg).then_with(|| a.seq.0.cmp(&b.seq.0)));
        let strs: Vec<String> = terms.iter().map(|m| m.to_string()).collect();
        write!(f, "{}", strs.join(" + "))
    }
}
