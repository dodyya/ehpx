#![allow(dead_code)]

pub mod curtis;

use std::fmt;
use std::ops::{Add, Mul};
use std::collections::HashSet;

#[derive(PartialEq, Eq, Hash, Clone, Debug, PartialOrd, Ord)]
pub struct Admissible(pub Vec<usize>);

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct Monomial {
    pub seq: Admissible,
    pub deg: usize,
}

#[derive(PartialEq, Eq, Clone)]
pub struct Element(pub HashSet<Monomial>);

impl Admissible {
    pub fn new(seq: &[usize]) -> Option<Self> {
        for (i, j) in seq.iter().zip(seq.iter().skip(1)) {
            if *j > 2 * i {
                return None;
            }
        }
        Some(Admissible(seq.to_vec()))
    }
}

impl Monomial {
    pub fn new(seq: &[usize]) -> Option<Self> {
        Some(Self {
            seq: Admissible::new(seq)?,
            deg: seq.into_iter().sum(),
        })
    }
}

impl Element {
    pub fn new(seq: &[usize]) -> Option<Self> {
        let mono = Monomial::new(seq)?;
        let mut set = HashSet::new();
        set.insert(mono);
        Some(Self(set))
    }

    pub fn zero() -> Self {
        Self(HashSet::new())
    }

    pub fn add_mono(&mut self, mono: Monomial) {
        if !self.0.remove(&mono) {
            self.0.insert(mono);
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

// Leibniz: d(xy) = x*dy + dx*y, splitting seq near the middle.
pub fn diff_mono(seq: &[usize]) -> Element {
    match seq.len() {
        0 => Element::zero(),
        1 => diff_generator(seq[0]),
        n => {
            let mid = n / 2;
            let (x_seq, y_seq) = seq.split_at(mid);
            let x = Element(HashSet::from([Monomial {
                seq: Admissible(x_seq.to_vec()),
                deg: x_seq.iter().sum(),
            }]));
            let y = Element(HashSet::from([Monomial {
                seq: Admissible(y_seq.to_vec()),
                deg: y_seq.iter().sum(),
            }]));
            let dx = diff_mono(x_seq);
            let dy = diff_mono(y_seq);
            x.clone() * dy + dx * y.clone()
        }
    }
}

impl Add for Element {
    type Output = Element;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        for item in rhs.0 {
            result.add_mono(item);
        }
        result
    }
}

impl Mul for Monomial {
    type Output = Element;

    fn mul(self, rhs: Self) -> Self::Output {
        let admissible_junction = match (self.seq.0.last(), rhs.seq.0.first()) {
            (Some(&l), Some(&r)) => r <= 2 * l,
            _ => true,
        };
        if admissible_junction {
            let seq: Vec<usize> = self.seq.0.iter().chain(&rhs.seq.0).copied().collect();
            let mono = Monomial { deg: self.deg + rhs.deg, seq: Admissible(seq) };
            let mut set = HashSet::new();
            set.insert(mono);
            Element(set)
        } else {
            let z = *self.seq.0.last().unwrap();
            let v = *rhs.seq.0.first().unwrap();
            let left_prefix = Monomial {
                seq: Admissible(self.seq.0[..self.seq.0.len() - 1].to_vec()),
                deg: self.deg - z,
            };
            let right_suffix = Monomial {
                seq: Admissible(rhs.seq.0[1..].to_vec()),
                deg: rhs.deg - v,
            };
            let middle = adem(z, v);
            Element(HashSet::from([left_prefix])) * middle * Element(HashSet::from([right_suffix]))
        }
    }
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
            let left = Monomial { seq: Admissible(vec![n + i - j]), deg: n + i - j };
            let right = Monomial { seq: Admissible(vec![2 * n + j + 1]), deg: 2 * n + j + 1 };
            result = result + left * right;
        }
    }
    result
}

impl Mul for Element {
    type Output = Element;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = Element::zero();
        for l in &self.0 {
            for r in &rhs.0 {
                result = result + l.clone() * r.clone();
            }
        }
        result
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
        let mut terms: Vec<&Monomial> = self.0.iter().collect();
        terms.sort_by(|a, b| a.deg.cmp(&b.deg).then_with(|| a.seq.0.cmp(&b.seq.0)));
        let strs: Vec<String> = terms.iter().map(|m| m.to_string()).collect();
        write!(f, "{}", strs.join(" + "))
    }
}
