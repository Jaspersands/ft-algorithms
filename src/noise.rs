//! Static noise sites and how faults are drawn.
//!
//! A *site* is one target group of a noise instruction outside every TABLE, counted in program
//! order with REPEAT bodies repeated: site i fires with probability q_i, independently. The tree
//! mirrors the program's REPEAT nesting, so a body executed a million times is stored once.
//!
//! - Plain sampling walks a unit-rate Poisson process along the cumulative hazard
//!   h_i = −ln(1 − q_i): site i fires iff a point lands in its interval, probability 1 − e^{−h_i} = q_i.
//! - Sampling *exactly k* faults draws k sites i.i.d. with weight w_i = q_i / (1 − q_i) and
//!   redraws if two coincide; the accepted sets have probability ∝ Π w_i, which is the law of the
//!   fired set given that exactly k sites fired.
//! - P(exactly k faults) is the coefficient of z^k in Π (1 − q_i + q_i z).

use crate::program::{Block, Op, Program};
use crate::rng::Xoshiro;

#[derive(Clone, Debug)]
enum Item {
    Site,
    Repeat(u64, Box<Node>),
}

#[derive(Clone, Debug, Default)]
struct Node {
    items: Vec<Item>,
    /// Inclusive prefix sums over items: sites, weights, hazards.
    cum_count: Vec<u64>,
    cum_w: Vec<f64>,
    cum_h: Vec<f64>,
    /// Per-site q (for Site items; 0 for Repeat items), used by the fault-count polynomial.
    q: Vec<f64>,
}

impl Node {
    fn count(&self) -> u64 {
        *self.cum_count.last().unwrap_or(&0)
    }
    fn weight(&self) -> f64 {
        *self.cum_w.last().unwrap_or(&0.0)
    }
    fn hazard(&self) -> f64 {
        *self.cum_h.last().unwrap_or(&0.0)
    }

    fn push(&mut self, item: Item, count: u64, w: f64, h: f64, q: f64) {
        let (c0, w0, h0) = (self.count(), self.weight(), self.hazard());
        self.items.push(item);
        self.cum_count.push(c0.saturating_add(count));
        self.cum_w.push(w0 + w);
        self.cum_h.push(h0 + h);
        self.q.push(q);
    }

    fn build(block: &Block) -> Node {
        let mut node = Node::default();
        for op in block {
            match op {
                Op::Noise(n, qs) => {
                    let q = n.total();
                    let w = q / (1.0 - q);
                    let h = -(1.0 - q).ln();
                    for _ in 0..qs.len() / n.arity() {
                        node.push(Item::Site, 1, w, h, q);
                    }
                }
                Op::Repeat(k, body) => {
                    let b = Node::build(body);
                    if b.count() > 0 {
                        let (c, w, h) = (b.count().saturating_mul(*k), b.weight() * *k as f64, b.hazard() * *k as f64);
                        node.push(Item::Repeat(*k, Box::new(b)), c, w, h, 0.0);
                    }
                }
                // Noise inside TABLE cases is dynamic: not a static site.
                _ => {}
            }
        }
        node
    }

    /// The site at weight coordinate x ∈ [0, weight()), as a global index from `base`.
    fn find_weight(&self, mut x: f64, base: u64) -> u64 {
        let i = self.cum_w.partition_point(|&c| c <= x).min(self.items.len() - 1);
        let prev_w = if i == 0 { 0.0 } else { self.cum_w[i - 1] };
        let prev_c = if i == 0 { 0 } else { self.cum_count[i - 1] };
        x -= prev_w;
        match &self.items[i] {
            Item::Site => base + prev_c,
            Item::Repeat(k, body) => {
                let bw = body.weight();
                let it = ((x / bw) as u64).min(k - 1);
                body.find_weight((x - it as f64 * bw).max(0.0), base + prev_c + it * body.count())
            }
        }
    }

    /// The first site whose cumulative hazard (inclusive) exceeds t, with that inclusive
    /// cumulative hazard; None if t ≥ hazard().
    fn find_hazard(&self, t: f64, base: u64, base_h: f64) -> Option<(u64, f64)> {
        if !(t < self.hazard()) {
            return None;
        }
        let i = self.cum_h.partition_point(|&c| c <= t).min(self.items.len() - 1);
        let prev_h = if i == 0 { 0.0 } else { self.cum_h[i - 1] };
        let prev_c = if i == 0 { 0 } else { self.cum_count[i - 1] };
        match &self.items[i] {
            Item::Site => Some((base + prev_c, base_h + self.cum_h[i])),
            Item::Repeat(k, body) => {
                let bh = body.hazard();
                let x = t - prev_h;
                let mut it = ((x / bh) as u64).min(k - 1);
                loop {
                    let r = body.find_hazard(x - it as f64 * bh, base + prev_c + it * body.count(), base_h + prev_h + it as f64 * bh);
                    match r {
                        Some(v) => return Some(v),
                        // Rounding put x at the very end of iteration `it`: the site is in the next one.
                        None if it + 1 < *k => it += 1,
                        None => return Some((base + prev_c + k * body.count() - 1, base_h + self.cum_h[i])),
                    }
                }
            }
        }
    }

    /// Π over sites of (1 − q + q z), truncated to degree `deg`.
    fn poly(&self, deg: usize) -> Vec<f64> {
        let mut p = vec![0.0; deg + 1];
        p[0] = 1.0;
        for (item, &q) in self.items.iter().zip(&self.q) {
            match item {
                Item::Site => {
                    for j in (0..=deg).rev() {
                        let below = if j > 0 { p[j - 1] } else { 0.0 };
                        p[j] = p[j] * (1.0 - q) + below * q;
                    }
                }
                Item::Repeat(k, body) => {
                    let b = poly_pow(&body.poly(deg), *k, deg);
                    p = poly_mul(&p, &b, deg);
                }
            }
        }
        p
    }
}

fn poly_mul(a: &[f64], b: &[f64], deg: usize) -> Vec<f64> {
    let mut c = vec![0.0; deg + 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0.0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            if i + j > deg {
                break;
            }
            c[i + j] += x * y;
        }
    }
    c
}

fn poly_pow(base: &[f64], mut e: u64, deg: usize) -> Vec<f64> {
    let mut result = vec![0.0; deg + 1];
    result[0] = 1.0;
    let mut b = base.to_vec();
    while e > 0 {
        if e & 1 == 1 {
            result = poly_mul(&result, &b, deg);
        }
        e >>= 1;
        if e > 0 {
            b = poly_mul(&b, &b, deg);
        }
    }
    result
}

/// Which faults a shot gets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FaultPlan {
    /// No noise at all (static or dynamic).
    None,
    /// Every site fires independently with its probability.
    Plain,
    /// Exactly k static sites fire (drawn from the conditional law); dynamic sites fire
    /// independently.
    Exactly(usize),
}

#[derive(Clone, Debug)]
pub struct SiteTree {
    root: Node,
}

/// Redraw attempts before `sample_exactly` gives up (only reachable when k is close to the
/// number of sites with non-zero probability).
const MAX_REDRAWS: usize = 1_000_000;

impl SiteTree {
    pub fn build(p: &Program) -> SiteTree {
        SiteTree { root: Node::build(&p.body) }
    }

    pub fn num_sites(&self) -> u64 {
        self.root.count()
    }

    /// Expected number of faults Σ q_i.
    pub fn expected_faults(&self) -> f64 {
        fn mean(n: &Node) -> f64 {
            n.items
                .iter()
                .zip(&n.q)
                .map(|(it, &q)| match it {
                    Item::Site => q,
                    Item::Repeat(k, b) => *k as f64 * mean(b),
                })
                .sum()
        }
        mean(&self.root)
    }

    /// [P(K = 0), …, P(K = kmax), P(K > kmax)].
    pub fn fault_count_distribution(&self, kmax: usize) -> Vec<f64> {
        let mut p = self.root.poly(kmax);
        let s: f64 = p.iter().sum();
        p.push((1.0 - s).max(0.0));
        p
    }

    /// Sorted indices of the sites that fire, each independently with its probability.
    pub fn sample_plain(&self, rng: &mut Xoshiro) -> Vec<u64> {
        let mut out = Vec::new();
        let mut t = rng.next_exp();
        while let Some((site, h_end)) = self.root.find_hazard(t, 0, 0.0) {
            out.push(site);
            t = h_end + rng.next_exp();
        }
        out
    }

    /// Sorted, distinct indices of exactly k sites, from the law of the fired set given K = k.
    pub fn sample_exactly(&self, k: usize, rng: &mut Xoshiro) -> Result<Vec<u64>, String> {
        if k == 0 {
            return Ok(vec![]);
        }
        let w = self.root.weight();
        if w <= 0.0 || (k as u64) > self.num_sites() {
            return Err(format!("cannot place {k} faults: the program has {} noise sites", self.num_sites()));
        }
        for _ in 0..MAX_REDRAWS {
            let mut s: Vec<u64> = (0..k).map(|_| self.root.find_weight(rng.next_f64() * w, 0)).collect();
            s.sort_unstable();
            if s.windows(2).all(|p| p[0] != p[1]) {
                return Ok(s);
            }
        }
        Err(format!("could not draw {k} distinct faults (too few sites with non-zero probability)"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::shot_rng;

    fn tree(text: &str) -> SiteTree {
        SiteTree::build(&Program::parse(text).unwrap())
    }

    fn binom(n: u64, k: u64, q: f64) -> f64 {
        let mut c = 1.0;
        for i in 0..k {
            c *= (n - i) as f64 / (i + 1) as f64;
        }
        c * q.powi(k as i32) * (1.0 - q).powi((n - k) as i32)
    }

    #[test]
    fn counts_sites() {
        assert_eq!(tree("REPEAT 4 {\nDEPOLARIZE1(0.1) 0 1\n}").num_sites(), 8);
        assert_eq!(tree("DEPOLARIZE2(0.1) 0 1 2 3\nZ_CHANNEL_3(0.01,0,0,0,0,0,0) 0 1 2").num_sites(), 3);
        // Noise inside TABLE cases is dynamic.
        assert_eq!(tree("M 0\nTABLE rec[-1] {X_ERROR(0.1) 0} {X_ERROR(0.1) 0}\nX_ERROR(0.1) 1").num_sites(), 1);
    }

    #[test]
    fn fault_count_distribution_is_binomial() {
        let t = tree("REPEAT 2 {\nREPEAT 2 {\nDEPOLARIZE1(0.1) 0 1\n}\n}");
        let p = t.fault_count_distribution(8);
        for k in 0..=8 {
            assert!((p[k] - binom(8, k as u64, 0.1)).abs() < 1e-12, "k={k}");
        }
        assert!(p[9].abs() < 1e-12);
        let p = t.fault_count_distribution(2);
        let tail: f64 = (3..=8).map(|k| binom(8, k, 0.1)).sum();
        assert!((p[3] - tail).abs() < 1e-12);
    }

    #[test]
    fn mixed_probabilities_match_direct_dp() {
        let t = tree("X_ERROR(0.1) 0\nREPEAT 3 {\nZ_ERROR(0.02) 0\nX_ERROR(0.3) 1\n}\nY_ERROR(0.5) 2");
        let qs = [0.1, 0.02, 0.3, 0.02, 0.3, 0.02, 0.3, 0.5];
        let mut dp = vec![1.0];
        for q in qs {
            let mut nd = vec![0.0; dp.len() + 1];
            for (j, &v) in dp.iter().enumerate() {
                nd[j] += v * (1.0 - q);
                nd[j + 1] += v * q;
            }
            dp = nd;
        }
        let p = t.fault_count_distribution(8);
        for k in 0..=8 {
            assert!((p[k] - dp[k]).abs() < 1e-12);
        }
        let mean: f64 = qs.iter().sum();
        assert!((t.expected_faults() - mean).abs() < 1e-12);
    }

    #[test]
    fn plain_sampling_per_site_frequencies() {
        // 3 + 2·(2) = 7 sites with mixed q.
        let t = tree("X_ERROR(0.05) 0 1 2\nREPEAT 2 {\nX_ERROR(0.2) 0\nX_ERROR(0.01) 1\n}");
        let qs = [0.05, 0.05, 0.05, 0.2, 0.01, 0.2, 0.01];
        let shots = 400_000u64;
        let mut hits = [0u64; 7];
        for s in 0..shots {
            for i in t.sample_plain(&mut shot_rng(3, s)) {
                hits[i as usize] += 1;
            }
        }
        for i in 0..7 {
            let m = shots as f64 * qs[i];
            let sd = (m * (1.0 - qs[i])).sqrt();
            assert!((hits[i] as f64 - m).abs() < 5.0 * sd, "site {i}: {} vs {m}", hits[i]);
        }
    }

    #[test]
    fn plain_sampling_on_a_long_repeat() {
        let t = tree("REPEAT 1000 {\nX_ERROR(0.001) 0\n}");
        let shots = 100_000u64;
        let mut total = 0u64;
        let mut first_half = 0u64;
        for s in 0..shots {
            let f = t.sample_plain(&mut shot_rng(4, s));
            assert!(f.windows(2).all(|w| w[0] < w[1]));
            total += f.len() as u64;
            first_half += f.iter().filter(|&&i| i < 500).count() as u64;
        }
        let mean = total as f64 / shots as f64;
        assert!((mean - 1.0).abs() < 5.0 * (1.0 / shots as f64).sqrt(), "{mean}");
        let frac = first_half as f64 / total as f64;
        assert!((frac - 0.5).abs() < 0.01, "{frac}");
    }

    #[test]
    fn exactly_k_follows_the_conditional_law() {
        // Sites with weights w = q/(1−q); pairs should appear ∝ w_i w_j.
        let t = tree("X_ERROR(0.1) 0\nX_ERROR(0.2) 0\nREPEAT 2 {\nX_ERROR(0.3) 0\n}");
        let q = [0.1, 0.2, 0.3, 0.3];
        let w: Vec<f64> = q.iter().map(|q| q / (1.0 - q)).collect();
        let mut expect = vec![];
        for i in 0..4 {
            for j in i + 1..4 {
                expect.push(((i, j), w[i] * w[j]));
            }
        }
        let z: f64 = expect.iter().map(|e| e.1).sum();
        let shots = 300_000u64;
        let mut counts = std::collections::HashMap::new();
        for s in 0..shots {
            let f = t.sample_exactly(2, &mut shot_rng(5, s)).unwrap();
            *counts.entry((f[0] as usize, f[1] as usize)).or_insert(0u64) += 1;
        }
        for ((i, j), wij) in expect {
            let m = shots as f64 * wij / z;
            let got = *counts.get(&(i, j)).unwrap_or(&0) as f64;
            assert!((got - m).abs() < 5.0 * m.sqrt(), "pair ({i},{j}): {got} vs {m}");
        }
        assert!(t.sample_exactly(5, &mut shot_rng(1, 1)).is_err());
    }
}
