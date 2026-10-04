//! Seeded random numbers: xoshiro256++ streams, one per shot, derived from (seed, shot) by
//! SplitMix64, so a batch's output never depends on how shots are spread over threads.

#[derive(Clone, Debug)]
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

#[derive(Clone, Debug)]
pub struct Xoshiro {
    s: [u64; 4],
}

impl Xoshiro {
    pub fn from_seed(seed: u64) -> Xoshiro {
        let mut sm = SplitMix64(seed);
        Xoshiro { s: [sm.next_u64(), sm.next_u64(), sm.next_u64(), sm.next_u64()] }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = (s[0].wrapping_add(s[3])).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Uniform in [0, 1) with 53 random bits.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform in (0, 1]: safe for logarithms.
    #[inline]
    pub fn next_f64_open0(&mut self) -> f64 {
        ((self.next_u64() >> 11) + 1) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// An exponential(1) variate.
    #[inline]
    pub fn next_exp(&mut self) -> f64 {
        -self.next_f64_open0().ln()
    }
}

/// The RNG stream of one shot.
pub fn shot_rng(seed: u64, shot: u64) -> Xoshiro {
    let mut sm = SplitMix64(seed);
    let a = sm.next_u64();
    let mut sm2 = SplitMix64(a ^ shot.wrapping_mul(0xD1B5_4A32_D192_ED03));
    Xoshiro::from_seed(sm2.next_u64())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn streams_differ_and_are_uniform() {
        let mut a = shot_rng(1, 0);
        let mut b = shot_rng(1, 1);
        assert_ne!(a.next_u64(), b.next_u64());
        let mut r = shot_rng(7, 3);
        let n = 100_000;
        let mean: f64 = (0..n).map(|_| r.next_f64()).sum::<f64>() / n as f64;
        assert!((mean - 0.5).abs() < 0.01, "{mean}");
        let e: f64 = (0..n).map(|_| r.next_exp()).sum::<f64>() / n as f64;
        assert!((e - 1.0).abs() < 0.02, "{e}");
    }
    #[test]
    fn deterministic() {
        assert_eq!(shot_rng(5, 9).next_u64(), shot_rng(5, 9).next_u64());
    }
}
