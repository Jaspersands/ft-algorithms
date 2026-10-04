//! A minimal complex number, so the crate has no numeric dependencies.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

impl C64 {
    pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };
    pub const ONE: C64 = C64 { re: 1.0, im: 0.0 };
    pub const I: C64 = C64 { re: 0.0, im: 1.0 };

    #[inline]
    pub const fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }
    #[inline]
    pub fn from_polar(r: f64, theta: f64) -> C64 {
        C64 { re: r * theta.cos(), im: r * theta.sin() }
    }
    #[inline]
    pub fn conj(self) -> C64 {
        C64 { re: self.re, im: -self.im }
    }
    #[inline]
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    #[inline]
    pub fn scale(self, s: f64) -> C64 {
        C64 { re: self.re * s, im: self.im * s }
    }
    /// Multiplication by i.
    #[inline]
    pub fn mul_i(self) -> C64 {
        C64 { re: -self.im, im: self.re }
    }
}

impl Add for C64 {
    type Output = C64;
    #[inline]
    fn add(self, o: C64) -> C64 {
        C64 { re: self.re + o.re, im: self.im + o.im }
    }
}
impl AddAssign for C64 {
    #[inline]
    fn add_assign(&mut self, o: C64) {
        self.re += o.re;
        self.im += o.im;
    }
}
impl Sub for C64 {
    type Output = C64;
    #[inline]
    fn sub(self, o: C64) -> C64 {
        C64 { re: self.re - o.re, im: self.im - o.im }
    }
}
impl Mul for C64 {
    type Output = C64;
    #[inline]
    fn mul(self, o: C64) -> C64 {
        C64 { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
}
impl MulAssign for C64 {
    #[inline]
    fn mul_assign(&mut self, o: C64) {
        *self = *self * o;
    }
}
impl Neg for C64 {
    type Output = C64;
    #[inline]
    fn neg(self) -> C64 {
        C64 { re: -self.re, im: -self.im }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn arithmetic() {
        let a = C64::new(1.0, 2.0);
        let b = C64::new(-3.0, 0.5);
        assert_eq!(a * b, C64::new(-3.0 - 1.0, 0.5 - 6.0));
        assert_eq!(a.mul_i(), a * C64::I);
        assert_eq!((a * a.conj()).im, 0.0);
        assert!((C64::from_polar(2.0, std::f64::consts::FRAC_PI_2).im - 2.0).abs() < 1e-15);
    }
}
