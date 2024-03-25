//! Provide methods to support modular arithmetic.
//!
//! Implemented features is still few and the API is subject to significant change.

mod barrett;
mod montgomery;

pub use barrett::*;
pub use montgomery::*;

/// Calculate `a^exp (mod m)`.
///
/// # Panics
/// - `m > 0` must be satisfied.
///
/// # Examples
/// ```rust
/// use modutil::pow_mod;
///
/// assert_eq!(pow_mod(2, 10, 2048), 1024);
/// assert_eq!(pow_mod(2, 10, 512), 0);
/// assert_eq!(pow_mod(3, 4, 7), 4);
/// assert_eq!(pow_mod(3, 4, 2), 1);
/// assert_eq!(pow_mod(3, 4, 14), 11);
/// ```
pub const fn pow_mod(a: u64, mut exp: u64, m: u64) -> u64 {
    assert!(m > 0);
    if m == 1 {
        return 0;
    }
    if exp == 0 {
        return 1;
    }

    let a = a % m;
    if exp == 1 {
        return a;
    }
    if a == 1 {
        return 1;
    }

    // If m is even, Montgomery multiplication cannot be used.
    // In this case, m=m' * 2^k (m' is odd), and the CRT restoration can be done as `pow_mod(a,exp,m')` and `pow_mod(a,exp,2^k)`.
    if m & 1 == 0 {
        let s = m.trailing_zeros();
        let r = pow_mod(a, exp, m >> s);
        let (mut res, mut val) = (1u64, a);
        while exp > 0 {
            if exp & 1 != 0 {
                res = res.wrapping_mul(val);
            }
            val = val.wrapping_mul(val);
            exp >>= 1;
        }
        res &= (1 << s) - 1;

        if m == 1 << s {
            return res;
        }

        let Some(res) = crt(r, m >> s, res, 1 << s) else {
            unreachable!()
        };
        return res;
    }

    let mont = Montgomery::<u64>::new(m);
    mont.reduce(mont.pow(mont.convert(a), exp))
}

/// Find `x` that satisfies `ax = 1 (mod m)`.  
/// If not found, return `None`.
///
/// # Panics
/// - `m > 0` must be satisfied.
///
/// # Examples
/// ```rust
/// use modutil::inverse_mod;
///
/// assert_eq!(inverse_mod(3, 8), Some(3));
/// assert_eq!(inverse_mod(2, 8), None);
/// assert_eq!(inverse_mod(4, 17), Some(13));
/// assert_eq!(inverse_mod(0, 17), None);
/// ```
pub const fn inverse_mod(a: u64, m: u64) -> Option<u64> {
    assert!(m > 0);
    if a == 0 {
        return None;
    }
    let (mut s, mut sy) = (m, 0u64);
    let (mut t, mut ty) = (a, 1u64);
    while t != 0 {
        let d = s / t;
        let u = s % t;

        let (v, f) = ty.overflowing_mul(d);
        let uy = if f || v >= m {
            sy.wrapping_add(ty.wrapping_neg().wrapping_mul(d))
        } else {
            sy.wrapping_sub(v)
        };

        (s, sy, t, ty) = (t, ty, u, uy);
    }

    if sy > m {
        sy = sy.wrapping_add(m);
    }

    if s == 1 {
        Some(sy)
    } else {
        None
    }
}

/// If x satisfies x = p (mod m) and x = q (mod n) is found, return Some(x).  
/// Otherwise, return None.
/// ```ignore
/// p + ma = x
/// q + nb = x
///     => p + ma = q + nb
///     => ma - nb = q - p
///     => ma = q - p (mod n)
/// ```
///
/// # Panics
/// - Both p < m and q < n must be satisfied.
///
/// # Examples
/// ```rust
/// use modutil::crt;
///
/// assert_eq!(crt(2, 3, 3, 5), Some(8));
/// assert_eq!(crt(3, 5, 2, 7), Some(23));
/// assert_eq!(crt(2, 7, 2, 3), Some(2));
/// assert_eq!(crt(2, 3, 23, 35), Some(23));
/// ```
pub const fn crt(p: u64, m: u64, q: u64, n: u64) -> Option<u64> {
    assert!(p < m && q < n);

    if q < p {
        return crt(q, n, p, m);
    }

    if p == q {
        return Some(p);
    }

    let w = q - p;
    let (mut s, mut ys) = (n, 0u64);
    let (mut t, mut yt) = (m, 1u64);
    while s % t != 0 {
        let tmp = s / t;
        let u = s % t;

        let (v, f) = yt.overflowing_mul(tmp);
        let yu = if f || v >= n {
            ys.wrapping_add(yt.wrapping_neg().wrapping_mul(tmp))
        } else {
            ys.wrapping_sub(v)
        };

        (s, ys, t, yt) = (t, yt, u, yu);
    }

    if yt > n {
        yt = yt.wrapping_add(n);
    }

    if w % t == 0 {
        let g = w / t;
        let lcm = m / t * n;
        let h = (g * yt) % (n / t);
        let res = (p + m * h) % lcm;
        debug_assert!(res % m == p);
        debug_assert!(res % n == q);
        Some(res)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crt_test() {
        assert_eq!(crt(2, 3, 3, 5), Some(8));
        assert_eq!(crt(3, 5, 2, 7), Some(23));
        assert_eq!(crt(2, 7, 2, 3), Some(2));
        assert_eq!(crt(2, 3, 23, 35), Some(23));
        assert_eq!(crt(4, 7, 1, 2), Some(11));
    }

    #[test]
    fn pow_mod_test() {
        assert_eq!(pow_mod(2, 10, 998244353), 1024);
        assert_eq!(pow_mod(2, 32, 998244353), 2u64.pow(32) % 998244353);
        assert_eq!(pow_mod(2, 10, 2048), 1024);
        assert_eq!(pow_mod(2, 10, 8), 0);
        assert_eq!(pow_mod(3, 4, 2), 1);
        assert_eq!(pow_mod(3, 4, 7), 4);
        assert_eq!(pow_mod(3, 4, 14), 11);
    }

    #[test]
    fn inverse_mod_test() {
        for m in 2..1000 {
            for i in 1..m {
                if let Some(res) = inverse_mod(i, m) {
                    assert_eq!(res * i % m, 1);
                }
            }
        }
    }
}
