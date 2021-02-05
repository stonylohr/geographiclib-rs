#![allow(non_snake_case)]

pub const DIGITS: u64 = 53;
pub const TWO: f64 = 2.0;

pub fn get_epsilon() -> f64 {
    TWO.powi(1 - DIGITS as i32)
}

pub fn get_min_val() -> f64 {
    TWO.powi(-1022)
}

// Square
pub fn sq(x: f64) -> f64 {
    x.powi(2)
}

// We use the built-in impl (f64::cbrt) rather than this.
// Real cube root
pub fn cbrt(x: f64) -> f64 {
    // y = math.pow(abs(x), 1/3.0)
    let y = x.abs().powf(1.0 / 3.0);

    // return y if x > 0 else (-y if x < 0 else x)
    if x > 0.0 {
        y
    } else if x < 0.0 {
        -y
    } else {
        x
    }
}

// Normalize a two-vector
pub fn norm(x: &mut f64, y: &mut f64) {
    let r = x.hypot(*y);
    *x /= r;
    *y /= r;
}

// Error free transformation of a sum
pub fn sum(u: f64, v: f64) -> (f64, f64) {
    let s = u + v;
    let up = s - v;
    let vpp = s - up;
    let up = up - u;
    let vpp = vpp - v;
    let t = -(up + vpp);
    (s, t)
}

// Evaluate a polynomial.
pub fn polyval(n: isize, p: &[f64], x: f64) -> f64 {
    if n < 0 {
        0.0
    } else {
        let mut y = p[0];
        for val in &p[1..=n as usize] {
            y = y * x + val;
        }
        y
    }
}

// Round an angle so taht small values underflow to 0
pub fn ang_round(x: f64) -> f64 {
    // The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
    // for reals = 0.7 pm on the earth if x is an angle in degrees.  (This
    // is about 1000 times more resolution than we get with angles around 90
    // degrees.)  We use this to avoid having to deal with near singular
    // cases when x is non-zero but tiny (e.g., 1.0e-200).
    let z = 1.0 / 16.0;
    let mut y = x.abs();
    // The compiler mustn't "simplify" z - (z - y) to y
    if y < z {
        y = z - (z - y);
    };
    if x == 0.0 {
        0.0
    } else {
        if x < 0.0 {
            -y
        } else {
            y
        }
    }
}

/// remainder of x/y in the range [-y/2, y/2]
fn remainder(x: f64, y: f64) -> f64 {
    // z = math.fmod(x, y) if Math.isfinite(x) else Math.nan
    let z = if x.is_finite() { x % y } else { std::f64::NAN };

    // # On Windows 32-bit with python 2.7, math.fmod(-0.0, 360) = +0.0
    // # This fixes this bug.  See also Math::AngNormalize in the C++ library.
    // # sincosd has a similar fix.
    // z = x if x == 0 else z
    let z = if x == 0.0 { x } else { z };

    // return (z + y if z < -y/2 else
    // (z if z < y/2 else z -y))
    if z < -y / 2.0 {
        z + y
    } else {
        if z < y / 2.0 {
            z
        } else {
            z - y
        }
    }
}

/// reduce angle to (-180,180]
pub fn ang_normalize(x: f64) -> f64 {
    // y = Math.remainder(x, 360)
    // return 180 if y == -180 else y
    let y = remainder(x, 360.0);
    if y == -180.0 {
        180.0
    } else {
        y
    }
}

// Replace angles outside [-90,90] with NaN
pub fn lat_fix(x: f64) -> f64 {
    if x.abs() > 90.0 {
        std::f64::NAN
    } else {
        x
    }
}

// compute y - x and reduce to [-180,180] accurately
pub fn ang_diff(x: f64, y: f64) -> (f64, f64) {
    let (d, t) = sum(ang_normalize(-x), ang_normalize(y));
    let d = ang_normalize(d);
    if d == 180.0 && t > 0.0 {
        sum(-180.0, t)
    } else {
        sum(d, t)
    }
}

pub fn fmod(x: f64, y: f64) -> f64 {
    x % y
}

/// Compute sine and cosine of x in degrees
pub fn sincosd(x: f64) -> (f64, f64) {
    // r = math.fmod(x, 360) if Math.isfinite(x) else Math.nan
    let mut r = if x.is_finite() {
        fmod(x, 360.0)
    } else {
        std::f64::NAN
    };

    // q = 0 if Math.isnan(r) else int(round(r / 90))
    let mut q = if r.is_nan() {
        0
    } else {
        (r / 90.0).round() as i32
    };

    // r -= 90 * q; r = math.radians(r)
    r -= 90.0 * q as f64;
    r = r.to_radians();

    // s = math.sin(r); c = math.cos(r)
    let s = r.sin();
    let c = r.cos();

    // q = q % 4
    q = q % 4;

    // if q == 1:
    //     s, c =  c, -s
    // elif q == 2:
    //     s, c = -s, -c
    // elif q == 3:
    //     s, c = -c,  s

    let q = if q < 0 { q + 4 } else { q };

    let (s, c) = if q == 1 {
        (c, -s)
    } else if q == 2 {
        (-s, -c)
    } else if q == 3 {
        (-c, s)
    } else {
        debug_assert_eq!(q, 0);
        (s, c)
    };

    // # Remove the minus sign on -0.0 except for sin(-0.0).
    // # On Windows 32-bit with python 2.7, math.fmod(-0.0, 360) = +0.0
    // # (x, c) here fixes this bug.  See also Math::sincosd in the C++ library.
    // # AngNormalize has a similar fix.
    //     s, c = (x, c) if x == 0 else (0.0+s, 0.0+c)
    // return s, c
    let (s, c) = if x == 0.0 { (x, c) } else { (0.0 + s, 0.0 + c) };

    (s, c)
}

// Compute atan2(y, x) with result in degrees
pub fn atan2d(y: f64, x: f64) -> f64 {
    let mut x = x;
    let mut y = y;
    let mut q = if y.abs() > x.abs() {
        let _x = x;
        x = y;
        y = _x;
        2.0
    } else {
        0.0
    };
    if x < 0.0 {
        q += 1.0;
        x = -x;
    }
    let mut ang = y.atan2(x).to_degrees();
    if q == 1.0 {
        ang = if y >= 0.0 { 180.0 - ang } else { -180.0 - ang };
    } else if q == 2.0 {
        ang = 90.0 - ang;
    } else if q == 3.0 {
        ang = -90.0 + ang;
    }
    ang
}

pub fn eatanhe(x: f64, es: f64) -> f64 {
    if es > 0.0 { es * (es * x).atanh() } else { -es * (es * x).atan() }
}

// Functions that used to be inside Geodesic
pub fn sin_cos_series(sinp: bool, sinx: f64, cosx: f64, c: &[f64]) -> f64 {
    let mut k = c.len();
    let mut n: i64 = k as i64 - if sinp { 1 } else { 0 };
    let ar: f64 = 2.0 * (cosx - sinx) * (cosx + sinx);
    let mut y1 = 0.0;
    let mut y0: f64 = if n & 1 != 0 {
        k -= 1;
        c[k]
    } else {
        0.0
    };
    n = n / 2;
    while n > 0 {
        n -= 1;
        k -= 1;
        y1 = ar * y0 - y1 + c[k];
        k -= 1;
        y0 = ar * y1 - y0 + c[k];
    }
    if sinp {
        2.0 * sinx * cosx * y0
    } else {
        cosx * (y0 - y1)
    }
}

// Solve astroid equation
pub fn astroid(x: f64, y: f64) -> f64 {
    let p = sq(x);
    let q = sq(y);
    let r = (p + q - 1.0) / 6.0;
    if !(q == 0.0 && r <= 0.0) {
        let s = p * q / 4.0;
        let r2 = sq(r);
        let r3 = r * r2;
        let disc = s * (s + 2.0 * r3);
        let mut u = r;
        if disc >= 0.0 {
            let mut t3 = s + r3;
            t3 += if t3 < 0.0 { -disc.sqrt() } else { disc.sqrt() };
            let t = cbrt(t3); // we could use built-in T.cbrt
            u += t + if t != 0.0 { r2 / t } else { 0.0 };
        } else {
            let ang = (-disc).sqrt().atan2(-(s + r3));
            u += 2.0 * r * (ang / 3.0).cos();
        }
        let v = (sq(u) + q).sqrt();
        let uv = if u < 0.0 { q / (v - u) } else { u + v };
        let w = (uv - q) / (2.0 * v);
        uv / ((uv + sq(w)).sqrt() + w)
    } else {
        0.0
    }
}

pub fn _A1m1f(eps: f64, geodesic_order: i64) -> f64 {
    const COEFF: [f64; 5] = [1.0, 4.0, 64.0, 0.0, 256.0];
    let m: i64 = geodesic_order / 2;
    let t = polyval(m as isize, &COEFF, sq(eps)) / COEFF[(m + 1) as usize] as f64;
    (t + eps) / (1.0 - eps)
}

pub fn _C1f(eps: f64, c: &mut [f64], geodesic_order: i64) {
    const COEFF: [f64; 18] = [
        -1.0, 6.0, -16.0, 32.0, -9.0, 64.0, -128.0, 2048.0, 9.0, -16.0, 768.0, 3.0, -5.0, 512.0,
        -7.0, 1280.0, -7.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    for l in 1..=geodesic_order {
        let m = ((geodesic_order - l) / 2) as i64;
        c[l as usize] =
            d * polyval(m as isize, &COEFF[o as usize..], eps2) / COEFF[(o + m + 1) as usize] as f64;
        o += m + 2;
        d *= eps;
    }
}

pub fn _C1pf(eps: f64, c: &mut [f64], geodesic_order: i64) {
    const COEFF: [f64; 18] = [
        205.0, -432.0, 768.0, 1536.0, 4005.0, -4736.0, 3840.0, 12288.0, -225.0, 116.0, 384.0,
        -7173.0, 2695.0, 7680.0, 3467.0, 7680.0, 38081.0, 61440.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    for l in 1..=geodesic_order {
        let m = (geodesic_order - l) / 2;
        c[l as usize] =
            d * polyval(m as isize, &COEFF[o as usize..], eps2) / COEFF[(o + m + 1) as usize] as f64;
        o += m + 2;
        d *= eps;
    }
}

pub fn _A2m1f(eps: f64, geodesic_order: i64) -> f64 {
    const COEFF: [f64; 5] = [-11.0, -28.0, -192.0, 0.0, 256.0];
    let m: i64 = geodesic_order / 2;
    let t = polyval(m as isize, &COEFF, sq(eps)) / COEFF[(m + 1) as usize] as f64;
    (t - eps) / (1.0 + eps)
}

pub fn _C2f(eps: f64, c: &mut [f64], geodesic_order: i64) {
    const COEFF: [f64; 18] = [
        1.0, 2.0, 16.0, 32.0, 35.0, 64.0, 384.0, 2048.0, 15.0, 80.0, 768.0, 7.0, 35.0, 512.0, 63.0,
        1280.0, 77.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    for l in 1..=geodesic_order {
        let m = (geodesic_order - l) / 2;
        c[l as usize] =
            d * polyval(m as isize, &COEFF[o as usize..], eps2) / COEFF[(o + m + 1) as usize] as f64;
        o += m + 2;
        d *= eps;
    }
}

// todo: remove
// pub fn fmt_hex(u: u64) -> String {
//     format!("{:x}", u)
// }
// pub fn fmt_i64(i: i64) -> String {
//     format!("{}", i)
// }
// pub fn fmt_f64(f: f64) -> String {
//     let sign_helper = if f.is_sign_negative() && (f.is_nan() || f == 0.0) {
//         "-"
//     } else {
//         ""
//     };
//     let root1 = format!("{}", f);
//     let root2 = format!("{:e}", f);
//     format!("{}{}", sign_helper, if root2.len() < root1.len() {root2} else {root1})
// }
// pub fn fmt_slice_f64(v: &[f64]) -> String {
//     let v2: Vec<String> = v.iter().map(|&o| fmt_f64(o)).collect();
//     v2.join(" ")
// }

#[cfg(test)]
mod tests {
    extern crate utilities;

    use super::*;
    use utilities::{assert_delta, util};
    use utilities::util::test_basic;
    use utilities::delta_entry::DeltaEntry;
    use std::sync::{Arc, Mutex};

    // Results for assertions in * tests are taken by running the python implementation

    #[test]
    fn test_sincosd() {
        let res = sincosd(-77.03196);
        assert_eq!(res.0, -0.9744953925159129);
        assert_eq!(res.1, 0.22440750870961693);

        let res = sincosd(69.48894);
        assert_eq!(res.0, 0.9366045700708676);
        assert_eq!(res.1, 0.3503881837653281);
        let res = sincosd(-1.0);
        assert_eq!(res.0, -0.01745240643728351);
        assert_eq!(res.1, 0.9998476951563913);
    }

    #[test]
    fn test__C2f() {
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        _C2f(0.12, &mut c, 6);
        assert_eq!(
            c,
            vec![
                1.0,
                0.0601087776,
                0.00270653103,
                0.000180486,
                1.4215824e-05,
                1.22472e-06,
                1.12266e-07
            ]
        )
    }

    #[test]
    fn test__A2m1f() {
        assert_eq!(_A2m1f(0.12, 6), -0.11680607884285714);
    }

    #[test]
    fn test__C1pf() {
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        _C1pf(0.12, &mut c, 6);
        assert_eq!(
            c,
            vec![
                1.0,
                0.059517321000000005,
                0.004421053215,
                0.0005074200000000001,
                6.997613759999999e-05,
                1.1233080000000001e-05,
                1.8507366e-06
            ]
        )
    }

    #[test]
    fn test__C1f() {
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        _C1f(0.12, &mut c, 6);
        assert_eq!(
            c,
            vec![
                1.0,
                -0.059676777599999994,
                -0.000893533122,
                -3.57084e-05,
                -2.007504e-06,
                -1.3607999999999999e-07,
                -1.0205999999999999e-08
            ]
        )
    }

    #[test]
    fn test__A1m1f() {
        assert_eq!(_A1m1f(0.12, 6), 0.1404582405272727);
    }

    #[test]
    fn test_astroid() {
        assert_eq!(astroid(21.0, 12.0), 23.44475767500982);
    }

    #[test]
    fn test_sin_cos_series() {
        assert_eq!(
            sin_cos_series(
                false,
                -0.8928657853278468,
                0.45032287238256896,
                &vec![
                    0.6660771734724675,
                    1.5757752625233906e-05,
                    3.8461688963148916e-09,
                    1.3040960748120204e-12,
                    5.252912023008548e-16,
                    2.367770858285795e-19
                ],
            ),
            0.29993425660538664
        );

        assert_eq!(
            sin_cos_series(
                false,
                -0.8928657853278468,
                0.45032287238256896,
                &vec![0., 1., 2., 3., 4., 5.],
            ),
            1.8998562852254026
        );
        assert_eq!(
            sin_cos_series(
                true,
                0.2969032234925426,
                0.9549075745221299,
                &vec![
                    0.0,
                    -0.0003561309485314716,
                    -3.170731714689771e-08,
                    -7.527972480734327e-12,
                    -2.5133854116682488e-15,
                    -1.0025061462383107e-18,
                    -4.462794158625518e-22
                ],
            ),
            -0.00020196665516199853
        );
        assert_eq!(
            sin_cos_series(
                true,
                -0.8928657853278468,
                0.45032287238256896,
                &vec![
                    0.0,
                    -0.0003561309485314716,
                    -3.170731714689771e-08,
                    -7.527972480734327e-12,
                    -2.5133854116682488e-15,
                    -1.0025061462383107e-18,
                    -4.462794158625518e-22
                ],
            ),
            0.00028635444718997857
        );

        assert_eq!(
            sin_cos_series(true, 0.12, 0.21, &vec![1.0, 2.0]),
            0.10079999999999999
        );
        assert_eq!(
            sin_cos_series(
                true,
                -0.024679833885152578,
                0.9996954065111039,
                &vec![
                    0.0,
                    -0.0008355098973052918,
                    -1.7444619952659748e-07,
                    -7.286557795511902e-11,
                    -3.80472772706481e-14,
                    -2.2251271876594078e-17,
                    1.2789961247944744e-20
                ],
            ),
            4.124513511893872e-05
        );
    }

    // *_vs_cpp_* tests are based on instrumented inputs and outputs from C++
    // They're flagged as "ignore" both because they're slow and not self-contained
    // (since they need to read data files), so they only run if specifically requested.

    #[test]
    #[ignore] // Fails current behavior. Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_ang_diff() {
        // Line format: x y result=z e-out
        test_basic("Math_AngDiff_x_y_e", 4, |line_num, items| {
            let result = ang_diff(items[0], items[1]);
            assert_delta!(items[2], result.0, 0.0, false, "result (result.0)", line_num);
            assert_delta!(items[3], result.1, 0.0, false, "e (result.1)", line_num);
        });
    }

    #[test]
    #[ignore] // Fails current behavior. Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_ang_normalize() {
        // Line format: x result
        test_basic("Math_AngNormalize", 2, |line_num, items| {
            let result = ang_normalize(items[0]);
            assert_delta!(items[1], result, 0.0, false, "result", line_num);
        });
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_ang_round() {
        // Line format: x result
        test_basic("Math_AngRound", 2, |line_num, items| {
            let result = ang_round(items[0]);
            assert_delta!(items[1], result, 0.0, false, "result", line_num);
        });
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_atan2d() {
        // Line format: y x result
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_atan2d ", &[
                // todo: this is a prime candidate for an ULP/step tolerance
                ("result", 3e-14, false, false),
            ])));
        test_basic("Math_atan2d", 3, |line_num, items| {
            let result = atan2d(items[0], items[1]);
            let mut entries = delta_entries.lock().unwrap();
            // For input of (-0, 1), the ideal result is +0, but C++ currently gets away with -0.
            // Ideally, we'd accept a value that is either positive 0,
            // or that matches the C++ result for this case.
            // For simplicity, require matching C++ for now.
            // let special_case = items[0] == 0.0 && items[0].is_sign_negative() && items[1] == 1.0 && result.is_sign_positive();
            entries[0].add(items[2], result, line_num);
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    // Placeholder: Math_atand

    #[test]
    #[ignore] // Relies on non-Karney outside files.
    fn test_vs_cpp_geomath_consts() {
        // Line format: digits digits10 extra_digits bigendian pi degree GEOGRAPHICLIB_PRECISION GEOGRAPHICLIB_WORDS_BIGENDIAN
        let items = util::read_consts_basic("Math_consts", 8);
        let line_num = 2;
        assert_delta!(items[0], DIGITS as f64, 0.0, false, "DIGITS", line_num);
        assert_delta!(items[4], std::f64::consts::PI, 0.0, false, "PI", line_num);
    }

    // Placeholder: Math_cosd

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_eatanhe() {
        // Line format: x es result
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_eatanhe ", &[
                ("result", 2e-18, false, false),
            ])));
        test_basic("Math_eatanhe", 3, |line_num, items| {
            let result = eatanhe(items[0], items[1]);
            let mut entries = delta_entries.lock().unwrap();
            entries[0].add(items[2], result, line_num);
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }
    
    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_lat_fix() {
        // Line format: x result
        test_basic("Math_LatFix", 2, |line_num, items| {
            let result = lat_fix(items[0]);
            assert_delta!(items[1], result, 0.0, false, "result", line_num);
        });
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_norm() {
        // Line format: x-in y-in x-out y-out
        test_basic("Math_norm", 4, |line_num, items| {
            let mut x = items[0];
            let mut y = items[1];
            norm(&mut x, &mut y);
            assert_delta!(items[2], x, 0.0, false, "x-out (result.0)", line_num);
            assert_delta!(items[3], y, 0.0, false, "y-out (result.1)", line_num);
        });
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_polyval() {
        // Line format: N p(N+1) x result
        test_basic("Math_polyval", -1, |line_num, items| {
            assert!(items.len() > 2, "Expected a minimum of 3 items per line. Line {} had {}.", line_num, items.len());
            let n = items[0] as isize;
            let p_len = if n < 0 { 0 } else { n as usize + 1 };
            let arg_count = p_len + 3;
            assert!(items.len() == arg_count, "Expected {} items on line {}, based on first item, but found {}.", arg_count, line_num, items.len());
            let p = &items[1..arg_count-2];
            let x = items[arg_count-2];
            assert_eq!(p_len, p.len(), "Internal error: On line {}, tried to construct slice size {} but found {}", line_num, p_len, p.len());
            let result = polyval(n, p, x);
            assert_delta!(items[arg_count - 1], result, 0.0, false, "result", line_num);
        });
    }

    #[test]
    #[ignore] // Fails current behavior. Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_sincosd() {
        // Line format: x sinx-out cosx-out
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_sincosd ", &[
                ("sinx-out (result.0)", 2e-16, false, false),
                ("cosx-out (result.1)", 2e-16, false, false),
            ])));
        test_basic("Math_sincosd", 3, |line_num, items| {
            let result = sincosd(items[0]);
            let mut entries = delta_entries.lock().unwrap();
            entries[0].add(items[1], result.0, line_num);
            entries[1].add(items[2], result.1, line_num);
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    // Placeholder: Math_sind

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_sq() {
        // Line format: x result
        test_basic("Math_sq", 2, |line_num, items| {
            let result = sq(items[0]);
            assert_delta!(items[1], result, 0.0, false, "result", line_num);
        });
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_sum() {
        // Line format: u v result=s t-out
        test_basic("Math_sum", 4, |line_num, items| {
            let result = sum(items[0], items[1]);
            assert_delta!(items[2], result.0, 0.0, false, "result (result.0)", line_num);
            assert_delta!(items[3], result.1, 0.0, false, "t-out (result.1)", line_num);
        });
    }

    // Placeholder: Math_swab
    // Placeholder: Math_tand
    // Placeholder: Math_tauf
    // Placeholder: Math_taupf

    // All remaining tests are for operations that are defined on Geodesic rather than Math in the C++ version.

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_astroid() {
        // Line format: x y result
        // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_astroid ", &[
                ("result", 5e-16, false, false),
            ])));
        test_basic("Geodesic_Astroid", 3, |line_num, items| {
            let result = astroid(items[0], items[1]);
            let mut entries = delta_entries.lock().unwrap();
            entries[0].add(items[2], result, line_num);
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_a1m1f() {
        // Line format: eps result
        // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_a1m1f ", &[
                ("result", 0.0, false, false),
            ])));
        test_basic("Geodesic_A1m1f", 2, |line_num, items| {
            let result = _A1m1f(items[0], crate::geodesic::GEODESIC_ORDER);
            let mut entries = delta_entries.lock().unwrap();
            entries[0].add(items[1], result, line_num);
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_a2m1f() {
        // Line format: eps result
        // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_a2m1f ", &[
                ("result", 0.0, false, false),
            ])));
        test_basic("Geodesic_A2m1f", 2, |line_num, items| {
            let result = _A2m1f(items[0], crate::geodesic::GEODESIC_ORDER);
            let mut entries = delta_entries.lock().unwrap();
            entries[0].add(items[1], result, line_num);
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_c1f() {
        // Line format: eps c-out(nC1_+1)
        // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_c1f ", &[
                ("c item", 0.0, false, false),
            ])));
        test_basic("Geodesic_C1f", crate::geodesic::GEODESIC_ORDER as isize + 2, |line_num, items| {
            let mut c = [0.0; crate::geodesic::GEODESIC_ORDER as usize + 1];
            _C1f(items[0], &mut c, crate::geodesic::GEODESIC_ORDER);
            let mut entries = delta_entries.lock().unwrap();
            for i in 0..c.len() {
                entries[0].add(items[1 + i], c[i], line_num);
            }
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_c1pf() {
        // Line format: eps c-out(nC1p_+1)
        // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_c1pf ", &[
                ("c item", 0.0, false, false),
            ])));
        test_basic("Geodesic_C1pf", crate::geodesic::GEODESIC_ORDER as isize + 2, |line_num, items| {
            let mut c = [0.0; crate::geodesic::GEODESIC_ORDER as usize + 1];
            _C1pf(items[0], &mut c, crate::geodesic::GEODESIC_ORDER);
            let mut entries = delta_entries.lock().unwrap();
            // Element 0 of this array is unused, and not initialized in c++, so value is unpredictable
            for i in 1..c.len() {
                entries[0].add(items[i + 1], c[i], line_num);
            }
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_c2f() {
        // Line format: eps c-out(nC1p_+1)
        // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geomath_c2f ", &[
                ("c item", 0.0, false, false),
            ])));
        test_basic("Geodesic_C2f", crate::geodesic::GEODESIC_ORDER as isize + 2, |line_num, items| {
            let mut c = [0.0; crate::geodesic::GEODESIC_ORDER as usize + 1];
            _C2f(items[0], &mut c, crate::geodesic::GEODESIC_ORDER);
            let mut entries = delta_entries.lock().unwrap();
            // The first element isn't used, so isn't useful to compare.
            for i in 1..c.len() {
                entries[0].add(items[1+i], c[i], line_num);
            }
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geomath_sin_cos_series() {
        // Line format: sinp sinx cosx n c(n+sinp) result
        // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
        test_basic("Geodesic_SinCosSeries", -1, |line_num, items| {
            assert!(items.len() > 4, "Expected a minimum of 5 items per line. Line {} had {}.", line_num, items.len());
            let sinp = items[0] != 0f64;
            let n = items[3] as i64;
            let mut s = if n < 0 { 0 } else { n as usize };
            if sinp {
                s += 1;
            }
            let arg_count = s + 5;
            assert!(items.len() == arg_count, "Expected {} items on line {}, based on first item, but found {}.", arg_count, line_num, items.len());
            let c = &items[4..s+4];
            assert_eq!(s, c.len(), "Internal error: On line {}, tried to construct slice size {} but found {}", line_num, s, c.len());
            let result = sin_cos_series(sinp, items[1], items[2], c);
            assert_delta!(items[arg_count - 1], result, 0.0, false, "result", line_num);
        });
    }

}
