#![allow(non_snake_case)]

use crate::geodesic::{self, GEODESIC_ORDER};
use crate::geodesiccapability as caps;
use crate::geomath;
use std::collections::HashMap;

#[derive(Debug)]
pub struct GeodesicLine {
    tiny_: f64, // This should be moved to consts
    _A1m1: f64,
    _A2m1: f64,
    _A3c: f64,
    _A4: f64,
    _B11: f64,
    _B21: f64,
    _B31: f64,
    _B41: f64,
    _C1a: [f64; GEODESIC_ORDER as usize + 1],
    _C1pa: [f64; GEODESIC_ORDER as usize + 1],
    _C2a: [f64; GEODESIC_ORDER as usize + 1],
    _C3a: [f64; GEODESIC_ORDER as usize + 1],
    _C4a: [f64; GEODESIC_ORDER as usize + 1],
    _b: f64,
    _c2: f64,
    _calp0: f64,
    _csig1: f64,
    _comg1: f64,
    _ctau1: f64,
    _dn1: f64,
    _f1: f64,
    _k2: f64,
    _salp0: f64,
    _somg1: f64,
    _ssig1: f64,
    _stau1: f64,
    a13: f64,
    a: f64,
    azi1: f64,
    calp1: f64,
    caps: u64,
    f: f64,
    lat1: f64,
    lon1: f64,
    s13: f64,
    salp1: f64,
}

impl GeodesicLine {
    pub fn new(
        geod: &geodesic::Geodesic,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        caps: Option<u64>,
        salp1: Option<f64>,
        calp1: Option<f64>,
    ) -> Self {
        let caps = match caps {
            None => caps::STANDARD | caps::DISTANCE_IN,
            Some(caps) => caps,
        };
        let salp1 = match salp1 {
            None => std::f64::NAN,
            Some(salp1) => salp1,
        };
        let calp1 = match calp1 {
            None => std::f64::NAN,
            Some(calp1) => calp1,
        };

        // This was taken from geodesic, putting it here for convenience
        let tiny_ = geomath::get_min_val().sqrt();

        let a = geod.a;
        let f = geod.f;
        let _b = geod._b;
        let _c2 = geod._c2;
        let _f1 = geod._f1;
        let caps = caps | caps::LATITUDE | caps::AZIMUTH | caps::LONG_UNROLL;
        let lat1 = geomath::lat_fix(lat1);
        let lon1 = lon1;
        let (azi1, salp1, calp1) = if salp1.is_nan() || calp1.is_nan() {
            let (salp1, calp1) = geomath::sincosd(geomath::ang_round(azi1));
            let azi1 = geomath::ang_normalize(azi1);
            (azi1, salp1, calp1)
        } else {
            (azi1, salp1, calp1)
        };

        let (mut sbet1, cbet1) = geomath::sincosd(geomath::ang_round(lat1));
        sbet1 *= _f1;
        let (sbet1, mut cbet1) = geomath::norm(sbet1, cbet1);
        cbet1 = tiny_.max(cbet1);
        let _dn1 = (1.0 + geod._ep2 * geomath::sq(sbet1)).sqrt();
        let _salp0 = salp1 * cbet1;
        let _calp0 = calp1.hypot(salp1 * sbet1);
        let _ssig1 = sbet1;
        let _somg1 = _salp0 * sbet1;
        let _csig1 = if sbet1 != 0.0 || calp1 != 0.0 {
            cbet1 * calp1
        } else {
            1.0
        };
        let _comg1 = _csig1;
        let (_ssig1, _csig1) = geomath::norm(_ssig1, _csig1);
        let _k2 = geomath::sq(_calp0) * geod._ep2;
        let eps = _k2 / (2.0 * (1.0 + (1.0 + _k2).sqrt()) + _k2);

        let mut _A1m1 = 0.0;
        let mut _C1a: [f64; GEODESIC_ORDER as usize + 1] = [0.0; GEODESIC_ORDER as usize + 1];
        let mut _B11 = 0.0;
        let mut _stau1 = 0.0;
        let mut _ctau1 = 0.0;
        if caps & caps::CAP_C1 != 0 {
            _A1m1 = geomath::_A1m1f(eps, geod.GEODESIC_ORDER);
            geomath::_C1f(eps, &mut _C1a, geod.GEODESIC_ORDER);
            _B11 = geomath::sin_cos_series(true, _ssig1, _csig1, &_C1a);
            let s = _B11.sin();
            let c = _B11.cos();
            _stau1 = _ssig1 * c + _csig1 * s;
            _ctau1 = _csig1 * c - _ssig1 * s;
        }

        let mut _C1pa: [f64; GEODESIC_ORDER as usize + 1] = [0.0; GEODESIC_ORDER as usize + 1];
        if caps & caps::CAP_C1p != 0 {
            geomath::_C1pf(eps, &mut _C1pa, geod.GEODESIC_ORDER);
        }

        let mut _A2m1 = 0.0;
        let mut _C2a: [f64; GEODESIC_ORDER as usize + 1] = [0.0; GEODESIC_ORDER as usize + 1];
        let mut _B21 = 0.0;
        if caps & caps::CAP_C2 != 0 {
            _A2m1 = geomath::_A2m1f(eps, geod.GEODESIC_ORDER);
            geomath::_C2f(eps, &mut _C2a, geod.GEODESIC_ORDER);
            _B21 = geomath::sin_cos_series(true, _ssig1, _csig1, &_C2a);
        }

        let mut _C3a: [f64; GEODESIC_ORDER as usize + 1] = [0.0; GEODESIC_ORDER as usize + 1];
        let mut _A3c = 0.0;
        let mut _B31 = 0.0;
        if caps & caps::CAP_C3 != 0 {
            geod._C3f(eps, &mut _C3a);
            _A3c = -f * _salp0 * geod._A3f(eps);
            _B31 = geomath::sin_cos_series(true, _ssig1, _csig1, &_C3a);
        }

        let mut _C4a: [f64; GEODESIC_ORDER as usize + 1] = [0.0; GEODESIC_ORDER as usize + 1];
        let mut _A4 = 0.0;
        let mut _B41 = 0.0;
        if caps & caps::CAP_C4 != 0 {
            geod._C4f(eps, &mut _C4a);
            _A4 = geomath::sq(a) * _calp0 * _salp0 * geod._e2;
            _B41 = geomath::sin_cos_series(false, _ssig1, _csig1, &_C4a);
        }

        let s13 = std::f64::NAN;
        let a13 = std::f64::NAN;

        GeodesicLine {
            tiny_,
            _A1m1,
            _A2m1,
            _A3c,
            _A4,
            _B11,
            _B21,
            _B31,
            _B41,
            _C1a,
            _C1pa,
            _comg1,
            _C2a,
            _C3a,
            _C4a,
            _b,
            _c2,
            _calp0,
            _csig1,
            _ctau1,
            _dn1,
            _f1,
            _k2,
            _salp0,
            _somg1,
            _ssig1,
            _stau1,
            a,
            a13,
            azi1,
            calp1,
            caps,
            f,
            lat1,
            lon1,
            s13,
            salp1,
        }
    }

    /// returns (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)
    pub fn _gen_position(
        &self,
        arcmode: bool,
        s12_a12: f64,
        outmask: u64,
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let mut a12 = std::f64::NAN;
        let mut lat2 = std::f64::NAN;
        let mut lon2 = std::f64::NAN;
        let mut azi2 = std::f64::NAN;
        let mut s12 = std::f64::NAN;
        let mut m12 = std::f64::NAN;
        let mut M12 = std::f64::NAN;
        let mut M21 = std::f64::NAN;
        let mut S12 = std::f64::NAN;
        let outmask = outmask & (self.caps & caps::OUT_MASK);
        if !(arcmode || (self.caps & (caps::OUT_MASK & caps::DISTANCE_IN) != 0)) {
            return (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12);
        }

        let mut B12 = 0.0;
        let mut AB1 = 0.0;
        let mut sig12: f64;
        let mut ssig12: f64;
        let mut csig12: f64;
        let mut ssig2: f64;
        let mut csig2: f64;
        if arcmode {
            sig12 = s12_a12.to_radians();
            let res = geomath::sincosd(s12_a12);
            ssig12 = res.0;
            csig12 = res.0;
        } else {
            // tau12 = s12_a12 / (self._b * (1 + self._A1m1))
            // tau12 = tau12 if Math.isfinite(tau12) else Math.nan
            let tau12 = s12_a12 / (self._b * (1.0 + self._A1m1));
            let tau12 = if tau12.is_finite() {
                tau12
            } else {
                std::f64::NAN
            };

            let s = tau12.sin();
            let c = tau12.cos();

            B12 = -geomath::sin_cos_series(
                true,
                self._stau1 * c + self._ctau1 * s,
                self._ctau1 * c - self._stau1 * s,
                &self._C1pa,
            );
            sig12 = tau12 - (B12 - self._B11);
            ssig12 = sig12.sin();
            csig12 = sig12.cos();
            if self.f.abs() > 0.01 {
                ssig2 = self._ssig1 * csig12 + self._csig1 * ssig12;
                csig2 = self._csig1 * csig12 - self._ssig1 * ssig12;
                B12 = geomath::sin_cos_series(true, ssig2, csig2, &self._C1a);
                let serr = (1.0 + self._A1m1) * (sig12 + (B12 - self._B11)) - s12_a12 / self._b;
                sig12 = sig12 - serr / (1.0 + self._k2 * ssig2.sqrt()).sqrt();
                ssig12 = sig12.sin();
                csig12 = sig12.cos();
            }
        };
        ssig2 = self._ssig1 * csig12 + self._csig1 * ssig12;
        csig2 = self._csig1 * csig12 - self._ssig1 * ssig12;
        let dn2 = (1.0 + self._k2 * geomath::sq(ssig2)).sqrt();
        if outmask & (caps::DISTANCE | caps::REDUCEDLENGTH | caps::GEODESICSCALE) != 0 {
            if arcmode || self.f.abs() > 0.01 {
                B12 = geomath::sin_cos_series(true, ssig2, csig2, &self._C1a);
            }
            AB1 = (1.0 + self._A1m1) * (B12 - self._B11);
        }

        let sbet2 = self._calp0 * ssig2;
        let mut cbet2 = self._salp0.hypot(self._calp0 * csig2);
        if cbet2 == 0.0 {
            cbet2 = self.tiny_;
            csig2 = self.tiny_;
        }
        let salp2 = self._salp0;
        let calp2 = self._calp0 * csig2;
        if outmask & caps::DISTANCE != 0 {
            s12 = if arcmode {
                self._b * ((1.0 + self._A1m1) * sig12 + AB1)
            } else {
                s12_a12
            }
        }
        if outmask & caps::LONGITUDE != 0 {
            let somg2 = self._salp0 * ssig2;
            let comg2 = csig2;
            let E = (1.0 as f64).copysign(self._salp0);
            let omg12 = if outmask & caps::LONG_UNROLL != 0 {
                E * (sig12 - (ssig2.atan2(csig2) - self._ssig1.atan2(self._csig1))
                    + ((E * somg2).atan2(comg2) - (E * self._somg1).atan2(self._comg1)))
            } else {
                (somg2 * self._comg1 - comg2 * self._somg1)
                    .atan2(comg2 * self._comg1 + somg2 * self._somg1)
            };
            let lam12 = omg12
                + self._A3c
                    * (sig12
                        + (geomath::sin_cos_series(true, ssig2, csig2, &self._C3a) - self._B31));
            let lon12 = lam12.to_degrees();
            lon2 = if outmask & caps::LONG_UNROLL != 0 {
                self.lon1 + lon12
            } else {
                geomath::ang_normalize(
                    geomath::ang_normalize(self.lon1) + geomath::ang_normalize(lon12),
                )
            };
        };

        if outmask & caps::LATITUDE != 0 {
            lat2 = geomath::atan2d(sbet2, self._f1 * cbet2);
        }
        if outmask & caps::AZIMUTH != 0 {
            azi2 = geomath::atan2d(salp2, calp2);
        }
        if outmask & (caps::REDUCEDLENGTH | caps::GEODESICSCALE) != 0 {
            let B22 = geomath::sin_cos_series(true, ssig2, csig2, &self._C2a);
            let AB2 = (1.0 + self._A2m1) * (B22 - self._B21);
            let J12 = (self._A1m1 - self._A2m1) * sig12 + (AB1 - AB2);
            if outmask & caps::REDUCEDLENGTH != 0 {
                m12 = self._b
                    * ((dn2 * (self._csig1 * ssig2) - self._dn1 * (self._ssig1 * csig2))
                        - self._csig1 * csig2 * J12);
            }
            if outmask & caps::GEODESICSCALE != 0 {
                let t =
                    self._k2 * (ssig2 - self._ssig1) * (ssig2 + self._ssig1) / (self._dn1 + dn2);
                M12 = csig12 + (t * ssig2 - csig2 * J12) * self._ssig1 / self._dn1;
                M21 = csig12 - (t * self._ssig1 - self._csig1 * J12) * ssig2 / dn2;
            }
        }
        if outmask & caps::AREA != 0 {
            let B42 = geomath::sin_cos_series(false, ssig2, csig2, &self._C4a);
            let salp12: f64;
            let calp12: f64;
            if self._calp0 == 0.0 || self._salp0 == 0.0 {
                salp12 = salp2 * self.calp1 - calp2 * self.salp1;
                calp12 = calp2 * self.calp1 + salp2 * self.salp1;
            } else {
                salp12 = self._calp0
                    * self._salp0
                    * (if csig12 <= 0.0 {
                        self._csig1 * (1.0 - csig12) + ssig12 * self._ssig1
                    } else {
                        ssig12 * (self._csig1 * ssig12 / (1.0 + csig12) + self._ssig1)
                    });
                calp12 = geomath::sq(self._salp0) + geomath::sq(self._calp0) * self._csig1 * csig2;
            }
            S12 = self._c2 * salp12.atan2(calp12) + self._A4 * (B42 - self._B41);
        }
        a12 = if arcmode { s12_a12 } else { sig12.to_degrees() };
        (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)
    }

    // not currently used, but maybe some day
    #[allow(dead_code)]
    pub fn Position(&self, s12: f64, outmask: Option<u64>) -> HashMap<String, f64> {
        let outmask = match outmask {
            Some(outmask) => outmask,
            None => caps::STANDARD,
        };
        let mut result: HashMap<String, f64> = HashMap::new();
        result.insert("lat1".to_string(), self.lat1);
        result.insert("azi1".to_string(), self.azi1);
        result.insert("s12".to_string(), s12);
        let lon1 = if outmask & caps::LONG_UNROLL != 0 {
            self.lon1
        } else {
            geomath::ang_normalize(self.lon1)
        };
        result.insert("lon1".to_string(), lon1);

        let (a12, lat2, lon2, azi2, _s12, m12, M12, M21, S12) =
            self._gen_position(false, s12, outmask);
        let outmask = outmask & caps::OUT_MASK;
        result.insert("a12".to_string(), a12);
        if outmask & caps::LATITUDE != 0 {
            result.insert("lat2".to_string(), lat2);
        }
        if outmask & caps::LONGITUDE != 0 {
            result.insert("lon2".to_string(), lon2);
        }
        if outmask & caps::AZIMUTH != 0 {
            result.insert("azi2".to_string(), azi2);
        }
        if outmask & caps::REDUCEDLENGTH != 0 {
            result.insert("m12".to_string(), m12);
        }
        if outmask & caps::GEODESICSCALE != 0 {
            result.insert("M12".to_string(), M12);
            result.insert("M21".to_string(), M21);
        }
        if outmask & caps::AREA != 0 {
            result.insert("S12".to_string(), S12);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    extern crate utilities;

    use super::*;
    use geodesic::Geodesic;
    use std::sync::{Arc, Mutex};
    use utilities::{log_assert_delta, util};
    use utilities::util::test_basic;
    use utilities::delta_entry::DeltaEntry;

    #[test]
    fn test_gen_position() {
        let geod = Geodesic::wgs84();
        let gl = GeodesicLine::new(&geod, 0.0, 0.0, 10.0, None, None, None);
        let res = gl._gen_position(false, 150.0, 3979);
        assert_eq!(res.0, 0.0013520059461334633);
        assert_eq!(res.1, 0.0013359451088740494);
        assert_eq!(res.2, 0.00023398621812867812);
        assert_eq!(res.3, 10.000000002727887);
        assert_eq!(res.4, 150.0);
        assert_eq!(res.5.is_nan(), true);
        assert_eq!(res.6.is_nan(), true);
        assert_eq!(res.7.is_nan(), true);
        assert_eq!(res.8.is_nan(), true);
    }

    #[test]
    fn test_init() {
        let geod = Geodesic::wgs84();
        let gl = GeodesicLine::new(&geod, 0.0, 0.0, 0.0, None, None, None);
        assert_eq!(gl.a, 6378137.0);
        assert_eq!(gl.f, 0.0033528106647474805);
        assert_eq!(gl._b, 6356752.314245179);
        assert_eq!(gl._c2, 40589732499314.76);
        assert_eq!(gl._f1, 0.9966471893352525);
        assert_eq!(gl.caps, 36747);
        assert_eq!(gl.lat1, 0.0);
        assert_eq!(gl.lon1, 0.0);
        assert_eq!(gl.azi1, 0.0);
        assert_eq!(gl.salp1, 0.0);
        assert_eq!(gl.calp1, 1.0);
        assert_eq!(gl._dn1, 1.0);
        assert_eq!(gl._salp0, 0.0);
        assert_eq!(gl._calp0, 1.0);
        assert_eq!(gl._ssig1, 0.0);
        assert_eq!(gl._somg1, 0.0);
        assert_eq!(gl._csig1, 1.0);
        assert_eq!(gl._comg1, 1.0);
        assert_eq!(gl._k2, geod._ep2);
        assert_eq!(gl.s13.is_nan(), true);
        assert_eq!(gl.a13.is_nan(), true);
    }

    // *_vs_cpp_* tests are based on instrumented inputs and outputs from C++.
    // These tests are flagged as "ignore" because they're slow and not self-contained
    // (since they need to read data files), so they only run if specifically requested.

    // Set GeodesicLine properties based on an array of items.
    // Needed because instrumented log files sometimes don't know which constructor was used to reach
    // the current state, and using the wrong constructor can result in small precision differences.
    fn set_line(line: &mut GeodesicLine, items: &[f64]) {
        assert_eq!(line.a, items[0]);
        assert_eq!(line.f, items[1]);
        line.lat1 = items[2];
        line.lon1 = items[3];
        line.azi1 = items[4];
        line.a13 = items[5];
        line.s13 = items[6];
        line.caps = items[7] as u64;
        line._salp0 = items[8];
        line._calp0 = items[9];
        line.tiny_ = items[10];
        line._b = items[11];
        line._c2 = items[12];
        line._f1 = items[13];
        line._k2 = items[14];
        line.salp1 = items[15];
        line.calp1 = items[16];
        line._ssig1 = items[17];
        line._csig1 = items[18];
        line._dn1 = items[19];
        line._stau1 = items[20];
        line._ctau1 = items[21];
        line._somg1 = items[22];
        line._comg1 = items[23];
        line._A1m1 = items[24];
        line._A2m1 = items[25];
        line._A3c = items[26];
        line._B11 = items[27];
        line._B21 = items[28];
        line._B31 = items[29];
        line._A4 = items[30];
        line._B41 = items[31];
        // _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)

        let mut i = 31;
        assert_eq!(GEODESIC_ORDER as usize + 1, line._C1a.len(), "self._C1a size mismatch");
        for item in &mut line._C1a {
            i += 1;
            *item = items[i];
        }

        assert_eq!(GEODESIC_ORDER as usize + 1, line._C1pa.len(), "self._C1pa size mismatch");
        for item in &mut line._C1pa {
            i += 1;
            *item = items[i];
        }

        assert_eq!(GEODESIC_ORDER as usize + 1, line._C2a.len(), "self._C2a size mismatch");
        for item in &mut line._C2a {
            i += 1;
            *item = items[i];
        }

        assert_eq!(GEODESIC_ORDER as usize, line._C3a.len(), "self._C3a size mismatch");
        for item in &mut line._C3a {
            i += 1;
            *item = items[i];
        }

        assert_eq!(GEODESIC_ORDER as usize, line._C4a.len(), "self._C4a size mismatch");
        for item in &mut line._C4a {
            i += 1;
            *item = items[i];
        }
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files
    fn test_vs_cpp_geodesicline_consts() {
        // Format: nC1_ nC1p_ nC2_ nC3_ nC4_
        let items = util::read_consts_basic("GeodesicLine_consts", 5);
        log_assert_delta("test_vs_cpp_geodesicline_consts", "nC1_", items[0], GEODESIC_ORDER as f64, 0.0, false);
        log_assert_delta("test_vs_cpp_geodesicline_consts", "nC1p_", items[1], GEODESIC_ORDER as f64, 0.0, false);
        log_assert_delta("test_vs_cpp_geodesicline_consts", "nC2_", items[2], GEODESIC_ORDER as f64, 0.0, false);
        log_assert_delta("test_vs_cpp_geodesicline_consts", "nC3_", items[3], GEODESIC_ORDER as f64, 0.0, false);
        log_assert_delta("test_vs_cpp_geodesicline_consts", "nC4_", items[4], GEODESIC_ORDER as f64, 0.0, false);
    }

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geodesicline_gen_position() {
        // Format: this-in[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)] arcmode s12_a12 outmask result=a12 lat2-out lon2-out azi2-out s12-out m12-out M12-out M21-out S12-out
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geodesicline_gen_position ", &[
                ("result (a12 or result.0)", 0.0, true, false),
                ("lat2-out (result.1)", 0.0, false, false),
                ("lon2-out (result.2)", 0.0, false , false),
                ("azi2-out (result.3)", 0.0, false, false),
                ("s12-out (result.4)", 0.0, true, false),
                ("m12-out (result.5)", 0.0, true, false),
                ("M12-out (result.6)", 0.0, true, false),
                ("M21-out (result.7)", 0.0, true, false),
                ("S12-out (result.8)", 0.0, true, false),
            ])));
        test_basic("GeodesicLine_GenPosition", 44 + 3 + 5 * GEODESIC_ORDER as isize, |line_num, items| {
            let g = Geodesic::new(items[0], items[1]);
            let mut line = GeodesicLine::new(&g, items[2], items[3], items[4], None, None, None);
            set_line(&mut line, items);
            // non-construction: arcmode s12_a12 outmask result=a12 lat2-out lon2-out azi2-out s12-out m12-out M12-out M21-out S12-out
            // 77 total items, 12 non-construction, so first non-construction is at 65
            let outmask = items[67] as u64;
            let result = line._gen_position(0.0 != items[65], items[66], outmask);
            let mut entries = delta_entries.lock().unwrap();
            entries[0].add(items[68], result.0, line_num);
            if outmask & caps::LATITUDE != 0 {
                entries[1].add(items[69], result.1, line_num);
            }
            if outmask & caps::LONGITUDE != 0 {
                entries[2].add(items[70], result.2, line_num);
            }
            if outmask & caps::AZIMUTH != 0 {
                entries[3].add(items[71], result.3, line_num);
            }
            if outmask & caps::DISTANCE != 0 {
                entries[4].add(items[72], result.4, line_num);
            }
            if outmask & caps::REDUCEDLENGTH != 0 {
                entries[5].add(items[73], result.5, line_num);
            }
            if outmask & caps::GEODESICSCALE != 0 {
                entries[6].add(items[74], result.6, line_num);
                entries[7].add(items[75], result.7, line_num);
            }
            if outmask & caps::AREA != 0 {
                entries[8].add(items[76], result.8, line_num);
            }
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    // placeholder: GeodesicLine_GenSetDistance

    #[test]
    #[ignore] // Relies on non-Karney outside files. Slow.
    fn test_vs_cpp_geodesicline_new_a() {
        // Format: g[_a _f] lat1 lon1 azi1 caps this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)]
        let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
            "test_vs_cpp_geodesicline_new_a ", &[
                ("a", 0.0, false, false),
                ("f", 0.0, false, false),
                ("lat1", 0.0, false, false),
                ("lon1", 0.0, false, false),
                ("azi1", 0.0, false, false),
                ("a13", 0.0, false, false),
                ("s13", 0.0, false, false),
                ("caps", 0.0, false, false),
                ("salp0", 0.0, false, false),
                ("calp0", 0.0, false, false),
                ("tiny", 0.0, false, false),
                ("b", 0.0, false, false),
                ("c2", 0.0, false, false),
                ("f1", 0.0, false, false),
                ("k2", 0.0, false, false),
                ("salp1", 0.0, false, false),
                ("calp1", 0.0, false, false),
                ("ssig1", 0.0, false, false),
                ("csig1", 0.0, false, false),
                ("dn1", 0.0, false, false),
                ("stau1", 0.0, false, false),
                ("ctau1", 0.0, false, false),
                ("somg1", 0.0, false, false),
                ("comg1", 0.0, false, false),
                ("A1m1", 0.0, false, false),
                ("A2m1", 0.0, false, false),
                ("A3c", 0.0, false, false),
                ("B11", 0.0, false, false),
                ("B21", 0.0, false, false),
                ("B31", 0.0, false, false),
                ("A4", 0.0, false, false),
                ("B41", 0.0, false, false),
                ("C1a item", 0.0, false, false),
                ("C1pa item", 0.0, false, false),
                ("C2a item", 0.0, false, false),
                ("C3a item", 0.0, false, false),
                ("C4a item", 0.0, false, false),
            ])));
        // 38 +... _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)
        test_basic("GeodesicLine_GeodesicLine_5arg", 38 + 3 + 5 * (GEODESIC_ORDER as isize), |line_num, items| {
            let g = Geodesic::new(items[0], items[1]);
            let caps = items[5] as u64;
            let line = GeodesicLine::new(&g, items[2], items[3], items[4], Some(caps), None, None);
            let mut entries = delta_entries.lock().unwrap();
            entries[0].add(items[6], line.a, line_num);
            entries[1].add(items[7], line.f, line_num);
            entries[2].add(items[8], line.lat1, line_num);
            entries[3].add(items[9], line.lon1, line_num);
            entries[4].add(items[10], line.azi1, line_num);
            entries[5].add(items[11], line.a13, line_num);
            entries[6].add(items[12], line.s13, line_num);
            entries[7].add(items[13], line.caps as f64, line_num);
            entries[8].add(items[14], line._salp0, line_num);
            entries[9].add(items[15], line._calp0, line_num);
            entries[10].add(items[16], line.tiny_, line_num);
            entries[11].add(items[17], line._b, line_num);
            entries[12].add(items[18], line._c2, line_num);
            entries[13].add(items[19], line._f1, line_num);
            entries[14].add(items[20], line._k2, line_num);
            entries[15].add(items[21], line.salp1, line_num);
            entries[16].add(items[22], line.calp1, line_num);
            entries[17].add(items[23], line._ssig1, line_num);
            entries[18].add(items[24], line._csig1, line_num);
            entries[19].add(items[25], line._dn1, line_num);
            entries[20].add(items[26], line._stau1, line_num);
            entries[21].add(items[27], line._ctau1, line_num);
            entries[22].add(items[28], line._somg1, line_num);
            entries[23].add(items[29], line._comg1, line_num);
            entries[24].add(items[30], line._A1m1, line_num);
            entries[25].add(items[31], line._A2m1, line_num);
            entries[26].add(items[32], line._A3c, line_num);
            entries[27].add(items[33], line._B11, line_num);
            entries[28].add(items[34], line._B21, line_num);
            entries[29].add(items[35], line._B31, line_num);
            entries[30].add(items[36], line._A4, line_num);
            entries[31].add(items[37], line._B41, line_num);

            let mut i = 37;
            assert_eq!(GEODESIC_ORDER as usize + 1, line._C1a.len(), "self._C1a size mismatch");
            for item in &line._C1a {
                i += 1;
                entries[32].add(items[i], *item, line_num);
            }

            assert_eq!(GEODESIC_ORDER as usize + 1, line._C1pa.len(), "self._C1pa size mismatch");
            for item in &line._C1pa {
                i += 1;
                entries[33].add(items[i], *item, line_num);
            }

            assert_eq!(GEODESIC_ORDER as usize + 1, line._C2a.len(), "self._C2a size mismatch");
            for item in &line._C2a {
                i += 1;
                entries[34].add(items[i], *item, line_num);
            }

            assert_eq!(GEODESIC_ORDER as usize, line._C3a.len(), "self._C3a size mismatch");
            for item in &line._C3a {
                i += 1;
                entries[35].add(items[i], *item, line_num);
            }

            assert_eq!(GEODESIC_ORDER as usize, line._C4a.len(), "self._C4a size mismatch");
            for item in &line._C4a {
                i += 1;
                entries[36].add(items[i], *item, line_num);
            }
        });
        println!();
        delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
        delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    }

    // #[test]
    // #[ignore] // Relies on non-Karney outside files. Slow.
    // fn test_vs_cpp_geodesicline_new_b() {
    //     // Format: g[_a _f] lat1 lon1 azi1 salp1 calp1 caps arcmode s13_a13 this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)]
    //     let delta_entries = Arc::new(Mutex::new(DeltaEntry::new_vec(
    //         "test_vs_cpp_geodesicline_new_b ", &[
    //             ("a", 0.0, false, false),
    //             ("f", 0.0, false, false),
    //             ("lat1", 0.0, false, false),
    //             ("lon1", 0.0, false, false),
    //             ("azi1", 0.0, false, false),
    //             ("a13", 0.0, false, false),
    //             ("s13", 0.0, false, false),
    //             ("caps", 0.0, false, false),
    //             ("salp0", 0.0, false, false),
    //             ("calp0", 0.0, false, false),
    //             ("tiny", 0.0, false, false),
    //             ("b", 0.0, false, false),
    //             ("c2", 0.0, false, false),
    //             ("f1", 0.0, false, false),
    //             ("k2", 0.0, false, false),
    //             ("salp1", 0.0, false, false),
    //             ("calp1", 0.0, false, false),
    //             ("ssig1", 0.0, false, false),
    //             ("csig1", 0.0, false, false),
    //             ("dn1", 0.0, false, false),
    //             ("stau1", 0.0, false, false),
    //             ("ctau1", 0.0, false, false),
    //             ("somg1", 0.0, false, false),
    //             ("comg1", 0.0, false, false),
    //             ("A1m1", 0.0, false, false),
    //             ("A2m1", 0.0, false, false),
    //             ("A3c", 0.0, false, false),
    //             ("B11", 0.0, false, false),
    //             ("B21", 0.0, false, false),
    //             ("B31", 0.0, false, false),
    //             ("A4", 0.0, false, false),
    //             ("B41", 0.0, false, false),
    //             ("C1a item", 0.0, false, false),
    //             ("C1pa item", 0.0, false, false),
    //             ("C2a item", 0.0, false, false),
    //             ("C3a item", 0.0, false, false),
    //             ("C4a item", 0.0, false, false),
    //         ])));
    //     // 42 +... _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)
    //     test_basic("GeodesicLine_GeodesicLine_9arg", 42 + 3 + 5 * (GEODESIC_ORDER as isize), |line_num, items| {
    //         let g = Geodesic::new(items[0], items[1]);
    //         let caps = items[7] as u64;
    //         // todo: modify to pass arcmode and s12_a13, but we'll need changes to "new" or something
    //         let line = GeodesicLine::new(&g, items[2], items[3], items[4], Some(caps), Some(items[5]), Some(items[6]));
    //         let mut entries = delta_entries.lock().unwrap();
    //         entries[0].add(items[10], line.a, line_num);
    //         entries[1].add(items[11], line.f, line_num);
    //         entries[2].add(items[12], line.lat1, line_num);
    //         entries[3].add(items[13], line.lon1, line_num);
    //         entries[4].add(items[14], line.azi1, line_num);
    //         entries[5].add(items[15], line.a13, line_num);
    //         entries[6].add(items[16], line.s13, line_num);
    //         entries[7].add(items[17], line.caps as f64, line_num);
    //         entries[8].add(items[18], line._salp0, line_num);
    //         entries[9].add(items[19], line._calp0, line_num);
    //         entries[10].add(items[20], line.tiny_, line_num);
    //         entries[11].add(items[21], line._b, line_num);
    //         entries[12].add(items[22], line._c2, line_num);
    //         entries[13].add(items[23], line._f1, line_num);
    //         entries[14].add(items[24], line._k2, line_num);
    //         entries[15].add(items[25], line.salp1, line_num);
    //         entries[16].add(items[26], line.calp1, line_num);
    //         entries[17].add(items[27], line._ssig1, line_num);
    //         entries[18].add(items[28], line._csig1, line_num);
    //         entries[19].add(items[29], line._dn1, line_num);
    //         entries[20].add(items[30], line._stau1, line_num);
    //         entries[21].add(items[31], line._ctau1, line_num);
    //         entries[22].add(items[32], line._somg1, line_num);
    //         entries[23].add(items[33], line._comg1, line_num);
    //         entries[24].add(items[34], line._A1m1, line_num);
    //         entries[25].add(items[35], line._A2m1, line_num);
    //         entries[26].add(items[36], line._A3c, line_num);
    //         entries[27].add(items[37], line._B11, line_num);
    //         entries[28].add(items[38], line._B21, line_num);
    //         entries[29].add(items[39], line._B31, line_num);
    //         entries[30].add(items[40], line._A4, line_num);
    //         entries[31].add(items[41], line._B41, line_num);

    //         let mut i = 41;
    //         assert_eq!(GEODESIC_ORDER as usize + 1, line._C1a.len(), "self._C1a size mismatch");
    //         for item in &line._C1a {
    //             i += 1;
    //             entries[32].add(items[i], *item, line_num);
    //         }

    //         assert_eq!(GEODESIC_ORDER as usize + 1, line._C1pa.len(), "self._C1pa size mismatch");
    //         for item in &line._C1pa {
    //             i += 1;
    //             entries[33].add(items[i], *item, line_num);
    //         }

    //         assert_eq!(GEODESIC_ORDER as usize + 1, line._C2a.len(), "self._C2a size mismatch");
    //         for item in &line._C2a {
    //             i += 1;
    //             entries[34].add(items[i], *item, line_num);
    //         }

    //         assert_eq!(GEODESIC_ORDER as usize, line._C3a.len(), "self._C3a size mismatch");
    //         for item in &line._C3a {
    //             i += 1;
    //             entries[35].add(items[i], *item, line_num);
    //         }

    //         assert_eq!(GEODESIC_ORDER as usize, line._C4a.len(), "self._C4a size mismatch");
    //         for item in &line._C4a {
    //             i += 1;
    //             entries[36].add(items[i], *item, line_num);
    //         }
    //     });
    //     println!();
    //     delta_entries.lock().unwrap().iter().for_each(|entry| println!("{}", entry));
    //     delta_entries.lock().unwrap().iter().for_each(|entry| entry.assert());
    // }

    // placeholder: GeodesicLine_LineInit
    // placeholder: GeodesicLine_SetArc
    // placeholder: GeodesicLine_SetDistance

}
