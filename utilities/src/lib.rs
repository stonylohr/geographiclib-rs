// todo: add difference calculations from tests/GeodTest to rs utilities (and review vs geod_error from validate-geographiclib
// todo: modify rs cpplog reader code to check for geographiclib and baseline version numbers in header line, and fail if they don't match expected
// todo: add rs test cases based on tests/ClosestApproach
// todo: add rs test cases based on tests/ConicTest
// todo: add rs test cases based on tests/EllipticTest
// todo: add rs test cases based on tests/GeodExact
// todo: add rs test cases based on tests/GeodShort
// todo: add rs test cases based on tests/GeodTest
// todo: add rs test cases based on tests/HarmTest
// todo: add rs test cases based on tests/intersect
// todo: add rs test cases based on tests/LevelEllipsoid
// todo: add rs test cases based on tests/M12zero
// todo: add rs test cases based on tests/NaNTester
// todo: add rs test cases based on tests/NormalTest
// todo: add rs test cases based on tests/ProjTest
// todo: add rs test cases based on tests/reformat
// todo: add rs test cases based on tests/test-distribution
// todo: add rs test cases based on tests/TMTest
// todo: add rs test cases based on tests/example-Geodesic
// todo: add tests from geographiclib tests folder to instrumented-crude cpp variation, rerun cpp baseline, and upload zip files
// todo: update rs local copy of cpplog data files


// testing tip: If you want full test output in a predictable order, including ignored tests,
//              use "cargo test -- --include-ignored --nocapture --test-threads 1"
//              (Note: include-ignored is still only available in nightly as of 2021/01/15, but looks like it's coming to stable soon)
// testing tip: If you want full test output in a predictable order, including ONLY ignored tests,
//              use "cargo test -- --ignored --nocapture --test-threads 1"
// testing tip: If a new "*_vs_cpp_*" test fails especially badly,
//              review the C++ "instrumented-crude" logging for possible bugs there
// testing tip: While assert variations are valuable for catching regressions, they're
//              often of more limited value when trying to _improve_ result precision.
//              In such cases, use of something like DeltaEntry is often helpful,
//              but remember to use nocapture, and often directing output to a file, for
//              syntax like "cargo test -- --nocapture > ../test-out-a.txt".
// benchmarking reminder: Remember to shut down unneeded processes before benching.
// benchmarking reminder: In some cases, it may be desirable to capture multiple before and
//                        after bench runs, to get a clearer sense of natural variability.

pub mod util;
pub mod delta_entry;

#[allow(non_upper_case_globals)]
pub const nC_: usize = 7; // todo: define in geodesic.rs instead

// Check whether x and y are within delta of each other, considering both absolute and relative difference.
// Consider them close enough if absolute difference is within diff_allowed inclusive,
// or if relative difference is within diff_allowed exclusive.
// nan values are considered equal to each other, and outside delta from everything else.
#[macro_export]
macro_rules! assert_delta {
    ($x:expr, $y:expr, $diff_allowed:expr, $allow_sign_change:expr, $msg:expr, $line_num:expr) => {
            let (diff_abs, diff_rel) = util::calc_delta($x, $y);
            assert!(
                    diff_abs <= $diff_allowed || diff_rel < $diff_allowed,
                    "assert_delta failed line {}, {}: {}{:e} vs {}{:e} diff {:e} rel {:e} outside inclusive {:e}", $line_num, $msg, util::help_sign($x), $x, util::help_sign($y), $y, diff_abs, diff_rel, $diff_allowed
            );
            // For the sign change check, allow (NAN vs NAN), but not (0.0 vs -0.0) or (NAN vs -NAN).
            assert!(
                $allow_sign_change || $x.is_sign_negative() == $y.is_sign_negative(),
                "assert_delta failed line {}, {}: {}{:e} vs {}{:e} sign difference disallowed", $line_num, $msg, util::help_sign($x), $x, util::help_sign($y), $y
            );
    }
}

// Check whether x and y are close enough to each other, considering both absolute and relative difference.
// Consider them close enough if absolute difference is within allow_abs inclusive,
// or if relative difference is within allow_rel exclusive. Passing a negative value
// for either allow_* value will exclude that check.
// nan values are considered equal to each other, and outside allowed range from everything else.
#[macro_export]
macro_rules! assert_delta_abs {
    ($x:expr, $y:expr, $allow_abs:expr, $allow_sign_change:expr, $msg:expr, $line_num:expr) => {
            let diff_abs = util::calc_delta_abs($x, $y);
            assert!(
                diff_abs <= $allow_abs,
                "assert_delta_abs failed line {}, {}: {}{:e} vs {}{:e} diff abs {:e} outside inclusive {:e}", $line_num, $msg, util::help_sign($x), $x, util::help_sign($y), $y, diff_abs, $allow_abs
            );
            // For the sign change check, allow (NAN vs NAN), but not (0.0 vs -0.0) or (NAN vs -NAN).
            assert!(
                $allow_sign_change || $x.is_sign_negative() == $y.is_sign_negative(),
                "assert_delta_abs failed line {}, {}: {}{:e} vs {}{:e} sign difference disallowed", $line_num, $msg, util::help_sign($x), $x, util::help_sign($y), $y
            );
    }
}

// Log a single comparison, using logic similar to DeltaEntry's
// logic for larger sets of comparisons.
// Note that a call to this function can typically be replaced with a
// call to assert_approx_eq!(x, y, allow_diff) for a non-logged version.
pub fn log_assert_delta(name_base: &str, name_detail: &str, x: f64, y: f64, allow_diff: f64, allow_sign_change: bool) {
    let diff_abs = util::calc_delta_abs(x, y);
    println!(
        "{} {}: {}{:e} vs {}{:e} abs {:e}, sign diff {}",
        name_base,
        name_detail,
        util::help_sign(x),
        x,
        util::help_sign(y),
        y,
        diff_abs,
        x.is_sign_negative() != y.is_sign_negative()
    );
    assert_delta_abs!(x, y, allow_diff, allow_sign_change, name_detail, 0);
}

