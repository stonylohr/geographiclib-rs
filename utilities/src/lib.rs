// todo: add difference calculations from tests/GeodTest to rs utilities (and review vs geod_error from validate-geographiclib
// todo: experiment with using diff::diff_ulps in some places
// todo: check uncommitted changes to make sure I didn't accidentally update any tolerances through copy-paste error
// todo: check uncommitted changes to make sure I didn't change any param names by copy-paste error
// todo: check uncommitted changes to make sure I didn't change any param indices by copy-paste error
// todo: change dependency on float-diff to specify version
// todo: change float-diff dependency back to github-based dependency

// testing tip: If you want full test output in a predictable order, including ignored tests,
//              use "cargo test -- --include-ignored --nocapture --test-threads 1"
//              (Note: include-ignored is still only available in nightly as of 2021/01/15, but looks like it's coming to stable soon)
// testing tip: If you want full test output in a predictable order, including ONLY ignored tests,
//              use "cargo test -- --ignored --nocapture --test-threads 1"
// testing tip: If a new "*_vs_cpp_*" test fails especially badly,
//              review the C++ "instrumented-crude" logging for possible bugs there
// testing tip: While assert variations are valuable for catching regressions, they're
//              often of more limited value when trying to _improve_ result precision.
//              In such cases, use of something like DiffSummary is often helpful,
//              but remember to use nocapture, and often directing output to a file, for
//              syntax like "cargo test -- --nocapture > ../test-out-a.txt".
// benchmarking reminder: Remember to shut down unneeded processes before benching.
// benchmarking reminder: In some cases, it may be desirable to capture multiple before and
//                        after bench runs, to get a clearer sense of natural variability.

pub mod util;

#[allow(non_upper_case_globals)]
pub const nC_: usize = 7; // todo: define in geodesic.rs instead
