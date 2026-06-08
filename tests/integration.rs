use googletest::prelude::*;
use test_case::test_case;

use bhtt::Histogram;

// dataset, histogram size, expected diff between true and approximated quantiles
#[test_case("utilities/testdata/pings.txt", 32, 0.5)]
#[test_case("utilities/testdata/pings.txt", 64, 0.25)]
#[test_case("utilities/testdata/pings.txt", 128, 0.035)]
#[test_case("utilities/testdata/pings.txt", 256, 0.005)]
fn quantile(filename: &str, histogram_size: usize, max_error_ratio: f64) {
    let dataset = utilities::Dataset::from_file(filename).unwrap();

    let h = Histogram::from_iter(histogram_size, dataset.values());
    for (q, expected_value) in dataset.quantiles() {
        let actual = h.quantile(**q).unwrap();
        assert_that!(actual, not(is_nan()));
        assert_that!(actual, is_finite());

        let tolerance = max_error_ratio * expected_value.abs().max(actual.abs());
        assert_that!(actual, near(*expected_value, tolerance));
    }
}
