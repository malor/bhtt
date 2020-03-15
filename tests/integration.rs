#[macro_use]
extern crate approx;

use test_case::test_case;

use bhtt::Histogram;

// dataset, histogram size, expected diff between true and approximated quantiles
#[test_case("utilities/testdata/pings.txt", 32, 0.5)]
#[test_case("utilities/testdata/pings.txt", 64, 0.25)]
#[test_case("utilities/testdata/pings.txt", 128, 0.035)]
#[test_case("utilities/testdata/pings.txt", 256, 0.005)]
fn quantile(filename: &str, histogram_size: usize, max_error_pct: f64) {
    let dataset = utilities::Dataset::from_file(filename).unwrap();

    let h = Histogram::from_iter(histogram_size, dataset.values());
    for (q, expected_value) in dataset.quantiles() {
        assert_relative_eq!(
            h.quantile(**q).unwrap(),
            expected_value,
            max_relative = max_error_pct
        );
    }
}
