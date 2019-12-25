# bhtt

The Ben-Haim/Tom-Tov (BHTT) streaming histogram sketch implementation
(http://www.jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf).

This implementation was heavily inspired by https://github.com/aaw/histk
and mainly is a fun way to experiment with Rust.

A BHTT histogram is an ordered list of bins, that is a compact approximate
representation of numerical data distribution.

The desired size of the list is chosen at histogram creation time.
A bin is a (value, count) pair, where values are weighted averages
of numbers inserted to the histogram, and counts show how many numbers
have been merged together to produce a bin.

To insert a new value to the histogram, a new (value, 1) bin is added first,
preserving the ascending order of bins in the list. If the total number of
bins now exceeds the configured size value, then two closest bins are merged
together to restore the invariant.

Unlike traditional histograms, which require one to specify both the
histogram size and the bin boundaries configuration at creation time,
BHTT streaming histograms automatically adjust the bins as new values
are inserted.

Typical operations on the constructed histograms include approximations
of quantiles and counts.

BHTT histograms are mergeable, but it's *not* possible to subtract one
histogram from another bin by bin to get a difference between two.

## Examples

```rust
use bhtt::histogram::Histogram;

let values = vec![1.0, 0.0, -5.4, -2.1, 8.5, 10.0, 8.6, 4.3, 7.8, 5.2];

// create a new Histogram and insert some values
let mut h = Histogram::new(5);
for value in values {
    h.insert(value);
}

// histogram keeps track of the total count of inserted values,
// as well as the minimum/maximum inserted value
assert_eq!(h.size(), 5);
assert_eq!(h.count(), 10);
assert_eq!(h.min(), Some(-5.4));
assert_eq!(h.max(), Some(10.0));

// histogram allows to get approximations of arbitrary quantiles
assert_eq!(h.quantile(0.0), Some(-5.4));
assert_eq!(h.quantile(0.5), Some(4.75));
assert_eq!(h.quantile(1.0), Some(10.0));

// or estimated counts of values less than or equal to the given value
assert_eq!(h.count_less_than_or_equal_to(-7.4), 0);
assert_eq!(h.count_less_than_or_equal_to(5.0), 5);
assert_eq!(h.count_less_than_or_equal_to(13.0), 10);
```

## Development

### Running tests

Unit tests are put directly to the code modules in `src/`.

Integration tests are in `tests/` directory with input data and
utility functions stored in `utilities/` (so that they can be shared
between integration tests and benchmarks).

Run the full test suite:

```shell
$ cargo test
```

### Running benchmarks

Run all benchmarks:

```shell
$ cargo bench
```
