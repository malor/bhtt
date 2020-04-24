use std::borrow::Borrow;

use ordered_float::OrderedFloat;
use superslice::*;

use super::bin::Bin;

/// The Ben-Haim/Tom-Tov (BHTT) streaming histogram sketch implementation
/// (http://www.jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf).
///
/// A BHTT histogram is an ordered list of bins, that is a compact approximate
/// representation of numerical data distribution.
///
/// The desired size of the list is chosen at histogram creation time.
/// A bin is a (value, count) pair, where values are weighted averages
/// of numbers inserted to the histogram, and counts show how many numbers
/// have been merged together to produce a bin.
///
/// To insert a new value to the histogram, a new (value, 1) bin is added first,
/// preserving the ascending order of bins in the list. If the total number of
/// bins now exceeds the configured size value, then two closest bins are merged
/// together to restore the invariant.
///
/// Unlike traditional histograms, which require one to specify both the
/// histogram size and the bin boundaries configuration at creation time,
/// BHTT streaming histograms automatically adjust the bins as new values
/// are inserted.
///
/// Typical operations on the constructed histograms include approximations
/// of quantiles and counts.
///
/// BHTT histograms are mergeable, but it's *not* possible to subtract one
/// histogram from another bin by bin to get a difference between two.
///
/// Examples:
///
/// ```
/// use bhtt::Histogram;
///
/// let values = vec![1.0, 0.0, -5.4, -2.1, 8.5, 10.0, 8.6, 4.3, 7.8, 5.2];
///
/// // create a new Histogram and insert some values
/// let mut h = Histogram::new(5);
/// for value in values {
///     h.insert(value);
/// }
///
/// // histogram keeps track of the total count of inserted values,
/// // as well as the minimum/maximum inserted value
/// assert_eq!(h.size(), 5);
/// assert_eq!(h.count(), 10);
/// assert_eq!(h.min(), Some(-5.4));
/// assert_eq!(h.max(), Some(10.0));
///
/// // histogram allows to get approximations of arbitrary quantiles
/// assert_eq!(h.quantile(0.0), Some(-5.4));
/// assert_eq!(h.quantile(0.5), Some(4.75));
/// assert_eq!(h.quantile(1.0), Some(10.0));
///
/// // or estimate counts of values less than or equal to the given value
/// assert_eq!(h.count_less_than_or_equal_to(-7.4), 0);
/// assert_eq!(h.count_less_than_or_equal_to(5.0), 5);
/// assert_eq!(h.count_less_than_or_equal_to(13.0), 10);
///
/// // it's possible to merge histograms that were built separately
/// h.merge(&Histogram::from_iter(5, &[1.0, -7.6, 0.0, 5.8, 4.3, 2.1, 11.6]));
/// assert_eq!(h.count(), 17);
/// assert_eq!(h.min(), Some(-7.6));
/// assert_eq!(h.max(), Some(11.6));
/// ```
#[derive(Debug)]
pub struct Histogram {
    size: usize,
    bins: Vec<Bin>,
    min_value: Option<f64>,
    max_value: Option<f64>,
}

impl Histogram {
    /// Create a new Histogram with the given number of bins.
    ///
    /// The larger the size of the histogram, the more accurate approximations
    /// of quantiles one can get from it. But updates will also be slower and
    /// the histogram will consume more byte space.
    pub fn new(size: usize) -> Histogram {
        assert!(size > 0, "histogram size must be greater than 0");

        Histogram {
            size: size,
            // reserve one extra slot for bins, which are temporarily added during
            // histogram updates. This will allow us to avoid unnecessary memory
            // allocations
            bins: Vec::with_capacity(size + 1),
            min_value: None,
            max_value: None,
        }
    }

    /// Create a new Histogram from components of another histogram.
    pub fn from_parts(
        size: usize,
        bins: Vec<Bin>,
        min_value: Option<f64>,
        max_value: Option<f64>,
    ) -> Histogram {
        assert!(size > 0, "histogram size must be greater than 0");

        let mut h = Histogram {
            size: size,
            bins: bins,
            min_value: min_value,
            max_value: max_value,
        };
        h.shrink();
        h.bins.shrink_to_fit();

        h
    }

    /// Create a new Histogram of the given size from an iterator.
    pub fn from_iter<T>(size: usize, iter: T) -> Histogram
    where
        T: IntoIterator,
        T::Item: Borrow<f64>,
    {
        let mut h = Histogram::new(size);

        for v in iter.into_iter() {
            h.insert(*v.borrow());
        }

        h
    }

    /// Returns the size of the histogram.
    pub fn size(&self) -> usize {
        self.size
    }

    /// Returns the bins of the histogram.
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }

    /// Returns the total number of values inserted to the histogram.
    pub fn count(&self) -> u64 {
        self.bins.iter().map(|bin| bin.count()).sum()
    }

    /// Returns the minimum inserted value or None, if the histogram is empty.
    pub fn min(&self) -> Option<f64> {
        self.min_value
    }

    /// Returns the maximum inserted value or None, if the histogram is empty.
    pub fn max(&self) -> Option<f64> {
        self.max_value
    }

    /// Returns an approximated value of the `q`'th quantile of the inserted values
    /// or None, if the histogram is empty. `q` must be in the range [0.0; 1.0].
    pub fn quantile(&self, q: f64) -> Option<f64> {
        assert!(q >= 0.0 && q <= 1.0, "q must be in the range [0.0; 1.0]");

        if q == 0.0 {
            self.min()
        } else if q == 1.0 {
            self.max()
        } else {
            let total_count = self.count();
            match total_count {
                0 => None,
                _ => {
                    // Algorithm 4: Uniform procedure (from the paper mentioned in the description)
                    //
                    // To estimate the value of a given quantile we first find a pair of bins in
                    // the histogram that enclose the target cumulative count. Then we use the
                    // difference between the target count and the cumulative count up to the left
                    // bin in the pair to form a (quadratic) equation, that allows us to calculate
                    // the target value based on its proximity to the right bin.

                    let qth_count = self.count() as f64 * q;
                    let (i, up_to_qth_count) = self.index_of_cumulative_count_less_than(qth_count);

                    let (left_bin, right_bin) = self.get_bordering_bins(i);
                    let (left_value, left_count) = (left_bin.value(), left_bin.count() as f64);
                    let (right_value, right_count) = (right_bin.value(), right_bin.count() as f64);

                    let d = qth_count - up_to_qth_count;
                    let a = right_count - left_count;
                    if a == 0.0 {
                        Some(left_value + (right_value - left_value) * d / left_count)
                    } else {
                        let b = 2.0 * left_count;
                        let c = -2.0 * d;
                        let z = (-b + (b.powi(2) - 4.0 * a * c).sqrt()) / (2.0 * a);

                        Some(left_value + (right_value - left_value) * z)
                    }
                }
            }
        }
    }

    /// Returns an estimate of the number of values in the histogram that are less than or equal to
    /// `value`.
    pub fn count_less_than_or_equal_to(&self, value: f64) -> u64 {
        assert!(!value.is_nan(), "value must not be NaN");

        let total_count = self.count();
        if total_count == 0 || value < self.min().unwrap_or(std::f64::NAN) {
            // histogram is empty, or the interval (-inf; value] does not intersect
            // with the interval [min; max]
            0
        } else if value >= self.max().unwrap_or(std::f64::NAN) {
            // the interval (-inf; value] includes all the values in the histogram
            total_count
        } else {
            // Algorithm 3: Sum Procedure (from the paper mentioned in the description)
            //
            // In order to estimate the number of values in the histogram that are less than or
            // equal to the given value we need to find a pair of bins, which would be adjacent to
            // the (value, count) bin if we were to insert it to the histogram. The resulting count
            // will be equal to the sum of the following components:
            //
            // 1) sum of counts of the bins preceeding the left neighbour
            // 2) one half of left neighbour's count
            // 3) count of values between the left neighbour and the (value, count) bin

            // find the position of the bin if we were to insert it to the histogram
            let pos = self.bins.upper_bound(&Bin::new(value, 0));

            // calculate the sum of counts of the bins preceeding the left neighbour of that bin
            let left = pos.checked_sub(1).unwrap_or(0);
            let count_up_to_left: u64 = self.bins[..left].iter().map(|bin| bin.count()).sum();

            // determine the bordering bins
            let (left_bin, right_bin) = self.get_bordering_bins(pos);
            let (left_value, left_count) = (left_bin.value(), left_bin.count() as f64);
            let (right_value, right_count) = (right_bin.value(), right_bin.count() as f64);

            // estimate the count of values between the left neighbour and the (value, count) bins
            let count_left_to_value = if right_value - left_value <= 0.0 {
                0.0
            } else {
                let proximity_to_right = (value - left_value) / (right_value - left_value);
                let count = left_count + (right_count - left_count) * proximity_to_right;

                (left_count + count) / 2.0 * proximity_to_right
            };

            // add up all partial counts and round to the nearest integer number
            (count_up_to_left as f64 + left_count / 2.0 + count_left_to_value).round() as u64
        }
    }

    /// Update the histogram by inserting a new value.
    pub fn insert(&mut self, value: f64) {
        self.insert_bin(&Bin::new(value, 1));
    }

    /// Update the histogram by inserting a new bin.
    pub fn insert_bin(&mut self, bin: &Bin) {
        // insert the new bin preserving the ascending order. If the total number of bins exceeds
        // the configured size, the histogram is shrinked by merging two closest bins to restore
        // the invariant
        self.bins.insert(self.bins.upper_bound(bin), *bin);
        self.shrink();
        self.track_min_max(bin.value());
    }

    /// Merge the histogram with another one (in-place).
    pub fn merge(&mut self, other: &Histogram) {
        for bin in other.bins() {
            self.insert_bin(bin);
        }

        if let Some(min_value) = other.min() {
            self.track_min_max(min_value);
        }
        if let Some(max_value) = other.max() {
            self.track_min_max(max_value);
        }
    }

    /// Keep track of the minimum and the maximum inserted values
    /// (this will allow us to have more accurate quantile approximations).
    fn track_min_max(&mut self, value: f64) {
        self.min_value
            .replace(self.min_value.map_or(
                value,
                |current| if value < current { value } else { current },
            ));
        self.max_value
            .replace(self.max_value.map_or(
                value,
                |current| if value > current { value } else { current },
            ));
    }

    /// Merge two closest bins until the histogram shrinks back to the fixed size.
    fn shrink(&mut self) {
        while self.bins.len() > self.size {
            let (left, right) = self.find_closest_bins();
            self.bins[left] = Bin::merge(&self.bins[left], &self.bins[right]);
            self.bins.remove(right);
        }
    }

    /// Find a pair of bins that are closest to each other.
    fn find_closest_bins(&self) -> (usize, usize) {
        let right_index = (1..self.bins.len())
            .min_by_key(|i| {
                (
                    // distance between values is considered first
                    OrderedFloat((self.bins[*i].value() - self.bins[*i - 1].value()).abs()),
                    // if distances are equal, a pair of bins with smaller total count is preferred
                    self.bins[i - 1].count() + self.bins[*i].count(),
                )
            })
            .unwrap_or(self.bins.len() - 1);

        (right_index - 1, right_index)
    }

    fn index_of_cumulative_count_less_than(&self, target_count: f64) -> (usize, f64) {
        self.bins
            .iter()
            .zip(std::iter::once(&Bin::new(0.0, 0)).chain(&self.bins))
            .map(|(l, r)| (l.count() + r.count()) as f64 / 2.0)
            .scan(0.0, |partial_count, next_count| {
                *partial_count += next_count;
                Some(*partial_count)
            })
            .enumerate()
            .take_while(|(_, partial_count)| target_count > *partial_count)
            .last()
            .map_or((0, 0.0), |(i, sum)| (i + 1, sum))
    }

    fn get_bordering_bins(&self, i: usize) -> (Bin, Bin) {
        if i == 0 {
            (
                Bin::new(self.min_value.unwrap(), 0),
                *self.bins.first().unwrap(),
            )
        } else if i == self.bins.len() {
            (
                *self.bins.last().unwrap(),
                Bin::new(self.max_value.unwrap(), 0),
            )
        } else {
            (self.bins[i - 1], self.bins[i])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let h = Histogram::new(5);
        assert_eq!(h.size(), 5);
        assert_eq!(h.count(), 0);
        assert_eq!(h.min(), None);
        assert_eq!(h.max(), None);
        assert_eq!(h.bins(), &[]);
    }

    #[test]
    #[should_panic(expected = "histogram size must be greater than 0")]
    fn new_invalid_size() {
        Histogram::new(0);
    }

    #[test]
    fn from_parts() {
        let size = 5;
        let min_value = Some(-0.5);
        let max_value = Some(85.0);
        let bins = vec![Bin::new(0.0, 4), Bin::new(42.0, 8), Bin::new(84.0, 2)];

        let h = Histogram::from_parts(size, bins.clone(), min_value, max_value);
        assert_eq!(h.size(), size);
        assert_eq!(h.count(), 14);
        assert_eq!(h.min(), min_value);
        assert_eq!(h.max(), max_value);
        assert_eq!(h.bins(), bins.as_slice());
    }

    #[test]
    #[should_panic(expected = "histogram size must be greater than 0")]
    fn from_parts_invalid_size() {
        let size = 0;
        let bins = vec![Bin::new(0.0, 4), Bin::new(42.0, 8), Bin::new(84.0, 2)];

        Histogram::from_parts(size, bins.clone(), Some(0.0), Some(100.0));
    }

    #[test]
    fn insert() {
        let values = vec![
            1.0, 0.0, -5.4, -2.1, 8.5, 10.0, 8.6, 4.3, 7.8, 5.2, -6.0, -6.6, 0.5, 0.5, 2.625,
        ];
        let expected_bins = vec![
            Bin::new(-6.0, 3),
            Bin::new(-2.1, 1),
            Bin::new(0.5, 4),
            Bin::new(4.041666666666667, 3),
            Bin::new(8.725, 4),
        ];

        let mut h = Histogram::new(5);
        for v in &values {
            h.insert(*v);
        }

        assert_eq!(h.count(), values.len() as u64);
        assert_eq!(h.size(), 5);
        assert_eq!(h.min(), Some(-6.6));
        assert_eq!(h.max(), Some(10.0));
        assert_eq!(h.bins(), expected_bins.as_slice());
    }

    #[test]
    fn insert_single_value() {
        let mut h = Histogram::new(5);
        h.insert(42.0);

        assert_eq!(h.count(), 1);
        assert_eq!(h.size(), 5);
        assert_eq!(h.min(), Some(42.0));
        assert_eq!(h.max(), Some(42.0));
        assert_eq!(h.bins(), &[Bin::new(42.0, 1)]);
    }

    #[test]
    fn insert_bin() {
        let bins = vec![
            Bin::new(4.9, 6),
            Bin::new(5.0, 8),
            Bin::new(20.1, 7),
            Bin::new(4.0, 8),
            Bin::new(42.0, 14),
            Bin::new(17.4, 4),
            Bin::new(-10.0, 1),
        ];
        let expected_bins = vec![
            Bin::new(-10.0, 1),
            Bin::new(4.609090909090909, 22),
            Bin::new(17.4, 4),
            Bin::new(20.1, 7),
            Bin::new(42.0, 14),
        ];

        let mut h = Histogram::new(5);
        for bin in &bins {
            h.insert_bin(bin);
        }

        assert_eq!(h.count(), 48);
        assert_eq!(h.size(), 5);
        assert_eq!(h.min(), Some(-10.0));
        assert_eq!(h.max(), Some(42.0));
        assert_eq!(h.bins(), expected_bins.as_slice());
    }

    #[test]
    fn insert_single_bin() {
        let mut h = Histogram::new(5);
        h.insert_bin(&Bin::new(42.0, 84));

        assert_eq!(h.count(), 84);
        assert_eq!(h.size(), 5);
        assert_eq!(h.min(), Some(42.0));
        assert_eq!(h.max(), Some(42.0));
        assert_eq!(h.bins(), &[Bin::new(42.0, 84)]);
    }

    #[test]
    fn merge() {
        let bins1 = vec![
            Bin::new(-6.0, 3),
            Bin::new(-2.1, 1),
            Bin::new(0.5, 4),
            Bin::new(4.041666666666667, 3),
            Bin::new(8.725, 4),
        ];
        let mut h1 = Histogram::from_parts(5, bins1, Some(-6.6), Some(10.0));

        let bins2 = vec![
            Bin::new(33.32588794226721, 9977),
            Bin::new(1255.8137647058825, 17),
            Bin::new(3364.983, 2),
            Bin::new(5361.3435, 2),
            Bin::new(7349.9465, 2),
        ];
        let h2 = Histogram::from_parts(5, bins2, Some(9.48), Some(7829.851));

        let expected_bins = vec![
            Bin::new(33.27875390312249, 9992),
            Bin::new(1255.8137647058825, 17),
            Bin::new(3364.983, 2),
            Bin::new(5361.3435, 2),
            Bin::new(7349.9465, 2),
        ];

        h1.merge(&h2);

        assert_eq!(h1.count(), 10015);
        assert_eq!(h1.size(), 5);
        assert_eq!(h1.min(), Some(-6.6));
        assert_eq!(h1.max(), Some(7829.851));
        assert_eq!(h1.bins(), expected_bins.as_slice());
    }

    #[test]
    fn merge_empty() {
        let mut h1 = Histogram::new(5);
        let h2 = Histogram::new(10);

        h1.merge(&h2);

        assert_eq!(h1.count(), 0);
        assert_eq!(h1.size(), 5);
        assert_eq!(h1.min(), None);
        assert_eq!(h1.max(), None);
        assert_eq!(h1.bins(), &[]);
    }

    #[test]
    fn from_iter() {
        let values = vec![
            1.0, 0.0, -5.4, -2.1, 8.5, 10.0, 8.6, 4.3, 7.8, 5.2, -6.0, -6.6, 0.5, 0.5, 2.625,
        ];
        let expected_bins = vec![
            Bin::new(-6.0, 3),
            Bin::new(-2.1, 1),
            Bin::new(0.5, 4),
            Bin::new(4.041666666666667, 3),
            Bin::new(8.725, 4),
        ];

        let h = Histogram::from_iter(5, &values);

        assert_eq!(h.count(), values.len() as u64);
        assert_eq!(h.size(), 5);
        assert_eq!(h.min(), Some(-6.6));
        assert_eq!(h.max(), Some(10.0));
        assert_eq!(h.bins(), expected_bins.as_slice());
    }

    #[test]
    fn find_closest_bins_distance() {
        // proximity of bins is defined by the absolute distance between their values
        let bins = vec![
            Bin::new(-10.0, 1),
            Bin::new(4.9, 1),
            Bin::new(5.0, 1),
            Bin::new(17.4, 1),
            Bin::new(20.1, 1),
        ];
        let h = Histogram::from_parts(5, bins, Some(-10.0), Some(100.0));

        assert_eq!(h.find_closest_bins(), (1, 2));
    }

    #[test]
    fn find_closest_bins_count() {
        // distances between bin values are equal. Therefore, the pair of bins
        // with a smaller total count is considered to be the closest.
        let bins = vec![
            Bin::new(1.0, 1),
            Bin::new(2.0, 2),
            Bin::new(3.0, 3),
            Bin::new(4.0, 4),
            Bin::new(5.0, 5),
        ];
        let h = Histogram::from_parts(5, bins, Some(-10.0), Some(100.0));

        assert_eq!(h.find_closest_bins(), (0, 1));
    }

    #[test]
    fn find_closest_bins_index() {
        // distances between bin values are equal, and so are bin counts.
        // Therefore, the leftmost pair of bins is considered to be the closest.
        let bins = vec![
            Bin::new(1.0, 1),
            Bin::new(2.0, 1),
            Bin::new(3.0, 1),
            Bin::new(4.0, 1),
            Bin::new(5.0, 1),
        ];
        let h = Histogram::from_parts(5, bins, Some(-10.0), Some(100.0));

        assert_eq!(h.find_closest_bins(), (0, 1));
    }

    #[test]
    fn index_of_cumulative_count_less_than() {
        let bins = vec![
            Bin::new(1.0, 10),
            Bin::new(2.0, 8),
            Bin::new(3.0, 7),
            Bin::new(4.0, 15),
            Bin::new(5.0, 20),
        ];
        let h = Histogram::from_parts(5, bins, Some(0.0), Some(6.0));

        // [0] [ 5 : 5 ] [ 4 : 4 ] [ 3.5 : 3.5 ] [ 7.5 : 7.5 ] [ 10 : 10 ] [0]
        //    0         1         2             3             4          5

        assert_eq!(h.index_of_cumulative_count_less_than(-10.0), (0, 0.0));
        assert_eq!(h.index_of_cumulative_count_less_than(0.0), (0, 0.0));
        assert_eq!(h.index_of_cumulative_count_less_than(1.0), (0, 0.0));
        assert_eq!(h.index_of_cumulative_count_less_than(5.0), (0, 0.0));
        assert_eq!(h.index_of_cumulative_count_less_than(6.0), (1, 5.0));
        assert_eq!(h.index_of_cumulative_count_less_than(11.0), (1, 5.0));
        assert_eq!(h.index_of_cumulative_count_less_than(14.0), (1, 5.0));
        assert_eq!(h.index_of_cumulative_count_less_than(15.0), (2, 14.0));
        assert_eq!(h.index_of_cumulative_count_less_than(18.0), (2, 14.0));
        assert_eq!(h.index_of_cumulative_count_less_than(21.5), (2, 14.0));
        assert_eq!(h.index_of_cumulative_count_less_than(22.0), (3, 21.5));
        assert_eq!(h.index_of_cumulative_count_less_than(60.0), (5, 50.0));
        assert_eq!(h.index_of_cumulative_count_less_than(70.0), (5, 50.0));
    }

    #[test]
    fn quantile_empty() {
        let h = Histogram::new(5);

        assert_eq!(h.quantile(0.0), None);
        assert_eq!(h.quantile(0.5), None);
        assert_eq!(h.quantile(1.0), None);
    }

    #[test]
    #[should_panic(expected = "q must be in the range [0.0; 1.0]")]
    fn quantile_nan() {
        let h = Histogram::new(5);
        h.quantile(std::f64::NAN);
    }

    #[test]
    #[should_panic(expected = "q must be in the range [0.0; 1.0]")]
    fn quantile_not_in_range_high() {
        let h = Histogram::new(5);
        h.quantile(1.1);
    }

    #[test]
    #[should_panic(expected = "q must be in the range [0.0; 1.0]")]
    fn quantile_not_in_range_low() {
        let h = Histogram::new(5);
        h.quantile(-0.1);
    }

    #[test]
    fn quantile() {
        let bins = vec![
            Bin::new(2.0, 1),
            Bin::new(9.5, 2),
            Bin::new(19.33, 3),
            Bin::new(32.67, 3),
            Bin::new(45.0, 1),
        ];
        let h = Histogram::from_parts(5, bins, Some(2.0), Some(45.0));

        let expected = vec![
            // quantile, expected value, max relative difference
            (0.0, 2.0, 0.0),
            (0.1, 8.3, 0.4),
            (0.25, 11.5, 0.1),
            (0.5, 21.0, 0.1),
            (0.75, 31.5, 0.1),
            (0.9, 36.9, 0.1),
            (0.99, 44.19, 0.1),
            (1.0, 45.0, 0.0),
        ];
        for (q, expected_value, max_relative) in expected {
            assert_relative_eq!(
                h.quantile(q).unwrap(),
                expected_value,
                max_relative = max_relative
            );
        }
    }

    #[test]
    fn count_less_than_or_equal_to_empty() {
        let h = Histogram::new(5);

        assert_eq!(h.count_less_than_or_equal_to(-42.0), 0);
        assert_eq!(h.count_less_than_or_equal_to(0.0), 0);
        assert_eq!(h.count_less_than_or_equal_to(42.0), 0);
        assert_eq!(h.count_less_than_or_equal_to(std::f64::NEG_INFINITY), 0);
        assert_eq!(h.count_less_than_or_equal_to(std::f64::INFINITY), 0);
    }

    #[test]
    #[should_panic(expected = "value must not be NaN")]
    fn count_less_than_or_equal_to_nan() {
        let h = Histogram::new(5);
        h.count_less_than_or_equal_to(std::f64::NAN);
    }

    #[test]
    fn count_less_than_or_equal_to() {
        let bins = vec![
            Bin::new(2.0, 1),
            Bin::new(9.5, 2),
            Bin::new(19.33, 3),
            Bin::new(32.67, 3),
            Bin::new(45.0, 1),
        ];
        let h = Histogram::from_parts(5, bins, Some(2.0), Some(45.0));

        assert_eq!(h.count_less_than_or_equal_to(std::f64::NEG_INFINITY), 0);
        assert_eq!(h.count_less_than_or_equal_to(-42.0), 0);
        assert_eq!(h.count_less_than_or_equal_to(0.0), 0);
        assert_eq!(h.count_less_than_or_equal_to(2.1), 1);
        assert_eq!(h.count_less_than_or_equal_to(10.0), 2);
        assert_eq!(h.count_less_than_or_equal_to(15.0), 3);
        assert_eq!(h.count_less_than_or_equal_to(25.0), 6);
        assert_eq!(h.count_less_than_or_equal_to(38.0), 9);
        assert_eq!(h.count_less_than_or_equal_to(45.0), 10);
        assert_eq!(h.count_less_than_or_equal_to(std::f64::INFINITY), 10);
    }
}
