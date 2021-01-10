use ordered_float::NotNan;

/// Histogram bin stored as a (value, count) pair.
///
/// ```
/// use bhtt::Bin;
///
/// // Create a Bin with the value of 42.0 and the count of 84
/// let b1 = Bin::new(42.0, 84);
/// assert_eq!(b1.value(), 42.0);
/// assert_eq!(b1.count(), 84);
///
/// // Bins can be merged together. The value of the new Bin
/// // will be equal to the weighted average of the two being merged
/// let b2 = Bin::new(84.0, 42);
/// let b3 = Bin::merge(&b1, &b2);
/// assert_eq!(b3.value(), 56.0);
/// assert_eq!(b3.count(), 126);
///
/// // Bins have natural ordering: values are compared first,
/// // and counts are used to resolve the ties. Two bins are
/// // equal when both their values and their counts are equal
/// let reference = Bin::new(42.0, 84);
/// let gt_by_value = Bin::new(42.1, 84);
/// let gt_by_count = Bin::new(42.0, 85);
/// let equal = Bin::new(42.0, 84);
/// assert!(gt_by_value > reference);
/// assert!(reference < gt_by_value);
/// assert!(gt_by_count > reference);
/// assert!(reference < gt_by_count);
/// assert_eq!(reference, equal);
/// assert_eq!(equal, reference);
/// ```
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Eq, Ord)]
pub struct Bin {
    value: NotNan<f64>,
    count: u64,
}

impl Bin {
    /// Returns a new Bin with the given value and count
    ///
    /// # Arguments
    ///
    /// * `value` - Bin value. Must be a finite non-NaN value; otherwise a panic will be triggered.
    /// * `count` - Bin count. Must be greater than or equal to zero; otherwise a panic will be triggered.
    pub fn new(value: f64, count: u64) -> Bin {
        assert!(!value.is_nan(), "value must not be NaN");
        assert!(value.is_finite(), "value must be finite");

        Bin {
            value: NotNan::new(value).unwrap(),
            count: count,
        }
    }

    /// Returns a new Bin that is an approximation of two bins merged together.
    pub fn merge(left: &Bin, right: &Bin) -> Bin {
        let count = left.count() + right.count();
        assert!(count > 0, "count must be greater than zero");

        let value = (left.value() * left.count() as f64 + right.value() * right.count() as f64)
            / count as f64;

        Bin::new(value, count)
    }

    /// Returns the value of the bin
    pub fn value(&self) -> f64 {
        self.value.into_inner()
    }

    /// Returns the count of the bin
    pub fn count(&self) -> u64 {
        self.count
    }
}

impl From<f32> for Bin {
    fn from(value: f32) -> Self {
        Bin::new(value as f64, 1)
    }
}
impl From<f64> for Bin {
    fn from(value: f64) -> Self {
        Bin::new(value, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let b1 = Bin::new(42.0, 84);
        assert_eq!(b1.value(), 42.0);
        assert_eq!(b1.count(), 84);

        let b2 = Bin::new(-42.0, 0);
        assert_eq!(b2.value(), -42.0);
        assert_eq!(b2.count(), 0);
    }

    #[test]
    #[should_panic(expected = "value must not be NaN")]
    fn new_invalid_value_nan() {
        Bin::new(std::f64::NAN, 84);
    }

    #[test]
    #[should_panic(expected = "value must be finite")]
    fn new_invalid_value_positive_infinity() {
        Bin::new(std::f64::INFINITY, 84);
    }

    #[test]
    #[should_panic(expected = "value must be finite")]
    fn new_invalid_value_negative_infinity() {
        Bin::new(std::f64::NEG_INFINITY, 84);
    }

    #[test]
    fn ordering() {
        let reference = Bin::new(42.0, 84);

        let gt_by_value = Bin::new(42.1, 84);
        assert!(gt_by_value > reference);
        assert!(gt_by_value >= reference);
        assert!(reference < gt_by_value);
        assert!(reference <= gt_by_value);

        let gt_by_count = Bin::new(42.0, 85);
        assert!(gt_by_count > reference);
        assert!(gt_by_count >= reference);
        assert!(reference < gt_by_count);
        assert!(reference <= gt_by_count);

        let equal = Bin::new(42.0, 84);
        assert!(!(reference < equal));
        assert!(!(reference > equal));
        assert!(reference >= equal);
        assert!(reference <= equal);
    }

    #[test]
    fn equality() {
        let reference = Bin::new(42.0, 84);
        let equal = Bin::new(42.0, 84);
        let not_equal_by_value = Bin::new(42.1, 84);
        let not_equal_by_count = Bin::new(42.0, 85);

        assert_eq!(reference, equal);
        assert_eq!(equal, reference);
        assert_ne!(reference, not_equal_by_value);
        assert_ne!(not_equal_by_value, reference);
        assert_ne!(reference, not_equal_by_count);
        assert_ne!(not_equal_by_count, reference);
    }

    #[test]
    fn merge() {
        let left = Bin::new(42.0, 84);
        let right = Bin::new(84.0, 42);
        let expected = Bin::new(
            (42.0 * 84 as f64 + 84.0 * 42 as f64) / (84 + 42) as f64,
            84 + 42,
        );

        let actual = Bin::merge(&left, &right);
        assert_eq!(actual, expected);
    }

    #[test]
    #[should_panic(expected = "count must be greater than zero")]
    fn merge_invalid_count() {
        let left = Bin::new(42.0, 0);
        let right = Bin::new(84.0, 0);

        Bin::merge(&left, &right);
    }

    #[test]
    fn from() {
        let b1 = Bin::from(42.0f32);
        assert_eq!(b1.value(), 42.0);
        assert_eq!(b1.count(), 1);

        let b2 = Bin::from(-7.5f32);
        assert_eq!(b2.value(), -7.5);
        assert_eq!(b2.count(), 1);
    }

    #[test]
    #[should_panic(expected = "value must not be NaN")]
    fn from_nan() {
        Bin::from(std::f64::NAN);
    }

    #[test]
    #[should_panic(expected = "value must be finite")]
    fn from_infinite() {
        Bin::from(std::f64::INFINITY);
    }
}
