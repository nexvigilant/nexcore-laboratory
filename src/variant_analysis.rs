//! Variant frequency distribution analysis — "Analyze" step of the lab pipeline.
//!
//! Given a population of labelled variants, this module computes the full
//! frequency distribution and a suite of information-theoretic statistics:
//! Shannon entropy, uniformity, chi-square goodness-of-fit, and a ranked
//! dispatch table with cumulative proportions.
//!
//! ## T1 Grounding
//!
//! | Primitive | Role |
//! |-----------|------|
//! | Σ (Sum) | accumulating counts and proportions |
//! | N (Quantity) | raw variant counts |
//! | κ (Comparison) | uniformity and chi-square ranking |
//! | σ (Sequence) | ordered dispatch table |
//!
//! ## PV Transfer
//!
//! ADR outcome type distribution — which reaction types dominate? Is the
//! distribution shifting across reporting periods? A heavily skewed
//! distribution (low uniformity, high chi-square) signals an emerging safety
//! signal concentrated in one outcome category.
//!
//! ## Example
//!
//! ```rust
//! use nexcore_laboratory::variant_analysis::{VariantDistribution, variant_entropy};
//!
//! let dist = VariantDistribution::from_counts(&[
//!     ("Hepatotoxicity", 50),
//!     ("Nephrotoxicity",  30),
//!     ("Cardiotoxicity",  20),
//! ])
//! .unwrap();
//!
//! let h = variant_entropy(&dist).unwrap();
//! assert!(h > 0.0 && h < 2.0);
//! ```
//!
//! ## Tier: T2-C (μ + Σ + κ + × + σ)

use serde::{Deserialize, Serialize};

// ============================================================================
// Error type
// ============================================================================

/// Errors that can occur during variant distribution construction or analysis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantError {
    /// The distribution contains no variants.
    EmptyDistribution,
    /// The total count supplied (or derived) is zero.
    ZeroTotal,
    /// The sum of individual frequencies does not match the recorded total.
    InconsistentTotal,
}

impl std::fmt::Display for VariantError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyDistribution => write!(f, "variant distribution is empty"),
            Self::ZeroTotal => write!(f, "total count is zero"),
            Self::InconsistentTotal => {
                write!(f, "sum of frequencies does not match total")
            }
        }
    }
}

impl std::error::Error for VariantError {}

// ============================================================================
// Core types
// ============================================================================

/// A single variant with its absolute count and relative proportion.
///
/// Tier: T2-P (N Quantity + μ Mapping)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VariantFrequency {
    /// Human-readable variant name (e.g. `"Hepatotoxicity"`).
    pub name: String,
    /// Absolute observed count.
    pub count: u64,
    /// Relative proportion in `[0.0, 1.0]`; `count / total`.
    pub proportion: f64,
}

/// Full frequency distribution over a set of named variants.
///
/// Constructed via [`VariantDistribution::from_counts`] which computes
/// proportions automatically.
///
/// Tier: T2-C (Σ Sum + N Quantity + μ Mapping)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VariantDistribution {
    /// Individual variant frequencies, in insertion order.
    pub frequencies: Vec<VariantFrequency>,
    /// Total count across all variants.
    pub total: u64,
}

/// A single row of the ranked dispatch table.
///
/// Tier: T2-P (σ Sequence + κ Comparison + N Quantity)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DispatchEntry {
    /// Variant name.
    pub name: String,
    /// Absolute count.
    pub count: u64,
    /// Relative proportion.
    pub proportion: f64,
    /// Cumulative proportion up to and including this entry
    /// (sum of proportions of all entries ranked at or above this one).
    pub cumulative: f64,
}

/// Ranked dispatch table derived from a [`VariantDistribution`].
///
/// Entries are sorted by count (descending). Includes aggregate
/// information-theoretic statistics.
///
/// Tier: T2-C (σ + Σ + κ + N)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DispatchTable {
    /// Ranked entries, highest count first.
    pub entries: Vec<DispatchEntry>,
    /// Shannon entropy of the distribution in bits.
    pub entropy: f64,
    /// Uniformity in `[0.0, 1.0]`; 1.0 = perfectly uniform.
    pub uniformity: f64,
}

// ============================================================================
// VariantDistribution construction
// ============================================================================

impl VariantDistribution {
    /// Build a distribution from a slice of `(name, count)` pairs.
    ///
    /// Proportions are computed as `count / total` where `total` is the sum
    /// of all counts.  Returns [`VariantError::EmptyDistribution`] when the
    /// slice is empty and [`VariantError::ZeroTotal`] when every count is
    /// zero.
    ///
    /// # Errors
    ///
    /// - [`VariantError::EmptyDistribution`] — `counts` slice is empty.
    /// - [`VariantError::ZeroTotal`] — all counts are zero (total = 0).
    ///
    /// # Example
    ///
    /// ```rust
    /// use nexcore_laboratory::variant_analysis::VariantDistribution;
    ///
    /// let dist = VariantDistribution::from_counts(&[
    ///     ("Alpha", 3),
    ///     ("Beta",  1),
    /// ])
    /// .unwrap();
    ///
    /// assert_eq!(dist.total, 4);
    /// assert!((dist.frequencies[0].proportion - 0.75).abs() < 1e-9);
    /// ```
    pub fn from_counts(counts: &[(&str, u64)]) -> Result<Self, VariantError> {
        if counts.is_empty() {
            return Err(VariantError::EmptyDistribution);
        }

        let total: u64 = counts.iter().map(|&(_, c)| c).sum();

        if total == 0 {
            return Err(VariantError::ZeroTotal);
        }

        let total_f = total as f64;
        let frequencies = counts
            .iter()
            .map(|&(name, count)| VariantFrequency {
                name: name.to_string(),
                count,
                proportion: count as f64 / total_f,
            })
            .collect();

        Ok(Self { frequencies, total })
    }
}

// ============================================================================
// Analysis functions
// ============================================================================

/// Compute the Shannon entropy of a variant distribution in bits.
///
/// `H = -Σ p(x) * log₂(p(x))` for each variant where `p(x) > 0`.
/// Zero-proportion variants contribute 0 (by the convention `0 * log(0) = 0`).
///
/// # Errors
///
/// - [`VariantError::EmptyDistribution`] — distribution has no variants.
/// - [`VariantError::ZeroTotal`] — total count is zero.
///
/// # Example
///
/// ```rust
/// use nexcore_laboratory::variant_analysis::{VariantDistribution, variant_entropy};
///
/// // Two equally probable variants → entropy = 1 bit
/// let dist = VariantDistribution::from_counts(&[("A", 1), ("B", 1)]).unwrap();
/// let h = variant_entropy(&dist).unwrap();
/// assert!((h - 1.0).abs() < 1e-9);
/// ```
pub fn variant_entropy(dist: &VariantDistribution) -> Result<f64, VariantError> {
    validate_non_empty(dist)?;

    let h = dist
        .frequencies
        .iter()
        .filter(|vf| vf.proportion > 0.0)
        .map(|vf| -vf.proportion * vf.proportion.log2())
        .sum();

    Ok(h)
}

/// Compute the uniformity of a variant distribution.
///
/// `uniformity = H / H_max` where `H_max = log₂(n)` for `n` variants.
/// Returns a value in `[0.0, 1.0]` where `1.0` means perfectly uniform.
///
/// A single-variant distribution is trivially uniform (`uniformity = 1.0`).
///
/// # Errors
///
/// - [`VariantError::EmptyDistribution`] — distribution has no variants.
/// - [`VariantError::ZeroTotal`] — total count is zero.
///
/// # Example
///
/// ```rust
/// use nexcore_laboratory::variant_analysis::{VariantDistribution, variant_uniformity};
///
/// // Single variant is trivially uniform
/// let dist = VariantDistribution::from_counts(&[("Only", 7)]).unwrap();
/// let u = variant_uniformity(&dist).unwrap();
/// assert!((u - 1.0).abs() < 1e-9);
/// ```
pub fn variant_uniformity(dist: &VariantDistribution) -> Result<f64, VariantError> {
    validate_non_empty(dist)?;

    let n = dist.frequencies.len();

    // Single variant is trivially uniform.
    if n == 1 {
        return Ok(1.0);
    }

    let h_max = (n as f64).log2();
    let h = variant_entropy(dist)?;

    Ok(h / h_max)
}

/// Return a reference to the variant with the highest count.
///
/// If multiple variants share the maximum count the first one encountered
/// (insertion order) is returned.  Returns `None` when the distribution
/// contains no variants.
///
/// # Example
///
/// ```rust
/// use nexcore_laboratory::variant_analysis::{VariantDistribution, dominant_variant};
///
/// let dist = VariantDistribution::from_counts(&[
///     ("Minor", 10),
///     ("Major", 90),
/// ])
/// .unwrap();
///
/// let dom = dominant_variant(&dist).unwrap();
/// assert_eq!(dom.name, "Major");
/// ```
#[must_use]
pub fn dominant_variant(dist: &VariantDistribution) -> Option<&VariantFrequency> {
    // Use a manual fold so that ties are broken by insertion order (first wins).
    // `Iterator::max_by_key` returns the *last* maximum on ties, which is not
    // the desired behaviour.
    dist.frequencies.iter().reduce(|best, current| {
        if current.count > best.count {
            current
        } else {
            best
        }
    })
}

/// Build a ranked dispatch table from a variant distribution.
///
/// Entries are sorted by count (descending, stable).  Cumulative proportions
/// are computed in that sorted order so the final entry's `cumulative` reaches
/// `1.0` (within floating-point precision).
///
/// # Errors
///
/// - [`VariantError::EmptyDistribution`] — distribution has no variants.
/// - [`VariantError::ZeroTotal`] — total count is zero.
///
/// # Example
///
/// ```rust
/// use nexcore_laboratory::variant_analysis::{VariantDistribution, dispatch_table};
///
/// let dist = VariantDistribution::from_counts(&[
///     ("Rare",   5),
///     ("Common", 95),
/// ])
/// .unwrap();
///
/// let table = dispatch_table(&dist).unwrap();
/// assert_eq!(table.entries[0].name, "Common");
/// assert!((table.entries.last().unwrap().cumulative - 1.0).abs() < 1e-9);
/// ```
pub fn dispatch_table(dist: &VariantDistribution) -> Result<DispatchTable, VariantError> {
    validate_non_empty(dist)?;

    let entropy = variant_entropy(dist)?;
    let uniformity = variant_uniformity(dist)?;

    // Clone and sort by count descending (stable preserves tie-break order).
    let mut sorted: Vec<&VariantFrequency> = dist.frequencies.iter().collect();
    sorted.sort_by(|a, b| b.count.cmp(&a.count));

    let mut cumulative = 0.0_f64;
    let entries = sorted
        .iter()
        .map(|vf| {
            cumulative += vf.proportion;
            DispatchEntry {
                name: vf.name.clone(),
                count: vf.count,
                proportion: vf.proportion,
                cumulative,
            }
        })
        .collect();

    Ok(DispatchTable {
        entries,
        entropy,
        uniformity,
    })
}

/// Compute the chi-square statistic for a uniformity test.
///
/// Compares the observed distribution to the expected uniform distribution:
///
/// `expected = total / n_variants` (identical for every variant)
///
/// `χ² = Σ (observed - expected)² / expected`
///
/// A value near zero indicates the distribution is close to uniform.
/// Large values indicate significant skew.
///
/// # Errors
///
/// - [`VariantError::EmptyDistribution`] — distribution has no variants.
/// - [`VariantError::ZeroTotal`] — total count is zero.
///
/// # Example
///
/// ```rust
/// use nexcore_laboratory::variant_analysis::{VariantDistribution, chi_square_uniformity};
///
/// // Perfectly uniform — chi-square should be ~0.
/// let dist = VariantDistribution::from_counts(&[
///     ("A", 25), ("B", 25), ("C", 25), ("D", 25),
/// ])
/// .unwrap();
///
/// let chi2 = chi_square_uniformity(&dist).unwrap();
/// assert!(chi2 < 1e-9);
/// ```
pub fn chi_square_uniformity(dist: &VariantDistribution) -> Result<f64, VariantError> {
    validate_non_empty(dist)?;

    let n = dist.frequencies.len() as f64;
    let expected = dist.total as f64 / n;

    let chi2 = dist
        .frequencies
        .iter()
        .map(|vf| {
            let diff = vf.count as f64 - expected;
            diff * diff / expected
        })
        .sum();

    Ok(chi2)
}

// ============================================================================
// Internal helpers
// ============================================================================

/// Validate that a distribution is non-empty and has a positive total.
fn validate_non_empty(dist: &VariantDistribution) -> Result<(), VariantError> {
    if dist.frequencies.is_empty() {
        return Err(VariantError::EmptyDistribution);
    }
    if dist.total == 0 {
        return Err(VariantError::ZeroTotal);
    }
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    /// Floating-point equality with a generous tolerance for f64 arithmetic.
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ------------------------------------------------------------------
    // VariantDistribution::from_counts — construction
    // ------------------------------------------------------------------

    #[test]
    fn from_counts_basic_construction() {
        let dist = VariantDistribution::from_counts(&[("Alpha", 3), ("Beta", 1)]).unwrap();
        assert_eq!(dist.total, 4);
        assert_eq!(dist.frequencies.len(), 2);
        assert_eq!(dist.frequencies[0].name, "Alpha");
        assert_eq!(dist.frequencies[0].count, 3);
        assert_eq!(dist.frequencies[1].name, "Beta");
        assert_eq!(dist.frequencies[1].count, 1);
    }

    #[test]
    fn from_counts_proportions_are_correct() {
        let dist = VariantDistribution::from_counts(&[("X", 75), ("Y", 25)]).unwrap();
        assert!(approx_eq(dist.frequencies[0].proportion, 0.75, 1e-9));
        assert!(approx_eq(dist.frequencies[1].proportion, 0.25, 1e-9));
    }

    #[test]
    fn from_counts_proportions_sum_to_one() {
        let dist = VariantDistribution::from_counts(&[("A", 10), ("B", 30), ("C", 60)]).unwrap();
        let sum: f64 = dist.frequencies.iter().map(|vf| vf.proportion).sum();
        assert!(approx_eq(sum, 1.0, 1e-9));
    }

    #[test]
    fn from_counts_empty_returns_error() {
        let result = VariantDistribution::from_counts(&[]);
        assert_eq!(result, Err(VariantError::EmptyDistribution));
    }

    #[test]
    fn from_counts_all_zero_returns_error() {
        let result = VariantDistribution::from_counts(&[("A", 0), ("B", 0)]);
        assert_eq!(result, Err(VariantError::ZeroTotal));
    }

    #[test]
    fn from_counts_single_variant_proportion_is_one() {
        let dist = VariantDistribution::from_counts(&[("Solo", 42)]).unwrap();
        assert!(approx_eq(dist.frequencies[0].proportion, 1.0, 1e-9));
    }

    // ------------------------------------------------------------------
    // variant_entropy
    // ------------------------------------------------------------------

    #[test]
    fn entropy_single_variant_is_zero() {
        // p=1 → H = -(1 * log2(1)) = 0
        let dist = VariantDistribution::from_counts(&[("Only", 100)]).unwrap();
        let h = variant_entropy(&dist).unwrap();
        assert!(approx_eq(h, 0.0, 1e-9));
    }

    #[test]
    fn entropy_two_equal_variants_is_one_bit() {
        let dist = VariantDistribution::from_counts(&[("A", 50), ("B", 50)]).unwrap();
        let h = variant_entropy(&dist).unwrap();
        assert!(approx_eq(h, 1.0, 1e-9));
    }

    #[test]
    fn entropy_two_unequal_variants_less_than_one() {
        // 75/25 split — H < 1 bit
        let dist = VariantDistribution::from_counts(&[("A", 75), ("B", 25)]).unwrap();
        let h = variant_entropy(&dist).unwrap();
        // Hand-calculated: -(0.75*log2(0.75) + 0.25*log2(0.25))
        // = -(0.75*-0.41504 + 0.25*-2.0) = 0.81128 bits
        assert!(h < 1.0);
        assert!(h > 0.0);
        assert!(approx_eq(h, 0.811_278_124_459_133_6, 1e-9));
    }

    #[test]
    fn entropy_four_uniform_variants_is_two_bits() {
        // log2(4) = 2
        let dist = VariantDistribution::from_counts(&[("A", 25), ("B", 25), ("C", 25), ("D", 25)])
            .unwrap();
        let h = variant_entropy(&dist).unwrap();
        assert!(approx_eq(h, 2.0, 1e-9));
    }

    #[test]
    fn entropy_heavily_skewed_distribution() {
        // 97 in one bucket, 1 each in three others → near-zero entropy
        let dist = VariantDistribution::from_counts(&[
            ("Major", 97),
            ("Minor1", 1),
            ("Minor2", 1),
            ("Minor3", 1),
        ])
        .unwrap();
        let h = variant_entropy(&dist).unwrap();
        // Entropy exists but is much less than log2(4)=2
        assert!(h > 0.0);
        assert!(h < 0.5);
    }

    #[test]
    fn entropy_empty_distribution_returns_error() {
        // Manually construct an empty distribution to exercise the guard.
        let dist = VariantDistribution {
            frequencies: vec![],
            total: 0,
        };
        assert_eq!(variant_entropy(&dist), Err(VariantError::EmptyDistribution));
    }

    // ------------------------------------------------------------------
    // variant_uniformity
    // ------------------------------------------------------------------

    #[test]
    fn uniformity_single_variant_is_one() {
        let dist = VariantDistribution::from_counts(&[("X", 10)]).unwrap();
        let u = variant_uniformity(&dist).unwrap();
        assert!(approx_eq(u, 1.0, 1e-9));
    }

    #[test]
    fn uniformity_two_equal_variants_is_one() {
        let dist = VariantDistribution::from_counts(&[("A", 1), ("B", 1)]).unwrap();
        let u = variant_uniformity(&dist).unwrap();
        assert!(approx_eq(u, 1.0, 1e-9));
    }

    #[test]
    fn uniformity_four_equal_variants_is_one() {
        let dist = VariantDistribution::from_counts(&[("A", 25), ("B", 25), ("C", 25), ("D", 25)])
            .unwrap();
        let u = variant_uniformity(&dist).unwrap();
        assert!(approx_eq(u, 1.0, 1e-9));
    }

    #[test]
    fn uniformity_skewed_distribution_is_less_than_one() {
        let dist = VariantDistribution::from_counts(&[("A", 75), ("B", 25)]).unwrap();
        let u = variant_uniformity(&dist).unwrap();
        assert!(u < 1.0);
        assert!(u > 0.0);
    }

    #[test]
    fn uniformity_heavily_skewed_is_near_zero() {
        let dist = VariantDistribution::from_counts(&[
            ("Dom", 990),
            ("Minor1", 4),
            ("Minor2", 4),
            ("Minor3", 2),
        ])
        .unwrap();
        let u = variant_uniformity(&dist).unwrap();
        // Not perfectly zero because minor variants have non-zero proportion,
        // but should be well below 0.5.
        assert!(u < 0.5);
    }

    // ------------------------------------------------------------------
    // dominant_variant
    // ------------------------------------------------------------------

    #[test]
    fn dominant_variant_returns_highest_count() {
        let dist = VariantDistribution::from_counts(&[("Minor", 10), ("Major", 90)]).unwrap();
        let dom = dominant_variant(&dist).unwrap();
        assert_eq!(dom.name, "Major");
    }

    #[test]
    fn dominant_variant_empty_returns_none() {
        let dist = VariantDistribution {
            frequencies: vec![],
            total: 0,
        };
        assert!(dominant_variant(&dist).is_none());
    }

    #[test]
    fn dominant_variant_tie_returns_first() {
        // Both have count 50 — insertion-order first should win.
        let dist = VariantDistribution::from_counts(&[("First", 50), ("Second", 50)]).unwrap();
        let dom = dominant_variant(&dist).unwrap();
        assert_eq!(dom.name, "First");
    }

    #[test]
    fn dominant_variant_single_returns_that_variant() {
        let dist = VariantDistribution::from_counts(&[("Sole", 7)]).unwrap();
        let dom = dominant_variant(&dist).unwrap();
        assert_eq!(dom.name, "Sole");
    }

    // ------------------------------------------------------------------
    // dispatch_table
    // ------------------------------------------------------------------

    #[test]
    fn dispatch_table_entries_sorted_descending() {
        let dist =
            VariantDistribution::from_counts(&[("Rare", 5), ("Common", 80), ("Mid", 15)]).unwrap();
        let table = dispatch_table(&dist).unwrap();
        assert_eq!(table.entries[0].name, "Common");
        assert_eq!(table.entries[1].name, "Mid");
        assert_eq!(table.entries[2].name, "Rare");
    }

    #[test]
    fn dispatch_table_cumulative_reaches_one() {
        let dist = VariantDistribution::from_counts(&[("A", 10), ("B", 40), ("C", 50)]).unwrap();
        let table = dispatch_table(&dist).unwrap();
        let last = table.entries.last().unwrap();
        assert!(approx_eq(last.cumulative, 1.0, 1e-9));
    }

    #[test]
    fn dispatch_table_cumulative_is_monotone() {
        let dist = VariantDistribution::from_counts(&[("A", 10), ("B", 40), ("C", 50)]).unwrap();
        let table = dispatch_table(&dist).unwrap();
        for window in table.entries.windows(2) {
            assert!(window[1].cumulative >= window[0].cumulative);
        }
    }

    #[test]
    fn dispatch_table_contains_entropy_and_uniformity() {
        let dist = VariantDistribution::from_counts(&[("A", 1), ("B", 1)]).unwrap();
        let table = dispatch_table(&dist).unwrap();
        assert!(approx_eq(table.entropy, 1.0, 1e-9));
        assert!(approx_eq(table.uniformity, 1.0, 1e-9));
    }

    #[test]
    fn dispatch_table_empty_returns_error() {
        let dist = VariantDistribution {
            frequencies: vec![],
            total: 0,
        };
        assert_eq!(dispatch_table(&dist), Err(VariantError::EmptyDistribution));
    }

    // ------------------------------------------------------------------
    // chi_square_uniformity
    // ------------------------------------------------------------------

    #[test]
    fn chi_square_perfectly_uniform_is_zero() {
        let dist = VariantDistribution::from_counts(&[("A", 25), ("B", 25), ("C", 25), ("D", 25)])
            .unwrap();
        let chi2 = chi_square_uniformity(&dist).unwrap();
        assert!(chi2 < 1e-9);
    }

    #[test]
    fn chi_square_skewed_distribution_gives_large_value() {
        // Heavy skew → large chi-square
        let dist = VariantDistribution::from_counts(&[
            ("Dom", 970),
            ("Minor1", 10),
            ("Minor2", 10),
            ("Minor3", 10),
        ])
        .unwrap();
        let chi2 = chi_square_uniformity(&dist).unwrap();
        // Critical value for df=3, p=0.001 is ~16.3
        assert!(chi2 > 100.0);
    }

    #[test]
    fn chi_square_two_variants_symmetric() {
        // 75/25 → χ² = (75-50)²/50 + (25-50)²/50 = 25 + 25 = 25 (total=100, expected=50)
        let dist = VariantDistribution::from_counts(&[("A", 75), ("B", 25)]).unwrap();
        let chi2 = chi_square_uniformity(&dist).unwrap();
        assert!(approx_eq(chi2, 25.0, 1e-6));
    }

    #[test]
    fn chi_square_empty_returns_error() {
        let dist = VariantDistribution {
            frequencies: vec![],
            total: 0,
        };
        assert_eq!(
            chi_square_uniformity(&dist),
            Err(VariantError::EmptyDistribution)
        );
    }

    // ------------------------------------------------------------------
    // Known value cross-checks
    // ------------------------------------------------------------------

    #[test]
    fn known_entropy_three_variants() {
        // p = [0.5, 0.3, 0.2]
        // H = -(0.5*log2(0.5) + 0.3*log2(0.3) + 0.2*log2(0.2))
        //   = -(0.5*(-1) + 0.3*(-1.73697) + 0.2*(-2.32193))
        //   = -(−0.5 − 0.52109 − 0.46439)
        //   = 1.48548 bits
        let dist = VariantDistribution::from_counts(&[("A", 50), ("B", 30), ("C", 20)]).unwrap();
        let h = variant_entropy(&dist).unwrap();
        let expected =
            -(0.5_f64 * 0.5_f64.log2() + 0.3_f64 * 0.3_f64.log2() + 0.2_f64 * 0.2_f64.log2());
        assert!(approx_eq(h, expected, 1e-12));
    }

    #[test]
    fn error_display_messages_are_correct() {
        assert_eq!(
            VariantError::EmptyDistribution.to_string(),
            "variant distribution is empty"
        );
        assert_eq!(VariantError::ZeroTotal.to_string(), "total count is zero");
        assert_eq!(
            VariantError::InconsistentTotal.to_string(),
            "sum of frequencies does not match total"
        );
    }

    // ------------------------------------------------------------------
    // Serde round-trip
    // ------------------------------------------------------------------

    #[test]
    fn variant_distribution_serde_round_trip() {
        let dist = VariantDistribution::from_counts(&[("X", 3), ("Y", 7)]).unwrap();
        let json = serde_json::to_string(&dist).unwrap();
        let restored: VariantDistribution = serde_json::from_str(&json).unwrap();
        assert_eq!(dist, restored);
    }

    #[test]
    fn dispatch_table_serde_round_trip() {
        let dist = VariantDistribution::from_counts(&[("A", 10), ("B", 90)]).unwrap();
        let table = dispatch_table(&dist).unwrap();
        let json = serde_json::to_string(&table).unwrap();
        let restored: DispatchTable = serde_json::from_str(&json).unwrap();
        assert_eq!(table, restored);
    }
}
