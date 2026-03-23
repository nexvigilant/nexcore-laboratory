//! Chemical reaction metaphor for concept combination.
//!
//! "Reacting" two concepts means combining their primitive compositions
//! and analyzing the resulting product through the molecular weight lens.
//!
//! | Chemistry | Laboratory |
//! |-----------|------------|
//! | Reactant A | Concept A primitives |
//! | Reactant B | Concept B primitives |
//! | Catalyst | Shared primitives (lower activation energy) |
//! | Product | Union of all primitives |
//! | Enthalpy | Weight delta (exothermic if product lighter) |
//!
//! ## Tier: T2-C (× Product + κ Comparison + Σ Sum)

use crate::experiment::{ExperimentResult, Specimen, run_experiment};
use nexcore_lex_primitiva::primitiva::LexPrimitiva;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

// ============================================================================
// Reaction Types
// ============================================================================

/// Result of "reacting" two concepts.
///
/// Tier: T2-C (× + κ + Σ + μ)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionResult {
    /// Reactant A experiment result
    pub reactant_a: ExperimentResult,
    /// Reactant B experiment result
    pub reactant_b: ExperimentResult,
    /// Product experiment result (combined)
    pub product: ExperimentResult,
    /// Shared primitives (catalyst)
    pub catalyst: Vec<String>,
    /// Primitives unique to A
    pub unique_a: Vec<String>,
    /// Primitives unique to B
    pub unique_b: Vec<String>,
    /// Jaccard similarity of primitive sets
    pub jaccard_similarity: f64,
    /// Weight delta: product - (A + B). Negative = exothermic (more compact)
    pub enthalpy: f64,
    /// Whether the reaction is exothermic (product lighter than sum of parts)
    pub exothermic: bool,
    /// Human-readable reaction equation
    pub equation: String,
    /// Interpretation of the reaction
    pub interpretation: String,
}

// ============================================================================
// react()
// ============================================================================

/// React two specimens — combine their primitive compositions.
///
/// The product is the *unique* union of all primitives from both reactants.
/// Shared primitives act as catalysts (already present in both, reduce novelty).
#[must_use]
pub fn react(a: &Specimen, b: &Specimen) -> ReactionResult {
    let set_a: HashSet<LexPrimitiva> = a.primitives.iter().copied().collect();
    let set_b: HashSet<LexPrimitiva> = b.primitives.iter().copied().collect();

    let catalyst: Vec<LexPrimitiva> = set_a.intersection(&set_b).copied().collect();
    let unique_a: Vec<LexPrimitiva> = set_a.difference(&set_b).copied().collect();
    let unique_b: Vec<LexPrimitiva> = set_b.difference(&set_a).copied().collect();
    let product_prims: Vec<LexPrimitiva> = set_a.union(&set_b).copied().collect();

    let jaccard = if set_a.is_empty() && set_b.is_empty() {
        0.0
    } else {
        catalyst.len() as f64 / product_prims.len() as f64
    };

    let product_name = format!("{}+{}", a.name, b.name);
    let product_specimen = Specimen::new(&product_name, product_prims);

    let result_a = run_experiment(a);
    let result_b = run_experiment(b);
    let result_product = run_experiment(&product_specimen);

    let sum_weight = result_a.molecular_weight + result_b.molecular_weight;
    let enthalpy = round3(result_product.molecular_weight - sum_weight);
    let exothermic = enthalpy < 0.0;

    let equation = format!(
        "{} [{}] + {} [{}] → {} [{}]",
        a.name, result_a.formula, b.name, result_b.formula, product_name, result_product.formula,
    );

    let interpretation = if exothermic {
        format!(
            "Exothermic reaction (ΔH={:.2} Da): {} shared primitives reduce redundancy. Product is more compact than sum of parts.",
            enthalpy,
            catalyst.len()
        )
    } else if enthalpy.abs() < 0.001 {
        format!(
            "Neutral reaction: no shared primitives. Product is exactly the sum of parts ({:.2} Da).",
            result_product.molecular_weight
        )
    } else {
        format!(
            "Endothermic reaction (ΔH=+{:.2} Da): product is heavier than expected.",
            enthalpy
        )
    };

    ReactionResult {
        reactant_a: result_a,
        reactant_b: result_b,
        product: result_product,
        catalyst: catalyst.iter().map(|p| p.symbol().to_string()).collect(),
        unique_a: unique_a.iter().map(|p| p.symbol().to_string()).collect(),
        unique_b: unique_b.iter().map(|p| p.symbol().to_string()).collect(),
        jaccard_similarity: round3(jaccard),
        enthalpy,
        exothermic,
        equation,
        interpretation,
    }
}

fn round3(v: f64) -> f64 {
    (v * 1000.0).round() / 1000.0
}

// ============================================================================
// Batch Experiment
// ============================================================================

/// Result of running experiments on multiple specimens.
///
/// Tier: T2-C (σ Sequence + κ Comparison + Σ Sum)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchResult {
    /// Individual experiment results
    pub experiments: Vec<ExperimentResult>,
    /// Lightest concept
    pub lightest: String,
    /// Heaviest concept
    pub heaviest: String,
    /// Average molecular weight across all specimens
    pub average_weight: f64,
    /// Weight standard deviation
    pub weight_std_dev: f64,
    /// Count by transfer class
    pub class_distribution: ClassDistribution,
}

/// Distribution of specimens across transfer classes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassDistribution {
    pub light: usize,
    pub medium: usize,
    pub heavy: usize,
}

/// Run experiments on a batch of specimens.
#[must_use]
pub fn run_batch(specimens: &[Specimen]) -> BatchResult {
    let experiments: Vec<ExperimentResult> = specimens.iter().map(run_experiment).collect();

    if experiments.is_empty() {
        return BatchResult {
            experiments: vec![],
            lightest: String::new(),
            heaviest: String::new(),
            average_weight: 0.0,
            weight_std_dev: 0.0,
            class_distribution: ClassDistribution {
                light: 0,
                medium: 0,
                heavy: 0,
            },
        };
    }

    let weights: Vec<f64> = experiments.iter().map(|e| e.molecular_weight).collect();
    let avg = weights.iter().sum::<f64>() / weights.len() as f64;
    let variance = weights.iter().map(|w| (w - avg).powi(2)).sum::<f64>() / weights.len() as f64;
    let std_dev = variance.sqrt();

    let lightest = experiments
        .iter()
        .min_by(|a, b| {
            a.molecular_weight
                .partial_cmp(&b.molecular_weight)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|e| e.name.clone())
        .unwrap_or_default();

    let heaviest = experiments
        .iter()
        .max_by(|a, b| {
            a.molecular_weight
                .partial_cmp(&b.molecular_weight)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|e| e.name.clone())
        .unwrap_or_default();

    let mut light = 0;
    let mut medium = 0;
    let mut heavy = 0;
    for e in &experiments {
        if e.transfer_class.contains("Light") {
            light += 1;
        } else if e.transfer_class.contains("Medium") {
            medium += 1;
        } else {
            heavy += 1;
        }
    }

    BatchResult {
        experiments,
        lightest,
        heaviest,
        average_weight: round3(avg),
        weight_std_dev: round3(std_dev),
        class_distribution: ClassDistribution {
            light,
            medium,
            heavy,
        },
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn specimen_guardian() -> Specimen {
        Specimen::new(
            "Guardian",
            vec![
                LexPrimitiva::State,
                LexPrimitiva::Boundary,
                LexPrimitiva::Comparison,
            ],
        )
    }

    fn specimen_signal() -> Specimen {
        Specimen::new(
            "Signal",
            vec![LexPrimitiva::Boundary, LexPrimitiva::Quantity],
        )
    }

    #[test]
    fn test_react_produces_product() {
        let result = react(&specimen_guardian(), &specimen_signal());
        assert!(result.product.name.contains("Guardian"));
        assert!(result.product.name.contains("Signal"));
    }

    #[test]
    fn test_react_shared_primitives_are_catalyst() {
        let result = react(&specimen_guardian(), &specimen_signal());
        // Both share ∂ (Boundary)
        assert!(result.catalyst.contains(&"∂".to_string()));
    }

    #[test]
    fn test_react_exothermic_with_overlap() {
        let result = react(&specimen_guardian(), &specimen_signal());
        // Product is union (4 unique) but sum of parts is 3+2=5 primitives
        // Product MW < sum of parts MW → exothermic
        assert!(result.exothermic);
        assert!(result.enthalpy < 0.0);
    }

    #[test]
    fn test_react_no_overlap_neutral() {
        let a = Specimen::new("A", vec![LexPrimitiva::State]);
        let b = Specimen::new("B", vec![LexPrimitiva::Quantity]);
        let result = react(&a, &b);
        // No shared primitives → product = sum of parts
        assert!(result.catalyst.is_empty());
        // Product weight should equal sum (within rounding)
        let sum = result.reactant_a.molecular_weight + result.reactant_b.molecular_weight;
        let delta = (result.product.molecular_weight - sum).abs();
        assert!(delta < 0.01, "Expected neutral, got delta={:.3}", delta);
    }

    #[test]
    fn test_react_jaccard_similarity() {
        let result = react(&specimen_guardian(), &specimen_signal());
        // Guardian: {ς, ∂, κ}, Signal: {∂, N} → shared: {∂}, union: {ς, ∂, κ, N}
        // Jaccard = 1/4 = 0.25
        assert!(
            (result.jaccard_similarity - 0.25).abs() < 0.01,
            "Expected Jaccard ~0.25, got {:.3}",
            result.jaccard_similarity
        );
    }

    #[test]
    fn test_react_equation_format() {
        let result = react(&specimen_guardian(), &specimen_signal());
        assert!(result.equation.contains("→"));
        assert!(result.equation.contains("+"));
    }

    #[test]
    fn test_batch_experiment() {
        let specimens = vec![
            specimen_guardian(),
            specimen_signal(),
            Specimen::new(
                "Cascade",
                vec![LexPrimitiva::Sequence, LexPrimitiva::Causality],
            ),
        ];
        let batch = run_batch(&specimens);
        assert_eq!(batch.experiments.len(), 3);
        assert!(batch.average_weight > 0.0);
        assert!(!batch.lightest.is_empty());
        assert!(!batch.heaviest.is_empty());
    }

    #[test]
    fn test_batch_empty() {
        let batch = run_batch(&[]);
        assert!(batch.experiments.is_empty());
        assert_eq!(batch.average_weight, 0.0);
    }

    #[test]
    fn test_batch_class_distribution() {
        let specimens = vec![
            Specimen::new("Light1", vec![LexPrimitiva::Quantity]),
            Specimen::new("Light2", vec![LexPrimitiva::Boundary]),
            Specimen::new(
                "Heavy1",
                vec![
                    LexPrimitiva::State,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Causality,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Persistence,
                ],
            ),
        ];
        let batch = run_batch(&specimens);
        assert_eq!(batch.class_distribution.light, 2);
        assert!(batch.class_distribution.heavy >= 1 || batch.class_distribution.medium >= 1);
    }

    #[test]
    fn test_batch_std_dev_positive_for_varied() {
        let specimens = vec![
            Specimen::new("Tiny", vec![LexPrimitiva::Quantity]),
            Specimen::new(
                "Huge",
                vec![
                    LexPrimitiva::State,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Causality,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Persistence,
                    LexPrimitiva::Product,
                ],
            ),
        ];
        let batch = run_batch(&specimens);
        assert!(batch.weight_std_dev > 0.0);
    }
}
