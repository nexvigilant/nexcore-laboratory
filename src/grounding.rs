//! # GroundsTo implementations for nexcore-laboratory types
//!
//! Connects virtual word/concept experiment types to the Lex Primitiva type system.
//!
//! ## Mapping (mu) Focus
//!
//! The laboratory decomposes words into primitive compositions (mu: mapping),
//! runs spectral analysis (sigma: sequence), and reacts concepts (Sigma: sum).

use nexcore_lex_primitiva::grounding::GroundsTo;
use nexcore_lex_primitiva::primitiva::{LexPrimitiva, PrimitiveComposition};

use crate::{
    BatchResult, ClassDistribution, ExperimentResult, ReactionResult, Specimen, SpectralLine,
};

// ---------------------------------------------------------------------------
// Input types -- mu (Mapping) dominant
// ---------------------------------------------------------------------------

/// Specimen: T2-P (mu + sigma), dominant mu
///
/// A word/concept to be experimentally decomposed.
/// Mapping-dominant: the specimen maps a word to its primitive composition.
/// Sequence is secondary (the word is a character sequence).
impl GroundsTo for Specimen {
    fn primitive_composition() -> PrimitiveComposition {
        PrimitiveComposition::new(vec![
            LexPrimitiva::Mapping,  // mu -- word -> primitive mapping
            LexPrimitiva::Sequence, // sigma -- character sequence
        ])
        .with_dominant(LexPrimitiva::Mapping, 0.85)
    }
}

// ---------------------------------------------------------------------------
// Result types -- multi-primitive composites
// ---------------------------------------------------------------------------

/// SpectralLine: T2-P (N + Sigma), dominant N
///
/// A single line in the spectral output: primitive name + weight.
/// Quantity-dominant: the weight IS a numeric measurement.
impl GroundsTo for SpectralLine {
    fn primitive_composition() -> PrimitiveComposition {
        PrimitiveComposition::new(vec![
            LexPrimitiva::Quantity, // N -- spectral weight
            LexPrimitiva::Sum,      // Sigma -- primitive classification
        ])
        .with_dominant(LexPrimitiva::Quantity, 0.85)
    }
}

/// ExperimentResult: T2-C (mu + N + sigma + kappa), dominant mu
///
/// Full result of decomposing a specimen: composition, molecular weight, tier.
/// Mapping-dominant: the result maps a word to a complete analysis.
/// Quantity is secondary (molecular weight, spectral weights).
/// Sequence is tertiary (spectral line ordering).
/// Comparison is quaternary (tier classification).
impl GroundsTo for ExperimentResult {
    fn primitive_composition() -> PrimitiveComposition {
        PrimitiveComposition::new(vec![
            LexPrimitiva::Mapping,    // mu -- word -> analysis mapping
            LexPrimitiva::Quantity,   // N -- molecular weight
            LexPrimitiva::Sequence,   // sigma -- spectral ordering
            LexPrimitiva::Comparison, // kappa -- tier classification
        ])
        .with_dominant(LexPrimitiva::Mapping, 0.80)
    }
}

/// ReactionResult: T2-C (mu + Sigma + N + causality), dominant mu
///
/// Result of reacting two concepts: product composition, shared catalysts.
/// Mapping-dominant: the reaction maps two inputs to a combined product.
/// Sum is secondary (combining compositions).
/// Quantity is tertiary (combined weight).
/// Causality is quaternary (reaction direction).
impl GroundsTo for ReactionResult {
    fn primitive_composition() -> PrimitiveComposition {
        PrimitiveComposition::new(vec![
            LexPrimitiva::Mapping,   // mu -- two inputs -> product
            LexPrimitiva::Sum,       // Sigma -- composition combination
            LexPrimitiva::Quantity,  // N -- combined weight
            LexPrimitiva::Causality, // causality -- reaction direction
        ])
        .with_dominant(LexPrimitiva::Mapping, 0.80)
    }
}

/// BatchResult: T2-C (sigma + mu + N + kappa), dominant sigma
///
/// Results from running multiple reactions in batch.
/// Sequence-dominant: the batch IS an ordered collection of results.
impl GroundsTo for BatchResult {
    fn primitive_composition() -> PrimitiveComposition {
        PrimitiveComposition::new(vec![
            LexPrimitiva::Sequence,   // sigma -- ordered batch
            LexPrimitiva::Mapping,    // mu -- individual reactions
            LexPrimitiva::Quantity,   // N -- counts, statistics
            LexPrimitiva::Comparison, // kappa -- classification
        ])
        .with_dominant(LexPrimitiva::Sequence, 0.80)
    }
}

/// ClassDistribution: T2-P (N + Sigma + kappa), dominant N
///
/// Distribution of tier classifications across batch results.
/// Quantity-dominant: counts per tier class.
impl GroundsTo for ClassDistribution {
    fn primitive_composition() -> PrimitiveComposition {
        PrimitiveComposition::new(vec![
            LexPrimitiva::Quantity,   // N -- counts per class
            LexPrimitiva::Sum,        // Sigma -- tier classes
            LexPrimitiva::Comparison, // kappa -- classification
        ])
        .with_dominant(LexPrimitiva::Quantity, 0.85)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use nexcore_lex_primitiva::tier::Tier;

    #[test]
    fn specimen_is_t2p() {
        assert_eq!(Specimen::tier(), Tier::T2Primitive);
        assert_eq!(Specimen::dominant_primitive(), Some(LexPrimitiva::Mapping));
    }

    #[test]
    fn experiment_result_is_t2c() {
        assert_eq!(ExperimentResult::tier(), Tier::T2Composite);
        assert_eq!(
            ExperimentResult::dominant_primitive(),
            Some(LexPrimitiva::Mapping)
        );
    }

    #[test]
    fn reaction_result_is_t2c() {
        assert_eq!(ReactionResult::tier(), Tier::T2Composite);
    }

    #[test]
    fn batch_result_is_t2c() {
        assert_eq!(BatchResult::tier(), Tier::T2Composite);
        assert_eq!(
            BatchResult::dominant_primitive(),
            Some(LexPrimitiva::Sequence)
        );
    }

    #[test]
    fn all_confidences_valid() {
        let compositions = [
            Specimen::primitive_composition(),
            SpectralLine::primitive_composition(),
            ExperimentResult::primitive_composition(),
            ReactionResult::primitive_composition(),
            BatchResult::primitive_composition(),
            ClassDistribution::primitive_composition(),
        ];
        for comp in &compositions {
            assert!(comp.confidence >= 0.80);
            assert!(comp.confidence <= 1.0);
        }
    }
}
