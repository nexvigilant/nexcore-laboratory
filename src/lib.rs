//! # NexVigilant Core — Laboratory — Virtual Word/Concept Experiment Engine
//!
//! A laboratory for running experiments on words and concepts through the
//! lens of Lex Primitiva molecular weight and T1 primitive decomposition.
//!
//! ## Experiment Pipeline
//!
//! ```text
//! Specimen → Decompose → Weigh → Classify → Analyze → Report
//!    (word)    (μ)        (Σ+N)    (κ+∂)      (σ)      (→)
//! ```
//!
//! ## Chemical Reaction Metaphor
//!
//! Two concepts can be "reacted" — combining their primitive compositions
//! to produce a new compound. Shared primitives act as catalysts.
//!
//! ## Tier: T2-C (μ + Σ + κ + × + σ)

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![cfg_attr(
    not(test),
    deny(clippy::unwrap_used, clippy::expect_used, clippy::panic)
)]

pub mod experiment;
pub mod grounding;
pub mod reaction;
pub mod variant_analysis;

// Re-exports for ergonomic API
pub use experiment::{ExperimentResult, Specimen, SpectralLine, run_experiment};
pub use reaction::{BatchResult, ClassDistribution, ReactionResult, react, run_batch};
