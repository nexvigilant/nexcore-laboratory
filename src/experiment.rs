//! Experiment types and pipeline stages.
//!
//! A laboratory experiment takes a word/concept through 5 stages:
//! 1. **Decompose** — resolve to T1 primitives (μ Mapping)
//! 2. **Weigh** — compute molecular weight in daltons (Σ Sum + N Quantity)
//! 3. **Classify** — determine tier and transfer confidence (κ Comparison)
//! 4. **Analyze** — spectral analysis of constituent masses (∂ Boundary)
//! 5. **Report** — structured result with interpretation (σ Sequence)
//!
//! ## Tier: T2-C (μ + Σ + κ + ∂ + σ)

use nexcore_lex_primitiva::molecular_weight::{
    AtomicMass, MolecularFormula, MolecularWeight, TierPrediction, TransferClass,
};
use nexcore_lex_primitiva::primitiva::LexPrimitiva;
use serde::{Deserialize, Serialize};

// ============================================================================
// Experiment Input
// ============================================================================

/// A word or concept to experiment on.
///
/// Tier: T2-P (μ Mapping + ∃ Existence)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Specimen {
    /// Human-readable name of the concept
    pub name: String,
    /// T1 primitives that compose this concept
    pub primitives: Vec<LexPrimitiva>,
}

impl Specimen {
    /// Create a new specimen from name and primitive list.
    #[must_use]
    pub fn new(name: &str, primitives: Vec<LexPrimitiva>) -> Self {
        Self {
            name: name.to_string(),
            primitives,
        }
    }

    /// Compute molecular formula for this specimen.
    #[must_use]
    pub fn formula(&self) -> MolecularFormula {
        MolecularFormula::new(&self.name).with_all(&self.primitives)
    }

    /// Compute molecular weight.
    #[must_use]
    pub fn weight(&self) -> MolecularWeight {
        self.formula().weight()
    }
}

// ============================================================================
// Experiment Result
// ============================================================================

/// Complete experiment result for a single specimen.
///
/// Tier: T2-C (μ + Σ + κ + N + σ)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExperimentResult {
    /// Specimen name
    pub name: String,
    /// Chemical formula string (e.g., "ς∂κ")
    pub formula: String,
    /// Molecular weight in daltons
    pub molecular_weight: f64,
    /// Number of constituent primitives
    pub primitive_count: usize,
    /// Average mass per primitive
    pub average_mass: f64,
    /// Weight-based transfer class
    pub transfer_class: String,
    /// Predicted transfer confidence (0.0-1.0)
    pub transfer_confidence: f64,
    /// Hybrid tier-aware prediction
    pub tier_prediction: String,
    /// Hybrid transfer confidence
    pub hybrid_transfer: f64,
    /// Individual constituent masses
    pub spectrum: Vec<SpectralLine>,
    /// Human-readable interpretation
    pub interpretation: String,
}

/// A single line in the mass spectrum — one constituent primitive.
///
/// Tier: T2-P (N Quantity + μ Mapping)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralLine {
    /// Primitive name
    pub primitive: String,
    /// Unicode symbol
    pub symbol: String,
    /// Mass in bits
    pub mass_bits: f64,
    /// Raw frequency in codebase
    pub frequency: u32,
    /// Probability (freq / total)
    pub probability: f64,
}

// ============================================================================
// Pipeline: run_experiment
// ============================================================================

/// Run a complete experiment on a specimen.
///
/// Pipeline: Decompose → Weigh → Classify → Analyze → Report
#[must_use]
pub fn run_experiment(specimen: &Specimen) -> ExperimentResult {
    let formula = specimen.formula();
    let weight = formula.weight();

    // Stage 3: Spectral analysis — individual masses
    let spectrum: Vec<SpectralLine> = formula
        .atomic_masses()
        .iter()
        .map(|am| SpectralLine {
            primitive: am.primitive().name().to_string(),
            symbol: am.primitive().symbol().to_string(),
            mass_bits: round3(am.bits()),
            frequency: am.frequency(),
            probability: round3(am.probability()),
        })
        .collect();

    // Stage 4: Classification
    let tc = weight.transfer_class();
    let tp = weight.tier_aware_class();
    let transfer = weight.predicted_transfer();
    let hybrid = weight.predicted_transfer_hybrid();

    // Stage 5: Interpretation
    let interpretation = format!(
        "{} [{:.2} Da, {} primitives] — {} ({:.0}% transfer, hybrid {:.0}%)",
        specimen.name,
        weight.daltons(),
        weight.primitive_count(),
        tc,
        transfer * 100.0,
        hybrid * 100.0,
    );

    ExperimentResult {
        name: specimen.name.clone(),
        formula: formula.formula_string(),
        molecular_weight: round3(weight.daltons()),
        primitive_count: weight.primitive_count(),
        average_mass: round3(weight.average_mass()),
        transfer_class: format!("{tc}"),
        transfer_confidence: round3(transfer),
        tier_prediction: format!("{tp}"),
        hybrid_transfer: round3(hybrid),
        spectrum,
        interpretation,
    }
}

fn round3(v: f64) -> f64 {
    (v * 1000.0).round() / 1000.0
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_specimen_creation() {
        let s = Specimen::new(
            "Vigilance",
            vec![
                LexPrimitiva::State,
                LexPrimitiva::Boundary,
                LexPrimitiva::Frequency,
                LexPrimitiva::Causality,
            ],
        );
        assert_eq!(s.name, "Vigilance");
        assert_eq!(s.primitives.len(), 4);
    }

    #[test]
    fn test_specimen_weight_is_positive() {
        let s = Specimen::new(
            "Signal",
            vec![LexPrimitiva::Boundary, LexPrimitiva::Quantity],
        );
        assert!(s.weight().daltons() > 0.0);
    }

    #[test]
    fn test_run_experiment_produces_result() {
        let s = Specimen::new(
            "Guardian",
            vec![
                LexPrimitiva::State,
                LexPrimitiva::Boundary,
                LexPrimitiva::Comparison,
            ],
        );
        let result = run_experiment(&s);
        assert_eq!(result.name, "Guardian");
        assert_eq!(result.primitive_count, 3);
        assert!(result.molecular_weight > 0.0);
        assert!(!result.formula.is_empty());
        assert_eq!(result.spectrum.len(), 3);
    }

    #[test]
    fn test_experiment_spectrum_has_correct_primitives() {
        let s = Specimen::new(
            "Cascade",
            vec![LexPrimitiva::Sequence, LexPrimitiva::Causality],
        );
        let result = run_experiment(&s);
        let symbols: Vec<&str> = result.spectrum.iter().map(|l| l.symbol.as_str()).collect();
        assert!(symbols.contains(&"σ"));
        assert!(symbols.contains(&"→"));
    }

    #[test]
    fn test_experiment_transfer_in_bounds() {
        let s = Specimen::new("Test", vec![LexPrimitiva::Quantity]);
        let result = run_experiment(&s);
        assert!(result.transfer_confidence >= 0.05);
        assert!(result.transfer_confidence <= 0.98);
    }

    #[test]
    fn test_experiment_interpretation_contains_name() {
        let s = Specimen::new(
            "Polypharmacy",
            vec![
                LexPrimitiva::Product,
                LexPrimitiva::Boundary,
                LexPrimitiva::Comparison,
                LexPrimitiva::Quantity,
            ],
        );
        let result = run_experiment(&s);
        assert!(result.interpretation.contains("Polypharmacy"));
    }

    #[test]
    fn test_heavy_concept_classified_correctly() {
        let s = Specimen::new(
            "ComplexDomain",
            vec![
                LexPrimitiva::State,
                LexPrimitiva::Boundary,
                LexPrimitiva::Causality,
                LexPrimitiva::Sequence,
                LexPrimitiva::Mapping,
                LexPrimitiva::Persistence,
                LexPrimitiva::Product,
            ],
        );
        let result = run_experiment(&s);
        // 7 primitives → heavy, domain-locked
        assert!(result.molecular_weight > 18.0);
        assert!(result.transfer_class.contains("Heavy"));
    }

    #[test]
    fn test_light_concept_high_transfer() {
        let s = Specimen::new("Count", vec![LexPrimitiva::Quantity]);
        let result = run_experiment(&s);
        assert!(result.transfer_confidence > 0.85);
    }

    /// Demo experiment — prints results for 5 domain words.
    /// Run with: cargo test -p nexcore-laboratory -- --nocapture demo_word_spectroscopy
    #[test]
    fn demo_word_spectroscopy() {
        let specimens = vec![
            Specimen::new(
                "Vigilance",
                vec![
                    LexPrimitiva::State,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Causality,
                ],
            ),
            Specimen::new(
                "Signal",
                vec![LexPrimitiva::Boundary, LexPrimitiva::Quantity],
            ),
            Specimen::new(
                "Guardian",
                vec![
                    LexPrimitiva::State,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Comparison,
                ],
            ),
            Specimen::new(
                "Cascade",
                vec![
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Causality,
                    LexPrimitiva::Irreversibility,
                ],
            ),
            Specimen::new(
                "Polypharmacy",
                vec![
                    LexPrimitiva::Product,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Quantity,
                ],
            ),
        ];

        eprintln!("\n--- WORD SPECTROSCOPY ---");
        for s in &specimens {
            let r = run_experiment(s);
            eprintln!(
                "  {} [{}] = {:.2} Da | {} | tier: {} | transfer: {:.0}%",
                r.name,
                r.formula,
                r.molecular_weight,
                r.transfer_class,
                r.tier_prediction,
                r.hybrid_transfer * 100.0,
            );
            for line in &r.spectrum {
                eprintln!(
                    "    {} {} = {:.2} bits (freq={}, p={:.3})",
                    line.symbol, line.primitive, line.mass_bits, line.frequency, line.probability,
                );
            }
        }

        let batch = crate::run_batch(&specimens);
        eprintln!("\n--- BATCH SUMMARY ---");
        eprintln!("  Lightest: {}", batch.lightest);
        eprintln!("  Heaviest: {}", batch.heaviest);
        eprintln!(
            "  Avg MW: {:.2} Da (σ={:.2})",
            batch.average_weight, batch.weight_std_dev
        );
        eprintln!(
            "  Distribution: {} light, {} medium, {} heavy",
            batch.class_distribution.light,
            batch.class_distribution.medium,
            batch.class_distribution.heavy,
        );

        // React Guardian + Signal
        let guardian = Specimen::new(
            "Guardian",
            vec![
                LexPrimitiva::State,
                LexPrimitiva::Boundary,
                LexPrimitiva::Comparison,
            ],
        );
        let signal = Specimen::new(
            "Signal",
            vec![LexPrimitiva::Boundary, LexPrimitiva::Quantity],
        );
        let rxn = crate::react(&guardian, &signal);
        eprintln!("\n--- REACTION ---");
        eprintln!("  {}", rxn.equation);
        eprintln!("  Catalyst: {:?}", rxn.catalyst);
        eprintln!(
            "  ΔH = {:.2} Da ({})",
            rxn.enthalpy,
            if rxn.exothermic {
                "exothermic"
            } else {
                "endothermic"
            }
        );
        eprintln!("  Jaccard: {:.2}", rxn.jaccard_similarity);
        eprintln!("  {}", rxn.interpretation);

        // Assertions to keep it a real test
        assert_eq!(batch.experiments.len(), 5);
        assert!(batch.average_weight > 0.0);
        assert!(rxn.exothermic); // Shared ∂ should make it exothermic
    }

    /// PSAS decomposition — Preemptive State Analytical Services
    /// Run with: cargo test -p nexcore-laboratory -- --nocapture demo_psas
    #[test]
    fn demo_psas() {
        let preemptive = Specimen::new(
            "Preemptive",
            vec![
                LexPrimitiva::Causality,
                LexPrimitiva::Frequency,
                LexPrimitiva::Irreversibility,
                LexPrimitiva::Boundary,
            ],
        );
        let state_word = Specimen::new("State", vec![LexPrimitiva::State]);
        let analytical = Specimen::new(
            "Analytical",
            vec![
                LexPrimitiva::Comparison,
                LexPrimitiva::Mapping,
                LexPrimitiva::Quantity,
            ],
        );
        let services = Specimen::new(
            "Services",
            vec![LexPrimitiva::Sequence, LexPrimitiva::Product],
        );
        let psas = Specimen::new(
            "PSAS",
            vec![
                LexPrimitiva::Causality,
                LexPrimitiva::Frequency,
                LexPrimitiva::Irreversibility,
                LexPrimitiva::Boundary,
                LexPrimitiva::State,
                LexPrimitiva::Comparison,
                LexPrimitiva::Mapping,
                LexPrimitiva::Quantity,
                LexPrimitiva::Sequence,
                LexPrimitiva::Product,
            ],
        );

        // Individual spectroscopy
        let components = [&preemptive, &state_word, &analytical, &services, &psas];
        eprintln!("\n{}", "=".repeat(70));
        eprintln!(" PSAS DECOMPOSITION: Preemptive State Analytical Services");
        eprintln!("{}\n", "=".repeat(70));

        for s in &components {
            let r = run_experiment(s);
            eprintln!(
                "  {:<14} [{:<12}] = {:>6.2} Da | {:<10} | {} | transfer: {:>3.0}%",
                r.name,
                r.formula,
                r.molecular_weight,
                r.transfer_class,
                r.tier_prediction,
                r.hybrid_transfer * 100.0,
            );
        }

        // Sum-of-parts analysis
        let mw_pre = run_experiment(&preemptive).molecular_weight;
        let mw_st = run_experiment(&state_word).molecular_weight;
        let mw_an = run_experiment(&analytical).molecular_weight;
        let mw_sv = run_experiment(&services).molecular_weight;
        let mw_psas = run_experiment(&psas).molecular_weight;
        let sum_parts = mw_pre + mw_st + mw_an + mw_sv;

        eprintln!("\n  --- SYNTHESIS ---");
        eprintln!(
            "  Sum of parts:   {:.2} Da ({:.2} + {:.2} + {:.2} + {:.2})",
            sum_parts, mw_pre, mw_st, mw_an, mw_sv
        );
        eprintln!(
            "  Actual PSAS:    {:.2} Da (union of 10 primitives)",
            mw_psas
        );
        eprintln!("  DeltaH synth:   {:+.2} Da", mw_psas - sum_parts);
        eprintln!(
            "  Compression:    {:.1}%  (no shared primitives = 0% overlap)",
            (1.0 - mw_psas / sum_parts) * 100.0
        );

        // Component reactions
        let pairs: Vec<(&Specimen, &Specimen, &str)> = vec![
            (&preemptive, &state_word, "Preemptive + State"),
            (&preemptive, &analytical, "Preemptive + Analytical"),
            (&preemptive, &services, "Preemptive + Services"),
            (&analytical, &services, "Analytical + Services"),
            (&state_word, &analytical, "State + Analytical"),
        ];

        eprintln!("\n  --- COMPONENT REACTIONS ---");
        for (a, b, label) in &pairs {
            let rxn = crate::react(a, b);
            eprintln!(
                "  {:<26} Catalyst: {:<8} Jaccard: {:.2}  DeltaH: {:+.2} Da ({})",
                label,
                format!("{:?}", rxn.catalyst),
                rxn.jaccard_similarity,
                rxn.enthalpy,
                if rxn.exothermic {
                    "exo"
                } else if rxn.enthalpy.abs() < 0.001 {
                    "neutral"
                } else {
                    "endo"
                },
            );
        }

        // React PSAS with NexVigilant
        let nexvigilant = Specimen::new(
            "NexVigilant",
            vec![
                LexPrimitiva::State,
                LexPrimitiva::Boundary,
                LexPrimitiva::Causality,
                LexPrimitiva::Frequency,
                LexPrimitiva::Comparison,
                LexPrimitiva::Persistence,
                LexPrimitiva::Sequence,
            ],
        );
        let rxn_nv = crate::react(&psas, &nexvigilant);

        eprintln!("\n  --- PSAS vs NEXVIGILANT ---");
        eprintln!("  {}", rxn_nv.equation);
        eprintln!("  Catalyst (shared):   {:?}", rxn_nv.catalyst);
        eprintln!("  Only in PSAS:        {:?}", rxn_nv.unique_a);
        eprintln!("  Only in NexVigilant: {:?}", rxn_nv.unique_b);
        eprintln!(
            "  Jaccard: {:.2}  DeltaH: {:+.2} Da ({})",
            rxn_nv.jaccard_similarity,
            rxn_nv.enthalpy,
            if rxn_nv.exothermic {
                "highly exothermic"
            } else {
                "endothermic"
            },
        );
        eprintln!(
            "  Overlap: {}/{} primitives shared",
            rxn_nv.catalyst.len(),
            rxn_nv.product.primitive_count,
        );

        // Primitive coverage
        eprintln!("\n  --- PRIMITIVE COVERAGE (PSAS = 10/16 T1) ---");
        let covered = [
            ("→", "Causality", "Preemptive"),
            ("ν", "Frequency", "Preemptive"),
            ("∝", "Irreversibility", "Preemptive"),
            ("∂", "Boundary", "Preemptive"),
            ("ς", "State", "State"),
            ("κ", "Comparison", "Analytical"),
            ("μ", "Mapping", "Analytical"),
            ("N", "Quantity", "Analytical"),
            ("σ", "Sequence", "Services"),
            ("×", "Product", "Services"),
        ];
        for (sym, name, source) in &covered {
            eprintln!("  [x] {} {:<16} <- {}", sym, name, source);
        }
        let missing = [
            ("ρ", "Recursion"),
            ("∅", "Void"),
            ("∃", "Existence"),
            ("π", "Persistence"),
            ("λ", "Location"),
            ("Σ", "Sum"),
        ];
        eprintln!();
        for (sym, name) in &missing {
            eprintln!("  [ ] {} {:<16}    (not in PSAS)", sym, name);
        }

        // Assertions
        assert!(mw_psas > 35.0);
        assert_eq!(psas.primitives.len(), 10);
        assert!(!rxn_nv.catalyst.is_empty());
    }

    /// Wordsmith engine — explore naming candidates for maximum primitive coverage.
    /// Run with: cargo test -p nexcore-laboratory -- --nocapture demo_wordsmith
    #[test]
    fn demo_wordsmith() {
        // Current: PSAS = Preemptive State Analytical Services (10/16, 39.80 Da)
        // Missing: ρ Recursion, ∅ Void, ∃ Existence, π Persistence, λ Location, Σ Sum
        //
        // Goal: Find words that fill the gaps or reframe with higher coverage.

        eprintln!("\n{}", "=".repeat(70));
        eprintln!(" WORDSMITH ENGINE: Naming Candidates");
        eprintln!("{}\n", "=".repeat(70));

        // ── Candidate words for the missing primitives ──
        let gap_fillers: Vec<(&str, Vec<LexPrimitiva>)> = vec![
            // Words that contain missing primitives
            (
                "Recursive",
                vec![LexPrimitiva::Recursion, LexPrimitiva::Sequence],
            ),
            (
                "Persistent",
                vec![LexPrimitiva::Persistence, LexPrimitiva::State],
            ),
            (
                "Validated",
                vec![LexPrimitiva::Existence, LexPrimitiva::Comparison],
            ),
            (
                "Void-Aware",
                vec![LexPrimitiva::Void, LexPrimitiva::Boundary],
            ),
            (
                "Localized",
                vec![LexPrimitiva::Location, LexPrimitiva::Boundary],
            ),
            ("Aggregate", vec![LexPrimitiva::Sum, LexPrimitiva::Quantity]),
            (
                "Comprehensive",
                vec![
                    LexPrimitiva::Recursion,
                    LexPrimitiva::Existence,
                    LexPrimitiva::Persistence,
                    LexPrimitiva::Sum,
                ],
            ),
            (
                "Vigilant",
                vec![
                    LexPrimitiva::State,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Causality,
                ],
            ),
            (
                "Sentinel",
                vec![
                    LexPrimitiva::Existence,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Persistence,
                    LexPrimitiva::State,
                ],
            ),
            (
                "Observatory",
                vec![
                    LexPrimitiva::Existence,
                    LexPrimitiva::Location,
                    LexPrimitiva::Persistence,
                    LexPrimitiva::Frequency,
                ],
            ),
        ];

        eprintln!("  Gap-filling words (missing: ρ ∅ ∃ π λ Σ):\n");
        for (name, prims) in &gap_fillers {
            let s = Specimen::new(name, prims.clone());
            let r = run_experiment(&s);
            let missing_hit: Vec<&str> = prims
                .iter()
                .filter(|p| {
                    matches!(
                        p,
                        LexPrimitiva::Recursion
                            | LexPrimitiva::Void
                            | LexPrimitiva::Existence
                            | LexPrimitiva::Persistence
                            | LexPrimitiva::Location
                            | LexPrimitiva::Sum
                    )
                })
                .map(|p| p.symbol())
                .collect();
            eprintln!(
                "    {:<16} [{:<6}] {:>5.2} Da | fills: {:?}",
                name, r.formula, r.molecular_weight, missing_hit,
            );
        }

        // ── Full naming candidates (4-5 word acronyms) ──
        eprintln!("\n  --- FULL NAME CANDIDATES ---\n");

        let candidates: Vec<(&str, &str, Vec<LexPrimitiva>)> = vec![
            (
                "PSAS",
                "Preemptive State Analytical Services",
                vec![
                    LexPrimitiva::Causality,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::State,
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Quantity,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Product,
                ],
            ),
            (
                "PRISM",
                "Preemptive Recursive Intelligence State Monitor",
                vec![
                    LexPrimitiva::Causality,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Recursion,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Comparison,
                    LexPrimitiva::State,
                    LexPrimitiva::Existence,
                ],
            ),
            (
                "SPEAR",
                "State Preemptive Existence Analytical Registry",
                vec![
                    LexPrimitiva::State,
                    LexPrimitiva::Causality,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Existence,
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Quantity,
                    LexPrimitiva::Persistence,
                ],
            ),
            (
                "AEGIS",
                "Analytical Existence Guard Intelligence Service",
                vec![
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Quantity,
                    LexPrimitiva::Existence,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::State,
                    LexPrimitiva::Causality,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Product,
                ],
            ),
            (
                "ATLAS",
                "Analytical Tracking Localized Alert System",
                vec![
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Quantity,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Location,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Product,
                ],
            ),
            (
                "VIGIL",
                "Validated Intelligence Guard Irreversibility Lattice",
                vec![
                    LexPrimitiva::Existence,
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::State,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Location,
                    LexPrimitiva::Recursion,
                ],
            ),
            (
                "SHIELD",
                "State Harm Intelligence Existence Localized Defense",
                vec![
                    LexPrimitiva::State,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Causality,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Existence,
                    LexPrimitiva::Location,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Persistence,
                ],
            ),
            (
                "ORACLE",
                "Observable Recursive Analytical Causal Lattice Engine",
                vec![
                    LexPrimitiva::Existence,
                    LexPrimitiva::Recursion,
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Quantity,
                    LexPrimitiva::Causality,
                    LexPrimitiva::Location,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Product,
                ],
            ),
            (
                "PHALANX",
                "Preemptive Harm Analysis Lattice Aggregate Nexus",
                vec![
                    LexPrimitiva::Causality,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Comparison,
                    LexPrimitiva::Mapping,
                    LexPrimitiva::Location,
                    LexPrimitiva::Sum,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Product,
                ],
            ),
            (
                "PRECEPT",
                "Preemptive Recursive Existence Causal Persistent Tracker",
                vec![
                    LexPrimitiva::Causality,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Boundary,
                    LexPrimitiva::Recursion,
                    LexPrimitiva::Existence,
                    LexPrimitiva::Persistence,
                    LexPrimitiva::State,
                    LexPrimitiva::Sequence,
                    LexPrimitiva::Mapping,
                ],
            ),
        ];

        for (acronym, full_name, prims) in &candidates {
            let s = Specimen::new(acronym, prims.clone());
            let r = run_experiment(&s);
            let coverage = prims.len();
            eprintln!(
                "  {:<8} {:>5.2} Da | {}/16 T1 | {:<8} | {:>2.0}% xfer | {}",
                acronym,
                r.molecular_weight,
                coverage,
                r.transfer_class,
                r.hybrid_transfer * 100.0,
                full_name,
            );
        }

        // ── Orthogonality check for top candidates ──
        eprintln!("\n  --- ORTHOGONALITY CHECK (word-by-word) ---\n");

        // PRECEPT decomposition
        let precept_words: Vec<(&str, Vec<LexPrimitiva>)> = vec![
            (
                "Preemptive",
                vec![
                    LexPrimitiva::Causality,
                    LexPrimitiva::Frequency,
                    LexPrimitiva::Irreversibility,
                    LexPrimitiva::Boundary,
                ],
            ),
            ("Recursive", vec![LexPrimitiva::Recursion]),
            ("Existence", vec![LexPrimitiva::Existence]),
            (
                "Causal",
                vec![LexPrimitiva::Causality, LexPrimitiva::Sequence],
            ),
            (
                "Persistent",
                vec![LexPrimitiva::Persistence, LexPrimitiva::State],
            ),
            (
                "Tracker",
                vec![LexPrimitiva::Sequence, LexPrimitiva::Mapping],
            ),
        ];

        eprintln!("  PRECEPT word-by-word:");
        let mut overlaps = 0u32;
        let mut seen = std::collections::HashSet::new();
        for (word, prims) in &precept_words {
            let new: Vec<&str> = prims
                .iter()
                .filter(|p| seen.insert(p.symbol()))
                .map(|p| p.symbol())
                .collect();
            let dup: Vec<&str> = prims
                .iter()
                .filter(|p| !new.contains(&p.symbol()))
                .map(|p| p.symbol())
                .collect();
            if !dup.is_empty() {
                overlaps += dup.len() as u32;
            }
            eprintln!(
                "    {:<14} new: {:<16} overlap: {:?}",
                word,
                format!("{:?}", new),
                dup,
            );
        }
        eprintln!("    Total overlaps: {} (PSAS has 0)", overlaps);

        // Assert basic properties
        assert!(candidates.len() >= 5);
    }
}
