# Evaluation protocol: regime-selection rule

**Fixed:** 2026-10-19T17:46+11:00, before any evaluation draw was generated.

**Hashes.** The SHA-256 of this file, of `THEOREM.md` and of `code/run_study.py` are recorded in `PROTOCOL.sha256`, written at the same time.

**Changes after freezing.** Any change to the rule, the constants or the design after this point is a deviation. It is reported in the results as exploratory and evaluated on separate draws with a different master seed.

## 1. Rules compared (all on the same draws)

| Label | Rule |
|---|---|
| **S** (primary) | Regime selection. C is selected iff Â > 0, m̂₀ > 0 and ν̂ = σ̂/(c_ν Â) ≤ (Â/m̂₀)^γ; otherwise L. Then δ_S = δ_{selected}. |
| I | The original agreement rule: declare only if L and C declare the same branch. |
| C | The cost test alone, at level α. |
| L | The linear test alone, at level α. |
| C½ | The cost test at the fixed stricter level α/2 (q = t_{21, 0.975}). |

**Constants:**
- α = 0.05, q = t_{21, 0.95}.
- γ = 1/(4r), that is 1/4 for K = 2 and 1/8 for K = 3. This is the midpoint of the admissible range (0, 1/(2r)) of `THEOREM.md`.
- c_ν = ‖P q_r‖/√N, the noise unit already used to define ν: 0.0370563141 (K = 2) and 1.916883326·10⁻⁴ (K = 3).
- λ_* is the certified value. The arcs are as before.

**How the constants were chosen.**
- No constant was chosen from simulation output.
- The nominal switch points ν ≈ (A/μ₀)^γ were computed from the design parameters only, to place the transition cells: 0.55 (pair, z = 0.3), 0.45 (symmetric triple, z = 0.3), 0.51 (unequal triple, z = 0.3), 0.36 (unequal triple, z = 0.15).
- The earlier records (`reliability_followup/results`) motivated the change of rule. They are not used for evaluation.

**Implementation.**
- `code/run_study.py` reuses `reliability_followup/code/reliability.py` without modification:
  - restricted Hermitian fits by arc profiling plus Levenberg–Marquardt;
  - the complex global fit by full-circle profile plus multistart.
- The fits are numerical optima, not certified global optima.
- The noiseless branch L(x) is the label of the smaller numerical restricted noiseless cost.
- λ_c(z) is taken from the earlier deterministic bisection (`validation_raw.json`, field `lambda_c`).

## 2. Design

**Model.** Model (1), with N = 21 and b = 10.

**Cluster cases** (unchanged families):
- equal pair at ±1 (K = 2, z = 0.3);
- equal triple at {−1, 0, 1} (K = 3, z = 0.3);
- triple with weights (1, 2, ½) at (−1, 0.3, 1.5) (K = 3, z = 0.3 and z = 0.15).

**Noise levels:** ν ∈ {0.005, 0.02, 0.1, 0.35, 0.5, 0.7, 1.0}.
- weak noise: 0.005, 0.02, 0.1;
- transition, where selection is expected to be random: 0.35, 0.5, 0.7;
- critical noise: 1.0.

**Amplitude positions**, with w = ρ√2‖d₀‖/Φ′:
- λ = λ_c(z) + x·w for x ∈ {−2, −0.5, −0.1, +0.1, +0.5, +2}, which covers both sides of the crossing;
- λ = (λ_* + λ_c)/2.

At x = ±2 and ν = 1, λ lies up to about 50% from λ_*. The truth label there is still the numerical noiseless label.

**Size.** 4 × 7 × 7 = 196 cells, with 1000 independent draws per cell (196,000 draws in total).

**Seeds.**
- Each cell is split into 20 chunks of 50 draws.
- Chunk k of cell i (cells in the order case × ν × position) uses `np.random.default_rng(SeedSequence([2026101902, i, k]))`.
- A smoke test (`--smoke`: 2 cells, seed word 1) may be run first, to check execution only. Its output file is labelled `SMOKE_TEST_DISCARD`, is not analysed, and its decisions are not inspected.

## 3. Outcomes and reporting

All outcomes are computed per cell. Cell-specific results are primary. Pooled summaries are secondary and labelled as pooled.

**Calibration (primary).** For each rule:
- the wrong-declaration count and rate, with a Wilson 95% interval;
- the cell is flagged if the lower bound exceeds α;
- with 196 cells, a rule whose rate is exactly α near the tie is expected to be flagged in a few cells by chance, so a Bonferroni-adjusted flag is also reported (Wilson interval at level 1 − 0.05/196);
- results are reported separately by side of the crossing (truth A or B).

**Declarations.** For each rule:
- correct declarations, wrong declarations, abstentions;
- the declaration rate;
- wrong declarations among declarations, with the denominator and a Clopper–Pearson 95% interval; it is reported as undefined when there are no declarations.

**Selection.** The number of draws in which C or L is selected, and the number of exceptional selections (Â ≤ 0 or m̂₀ ≤ 0).

**Paired comparisons (same draws).**
- **S against I.**
  - Draws where S declares correctly and I abstains (restored correct declarations).
  - Draws where S declares wrongly and I abstains (added wrong declarations).
  - A check that I ⊂ S holds draw by draw, as Proposition 3(a) states.
- **S against C.** The four discordance counts: correct only under S, correct only under C, wrong only under S, wrong only under C.
- **S against C½ and L.** Correct and wrong counts.
- Each paired count is reported with its Wilson 95% interval as a fraction of the 1000 draws.

**Usefulness (primary).**
- Correct declarations of S against I and against C in every cell.
- Pre-specified focus cells:
  - weak noise (ν ≤ 0.1) on the side of the crossing opposite to the linear test's finite-spacing offset (for the symmetric triple, truth A: x = −0.5, −0.1 and the between position);
  - the transition cells (ν = 0.35, 0.5, 0.7).

**Fitted optimizer (separate).** For each rule:
- draws where it declares and the complex global fit lies in the declared arc;
- the number of those where the fitted branch is not L(x);
- the fraction of fits outside both arcs.

These are not branch declarations and are reported separately.

**Diagnostics.**
- The agreement between Â from the selector and from the linear statistic.
- The distribution of ν̂/ν.

## 4. What will not be done

- No change of γ, c_ν, q or the design in response to these results within the primary analysis.
- No conditional-correctness claims.
- No new noise models, large N, or unknown b.
