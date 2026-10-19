# Regime-adaptive branch decisions for deficient-order fits of spectral clusters

Code and simulation data for the statistical-method paper by Jinghang Mei, prepared for *Signal Processing*.

The method returns branch A, branch B or an inconclusive decision. It selects between a linear contrast test and a restricted fitting-cost test using estimated noise and cluster defect. The theoretical target is a uniform asymptotic bound on the unconditional probability of a wrong declaration. The numerical study compares correct declarations, wrong declarations and abstentions across noise regimes.

**Version 1.0.0 archive: [DOI 10.5281/zenodo.23474289](https://doi.org/10.5281/zenodo.23474289)**  
[GitHub Release v1.0.0](https://github.com/surnamemei/regime-adaptive-spectral-decisions/releases/tag/v1.0.0)

## This archive and the deterministic companion

This repository contains the **statistical decision study**: the tests, adaptive selector, simulation records, operating comparisons, finite-spacing checks and figure generators. Its first software release is **v1.0.0**.

The deterministic localization theorems and interval certificates belong to the separate ACHA companion, [Global Representation Exchange in Best Approximation of Positive Spectral Clusters](https://github.com/surnamemei/positive-spectral-cluster-exchange), archived as [DOI 10.5281/zenodo.23451747](https://doi.org/10.5281/zenodo.23451747). That DOI identifies the companion's material. Statistical simulations here are numerical experiments, rather than interval-certified parameter regions.

## Contents

| Location | Scientific content |
|---|---|
| `noise_stochastic_transfer/adaptive_selection/` | Original 196-cell study: 196,000 per-draw records, study generator and protocol |
| `noise_stochastic_transfer/reliability_followup/` | Preliminary validation records supplying finite-spacing thresholds and comparison calculations |
| `noise_stochastic_transfer/targeted_validation/` | Fixed follow-up design: 72 cells and 36,000 draws; matched-error curves; spacing and tuning checks; high-precision thresholds and labels |
| `paper_figures/` | Paper figure/table inputs and the original-study figure generator |
| `tools/` | Integrity and stored-data checks |

The follow-up includes three spacings, three noise levels, four positions relative to the crossing and two cluster families. Correct/wrong counts are compared at common design-average error budgets. The fixed protocols, seeds, design weights and computational amendments are included. `PROTOCOL.sha256` records file integrity; it is not a third-party timestamp.

## Quick check

Python 3.12 was used. Create an environment and install the pinned dependencies:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python tools/verify_release.py
```

The check verifies the release manifest, record counts and default fresh-study totals. It does not run new simulations.

## Reproduce the paper displays

Use a working copy so regenerated outputs can be compared with the supplied files:

```sh
export OPENBLAS_NUM_THREADS=1
python paper_figures/make_figures.py
python noise_stochastic_transfer/targeted_validation/code/make_outputs.py
```

The first script produces the original-study displays; the second produces the matched-error, spacing and tuning displays. Figure metadata may include creation timestamps; scientific tables can be compared directly.

For the operating curves and paired bootstrap summaries, run from `noise_stochastic_transfer/targeted_validation/`:

```sh
python code/operating.py old
python code/operating.py fresh
python code/matched_cost.py old
python code/matched_cost.py fresh
```

`README.md` in that directory lists the full fresh-study, high-precision and selected-record replay commands. The raw original-study JSON is about 65 MB; the release ZIP compresses the data for download.

## Model and interpretation

The study uses N=21, a known remote frequency b=10, a known phase reference, a phase-opposed remote tone and circular Gaussian noise. The nominal error target is unconditional across repeated records; error among declarations is reported separately. Matched-error comparisons are design averages. Finite-spacing results retain both gains and losses against the agreement rule.

## Licence and citation

Code: MIT. Data, figures and documentation: CC BY 4.0, as specified in `LICENSE-DATA.md`. Release citation metadata are in `CITATION.cff` and `.zenodo.json`.
