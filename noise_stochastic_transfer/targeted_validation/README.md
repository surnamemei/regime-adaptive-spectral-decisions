# Targeted validation

Read REPORT.md, then PROTOCOL.md/config.json. AMENDMENTS.md records computational repairs and the added S/C matched-error summary. The stored roots, label checks and selected-record replay results are in results/.

## Reproduction

Python 3.12; numpy, scipy, matplotlib, mpmath. Executed versions: NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, mpmath 1.3.0. Install mpmath with the other pinned dependencies.

From this directory, set OPENBLAS_NUM_THREADS=1 and run:

```
python code/high_precision.py
python code/audit_records.py
python code/fresh_study.py
python code/check_design_means.py
python code/operating.py old
python code/operating.py fresh
python code/matched_cost.py old
python code/matched_cost.py fresh
python code/make_outputs.py
python code/replay_selected.py
```

The fresh study uses 12 workers and took about 214 seconds. High-precision checks take longer. Scripts write into results/ and figures/; run in a copy to compare with supplied outputs. Preserve the relative noise_stochastic_transfer layout. Inputs are adaptive_selection/results/raw.json and reliability_followup/results/validation_raw.json; both are included in the local archive. Their original protocol and source files are included, with relative dependencies preserved.

- code/: direct high-precision regression, stable restricted fitting, fresh draws, curves, plotting and a 100-record replay.
- results/: valid outputs; invalid_generator_run/ is an invalid internal diagnostic run, excluded from the scientific archive.
- figures/: generated figures and TeX tables.
- manuscript_before/: manuscript snapshots preceding this pass, excluded from the scientific archive.

Unattainable matched-error targets are NaN, with no extrapolation. Wilson intervals are per-cell, curve intervals are paired bootstrap intervals, and spacing-figure bars use within-cell Monte Carlo variances. Seeds, design weights and all nine tuning settings are reported. The repository and its versioned release contain these code and data. Typesetting sources for the article are separate.
