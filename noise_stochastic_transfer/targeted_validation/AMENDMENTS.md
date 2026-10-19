# Computational amendments

No computational amendments at protocol creation.

## High-precision root preparation

The first deterministic pass solved the six pair/symmetric roots and stalled on the unequal triple. Before any fresh draws, the implementation was repaired to derive h_r from the exact decimal weights/offsets rather than a binary-float normalization, use the existing unequal-triple crossing as a numerical starting value, and allow more Newton steps. The initial partial roots are retained in roots_first_attempt.json. No design points or decision-rule parameters changed.

## Fresh generator coordinate repair

The first 36,000-draw run used iota(x) as the complex mean. That is not the specified Hermitian observation model and yielded spurious signal in the anti-Hermitian noise estimator. The entire invalid run and its summaries are retained in results/invalid_generator_run. The generator now applies the inverse Hermitian coordinate map before adding circular complex noise and asserts Hermitian symmetry. All cells are rerun with the same fixed design and seeds. The old-record operating curves and independent root/record checks did not use this generator and are unchanged.

## Direct cost-test matched-error summary

The fixed protocol includes C in every operating curve but specifies paired matched-error intervals for S versus I. During manuscript assembly, the same five target budgets and interpolation were applied to S versus C as an additional secondary summary of the already planned curves. A separate paired within-cell bootstrap uses master seed 2026101910 (plus dataset index), 399 resamples. This addition was made after observing outcomes, adds no draws, changes no rule or default, and is not described as pre-specified.
