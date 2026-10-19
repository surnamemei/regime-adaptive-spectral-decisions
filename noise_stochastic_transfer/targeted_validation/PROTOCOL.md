# Targeted validation protocol

Fixed at 2026-10-19T21:28:40.973361+11:00 (Australia/Sydney), before the new analyses and draws described below. This is a follow-up to the existing 196-cell study; the earlier results and reviewer concerns informed the questions. It is not a retrospective preregistration of that earlier study.

## Model and fixed decision rules

Keep N=21, b=10, positive common-phase clusters, phase-pi remote tone, centered grid and circular white Gaussian noise. Use the existing equal pair, symmetric equal triple, arcs, lambda-star values and Gauss-defect normalization. No larger N, K, phase uncertainty or new model is introduced. The existing theoretical statements remain unchanged.

Default selector: gamma=1/(4r), multiplier 1. Sensitivity: gamma in {1/(8r),1/(4r),3/(8r)} and multiplier in {0.5,1,2}. All nine settings use the same generated records. No setting will be chosen as a new default from its performance.

## A. Fair comparison

First re-evaluate all stored 196,000 records on the fixed nominal-level grid in config.json, using the saved t_L, T_C and selector. This is a secondary analysis of existing draws. Display probability of correct declaration versus probability of wrong declaration for S,I,C,L. Also report abstention and conditional wrong/declaration rates. Groups: all cells, and each of the seven existing noise levels, with equal weight for each cell in its group.

At target wrong-declaration rates {0.005,0.01,0.02,0.03,0.05}, compare S and I by piecewise-linear interpolation of adjacent observed curve points; retain the highest correct rate at repeated error coordinates. Do not extrapolate outside observed ranges. These are design-average operating-curve comparisons, not a new calibrated procedure or a bound uniform over cells. Also report correct/wrong/abstention at the original alpha=0.05.

Use 399 paired stratified bootstrap resamples, resampling draws within each fixed cell and preserving the joint S/I outcomes over all cutoffs. Report percentile intervals on the matched-error difference when at least 95% of bootstrap curves reach the target; otherwise mark interval unavailable. Fixed bootstrap seed is in config.json. Apply the same analysis to the fresh study.

## B. Difficult thresholds and labels

Recompute the noiseless crossing for eight cluster/spacing pairs: pair and symmetric triple at z=0.3,0.15,0.075, plus the existing unequal triple at z=0.3,0.15. Obtain stationary wells and equal costs by direct sampled regression in mpmath at 70 and 100 digits, independent of the double-precision profiling implementation and without Fourier-kernel formulas. Report root agreement, residual-gradient norms, Hessian signs and the scaled offset (lambda_c-lambda_L)/z^kappa. Obtain lambda_L as the zero of the noiseless regression contrast. Check both signs of the cost difference around each root.

For the existing four cases, audit the two closest nonzero standardized positions (-0.1,+0.1 widths) at nu=0.005 and the between-threshold position at nu=0.5: twelve preselected cells. Reconstruct old draw indices 0 and 1 from their recorded seed construction, giving 24 noisy records. Refit both restricted wells at high precision, compare cost statistics and declarations, and report discrepancies. These checks establish numerical agreement and stationary minima, not interval-certified global optimality. Use the existing symmetric-triple bridge only where its stated band and normalization apply; do not infer coverage of all simulation cells from it.

## C. Fresh small-spacing and tuning study

Two existing cluster families x three spacings (0.3,0.15,0.075) x three relative noise levels (0.005,0.5,1) x four positions (-2,-0.1,+0.1,+2 noise widths from the newly recomputed crossing). This is 72 cells with 500 draws each (36,000 fresh draws). Width is the same reduced-contrast width as in the earlier study. Truth is the side of the independently recomputed noiseless crossing, cross-checked against restricted costs. Absolute sigma is nu*c_nu*h_r*z^(2r).

Use master seed 2026101907, independent of earlier seeds, with SeedSequence([master,cell_index,chunk_index]); chunks have 50 draws. Compute only the two restricted Hermitian fits needed for the branch tests, not the auxiliary complex global-fit benchmark. Scale the residual objective by a positive data-derived norm to avoid absolute optimizer stopping at a tiny residual. This scaling leaves its minimizers unchanged. Record solver failures and convergence quantities; retain all draws, marking failed numerical decisions inconclusive and counting failures separately.

At each spacing/noise/family report correct, wrong and abstaining probabilities, conditional error among declarations, cost-test selection frequency and Wilson 95% intervals for individual-cell proportions. The fresh study also supplies the fixed operating curves in A. Sensitivity results cover all nine fixed settings, including negative or null findings.

## Reporting and stopping

Do not add cells, tune defaults, replace inconvenient results or extend the model after inspecting outcomes. Computational bugs may be repaired and recorded in AMENDMENTS.md with the effect on data; rerun affected outputs. Numerical failure is an outcome, not an exclusion criterion.

Main text: the operating comparison and the spacing/tuning findings that affect interpretation. Supplement: full cells, tolerances, numerical checks and implementation details. Keep mathematical scope conditions and statistical definitions once where they apply. Delete sentences whose sole job is rebutting an imagined objection.
