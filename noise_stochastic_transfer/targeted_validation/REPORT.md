# Targeted validation results

## Design

The protocol/configuration were fixed before the follow-up computations and fresh draws. Keep N=21, b=10, known phase reference, phase-opposed remote tone, positive clusters and circular Gaussian noise. Analyze the existing 196,000 draws separately from 36,000 new draws. The theorem, default rule and model are unchanged.

## Fair operating comparisons

Compare S (adaptive selection), C (cost only), I (agreement) and L (linear only) across ten nominal levels. Match by linear interpolation at design-average wrong-declaration budgets, with no extrapolation. Average within each stated design. Paired bootstrap intervals use 399 resamples within each cell. S/I matching was specified in the protocol; direct S/C scalar matching was added during assembly of the planned curves and is documented in AMENDMENTS.md.

At 2% average wrong declarations:

| Design/noise | S correct | C correct | S minus C (percentage points), 95% interval | I correct |
|---|---:|---:|---:|---:|
| Original, all | 26.01% | 22.85% | +3.16 [+2.96, +3.40] | 26.87% |
| Original, nu=0.005 | 29.59% | 29.59% | +0.00 [+0.00, +0.00] | 30.63% |
| Original, nu=0.5 | 26.37% | 21.32% | +5.05 [+4.34, +5.70] | 26.19% |
| Original, nu=1.0 | 29.46% | 20.72% | +8.75 [+7.93, +9.59] | 24.90% |
| Fresh, all | 38.95% | 31.95% | +7.00 [+6.30, +7.65] | 36.56% |
| Fresh, nu=0.005 | 32.45% | 32.45% | +0.00 [+0.00, +0.00] | 34.64% |
| Fresh, nu=0.5 | 41.67% | 32.69% | +8.99 [+7.65, +10.53] | 39.24% |
| Fresh, nu=1.0 | 42.56% | 30.77% | +11.79 [+10.41, +13.25] | 36.02% |

S improves on C at transition/critical noise and coincides with C at weak noise in these designs. Relative to I, the pooled difference is -0.86 points originally and +2.39 in the fresh study; critical-noise gains accompany weak-noise losses. The fresh design assigns 1/3 of cells to critical noise and 1/2 to far positions, versus 1/7 and 2/7 originally. The pooled ranking against I therefore depends on design weights.

## Difficult thresholds and labels

Eight equal-cost stationary roots were recomputed from direct 21-sample regression at 70/100 digits. Their maximum difference is below 2e-55. Both profiled Hessians are positive at each root and the cost difference has opposite signs on its two sides. All 12 preselected truth labels agree. All 24 reconstructed records retain S and I decisions; maximum statistic changes are 5.68e-5 (tL) and 1.92e-5 (TC). The three symmetric-triple crossings lie in the existing bridge tube, with predictor discrepancies below 1.54e-15. These are restricted-fit stationary/numerical checks, not new interval certificates for noisy global minima.

## Spacing and tuning

Two families times three spacings (.3,.15,.075), three noise levels (.005,.5,1), four positions (-2,-.1,.1,2 widths), 500 draws/cell: 72 cells and 36,000 draws. There were no numerical failures. Default counts: 14,149 correct, 741 wrong, 21,110 inconclusive. Maximum default cell error: 29/500=.058, Wilson 95% interval [.0407,.0821].

Nine pre-fixed settings on the same records give 12,010-14,361 correct and 709-759 wrong declarations; maximum individual wrong rate .060. The default was not tuned from outcomes. At transition noise, C-selection falls from 52.55% to 5.825% to .10% as spacing decreases. The scaled noiseless threshold displacement approaches a nonzero constant in both symmetric families.

All 36 new triple cells, like the 147 original triple cells, are outside the stated reduced lambda box. Pair interval membership is undetermined and the spacing radius remains existential. The study evaluates finite-spacing numerical performance separately from the asymptotic theorem.

## Reproduction and manuscript

All 72 mean-generation identities/noiseless labels pass; two fixed 50-record chunks replay exactly. FINAL_VERIFICATION.json records the checks. The author AI declaration is unchanged, as are the ten formal statements after deleting one redundant sentence on sufficiency. The four ACHA source/PDF hashes match the baseline.

The main article contains operating curves, a matched-error table, a spacing figure and concise results. Supplement S7 gives high-precision calculations, tolerances, all tuning settings and every fresh cell. AMENDMENTS.md records root-preparation and generator repairs. The invalid first generator run is retained internally and excluded from scientific summaries.
