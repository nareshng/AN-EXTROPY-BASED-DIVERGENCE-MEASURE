# What changed, and why

Everything here is a change to `Uncensored Simulations/`. `Censored
Simulations/` and `Real Data Analysis/` are untouched.

## 1. The JEL constraint is now centred at E[V_i] = D

`EV_vec()` returned group-specific constants

```r
cx <- (n / (n - 2)) * ((n2 - 1) * (2 / n1) - 1)
cy <- (n / (n - 2)) * ((n1 - 1) * (2 / n2) - 1)
```

matching equation (3.8) of the manuscript. Those constants do not describe the
pseudo-values the code actually builds. By (3.5),

```text
V_i = n T_n - (n - 1) T_{n-1}^{(-i)},
```

and both `T_n` and `T_{n-1}^{(-i)}` are two-sample U-statistics with the same
kernel, unbiased for D on (n1, n2) and on (n1 - 1, n2) or (n1, n2 - 1)
observations. So `E[V_i] = n D - (n - 1) D = D`, for every i, in both blocks.
Simulating the code's own pseudo-values (B = 3000, exponential (0.5, 1)):

| (n1, n2) | E[V_X]/D | c_x from (3.8) | E[V_Y]/D | c_y from (3.8) |
|---|---|---|---|---|
| (20, 20) | 0.9872 | 0.9474 | 0.9872 | 0.9474 |
| (50, 50) | 0.9882 | 0.9796 | 0.9882 | 0.9796 |
| (30, 40) | 1.0001 | 1.6471 | 1.0001 | 0.4632 |
| (10, 30) | 0.9909 | 5.0526 | 0.9909 | -0.4211 |
| (40, 10) | 0.9757 | -0.5729 | 0.9757 | 7.0833 |

`EV_vec()` now returns `rep(theta, n1 + n2)`. Tables 5 and 6 must be replaced
with the regenerated values in `expected_outputs/`; equation (3.8) of the
manuscript must be replaced by `E[V_k] = theta` (draft wording in
`Section3_equation_3_8.md`).

Effect on the tables: the corrected constraint gives slightly shorter intervals
and slightly lower coverage at the small and unequal sample sizes, and the same
values at (100, 100). It also brings the JEL and normal-approximation columns
close together, which is what first-order equivalence predicts; under the old
constants JEL sat 1 to 4 points above NA throughout, which was an artefact of
the inflated constraint rather than a property of the method.

## 2. `emplik` is no longer required

`el_functions.R` computes the empirical-likelihood ratio in base R through the
Lagrange dual: with `lambda` solving `sum_i w_i / (1 + lambda w_i) = 0`, the
maximiser is `p_i = 1 / {n (1 + lambda w_i)}` and `-2 log R = 2 sum_i log(1 +
lambda w_i)`. `el_self_test()` checks that the dual solution satisfies both
constraints exactly and attains the maximum, and, when `emplik` happens to be
installed, that the two agree to 1e-6. Each confidence-interval script runs
that self-test before simulating.

One detail matters in practice. `-2 log R` is nonnegative in theory but comes
back as roughly `-1e-16` when the constraint is already satisfied at
`lambda = 0`, which is exactly what happens at the starting point of the
interval search. A guard that treats a negative value as a failure makes the
search report the interval as unavailable about half the time. `el_m2llr()`
clamps at zero.

## 3. Section 5.3 is seeded per replication

Each replication is now seeded from `(base seed, cell id, replication index)`
rather than drawing from one long stream. A cell's result therefore does not
depend on the order in which cells run or on how many cores are used:
`--cores=1` and `--cores=2` give bit-identical CSVs, and a run split into
`--cells=1-6` then `--cells=7-20` reproduces the single-command run exactly.
Cells are checkpointed as they finish, so a long run can be stopped and resumed.

## 4. The kernel bootstrap interval is opt-in

Tables 5 and 6 report three intervals. The code also computed a fourth, a
percentile bootstrap around the kernel estimator, which costs roughly twenty
times the other three put together and appears in no table. It now runs only
under `--with-kernel`. Tables 5 and 6 take about five minutes each on two cores
instead of about two hours.

## 5. Fixes and housekeeping

- `Sim_RelativeMSE_comparison_Exponential_dist.R` ended by comparing
  `results_100` with `results_150`, which the script never creates, so it always
  exited with an error after writing its table. The block now compares
  `n_quad = 200` against `n_quad = 300`, and is opt-in with `--quadrature-check`.
- ` Sim_Confidence_Intvals_Exponential_dist.R` (leading space, "Intvals") is
  renamed `Sim_Confidence_Intervals_Exponential_dist.R`, which is the name the
  root README already used.
- Every script takes `--output-dir=` and writes `TableN_...csv`; the two
  confidence-interval scripts also write a `TableN_run_settings.txt` recording
  replications, seed, cores, elapsed time, R version and RNG kinds.
- `run_uncensored_tables.R` runs all six into one directory and fails loudly if
  any output is missing or empty.
- `--quick` everywhere for a software check; the whole suite runs in about a
  minute.

## 6. Tables 1 to 4 are unchanged

The single-stream seeding of Sections 5.1 and 5.2 was deliberately left alone,
so Tables 1, 2 and 4 keep the values printed in the manuscript. Verified by
running before and after: all three CSVs are bit-identical.

Table 3 is the exception, and it was already the exception before these
changes: `Sim_MSE_Estimators_Exponential_dist.R` run unmodified does not
reproduce the Table 3 in the manuscript, differing by up to 7.7% in both
directions, for instance 3.9122 against 3.7562 and 0.9883 against 1.0711 in the
kernel row of (0.5, 0.1). The regenerated table is in
`expected_outputs/Table3_MSE_Exponential.csv`.
