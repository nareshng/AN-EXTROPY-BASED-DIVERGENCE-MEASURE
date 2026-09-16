# Uncensored simulations (Sections 5.1-5.3, Tables 1-6)

Base R only. No contributed package is required: the empirical-likelihood ratio
used by the JEL interval is computed in `el_functions.R` rather than by
`emplik`, and `el_self_test()` cross-checks against `emplik::el.test()` when
that package happens to be installed.

Tested with R 4.3.3 on Linux. `parallel` (shipped with R) is used when
`--cores` is greater than 1.

## One command per table

```bash
Rscript Sim_RelativeMSE_comparison_Exponential_dist.R --output-dir=results   # Table 1
Rscript Sim_RelativeMSE_comparison_Weibull_dist.R     --output-dir=results   # Table 2
Rscript Sim_MSE_Estimators_Exponential_dist.R         --output-dir=results   # Table 3
Rscript Sim_MSE_Estimators_Weibull_dist.R             --output-dir=results   # Table 4
Rscript Sim_Confidence_Intervals_Exponential_dist.R   --output-dir=results --cores=4   # Table 5
Rscript Sim_Confidence_Intervals_Weibull_dist.R       --output-dir=results --cores=4   # Table 6
```

or all six at once:

```bash
Rscript run_uncensored_tables.R --output-dir=results --cores=4
Rscript run_uncensored_tables.R --output-dir=check --cores=4 --quick   # software check
```

Add `--quick` to any script for a fast software check; quick-mode numbers must
not be reported in the paper.

Options accepted by the two confidence-interval scripts:

| Option | Meaning |
|---|---|
| `--iterations=N` | Monte Carlo replications per cell (default 2000) |
| `--seed=N` | base seed (default 2026) |
| `--cores=N` | forked workers over simulation cells; results do not depend on this |
| `--cells=1,4,9-12` | run only these cells, checkpoint them, and write the table once every cell exists |
| `--no-resume` | recompute selected cells instead of reusing their checkpoints |
| `--with-kernel` | also compute the kernel bootstrap interval (see below) |
| `--B-boot=N` | kernel bootstrap resamples, used only with `--with-kernel` |
| `--output-dir=DIR`, `--quick`, `--quiet` | as above |

### The kernel bootstrap interval is off by default

Tables 5 and 6 report the JEL, normal-approximation and empirical intervals.
The fourth interval in the code, a percentile bootstrap around the kernel
estimator, costs roughly twenty times the other three put together and appears
in no table, so it is computed only under `--with-kernel`; without it the
`Ker_CP` and `Ker_AL` columns are `NA`. With the default settings Tables 5 and 6
each take about five minutes on two cores; with `--with-kernel` they take
closer to two hours.

### Long runs

Every cell is written to `intermediate_exp/` or `intermediate_wei/` inside the
output directory as soon as it finishes, keyed by everything that could change
its numbers (cell, replications, bootstrap size, quadrature settings, alpha,
seed, RNG kinds). Re-running the same command reuses those checkpoints, and a
run split into `--cells=1-6` then `--cells=7-20` produces a table bit-identical
to the single-command run. A partial run writes no table and reports which
cells are still missing.

| Table | Script | Output file | Settings |
|---|---|---|---|
| 1 | `Sim_RelativeMSE_comparison_Exponential_dist.R` | `Table1_RelMSE_Exponential.csv` | B = 2000, seed 2024, `n_quad` = 100 |
| 2 | `Sim_RelativeMSE_comparison_Weibull_dist.R` | `Table2_RelMSE_Weibull.csv` | B = 2000, seed 2024, `n_quad` = 100 |
| 3 | `Sim_MSE_Estimators_Exponential_dist.R` | `Table3_MSE_Exponential.csv` | B = 2000, seed 2026, `n_quad` = 80 |
| 4 | `Sim_MSE_Estimators_Weibull_dist.R` | `Table4_MSE_Weibull.csv` | B = 2000, seed 2024, `n_quad` = 80 |
| 5 | `Sim_Confidence_Intervals_Exponential_dist.R` | `Table5_CI_Exponential.csv` | B = 2000, seed 2026, 499 kernel bootstrap resamples |
| 6 | `Sim_Confidence_Intervals_Weibull_dist.R` | `Table6_CI_Weibull.csv` | B = 2000, seed 2026, 199 kernel bootstrap resamples |

Tables 1-4 run in a few minutes each on one core. Tables 5-6 are the long ones;
budget roughly two hours each on one core, proportionally less with `--cores`.
Each confidence-interval run also writes a `TableN_run_settings.txt` recording
the replication count, seed, core count, elapsed time, R version and RNG kinds.

## The JEL constraint

Equation (3.5) defines the jackknife pseudo-values as

```text
V_i = n T_n - (n - 1) T_{n-1}^{(-i)},        n = n1 + n2.
```

`T_n` and `T_{n-1}^{(-i)}` are two-sample U-statistics with the same kernel,
evaluated on (n1, n2) and on (n1 - 1, n2) or (n1, n2 - 1) observations, so both
are unbiased for D. Hence

```text
E[V_i] = n D - (n - 1) D = D           for every i, in both blocks,
```

with no dependence on n1, n2 or on which sample observation i came from. The
constraint in the empirical likelihood (3.7) is therefore

```text
sum_i p_i (V_i - theta) = 0,
```

which is what `EV_vec()` returns. Any expression for `E[V_i]` that depends on
n1 and n2 separately does not describe the pseudo-values of (3.5): for
(n1, n2) = (30, 40) a group-specific rule of the form
`(n/(n-2)){(n2-1)(2/n1) - 1}` gives 1.647 and 0.463 for the two blocks, whereas
simulation of the pseudo-values themselves returns 1.000 for both.

Note also that the variance estimator (3.2) uses the *group-specific*
pseudo-values `n_r U - (n_r - 1) U^{(-i)}`, which is a different and equally
standard construction; the code keeps the two separate, as the paper does.

## Reproducibility notes

- Every Section 5.3 replication is seeded from `(base seed, cell id,
  replication index)`, so a cell's result does not depend on the order in which
  cells are run or on `--cores`. `--cores=1` and `--cores=4` give bit-identical
  CSVs.
- Sections 5.1 and 5.2 keep the original single-stream seeding, so Tables 1-4
  are unchanged by this reorganisation; only the output path and file name moved.
- `Sim_RelativeMSE_comparison_Exponential_dist.R` carries an optional quadrature
  sensitivity check. It is off by default and is enabled with
  `--quadrature-check`; it roughly triples the runtime.

## Files

- `el_functions.R` - empirical-likelihood ratio in base R, with self-tests.
- `unc_simulation_engine.R` - per-replication seeding, cell grid, optional
  forked parallelism, command-line parsing and atomic CSV writing for Tables 5-6.
- `run_uncensored_tables.R` - runs all six scripts into one output directory.
- the six table scripts listed above.
