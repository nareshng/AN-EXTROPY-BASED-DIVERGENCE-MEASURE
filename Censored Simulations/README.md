# Corrected code for Tables 7--10

This folder contains only the corrected right-censoring simulation code needed
to generate Tables 7--10 of the paper.

## Main command

Run a short software check first:

```bash
Rscript run_tables_7_10.R --mode=quick --cores=4 --tidy
```

Generate the manuscript-scale tables with:

```bash
Rscript run_tables_7_10.R --mode=full --cores=4 --tidy
```

Full mode uses 2,000 Monte Carlo replications in every simulation cell, 1,000
bootstrap resamples for Tables 9--10, point-estimation seed 5401, and
confidence-interval seed 5402. Quick-mode numbers are only for checking the
software and must not be reported in the paper.

With `--tidy`, the four manuscript CSVs are placed in `paper_tables/`; all
replication-level and reproducibility records are retained under `details/`.

No external R package is required. The code uses base R and can be run from any
working directory.

The exponential parameter pairs are rates, matching `rexp(rate = lambda)`.
Each Weibull tuple is `(shape, scale)`, matching R's `rweibull()` convention.

## Table mapping

| Table | Output file | Contents |
|---|---|---|
| 7 | `Table7_RelMSE_Exponential.csv` | Relative MSEs of the fixed-horizon KM estimators of `D_tau` and `D_CC,tau` for exponential lifetimes |
| 8 | `Table8_RelMSE_Weibull.csv` | Relative MSEs of the fixed-horizon KM estimators of `D_tau` and `D_CC,tau` for Weibull lifetimes |
| 9 | `Table9_CI_Exponential.csv` | Greenwood normal and diagnostic percentile-bootstrap intervals |
| 10 | `Table10_CI_Weibull.csv` | Greenwood normal and diagnostic percentile-bootstrap intervals |

Detailed long-format, replication-level, task-manifest, and checkpoint-manifest
files are also written. Tables 7--8 report relative mean squared error,
calculated separately against the target of each estimator:

```text
RelMSE_D   = mean{(D_hat - D_tau)^2} / D_tau^2
RelMSE_DCC = mean{(DCC_hat - D_CC,tau)^2} / D_CC,tau^2.
```
We deliberately removed `RelMSE_DCC` from the paper because the two measures are not directly comparable.



Each sample-size block therefore contains two columns, `D_KM` and `DCC_KM`.
Because the two estimators target different divergence functionals, their
relative-MSE values should not be interpreted as a direct efficiency comparison
between the two measures.






Tables 7--8 use the four sample-size pairs `(20,20)`, `(50,50)`, `(100,100)`,
and `(200,200)`. Tables 9--10 use the three sample-size pairs `(30,40)`,
`(70,50)`, and `(100,100)`.

## Restarting the long run

Every completed simulation cell is atomically saved under
`intermediate_results/point/` or `intermediate_results/ci/`. If a run stops,
repeat the same command with the same output directory and `--force`; compatible
cells will be reused.

To distribute individual cells, use the lower-level scripts. First inspect the
stable task map. Run these commands from this folder:

```bash
Rscript Point_estimation_right_censoring.R \
  --mode=full --B=2000 --seed=5401 --list-tasks

Rscript Confidence_Intervals_right_censoring.R \
  --mode=full --B=2000 --R-boot=1000 --seed=5402 --list-tasks
```

Then run selected task IDs into a common output directory, for example:

```bash
Rscript Point_estimation_right_censoring.R \
  --mode=full --B=2000 --seed=5401 --task-ids=1,4,9-12 \
  --output-dir=results/cells

Rscript Confidence_Intervals_right_censoring.R \
  --mode=full --B=2000 --R-boot=1000 --seed=5402 --task-ids=1 \
  --output-dir=results/ci_cells
```

After all tasks have been completed, construct the tables only from the complete
validated checkpoint cohort:

```bash
Rscript Aggregate_right_censoring_results.R \
  --analysis=point --mode=full --B=2000 --seed=5401 \
  --output-dir=results/cells

Rscript Aggregate_right_censoring_results.R \
  --analysis=ci --mode=full --B=2000 --R-boot=1000 --seed=5402 \
  --output-dir=results/ci_cells
```

Aggregation stops if any expected cell is missing, corrupt, or incompatible.
After a selected-task run, do not use any table CSV already present in that
output directory until the corresponding `--aggregate-only` command succeeds.

## Files

- `run_tables_7_10.R`: one-command master runner.
- `Point_estimation_right_censoring.R`: Tables 7--8 driver.
- `Confidence_Intervals_right_censoring.R`: Tables 9--10 driver.
- `Aggregate_right_censoring_results.R`: strict checkpoint-only aggregation.
- `km_functions.R`: Kaplan--Meier estimator, Greenwood variance, bootstrap, and
  censoring calibration.
- `censored_simulation_helpers.R`: scenario grids, simulations, diagnostics,
  tables, task selection, and checkpoints.
- `diagnose_tables.R`: standalone contract check for a completed output folder.
- `atomic_io_functions.R`: safe atomic CSV and RDS writing.


