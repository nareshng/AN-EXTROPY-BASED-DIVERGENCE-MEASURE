# Reproducing Tables 1--6

This folder contains the six standalone simulation scripts used for manuscript
Tables 1--6 and a single deterministic runner. The statistical functions remain
inside their scripts; the runner only supplies replication count, seed,
verbosity, and output location through environment variables.

## Requirements

- R >= 4.1.0 (run with `Rscript --vanilla`)
- the pinned R package `emplik` 1.3-3 for Tables 5--6

The archival reference environment is R 4.6.1 with `emplik` 1.3-3. Install the
pinned package and its declared dependencies with:

```r
install.packages("remotes")
remotes::install_version(
  "emplik",
  version = "1.3-3",
  repos = "https://cloud.r-project.org",
  upgrade = "never"
)
```

The runner stops if a different `emplik` version is used and warns if R differs
from 4.6.1. It records the actual versions in `run_manifest.csv` and
`sessionInfo.txt`. The validator, all six quick simulations, output checks, and
the two-run determinism test were executed successfully with R 4.6.1 and
`emplik` 1.3-3.

## Validate the source files

From this directory, run:

```bash
Rscript validate_tables_1_6.R
```

The validator checks that every script parses, uses the shared runner interface,
contains no machine-specific working directory, and names the expected CSV. For

## Quick smoke test

```bash
Rscript run_tables_1_6.R --mode=quick
```

Quick mode uses only 20 Monte Carlo replications per scenario. It tests execution
and output structure.
Verify that two independent same-seed quick runs are identical with:

```bash
Rscript test_tables_1_6_determinism.R
```

Use `--tables=5,6` to restrict this regression test to selected tables.

## Full manuscript run

```bash
Rscript run_tables_1_6.R --mode=full
```

Full mode uses 2,000 Monte Carlo replications per scenario and the table-specific
seeds shown below. Tables 1--2 use 200 Gauss--Legendre nodes; Tables 3--4 use 80.
These values are set explicitly by the runner. Full mode also reruns Tables 1--2
with a finer quadrature rule on the identical Monte Carlo samples and archives
the cellwise differences. To run a subset, for example Tables 5--6, use:

```bash
Rscript run_tables_1_6.R --mode=full --tables=5,6
```

Existing raw CSVs, manuscript CSVs, logs, and quadrature checks are not replaced
unless `--overwrite` is supplied. A custom output root can be selected with
`--output-dir=PATH`. The runner keeps quick and full results in separate `quick/`
and `full/` subdirectories. Subset runs merge their rows into the cumulative
manifest rather than deleting metadata from earlier table runs.

## Table-to-file map

| Table | Distribution and comparison | Script | Seed | Output CSV |
|---:|---|---|---:|---|
| 1 | Kernel divergences, exponential | `Sim_RelativeMSE_comparison_Exponential_dist.R` | 2024 | `MSE_and_Relative_MSE_Section_5_1.csv` |
| 2 | Kernel divergences, Weibull | `Sim_RelativeMSE_comparison_Weibull_dist.R` | 2024 | `MSE_and_Relative_MSE_Section_5_1_Weibull.csv` |
| 3 | Three estimators, exponential | `Sim_MSE_Estimators_Exponential_dist.R` | 2026 | `MSE_and_Relative_MSE_results.csv` |
| 4 | Three estimators, Weibull | `Sim_MSE_Estimators_Weibull_dist.R` | 2024 | `MSE_and_Relative_MSE_results_Weibull.csv` |
| 5 | Three confidence intervals, exponential | `Sim_Confidence_Intervals_Exponential_dist.R` | 2026 | `Coverage_and_Average_Length_Section_5_3_Exponential_Efficient_Corrected.csv` |
| 6 | Three confidence intervals, Weibull | `Sim_Confidence_Intervals_Weibull_dist.R` | 2026 | `Coverage_and_Average_Length_Section_5_3_Weibull.csv` |

Tables 5--6 contain JEL, U-statistic normal, and empirical-estimator normal
intervals.

## Output layout and audit files

For a full run, results are written under:

```text
results/tables_1_6/full/
  table_1/
  table_2/
  table_3/
  table_4/
  table_5/
  table_6/
  run_manifest.csv
  sessionInfo.txt
```

Each table directory contains the unrounded raw CSV, a rounded
`Table_N_manuscript.csv`, and `run.log`. 
Full Table 1--2 directories additionally contain the same-sample quadrature-check
CSV and its checksum is recorded in the manifest.
`sessionInfo.txt` records the R and package environment. Validate a
completed full output tree with:

```bash
Rscript validate_tables_1_6.R \
  --output-dir=results/tables_1_6/full
```

Run the full workflow on a clean R installation. Commit the full CSVs, manifest, logs, and `sessionInfo.txt`
with the revision archive.
