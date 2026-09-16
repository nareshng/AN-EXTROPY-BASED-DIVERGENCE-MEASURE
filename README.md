# Nonparametric Inference for an Extropy-Based Divergence Measure

R code for reproducing every table and figure in

> Garg, N., Dewan, I. and Kattumannil, S. K. *Nonparametric inference for an extropy-based divergence measure.* (Manuscript under revision, *Biometrical Journal*.)

## Abstract

Survival extropy, which quantifies the uncertainty associated with the remaining lifetime distribution, provides an information-theoretic perspective on survival behaviour. We consider a divergence measure based on survival extropy and derive its nonparametric estimators based on U-statistics, empirical distribution functions and kernel density estimation. We construct confidence intervals for the divergence measure using the jackknife empirical likelihood (JEL) method and the normal approximation method with a jackknife pseudo-value-based variance estimator. The divergence measure is also extended to randomly right-censored data, covering both point estimation and confidence interval construction. A Monte Carlo simulation study compares the proposed estimators with estimators of other divergence measures and evaluates the finite-sample performance of the proposed estimators and intervals for uncensored and censored data. The measure is illustrated on real survival data sets and on MRI image data.

**Keywords:** Extropy; Jackknife empirical likelihood; Measure of divergence; U-statistics.

---

## Repository structure

```text
.
├── README.md
├── Uncensored Simulations/            Section 5.1–5.3 (Tables 1–6)
│   ├── run_tables_1_6.R               one-command runner for Tables 1–6
│   ├── validate_tables_1_6.R          source and output checks
│   ├── test_tables_1_6_determinism.R  same-seed determinism test
│   ├── Sim_RelativeMSE_comparison_Exponential_dist.R   Table 1
│   ├── Sim_RelativeMSE_comparison_Weibull_dist.R       Table 2
│   ├── Sim_MSE_Estimators_Exponential_dist.R           Table 3
│   ├── Sim_MSE_Estimators_Weibull_dist.R               Table 4
│   ├── Sim_Confidence_Intervals_Exponential_dist.R     Table 5
│   ├── Sim_Confidence_Intervals_Weibull_dist.R         Table 6
│   └── README.md
├── Censored Simulations/              Section 5.4 (Tables 7–10)
│   ├── run_tables_7_10.R              one-command runner for Tables 7–10
│   ├── Point_estimation_right_censoring.R              Tables 7–8
│   ├── Confidence_Intervals_right_censoring.R          Tables 9–10
│   ├── Aggregate_right_censoring_results.R
│   ├── km_functions.R                 Kaplan–Meier estimator, Greenwood variance, bootstrap
│   ├── censored_simulation_helpers.R
│   ├── diagnose_tables.R
│   ├── atomic_io_functions.R
│   └── README.md
└── Real Data Analysis/                Section 6 (Table 11–14, Figures 1–4)
    ├── Censored_Real_Data_Analysis.R  Table 11, Figures 1–3
    ├── Image_based_Real_Data_Analysis.R               Tables 12–14, Figure 4
    ├── real_data_functions.R
    ├── image_functions.R
    ├── figures_paper_style.R
    ├── Images/                        NT1–3, BT1–3, MT1–3 (.jpg)
    └── README.md
```

## Requirements

| Software | Needed for | Version used |
|---|---|---|
| R (≥ 4.1) | everything | 4.3.3 and 4.6.1 |
| `emplik` **1.3-3** | Tables 5–6 (JEL) | 1.3-3 (exact version is checked) |
| `survival` | Table 11, Figures 1–3 | 3.5-8 |
| `TH.data` | Table 11 (GBSG2 data) | 1.1-2 |
| `jpeg` | Tables 12–14, Figure 4 | 0.1-10 |
| `ggplot2`, `survminer` | paper-style Figures 1–3 | 3.4.4, 0.4.9 |

Tables 1–4 and 7–10 use base R only.

```r
install.packages(c("survival", "TH.data", "jpeg", "ggplot2", "survminer", "remotes"))
remotes::install_version("emplik", version = "1.3-3",
                         repos = "https://cloud.r-project.org", upgrade = "never")
```

`TH.data` is strongly recommended: without it the script falls back to `survival::gbsg` (same 686 patients, different row order), which changes the GBSG2 bootstrap interval in Table 11 in the second decimal.

## Quick check (a few minutes)

Run from inside each folder. These use very few replications and only test that the code runs; their numbers are **not** the paper's results.

```bash
cd "Uncensored Simulations"
Rscript validate_tables_1_6.R
Rscript run_tables_1_6.R --mode=quick

cd "../Censored Simulations"
Rscript run_tables_7_10.R --mode=quick --cores=2 --tidy

cd "../Real Data Analysis"
Rscript Censored_Real_Data_Analysis.R --self-test
Rscript Image_based_Real_Data_Analysis.R --self-test
```

## Reproducing the paper

### Tables 1–6 (Sections 5.1–5.3)

```bash
cd "Uncensored Simulations"
Rscript run_tables_1_6.R --mode=full                 # all six tables
Rscript run_tables_1_6.R --mode=full --tables=5,6    # a subset
```

Each table runs in a clean R process with 2000 Monte Carlo replications. Outputs are written to `results/tables_1_6/full/table_N/`, one folder per table: the raw CSV, `Table_N_manuscript.csv` (rounded as printed), the run log, `run_manifest.csv` (seeds, versions, MD5 checksums) and `sessionInfo.txt`. Existing results are not overwritten unless `--overwrite` is given.

| Table | Content | Script | Seed | Quadrature nodes |
|---|---|---|---|---|
| 1 | Relative MSE of kernel estimators of D, D_CC, KL (exponential) | `Sim_RelativeMSE_comparison_Exponential_dist.R` | 2024 | 200 |
| 2 | Same, Weibull | `Sim_RelativeMSE_comparison_Weibull_dist.R` | 2024 | 200 |
| 3 | MSE of kernel, empirical and U-statistic estimators (exponential) | `Sim_MSE_Estimators_Exponential_dist.R` | 2026 | 80 |
| 4 | Same, Weibull | `Sim_MSE_Estimators_Weibull_dist.R` | 2024 | 80 |
| 5 | CP and AL of JEL, NA and empirical intervals (exponential) | `Sim_Confidence_Intervals_Exponential_dist.R` | 2026 | – |
| 6 | Same, Weibull | `Sim_Confidence_Intervals_Weibull_dist.R` | 2026 | – |

Notes:

* The JEL interval uses the jackknife pseudo-values of Jing, Yuan and Zhou (2009, *JASA* 104, 1224–1232) for two-sample U-statistics, with the leave-one-out statistic evaluated with the full-sample normalising constants and centred at E[V_k] from their eq. (14) (manuscript eq. 3.8).
* The KL kernel estimator bounds the density estimates below by 10⁻¹², as stated in Section 5.1.
* Parameterisation: exponential λ is a **rate** (`rexp(rate = λ)`); Weibull (k, λ) is (**shape**, **scale**) (`rweibull(shape = k, scale = λ)`).
* On Windows, run the individual scripts with the environment variables `EXTROPY_MC_REPS`, `EXTROPY_SEED` and `EXTROPY_OUTPUT_DIR` set, because `system2(env = …)` in the runner is not supported for `Rscript` on Windows. Example (PowerShell):
  `$env:EXTROPY_MC_REPS=2000; $env:EXTROPY_SEED=2026; Rscript Sim_Confidence_Intervals_Weibull_dist.R`

### Tables 7–10 (Section 5.4, right censoring)

```bash
cd "Censored Simulations"
Rscript run_tables_7_10.R --mode=full --cores=4 --tidy
```

Full mode uses 2000 Monte Carlo replications per cell, 1000 bootstrap resamples for Tables 9–10, point-estimation seed 5401 and confidence-interval seed 5402. Each cell has its own seeds, so results are identical for any `--cores` value and for runs split across machines. With `--tidy`, the four tables are written to `paper_tables/`:

| Table | Output file |
|---|---|
| 7 | `Table7_RelMSE_Exponential.csv` (column `D_KM`) |
| 8 | `Table8_RelMSE_Weibull.csv` (column `D_KM`) |
| 9 | `Table9_CI_Exponential.csv` |
| 10 | `Table10_CI_Weibull.csv` |

The CSVs contain more columns than the printed tables:

* `DCC_KM` is a Cox–Czanner comparator that is not reported in the paper.
* `Valid_%` is the percentage of replications in which an interval could be computed.
* `*_MCSE` columns give Monte Carlo standard errors.

Coverage in Tables 9–10 is computed over all 2000 replications. A replication in which no interval is available counts as non-coverage; this happens when a group has no subject at risk at τ or has a terminal event before τ.

Long runs can be split and resumed; see `Censored Simulations/README.md` (`--task-ids`, `--aggregate-only`).

### Table 11 and Figures 1–3 (Section 6)

```bash
cd "Real Data Analysis"
Rscript Censored_Real_Data_Analysis.R
```

* Defaults: τ is the smaller of the two 80th percentiles of the observed event times, 5000 bootstrap resamples, seed 2026.
* Outputs go to `real_data_outputs/`:
  * `Table11_formatted.csv`: the printed table, except that the normal limits are not truncated at 0 in the CSV.
  * `Table11_real_data_full.csv`
  * `Table11_tau_sensitivity.csv`
  * Figure PDFs
  * `sessionInfo.txt`

Data sources, all loaded from R packages:

* Veteran lung cancer: `survival::veteran` (Kalbfleisch and Prentice, 1980).
* Lung cancer (NCCTG): `survival::lung` (Loprinzi et al., 1994).
* GBSG2: `TH.data::GBSG2` (Schumacher et al., 1994).

### Tables 12–14 and Figure 4 (Section 6.1)

```bash
cd "Real Data Analysis"
Rscript Image_based_Real_Data_Analysis.R
```

* The script reads the nine images in `Images/`: no tumour (NT), benign tumour (BT) and malignant tumour (MT). They are taken from the Kaggle brain-tumour MRI dataset cited in the paper.
* It converts each image to grayscale in [0, 1] and computes the tie-safe empirical estimate for each pair.
* Image indices 1, 2 and 3 correspond to Tables 12, 13 and 14.
* Outputs are written to `image_outputs/`.

## Approximate run times

Measured on a single core (R 4.3.3, Linux). Other machines will differ.

| Step | Time |
|---|---|
| Tables 1 and 2 (each, including the quadrature check) | about 5–10 min |
| Tables 5 and 6 (each) | about 15–25 min |
| Tables 7–8 | about 5 min |
| Tables 9–10 | about 2–3 CPU-hours (use `--cores`) |
| Table 11 and Figures 1–3 | under 1 min |
| Tables 12–14 and Figure 4 | a few seconds |

## Reproducibility checks

The following outputs were regenerated from this repository and compared with the printed tables value by value:

* Tables 1–8
* Table 11
* Tables 12–14
* A subset of cells of Tables 9–10

To verify determinism of Tables 1–6:

```bash
cd "Uncensored Simulations"
Rscript test_tables_1_6_determinism.R --tables=3,4
```

Random-number settings are fixed explicitly (`RNGkind("Mersenne-Twister", "Inversion", "Rejection")`), so results do not depend on the R version's default sampler.

## Citation

```bibtex
@article{GargDewanKattumannil,
  author  = {Garg, Naresh and Dewan, Isha and Kattumannil, Sudheesh Kumar},
  title   = {Nonparametric Inference for an Extropy-Based Divergence Measure},
  journal = {Biometrical Journal},
  note    = {Under revision}
}
```

## Contact

Naresh Garg — garg.naresh22@gmail.com
