# Section 6 — real-data analysis

Two independent analyses, each runnable with one command from this folder.

| Paper output | Script | Command |
|---|---|---|
| Table 11, Figures 1–3 | `Censored_Real_Data_Analysis.R` | `Rscript Censored_Real_Data_Analysis.R` |
| Tables 12–14, Figure 4 | `Image_based_Real_Data_Analysis.R` | `Rscript Image_based_Real_Data_Analysis.R` |

Shared code: `real_data_functions.R` (data loading, analysis, figures, checks) and
`image_functions.R` (image reading and the divergence estimator). The censored
analysis reads the estimator, the Greenwood variance and the bootstrap from
`../Censored Simulations/km_functions.R`, so Section 6 and the Section 5.4
simulations use one implementation.

## Software

R (>= 4.1) and `survival` for Table 11; `jpeg` for Tables 12–14 (or `png` with
`--extension=png`). `TH.data` is optional: `survival::gbsg` holds the same 686
GBSG2 patients and is used automatically when `TH.data` is absent. Both scripts
run in a few seconds; the default 5000 bootstrap resamples for Table 11 take
about a minute in total.

Check the installation before running anything:

```bash
Rscript Censored_Real_Data_Analysis.R --self-test   # 5 checks
Rscript Image_based_Real_Data_Analysis.R --self-test # 8 checks
```

The censored self-test recomputes D and its Greenwood variance on the veteran
data from `survival::survfit` output, and recomputes bootstrap draws one by one;
the image self-test verifies the estimator against its V-statistic form, its
behaviour under ties, symmetry, scale equivariance and the grayscale conversion.

## Table 11

```bash
Rscript Censored_Real_Data_Analysis.R                       # defaults
Rscript Censored_Real_Data_Analysis.R --bootstrap-reps=1000 --output-dir=out
```

Outputs in `real_data_outputs/`: `Table11_formatted.csv` (the printed table),
`Table11_real_data_full.csv` (all quantities, including the bias estimate, the
bootstrap SD and the basic interval), `Table11_tau_sensitivity.csv`,
`<dataset>_cumulative_divergence.csv`, `<dataset>_figure.pdf`,
`Section6_run_settings.txt` and `sessionInfo.txt`.

The truncation point is stated explicitly rather than implied: the primary table
uses the **smaller** of the two 80th percentiles of the observed event times
(`--tau-rule=min`, the rule stated in Section 6). Because the choice matters, the
script also reports the larger percentile and a fixed calendar horizon for every
data set; with 1000 resamples the three rules give

| Data set | τ rule | τ | at risk at τ | D̂ | SE |
|---|---|---|---|---|---|
| Veteran | min / max / 180 days | 168.0 / 172.8 / 180 | (14,14) / (14,14) / (13,14) | 1.700 / 1.701 / 1.702 | 2.023 / 2.021 / 2.013 |
| Lung | min / max / 365 days | 421.2 / 501.2 / 365 | (29,26) / (20,21) / (35,30) | 12.679 / 15.583 / 10.442 | 6.829 / 8.499 / 5.585 |
| GBSG2 | min / max / 1825 days | 1165.2 / 1356.4 / 1825 | (185,126) / (141,109) / (63,60) | 5.365 / 8.023 / 15.536 | 3.990 / 5.611 / 10.088 |

The quantile rule is `stats::quantile` type 7 (R's default). The type matters —
for the veteran data the 80th percentile of the observed event times is 165, 168,
174 or 177 under types 4, 7, 8 and 1 — so it is fixed in the code and recorded in
`Section6_run_settings.txt`. Within a data set the three tau rules share the same
bootstrap resamples, so the sensitivity comparison is paired.

Section 6 of the paper must name the rule it reports. Standard errors follow
Section 4.3 with the factors n₁ and n₂, i.e. Var(D̂) = 4σ̂²_F/n₁ + 4σ̂²_G/n₂; the
self-test compares them with the direct double sum.

### Figures

`--figure-style=paper` is the default and reproduces the figures of the submitted
manuscript with the same code as the original repository: `survminer::ggsurvplot`
for the Kaplan–Meier panel (`conf.int = TRUE`, `pval = TRUE`, `theme_bw()`),
`ggplot2::geom_step` for the squared survival difference and for the cumulative
divergence, each with the dashed line at τ and the τ / D̂ subtitle, plus the
combined `facet_wrap(~ Dataset, scales = "free")` figure. File names, page sizes
(7.2 × 4.8 in, and 10 × 5.8 in for the combined figure) and axis expressions are
those of the original script, and everything is written to
`<output-dir>/figures/`:

```text
<Dataset>_<group1>_vs_<group2>_KM_curve.pdf
<Dataset>_<group1>_vs_<group2>_squared_difference.pdf
<Dataset>_<group1>_vs_<group2>_cumulative_divergence.pdf
three_cases_combined_cumulative_divergence.pdf
three_cases_stepwise_divergence_values.csv
```

This style needs `ggplot2` and `survminer`; without `survminer` the Kaplan–Meier
panel is drawn with `ggplot2` alone and the script says so. `--figure-style=base`
draws everything with base graphics and needs no extra package.

`--cumulative-geom` controls the cumulative panels: `step` (default) matches the
submitted figures exactly, `line` is the exact rendering, since the integral of a
squared step function is continuous and piecewise linear between the pooled jump
times (the staircase misstates it between knots by about 4% of D̂ on the veteran
data; endpoints and all reported numbers are unaffected).

The numbers in the figures change from the submitted version, because τ now
follows the rule stated in Section 6 (the smaller 80th percentile) and the
standard error uses the corrected Section 4.3 scaling. For the lung comparison
the subtitle becomes τ = 421.2, D̂ = 12.6792 in place of τ = 501.2, D̂ = 15.5826.

## Tables 12–14

```bash
Rscript Image_based_Real_Data_Analysis.R                    # every complete triple in Images/
Rscript Image_based_Real_Data_Analysis.R --indices=1,2,4
```

Outputs in `image_outputs/`: `Divergence_matrix_index<k>.csv` (one 3 × 3 matrix
per index), `Image_divergences_long.csv` (also carrying the equation-(2.5)
value, which subtracts the diagonal terms and can therefore be slightly negative
for near-identical images), `Image_inventory.csv` (size, pixel count and number of distinct
intensities per file), `Figure4_MRI_images.pdf` and the run settings.

Two points to reconcile with the manuscript. The repository ships images with
indices 1, 2, 4 and 5, while the paper labels the three panels of Tables 12–14 as
NT1/NT2/NT3; the third table corresponds to the files with index 4, so either
rename the files or relabel the tables. And the estimator is evaluated exactly on
the pooled distinct intensities: an 8-bit image has at most 256 distinct values,
and a rank formula with `ties.method = "max"` does not evaluate equation (2.5)
under that many ties. The exact values are

| Index | NT–BT | NT–MT | BT–MT |
|---|---|---|---|
| 1 | 0.0358 | 0.0485 | 0.0118 |
| 2 | 0.0239 | 0.0093 | 0.0460 |
| 4 | 0.0229 | 0.0317 | 0.0232 |
| 5 | 0.0162 | 0.0095 | 0.0083 |

against 0.041 / 0.052 / 0.017, 0.027 / 0.014 / 0.050 and 0.026 / 0.039 / 0.031 in
the submitted Tables 12–14. The ordering of the pairs is unchanged, so the
conclusions of Section 6.1 stand, but the printed numbers need updating.
