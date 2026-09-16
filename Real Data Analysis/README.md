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
Rscript Image_based_Real_Data_Analysis.R --self-test # 10 checks
```

The censored self-test recomputes D and its Greenwood variance on the veteran
data from `survival::survfit` output, and recomputes bootstrap draws one by one.
The image self-test verifies the plug-in integral against its V-statistic form,
the equation-(2.5) correction, behaviour under ties, symmetry, scale
equivariance, grayscale conversion and deterministic resizing.

## Table 11

```bash
Rscript Censored_Real_Data_Analysis.R                       # defaults
Rscript Censored_Real_Data_Analysis.R --bootstrap-reps=1000 --output-dir=out
```

Outputs in `real_data_outputs/`: `Table11_formatted.csv` (the printed table),
`Table11_real_data_full.csv` (all quantities, including the bias estimate, the
bootstrap SD and the basic interval), `Table11_tau_sensitivity.csv`,
`<stem>_cumulative_divergence.csv`, the six figure PDFs above,
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

Each run writes exactly six PDFs, named for the manuscript, in the output
directory (no other figure files, no subfolder):

```text
Veteran_KM_curve.pdf        Veteran_cumulative_divergence.pdf
Lung_cancer_KM_curve.pdf    Lung_cancer_cumulative_divergence.pdf
GBSG2_KM_curve.pdf          GBSG2_cumulative_divergence.pdf
```

`--figure-style=paper` (default) draws them as in the submitted manuscript:
`survminer::ggsurvplot(conf.int = TRUE, pval = TRUE, ggtheme = theme_bw())` for
the Kaplan–Meier panel and `ggplot2::geom_step` with a dashed line at τ and the
τ / D̂ subtitle for the cumulative divergence, both 7.2 × 4.8 in. It needs
`ggplot2` and `survminer`; without `survminer` the Kaplan–Meier panel is drawn
with `ggplot2` alone and the script says so. `--figure-style=base` produces the
same two files per data set with base graphics and no extra package, so the file
names in the manuscript never change. `--no-figures` skips them.

`--cumulative-geom` controls the cumulative panels: `step` (default) matches the
submitted figures exactly; `line` is the exact rendering, since the integral of a
squared step function is continuous and piecewise linear between the pooled jump
times (the staircase misstates it between knots by about 4% of D̂ on the veteran
data; endpoints and all reported numbers are unaffected).

The numbers in the figures change from the submitted version, because τ now
follows the rule stated in Section 6 (the smaller 80th percentile) and the
standard error uses the corrected Section 4.3 scaling. For the lung comparison
the subtitle becomes τ = 421.2, D̂ = 12.6792 in place of τ = 501.2, D̂ = 15.5826.

## Tables 12–14

```bash
Rscript Image_based_Real_Data_Analysis.R
Rscript Image_based_Real_Data_Analysis.R --indices=1,2,3
Rscript Image_based_Real_Data_Analysis.R --resize=none       # sensitivity analysis
Rscript Image_based_Real_Data_Analysis.R --estimator=plugin  # diagnostic only
```

The default analysis now implements the preprocessing and estimator that must
be stated in the manuscript:

- every image is converted to grayscale and resized to **227 × 255 pixels
  (width × height)** by bilinear interpolation;
- the primary table value is the empirical estimator in equation (2.5),
  evaluated through the tie-safe identity
  `D_Emp = D_plugin - mean(x)/n1 - mean(y)/n2`;
- the nonnegative plug-in integral is retained only as an explicitly labelled
  diagnostic; and
- Figure 4 is written with rows `NT`, `BT`, `MT` and columns equal to the image
  indices, matching the manuscript ordering.

The equation-(2.5) estimator can be slightly negative in finite samples. The
script reports such values rather than truncating or silently replacing them by
the plug-in estimator.

Outputs in `image_outputs/` are:

- `Divergence_matrix_index<k>.csv`: the selected primary estimator; under the
  default settings this is equation (2.5);
- `Divergence_matrix_index<k>_plugin.csv`: the plug-in diagnostic;
- `Image_divergences_long.csv`: both labelled estimators at full precision;
- `Image_inventory.csv`: input MD5 checksum, original and analysed dimensions,
  resize status, pixel counts and intensity summaries;
- `Figure4_MRI_images.pdf`, `Section6_1_run_settings.txt` and `sessionInfo.txt`.

Do not copy the old values in Tables 12–14 after this preprocessing correction.
Run the script and update the manuscript from `Divergence_matrix_index<k>.csv`.
State the 227 × 255 target and bilinear interpolation explicitly in Section 6.1.
The run settings and inventory should be archived with those tables.

### Image provenance

For a reproducible paper, add `Images/image_manifest.csv` with one row for each
image. The required column is `file`; the recommended fields are
`analysis_group`, `analysis_index`, `source_dataset_class`, `source_record_id`,
`source_url` and `selection_rule`. For example, the header can be

```text
file,analysis_group,analysis_index,source_dataset_class,source_record_id,source_url,selection_rule
```

The script checks the filename/group/index mapping and appends all manifest
fields to `Image_inventory.csv`. If no manifest is supplied, it records file
checksums but warns that the dataset labels and source mapping cannot be
verified. Do not describe `BT` and `MT` as benign and malignant unless that
mapping is supported by the original dataset metadata.

### Interpretation limitation

This is a descriptive comparison of marginal pixel-intensity distributions.
Pixels from one image are spatially dependent, flattening discards tumour
location and texture, and three manually selected images per label are not a
representative sample. Consequently, these tables do not establish diagnostic
accuracy, class separation or a biological difference. In particular, no text
should claim that every normal–tumour divergence exceeds every tumour–tumour
divergence unless the updated tables actually show that ordering.
