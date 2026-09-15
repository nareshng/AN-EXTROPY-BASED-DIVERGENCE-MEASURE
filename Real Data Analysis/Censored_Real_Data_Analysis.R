#!/usr/bin/env Rscript
## =============================================================================
## Section 6 -- Table 11 and Figures 1-3: the Kaplan-Meier estimator D_tau^KM on
## three censored data sets.
##
##   Veteran lung cancer  standard vs test treatment   (survival::veteran)
##   Lung cancer          male vs female               (survival::lung)
##   GBSG2 breast cancer  no vs yes hormonal therapy   (TH.data::GBSG2 or survival::gbsg)
##
## For every comparison the script reports tau, the number still at risk at tau,
## D_tau^KM (eq. 4.5), the plug-in bias b_tau, the Greenwood standard error of
## Section 4.3, the normal interval (4.15), a bias-corrected normal interval, a
## nonparametric bootstrap of the (T, delta) pairs (SD, percentile and basic
## intervals, tau held fixed) and the log-rank p-value.  The primary table uses
## tau = smaller 80th percentile of the observed event times in the two groups,
## as stated in Section 6; a sensitivity table repeats the analysis for the larger
## 80th percentile and for fixed calendar horizons.
##
## Usage (run from this folder, or give the path to the script):
##   Rscript Censored_Real_Data_Analysis.R
##   Rscript Censored_Real_Data_Analysis.R --bootstrap-reps=1000 --output-dir=out
##   Rscript Censored_Real_Data_Analysis.R --self-test
##
## Options
##   --tau-rule=min|max        rule for the primary table (default min)
##   --tau-quantile=Q          quantile of the observed event times (default 0.80)
##   --fixed-horizons=A,B,C    horizons for the sensitivity table, one per data set
##                             in the order veteran, lung, GBSG2 (default 180,365,1825)
##   --bootstrap-reps=N        bootstrap resamples per comparison (default 5000)
##   --seed=N                  base seed (default 2026)
##   --alpha=A                 1 - confidence level (default 0.05)
##   --gbsg2-source=auto|TH.data|survival   GBSG2 provider (default auto)
##   --no-sensitivity          primary table only
##   --figure-style=paper|base  figure renderer. "paper" (default) reproduces the
##                             figures of the submitted manuscript with
##                             survminer::ggsurvplot and ggplot2 (three PDFs per
##                             data set plus the combined faceted figure, in
##                             <output-dir>/figures/); "base" uses base graphics
##                             and needs no extra package.
##   --cumulative-geom=step|line  geom for the cumulative-divergence panels
##                             (default step, as in the submitted figures; the
##                             integral is in fact piecewise linear, so "line"
##                             is the exact rendering)
##   --no-figures              skip the PDF figures
##   --output-dir=PATH         output directory (default real_data_outputs)
##   --self-test               run the deterministic checks and exit
##   --help
##
## Requires: survival.  No other package is needed; TH.data is optional.
## =============================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  ## Rscript encodes spaces in the script path as "~+~".
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE),
                        winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(script_dir, "real_data_functions.R"))
source(rd_find_km_functions(script_dir))
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required.", call. = FALSE)
}

## ---- command line ------------------------------------------------------------
CFG <- list(tau_rule = "min", tau_quantile = 0.80, fixed_horizons = c(180, 365, 1825),
            bootstrap_reps = 5000L, seed = 2026L, alpha = 0.05, gbsg2_source = "auto",
            sensitivity = TRUE, figures = TRUE, figure_style = "paper",
            cumulative_geom = "step", output_dir = "real_data_outputs",
            self_test = FALSE)
usage <- function() {
  cat(paste(readLines(file.path(script_dir, "Censored_Real_Data_Analysis.R"), n = 40)[-1L],
            collapse = "\n"), "\n")
}
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (identical(arg, "--help")) { usage(); quit(save = "no", status = 0L) }
  else if (identical(arg, "--self-test")) CFG$self_test <- TRUE
  else if (identical(arg, "--no-sensitivity")) CFG$sensitivity <- FALSE
  else if (identical(arg, "--no-figures")) CFG$figures <- FALSE
  else if (startsWith(arg, "--tau-rule=")) CFG$tau_rule <- sub("^--tau-rule=", "", arg)
  else if (startsWith(arg, "--tau-quantile=")) CFG$tau_quantile <- as.numeric(sub("^--tau-quantile=", "", arg))
  else if (startsWith(arg, "--fixed-horizons=")) CFG$fixed_horizons <- as.numeric(strsplit(sub("^--fixed-horizons=", "", arg), ",")[[1L]])
  else if (startsWith(arg, "--bootstrap-reps=")) CFG$bootstrap_reps <- as.integer(sub("^--bootstrap-reps=", "", arg))
  else if (startsWith(arg, "--seed=")) CFG$seed <- as.integer(sub("^--seed=", "", arg))
  else if (startsWith(arg, "--alpha=")) CFG$alpha <- as.numeric(sub("^--alpha=", "", arg))
  else if (startsWith(arg, "--gbsg2-source=")) CFG$gbsg2_source <- sub("^--gbsg2-source=", "", arg)
  else if (startsWith(arg, "--figure-style=")) CFG$figure_style <- sub("^--figure-style=", "", arg)
  else if (startsWith(arg, "--cumulative-geom=")) CFG$cumulative_geom <- sub("^--cumulative-geom=", "", arg)
  else if (startsWith(arg, "--output-dir=")) CFG$output_dir <- sub("^--output-dir=", "", arg)
  else stop("Unknown option: ", arg, call. = FALSE)
}
if (!CFG$tau_rule %in% c("min", "max")) stop("--tau-rule must be min or max.", call. = FALSE)
if (!CFG$figure_style %in% c("paper", "base")) {
  stop("--figure-style must be paper or base.", call. = FALSE)
}
if (!CFG$cumulative_geom %in% c("step", "line")) {
  stop("--cumulative-geom must be step or line.", call. = FALSE)
}
if (!is.finite(CFG$tau_quantile) || CFG$tau_quantile <= 0 || CFG$tau_quantile >= 1) {
  stop("--tau-quantile must lie strictly between 0 and 1.", call. = FALSE)
}
if (!is.finite(CFG$bootstrap_reps) || CFG$bootstrap_reps < 19L) {
  stop("--bootstrap-reps must be at least 19.", call. = FALSE)
}
if (length(CFG$fixed_horizons) != 3L || any(!is.finite(CFG$fixed_horizons)) ||
    any(CFG$fixed_horizons <= 0)) {
  stop("--fixed-horizons needs three positive values.", call. = FALSE)
}

if (CFG$self_test) {
  ok <- rd_self_test()
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

## The figure renderer is only needed for an actual run, never for --self-test.
if (CFG$figures && CFG$figure_style == "paper") {
  source(file.path(script_dir, "figures_paper_style.R"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("--figure-style=paper needs ggplot2 (and survminer for the ",
         "Kaplan-Meier panel). Install them, or use --figure-style=base.",
         call. = FALSE)
  }
}

## ---- data --------------------------------------------------------------------
veteran <- rd_load_data("veteran", "survival")   # status: 1 = death, 0 = censored
lung <- rd_load_data("lung", "survival")         # status: 2 = death, 1 = censored
gb <- rd_gbsg2(CFG$gbsg2_source)
lung_event <- as.integer(lung$status == 2L)
stopifnot(all(veteran$status %in% 0:1), all(lung$status %in% 1:2), all(gb$event %in% 0:1))

cases <- list(
  list(name = "Veteran lung cancer", comparison = "Standard treatment vs Test treatment",
       t1 = veteran$time[veteran$trt == 1], s1 = veteran$status[veteran$trt == 1],
       t2 = veteran$time[veteran$trt == 2], s2 = veteran$status[veteran$trt == 2],
       labels = c("Standard treatment", "Test treatment"),
       source = paste0("survival::veteran (survival ", utils::packageVersion("survival"), ")"),
       horizon = CFG$fixed_horizons[1], seed = CFG$seed + 1L),
  list(name = "Lung cancer", comparison = "Male vs Female",
       t1 = lung$time[lung$sex == 1], s1 = lung_event[lung$sex == 1],
       t2 = lung$time[lung$sex == 2], s2 = lung_event[lung$sex == 2],
       labels = c("Male", "Female"),
       source = paste0("survival::lung (survival ", utils::packageVersion("survival"), ")"),
       horizon = CFG$fixed_horizons[2], seed = CFG$seed + 2L),
  list(name = "GBSG2 breast cancer", comparison = "No hormonal therapy vs Hormonal therapy",
       t1 = gb$time[gb$hormone == 0], s1 = gb$event[gb$hormone == 0],
       t2 = gb$time[gb$hormone == 1], s2 = gb$event[gb$hormone == 1],
       labels = c("No hormonal therapy", "Hormonal therapy"),
       source = gb$source, horizon = CFG$fixed_horizons[3], seed = CFG$seed + 3L)
)

out_dir <- CFG$output_dir
if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create the output directory: ", out_dir, call. = FALSE)
}
slug <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
fig_dir <- file.path(out_dir, "figures")
if (CFG$figures && CFG$figure_style == "paper") {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}
all_step_dfs <- list()

## ---- primary table (Table 11) -------------------------------------------------
primary <- list()
for (k in seq_along(cases)) {
  cs <- cases[[k]]
  tau <- rd_choose_tau(cs$t1, cs$s1, cs$t2, cs$s2, CFG$tau_rule, CFG$tau_quantile)
  label <- sprintf("%s of the two %g quantiles of the observed event times",
                   CFG$tau_rule, CFG$tau_quantile)
  row <- rd_analyse(cs$name, cs$comparison, cs$t1, cs$s1, cs$t2, cs$s2, tau, label,
                    R_boot = CFG$bootstrap_reps, alpha = CFG$alpha, seed = cs$seed)
  row$Data_source <- cs$source
  steps <- rd_cumulative(cs$t1, cs$s1, cs$t2, cs$s2, tau)
  utils::write.csv(steps, file.path(out_dir, paste0(slug(cs$name), "_cumulative_divergence.csv")),
                   row.names = FALSE)
  if (CFG$figures) {
    if (CFG$figure_style == "base") {
      rd_figure(row, cs$t1, cs$s1, cs$t2, cs$s2, steps, cs$labels,
                file.path(out_dir, paste0(slug(cs$name), "_figure.pdf")))
    } else {
      sdf <- fig_step_df(steps, cs$name, cs$comparison)
      all_step_dfs[[length(all_step_dfs) + 1L]] <- sdf
      fig_dataset_figures(row, cs$t1, cs$s1, cs$t2, cs$s2, sdf, cs$labels, fig_dir,
                          cumulative_geom = CFG$cumulative_geom)
    }
  }
  primary[[k]] <- row
  cat(sprintf("%-22s tau = %8.1f  D = %.4f  SE = %.4f  at risk = (%d, %d)\n",
              cs$name, tau, row$D_tau_KM, row$SE_Greenwood, row$At_risk_tau_1, row$At_risk_tau_2))
}
primary <- do.call(rbind, primary)
utils::write.csv(primary, file.path(out_dir, "Table11_real_data_full.csv"), row.names = FALSE)

if (length(all_step_dfs)) {
  stacked <- do.call(rbind, all_step_dfs)
  utils::write.csv(stacked, file.path(out_dir, "three_cases_stepwise_divergence_values.csv"),
                   row.names = FALSE)
  fig_combined_figure(stacked, fig_dir, cumulative_geom = CFG$cumulative_geom)
}

printed <- data.frame(
  Dataset = primary$Dataset, Comparison = primary$Comparison,
  n1 = primary$n1, n2 = primary$n2,
  `Cens. 1 (%)` = round(primary$Censoring_1_pct, 1),
  `Cens. 2 (%)` = round(primary$Censoring_2_pct, 1),
  tau = round(primary$tau, 1),
  `At risk at tau` = paste0("(", primary$At_risk_tau_1, ", ", primary$At_risk_tau_2, ")"),
  `D_tau^KM` = round(primary$D_tau_KM, 4),
  `b_tau` = round(primary$bias_hat, 4),
  SE = round(primary$SE_Greenwood, 4),
  `95% normal CI` = sprintf("(%.3f, %.3f)", primary$Normal_lo, primary$Normal_hi),
  `95% bias-corrected CI` = sprintf("(%.3f, %.3f)", primary$BiasCorr_lo, primary$BiasCorr_hi),
  `95% percentile bootstrap CI` = sprintf("(%.3f, %.3f)", primary$Percentile_lo, primary$Percentile_hi),
  `Log-rank p` = signif(primary$Logrank_p, 3),
  check.names = FALSE, stringsAsFactors = FALSE
)
utils::write.csv(printed, file.path(out_dir, "Table11_formatted.csv"), row.names = FALSE)
cat("\n"); print(printed, row.names = FALSE)

## ---- sensitivity in tau -------------------------------------------------------
if (CFG$sensitivity) {
  rows <- list()
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    taus <- list(
      list(tau = rd_choose_tau(cs$t1, cs$s1, cs$t2, cs$s2, "min", CFG$tau_quantile),
           label = sprintf("min of the two %g quantiles", CFG$tau_quantile)),
      list(tau = rd_choose_tau(cs$t1, cs$s1, cs$t2, cs$s2, "max", CFG$tau_quantile),
           label = sprintf("max of the two %g quantiles", CFG$tau_quantile)),
      list(tau = cs$horizon, label = sprintf("fixed horizon (%g time units)", cs$horizon))
    )
    for (tt in taus) {
      rows[[length(rows) + 1L]] <- rd_analyse(
        cs$name, cs$comparison, cs$t1, cs$s1, cs$t2, cs$s2, tt$tau, tt$label,
        R_boot = CFG$bootstrap_reps, alpha = CFG$alpha, seed = cs$seed)
    }
  }
  sens <- do.call(rbind, rows)
  utils::write.csv(sens, file.path(out_dir, "Table11_tau_sensitivity.csv"), row.names = FALSE)
  cat("\nSensitivity of Table 11 to the truncation point:\n")
  print(sens[, c("Dataset", "tau_rule", "tau", "At_risk_tau_1", "At_risk_tau_2",
                 "D_tau_KM", "SE_Greenwood", "Normal_lo", "Normal_hi")],
        row.names = FALSE, digits = 4)
}

## ---- provenance ---------------------------------------------------------------
settings <- c(
  paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("tau rule:", CFG$tau_rule, "of the two", CFG$tau_quantile,
        "quantiles of the observed event times (stats::quantile type 7)"),
  "Bootstrap draws are shared across the tau rules of one data set (paired sensitivity)",
  paste("Figure renderer:", CFG$figure_style,
        if (CFG$figure_style == "paper") {
          paste0("(survminer ", if (requireNamespace("survminer", quietly = TRUE))
            as.character(utils::packageVersion("survminer")) else "not installed",
            "; cumulative geom: ", CFG$cumulative_geom, ")")
        } else ""),
  paste("Bootstrap resamples:", CFG$bootstrap_reps),
  paste("Base seed:", CFG$seed, "(dataset seeds: base + 1, 2, 3)"),
  paste("Confidence level:", 1 - CFG$alpha),
  paste("GBSG2 source:", gb$source),
  paste("R version:", R.version.string),
  paste("survival:", as.character(utils::packageVersion("survival"))),
  "RNG: Mersenne-Twister; Inversion; Rejection"
)
writeLines(settings, file.path(out_dir, "Section6_run_settings.txt"))
utils::capture.output(utils::sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("\nOutputs written to ", normalizePath(out_dir), "\n", sep = "")
