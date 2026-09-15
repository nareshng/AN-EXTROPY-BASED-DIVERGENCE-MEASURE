## =============================================================================
## Section 6 -- shared functions for the censored real-data analysis.
##
## The estimator, its Greenwood variance, the plug-in bias and the bootstrap all
## come from "Censored Simulations/km_functions.R", so the real-data analysis and
## the simulations of Section 5.4 use one implementation.
##
## Requires: survival (data sets and the log-rank test).  TH.data is used for
## GBSG2 when it is installed; survival::gbsg holds the same 686 patients.
## =============================================================================

## Input: the directory that contains this file.
## Output: the path of km_functions.R, or an informative error.
rd_find_km_functions <- function(script_dir) {
  candidates <- c(
    file.path(dirname(script_dir), "Censored Simulations", "km_functions.R"),
    file.path(script_dir, "km_functions.R")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Cannot find km_functions.R in '../Censored Simulations/' or beside this script.",
         call. = FALSE)
  }
  hit[1L]
}

## Input: a data set name and the package that provides it.
## Output: the data set. Lazy-loaded data (survival >= 3.5) is also handled.
rd_load_data <- function(object, package) {
  env <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(list = object, package = package, envir = env))
  if (exists(object, envir = env, inherits = FALSE)) {
    return(get(object, envir = env, inherits = FALSE))
  }
  value <- tryCatch(getExportedValue(package, object), error = function(e) NULL)
  if (is.null(value)) {
    stop("Data set '", object, "' was not found in package '", package, "'.", call. = FALSE)
  }
  value
}

## Input: the requested GBSG2 source ("auto", "TH.data" or "survival").
## Output: list(time, event, hormone, source); the two sources are the same study.
rd_gbsg2 <- function(source = c("auto", "TH.data", "survival")) {
  source <- match.arg(source)
  use_th <- switch(source,
    auto = requireNamespace("TH.data", quietly = TRUE),
    TH.data = TRUE,
    survival = FALSE
  )
  if (use_th) {
    if (!requireNamespace("TH.data", quietly = TRUE)) {
      stop("Package 'TH.data' was requested but is not installed.", call. = FALSE)
    }
    gb <- rd_load_data("GBSG2", "TH.data")
    out <- list(
      time = as.numeric(gb$time),
      event = as.integer(as.character(gb$cens)),
      hormone = as.integer(tolower(as.character(gb$horTh)) == "yes"),
      source = paste0("TH.data::GBSG2 (TH.data ", utils::packageVersion("TH.data"), ")")
    )
  } else {
    gb <- rd_load_data("gbsg", "survival")
    out <- list(
      time = as.numeric(gb$rfstime),
      event = as.integer(gb$status),
      hormone = as.integer(gb$hormon),
      source = paste0("survival::gbsg (survival ", utils::packageVersion("survival"), ")")
    )
  }
  ## Both sources are the German Breast Cancer Study Group 2 trial: 686 patients,
  ## 299 recurrence-or-death events, 246 on hormonal therapy.
  stopifnot(length(out$time) == 686L, sum(out$event) == 299L, sum(out$hormone) == 246L)
  out
}

## Input: two censored samples, a rule and a quantile.
## Output: the truncation point tau on the observed event-time scale.
rd_choose_tau <- function(t1, s1, t2, s2, rule = c("min", "max"), q = 0.80) {
  rule <- match.arg(rule)
  e1 <- t1[s1 == 1L]; e2 <- t2[s2 == 1L]
  if (!length(e1) || !length(e2)) stop("Both groups need at least one event.", call. = FALSE)
  ## type = 7 is R's default; the choice matters (for the veteran data the
  ## 80th percentile rule gives tau = 165, 168, 174 or 177 for types 4, 7, 8, 1),
  ## so it is fixed here and recorded in the run settings.
  q1 <- stats::quantile(e1, q, names = FALSE, type = 7)
  q2 <- stats::quantile(e2, q, names = FALSE, type = 7)
  if (rule == "min") min(q1, q2) else max(q1, q2)
}

## Input: two censored samples, a fixed tau, bootstrap settings.
## Output: one row of estimates, standard errors, intervals and diagnostics.
rd_analyse <- function(name, comparison, t1, s1, t2, s2, tau, tau_label,
                       R_boot = 5000L, alpha = 0.05, seed = 2026L) {
  fit <- D_km(t1, s1, t2, s2, tau)                       # eq. (4.5) + Section 4.3
  z <- stats::qnorm(1 - alpha / 2)
  set.seed(seed)
  bt <- boot_engine(t1, s1, t2, s2, tau, fit$D, R = R_boot, alpha = alpha)
  group <- rep(1:2, c(length(t1), length(t2)))
  lr <- survival::survdiff(survival::Surv(c(t1, t2), c(s1, s2)) ~ group)
  data.frame(
    Dataset = name, Comparison = comparison, tau_rule = tau_label, tau = tau,
    n1 = length(t1), n2 = length(t2),
    Events_1 = sum(s1 == 1L), Events_2 = sum(s2 == 1L),
    Censoring_1_pct = 100 * mean(s1 == 0L), Censoring_2_pct = 100 * mean(s2 == 0L),
    At_risk_tau_1 = fit$Y1tau, At_risk_tau_2 = fit$Y2tau,
    inference_regular = fit$inference_regular,
    nonregular_reason = fit$nonregular_reason,
    D_tau_KM = fit$D, bias_hat = fit$bias_hat, D_bias_corrected = fit$D - fit$bias_hat,
    SE_Greenwood = fit$se,
    Normal_lo = fit$D - z * fit$se, Normal_hi = fit$D + z * fit$se,
    BiasCorr_lo = fit$D - fit$bias_hat - z * fit$se,
    BiasCorr_hi = fit$D - fit$bias_hat + z * fit$se,
    Bootstrap_SD = unname(bt["boot_sd"]), Bootstrap_mean = unname(bt["boot_mean"]),
    Percentile_lo = unname(bt["perc_lo"]), Percentile_hi = unname(bt["perc_hi"]),
    Basic_lo = unname(bt["basic_lo"]), Basic_hi = unname(bt["basic_hi"]),
    Bootstrap_nonregular_rate = unname(bt["boot_nonregular_rate"]),
    R_boot = R_boot, seed = seed,
    Logrank_p = 1 - stats::pchisq(lr$chisq, df = 1),
    stringsAsFactors = FALSE
  )
}

## Input: two censored samples and a fixed tau.
## Output: the cumulative divergence on the pooled Kaplan-Meier jump grid.
rd_cumulative <- function(t1, s1, t2, s2, tau) {
  fit <- D_km(t1, s1, t2, s2, tau, want_var = FALSE)
  out <- data.frame(
    left_time = fit$left, right_time = fit$grid[-1], width = fit$w,
    S_group1 = fit$S1, S_group2 = fit$S2,
    squared_difference = (fit$S1 - fit$S2)^2,
    contribution = (fit$S1 - fit$S2)^2 * fit$w
  )
  out$cumulative_D_tau <- cumsum(out$contribution)
  out
}

## Input: one analysis row, the samples, the step table, labels, the output
##        directory and the file stem.
## Output: the two file paths invisibly; writes <stem>_KM_curve.pdf and
##         <stem>_cumulative_divergence.pdf with base graphics (same names as the
##         ggplot renderer, so the manuscript never has to change).
rd_figure <- function(row, t1, s1, t2, s2, steps, labels, out_dir, stem) {
  km_file <- file.path(out_dir, paste0(stem, "_KM_curve.pdf"))
  cum_file <- file.path(out_dir, paste0(stem, "_cumulative_divergence.pdf"))
  group <- rep(1:2, c(length(t1), length(t2)))

  grDevices::pdf(km_file, width = 7.2, height = 4.8)
  graphics::par(mar = c(4.2, 4.4, 3, 1))
  fit <- survival::survfit(survival::Surv(c(t1, t2), c(s1, s2)) ~ group)
  plot(fit, col = c("firebrick", "steelblue"), lwd = 2, conf.int = TRUE,
       xlab = "Time", ylab = "Survival probability",
       main = paste0(row$Dataset, ": Kaplan-Meier curves"))
  graphics::abline(v = row$tau, lty = 2)
  graphics::legend("topright", legend = labels, col = c("firebrick", "steelblue"),
                   lwd = 2, bty = "n")
  graphics::mtext(sprintf("log-rank p = %.3g", row$Logrank_p), side = 3, line = 0.2,
                  adj = 0, cex = 0.85)
  grDevices::dev.off()

  grDevices::pdf(cum_file, width = 7.2, height = 4.8)
  graphics::par(mar = c(4.2, 4.4, 3, 1))
  ## The cumulative integral is piecewise LINEAR between the pooled knots.
  plot(c(0, steps$right_time), c(0, steps$cumulative_D_tau), type = "l", lwd = 2,
       xlab = "Time",
       ylab = expression(integral((hat(bar(F))[n[1]](u) - hat(bar(G))[n[2]](u))^2 * du, 0, t)),
       main = paste0(row$Dataset, ": cumulative divergence"))
  graphics::abline(v = row$tau, lty = 2)
  graphics::mtext(sprintf("tau = %.1f,  D = %.4f,  SE = %.4f", row$tau, row$D_tau_KM,
                          row$SE_Greenwood), side = 3, line = 0.2, adj = 0, cex = 0.85)
  grDevices::dev.off()
  invisible(c(km = km_file, cumulative = cum_file))
}

## Input: none. Output: TRUE invisibly; stops if any deterministic check fails.
## The checks use the veteran data, so they exercise the real-data code path.
rd_self_test <- function() {
  checks <- logical(0)
  checks["km engine self-checks"] <- isTRUE(verify_km_functions(verbose = FALSE))

  vet <- rd_load_data("veteran", "survival")
  t1 <- vet$time[vet$trt == 1]; s1 <- vet$status[vet$trt == 1]
  t2 <- vet$time[vet$trt == 2]; s2 <- vet$status[vet$trt == 2]
  tau <- rd_choose_tau(t1, s1, t2, s2, "min")
  fit <- D_km(t1, s1, t2, s2, tau)

  ## (a) D equals a direct integral built from survival::survfit.
  f1 <- survival::survfit(survival::Surv(t1, s1) ~ 1)
  f2 <- survival::survfit(survival::Surv(t2, s2) ~ 1)
  step <- function(f, u) {
    ev <- f$n.event > 0
    k <- findInterval(u, f$time[ev]); out <- rep(1, length(u))
    out[k > 0] <- f$surv[ev][k[k > 0]]; out
  }
  grid <- sort(unique(c(0, f1$time[f1$n.event > 0 & f1$time < tau],
                        f2$time[f2$n.event > 0 & f2$time < tau], tau)))
  left <- grid[-length(grid)]; w <- diff(grid)
  S1 <- step(f1, left); S2 <- step(f2, left)
  checks["D equals survfit-based integral"] <- isTRUE(all.equal(fit$D, sum((S1 - S2)^2 * w)))

  ## (b) The Greenwood variance equals the direct double sum of Section 4.3,
  ##     which is the check that guards the n1, n2 factors.
  gw <- function(f, u) {
    ev <- f$n.event > 0
    d <- f$n.event[ev]; Y <- f$n.risk[ev]
    inc <- ifelse(Y * (Y - d) > 0, d / (Y * (Y - d)), 0)
    k <- findInterval(u, f$time[ev]); out <- rep(0, length(u))
    out[k > 0] <- cumsum(inc)[k[k > 0]]; out
  }
  mn <- outer(left, left, pmin)
  V1 <- outer(S1, S1) * matrix(gw(f1, as.vector(mn)), length(left))
  V2 <- outer(S2, S2) * matrix(gw(f2, as.vector(mn)), length(left))
  direct <- 4 * sum(outer((S1 - S2) * w, (S1 - S2) * w) * (V1 + V2))
  checks["Greenwood variance equals direct double sum"] <-
    isTRUE(all.equal(fit$var, direct, tolerance = 1e-8))

  ## (c) The vectorised bootstrap equals explicit recomputation on the same draws.
  set.seed(11)
  bt <- boot_engine(t1, s1, t2, s2, tau, fit$D, R = 25L, return_draws = TRUE)
  i1 <- attr(bt, "indices_1"); i2 <- attr(bt, "indices_2")
  o1 <- order(t1); o2 <- order(t2)                    # boot_engine sorts each group
  slow <- vapply(seq_len(25L), function(r) {
    a <- i1[((r - 1L) * length(t1) + 1L):(r * length(t1))]
    b <- i2[((r - 1L) * length(t2) + 1L):(r * length(t2))]
    D_km(t1[o1][a], s1[o1][a], t2[o2][b], s2[o2][b], tau, want_var = FALSE)$D
  }, numeric(1))
  checks["vectorised bootstrap equals recomputation"] <-
    isTRUE(all.equal(attr(bt, "draws"), slow, tolerance = 1e-12))

  ## (d) The estimator is symmetric in the two groups.
  checks["estimator is symmetric"] <-
    isTRUE(all.equal(fit$D, D_km(t2, s2, t1, s1, tau, want_var = FALSE)$D))

  message("Real-data self-checks: ", sum(checks), "/", length(checks), " passed.")
  if (!all(checks)) {
    message("Failed: ", paste(names(checks)[!checks], collapse = ", "))
  }
  invisible(all(checks))
}
