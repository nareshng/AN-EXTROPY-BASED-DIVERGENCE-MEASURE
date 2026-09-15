## =============================================================================
## Shared functions for inference on the fixed, truncated divergence
##
##   D_tau = integral_0^tau {S_1(t) - S_2(t)}^2 dt
##
## under independent right censoring. The implementation uses the exact
## Kaplan--Meier step functions on [0, tau]. Inference requires tau to be fixed
## before a sample is analysed and observable follow-up through tau in both
## groups. Every function validates its inputs and explicitly reports samples or
## resamples for which these support conditions fail.
##
## The Greenwood calculation returned by D_km() is the direct variance estimate
##
##   4 sum_l d_1l A_1l^2 / {Y_1l (Y_1l-d_1l)}
## + 4 sum_l d_2l A_2l^2 / {Y_2l (Y_2l-d_2l)}.
##
## No additional division by n_1, n_2, or n is applied by calling scripts.
##
## This file also implements the truncated Cox--Czanner divergence used as the
## comparator in Tables 7--8.  Its Kaplan--Meier estimator uses linearly
## interpolated KM curves and a deterministic one-crossing total-variation
## rule.  Taking an absolute value separately on every pooled interval is not
## used: that construction accumulates sample-level sign noise and is
## inconsistent even for two exponential distributions.
## =============================================================================

## Inputs: value, label, bounds, openness, and integer flag.
## Returns: TRUE invisibly; stops when the scalar violates its contract.
assert_scalar_numeric <- function(x, name, lower = -Inf, upper = Inf,
                                  lower_open = FALSE, upper_open = FALSE,
                                  integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("%s must be one finite numeric value.", name), call. = FALSE)
  }
  lower_bad <- if (lower_open) x <= lower else x < lower
  upper_bad <- if (upper_open) x >= upper else x > upper
  if (lower_bad || upper_bad) {
    brackets <- paste0(if (lower_open) "(" else "[", lower, ", ", upper,
                       if (upper_open) ")" else "]")
    stop(sprintf("%s must lie in %s.", name, brackets), call. = FALSE)
  }
  if (integer && x != as.integer(x)) {
    stop(sprintf("%s must be an integer.", name), call. = FALSE)
  }
  invisible(TRUE)
}

## Inputs: follow-up times, 0/1 event indicators, and sample label.
## Returns: validated numeric time and integer status vectors.
validate_survival_sample <- function(time, status, name = "sample") {
  if (!is.numeric(time) || length(time) < 1L || anyNA(time) ||
      any(!is.finite(time)) || any(time < 0)) {
    stop(sprintf("%s times must be a nonempty finite, nonnegative numeric vector.", name),
         call. = FALSE)
  }
  if (!(is.numeric(status) || is.logical(status)) || length(status) != length(time) ||
      anyNA(status) || any(!status %in% c(0, 1))) {
    stop(sprintf("%s status must have the same length as time and contain only 0/1.", name),
         call. = FALSE)
  }
  list(time = as.numeric(time), status = as.integer(status))
}

## Input: diagnostic strings. Returns: unique nonempty reasons joined by semicolons.
collapse_reasons <- function(...) {
  x <- unlist(list(...), use.names = FALSE)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) paste(unique(x), collapse = "; ") else ""
}

## ---- Kaplan--Meier ------------------------------------------------------------
## Inputs: follow-up times and 0/1 event indicators.
## Returns: KM event grid, survival, risk, event, and Greenwood components.
km_fast <- function(time, status) {
  dat <- validate_survival_sample(time, status)
  time <- dat$time
  status <- dat$status
  n <- length(time)
  ev <- time[status == 1L]
  if (!length(ev)) {
    return(list(n = n, t = numeric(0), S = numeric(0), Y = numeric(0),
                d = numeric(0), gw = numeric(0), inc = numeric(0),
                terminal_event = FALSE, terminal_event_time = NA_real_))
  }

  ut <- sort(unique(ev))
  st <- sort(time)
  Y <- n - findInterval(ut, st, left.open = TRUE)  # #{i: time_i >= ut}
  d <- tabulate(match(ev, ut), nbins = length(ut))
  if (any(d <= 0L) || any(Y <= 0L) || any(d > Y)) {
    stop("Internal Kaplan--Meier risk-set inconsistency.", call. = FALSE)
  }

  S <- cumprod(1 - d / Y)
  den <- as.numeric(Y) * (Y - d)
  inc <- rep(NA_real_, length(ut))
  regular_jump <- den > 0
  inc[regular_jump] <- d[regular_jump] / den[regular_jump]

  ## A terminal event (Y=d) makes the usual Greenwood representation singular.
  ## Retain NA from that jump onward instead of silently replacing infinity by 0.
  gw <- cumsum(replace(inc, !is.finite(inc), 0))
  first_terminal <- which(!regular_jump)[1]
  if (!is.na(first_terminal)) gw[first_terminal:length(gw)] <- NA_real_

  list(n = n, t = ut, S = S, Y = Y, d = d, gw = gw, inc = inc,
       terminal_event = any(!regular_jump),
       terminal_event_time = if (any(!regular_jump)) ut[first_terminal] else NA_real_)
}

## Inputs: jump times/values, evaluation points, and pre-jump value.
## Returns: right-continuous step-function values at all evaluation points.
## Explicit positive indexing avoids R's zero-index dropping/recycling behaviour.
step_eval <- function(t, vals, u, init = 1) {
  if (!is.numeric(t) || !is.numeric(vals) || length(t) != length(vals) ||
      anyNA(t) || any(!is.finite(t)) || is.unsorted(t, strictly = TRUE)) {
    stop("t and vals must be numeric vectors of equal length and t must be strictly increasing.",
         call. = FALSE)
  }
  if (!is.numeric(u) || anyNA(u) || any(!is.finite(u))) {
    stop("u must contain only finite numeric values.", call. = FALSE)
  }
  assert_scalar_numeric(init, "init")
  k <- findInterval(u, t)
  out <- rep(init, length(u))
  pos <- k > 0L
  out[pos] <- vals[k[pos]]
  out
}

## ---- D_tau^KM, direct Greenwood variance, and plug-in bias --------------------
## Inputs: two censored samples, fixed horizon, and variance-computation flag.
## Returns: exact KM divergence plus support, variance, bias, and grid diagnostics.
D_km <- function(t1, s1, t2, s2, tau, want_var = TRUE) {
  d1 <- validate_survival_sample(t1, s1, "group 1")
  d2 <- validate_survival_sample(t2, s2, "group 2")
  assert_scalar_numeric(tau, "tau", lower = 0, lower_open = TRUE)
  if (!is.logical(want_var) || length(want_var) != 1L || is.na(want_var)) {
    stop("want_var must be TRUE or FALSE.", call. = FALSE)
  }

  t1 <- d1$time; s1 <- d1$status
  t2 <- d2$time; s2 <- d2$status
  k1 <- km_fast(t1, s1)
  k2 <- km_fast(t2, s2)
  grid <- sort(unique(c(0, k1$t[k1$t < tau], k2$t[k2$t < tau], tau)))
  left <- grid[-length(grid)]
  w <- diff(grid)
  if (!length(w) || any(w <= 0)) stop("Failed to construct the integration grid.", call. = FALSE)

  S1 <- step_eval(k1$t, k1$S, left)
  S2 <- step_eval(k2$t, k2$S, left)
  dif <- S1 - S2
  D <- sum(dif^2 * w)

  Y1tau <- sum(t1 >= tau)
  Y2tau <- sum(t2 >= tau)
  terminal1 <- any(k1$t < tau & k1$Y == k1$d)
  terminal2 <- any(k2$t < tau & k2$Y == k2$d)
  support_regular <- Y1tau > 0L && Y2tau > 0L && !terminal1 && !terminal2
  reason <- collapse_reasons(
    if (Y1tau == 0L) "group 1 has an empty risk set at tau",
    if (Y2tau == 0L) "group 2 has an empty risk set at tau",
    if (terminal1) "group 1 has a terminal event before tau",
    if (terminal2) "group 2 has a terminal event before tau"
  )

  common <- list(
    D = D, Y1tau = Y1tau, Y2tau = Y2tau,
    zero_risk_1 = Y1tau == 0L, zero_risk_2 = Y2tau == 0L,
    terminal_before_tau_1 = terminal1, terminal_before_tau_2 = terminal2,
    support_regular = support_regular, nonregular_reason = reason,
    last_observed_1 = max(t1), last_observed_2 = max(t2),
    grid = grid, left = left, w = w, S1 = S1, S2 = S2
  )
  if (!want_var) return(common)

  if (!support_regular) {
    return(c(common, list(se = NA_real_, var = NA_real_, bias_hat = NA_real_,
                          variance_positive = FALSE, inference_regular = FALSE)))
  }

  ## A_r(t_l) = integral_{t_l}^tau (S1-S2) S_r du, exactly on the pooled grid.
  tailF <- rev(cumsum(rev(dif * S1 * w)))
  tailG <- rev(cumsum(rev(dif * S2 * w)))
  gcomp <- function(k, tailv) {
    sel <- k$t < tau
    if (!any(sel)) return(0)
    pos <- match(k$t[sel], grid)
    increments <- k$inc[sel]
    if (anyNA(pos) || any(!is.finite(increments))) return(NA_real_)
    sum(increments * tailv[pos]^2)
  }
  gF <- gcomp(k1, tailF)
  gG <- gcomp(k2, tailG)
  var_D <- 4 * gF + 4 * gG

  G1 <- step_eval(k1$t, k1$gw, left, init = 0)
  G2 <- step_eval(k2$t, k2$gw, left, init = 0)
  bias_hat <- sum((S1^2 * G1 + S2^2 * G2) * w)
  variance_positive <- is.finite(var_D) && var_D > 0
  inference_regular <- support_regular && variance_positive && is.finite(bias_hat)
  if (!variance_positive) {
    reason <- collapse_reasons(reason, "estimated first-order variance is zero or non-finite")
  }
  if (!is.finite(bias_hat)) {
    reason <- collapse_reasons(reason, "plug-in bias estimate is non-finite")
  }

  common$nonregular_reason <- reason
  c(common, list(se = if (is.finite(var_D) && var_D >= 0) sqrt(var_D) else NA_real_,
                 var = var_D, bias_hat = bias_hat,
                 variance_positive = variance_positive,
                 inference_regular = inference_regular))
}

## ---- truncated Cox--Czanner comparator --------------------------------------
## Input: a km_fast() result and evaluation points.
## Output: the linearly interpolated KM curve used by Cox and Czanner (2016).
km_linear_eval <- function(km, u) {
  if (!is.list(km) || !all(c("t", "S") %in% names(km))) {
    stop("km must be a result returned by km_fast().", call. = FALSE)
  }
  if (!is.numeric(u) || anyNA(u) || any(!is.finite(u)) || any(u < 0)) {
    stop("u must contain finite nonnegative values.", call. = FALSE)
  }
  if (!length(km$t)) return(rep(1, length(u)))
  if (km$t[1L] <= 0) {
    stop("The linearly interpolated Cox--Czanner estimator requires positive event times.",
         call. = FALSE)
  }
  stats::approx(
    x = c(0, km$t), y = c(1, km$S), xout = u,
    method = "linear", rule = 2, ties = "ordered"
  )$y
}

## Inputs: values of two linearly interpolated survival curves on one grid.
## Output: estimated total variation when the hazard difference changes sign at
## most once.  The maximum selects the deterministic global extremum suggested
## by Cox and Czanner for a single crossing.  With no crossing it reduces
## asymptotically to the absolute value of the global signed sum.
dcc_one_crossing_value <- function(F_values, G_values) {
  if (!is.numeric(F_values) || !is.numeric(G_values) ||
      length(F_values) != length(G_values) || length(F_values) < 2L ||
      anyNA(F_values) || anyNA(G_values) ||
      any(!is.finite(F_values)) || any(!is.finite(G_values))) {
    stop("F_values and G_values must be equal-length finite numeric vectors.",
         call. = FALSE)
  }
  increments <-
    F_values[-1L] * G_values[-length(G_values)] -
    F_values[-length(F_values)] * G_values[-1L]
  cumulative <- c(0, cumsum(increments))
  total <- cumulative[length(cumulative)]
  max(abs(cumulative) + abs(total - cumulative))
}

## Inputs: two censored samples and a fixed truncation horizon.
## Output: the Cox--Czanner KM estimator on [0,tau].  The simulation families
## used in Tables 7--8 have at most one crossing of their hazard functions.
DCC_km <- function(t1, s1, t2, s2, tau) {
  d1 <- validate_survival_sample(t1, s1, "group 1")
  d2 <- validate_survival_sample(t2, s2, "group 2")
  assert_scalar_numeric(tau, "tau", lower = 0, lower_open = TRUE)
  k1 <- km_fast(d1$time, d1$status)
  k2 <- km_fast(d2$time, d2$status)
  grid <- sort(unique(c(0, k1$t[k1$t <= tau], k2$t[k2$t <= tau], tau)))
  F_values <- km_linear_eval(k1, grid)
  G_values <- km_linear_eval(k2, grid)
  value <- dcc_one_crossing_value(F_values, G_values)
  if (!is.finite(value) || value < -sqrt(.Machine$double.eps)) {
    stop("The Cox--Czanner estimator is non-finite or negative.", call. = FALSE)
  }
  max(value, 0)
}

## ---- vectorised within-group pairs bootstrap ---------------------------------
## A resample is represented by multinomial row counts. The resulting KM curve is
## algebraically identical to recomputing the estimator after separately resampling
## the observed (T, delta) pairs in each group.
## Inputs: one censored sample. Returns: sorted matrices for bootstrap KM curves.
prep_group <- function(time, status) {
  dat <- validate_survival_sample(time, status)
  o <- order(dat$time)
  time <- dat$time[o]
  status <- dat$status[o]
  n <- length(time)
  ev <- status == 1L
  ut <- sort(unique(time[ev]))
  K <- length(ut)
  U <- matrix(0, n, n)
  U[lower.tri(U, diag = TRUE)] <- 1
  first <- if (K) match(ut, time) else integer(0)
  Emat <- matrix(0, n, K)
  if (K) {
    for (k in seq_len(K)) Emat[time == ut[k] & ev, k] <- 1
  }
  list(n = n, time = time, status = status, ut = ut, K = K,
       U = U, first = first, Emat = Emat)
}

## Inputs: two censored samples, fixed horizon/estimate, resample count, and alpha.
## Returns: bootstrap intervals, scale/support diagnostics, and optional draws.
boot_engine <- function(t1, s1, t2, s2, tau, D_hat, R = 1000,
                        alpha = 0.05, return_draws = FALSE) {
  assert_scalar_numeric(tau, "tau", lower = 0, lower_open = TRUE)
  assert_scalar_numeric(D_hat, "D_hat", lower = 0)
  assert_scalar_numeric(R, "R", lower = 2, integer = TRUE)
  assert_scalar_numeric(alpha, "alpha", lower = 0, upper = 1,
                        lower_open = TRUE, upper_open = TRUE)
  if (!is.logical(return_draws) || length(return_draws) != 1L || is.na(return_draws)) {
    stop("return_draws must be TRUE or FALSE.", call. = FALSE)
  }

  R <- as.integer(R)
  p1 <- prep_group(t1, s1)
  p2 <- prep_group(t2, s2)
  grid <- sort(unique(c(0, p1$ut[p1$ut < tau], p2$ut[p2$ut < tau], tau)))
  left <- grid[-length(grid)]
  w <- diff(grid)

  surv_mat <- function(p) {
    idx <- sample.int(p$n, p$n * R, replace = TRUE)
    encoded <- idx + p$n * rep(0:(R - 1L), each = p$n)
    M <- matrix(tabulate(encoded, nbins = p$n * R), nrow = R, byrow = TRUE)
    Ytau <- if (any(p$time >= tau)) {
      rowSums(M[, p$time >= tau, drop = FALSE])
    } else {
      rep(0, R)
    }

    if (p$K == 0L) {
      return(list(S = matrix(1, R, length(left)), Ytau = Ytau,
                  terminal = rep(FALSE, R), indices = idx))
    }

    Y <- (M %*% p$U)[, p$first, drop = FALSE]
    d <- M %*% p$Emat
    ratio <- matrix(0, nrow = R, ncol = p$K)
    positive_Y <- Y > 0
    ratio[positive_Y] <- d[positive_Y] / Y[positive_Y]
    if (any(ratio < -sqrt(.Machine$double.eps) |
            ratio > 1 + sqrt(.Machine$double.eps), na.rm = TRUE)) {
      stop("Internal bootstrap risk-set inconsistency.", call. = FALSE)
    }
    ratio <- pmin(pmax(ratio, 0), 1)  ## FIX: scalar-first pmin/pmax dropped dim()

    ## Row-wise cumulative product, preserving an exact zero after a terminal jump.
    S_event <- matrix(1, nrow = R, ncol = p$K)
    running <- rep(1, R)
    for (k in seq_len(p$K)) {
      running <- running * (1 - ratio[, k])
      S_event[, k] <- running
    }
    kk <- findInterval(left, p$ut)
    S_left <- cbind(1, S_event)[, kk + 1L, drop = FALSE]
    before_tau <- p$ut < tau
    terminal <- if (any(before_tau)) {
      rowSums((d[, before_tau, drop = FALSE] > 0) &
                (d[, before_tau, drop = FALSE] == Y[, before_tau, drop = FALSE])) > 0
    } else {
      rep(FALSE, R)
    }
    list(S = S_left, Ytau = Ytau, terminal = terminal, indices = idx)
  }

  b1 <- surv_mat(p1)
  b2 <- surv_mat(p2)
  Db <- as.vector(((b1$S - b2$S)^2) %*% w)
  if (any(!is.finite(Db)) || any(Db < -sqrt(.Machine$double.eps))) {
    stop("Bootstrap produced a non-finite or negative divergence.", call. = FALSE)
  }
  Db <- pmax(Db, 0)
  q <- stats::quantile(Db, c(alpha / 2, 1 - alpha / 2),
                       names = FALSE, type = 8)
  boot_regular <- b1$Ytau > 0 & b2$Ytau > 0 & !b1$terminal & !b2$terminal
  ans <- c(
    perc_lo = q[1], perc_hi = q[2],
    basic_lo = 2 * D_hat - q[2], basic_hi = 2 * D_hat - q[1],
    boot_sd = stats::sd(Db), boot_mean = mean(Db),
    boot_zero_risk_rate_1 = mean(b1$Ytau == 0),
    boot_zero_risk_rate_2 = mean(b2$Ytau == 0),
    boot_terminal_rate_1 = mean(b1$terminal),
    boot_terminal_rate_2 = mean(b2$terminal),
    boot_nonregular_rate = mean(!boot_regular),
    boot_regular_n = sum(boot_regular)
  )
  if (return_draws) {
    attr(ans, "draws") <- Db
    attr(ans, "regular") <- boot_regular
    attr(ans, "indices_1") <- b1$indices
    attr(ans, "indices_2") <- b2$indices
  }
  ans
}

## ---- censoring calibration ----------------------------------------------------
## Here lambda is an exponential rate, matching rexp(..., rate = lambda).
## P(censored) = P(C < X) = 1 - E{exp(-gamma X)}.
## Inputs: lifetime rate and censoring probability. Returns: censoring rate.
exp_censor_rate <- function(lambda, p) {
  assert_scalar_numeric(lambda, "lambda", lower = 0, lower_open = TRUE)
  assert_scalar_numeric(p, "p", lower = 0, upper = 1,
                        lower_open = TRUE, upper_open = TRUE)
  p * lambda / (1 - p)
}

## Inputs: censoring rate and Weibull parameters. Returns: P(C < X).
weibull_censor_probability <- function(gamma, shape, scale) {
  assert_scalar_numeric(gamma, "gamma", lower = 0)
  assert_scalar_numeric(shape, "shape", lower = 0, lower_open = TRUE)
  assert_scalar_numeric(scale, "scale", lower = 0, lower_open = TRUE)
  if (gamma == 0) return(0)
  ## P(C < X) = integral gamma exp(-gamma t) S_X(t) dt. This survival-form
  ## expression stays finite at zero even when a Weibull density with shape
  ## below one is unbounded there.
  stats::integrate(
    function(x) gamma * exp(-gamma * x - (x / scale)^shape),
    lower = 0, upper = Inf, subdivisions = 2000,
    rel.tol = 1e-10, abs.tol = 1e-12, stop.on.error = TRUE
  )$value
}

## Inputs: Weibull parameters and target P(C < X). Returns: calibrated rate.
weibull_censor_rate <- function(shape, scale, p) {
  assert_scalar_numeric(shape, "shape", lower = 0, lower_open = TRUE)
  assert_scalar_numeric(scale, "scale", lower = 0, lower_open = TRUE)
  assert_scalar_numeric(p, "p", lower = 0, upper = 1,
                        lower_open = TRUE, upper_open = TRUE)
  objective <- function(gamma) weibull_censor_probability(gamma, shape, scale) - p
  upper <- 1
  iter <- 0L
  while (objective(upper) < 0) {
    upper <- 2 * upper
    iter <- iter + 1L
    if (!is.finite(upper) || iter > 100L) {
      stop("Failed to bracket the Weibull censoring rate.", call. = FALSE)
    }
  }
  stats::uniroot(objective, c(0, upper), tol = 1e-10)$root
}

## ---- population targets -------------------------------------------------------
## Inputs: two survival functions and fixed horizon. Returns: population D_tau.
true_D_tau <- function(S1, S2, tau) {
  if (!is.function(S1) || !is.function(S2)) stop("S1 and S2 must be functions.", call. = FALSE)
  assert_scalar_numeric(tau, "tau", lower = 0, lower_open = TRUE)
  stats::integrate(function(t) (S1(t) - S2(t))^2, 0, tau,
                   subdivisions = 2000, rel.tol = 1e-10,
                   abs.tol = 1e-12, stop.on.error = TRUE)$value
}

## Inputs: two survival functions, their densities, a fixed horizon, and a
## positive endpoint-transformation power. Output: population Cox--Czanner
## divergence on [0,tau]. A shape-adaptive power is supplied for Weibull
## scenarios so density singularities at zero are integrated stably.
true_DCC_tau <- function(S1, S2, f1, f2, tau, transform_power = 1) {
  if (!is.function(S1) || !is.function(S2) ||
      !is.function(f1) || !is.function(f2)) {
    stop("S1, S2, f1 and f2 must be functions.", call. = FALSE)
  }
  assert_scalar_numeric(tau, "tau", lower = 0, lower_open = TRUE)
  assert_scalar_numeric(
    transform_power, "transform_power", lower = 0, lower_open = TRUE
  )
  transformed_integrand <- function(u) {
    t <- tau * u^transform_power
    jacobian <- transform_power * tau * u^(transform_power - 1)
    value <- abs(S1(t) * f2(t) - S2(t) * f1(t)) * jacobian
    value[!is.finite(value) & u == 0] <- 0
    value
  }
  stats::integrate(
    transformed_integrand, 0, 1, subdivisions = 2000,
    rel.tol = 1e-10, abs.tol = 1e-12, stop.on.error = TRUE
  )$value
}

## family = "exp": cfg = c(lambda1, lambda2), with rates lambda_r.
## family = "weibull": cfg = c(shape1, scale1, shape2, scale2).
## Inputs: family, parameters, target censoring, and population tau quantile.
## Returns: generators, survival functions, target, tau, and calibration details.
make_scenario <- function(family, cfg, cens, tau_quantile = 0.80) {
  family <- match.arg(family, c("exp", "weibull"))
  assert_scalar_numeric(cens, "cens", lower = 0, upper = 1,
                        lower_open = TRUE, upper_open = TRUE)
  assert_scalar_numeric(tau_quantile, "tau_quantile", lower = 0, upper = 1,
                        lower_open = TRUE, upper_open = TRUE)
  if (!is.numeric(cfg) || anyNA(cfg) || any(!is.finite(cfg)) || any(cfg <= 0)) {
    stop("cfg must contain only finite positive values.", call. = FALSE)
  }

  if (family == "exp") {
    if (length(cfg) != 2L) stop("An exponential cfg must be c(lambda1, lambda2).", call. = FALSE)
    l1 <- cfg[1]; l2 <- cfg[2]
    g1 <- exp_censor_rate(l1, cens); g2 <- exp_censor_rate(l2, cens)
    tau <- min(stats::qexp(tau_quantile, rate = l1),
               stats::qexp(tau_quantile, rate = l2))
    S1 <- function(t) exp(-l1 * t); S2 <- function(t) exp(-l2 * t)
    f1 <- function(t) stats::dexp(t, rate = l1)
    f2 <- function(t) stats::dexp(t, rate = l2)
    gen <- function(n1, n2) {
      assert_scalar_numeric(n1, "n1", lower = 2, integer = TRUE)
      assert_scalar_numeric(n2, "n2", lower = 2, integer = TRUE)
      list(X = stats::rexp(n1, rate = l1), Y = stats::rexp(n2, rate = l2),
           C = stats::rexp(n1, rate = g1), E = stats::rexp(n2, rate = g2))
    }
    label <- sprintf("(%.1f, %.1f)", l1, l2)
    cens1 <- g1 / (l1 + g1)
    cens2 <- g2 / (l2 + g2)
    dcc_transform_power <- 1
  } else {
    if (length(cfg) != 4L) {
      stop("A Weibull cfg must be c(shape1, scale1, shape2, scale2).", call. = FALSE)
    }
    k1 <- cfg[1]; scale1 <- cfg[2]; k2 <- cfg[3]; scale2 <- cfg[4]
    g1 <- weibull_censor_rate(k1, scale1, cens)
    g2 <- weibull_censor_rate(k2, scale2, cens)
    tau <- min(stats::qweibull(tau_quantile, shape = k1, scale = scale1),
               stats::qweibull(tau_quantile, shape = k2, scale = scale2))
    S1 <- function(t) exp(-(t / scale1)^k1)
    S2 <- function(t) exp(-(t / scale2)^k2)
    f1 <- function(t) stats::dweibull(t, shape = k1, scale = scale1)
    f2 <- function(t) stats::dweibull(t, shape = k2, scale = scale2)
    gen <- function(n1, n2) {
      assert_scalar_numeric(n1, "n1", lower = 2, integer = TRUE)
      assert_scalar_numeric(n2, "n2", lower = 2, integer = TRUE)
      list(X = stats::rweibull(n1, shape = k1, scale = scale1),
           Y = stats::rweibull(n2, shape = k2, scale = scale2),
           C = stats::rexp(n1, rate = g1), E = stats::rexp(n2, rate = g2))
    }
    label <- sprintf("((%.1f, %.1f), (%.1f, %.1f))", k1, scale1, k2, scale2)
    cens1 <- weibull_censor_probability(g1, k1, scale1)
    cens2 <- weibull_censor_probability(g2, k2, scale2)
    dcc_transform_power <- max(1, 2 / min(k1, k2))
  }

  if (max(abs(c(cens1, cens2) - cens)) > 1e-7) {
    stop("Censoring-rate calibration self-check failed.", call. = FALSE)
  }
  list(
    gen = gen, tau = tau,
    D = true_D_tau(S1, S2, tau),
    DCC = true_DCC_tau(
      S1, S2, f1, f2, tau,
      transform_power = dcc_transform_power
    ),
    S1 = S1, S2 = S2, f1 = f1, f2 = f2,
    label = label, g1 = g1, g2 = g2,
    cens1 = cens1, cens2 = cens2,
    P_obs_gt_tau_1 = S1(tau) * exp(-g1 * tau),
    P_obs_gt_tau_2 = S2(tau) * exp(-g2 * tau)
  )
}

## ---- deterministic implementation checks -------------------------------------
## Input: verbosity flag. Returns: TRUE invisibly; stops after a failed check.
verify_km_functions <- function(verbose = TRUE) {
  checks <- logical(0)

  k <- km_fast(c(1, 2, 3, 4), c(1, 0, 1, 0))
  checks["known KM times"] <- isTRUE(all.equal(k$t, c(1, 3)))
  checks["known KM risks"] <- isTRUE(all.equal(k$Y, c(4, 2)))
  checks["known KM survival"] <- isTRUE(all.equal(k$S, c(0.75, 0.375)))
  checks["step zero-index guard"] <- isTRUE(all.equal(
    step_eval(c(1, 3), c(0.75, 0.375), c(0, 1, 2, 3, 4)),
    c(1, 0.75, 0.75, 0.375, 0.375)
  ))

  same <- D_km(c(1, 2, 3, 4), c(1, 0, 1, 0),
               c(1, 2, 3, 4), c(1, 0, 1, 0), tau = 2.5)
  checks["identical KM divergence"] <- isTRUE(all.equal(same$D, 0))
  checks["degenerate variance flagged"] <- !same$inference_regular && same$var == 0
  checks["identical KM Cox--Czanner divergence"] <- isTRUE(all.equal(
    DCC_km(c(1, 2, 3, 4), c(1, 0, 1, 0),
           c(1, 2, 3, 4), c(1, 0, 1, 0), tau = 2.5),
    0
  ))

  dcc_forward <- DCC_km(
    c(0.5, 1, 1.5, 2, 3), c(1, 0, 1, 0, 1),
    c(0.4, 1.2, 1.8, 2.5, 3.2, 4), c(1, 1, 0, 1, 0, 1),
    tau = 2.4
  )
  dcc_reverse <- DCC_km(
    c(0.4, 1.2, 1.8, 2.5, 3.2, 4), c(1, 1, 0, 1, 0, 1),
    c(0.5, 1, 1.5, 2, 3), c(1, 0, 1, 0, 1),
    tau = 2.4
  )
  checks["Cox--Czanner estimator is symmetric"] <-
    isTRUE(all.equal(dcc_forward, dcc_reverse, tolerance = 1e-14))
  regression_F <- c(1, 0.9, 0.4, 0.3, 0.2)
  regression_G <- c(1, 0.7, 0.65, 0.2, 0.15)
  checks["one-crossing rule is not intervalwise variation"] <-
    isTRUE(all.equal(
      dcc_one_crossing_value(regression_F, regression_G),
      0.395,
      tolerance = 1e-14
    ))
  checks["known linearly interpolated Cox--Czanner value"] <-
    isTRUE(all.equal(
      DCC_km(c(1, 3), c(1, 1), c(2, 4), c(1, 1), tau = 2.5),
      0.34375,
      tolerance = 1e-14
    ))
  dcc_ties <- DCC_km(
    c(1, 1, 3), c(1, 1, 1),
    c(1, 2, 2), c(1, 1, 1),
    tau = 2.5
  )
  checks["tied and simultaneous failures are finite"] <-
    is.finite(dcc_ties) && dcc_ties >= 0
  checks["no-event Cox--Czanner value"] <- isTRUE(all.equal(
    DCC_km(c(1, 2), c(0, 0), c(1.5, 3), c(0, 0), tau = 2.5),
    0
  ))

  terminal <- D_km(c(1, 2), c(0, 1), c(1, 3), c(1, 0), tau = 2.5)
  checks["empty risk set flagged"] <- terminal$zero_risk_1 && !terminal$support_regular &&
    is.na(terminal$se)

  ## Verify the collapsed Greenwood calculation against the direct double-sum
  ## covariance integral on a small pooled step grid. This specifically guards
  ## against adding or omitting an extra n, n1, or n2 factor.
  vt1 <- c(0.5, 1, 1.5, 2, 3)
  vs1 <- c(1, 0, 1, 0, 0)
  vt2 <- c(0.4, 1.2, 1.8, 2.5, 3.2, 4)
  vs2 <- c(1, 1, 0, 1, 0, 1)
  vfit <- D_km(vt1, vs1, vt2, vs2, tau = 2.4, want_var = TRUE)
  vk1 <- km_fast(vt1, vs1)
  vk2 <- km_fast(vt2, vs2)
  vmin <- outer(vfit$left, vfit$left, pmin)
  vG1 <- matrix(step_eval(vk1$t, vk1$gw, as.vector(vmin), init = 0),
                 nrow = length(vfit$left))
  vG2 <- matrix(step_eval(vk2$t, vk2$gw, as.vector(vmin), init = 0),
                 nrow = length(vfit$left))
  vcov1 <- outer(vfit$S1, vfit$S1) * vG1
  vcov2 <- outer(vfit$S2, vfit$S2) * vG2
  vweighted_diff <- (vfit$S1 - vfit$S2) * vfit$w
  vdirect <- 4 * sum(outer(vweighted_diff, vweighted_diff) * (vcov1 + vcov2))
  checks["direct Greenwood variance scaling"] <-
    isTRUE(all.equal(vfit$var, vdirect, tolerance = 1e-12))

  sc_exp <- make_scenario("exp", c(0.5, 1), cens = 0.3)
  sc_wei <- make_scenario("weibull", c(1.3, 2, 0.8, 1), cens = 0.3)
  sc_wei_singular <- make_scenario("weibull", c(0.5, 1, 1, 1), cens = 0.3)
  checks["exponential censor calibration"] <- max(abs(c(sc_exp$cens1, sc_exp$cens2) - 0.3)) < 1e-10
  dcc_exp_closed <- abs(1 - 0.5) / (1 + 0.5) *
    (1 - exp(-(1 + 0.5) * sc_exp$tau))
  checks["exponential Cox--Czanner target"] <-
    isTRUE(all.equal(sc_exp$DCC, dcc_exp_closed, tolerance = 1e-9))
  sc_wei_exp <- make_scenario("weibull", c(1, 2, 1, 1), cens = 0.3)
  checks["Weibull shape-one target equals exponential target"] <-
    isTRUE(all.equal(sc_wei_exp$DCC, sc_exp$DCC, tolerance = 1e-9))
  dense_t <- seq(0, sc_exp$tau, length.out = 10001L)
  dense_dcc <- dcc_one_crossing_value(sc_exp$S1(dense_t), sc_exp$S2(dense_t))
  checks["one-crossing aggregation matches exponential target"] <-
    isTRUE(all.equal(dense_dcc, sc_exp$DCC, tolerance = 1e-4))
  checks["Weibull censor calibration"] <- max(abs(c(sc_wei$cens1, sc_wei$cens2) - 0.3)) < 1e-7
  checks["Weibull shape-below-one calibration"] <-
    max(abs(c(sc_wei_singular$cens1, sc_wei_singular$cens2) - 0.3)) < 1e-7

  ## Exact vectorised bootstrap check against explicit recomputation using the
  ## same resampled row indices.
  t1 <- c(0.5, 1, 1.5, 2, 3); s1 <- c(1, 0, 1, 0, 1)
  t2 <- c(0.4, 1.2, 1.8, 2.5, 3.2, 4); s2 <- c(1, 1, 0, 1, 0, 1)
  D0 <- D_km(t1, s1, t2, s2, tau = 2.4, want_var = FALSE)$D
  set.seed(917)
  bt <- boot_engine(t1, s1, t2, s2, tau = 2.4, D_hat = D0,
                    R = 7, return_draws = TRUE)
  i1 <- attr(bt, "indices_1"); i2 <- attr(bt, "indices_2")
  slow <- vapply(seq_len(7), function(r) {
    a <- i1[((r - 1L) * length(t1) + 1L):(r * length(t1))]
    b <- i2[((r - 1L) * length(t2) + 1L):(r * length(t2))]
    D_km(t1[a], s1[a], t2[b], s2[b], tau = 2.4, want_var = FALSE)$D
  }, numeric(1))
  checks["vectorised bootstrap"] <- isTRUE(all.equal(attr(bt, "draws"), slow,
                                                       tolerance = 1e-13))

  ## Optional comparison with survival::survfit. Its absence is not a failure
  ## because the simulations themselves use only base R.
  if (requireNamespace("survival", quietly = TRUE)) {
    set.seed(41)
    ok_survival <- TRUE
    for (r in seq_len(10)) {
      n <- 25L + r
      x <- round(stats::rweibull(n, shape = 0.8, scale = 1), 1)
      cc <- stats::rexp(n, rate = 0.5)
      tt <- pmin(x, cc); ss <- as.integer(x <= cc)
      fit <- survival::survfit(survival::Surv(tt, ss) ~ 1)
      kk <- km_fast(tt, ss); event_rows <- fit$n.event > 0
      ok_survival <- ok_survival &&
        isTRUE(all.equal(fit$time[event_rows], kk$t)) &&
        isTRUE(all.equal(fit$surv[event_rows], kk$S)) &&
        isTRUE(all.equal(fit$n.risk[event_rows], kk$Y)) &&
        isTRUE(all.equal(fit$n.event[event_rows], kk$d))
    }
    checks["survfit agreement"] <- ok_survival
  }

  passed <- all(checks)
  if (verbose) {
    message("KM self-checks: ", sum(checks), "/", length(checks), " passed.")
    if (!passed) message("Failed checks: ", paste(names(checks)[!checks], collapse = ", "))
  }
  isTRUE(passed)
}
