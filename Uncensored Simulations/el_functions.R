## =============================================================================
## Empirical likelihood ratio for a single mean constraint, in base R.
##
## Given w_1, ..., w_n, this returns the value of
##
##     -2 log R = -2 log max { prod_{i} (n p_i) : p_i > 0, sum p_i = 1,
##                             sum p_i w_i = 0 },
##
## which is the quantity emplik::el.test(w, mu = 0)$"-2LLR" reports.  The
## maximiser is available in closed form through the Lagrange dual: with lambda
## solving sum_i w_i / (1 + lambda w_i) = 0 we have p_i = 1 / {n (1 + lambda
## w_i)}, and then -2 log R = 2 sum_i log(1 + lambda w_i).  Both constraints hold
## identically at that lambda, so no numerical optimisation over the simplex is
## needed.
##
## Providing this here removes the hard dependency on the emplik package, and
## lets el_self_test() cross-check against emplik when it happens to be
## installed.  Requires: base R only.
## =============================================================================

## Input: a numeric vector of centred observations.
## Output: -2 log R at mu = 0; Inf when 0 is outside the convex hull of w.
el_m2llr <- function(w) {
  if (!is.numeric(w) || !length(w) || anyNA(w) || !all(is.finite(w))) return(Inf)
  if (min(w) >= 0 || max(w) <= 0) return(Inf)      # the constraint is infeasible
  lower <- (-1 / max(w)) * (1 - 1e-12)             # keeps every 1 + lambda w_i > 0
  upper <- (-1 / min(w)) * (1 - 1e-12)
  estimating_equation <- function(lambda) sum(w / (1 + lambda * w))
  lambda <- tryCatch(
    stats::uniroot(
      estimating_equation, c(lower, upper),
      tol = .Machine$double.eps^0.75
    )$root,
    error = function(e) NA_real_
  )
  if (!is.finite(lambda)) return(Inf)
  value <- 2 * sum(log1p(lambda * w))
  if (!is.finite(value)) return(Inf)
  ## -2 log R is nonnegative; clamp the O(1e-16) noise that appears when the
  ## constraint is already satisfied at lambda = 0.  Without this clamp the
  ## calling root search sees a spurious Inf at its own starting point and
  ## reports the interval as unavailable.
  max(value, 0)
}

## Input: verbosity flag. Output: TRUE invisibly; stops after a failed check.
el_self_test <- function(verbose = TRUE) {
  checks <- logical(0)

  checks["infeasible on one-sided data"] <-
    is.infinite(el_m2llr(c(1, 2, 3))) && is.infinite(el_m2llr(c(-1, -2, -3)))
  checks["zero at an exactly centred sample"] <-
    isTRUE(all.equal(el_m2llr(c(-1, 1)), 0, tolerance = 1e-12))
  checks["never negative"] <- {
    set.seed(11)
    all(vapply(seq_len(200), function(i) el_m2llr(stats::rnorm(30)) >= 0, logical(1)))
  }

  ## The dual solution must satisfy both constraints exactly.
  set.seed(12)
  feasible <- TRUE
  for (i in seq_len(50)) {
    w <- stats::rnorm(25)
    lower <- (-1 / max(w)) * (1 - 1e-12)
    upper <- (-1 / min(w)) * (1 - 1e-12)
    lambda <- stats::uniroot(function(l) sum(w / (1 + l * w)), c(lower, upper),
                             tol = .Machine$double.eps^0.75)$root
    p <- 1 / (length(w) * (1 + lambda * w))
    feasible <- feasible && all(p > 0) &&
      abs(sum(p) - 1) < 1e-10 && abs(sum(p * w)) < 1e-10 &&
      isTRUE(all.equal(-2 * sum(log(length(w) * p)), el_m2llr(w), tolerance = 1e-10))
  }
  checks["dual solution is feasible and attains the maximum"] <- feasible

  ## First-order agreement with the chi-square pivot for a well-separated mean.
  set.seed(13)
  w <- stats::rnorm(400) + 0.25
  checks["first-order agreement with n mean^2 / var"] <-
    abs(el_m2llr(w) - 400 * mean(w)^2 / stats::var(w)) < 0.5

  ## Optional agreement with emplik, when the package is installed.
  if (requireNamespace("emplik", quietly = TRUE)) {
    set.seed(14)
    agree <- TRUE
    for (i in seq_len(25)) {
      w <- stats::rnorm(40) + stats::runif(1, -0.3, 0.3)
      agree <- agree && isTRUE(all.equal(
        el_m2llr(w),
        as.numeric(emplik::el.test(w, mu = 0)$"-2LLR"),
        tolerance = 1e-6
      ))
    }
    checks["agreement with emplik::el.test"] <- agree
  }

  passed <- all(checks)
  if (verbose) {
    message("Empirical-likelihood self-checks: ", sum(checks), "/", length(checks),
            " passed.")
    if (!passed) {
      message("Failed checks: ", paste(names(checks)[!checks], collapse = ", "))
    }
  }
  if (!passed) stop("Empirical-likelihood self-check failed.", call. = FALSE)
  invisible(TRUE)
}
