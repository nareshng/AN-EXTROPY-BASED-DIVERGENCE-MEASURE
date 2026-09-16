# =============================================================================
# Section 5.3 – Coverage Probability and Average Length
# for THREE confidence interval methods under Weibull distributions:
#
#   1. U-statistic + JEL
#   2. U-statistic + Normal approximation
#   3. Empirical estimator + Normal approximation
#
# Model:
#   X ~ Weibull(shape1, scale1)
#   Y ~ Weibull(shape2, scale2)
#
#
#
# =============================================================================


# =============================================================================
# Package
# =============================================================================

if (!requireNamespace("emplik", quietly = TRUE)) {
  stop("Package 'emplik' 1.3-3 is required; see README.md for installation.")
}
if (utils::packageVersion("emplik") != "1.3.3") {
  stop("This reproducibility script requires package 'emplik' version 1.3-3.")
}

library(emplik)


# =============================================================================
# 1. True divergence for Weibull distributions
# =============================================================================

true_D_weibull <- function(shape1, scale1, shape2, scale2) {
  if (shape1 <= 0 || shape2 <= 0 || scale1 <= 0 || scale2 <= 0) {
    stop("All Weibull shape and scale parameters must be positive.")
  }
  
  # integral_0^inf S_X(x)^2 dx
  term1 <- scale1 * gamma(1 + 1 / shape1) / (2^(1 / shape1))
  
  # integral_0^inf S_Y(x)^2 dx
  term2 <- scale2 * gamma(1 + 1 / shape2) / (2^(1 / shape2))
  
  # integral_0^inf S_X(x) S_Y(x) dx
  # Closed form only when shape1 == shape2.
  if (abs(shape1 - shape2) < 1e-12) {
    shape <- shape1
    c_cross <- scale1^(-shape) + scale2^(-shape)
    cross <- gamma(1 + 1 / shape) / (c_cross^(1 / shape))
  } else {
    cross_integrand <- function(x) {
      exp(- (x / scale1)^shape1 - (x / scale2)^shape2)
    }
    
    cross <- integrate(
      cross_integrand,
      lower = 0,
      upper = Inf,
      rel.tol = 1e-10,
      subdivisions = 1000
    )$value
  }
  
  term1 + term2 - 2 * cross
}


# =============================================================================
# 2. Shared O(n log n) components
# =============================================================================

shared_components <- function(X, Y) {
  n1 <- length(X)
  n2 <- length(Y)
  
  if (n1 < 3L || n2 < 3L) {
    stop("Both samples must have size at least 3 for leave-one-out formulas.")
  }
  
  Xs <- sort(X)
  Ys <- sort(Y)
  
  cumX <- c(0, cumsum(Xs))
  cumY <- c(0, cumsum(Ys))
  
  sumXX <- sum(Xs * (n1 - seq_len(n1)))
  sumYY <- sum(Ys * (n2 - seq_len(n2)))
  
  # rowX[i] = sum_j min(X[i], Y[j])
  rowX <- vapply(X, function(xi) {
    k <- findInterval(xi, Ys)
    cumY[k + 1L] + xi * (n2 - k)
  }, numeric(1))
  
  # rowY[j] = sum_i min(X[i], Y[j])
  rowY <- vapply(Y, function(yj) {
    k <- findInterval(yj, Xs)
    cumX[k + 1L] + yj * (n1 - k)
  }, numeric(1))
  
  cross <- sum(rowX)
  
  UXX <- 2 * sumXX / (n1 * (n1 - 1))
  UYY <- 2 * sumYY / (n2 * (n2 - 1))
  UXY <- cross / (n1 * n2)
  
  # colX[i] = sum_{k != i} min(X[i], X[k])
  rX <- rank(X, ties.method = "first")
  rY <- rank(Y, ties.method = "first")
  
  colX <- cumX[rX] + X * (n1 - rX)
  colY <- cumY[rY] + Y * (n2 - rY)
  
  D_Ustat <- UXX + UYY - 2 * UXY
  
  # Empirical plug-in / V-statistic version
  D_Emp <- (n1 - 1) / n1 * UXX +
    (n2 - 1) / n2 * UYY -
    2 * UXY
  
  list(
    n1 = n1,
    n2 = n2,
    n = n1 + n2,
    sumXX = sumXX,
    sumYY = sumYY,
    cross = cross,
    rowX = rowX,
    rowY = rowY,
    colX = colX,
    colY = colY,
    UXX = UXX,
    UYY = UYY,
    UXY = UXY,
    D_Ustat = D_Ustat,
    D_Emp = D_Emp
  )
}


# =============================================================================
# 3. Leave-one-out estimates
# =============================================================================

loo_estimates <- function(sc) {
  n1 <- sc$n1
  n2 <- sc$n2
  
  # Remove X[i]
  UXX_loo_X <- 2 * (sc$sumXX - sc$colX) / ((n1 - 1) * (n1 - 2))
  UXY_loo_X <- (sc$cross - sc$rowX) / ((n1 - 1) * n2)
  
  D_Ustat_loo_X <- UXX_loo_X + sc$UYY - 2 * UXY_loo_X
  
  D_Emp_loo_X <- 2 * (sc$sumXX - sc$colX) / (n1 - 1)^2 +
    (n2 - 1) / n2 * sc$UYY -
    2 * UXY_loo_X
  
  # Remove Y[j]
  UYY_loo_Y <- 2 * (sc$sumYY - sc$colY) / ((n2 - 1) * (n2 - 2))
  UXY_loo_Y <- (sc$cross - sc$rowY) / (n1 * (n2 - 1))
  
  D_Ustat_loo_Y <- sc$UXX + UYY_loo_Y - 2 * UXY_loo_Y
  
  D_Emp_loo_Y <- (n1 - 1) / n1 * sc$UXX +
    2 * (sc$sumYY - sc$colY) / (n2 - 1)^2 -
    2 * UXY_loo_Y
  
  list(
    D_Ustat_loo_X = D_Ustat_loo_X,
    D_Ustat_loo_Y = D_Ustat_loo_Y,
    D_Emp_loo_X = D_Emp_loo_X,
    D_Emp_loo_Y = D_Emp_loo_Y
  )
}


# =============================================================================
# 4. Expected pseudo-values for JEL
# =============================================================================

# For a two-sample U-statistic with
# kernel degrees (m1, m2) = (2, 2), m = m1 + m2 = 4:
#   E[V_k] = theta * n/(n - m) * [(n2 - 1) m1/n1 - (m2 - 1)],  k = 1, ..., n1
#   E[V_k] = theta * n/(n - m) * [(n1 - 1) m2/n2 - (m1 - 1)],  k = n1 + 1, ..., n
# This holds for the JYZ pseudo-values built in jyz_pseudo_values() below.
EV_vec <- function(theta, n1, n2) {
  n <- n1 + n2
  
  cx <- (n / (n - 4)) * ((n2 - 1) * (2 / n1) - 1)
  cy <- (n / (n - 4)) * ((n1 - 1) * (2 / n2) - 1)
  
  c(rep(cx * theta, n1), rep(cy * theta, n2))
}



# the full-sample normalising constants kept fixed. Evaluating T at the n - 1
# remaining observations gives
#   T_{n-1}^{(-i)} = n/(n - 4) * (n1 - 2)/n1 * U^{(-i)},  if W_i is an X,
#   T_{n-1}^{(-i)} = n/(n - 4) * (n2 - 2)/n2 * U^{(-i)},  if W_i is a Y,
# with U^{(-i)} the two-sample U-statistic of the reduced samples.
# These pseudo-values satisfy mean(V) = D_Ustat exactly (eq. 3.6).
jyz_pseudo_values <- function(D_Ustat, loo_X, loo_Y, n1, n2) {
  n <- n1 + n2
  if (n <= 4L) {
    stop("JYZ pseudo-values require n1 + n2 > 4.")
  }
  
  T_loo_X <- (n / (n - 4)) * ((n1 - 2) / n1) * loo_X
  T_loo_Y <- (n / (n - 4)) * ((n2 - 2) / n2) * loo_Y
  
  c(n * D_Ustat - (n - 1) * T_loo_X,
    n * D_Ustat - (n - 1) * T_loo_Y)
}


neg2llr <- function(theta, V, n1, n2) {
  w <- V - EV_vec(theta, n1, n2)
  
  if (!all(is.finite(w))) {
    return(Inf)
  }
  
  # Empirical likelihood feasibility condition
  if (min(w) >= 0 || max(w) <= 0) {
    return(Inf)
  }
  
  tryCatch(
    emplik::el.test(w, mu = 0)$"-2LLR",
    error = function(e) Inf
  )
}


# =============================================================================
# 5. JEL confidence interval
# =============================================================================

jel_ci <- function(V, n1, n2, alpha = 0.05) {
  chi_crit <- stats::qchisq(1 - alpha, df = 1)
  D_hat <- mean(V)
  step <- stats::sd(V) / sqrt(n1 + n2)
  
  if (!is.finite(D_hat) || !is.finite(step) || step <= 0) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  # Upper endpoint
  u_hi <- max(D_hat + step, step)
  found_hi <- FALSE
  
  for (k in seq_len(80)) {
    val_hi <- neg2llr(u_hi, V, n1, n2)
    
    if (is.finite(val_hi) && val_hi >= chi_crit) {
      found_hi <- TRUE
      break
    }
    
    u_hi <- u_hi + step * 1.5^k
  }
  
  if (!found_hi) {
    ub <- NA_real_
  } else {
    ub <- tryCatch(
      stats::uniroot(
        function(t) neg2llr(t, V, n1, n2) - chi_crit,
        lower = D_hat,
        upper = u_hi,
        tol = 1e-5
      )$root,
      error = function(e) NA_real_
    )
  }
  
  # Lower endpoint
  # D is nonnegative. If 0 is inside the JEL confidence set, use 0.
  val_zero <- neg2llr(0, V, n1, n2)
  
  if (is.finite(val_zero) && val_zero <= chi_crit) {
    lb <- 0
  } else {
    u_lo <- max(D_hat - step, 0)
    found_lo <- FALSE
    
    for (k in seq_len(80)) {
      val_lo <- neg2llr(u_lo, V, n1, n2)
      
      if (u_lo <= 0 || (is.finite(val_lo) && val_lo >= chi_crit)) {
        found_lo <- TRUE
        break
      }
      
      u_lo <- max(u_lo - step * 1.5^k, 0)
    }
    
    if (!found_lo || u_lo >= D_hat) {
      lb <- 0
    } else {
      lb <- tryCatch(
        stats::uniroot(
          function(t) neg2llr(t, V, n1, n2) - chi_crit,
          lower = u_lo,
          upper = D_hat,
          tol = 1e-5
        )$root,
        error = function(e) 0
      )
    }
  }
  
  lb <- max(0, lb)
  
  if (!is.finite(ub) || ub < lb) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  c(lb = unname(lb), ub = unname(ub))
}


# =============================================================================
# 6. Confidence intervals for one sample pair
# =============================================================================

all_CIs <- function(
    X,
    Y,
    theta,
    alpha = 0.05
) {
  n1 <- length(X)
  n2 <- length(Y)
  n <- n1 + n2
  
  sc <- shared_components(X, Y)
  loo <- loo_estimates(sc)
  
  # U-statistic pseudo-values for JEL
  # (JYZ 2009 construction; see jyz_pseudo_values())
  V_JEL <- jyz_pseudo_values(
    sc$D_Ustat, loo$D_Ustat_loo_X, loo$D_Ustat_loo_Y, n1, n2
  )
  
  # U-statistic pseudo-values for normal approximation
  V_NA_X <- n1 * sc$D_Ustat - (n1 - 1) * loo$D_Ustat_loo_X
  V_NA_Y <- n2 * sc$D_Ustat - (n2 - 1) * loo$D_Ustat_loo_Y
  
  # Empirical estimator pseudo-values
  V_Emp_X <- n1 * sc$D_Emp - (n1 - 1) * loo$D_Emp_loo_X
  V_Emp_Y <- n2 * sc$D_Emp - (n2 - 1) * loo$D_Emp_loo_Y
  
  # CI 1: U-statistic + JEL
  ci_JEL <- jel_ci(V_JEL, n1, n2, alpha)
  
  if (all(is.finite(ci_JEL))) {
    lb <- unname(ci_JEL["lb"])
    ub <- unname(ci_JEL["ub"])
    
    cov_JEL <- theta >= lb && theta <= ub
    al_JEL <- ub - lb
  } else {
    cov_JEL <- NA
    al_JEL <- NA_real_
  }
  
  # CI 2: U-statistic + Normal approximation
  sigma2_NA <- stats::var(V_NA_X) / n1 + stats::var(V_NA_Y) / n2
  
  if (is.finite(sigma2_NA) && sigma2_NA >= 0) {
    hw_NA <- stats::qnorm(1 - alpha / 2) * sqrt(sigma2_NA)
    
    lb_NA <- max(0, sc$D_Ustat - hw_NA)
    ub_NA <- max(0, sc$D_Ustat + hw_NA)
    
    if (all(is.finite(c(lb_NA, ub_NA))) && ub_NA >= lb_NA) {
      cov_NA <- theta >= lb_NA && theta <= ub_NA
      al_NA <- ub_NA - lb_NA
    } else {
      cov_NA <- NA
      al_NA <- NA_real_
    }
  } else {
    cov_NA <- NA
    al_NA <- NA_real_
  }
  
  # CI 3: Empirical estimator + Normal approximation
  sigma2_Emp <- stats::var(V_Emp_X) / n1 + stats::var(V_Emp_Y) / n2
  
  if (is.finite(sigma2_Emp) && sigma2_Emp >= 0) {
    hw_Emp <- stats::qnorm(1 - alpha / 2) * sqrt(sigma2_Emp)
    
    lb_Emp <- max(0, sc$D_Emp - hw_Emp)
    ub_Emp <- max(0, sc$D_Emp + hw_Emp)
    
    if (all(is.finite(c(lb_Emp, ub_Emp))) && ub_Emp >= lb_Emp) {
      cov_Emp <- theta >= lb_Emp && theta <= ub_Emp
      al_Emp <- ub_Emp - lb_Emp
    } else {
      cov_Emp <- NA
      al_Emp <- NA_real_
    }
  } else {
    cov_Emp <- NA
    al_Emp <- NA_real_
  }
  
  list(
    cov = c(
      JEL = unname(cov_JEL),
      NA_ = unname(cov_NA),
      Emp = unname(cov_Emp)
    ),
    al = c(
      JEL = unname(al_JEL),
      NA_ = unname(al_NA),
      Emp = unname(al_Emp)
    )
  )
}


# =============================================================================
# 7. Functions
# =============================================================================

safe_mean <- function(z) {
  z <- z[is.finite(z)]
  
  if (length(z) == 0L) {
    return(NA_real_)
  }
  
  mean(z)
}


safe_cp <- function(success, planned) {
  success <- unname(success)
  planned <- unname(planned)
  
  if (!is.finite(planned) || planned <= 0) {
    return(NA_real_)
  }
  
  100 * success / planned
}


# =============================================================================
# 8. Simulation for Weibull distributions
# =============================================================================

run_coverage_simulation <- function(
    iterations = 2000,
    alpha = 0.05,
    seed = 2024,
    verbose = TRUE
) {
  RNGkind(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  set.seed(seed)
  
  # Parameter settings:
  # c(shape1, scale1, shape2, scale2)
  #
  # These include both shape and scale differences.
  # Do not use only scale differences; that would be a weak Weibull study.
  params_list <- list(
    c(0.5, 1.0, 1.0, 1.0),
    c(1.0, 1.0, 2.0, 1.0),
    c(2.0, 1.0, 1.0, 1.0),
    c(1.5, 1.0, 3.0, 1.5),
    c(0.7, 2.0, 1.5, 1.0),
    c(3.0, 1.0, 1.2, 2.0)
  )
  
  n1_vec <- c(10, 30, 70, 100)
  n2_vec <- c(10, 40, 50, 100)
  
  out <- vector("list", length(params_list) * length(n1_vec))
  idx <- 1L
  
  sep <- strrep("-", 110)
  
  for (params in params_list) {
    shape1 <- params[1]
    scale1 <- params[2]
    shape2 <- params[3]
    scale2 <- params[4]
    
    theta <- true_D_weibull(shape1, scale1, shape2, scale2)
    
    if (verbose) {
      cat(sprintf(
        "\n%s\nWeibull 1: shape = %.3f, scale = %.3f | Weibull 2: shape = %.3f, scale = %.3f | true D = %.8f\n%s\n",
        sep, shape1, scale1, shape2, scale2, theta, sep
      ))
      
      cat(sprintf(
        "%-12s  %8s %8s   %8s %8s   %8s %8s\n",
        "(n1,n2)",
        "JEL_CP", "JEL_AL",
        "NA_CP", "NA_AL",
        "Emp_CP", "Emp_AL"
      ))
    }
    
    for (pp in seq_along(n1_vec)) {
      n1 <- n1_vec[pp]
      n2 <- n2_vec[pp]
      
      acc_cov <- c(JEL = 0, NA_ = 0, Emp = 0)
      valid_cov <- c(JEL = 0, NA_ = 0, Emp = 0)
      
      acc_al <- list(
        JEL = rep(NA_real_, iterations),
        NA_ = rep(NA_real_, iterations),
        Emp = rep(NA_real_, iterations)
      )
      
      for (h in seq_len(iterations)) {
        X <- stats::rweibull(n1, shape = shape1, scale = scale1)
        Y <- stats::rweibull(n2, shape = shape2, scale = scale2)
        
        res <- tryCatch(
          all_CIs(
            X = X,
            Y = Y,
            theta = theta,
            alpha = alpha
          ),
          error = function(e) {
            list(
              cov = c(JEL = NA, NA_ = NA, Emp = NA),
              al = c(JEL = NA_real_, NA_ = NA_real_, Emp = NA_real_)
            )
          }
        )
        
        for (method in names(acc_cov)) {
          if (!is.na(res$cov[method])) {
            acc_cov[method] <- acc_cov[method] + as.numeric(res$cov[method])
            valid_cov[method] <- valid_cov[method] + 1
          }
          
          acc_al[[method]][h] <- unname(res$al[method])
        }
      }
      
      # Coverage uses all planned replications; failed intervals are noncoverage.
      # Average length is calculated over the intervals that were valid.
      cp <- c(
        JEL = unname(safe_cp(acc_cov["JEL"], iterations)),
        NA_ = unname(safe_cp(acc_cov["NA_"], iterations)),
        Emp = unname(safe_cp(acc_cov["Emp"], iterations))
      )
      
      al <- c(
        JEL = unname(safe_mean(acc_al$JEL)),
        NA_ = unname(safe_mean(acc_al$NA_)),
        Emp = unname(safe_mean(acc_al$Emp))
      )
      
      row <- data.frame(
        shape1 = shape1,
        scale1 = scale1,
        shape2 = shape2,
        scale2 = scale2,
        true_D = theta,
        n1 = n1,
        n2 = n2,
        
        JEL_CP = unname(cp["JEL"]),
        JEL_AL = unname(al["JEL"]),
        
        NA_CP = unname(cp["NA_"]),
        NA_AL = unname(al["NA_"]),
        
        Emp_CP = unname(cp["Emp"]),
        Emp_AL = unname(al["Emp"]),
        
        Covered_JEL = unname(acc_cov["JEL"]),
        Covered_NA = unname(acc_cov["NA_"]),
        Covered_Emp = unname(acc_cov["Emp"]),
        Valid_JEL = unname(valid_cov["JEL"]),
        Valid_NA = unname(valid_cov["NA_"]),
        Valid_Emp = unname(valid_cov["Emp"]),
        Failed_JEL = iterations - unname(valid_cov["JEL"]),
        Failed_NA = iterations - unname(valid_cov["NA_"]),
        Failed_Emp = iterations - unname(valid_cov["Emp"]),
        
        row.names = NULL
      )
      
      out[[idx]] <- row
      idx <- idx + 1L
      
      if (verbose) {
        cat(sprintf(
          "(%3d,%3d)     %8.2f %8.4f   %8.2f %8.4f   %8.2f %8.4f\n",
          n1, n2,
          cp["JEL"], al["JEL"],
          cp["NA_"], al["NA_"],
          cp["Emp"], al["Emp"]
        ))
      }
    }
  }
  
  results <- do.call(rbind, out)
  rownames(results) <- NULL
  
  results
}


# =============================================================================
# 9. Print results
# =============================================================================

print_results_by_parameter <- function(results, digits_cp = 2, digits_al = 4) {
  param_sets <- unique(results[, c("shape1", "scale1", "shape2", "scale2")])
  
  for (i in seq_len(nrow(param_sets))) {
    shape1 <- param_sets$shape1[i]
    scale1 <- param_sets$scale1[i]
    shape2 <- param_sets$shape2[i]
    scale2 <- param_sets$scale2[i]
    
    sub <- results[
      results$shape1 == shape1 &
        results$scale1 == scale1 &
        results$shape2 == shape2 &
        results$scale2 == scale2,
    ]
    
    cat("\n", strrep("-", 110), "\n", sep = "")
    cat(sprintf(
      "Weibull 1: shape = %.3f, scale = %.3f | Weibull 2: shape = %.3f, scale = %.3f | true D = %.8f\n",
      shape1, scale1, shape2, scale2, sub$true_D[1]
    ))
    cat(strrep("-", 110), "\n", sep = "")
    
    tab <- data.frame(
      n1 = sub$n1,
      n2 = sub$n2,
      
      JEL_CP = round(sub$JEL_CP, digits_cp),
      JEL_AL = round(sub$JEL_AL, digits_al),
      
      NA_CP = round(sub$NA_CP, digits_cp),
      NA_AL = round(sub$NA_AL, digits_al),
      
      Emp_CP = round(sub$Emp_CP, digits_cp),
      Emp_AL = round(sub$Emp_AL, digits_al),
      
      Covered_JEL = sub$Covered_JEL,
      Covered_NA = sub$Covered_NA,
      Covered_Emp = sub$Covered_Emp,
      Valid_JEL = sub$Valid_JEL,
      Valid_NA = sub$Valid_NA,
      Valid_Emp = sub$Valid_Emp,
      Failed_JEL = sub$Failed_JEL,
      Failed_NA = sub$Failed_NA,
      Failed_Emp = sub$Failed_Emp
    )
    
    print(tab, row.names = FALSE)
  }
}


# =============================================================================
# 10. Run Weibull simulation
# =============================================================================

read_env_integer <- function(name, default, minimum = 0L) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) {
    return(as.integer(default))
  }

  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || !is.finite(value) || value != floor(value) ||
      value < minimum || value > .Machine$integer.max) {
    stop(sprintf("%s must be an integer of at least %d.", name, minimum))
  }

  as.integer(value)
}


read_env_flag <- function(name, default = TRUE) {
  raw <- tolower(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(raw)) {
    return(default)
  }
  if (raw %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (raw %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop(sprintf("%s must be TRUE or FALSE.", name))
}


mc_reps <- read_env_integer("EXTROPY_MC_REPS", 2000L, minimum = 1L)
simulation_seed <- read_env_integer("EXTROPY_SEED", 2026L)
output_dir <- Sys.getenv("EXTROPY_OUTPUT_DIR", unset = ".")
verbose <- read_env_flag("EXTROPY_VERBOSE", TRUE)

if (!nzchar(output_dir)) {
  stop("EXTROPY_OUTPUT_DIR must not be empty.")
}
if (!dir.exists(output_dir) &&
    !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop(sprintf("Could not create output directory: %s", output_dir))
}

results_section_5_3_weibull <- run_coverage_simulation(
  iterations = mc_reps,
  alpha = 0.05,
  seed = simulation_seed,
  verbose = verbose
)

if (verbose) {
  print_results_by_parameter(results_section_5_3_weibull)
}

write.csv(
  results_section_5_3_weibull,
  file.path(
    output_dir,
    "Coverage_and_Average_Length_Section_5_3_Weibull.csv"
  ),
  row.names = FALSE
)
