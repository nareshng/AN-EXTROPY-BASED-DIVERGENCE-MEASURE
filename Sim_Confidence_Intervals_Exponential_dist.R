# =============================================================================
# Section 5.3 – Coverage Probability and Average Length
# for FOUR confidence interval methods under exponential distributions:
#
#   1. U-statistic + JEL
#   2. U-statistic + Normal approximation
#   3. Empirical estimator + Normal approximation
#   4. Kernel estimator + Bootstrap
#
# Model:
#   X ~ Exp(lambda1)
#   Y ~ Exp(lambda2)
#
#
#
# True value:
#   D = 1/(2 lambda1) + 1/(2 lambda2) - 2/(lambda1 + lambda2)
#
# 
# =============================================================================


# =============================================================================
# Package
# =============================================================================

if (!requireNamespace("emplik", quietly = TRUE)) {
  stop("Package 'emplik' is required. Install it using install.packages('emplik').")
}

library(emplik)


# =============================================================================
# 1. True divergence for exponential distributions
# =============================================================================

true_D_exp <- function(lambda1, lambda2) {
  if (lambda1 <= 0 || lambda2 <= 0) {
    stop("lambda1 and lambda2 must be positive.")
  }
  
  1 / (2 * lambda1) + 1 / (2 * lambda2) - 2 / (lambda1 + lambda2)
}


# =============================================================================
# 2. Bandwidth
# =============================================================================

silverman_bw <- function(z) {
  n <- length(z)
  
  if (n < 2L) {
    return(NA_real_)
  }
  
  h <- 0.9 * min(stats::sd(z), stats::IQR(z) / 1.34) * n^(-0.2)
  
  if (!is.finite(h) || h <= 0) {
    return(NA_real_)
  }
  
  h
}


# =============================================================================
# 3. Trapezoidal integration
# =============================================================================

trapz <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }
  
  if (length(x) < 2L) {
    return(NA_real_)
  }
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}


# =============================================================================
# 4. components
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
# 5. Leave-one-out estimates
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
# 6. Expected pseudo-values for JEL
# =============================================================================

EV_vec <- function(theta, n1, n2) {
  n <- n1 + n2
  
  cx <- (n / (n - 2)) * ((n2 - 1) * (2 / n1) - 1)
  cy <- (n / (n - 2)) * ((n1 - 1) * (2 / n2) - 1)
  
  c(rep(cx * theta, n1), rep(cy * theta, n2))
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
# 7. JEL confidence interval
# =============================================================================

jel_ci <- function(V, n1, n2, alpha = 0.05) {
  chi_crit <- stats::qchisq(1 - alpha, df = 1)
  D_hat <- mean(V)
  step <- stats::sd(V) / sqrt(n1 + n2)
  
  if (!is.finite(D_hat) || !is.finite(step) || step <= 0) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  
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
# 8. kernel estimator
# =============================================================================

make_kernel_object <- function(
    X,
    Y,
    grid_size = 400,
    tail_mult = 8
) {
  n1 <- length(X)
  n2 <- length(Y)
  
  h1 <- silverman_bw(X)
  h2 <- silverman_bw(Y)
  
  if (!is.finite(h1) || h1 <= 0 || !is.finite(h2) || h2 <= 0) {
    return(NULL)
  }
  
  upper <- max(c(X, Y)) + tail_mult * max(h1, h2)
  
  if (!is.finite(upper) || upper <= 0) {
    return(NULL)
  }
  
  xg <- seq(0, upper, length.out = grid_size)
  
  PhiX <- stats::pnorm(outer(X, xg, function(a, b) (a - b) / h1))
  PhiY <- stats::pnorm(outer(Y, xg, function(a, b) (a - b) / h2))
  
  list(
    n1 = n1,
    n2 = n2,
    xg = xg,
    PhiX = PhiX,
    PhiY = PhiY
  )
}


kernel_D_from_object <- function(obj, idxX = NULL, idxY = NULL) {
  if (is.null(obj)) {
    return(NA_real_)
  }
  
  if (is.null(idxX)) {
    idxX <- seq_len(obj$n1)
  }
  
  if (is.null(idxY)) {
    idxY <- seq_len(obj$n2)
  }
  
  Fb <- colMeans(obj$PhiX[idxX, , drop = FALSE])
  Gb <- colMeans(obj$PhiY[idxY, , drop = FALSE])
  
  trapz(obj$xg, (Fb - Gb)^2)
}


# =============================================================================
# 9. kernel bootstrap CI
# =============================================================================

kernel_bootstrap_ci_fast <- function(
    X,
    Y,
    B = 199,
    alpha = 0.05,
    grid_size = 400,
    tail_mult = 8,
    min_valid_frac = 0.8
) {
  obj <- make_kernel_object(
    X = X,
    Y = Y,
    grid_size = grid_size,
    tail_mult = tail_mult
  )
  
  if (is.null(obj)) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  D_hat <- kernel_D_from_object(obj)
  
  if (!is.finite(D_hat)) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  boot_D <- numeric(B)
  
  for (b in seq_len(B)) {
    idxX <- sample.int(obj$n1, obj$n1, replace = TRUE)
    idxY <- sample.int(obj$n2, obj$n2, replace = TRUE)
    
    boot_D[b] <- kernel_D_from_object(obj, idxX, idxY)
  }
  
  boot_D <- boot_D[is.finite(boot_D)]
  
  if (length(boot_D) < min_valid_frac * B) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  Q_lo <- as.numeric(stats::quantile(
    boot_D,
    probs = alpha / 2,
    names = FALSE,
    type = 8
  ))
  
  Q_hi <- as.numeric(stats::quantile(
    boot_D,
    probs = 1 - alpha / 2,
    names = FALSE,
    type = 8
  ))
  
  if (!is.finite(Q_lo) || !is.finite(Q_hi)) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  lb <- max(0, 2 * D_hat - Q_hi)
  ub <- max(0, 2 * D_hat - Q_lo)
  
  if (!is.finite(lb) || !is.finite(ub) || ub < lb) {
    return(c(lb = NA_real_, ub = NA_real_))
  }
  
  c(lb = unname(lb), ub = unname(ub))
}


# =============================================================================
# 10. Confidence intervals for one sample
# =============================================================================

all_CIs <- function(
    X,
    Y,
    theta,
    B_boot = 199,
    alpha = 0.05,
    grid_size = 400,
    tail_mult = 8
) {
  n1 <- length(X)
  n2 <- length(Y)
  n <- n1 + n2
  
  sc <- shared_components(X, Y)
  loo <- loo_estimates(sc)
  
  
  # Pseudo-values
  
  # U-statistic pseudo-values for JEL
  V_JEL_X <- n * sc$D_Ustat - (n - 1) * loo$D_Ustat_loo_X
  V_JEL_Y <- n * sc$D_Ustat - (n - 1) * loo$D_Ustat_loo_Y
  V_JEL <- c(V_JEL_X, V_JEL_Y)
  
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
  
  # ---------------------------------------------------------------------------
  # CI 2: U-statistic + Normal approximation
  # ---------------------------------------------------------------------------
  
  sigma2_NA <- stats::var(V_NA_X) / n1 + stats::var(V_NA_Y) / n2
  
  if (is.finite(sigma2_NA) && sigma2_NA >= 0) {
    hw_NA <- stats::qnorm(1 - alpha / 2) * sqrt(sigma2_NA)
    
    lb_NA <- max(0, sc$D_Ustat - hw_NA)
    ub_NA <- sc$D_Ustat + hw_NA
    
    cov_NA <- theta >= lb_NA && theta <= ub_NA
    al_NA <- ub_NA - lb_NA
  } else {
    cov_NA <- NA
    al_NA <- NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # CI 3: Empirical estimator + Normal approximation
  # ---------------------------------------------------------------------------
  
  sigma2_Emp <- stats::var(V_Emp_X) / n1 + stats::var(V_Emp_Y) / n2
  
  if (is.finite(sigma2_Emp) && sigma2_Emp >= 0) {
    hw_Emp <- stats::qnorm(1 - alpha / 2) * sqrt(sigma2_Emp)
    
    lb_Emp <- max(0, sc$D_Emp - hw_Emp)
    ub_Emp <- sc$D_Emp + hw_Emp
    
    cov_Emp <- theta >= lb_Emp && theta <= ub_Emp
    al_Emp <- ub_Emp - lb_Emp
  } else {
    cov_Emp <- NA
    al_Emp <- NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # CI 4: Kernel estimator + bootstrap
  # ---------------------------------------------------------------------------
  
  ci_Ker <- kernel_bootstrap_ci_fast(
    X = X,
    Y = Y,
    B = B_boot,
    alpha = alpha,
    grid_size = grid_size,
    tail_mult = tail_mult
  )
  
  if (all(is.finite(ci_Ker))) {
    lb <- unname(ci_Ker["lb"])
    ub <- unname(ci_Ker["ub"])
    
    cov_Ker <- theta >= lb && theta <= ub
    al_Ker <- ub - lb
  } else {
    cov_Ker <- NA
    al_Ker <- NA_real_
  }
  
  list(
    cov = c(
      JEL = unname(cov_JEL),
      NA_ = unname(cov_NA),
      Emp = unname(cov_Emp),
      Ker = unname(cov_Ker)
    ),
    al = c(
      JEL = unname(al_JEL),
      NA_ = unname(al_NA),
      Emp = unname(al_Emp),
      Ker = unname(al_Ker)
    )
  )
}


# =============================================================================
# 11. Functions
# =============================================================================

safe_mean <- function(z) {
  z <- z[is.finite(z)]
  
  if (length(z) == 0L) {
    return(NA_real_)
  }
  
  mean(z)
}


safe_cp <- function(success, valid) {
  success <- unname(success)
  valid <- unname(valid)
  
  if (!is.finite(valid) || valid <= 0) {
    return(NA_real_)
  }
  
  100 * success / valid
}


# =============================================================================
# 12. Simulation
# =============================================================================

run_coverage_simulation <- function(
    iterations = 1000,
    B_boot = 199,
    grid_size = 400,
    tail_mult = 8,
    alpha = 0.05,
    seed = 2024,
    verbose = TRUE
) {
  set.seed(seed)
  
  params_list <- list(
    c(0.2, 0.1),
    c(0.1, 0.5),
    c(1.0, 0.5),
    c(1.0, 2.0),
    c(10.0, 1.0)
  )
  
  n1_vec <- c(20, 30, 70, 100)
  n2_vec <- c(10, 40, 50, 100)
  
  out <- vector("list", length(params_list) * length(n1_vec))
  idx <- 1L
  
  sep <- strrep("-", 100)
  
  for (params in params_list) {
    lambda1 <- params[1]
    lambda2 <- params[2]
    theta <- true_D_exp(lambda1, lambda2)
    
    if (verbose) {
      cat(sprintf(
        "\n%s\n(lambda1, lambda2) = (%.3f, %.3f), true D = %.8f\n%s\n",
        sep, lambda1, lambda2, theta, sep
      ))
      
      cat(sprintf(
        "%-12s  %8s %8s   %8s %8s   %8s %8s   %8s %8s\n",
        "(n1,n2)",
        "JEL_CP", "JEL_AL",
        "NA_CP", "NA_AL",
        "Emp_CP", "Emp_AL",
        "Ker_CP", "Ker_AL"
      ))
    }
    
    for (pp in seq_along(n1_vec)) {
      n1 <- n1_vec[pp]
      n2 <- n2_vec[pp]
      
      acc_cov <- c(JEL = 0, NA_ = 0, Emp = 0, Ker = 0)
      valid_cov <- c(JEL = 0, NA_ = 0, Emp = 0, Ker = 0)
      
      acc_al <- list(
        JEL = rep(NA_real_, iterations),
        NA_ = rep(NA_real_, iterations),
        Emp = rep(NA_real_, iterations),
        Ker = rep(NA_real_, iterations)
      )
      
      for (h in seq_len(iterations)) {
        X <- stats::rexp(n1, rate = lambda1)
        Y <- stats::rexp(n2, rate = lambda2)
        
        res <- tryCatch(
          all_CIs(
            X = X,
            Y = Y,
            theta = theta,
            B_boot = B_boot,
            alpha = alpha,
            grid_size = grid_size,
            tail_mult = tail_mult
          ),
          error = function(e) {
            list(
              cov = c(JEL = NA, NA_ = NA, Emp = NA, Ker = NA),
              al = c(JEL = NA_real_, NA_ = NA_real_, Emp = NA_real_, Ker = NA_real_)
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
      
      cp <- c(
        JEL = unname(safe_cp(acc_cov["JEL"], valid_cov["JEL"])),
        NA_ = unname(safe_cp(acc_cov["NA_"], valid_cov["NA_"])),
        Emp = unname(safe_cp(acc_cov["Emp"], valid_cov["Emp"])),
        Ker = unname(safe_cp(acc_cov["Ker"], valid_cov["Ker"]))
      )
      
      al <- c(
        JEL = unname(safe_mean(acc_al$JEL)),
        NA_ = unname(safe_mean(acc_al$NA_)),
        Emp = unname(safe_mean(acc_al$Emp)),
        Ker = unname(safe_mean(acc_al$Ker))
      )
      
      row <- data.frame(
        lambda1 = lambda1,
        lambda2 = lambda2,
        true_D = theta,
        n1 = n1,
        n2 = n2,
        
        JEL_CP = unname(cp["JEL"]),
        JEL_AL = unname(al["JEL"]),
        
        NA_CP = unname(cp["NA_"]),
        NA_AL = unname(al["NA_"]),
        
        Emp_CP = unname(cp["Emp"]),
        Emp_AL = unname(al["Emp"]),
        
        Ker_CP = unname(cp["Ker"]),
        Ker_AL = unname(al["Ker"]),
        
        Valid_JEL = unname(valid_cov["JEL"]),
        Valid_NA = unname(valid_cov["NA_"]),
        Valid_Emp = unname(valid_cov["Emp"]),
        Valid_Ker = unname(valid_cov["Ker"]),
        
        row.names = NULL
      )
      
      out[[idx]] <- row
      idx <- idx + 1L
      
      if (verbose) {
        cat(sprintf(
          "(%3d,%3d)     %8.2f %8.4f   %8.2f %8.4f   %8.2f %8.4f   %8.2f %8.4f\n",
          n1, n2,
          cp["JEL"], al["JEL"],
          cp["NA_"], al["NA_"],
          cp["Emp"], al["Emp"],
          cp["Ker"], al["Ker"]
        ))
      }
    }
  }
  
  results <- do.call(rbind, out)
  rownames(results) <- NULL
  
  results
}


# =============================================================================
# 13. Print results
# =============================================================================

print_results_by_parameter <- function(results, digits_cp = 2, digits_al = 4) {
  param_pairs <- unique(results[, c("lambda1", "lambda2")])
  
  for (i in seq_len(nrow(param_pairs))) {
    lambda1 <- param_pairs$lambda1[i]
    lambda2 <- param_pairs$lambda2[i]
    
    sub <- results[
      results$lambda1 == lambda1 &
        results$lambda2 == lambda2,
    ]
    
    cat("\n", strrep("-", 100), "\n", sep = "")
    cat(sprintf(
      "(lambda1, lambda2) = (%.3f, %.3f), true D = %.8f\n",
      lambda1, lambda2, sub$true_D[1]
    ))
    cat(strrep("-", 100), "\n", sep = "")
    
    tab <- data.frame(
      n1 = sub$n1,
      n2 = sub$n2,
      
      JEL_CP = round(sub$JEL_CP, digits_cp),
      JEL_AL = round(sub$JEL_AL, digits_al),
      
      NA_CP = round(sub$NA_CP, digits_cp),
      NA_AL = round(sub$NA_AL, digits_al),
      
      Emp_CP = round(sub$Emp_CP, digits_cp),
      Emp_AL = round(sub$Emp_AL, digits_al),
      
      Ker_CP = round(sub$Ker_CP, digits_cp),
      Ker_AL = round(sub$Ker_AL, digits_al),
      
      Valid_JEL = sub$Valid_JEL,
      Valid_NA = sub$Valid_NA,
      Valid_Emp = sub$Valid_Emp,
      Valid_Ker = sub$Valid_Ker
    )
    
    print(tab, row.names = FALSE)
  }
}


# =============================================================================
# 14. Run simulation
# =============================================================================

results_section_5_3_exp <- run_coverage_simulation(
  iterations = 2000,
  B_boot = 499,
  grid_size = 400,
  tail_mult = 8,
  alpha = 0.05,
  seed = 2026,
  verbose = TRUE
)

print_results_by_parameter(results_section_5_3_exp)

write.csv(
  results_section_5_3_exp,
  "Coverage_and_Average_Length_Section_5_3_Exponential_Efficient_Corrected.csv",
  row.names = FALSE
)



