

# =============================================================================
# divergence for Weibull distributions
#
# For X ~ Weibull(shape1, scale1), Y ~ Weibull(shape2, scale2),
#
# # MSE and Relative MSE comparison of Kernel, Empirical, and U-statistic estimators
#
# D(F,G) = integral_0^inf [Fbar(x) - Gbar(x)]^2 dx
#
# Relative MSE is computed using U-statistic as baseline:
# RelMSE_Estimator = MSE_Estimator / MSE_Ustat
#
# Bandwidth is kept exactly as in the original code:
# h = 0.9 * min(sd(X), IQR(X) / 1.34) * n^(-1/5)
# =============================================================================


# =============================================================================
# 1. Gauss-Legendre quadrature
# =============================================================================

gauss_legendre <- function(n) {
  if (n < 2L) stop("n must be at least 2.")
  
  i <- seq_len(n - 1L)
  beta <- i / sqrt(4 * i^2 - 1)
  
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- beta
  J[cbind(i + 1L, i)] <- beta
  
  eig <- eigen(J, symmetric = TRUE)
  
  nodes <- eig$values
  weights <- 2 * eig$vectors[1, ]^2
  
  ord <- order(nodes)
  
  list(
    x = nodes[ord],
    w = weights[ord]
  )
}


# =============================================================================
# 2. Bandwidth rule
# =============================================================================

bw_original <- function(z) {
  n <- length(z)
  0.9 * min(stats::sd(z), stats::IQR(z) / 1.34) * n^(-1/5)
}


# =============================================================================
# 3. U-statistic estimator
# =============================================================================

calc_Ustat <- function(X, Y) {
  n1 <- length(X)
  n2 <- length(Y)
  
  if (n1 < 2L || n2 < 2L) {
    stop("Both samples must have size at least 2.")
  }
  
  Xs <- sort(X)
  Ys <- sort(Y)
  
  # E[min(X1, X2)]
  UXX <- 2 * sum(Xs * (n1 - seq_len(n1))) / (n1 * (n1 - 1))
  
  # E[min(Y1, Y2)]
  UYY <- 2 * sum(Ys * (n2 - seq_len(n2))) / (n2 * (n2 - 1))
  
  # E[min(X, Y)]
  csy <- c(0, cumsum(Ys))
  
  k <- findInterval(X, Ys)
  cross <- sum(csy[k + 1L] + X * (n2 - k))
  
  UXY <- cross / (n1 * n2)
  
  UXX + UYY - 2 * UXY
}


# =============================================================================
# 4. Empirical estimator
# =============================================================================

calc_Emp <- function(X, Y) {
  n1 <- length(X)
  n2 <- length(Y)
  
  if (n1 < 1L || n2 < 1L) {
    stop("Both samples must be nonempty.")
  }
  
  a <- 1 / n1 + 1 / n2
  
  combined <- c(X, Y)
  lab <- c(rep(1L, n1), rep(2L, n2))
  
  pooled_rank <- rank(combined, ties.method = "average")
  
  S <- pooled_rank[lab == 1L][order(X)]
  R <- pooled_rank[lab == 2L][order(Y)]
  
  Xs <- sort(X)
  Ys <- sort(Y)
  
  term_X <- (2 / n1) * sum(Xs * (S / n2 - seq_len(n1) * a))
  term_Y <- (2 / n2) * sum(Ys * (R / n1 - seq_len(n2) * a))
  
  term_X + term_Y
}


# =============================================================================
# 5. Kernel estimator using Gauss-Legendre quadrature
# =============================================================================

calc_Kernel_GL <- function(X, Y, gl, upper = NULL) {
  n1 <- length(X)
  n2 <- length(Y)
  
  h1 <- bw_original(X)
  h2 <- bw_original(Y)
  
  if (!is.finite(h1) || h1 <= 0 || !is.finite(h2) || h2 <= 0) {
    return(NA_real_)
  }
  
  if (is.null(upper)) {
    upper <- max(c(X, Y)) + 8 * max(h1, h2)
  }
  
  if (!is.finite(upper) || upper <= 0) {
    return(NA_real_)
  }
  
  # Transform Gauss-Legendre nodes from [-1, 1] to [0, upper]
  xg <- 0.5 * upper * (gl$x + 1)
  wg <- 0.5 * upper * gl$w
  
  Fbar <- colMeans(pnorm(outer(X, xg, function(xi, x) (xi - x) / h1)))
  Gbar <- colMeans(pnorm(outer(Y, xg, function(yi, x) (yi - x) / h2)))
  
  sum(wg * (Fbar - Gbar)^2)
}


# =============================================================================
# 6. True divergence for Weibull distributions
# =============================================================================

true_D_weibull <- function(shape1, scale1, shape2, scale2) {
  if (shape1 <= 0 || shape2 <= 0 || scale1 <= 0 || scale2 <= 0) {
    stop("All Weibull shape and scale parameters must be positive.")
  }
  
  # Integral of S_X(x)^2 dx
  term1 <- scale1 * gamma(1 + 1 / shape1) / (2^(1 / shape1))
  
  # Integral of S_Y(x)^2 dx
  term2 <- scale2 * gamma(1 + 1 / shape2) / (2^(1 / shape2))
  
  # Integral of S_X(x) S_Y(x) dx
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
# 7. One simulation cell
# =============================================================================

simulate_one_cell <- function(shape1, scale1, shape2, scale2,
                              n1, n2, iterations, gl) {
  
  true_D <- true_D_weibull(shape1, scale1, shape2, scale2)
  
  ker_est <- numeric(iterations)
  emp_est <- numeric(iterations)
  ust_est <- numeric(iterations)
  
  # Deterministic high quantile for tail coverage in kernel integration.
  q_upper <- max(
    qweibull(1 - 1e-12, shape = shape1, scale = scale1),
    qweibull(1 - 1e-12, shape = shape2, scale = scale2)
  )
  
  for (b in seq_len(iterations)) {
    X <- rweibull(n1, shape = shape1, scale = scale1)
    Y <- rweibull(n2, shape = shape2, scale = scale2)
    
    h1 <- bw_original(X)
    h2 <- bw_original(Y)
    
    upper_b <- max(
      q_upper,
      max(c(X, Y)) + 8 * max(h1, h2),
      na.rm = TRUE
    )
    
    ker_est[b] <- calc_Kernel_GL(X, Y, gl, upper = upper_b)
    emp_est[b] <- calc_Emp(X, Y)
    ust_est[b] <- calc_Ustat(X, Y)
  }
  
  valid_ker <- is.finite(ker_est)
  
  MSE_Kernel <- mean((ker_est[valid_ker] - true_D)^2)
  MSE_Emp    <- mean((emp_est - true_D)^2)
  MSE_Ustat  <- mean((ust_est - true_D)^2)
  
  RelMSE_Kernel <- MSE_Kernel / MSE_Ustat
  RelMSE_Emp    <- MSE_Emp / MSE_Ustat
  RelMSE_Ustat  <- 1
  
  c(
    MSE_Kernel = MSE_Kernel,
    MSE_Emp = MSE_Emp,
    MSE_Ustat = MSE_Ustat,
    RelMSE_Kernel = RelMSE_Kernel,
    RelMSE_Emp = RelMSE_Emp,
    RelMSE_Ustat = RelMSE_Ustat,
    Kernel_NA_Count = sum(!valid_ker)
  )
}


# =============================================================================
# 8.full Weibull simulation
# =============================================================================

run_simulation <- function(
    iterations = 2000,
    n_quad = 80,
    seed = 2024,
    verbose = TRUE
) {
  set.seed(seed)
  
  gl <- gauss_legendre(n_quad)
  
  # Parameter grid:
  # c(shape1, scale1, shape2, scale2)
  #
  # This grid includes both shape differences and scale differences.
  # That is important. A Weibull simulation that changes only scale is weak.
  param_grid <- list(
    c(0.5, 1.0, 1.0, 1.0),
    c(1.0, 1.0, 2.0, 1.0),
    c(1.5, 1.0, 3.0, 1.0),
    c(2.0, 1.0, 2.0, 2.0),
    c(0.7, 2.0, 1.5, 1.0),
    c(3.0, 1.0, 1.2, 2.0)
  )
  
  n1_vec <- c(5,  10, 20, 20, 40, 30, 50)
  n2_vec <- c(10, 10, 10, 20, 30, 40, 50)
  
  out <- vector("list", length(param_grid) * length(n1_vec))
  idx <- 1L
  
  for (p in seq_along(param_grid)) {
    shape1 <- param_grid[[p]][1]
    scale1 <- param_grid[[p]][2]
    shape2 <- param_grid[[p]][3]
    scale2 <- param_grid[[p]][4]
    
    true_D <- true_D_weibull(shape1, scale1, shape2, scale2)
    
    if (verbose) {
      cat("\n============================================================\n")
      cat(sprintf(
        "Weibull 1: shape = %.3f, scale = %.3f | Weibull 2: shape = %.3f, scale = %.3f | true D = %.8f\n",
        shape1, scale1, shape2, scale2, true_D
      ))
      cat("============================================================\n")
    }
    
    for (j in seq_along(n1_vec)) {
      n1 <- n1_vec[j]
      n2 <- n2_vec[j]
      
      if (verbose) {
        cat(sprintf("Running n1 = %d, n2 = %d ...\n", n1, n2))
      }
      
      sim <- simulate_one_cell(
        shape1 = shape1,
        scale1 = scale1,
        shape2 = shape2,
        scale2 = scale2,
        n1 = n1,
        n2 = n2,
        iterations = iterations,
        gl = gl
      )
      
      out[[idx]] <- data.frame(
        shape1 = shape1,
        scale1 = scale1,
        shape2 = shape2,
        scale2 = scale2,
        true_D = true_D,
        n1 = n1,
        n2 = n2,
        
        MSE_Kernel = sim["MSE_Kernel"],
        MSE_Emp = sim["MSE_Emp"],
        MSE_Ustat = sim["MSE_Ustat"],
        
        RelMSE_Kernel = sim["RelMSE_Kernel"],
        RelMSE_Emp = sim["RelMSE_Emp"],
        RelMSE_Ustat = sim["RelMSE_Ustat"],
        
        Kernel_NA_Count = sim["Kernel_NA_Count"],
        
        row.names = NULL
      )
      
      idx <- idx + 1L
    }
  }
  
  results <- do.call(rbind, out)
  rownames(results) <- NULL
  
  results
}


# =============================================================================
# 9. tables
# =============================================================================

print_results_by_parameter <- function(results, digits_mse = 4, digits_rel = 4) {
  param_sets <- unique(results[, c("shape1", "scale1", "shape2", "scale2")])
  
  for (i in seq_len(nrow(param_sets))) {
    sh1 <- param_sets$shape1[i]
    sc1 <- param_sets$scale1[i]
    sh2 <- param_sets$shape2[i]
    sc2 <- param_sets$scale2[i]
    
    sub <- results[
      results$shape1 == sh1 &
        results$scale1 == sc1 &
        results$shape2 == sh2 &
        results$scale2 == sc2,
    ]
    
    cat("\n============================================================\n")
    cat(sprintf(
      "Weibull 1: shape = %.3f, scale = %.3f | Weibull 2: shape = %.3f, scale = %.3f | true D = %.8f\n",
      sh1, sc1, sh2, sc2, sub$true_D[1]
    ))
    cat("============================================================\n")
    
    cat("\nMSE table:\n")
    mse_tab <- data.frame(
      n1 = sub$n1,
      n2 = sub$n2,
      Kernel = round(sub$MSE_Kernel, digits_mse),
      Empirical = round(sub$MSE_Emp, digits_mse),
      Ustat = round(sub$MSE_Ustat, digits_mse)
    )
    print(mse_tab, row.names = FALSE)
    
    cat("\nRelative MSE table, baseline = U-statistic:\n")
    rel_tab <- data.frame(
      n1 = sub$n1,
      n2 = sub$n2,
      Kernel = round(sub$RelMSE_Kernel, digits_rel),
      Empirical = round(sub$RelMSE_Emp, digits_rel),
      Ustat = round(sub$RelMSE_Ustat, digits_rel)
    )
    print(rel_tab, row.names = FALSE)
    
    if (any(sub$Kernel_NA_Count > 0)) {
      cat("\nWarning: Some kernel estimates returned NA.\n")
      print(
        data.frame(
          n1 = sub$n1,
          n2 = sub$n2,
          Kernel_NA_Count = sub$Kernel_NA_Count
        ),
        row.names = FALSE
      )
    }
  }
}


# =============================================================================
# 10. Run simulation
# =============================================================================

results_weibull <- run_simulation(
  iterations = 2000,
  n_quad = 80,
  seed = 2024,
  verbose = TRUE
)

print_results_by_parameter(results_weibull)

write.csv(
  results_weibull,
  "MSE_and_Relative_MSE_results_Weibull.csv",
  row.names = FALSE
)

