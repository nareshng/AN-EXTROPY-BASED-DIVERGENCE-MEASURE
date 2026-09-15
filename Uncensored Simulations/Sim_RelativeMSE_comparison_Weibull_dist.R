

# =============================================================================
# Section 5.1 – comparison of Relative MSEs under Weibull Distribution
# for kernel estimators of D, D_CC, and KL
#
# Weibull simulation version
#
# Bandwidth is kept exactly as original:
# h = 0.9 * min(sd(z), IQR(z)/1.34) * n^(-0.2)
#
# Relative MSE:
# RelMSE = MSE / true_value^2
#
# =============================================================================


# =============================================================================
# 1. Silverman bandwidth
# =============================================================================

silverman_bw <- function(z) {
  0.9 * min(sd(z), IQR(z) / 1.34) * length(z)^(-0.2)
}


# =============================================================================
# 2. Gauss-Legendre quadrature nodes and weights
# =============================================================================

gauss_legendre <- function(n) {
  if (n < 2L) {
    stop("n must be at least 2.")
  }
  
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
# 3. Kernel survival and density estimators evaluated on a grid
# =============================================================================

kern_surv_grid <- function(samples, x_grid, bw) {
  colMeans(
    pnorm(
      outer(samples, x_grid, function(s, x) (s - x) / bw)
    )
  )
}


kern_dens_grid <- function(samples, x_grid, bw) {
  colMeans(
    dnorm(
      outer(samples, x_grid, function(s, x) (x - s) / bw)
    )
  ) / bw
}


# =============================================================================
# 4. All three kernel estimators
# =============================================================================

calc_all_kernel <- function(
    X,
    Y,
    gl,
    upper = NULL,
    tail_mult = 10,
    eps = 1e-12
) {
  h1 <- silverman_bw(X)
  h2 <- silverman_bw(Y)
  
 
  if (!is.finite(h1) || h1 <= 0 || !is.finite(h2) || h2 <= 0) {
    return(c(D = NA_real_, DCC = NA_real_, KL = NA_real_))
  }
  
  if (is.null(upper)) {
    upper <- max(c(X, Y)) + tail_mult * max(h1, h2)
  }
  
  if (!is.finite(upper) || upper <= 0) {
    return(c(D = NA_real_, DCC = NA_real_, KL = NA_real_))
  }
  

  x_grid <- 0.5 * upper * (gl$x + 1)
  w_grid <- 0.5 * upper * gl$w
  
  Fb <- kern_surv_grid(X, x_grid, h1)
  Gb <- kern_surv_grid(Y, x_grid, h2)
  
  f_hat <- kern_dens_grid(X, x_grid, h1)
  g_hat <- kern_dens_grid(Y, x_grid, h2)
  
  # D estimator
  D_hat <- sum(w_grid * (Fb - Gb)^2)
  
  # D_CC estimator
  DCC_hat <- sum(w_grid * abs(Fb * g_hat - Gb * f_hat))
  
  # KL estimator
  # pmax avoids log(0).
  f_safe <- pmax(f_hat, eps)
  g_safe <- pmax(g_hat, eps)
  
  KL_hat <- sum(w_grid * f_safe * log(f_safe / g_safe))
  
  c(
    D = D_hat,
    DCC = DCC_hat,
    KL = KL_hat
  )
}


# =============================================================================
# 5. True values for Weibull(shape1, scale1) vs Weibull(shape2, scale2)
# =============================================================================

weibull_surv <- function(x, shape, scale) {
  exp(- (x / scale)^shape)
}


weibull_dens <- function(x, shape, scale) {
  dweibull(x, shape = shape, scale = scale)
}


true_D_weibull <- function(shape1, scale1, shape2, scale2) {
  if (shape1 <= 0 || shape2 <= 0 || scale1 <= 0 || scale2 <= 0) {
    stop("All Weibull shape and scale parameters must be positive.")
  }
  
  term1 <- scale1 * gamma(1 + 1 / shape1) / (2^(1 / shape1))
  term2 <- scale2 * gamma(1 + 1 / shape2) / (2^(1 / shape2))
  
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


# -----------------------------------------------------------------------------
# True D_CC:
#
# D_CC(F,G) = integral_0^inf |S_X(x) g(x) - S_Y(x) f(x)| dx.
#
# For general Weibull distributions this is safest to compute numerically.
# Do not fake a closed form unless you have proved one for your exact case.
# -----------------------------------------------------------------------------

true_DCC_weibull <- function(shape1, scale1, shape2, scale2) {
  if (shape1 <= 0 || shape2 <= 0 || scale1 <= 0 || scale2 <= 0) {
    stop("All Weibull shape and scale parameters must be positive.")
  }
  
  integrand <- function(x) {
    Sx <- weibull_surv(x, shape1, scale1)
    Sy <- weibull_surv(x, shape2, scale2)
    
    fx <- weibull_dens(x, shape1, scale1)
    gy <- weibull_dens(x, shape2, scale2)
    
    abs(Sx * gy - Sy * fx)
  }
  
  integrate(
    integrand,
    lower = 0,
    upper = Inf,
    rel.tol = 1e-10,
    subdivisions = 1000
  )$value
}


# -----------------------------------------------------------------------------
# True KL:
#
# KL(f_X || f_Y), where
#
# X ~ Weibull(shape1, scale1)
# Y ~ Weibull(shape2, scale2)
#
# Closed form:
#
# KL = log(shape1 / shape2)
#      + shape2 * log(scale2 / scale1)
#      - EulerGamma * (shape1 - shape2) / shape1
#      - 1
#      + (scale1 / scale2)^shape2 * Gamma(1 + shape2 / shape1)
# -----------------------------------------------------------------------------

true_KL_weibull <- function(shape1, scale1, shape2, scale2) {
  if (shape1 <= 0 || shape2 <= 0 || scale1 <= 0 || scale2 <= 0) {
    stop("All Weibull shape and scale parameters must be positive.")
  }
  
  euler_gamma <- -digamma(1)
  
  log(shape1 / shape2) +
    shape2 * log(scale2 / scale1) -
    euler_gamma * (shape1 - shape2) / shape1 -
    1 +
    (scale1 / scale2)^shape2 * gamma(1 + shape2 / shape1)
}


relative_mse <- function(mse, true_value) {
  if (!is.finite(true_value) || abs(true_value) < .Machine$double.eps) {
    return(NA_real_)
  }
  
  mse / true_value^2
}


# =============================================================================
# 6. Simulation cell
# =============================================================================

simulate_one_cell <- function(
    shape1,
    scale1,
    shape2,
    scale2,
    n1,
    n2,
    iterations,
    gl,
    tail_mult = 10,
    eps = 1e-12
) {
  tD <- true_D_weibull(shape1, scale1, shape2, scale2)
  tDCC <- true_DCC_weibull(shape1, scale1, shape2, scale2)
  tKL <- true_KL_weibull(shape1, scale1, shape2, scale2)
  
  D_est <- numeric(iterations)
  DCC_est <- numeric(iterations)
  KL_est <- numeric(iterations)
  
  # Deterministic high quantile for better tail coverage.
  # This matters more for Weibull than for exponential, especially when shape < 1.
  q_upper <- max(
    qweibull(1 - 1e-12, shape = shape1, scale = scale1),
    qweibull(1 - 1e-12, shape = shape2, scale = scale2)
  )
  
  for (b in seq_len(iterations)) {
    X <- rweibull(n1, shape = shape1, scale = scale1)
    Y <- rweibull(n2, shape = shape2, scale = scale2)
    
    h1 <- silverman_bw(X)
    h2 <- silverman_bw(Y)
    
    upper_b <- max(
      q_upper,
      max(c(X, Y)) + tail_mult * max(h1, h2),
      na.rm = TRUE
    )
    
    est <- calc_all_kernel(
      X = X,
      Y = Y,
      gl = gl,
      upper = upper_b,
      tail_mult = tail_mult,
      eps = eps
    )
    
    D_est[b] <- est["D"]
    DCC_est[b] <- est["DCC"]
    KL_est[b] <- est["KL"]
  }
  
  D_ok <- is.finite(D_est)
  DCC_ok <- is.finite(DCC_est)
  KL_ok <- is.finite(KL_est)
  
  MSE_D <- mean((D_est[D_ok] - tD)^2)
  MSE_DCC <- mean((DCC_est[DCC_ok] - tDCC)^2)
  MSE_KL <- mean((KL_est[KL_ok] - tKL)^2)
  
  c(
    true_D = tD,
    true_DCC = tDCC,
    true_KL = tKL,
    
    MSE_D = MSE_D,
    MSE_DCC = MSE_DCC,
    MSE_KL = MSE_KL,
    
    RelMSE_D = relative_mse(MSE_D, tD),
    RelMSE_DCC = relative_mse(MSE_DCC, tDCC),
    RelMSE_KL = relative_mse(MSE_KL, tKL),
    
    NA_D = sum(!D_ok),
    NA_DCC = sum(!DCC_ok),
    NA_KL = sum(!KL_ok)
  )
}


# =============================================================================
# 7. Full simulation
# =============================================================================

run_simulation <- function(
    iterations = 2000,
    n_quad = 100,
    tail_mult = 10,
    eps = 1e-12,
    seed = 2024,
    verbose = TRUE
) {
  set.seed(seed)
  
  gl <- gauss_legendre(n_quad)
  

  param_grid <- list(
    c(0.5, 1.0, 1.0, 1.0),
    c(1.0, 1.0, 2.0, 1.0),
    c(2.0, 1.0, 1.0, 1.0),
    c(1.5, 1.0, 3.0, 1.5),
    c(0.7, 2.0, 1.5, 1.0),
    c(3.0, 1.0, 1.2, 2.0)
  )
  
  n1_vec <- c(10, 50, 80, 200)
  n2_vec <- c(10, 40, 100, 200)
  
  out <- vector("list", length(param_grid) * length(n1_vec))
  idx <- 1L
  
  for (p in seq_along(param_grid)) {
    shape1 <- param_grid[[p]][1]
    scale1 <- param_grid[[p]][2]
    shape2 <- param_grid[[p]][3]
    scale2 <- param_grid[[p]][4]
    
    if (verbose) {
      tD <- true_D_weibull(shape1, scale1, shape2, scale2)
      tDCC <- true_DCC_weibull(shape1, scale1, shape2, scale2)
      tKL <- true_KL_weibull(shape1, scale1, shape2, scale2)
      
      cat("\n============================================================\n")
      cat(sprintf(
        "Weibull 1: shape = %.3f, scale = %.3f | Weibull 2: shape = %.3f, scale = %.3f\n",
        shape1, scale1, shape2, scale2
      ))
      cat(sprintf(
        "true D = %.8f, true D_CC = %.8f, true KL = %.8f\n",
        tD, tDCC, tKL
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
        gl = gl,
        tail_mult = tail_mult,
        eps = eps
      )
      
      out[[idx]] <- data.frame(
        shape1 = shape1,
        scale1 = scale1,
        shape2 = shape2,
        scale2 = scale2,
        n1 = n1,
        n2 = n2,
        
        true_D = sim["true_D"],
        true_DCC = sim["true_DCC"],
        true_KL = sim["true_KL"],
        
        MSE_D = sim["MSE_D"],
        MSE_DCC = sim["MSE_DCC"],
        MSE_KL = sim["MSE_KL"],
        
        RelMSE_D = sim["RelMSE_D"],
        RelMSE_DCC = sim["RelMSE_DCC"],
        RelMSE_KL = sim["RelMSE_KL"],
        
        NA_D = sim["NA_D"],
        NA_DCC = sim["NA_DCC"],
        NA_KL = sim["NA_KL"],
        
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
# 8. Tables
# =============================================================================

print_results_by_parameter <- function(
    results,
    digits_mse = 4,
    digits_rel = 4
) {
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
      "Weibull 1: shape = %.3f, scale = %.3f | Weibull 2: shape = %.3f, scale = %.3f\n",
      sh1, sc1, sh2, sc2
    ))
    cat(sprintf(
      "true D = %.8f, true D_CC = %.8f, true KL = %.8f\n",
      sub$true_D[1],
      sub$true_DCC[1],
      sub$true_KL[1]
    ))
    cat("============================================================\n")
    
    mse_tab <- data.frame(
      n1 = sub$n1,
      n2 = sub$n2,
      MSE_D = formatC(sub$MSE_D, digits = digits_mse, format = "e"),
      MSE_DCC = formatC(sub$MSE_DCC, digits = digits_mse, format = "e"),
      MSE_KL = formatC(sub$MSE_KL, digits = digits_mse, format = "e")
    )
    
    cat("\nMSE table:\n")
    print(mse_tab, row.names = FALSE)
    
    rel_tab <- data.frame(
      n1 = sub$n1,
      n2 = sub$n2,
      RelMSE_D = round(sub$RelMSE_D, digits_rel),
      RelMSE_DCC = round(sub$RelMSE_DCC, digits_rel),
      RelMSE_KL = round(sub$RelMSE_KL, digits_rel)
    )
    
    cat("\nRelative MSE table: MSE / true_value^2\n")
    print(rel_tab, row.names = FALSE)
    
    na_tab <- data.frame(
      n1 = sub$n1,
      n2 = sub$n2,
      NA_D = sub$NA_D,
      NA_DCC = sub$NA_DCC,
      NA_KL = sub$NA_KL
    )
    
    cat("\nNA counts:\n")
    print(na_tab, row.names = FALSE)
  }
}


# =============================================================================
# 9. Run main Weibull simulation
# =============================================================================

results_weibull <- run_simulation(
  iterations = 2000,
  n_quad = 100,
  tail_mult = 10,
  eps = 1e-12,
  seed = 2024,
  verbose = TRUE
)

print_results_by_parameter(results_weibull)

write.csv(
  results_weibull,
  "MSE_and_Relative_MSE_Section_5_1_Weibull.csv",
  row.names = FALSE
)


