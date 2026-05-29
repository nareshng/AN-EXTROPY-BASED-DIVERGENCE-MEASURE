# =============================================================================


# Divergence for Weibull distributions
#
# MSE and Relative MSE comparison of Kernel, Empirical, and U-statistic estimators
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
# 1. Code
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
# 2. bandwidth rule
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
  
  # Average ranks are safer for real data.
  # For continuous exponential simulation, ties occur with probability zero.
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
# 5. Kernel estimator 
# =============================================================================

calc_Kernel_GL <- function(X, Y, gl, upper = NULL) {
  n1 <- length(X)
  n2 <- length(Y)
  
  h1 <- bw_original(X)
  h2 <- bw_original(Y)
  
  # This keeps the original bandwidth formula.
  # The check below only prevents numerical crash in degenerate cases.
  # For exponential simulations, h1 and h2 should be positive almost surely.
  if (!is.finite(h1) || h1 <= 0 || !is.finite(h2) || h2 <= 0) {
    return(NA_real_)
  }
  
  # Finite upper bound approximation for integral_0^inf.
  # This is much faster than integrate(..., 0, Inf).
  if (is.null(upper)) {
    upper <- max(c(X, Y)) + 8 * max(h1, h2)
  }
  
  if (!is.finite(upper) || upper <= 0) {
    return(NA_real_)
  }
  
  # Transform Gauss-Legendre nodes from [-1,1] to [0, upper]
  xg <- 0.5 * upper * (gl$x + 1)
  wg <- 0.5 * upper * gl$w
  
  Fbar <- colMeans(pnorm(outer(X, xg, function(xi, x) (xi - x) / h1)))
  Gbar <- colMeans(pnorm(outer(Y, xg, function(yi, x) (yi - x) / h2)))
  
  sum(wg * (Fbar - Gbar)^2)
}


# =============================================================================
# 6. True divergence for exponential distributions
# =============================================================================

true_D_exp <- function(lambda1, lambda2) {
  1 / (2 * lambda1) + 1 / (2 * lambda2) - 2 / (lambda1 + lambda2)
}


# =============================================================================
# 7. simulation cell
# =============================================================================

simulate_one_cell <- function(lambda1, lambda2, n1, n2, iterations, gl) {
  true_D <- true_D_exp(lambda1, lambda2)
  
  ker_est <- numeric(iterations)
  emp_est <- numeric(iterations)
  ust_est <- numeric(iterations)
  
  for (b in seq_len(iterations)) {
    X <- rexp(n1, rate = lambda1)
    Y <- rexp(n2, rate = lambda2)
    
    ker_est[b] <- calc_Kernel_GL(X, Y, gl)
    emp_est[b] <- calc_Emp(X, Y)
    ust_est[b] <- calc_Ustat(X, Y)
  }
  
  # Remove possible NA kernel values.
  # For the exponential simulation this should almost never happen.
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
# 8. simulation
# =============================================================================

run_simulation <- function(
    iterations = 2000,
    n_quad = 80,
    seed = 2024,
    verbose = TRUE
) {
  set.seed(seed)
  
  gl <- gauss_legendre(n_quad)
  
  param_grid <- list(
    c(0.5, 0.1),
    c(0.1, 0.5),
    c(0.1, 1.0),
    c(1.0, 0.5),
    c(2.0, 5.0),
    c(10.0, 5.0)
  )
  
  n1_vec <- c(5,  10, 20, 20, 40, 30, 50)
  n2_vec <- c(10, 10, 10, 20, 30, 40, 50)
  
  out <- vector("list", length(param_grid) * length(n1_vec))
  idx <- 1L
  
  for (p in seq_along(param_grid)) {
    lambda1 <- param_grid[[p]][1]
    lambda2 <- param_grid[[p]][2]
    true_D <- true_D_exp(lambda1, lambda2)
    
    if (verbose) {
      cat("\n============================================================\n")
      cat(sprintf(
        "lambda1 = %.3f, lambda2 = %.3f, true D = %.8f\n",
        lambda1, lambda2, true_D
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
        lambda1 = lambda1,
        lambda2 = lambda2,
        n1 = n1,
        n2 = n2,
        iterations = iterations,
        gl = gl
      )
      
      out[[idx]] <- data.frame(
        lambda1 = lambda1,
        lambda2 = lambda2,
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
# 9.table output
# =============================================================================

print_results_by_parameter <- function(results, digits_mse = 4, digits_rel = 4) {
  param_pairs <- unique(results[, c("lambda1", "lambda2")])
  
  for (i in seq_len(nrow(param_pairs))) {
    l1 <- param_pairs$lambda1[i]
    l2 <- param_pairs$lambda2[i]
    
    sub <- results[results$lambda1 == l1 & results$lambda2 == l2, ]
    
    cat("\n============================================================\n")
    cat(sprintf(
      "(lambda1, lambda2) = (%.3f, %.3f), true D = %.8f\n",
      l1, l2, sub$true_D[1]
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
  }
}




# =============================================================================
# 10. Run  simulation
# =============================================================================

results <- run_simulation(
  iterations = 2000,
  n_quad = 80,
  seed = 2024,
  verbose = TRUE
)

print_results_by_parameter(results)

write.csv(results, "MSE_and_Relative_MSE_results.csv", row.names = FALSE)


