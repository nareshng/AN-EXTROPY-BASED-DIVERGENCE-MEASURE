# =============================================================================
# Section 5.1 – Comparison of Relative MSEs under Exponential Distribution
# for kernel estimators of D, D_CC, and KL
#
# Bandwidth
# h = 0.9 * min(sd(z), IQR(z)/1.34) * n^(-0.2)
#
# Relative MSE:
# RelMSE = MSE / true_value^2
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
    tail_mult = 10,
    eps = 1e-12
) {
  h1 <- silverman_bw(X)
  h2 <- silverman_bw(Y)
  
  
  if (!is.finite(h1) || h1 <= 0 || !is.finite(h2) || h2 <= 0) {
    return(c(D = NA_real_, DCC = NA_real_, KL = NA_real_))
  }
  
  upper <- max(c(X, Y)) + tail_mult * max(h1, h2)
  
  if (!is.finite(upper) || upper <= 0) {
    return(c(D = NA_real_, DCC = NA_real_, KL = NA_real_))
  }
  
  # Map Gauss-Legendre rule from [-1, 1] to [0, upper]
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
  # pmax avoids log(0), but does not use fake constant tail extrapolation.
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
# 5. True values for Exp(lambda1) vs Exp(lambda2)
# =============================================================================

true_D <- function(l1, l2) {
  1 / (2 * l1) + 1 / (2 * l2) - 2 / (l1 + l2)
}


true_DCC <- function(l1, l2) {
  abs(l1 - l2) / (l1 + l2)
}


true_KL <- function(l1, l2) {
  log(l1 / l2) + l2 / l1 - 1
}


relative_mse <- function(mse, true_value) {
  if (!is.finite(true_value) || abs(true_value) < .Machine$double.eps) {
    return(NA_real_)
  }
  
  mse / true_value^2
}


# =============================================================================
# 6. simulation cell
# =============================================================================

simulate_one_cell <- function(
    lambda1,
    lambda2,
    n1,
    n2,
    iterations,
    gl,
    tail_mult = 10,
    eps = 1e-12
) {
  tD <- true_D(lambda1, lambda2)
  tDCC <- true_DCC(lambda1, lambda2)
  tKL <- true_KL(lambda1, lambda2)
  
  D_est <- numeric(iterations)
  DCC_est <- numeric(iterations)
  KL_est <- numeric(iterations)
  
  for (b in seq_len(iterations)) {
    X <- rexp(n1, rate = lambda1)
    Y <- rexp(n2, rate = lambda2)
    
    est <- calc_all_kernel(
      X = X,
      Y = Y,
      gl = gl,
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
# 7. simulation
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
    c(0.1, 0.2),
    c(0.5, 1.0),
    c(2.0, 1.0),
    c(3.0, 5.0),
    c(7.0, 5.0),
    c(10.0, 7.0)
  )
  
  n1_vec <- c(10, 50, 80, 200)
  n2_vec <- c(10, 40, 100, 200)
  
  out <- vector("list", length(param_grid) * length(n1_vec))
  idx <- 1L
  
  for (p in seq_along(param_grid)) {
    lambda1 <- param_grid[[p]][1]
    lambda2 <- param_grid[[p]][2]
    
    if (verbose) {
      cat("\n============================================================\n")
      cat(sprintf("lambda1 = %.3f, lambda2 = %.3f\n", lambda1, lambda2))
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
        gl = gl,
        tail_mult = tail_mult,
        eps = eps
      )
      
      out[[idx]] <- data.frame(
        lambda1 = lambda1,
        lambda2 = lambda2,
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
  param_pairs <- unique(results[, c("lambda1", "lambda2")])
  
  for (i in seq_len(nrow(param_pairs))) {
    l1 <- param_pairs$lambda1[i]
    l2 <- param_pairs$lambda2[i]
    
    sub <- results[results$lambda1 == l1 & results$lambda2 == l2, ]
    
    cat("\n============================================================\n")
    cat(sprintf("(lambda1, lambda2) = (%.3f, %.3f)\n", l1, l2))
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
# 9. Run
# =============================================================================

results <- run_simulation(
  iterations = 2000,
  n_quad = 100,
  tail_mult = 10,
  eps = 1e-12,
  seed = 2024,
  verbose = TRUE
)

print_results_by_parameter(results)

write.csv(
  results,
  "MSE_and_Relative_MSE_Section_5_1.csv",
  row.names = FALSE
)






# =============================================================================
# 10.quadrature sensitivity check
# =============================================================================

# If n_quad = 200 and n_quad = 300 differ materially, use the larger value.

 results_100 <- run_simulation(
   iterations = 2000,
   n_quad = 200,
   tail_mult = 10,
   seed = 2024,
   verbose = FALSE
 )

 results_200 <- run_simulation(
   iterations = 2000,
   n_quad = 300,
   tail_mult = 10,
   seed = 2026,
   verbose = FALSE
 )

check_quad <- data.frame(
   lambda1 = results_100$lambda1,
   lambda2 = results_100$lambda2,
   n1 = results_100$n1,
   n2 = results_100$n2,

   Diff_MSE_D = abs(results_100$MSE_D - results_150$MSE_D),
   Diff_MSE_DCC = abs(results_100$MSE_DCC - results_150$MSE_DCC),
   Diff_MSE_KL = abs(results_100$MSE_KL - results_150$MSE_KL)
 )

 print(check_quad)
 
 
 





