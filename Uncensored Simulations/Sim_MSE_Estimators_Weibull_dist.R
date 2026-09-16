

# =============================================================================
# Divergence for Weibull distributions
#
# For X ~ Weibull(shape1, scale1), Y ~ Weibull(shape2, scale2),
#
# # MSE comparison of Kernel, Empirical, and U-statistic estimators
#
#
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
                              n1, n2, iterations, gl,
                              gl_check = NULL,
                              quadrature_tolerance = 0.01) {
  
  true_D <- true_D_weibull(shape1, scale1, shape2, scale2)
  
  ker_est <- numeric(iterations)
  emp_est <- numeric(iterations)
  ust_est <- numeric(iterations)
  
  for (b in seq_len(iterations)) {
    X <- rweibull(n1, shape = shape1, scale = scale1)
    Y <- rweibull(n2, shape = shape2, scale = scale2)
    
    h1 <- bw_original(X)
    h2 <- bw_original(Y)
    
    # The estimated survival curves are Gaussian-kernel mixtures, so their
    # relevant numerical support is determined by the observed samples and
    # bandwidths, not by an extreme quantile of the generating distribution.
    # In particular, qweibull(1 - 1e-12, shape = 0.5) is about 763; mapping
    # only 80 Gauss-Legendre nodes over that range gives poor resolution near
    # zero. Using eight bandwidths beyond the largest observation makes the
    # omitted Gaussian-kernel tail negligible while keeping the quadrature
    # stable.
    upper_b <- max(c(X, Y)) + 8 * max(h1, h2)
    
    ker_est[b] <- calc_Kernel_GL(X, Y, gl, upper = upper_b)

    # Same-data sensitivity check: changing the quadrature rule must not be
    # confounded with a new Monte Carlo sample or a different random seed.
    if (b == 1L && !is.null(gl_check)) {
      ker_check <- calc_Kernel_GL(X, Y, gl_check, upper = upper_b)
      check_scale <- max(abs(ker_check), abs(true_D), sqrt(.Machine$double.eps))
      check_difference <- abs(ker_est[b] - ker_check) / check_scale

      if (!is.finite(check_difference) ||
          check_difference > quadrature_tolerance) {
        stop(sprintf(
          paste0(
            "Gauss-Legendre sensitivity check failed (relative difference %.4g) ",
            "for shape1=%g, scale1=%g, shape2=%g, scale2=%g, n1=%d, n2=%d. ",
            "Increase EXTROPY_N_QUAD."
          ),
          check_difference,
          shape1, scale1, shape2, scale2, n1, n2
        ))
      }
    }

    emp_est[b] <- calc_Emp(X, Y)
    ust_est[b] <- calc_Ustat(X, Y)
  }
  
  valid_ker <- is.finite(ker_est)

  # All B planned Monte Carlo replications must contribute to each MSE.
  # Do not silently discard a failed estimate, since doing so changes the
  # Monte Carlo estimand and can make a table irreproducible.
  if (!all(valid_ker)) {
    stop(sprintf(
      paste0(
        "Kernel estimation failed in %d of %d replications for ",
        "shape1=%g, scale1=%g, shape2=%g, scale2=%g, n1=%d, n2=%d."
      ),
      sum(!valid_ker), iterations,
      shape1, scale1, shape2, scale2, n1, n2
    ))
  }

  if (!all(is.finite(emp_est)) || !all(is.finite(ust_est))) {
    stop(sprintf(
      paste0(
        "A non-finite empirical or U-statistic estimate was produced for ",
        "shape1=%g, scale1=%g, shape2=%g, scale2=%g, n1=%d, n2=%d."
      ),
      shape1, scale1, shape2, scale2, n1, n2
    ))
  }
  
  MSE_Kernel <- mean((ker_est - true_D)^2)
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
  if (length(iterations) != 1L || !is.numeric(iterations) ||
      !is.finite(iterations) || iterations < 1 ||
      iterations != floor(iterations) || iterations > .Machine$integer.max) {
    stop("iterations must be a positive integer.")
  }
  if (length(n_quad) != 1L || !is.numeric(n_quad) ||
      !is.finite(n_quad) || n_quad < 2 ||
      n_quad != floor(n_quad) ||
      n_quad > (.Machine$integer.max %/% 2L)) {
    stop("n_quad must be an integer of at least 2 that can safely be doubled.")
  }
  if (length(seed) != 1L || !is.numeric(seed) ||
      !is.finite(seed) || seed < 0 ||
      seed != floor(seed) || seed > .Machine$integer.max) {
    stop("seed must be a nonnegative integer.")
  }
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  iterations <- as.integer(iterations)
  n_quad <- as.integer(n_quad)
  seed <- as.integer(seed)

  RNGkind(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  set.seed(seed)
  
  gl <- gauss_legendre(n_quad)
  gl_check <- gauss_legendre(2L * n_quad)
  
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
        gl = gl,
        gl_check = gl_check
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
# 10. Reproducible standalone interface
# =============================================================================

read_env_integer <- function(name, default, minimum) {
  raw <- trimws(Sys.getenv(name, unset = ""))

  if (!nzchar(raw)) {
    return(as.integer(default))
  }

  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || !is.finite(value) ||
      value < minimum || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(sprintf(
      "%s must be an integer greater than or equal to %d; received '%s'.",
      name, minimum, raw
    ))
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

  stop(sprintf(
    "%s must be one of true/false, 1/0, yes/no; received '%s'.",
    name, raw
  ))
}


prepare_output_dir <- function(path) {
  if (!nzchar(path)) {
    stop("EXTROPY_OUTPUT_DIR must not be empty.")
  }
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) {
      stop(sprintf("Could not create output directory: %s", path))
    }
  }
  if (file.access(path, mode = 2L) != 0L) {
    stop(sprintf("Output directory is not writable: %s", path))
  }

  normalizePath(path, winslash = "/", mustWork = TRUE)
}


# Defaults reproduce the paper simulation. Environment variables permit a
# deterministic quick check without editing the statistical code. A rerun with
# EXTROPY_N_QUAD=160 provides a direct quadrature-sensitivity check.
mc_reps <- read_env_integer("EXTROPY_MC_REPS", 2000L, 1L)
quadrature_nodes <- read_env_integer("EXTROPY_N_QUAD", 80L, 2L)
simulation_seed <- read_env_integer("EXTROPY_SEED", 2024L, 0L)
show_progress <- read_env_flag("EXTROPY_VERBOSE", TRUE)
output_dir <- prepare_output_dir(
  Sys.getenv("EXTROPY_OUTPUT_DIR", unset = ".")
)


# =============================================================================
# 11. Run simulation
# =============================================================================

results_weibull <- run_simulation(
  iterations = mc_reps,
  n_quad = quadrature_nodes,
  seed = simulation_seed,
  verbose = show_progress
)

if (show_progress) {
  print_results_by_parameter(results_weibull)
}

write.csv(
  results_weibull,
  file.path(output_dir, "MSE_and_Relative_MSE_results_Weibull.csv"),
  row.names = FALSE
)

output_file <- file.path(output_dir, "MSE_and_Relative_MSE_results_Weibull.csv")
if (!file.exists(output_file)) {
  stop(sprintf("Expected output file was not created: %s", output_file))
}

message(sprintf("Wrote Table 4 simulation results to %s", output_file))
