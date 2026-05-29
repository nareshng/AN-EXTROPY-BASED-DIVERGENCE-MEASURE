# ================================================================
# RIGHT-CENSORED POWER SIMULATION
#
# Tests     : D_KM, DCC_KM, KLstar_IPCW, restricted Log-rank
# Scenarios : Persistent separation, Strong tail difference,
#             Non-proportional hazards (crossing curves)
# ================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(patchwork)
})

start_time <- Sys.time()

# ================================================================
# CONFIG — edit only here
# ================================================================

CFG <- list(
  scenarios     = c("Persistent separation",
                    "Strong tail difference",
                    "Non-proportional hazards (crossing curves)"),
  n             = 100,
  censor_levels = c(0.05, 0.20, 0.40),
  tau_levels    = c(0.80, 0.90),
  B_mc          = 1000,
  B_perm        = 1000,
  alpha         = 0.05,
  seed          = 2026,
  plot_q        = 0.995,
  ref_censor    = 0.20,
  ref_tau       = 0.80,          # used by Figure 2 (main power grid)
  out_dir       = "simulation_figures",
  csv_file      = "simulation_figures/power_grid_results_v6.csv"
)

stopifnot(
  "ref_tau must be in tau_levels"       = CFG$ref_tau    %in% CFG$tau_levels,
  "ref_censor must be in censor_levels" = CFG$ref_censor %in% CFG$censor_levels,
  "B_mc must be positive"               = CFG$B_mc > 0,
  "B_perm must be positive"             = CFG$B_perm > 0,
  "alpha must be in (0,1)"              = CFG$alpha > 0 && CFG$alpha < 1,
  "n must be positive"                  = CFG$n > 0
)

N_PER_GROUP <- CFG$n
B_MC        <- CFG$B_mc
B_PERM      <- CFG$B_perm
ALPHA       <- CFG$alpha

# Fixed method order throughout
METHOD_ORDER <- c("D_KM", "DCC_KM", "KLstar_IPCW", "Log-rank_tau")

# Paper-style axis labels
METHOD_LABELS <- c(
  "D_KM"         = expression(hat(D)[tau]^{KM}),
  "DCC_KM"       = expression(hat(D)[tau]^{CC}),
  "KLstar_IPCW"  = expression(hat(KL)[tau * "," ~ IPCW]^"*"),
  "Log-rank_tau" = expression("Log-rank"["[" * tau * "]"])
)

# Power plot colours
METHOD_COLS <- c(
  "D_KM"         = "#D95F02",
  "DCC_KM"       = "#1F78B4",
  "KLstar_IPCW"  = "#33A02C",
  "Log-rank_tau" = "#7B2D8B"
)

# Survival curve colours — survminer / ggsurvplot defaults
SURV_COLS <- c(
  "Group 1" = "#F8766D",   # salmon
  "Group 2" = "#00BFC4"    # teal
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================================================
# 1. GAUSS-LEGENDRE QUADRATURE
# ================================================================

gauss_legendre <- function(n = 120) {
  i    <- 1:(n - 1)
  beta <- i / sqrt(4 * i^2 - 1)
  J    <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- beta
  J[cbind(i + 1, i)] <- beta
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  list(nodes   = eig$values[ord],
       weights = 2 * eig$vectors[1, ord]^2)
}

GL <- gauss_legendre(120)

# ================================================================
# 2. KAPLAN-MEIER UTILITIES
# ================================================================

get_km_fit <- function(time, status) survfit(Surv(time, status) ~ 1)

km_step_eval <- function(fit, x) {
  tt <- fit$time[fit$n.event > 0]
  ss <- fit$surv[fit$n.event > 0]
  if (length(tt) == 0) return(rep(1, length(x)))
  as.numeric(ifelse(findInterval(x, tt) == 0, 1, ss[findInterval(x, tt)]))
}

get_linear_km_points <- function(fit, tau) {
  keep <- fit$n.event > 0 & fit$time <= tau
  tt   <- c(0, fit$time[keep])
  ss   <- c(1, fit$surv[keep])
  if (tt[length(tt)] < tau) {
    tt <- c(tt, tau)
    ss <- c(ss, km_step_eval(fit, tau))
  }
  u <- !duplicated(tt)
  list(time = tt[u], surv = ss[u])
}

# ================================================================
# 3. IPCW KERNEL DENSITY UTILITIES
# ================================================================

silverman_bw <- function(x) {
  x  <- x[is.finite(x)]; n <- length(x)
  if (n <= 1) return(1)
  sc <- min(sd(x), IQR(x) / 1.34)
  if (!is.finite(sc) || sc <= 0) sc <- sd(x)
  if (!is.finite(sc) || sc <= 0) sc <- max(abs(x), na.rm = TRUE) / 5
  if (!is.finite(sc) || sc <= 0) sc <- 1
  0.9 * sc * n^(-0.2)
}

km_censor_surv_left <- function(time, status, eval_time) {
  fit <- survfit(Surv(time, 1L - as.integer(status)) ~ 1)
  tt  <- fit$time[fit$n.event > 0]
  ss  <- fit$surv[fit$n.event > 0]
  if (length(tt) == 0) return(rep(1, length(eval_time)))
  idx <- findInterval(eval_time, tt, left.open = TRUE)
  pmax(ifelse(idx == 0, 1, ss[idx]), 1e-6)
}

ipcw_kernel_density <- function(xgrid, time, status, h) {
  if (!is.finite(h) || h <= 0) return(rep(NA_real_, length(xgrid)))
  wts  <- as.numeric(status) / km_censor_surv_left(time, status, time)
  z    <- outer(xgrid, time, "-") / h
  dens <- rowSums(matrix(wts, nrow = length(xgrid),
                         ncol = length(time), byrow = TRUE) *
                    dnorm(z)) / (length(time) * h)
  pmax(dens, 1e-12)
}

# ================================================================
# 4. TEST STATISTICS
# ================================================================

stat_D_KM <- function(time1, status1, time2, status2, tau) {
  fit1  <- get_km_fit(time1, status1)
  fit2  <- get_km_fit(time2, status2)
  knots <- sort(unique(c(0,
                         fit1$time[fit1$n.event > 0 & fit1$time < tau],
                         fit2$time[fit2$n.event > 0 & fit2$time < tau],
                         tau)))
  if (length(knots) < 2) return(0)
  lp <- knots[-length(knots)]
  sum((km_step_eval(fit1, lp) - km_step_eval(fit2, lp))^2 * diff(knots))
}

stat_DCC_KM <- function(time1, status1, time2, status2, tau) {
  p1 <- get_linear_km_points(get_km_fit(time1, status1), tau)
  p2 <- get_linear_km_points(get_km_fit(time2, status2), tau)
  tt <- sort(unique(c(p1$time, p2$time)))
  tt <- tt[tt >= 0 & tt <= tau]
  if (length(tt) < 2) return(0)
  Fv <- approx(p1$time, p1$surv, xout = tt, rule = 2)$y
  Gv <- approx(p2$time, p2$surv, xout = tt, rule = 2)$y
  N  <- length(tt)
  sum(abs(Fv[2:N] * Gv[1:(N-1)] - Fv[1:(N-1)] * Gv[2:N]))
}

stat_KLstar_IPCW <- function(time1, status1, time2, status2,
                             tau, eps = 1e-8) {
  if (sum(status1) < 3 || sum(status2) < 3) return(NA_real_)
  fit1 <- get_km_fit(time1, status1)
  fit2 <- get_km_fit(time2, status2)
  Ft   <- 1 - km_step_eval(fit1, tau)
  Gt   <- 1 - km_step_eval(fit2, tau)
  if (!is.finite(Ft) || !is.finite(Gt) ||
      Ft <= 1e-8 || Gt <= 1e-8) return(NA_real_)
  h1 <- silverman_bw(time1[status1 == 1])
  h2 <- silverman_bw(time2[status2 == 1])
  xg <- tau / 2 * (GL$nodes + 1)
  fc <- ipcw_kernel_density(xg, time1, status1, h1) / Ft
  gc <- ipcw_kernel_density(xg, time2, status2, h2) / Gt
  if (any(!is.finite(fc)) || any(!is.finite(gc))) return(NA_real_)
  tau / 2 * sum(GL$weights * fc * log((fc + eps) / (gc + eps)))
}

stat_KLstar_sym_IPCW <- function(time1, status1, time2, status2, tau) {
  a <- stat_KLstar_IPCW(time1, status1, time2, status2, tau)
  b <- stat_KLstar_IPCW(time2, status2, time1, status1, tau)
  if (!is.finite(a) || !is.finite(b)) return(NA_real_)
  a + b
}

stat_logrank_tau_pvalue <- function(time1, status1, time2, status2, tau) {
  t1 <- pmin(time1, tau); t2 <- pmin(time2, tau)
  s1 <- as.integer(status1 == 1 & time1 <= tau)
  s2 <- as.integer(status2 == 1 & time2 <= tau)
  if (sum(s1) + sum(s2) == 0) return(NA_real_)
  group <- factor(c(rep(1, length(t1)), rep(2, length(t2))))
  fit   <- tryCatch(
    survdiff(Surv(c(t1, t2), c(s1, s2)) ~ group),
    error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$chisq)) return(NA_real_)
  pchisq(fit$chisq, df = 1, lower.tail = FALSE)
}

# ================================================================
# 5. PERMUTATION TEST WRAPPER
# ================================================================

perm_test_stat <- function(time1, status1, time2, status2,
                           tau, stat_fun) {
  n1   <- length(time1); n <- n1 + length(time2)
  tall <- c(time1, time2); sall <- c(status1, status2)
  T0   <- tryCatch(stat_fun(time1, status1, time2, status2, tau),
                   error = function(e) NA_real_)
  if (!is.finite(T0)) return(NA_real_)
  Tv <- vapply(seq_len(B_PERM), function(b) {
    i1 <- sample.int(n, n1)
    tryCatch(stat_fun(tall[i1], sall[i1], tall[-i1], sall[-i1], tau),
             error = function(e) NA_real_)
  }, numeric(1))
  Tv <- Tv[is.finite(Tv)]
  if (length(Tv) < max(20, floor(0.8 * B_PERM))) return(NA_real_)
  (1 + sum(Tv >= T0)) / (length(Tv) + 1)
}

# ================================================================
# 6. PIECEWISE EXPONENTIAL HELPERS
# ================================================================

r_pw <- function(n, le = 1, ll = 0.20, m = log(2)) {
  U  <- runif(n); Sm <- exp(-le * m); T <- numeric(n)
  T[U >= Sm] <- -log(U[U >= Sm]) / le
  T[U <  Sm] <- m + (-log(U[U < Sm]) - le * m) / ll
  T
}
S_pw <- function(t, le = 1, ll = 0.20, m = log(2))
  ifelse(t <= m, exp(-le * t), exp(-le * m - ll * (t - m)))
q_pw <- function(p, le = 1, ll = 0.20, m = log(2)) {
  ts <- 1 - p; Sm <- exp(-le * m)
  ifelse(ts >= Sm, -log(ts) / le, m + (-log(ts) - le * m) / ll)
}

# ================================================================
# 7. SCENARIO DEFINITIONS
#    tau_fun = max of the two group quantiles
# ================================================================

make_scenario <- function(name) {
  sc <- switch(name,
               
               "Persistent separation" = {
                 sh <- 1.5; sc1 <- 1.0; sc2 <- 1.4
                 list(r1 = function(n) rweibull(n, sh, sc1),
                      r2 = function(n) rweibull(n, sh, sc2),
                      S1 = function(t) exp(-(t / sc1)^sh),
                      S2 = function(t) exp(-(t / sc2)^sh),
                      q1 = function(p) qweibull(p, sh, sc1),
                      q2 = function(p) qweibull(p, sh, sc2))
               },
               
               "Strong tail difference" = {
                 m <- log(2)
                 list(r1 = function(n) rexp(n, 1),
                      r2 = function(n) r_pw(n, m = m),
                      S1 = function(t) exp(-t),
                      S2 = function(t) S_pw(t, m = m),
                      q1 = function(p) qexp(p, 1),
                      q2 = function(p) q_pw(p, m = m))
               },
               
               "Non-proportional hazards (crossing curves)" = {
                 sh1 <- 0.6; sc1 <- 3.0; sh2 <- 2.5; sc2 <- 5.2
                 list(r1 = function(n) rweibull(n, sh1, sc1),
                      r2 = function(n) rweibull(n, sh2, sc2),
                      S1 = function(t) exp(-(t / sc1)^sh1),
                      S2 = function(t) exp(-(t / sc2)^sh2),
                      q1 = function(p) qweibull(p, sh1, sc1),
                      q2 = function(p) qweibull(p, sh2, sc2))
               },
               
               stop("Unknown scenario: ", name)
  )
  sc$tau_fun <- function(q) max(sc$q1(q), sc$q2(q))
  sc$name    <- name
  sc
}

# ================================================================
# 8. CENSORING CALIBRATION & DATA GENERATION
# ================================================================

calibrate_censor_rate <- function(r1, r2, target, Npilot = 1e5) {
  if (target <= 0) return(0)
  Tp  <- c(r1(Npilot / 2), r2(Npilot / 2))
  obj <- function(mu) mean(1 - exp(-mu * Tp)) - target
  upper <- 1
  while (obj(upper) < 0) {
    upper <- upper * 2
    if (upper > 1e6) stop("Censoring calibration failed.")
  }
  uniroot(obj, lower = 1e-10, upper = upper)$root
}

simulate_two_sample <- function(sc, tau_q, censor_rate) {
  X <- sc$r1(N_PER_GROUP); Y <- sc$r2(N_PER_GROUP)
  if (censor_rate > 0) {
    C <- rexp(N_PER_GROUP, censor_rate)
    E <- rexp(N_PER_GROUP, censor_rate)
  } else {
    C <- rep(Inf, N_PER_GROUP)
    E <- rep(Inf, N_PER_GROUP)
  }
  list(time1     = pmin(X, C), status1   = as.integer(X <= C),
       time2     = pmin(Y, E), status2   = as.integer(Y <= E),
       tau       = sc$tau_fun(tau_q),
       obs_cens1 = mean(X > C), obs_cens2 = mean(Y > E))
}

# ================================================================
# 9. POWER SIMULATION
# ================================================================

mc_ci <- function(x) {
  v <- is.finite(x); m <- sum(v)
  if (m == 0) return(c(power=NA, se=NA, lower=NA, upper=NA,
                       valid=0, invalid=length(x)))
  k <- sum(x[v]); p <- k / m; z <- 1.96
  denom  <- 1 + z^2 / m
  center <- (p + z^2 / (2 * m)) / denom
  half   <- z * sqrt((p * (1 - p) + z^2 / (4 * m)) / m) / denom
  c(power   = p,
    se      = sqrt(p * (1 - p) / m),
    lower   = max(0, center - half),
    upper   = min(1, center + half),
    valid   = m,
    invalid = length(x) - m)
}

run_one <- function(scenario_name, censor_prop, tau_q, seed) {
  set.seed(seed)
  sc    <- make_scenario(scenario_name)
  crate <- calibrate_censor_rate(sc$r1, sc$r2, censor_prop)
  
  rD <- rDCC <- rKL <- rLR <- rep(NA_real_, B_MC)
  c1 <- c2 <- numeric(B_MC)
  
  for (b in seq_len(B_MC)) {
    dat   <- simulate_two_sample(sc, tau_q, crate)
    c1[b] <- dat$obs_cens1; c2[b] <- dat$obs_cens2
    tau   <- dat$tau
    
    p <- perm_test_stat(dat$time1, dat$status1,
                        dat$time2, dat$status2, tau, stat_D_KM)
    rD[b] <- if (is.finite(p)) as.numeric(p < ALPHA) else NA_real_
    
    p <- perm_test_stat(dat$time1, dat$status1,
                        dat$time2, dat$status2, tau, stat_DCC_KM)
    rDCC[b] <- if (is.finite(p)) as.numeric(p < ALPHA) else NA_real_
    
    p <- perm_test_stat(dat$time1, dat$status1,
                        dat$time2, dat$status2, tau, stat_KLstar_sym_IPCW)
    rKL[b] <- if (is.finite(p)) as.numeric(p < ALPHA) else NA_real_
    
    p <- stat_logrank_tau_pvalue(dat$time1, dat$status1,
                                 dat$time2, dat$status2, tau)
    rLR[b] <- if (is.finite(p)) as.numeric(p < ALPHA) else NA_real_
    
    if (b %% 100 == 0)
      cat(sprintf("  %s | cens=%.2f tau_q=%.2f | %d/%d\n",
                  scenario_name, censor_prop, tau_q, b, B_MC))
  }
  
  rows <- Map(function(r, m) {
    ci <- mc_ci(r)
    data.frame(Scenario     = scenario_name,
               Method       = m,
               Power        = round(100 * ci["power"], 1),
               Lower        = round(100 * ci["lower"], 1),
               Upper        = round(100 * ci["upper"], 1),
               SE           = round(100 * ci["se"],    2),
               Valid        = as.integer(ci["valid"]),
               Invalid      = as.integer(ci["invalid"]),
               n            = N_PER_GROUP,
               CensorTarget = censor_prop,
               CensorObs    = round(mean(c(c1, c2)), 3),
               TauQ         = tau_q,
               B_mc         = B_MC,
               B_perm       = B_PERM,
               stringsAsFactors = FALSE)
  }, list(rD, rDCC, rKL, rLR),
  list("D_KM", "DCC_KM", "KLstar_IPCW", "Log-rank_tau"))
  
  do.call(rbind, rows)
}

# ================================================================
# 10. GRID RUNNER
# ================================================================

run_grid <- function() {
  results <- list(); counter <- 1
  for (sc in CFG$scenarios)
    for (cp in CFG$censor_levels)
      for (tq in CFG$tau_levels) {
        cat(sprintf("\n[%d] %s | censor=%.2f | tau_q=%.2f\n",
                    counter, sc, cp, tq))
        out <- tryCatch(
          run_one(sc, cp, tq, seed = CFG$seed + counter),
          error = function(e) { cat("ERROR:", e$message, "\n"); NULL })
        if (!is.null(out)) {
          out$SettingID     <- counter
          results[[length(results) + 1]] <- out
        }
        counter <- counter + 1
      }
  if (length(results) == 0) stop("No simulation results produced.")
  final <- do.call(rbind, results)
  rownames(final) <- NULL
  write.csv(final, CFG$csv_file, row.names = FALSE)
  cat("\nSaved:", CFG$csv_file, "\n")
  final
}

# ================================================================
# 11. FIGURE 1 — TRUE SURVIVAL CURVES
#     Clean: no tau line, no shading, no subtitle
#     Colours: survminer defaults (#F8766D salmon, #00BFC4 teal)
#     Output: PNG only
# ================================================================

make_survival_panel <- function(scenario_name) {
  sc    <- make_scenario(scenario_name)
  x_max <- max(sc$q1(CFG$plot_q), sc$q2(CFG$plot_q))
  t_all <- seq(0, x_max, length.out = 1000)
  s1    <- sc$S1(t_all); s2 <- sc$S2(t_all)
  
  both_low <- which(s1 < 0.01 & s2 < 0.01)
  if (length(both_low)) {
    x_max <- t_all[both_low[1]] * 1.06
    t_all  <- t_all[t_all <= x_max]
    s1     <- sc$S1(t_all); s2 <- sc$S2(t_all)
  }
  
  df_curves <- data.frame(
    t    = rep(t_all, 2),
    surv = c(s1, s2),
    grp  = rep(c("Group 1", "Group 2"), each = length(t_all))
  )
  
  ggplot(df_curves, aes(x = t, y = surv, colour = grp)) +
    geom_line(linewidth = 1.4) +
    scale_colour_manual(values = SURV_COLS) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01),
                       breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    labs(x      = "Time",
         y      = "Survival probability",
         colour = NULL,
         title  = scenario_name) +
    theme_classic(base_size = 12) +
    theme(
      plot.title        = element_text(face = "bold", size = 12),
      legend.position   = c(0.78, 0.88),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key.width  = unit(1.2, "cm"),
      legend.text       = element_text(size = 9),
      panel.grid.major  = element_line(colour = "grey93"),
      axis.line         = element_line(colour = "grey60")
    )
}

save_figure1 <- function() {
  panels <- lapply(CFG$scenarios, make_survival_panel)
  
  fig1 <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title = "True survival curves",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  
  # PNG only
  fname <- file.path(CFG$out_dir, "Figure1_survival_curves.png")
  ggsave(fname, fig1, width = 15, height = 5, dpi = 300)
  cat("Saved:", fname, "\n")
}

# ================================================================
# 12. FIGURE 2 — POWER GRID
#     Rows    = scenarios (3)
#     Columns = censoring levels (3)
#     One figure per tau_q level in CFG$tau_levels
#     Output: PNG only
# ================================================================

make_power_panel <- function(grid, scenario_name, censor_prop, tau_q) {
  df <- grid[grid$Scenario     == scenario_name &
               grid$CensorTarget == censor_prop   &
               grid$TauQ         == tau_q,        ]
  
  if (nrow(df) == 0)
    stop("No data: scenario=", scenario_name,
         " censor=", censor_prop, " tau_q=", tau_q)
  
  df$Method <- factor(df$Method, levels = rev(METHOD_ORDER))
  df        <- df[order(df$Method), ]
  df$ypos   <- seq_len(nrow(df))
  is_hero   <- df$Method == "D_KM"
  
  ggplot(df, aes(y = ypos)) +
    geom_segment(
      aes(x = Lower, xend = Upper, yend = ypos, colour = Method),
      linewidth = ifelse(is_hero, 2.4, 1.0)
    ) +
    geom_point(
      aes(x = Power, colour = Method),
      size  = ifelse(is_hero, 4.5, 2.8),
      shape = 21, fill = "white", stroke = 1.8
    ) +
    geom_text(
      aes(x        = Upper + 1.5,
          label    = sprintf("%.1f", Power),
          colour   = Method,
          fontface = ifelse(Method == "D_KM", "bold", "plain")),
      hjust = 0, size = 3.5
    ) +
    geom_vline(xintercept = 5, linetype = "dashed",
               linewidth = 0.5, colour = "grey60") +
    scale_colour_manual(values = METHOD_COLS, guide = "none") +
    scale_x_continuous(limits = c(0, 116), expand = c(0, 0),
                       breaks = seq(0, 100, 25)) +
    scale_y_continuous(
      breaks = df$ypos,
      labels = METHOD_LABELS[as.character(df$Method)]
    ) +
    labs(x     = "Power (%)",
         y     = NULL,
         title = sprintf("Censoring = %d%%", round(100 * censor_prop))) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 12),
      panel.grid.major.x = element_line(colour = "grey93"),
      axis.text.y        = element_text(size = 10.5),
      axis.line          = element_line(colour = "grey60")
    )
}

save_figure2 <- function(grid) {
  # One power grid figure per tau_q level
  for (tq in CFG$tau_levels) {
    all_panels <- list()
    for (sc in CFG$scenarios)
      for (cp in CFG$censor_levels)
        all_panels[[length(all_panels) + 1]] <- tryCatch(
          make_power_panel(grid, sc, cp, tq),
          error = function(e) {
            cat("Panel failed:", sc, cp, tq, "\n  ", e$message, "\n")
            NULL })
    
    all_panels <- Filter(Negate(is.null), all_panels)
    if (length(all_panels) == 0) {
      cat("No panels for tau_q =", tq, "— skipping.\n"); next
    }
    
    # Figure title uses Unicode escapes so math renders correctly
    fig_title <- sprintf(
      "Empirical power |  \u03C4_q = %.2f  |  \u03B1 = %.2f",
       tq, CFG$alpha)
    
    fig2 <- wrap_plots(all_panels,
                       ncol = length(CFG$censor_levels),
                       nrow = length(CFG$scenarios)) +
      plot_annotation(
        title = fig_title,
        theme = theme(plot.title = element_text(face = "bold", size = 12))
      )
    
    # PNG only; filename encodes tau level
    fname <- file.path(CFG$out_dir,
                       sprintf("Figure2_power_grid_tau%.2f.png", tq))
    ggsave(fname, fig2,
           width  = 5 * length(CFG$censor_levels),
           height = 4.8 * length(CFG$scenarios),
           dpi    = 300)
    cat("Saved:", fname, "\n")
  }
}

# ================================================================
# 13. FIGURE 3 — POWER VS CENSORING  (3 scenarios × 2 tau levels)
#     Layout : rows = scenarios, cols = tau_q levels
#     No CI bands.  Output: PNG only.
# ================================================================

save_figure3 <- function(grid) {
  
  # build one panel per (scenario, tau_q) combination
  # rows = scenarios, cols = tau_q levels
  all_panels <- list()
  for (sc in CFG$scenarios)
    for (tq in CFG$tau_levels) {
      df_sub <- grid[grid$Scenario == sc & grid$TauQ == tq, ]
      
      panel <- ggplot(df_sub,
                      aes(x = 100 * CensorTarget, y = Power,
                          colour = Method, group = Method)) +
        geom_line(linewidth = 1.1) +
        geom_point(size = 2.5) +
        scale_colour_manual(values = METHOD_COLS,
                            labels = METHOD_ORDER) +
        scale_x_continuous(breaks = 100 * CFG$censor_levels) +
        coord_cartesian(ylim = c(0, 100)) +
        labs(
          title    = sc,
          subtitle = sprintf("\u03C4_q = %.2f", tq),
          x        = "Target censoring (%)",
          y        = "Power (%)",
          colour   = "Method"
        ) +
        theme_bw(base_size = 11) +
        theme(
          plot.title       = element_text(face = "bold", size = 10),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      all_panels[[length(all_panels) + 1]] <- panel
    }
  
  # 3 rows (scenarios) × 2 cols (tau levels)
  # patchwork collects legends automatically with guides = "collect"
  fig3 <- wrap_plots(all_panels,
                     ncol = length(CFG$tau_levels),
                     nrow = length(CFG$scenarios),
                     guides = "collect") +
    plot_annotation(
      title = sprintf(
        "Power vs. censoring  |  \u03B1 = %.2f",
         CFG$alpha),
      theme = theme(
        plot.title      = element_text(face = "bold", size = 12),
        legend.position = "bottom")
    ) &
    theme(legend.position = "bottom")
  
  fname <- file.path(CFG$out_dir, "Figure3_power_by_censoring.png")
  ggsave(fname, fig3,
         width  = 5 * length(CFG$tau_levels),
         height = 4.5 * length(CFG$scenarios),
         dpi    = 300)
  cat("Saved:", fname, "\n")
}

# ================================================================
# 14. RUN EVERYTHING
# ================================================================

cat("\n=== CONFIG ===\n")
cat(sprintf("  Scenarios     : %s\n",
            paste(CFG$scenarios, collapse = "\n                  ")))
cat(sprintf("  n per group   : %d\n",   CFG$n))
cat(sprintf("  tau_levels    : %s\n",   paste(CFG$tau_levels,    collapse = ", ")))
cat(sprintf("  censor_levels : %s\n",   paste(CFG$censor_levels, collapse = ", ")))
cat(sprintf("  B_mc / B_perm : %d / %d\n", CFG$B_mc, CFG$B_perm))
cat(sprintf("  alpha         : %.2f\n", CFG$alpha))
cat(sprintf("  ref_tau       : %.2f\n", CFG$ref_tau))
cat(sprintf("  ref_censor    : %.2f\n", CFG$ref_censor))
cat(sprintf("  Output dir    : %s\n",   CFG$out_dir))
cat("==============\n\n")

grid <- run_grid()

cat("\n--- Results (censor = 0.20, tau_q =", CFG$ref_tau, ") ---\n")
print(grid[grid$CensorTarget == CFG$ref_censor & grid$TauQ == CFG$ref_tau,
           c("Scenario", "Method", "Power", "Lower", "Upper",
             "Valid", "Invalid")])

save_figure1()
save_figure2(grid)
save_figure3(grid)

cat("\n=== DONE ===\n")
cat("All output saved to:", CFG$out_dir, "\n")
cat(" ", basename(CFG$csv_file), "\n")
cat("  Figure1_survival_curves.png\n")
cat("  Figure2_power_grid_tau0.80.png\n")
cat("  Figure2_power_grid_tau0.90.png\n")
cat("  Figure3_power_by_censoring.png\n")

end_time   <- Sys.time()
time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat(sprintf("Time taken: %.1f seconds\n", time_taken))