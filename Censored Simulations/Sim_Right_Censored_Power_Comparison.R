# ================================================================
# Section 5.4 : RIGHT-CENSORED POWER SIMULATION
#
# Tests:
#   D_KM        : unweighted Kaplan--Meier divergence
#   DCC_KM      : Cox--Czanner type KM statistic
#   Log-rank_tau: restricted log-rank test
#
# Scenarios:
#   1. Persistent separation
#   2. Non-proportional hazards / crossing curves
#   3. Equal median different tail
#
# Sample-size settings:
#   1. n1 = 60, n2 = 40
#   2. n1 = 80, n2 = 120
# ================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(patchwork)
})

start_time <- Sys.time()

# ================================================================
# CONFIG
# ================================================================

CFG <- list(
  scenarios = c(
    "Persistent separation",
    "Non-proportional hazards (crossing curves)",
    "Equal median different tail"
  ),
  
  n_configs = data.frame(
    n1 = c(60, 80),
    n2 = c(40, 120)
  ),
  
  censor_levels = c(0.05, 0.20, 0.40),
  tau_levels    = c(0.80, 0.90),
  
  B_mc    = 1000,
  B_perm  = 1000,
  alpha   = 0.05,
  seed    = 2026,
  plot_q  = 0.995,
  
  ref_censor = 0.20,
  ref_tau    = 0.80,
  
  out_dir  = "simulation_figures_final",
  csv_file = "simulation_figures_final/power_grid_results_final.csv"
)

stopifnot(
  "ref_tau must be in tau_levels"       = CFG$ref_tau    %in% CFG$tau_levels,
  "ref_censor must be in censor_levels" = CFG$ref_censor %in% CFG$censor_levels,
  "B_mc must be positive"               = CFG$B_mc > 0,
  "B_perm must be positive"             = CFG$B_perm > 0,
  "alpha must be in (0,1)"              = CFG$alpha > 0 && CFG$alpha < 1,
  "n_configs must contain n1"           = "n1" %in% names(CFG$n_configs),
  "n_configs must contain n2"           = "n2" %in% names(CFG$n_configs)
)

B_MC   <- CFG$B_mc
B_PERM <- CFG$B_perm
ALPHA  <- CFG$alpha

METHOD_ORDER <- c(
  "D_KM",
  "DCC_KM",
  "Log-rank_tau"
)

METHOD_LABELS <- c(
  "D_KM"         = expression(hat(D)[tau]^{KM}),
  "DCC_KM"       = expression(hat(D)[tau]^{CC}),
  "Log-rank_tau" = expression("Log-rank"["[" * tau * "]"])
)

METHOD_COLS <- c(
  "D_KM"         = "#D95F02",
  "DCC_KM"       = "#1F78B4",
  "Log-rank_tau" = "#7B2D8B"
)

SURV_COLS <- c(
  "Group 1" = "#F8766D",
  "Group 2" = "#00BFC4"
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================================================
# 1. KAPLAN--MEIER UTILITIES
# ================================================================

get_km_fit <- function(time, status) {
  survival::survfit(survival::Surv(time, status) ~ 1)
}

km_step_eval <- function(fit, x) {
  tt <- fit$time[fit$n.event > 0]
  ss <- fit$surv[fit$n.event > 0]
  
  if (length(tt) == 0) {
    return(rep(1, length(x)))
  }
  
  idx <- findInterval(x, tt)
  as.numeric(ifelse(idx == 0, 1, ss[idx]))
}

get_linear_km_points <- function(fit, tau) {
  keep <- fit$n.event > 0 & fit$time <= tau
  
  tt <- c(0, fit$time[keep])
  ss <- c(1, fit$surv[keep])
  
  if (tt[length(tt)] < tau) {
    tt <- c(tt, tau)
    ss <- c(ss, km_step_eval(fit, tau))
  }
  
  u <- !duplicated(tt)
  list(time = tt[u], surv = ss[u])
}

# ================================================================
# 2. TEST STATISTICS
# ================================================================

stat_D_KM <- function(time1, status1, time2, status2, tau) {
  fit1 <- get_km_fit(time1, status1)
  fit2 <- get_km_fit(time2, status2)
  
  knots <- sort(unique(c(
    0,
    fit1$time[fit1$n.event > 0 & fit1$time < tau],
    fit2$time[fit2$n.event > 0 & fit2$time < tau],
    tau
  )))
  
  if (length(knots) < 2) return(0)
  
  left  <- knots[-length(knots)]
  width <- diff(knots)
  
  sum(
    (km_step_eval(fit1, left) - km_step_eval(fit2, left))^2 * width
  )
}

stat_DCC_KM <- function(time1, status1, time2, status2, tau) {
  p1 <- get_linear_km_points(get_km_fit(time1, status1), tau)
  p2 <- get_linear_km_points(get_km_fit(time2, status2), tau)
  
  tt <- sort(unique(c(p1$time, p2$time)))
  tt <- tt[tt >= 0 & tt <= tau]
  
  if (length(tt) < 2) return(0)
  
  Fv <- approx(p1$time, p1$surv, xout = tt, rule = 2)$y
  Gv <- approx(p2$time, p2$surv, xout = tt, rule = 2)$y
  
  N <- length(tt)
  
  sum(abs(Fv[2:N] * Gv[1:(N - 1)] - Fv[1:(N - 1)] * Gv[2:N]))
}

stat_logrank_tau_pvalue <- function(time1, status1, time2, status2, tau) {
  t1 <- pmin(time1, tau)
  t2 <- pmin(time2, tau)
  
  s1 <- as.integer(status1 == 1 & time1 <= tau)
  s2 <- as.integer(status2 == 1 & time2 <= tau)
  
  if (sum(s1) + sum(s2) == 0) {
    return(NA_real_)
  }
  
  group <- factor(c(rep(1, length(t1)), rep(2, length(t2))))
  
  fit <- tryCatch(
    survival::survdiff(survival::Surv(c(t1, t2), c(s1, s2)) ~ group),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !is.finite(fit$chisq)) {
    return(NA_real_)
  }
  
  stats::pchisq(fit$chisq, df = 1, lower.tail = FALSE)
}

# ================================================================
# 3. PERMUTATION TEST WRAPPER
# ================================================================

perm_test_stat <- function(time1, status1, time2, status2,
                           tau, stat_fun) {
  n1 <- length(time1)
  n  <- n1 + length(time2)
  
  tall <- c(time1, time2)
  sall <- c(status1, status2)
  
  T0 <- tryCatch(
    stat_fun(time1, status1, time2, status2, tau),
    error = function(e) NA_real_
  )
  
  if (!is.finite(T0)) {
    return(NA_real_)
  }
  
  Tv <- vapply(seq_len(B_PERM), function(b) {
    i1 <- sample.int(n, n1)
    
    tryCatch(
      stat_fun(tall[i1], sall[i1], tall[-i1], sall[-i1], tau),
      error = function(e) NA_real_
    )
  }, numeric(1))
  
  Tv <- Tv[is.finite(Tv)]
  
  if (length(Tv) < max(20, floor(0.8 * B_PERM))) {
    return(NA_real_)
  }
  
  (1 + sum(Tv >= T0)) / (length(Tv) + 1)
}

# ================================================================
# 4. SCENARIO DEFINITIONS
# ================================================================

make_scenario <- function(name) {
  sc <- switch(
    name,
    
    "Persistent separation" = {
      sh  <- 1.5
      sc1 <- 1.0
      sc2 <- 1.4
      
      list(
        r1 = function(n) rweibull(n, shape = sh, scale = sc1),
        r2 = function(n) rweibull(n, shape = sh, scale = sc2),
        S1 = function(t) exp(-(t / sc1)^sh),
        S2 = function(t) exp(-(t / sc2)^sh),
        q1 = function(p) qweibull(p, shape = sh, scale = sc1),
        q2 = function(p) qweibull(p, shape = sh, scale = sc2)
      )
    },
    
    "Non-proportional hazards (crossing curves)" = {
      sh1 <- 0.6
      sc1 <- 3.0
      sh2 <- 2.5
      sc2 <- 5.2
      
      list(
        r1 = function(n) rweibull(n, shape = sh1, scale = sc1),
        r2 = function(n) rweibull(n, shape = sh2, scale = sc2),
        S1 = function(t) exp(-(t / sc1)^sh1),
        S2 = function(t) exp(-(t / sc2)^sh2),
        q1 = function(p) qweibull(p, shape = sh1, scale = sc1),
        q2 = function(p) qweibull(p, shape = sh2, scale = sc2)
      )
    },
    
    "Equal median different tail" = {
      med <- 1
      
      sh1 <- 2.5
      sh2 <- 0.7
      
      sc1 <- med / (log(2))^(1 / sh1)
      sc2 <- med / (log(2))^(1 / sh2)
      
      list(
        r1 = function(n) rweibull(n, shape = sh1, scale = sc1),
        r2 = function(n) rweibull(n, shape = sh2, scale = sc2),
        S1 = function(t) exp(-(t / sc1)^sh1),
        S2 = function(t) exp(-(t / sc2)^sh2),
        q1 = function(p) qweibull(p, shape = sh1, scale = sc1),
        q2 = function(p) qweibull(p, shape = sh2, scale = sc2)
      )
    },
    
    stop("Unknown scenario: ", name)
  )
  
  sc$tau_fun <- function(q) max(sc$q1(q), sc$q2(q))
  sc$name <- name
  
  sc
}

# ================================================================
# 5. CENSORING CALIBRATION AND DATA GENERATION
# ================================================================

calibrate_censor_rate <- function(r1, r2, target, Npilot = 1e5) {
  if (target <= 0) {
    return(0)
  }
  
  Tp <- c(r1(Npilot / 2), r2(Npilot / 2))
  
  obj <- function(mu) mean(1 - exp(-mu * Tp)) - target
  
  upper <- 1
  
  while (obj(upper) < 0) {
    upper <- upper * 2
    
    if (upper > 1e6) {
      stop("Censoring calibration failed.")
    }
  }
  
  uniroot(obj, lower = 1e-10, upper = upper)$root
}

simulate_two_sample <- function(sc, n1, n2, tau_q, censor_rate) {
  X <- sc$r1(n1)
  Y <- sc$r2(n2)
  
  if (censor_rate > 0) {
    C <- rexp(n1, rate = censor_rate)
    E <- rexp(n2, rate = censor_rate)
  } else {
    C <- rep(Inf, n1)
    E <- rep(Inf, n2)
  }
  
  list(
    time1     = pmin(X, C),
    status1   = as.integer(X <= C),
    time2     = pmin(Y, E),
    status2   = as.integer(Y <= E),
    tau       = sc$tau_fun(tau_q),
    obs_cens1 = mean(X > C),
    obs_cens2 = mean(Y > E)
  )
}

# ================================================================
# 6. POWER SIMULATION
# ================================================================

mc_ci <- function(x) {
  v <- is.finite(x)
  m <- sum(v)
  
  if (m == 0) {
    return(c(
      power = NA,
      se = NA,
      lower = NA,
      upper = NA,
      valid = 0,
      invalid = length(x)
    ))
  }
  
  k <- sum(x[v])
  p <- k / m
  z <- 1.96
  
  denom  <- 1 + z^2 / m
  center <- (p + z^2 / (2 * m)) / denom
  half   <- z * sqrt((p * (1 - p) + z^2 / (4 * m)) / m) / denom
  
  c(
    power = p,
    se = sqrt(p * (1 - p) / m),
    lower = max(0, center - half),
    upper = min(1, center + half),
    valid = m,
    invalid = length(x) - m
  )
}

run_one <- function(scenario_name, n1, n2, censor_prop, tau_q, seed) {
  set.seed(seed)
  
  sc    <- make_scenario(scenario_name)
  crate <- calibrate_censor_rate(sc$r1, sc$r2, censor_prop)
  
  rD  <- rep(NA_real_, B_MC)
  rCC <- rep(NA_real_, B_MC)
  rLR <- rep(NA_real_, B_MC)
  
  c1 <- numeric(B_MC)
  c2 <- numeric(B_MC)
  
  for (b in seq_len(B_MC)) {
    dat <- simulate_two_sample(sc, n1, n2, tau_q, crate)
    
    c1[b] <- dat$obs_cens1
    c2[b] <- dat$obs_cens2
    
    tau <- dat$tau
    
    p <- perm_test_stat(
      dat$time1, dat$status1,
      dat$time2, dat$status2,
      tau,
      stat_D_KM
    )
    rD[b] <- if (is.finite(p)) as.numeric(p < ALPHA) else NA_real_
    
    p <- perm_test_stat(
      dat$time1, dat$status1,
      dat$time2, dat$status2,
      tau,
      stat_DCC_KM
    )
    rCC[b] <- if (is.finite(p)) as.numeric(p < ALPHA) else NA_real_
    
    p <- stat_logrank_tau_pvalue(
      dat$time1, dat$status1,
      dat$time2, dat$status2,
      tau
    )
    rLR[b] <- if (is.finite(p)) as.numeric(p < ALPHA) else NA_real_
    
    if (b %% 100 == 0) {
      cat(sprintf(
        "  %s | n1=%d n2=%d | cens=%.2f tau_q=%.2f | %d/%d\n",
        scenario_name, n1, n2, censor_prop, tau_q, b, B_MC
      ))
    }
  }
  
  rows <- Map(
    function(r, m) {
      ci <- mc_ci(r)
      
      data.frame(
        Scenario     = scenario_name,
        Method       = m,
        Power        = round(100 * ci["power"], 1),
        Lower        = round(100 * ci["lower"], 1),
        Upper        = round(100 * ci["upper"], 1),
        SE           = round(100 * ci["se"], 2),
        Valid        = as.integer(ci["valid"]),
        Invalid      = as.integer(ci["invalid"]),
        n1           = n1,
        n2           = n2,
        CensorTarget = censor_prop,
        CensorObs    = round(mean(c(c1, c2)), 3),
        TauQ         = tau_q,
        B_mc         = B_MC,
        B_perm       = B_PERM,
        stringsAsFactors = FALSE
      )
    },
    list(rD, rCC, rLR),
    METHOD_ORDER
  )
  
  do.call(rbind, rows)
}

# ================================================================
# 7. GRID RUNNER
# ================================================================

run_grid <- function() {
  results <- list()
  counter <- 1
  
  for (i in seq_len(nrow(CFG$n_configs))) {
    n1 <- CFG$n_configs$n1[i]
    n2 <- CFG$n_configs$n2[i]
    
    for (sc in CFG$scenarios) {
      for (cp in CFG$censor_levels) {
        for (tq in CFG$tau_levels) {
          cat(sprintf(
            "\n[%d] %s | n1=%d n2=%d | censor=%.2f | tau_q=%.2f\n",
            counter, sc, n1, n2, cp, tq
          ))
          
          out <- tryCatch(
            run_one(sc, n1, n2, cp, tq, seed = CFG$seed + counter),
            error = function(e) {
              cat("ERROR:", e$message, "\n")
              NULL
            }
          )
          
          if (!is.null(out)) {
            out$SettingID <- counter
            results[[length(results) + 1]] <- out
          }
          
          counter <- counter + 1
        }
      }
    }
  }
  
  if (length(results) == 0) {
    stop("No simulation results produced.")
  }
  
  final <- do.call(rbind, results)
  rownames(final) <- NULL
  
  write.csv(final, CFG$csv_file, row.names = FALSE)
  cat("\nSaved:", CFG$csv_file, "\n")
  
  final
}

# ================================================================
# 8. FIGURE 1 — TRUE SURVIVAL CURVES
# ================================================================

make_survival_panel <- function(scenario_name) {
  sc <- make_scenario(scenario_name)
  
  x_max <- max(sc$q1(CFG$plot_q), sc$q2(CFG$plot_q))
  t_all <- seq(0, x_max, length.out = 1000)
  
  s1 <- sc$S1(t_all)
  s2 <- sc$S2(t_all)
  
  both_low <- which(s1 < 0.01 & s2 < 0.01)
  
  if (length(both_low)) {
    x_max <- t_all[both_low[1]] * 1.06
    t_all <- t_all[t_all <= x_max]
    s1 <- sc$S1(t_all)
    s2 <- sc$S2(t_all)
  }
  
  df_curves <- data.frame(
    t    = rep(t_all, 2),
    surv = c(s1, s2),
    grp  = rep(c("Group 1", "Group 2"), each = length(t_all))
  )
  
  ggplot(df_curves, aes(x = t, y = surv, colour = grp)) +
    geom_line(linewidth = 1.4) +
    scale_colour_manual(values = SURV_COLS) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0.01, 0.01),
      breaks = seq(0, 1, 0.25)
    ) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    labs(
      x = "Time",
      y = "Survival probability",
      colour = NULL,
      title = scenario_name
    ) +
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
  
  fig1 <- patchwork::wrap_plots(panels, ncol = 3) +
    patchwork::plot_annotation(
      title = "True survival curves",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  
  fname <- file.path(CFG$out_dir, "Figure1_survival_curves.png")
  
  ggsave(fname, fig1, width = 15, height = 5, dpi = 300)
  cat("Saved:", fname, "\n")
}

# ================================================================
# 9. FIGURE 2 — POWER GRID FOR EACH SAMPLE-SIZE SETTING
# ================================================================

make_power_panel <- function(grid, scenario_name, n1, n2, censor_prop, tau_q) {
  df <- grid[
    grid$Scenario == scenario_name &
      grid$n1 == n1 &
      grid$n2 == n2 &
      grid$CensorTarget == censor_prop &
      grid$TauQ == tau_q,
  ]
  
  if (nrow(df) == 0) {
    stop(
      "No data: scenario=", scenario_name,
      " n1=", n1,
      " n2=", n2,
      " censor=", censor_prop,
      " tau_q=", tau_q
    )
  }
  
  df$Method <- factor(df$Method, levels = rev(METHOD_ORDER))
  df <- df[order(df$Method), ]
  df$ypos <- seq_len(nrow(df))
  
  ggplot(df, aes(y = ypos)) +
    geom_segment(
      aes(x = Lower, xend = Upper, yend = ypos, colour = Method),
      linewidth = 1.2
    ) +
    geom_point(
      aes(x = Power, colour = Method),
      size = 3.2,
      shape = 21,
      fill = "white",
      stroke = 1.5
    ) +
    geom_text(
      aes(
        x = Upper + 1.5,
        label = sprintf("%.1f", Power),
        colour = Method
      ),
      hjust = 0,
      size = 3.3
    ) +
    geom_vline(
      xintercept = 5,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey60"
    ) +
    scale_colour_manual(values = METHOD_COLS, guide = "none") +
    scale_x_continuous(
      limits = c(0, 116),
      expand = c(0, 0),
      breaks = seq(0, 100, 25)
    ) +
    scale_y_continuous(
      breaks = df$ypos,
      labels = METHOD_LABELS[as.character(df$Method)]
    ) +
    labs(
      x = "Power (%)",
      y = NULL,
      title = sprintf("Censoring = %d%%", round(100 * censor_prop))
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 12),
      panel.grid.major.x = element_line(colour = "grey93"),
      axis.text.y        = element_text(size = 10.5),
      axis.line          = element_line(colour = "grey60")
    )
}

save_figure2 <- function(grid) {
  for (i in seq_len(nrow(CFG$n_configs))) {
    n1 <- CFG$n_configs$n1[i]
    n2 <- CFG$n_configs$n2[i]
    
    for (tq in CFG$tau_levels) {
      all_panels <- list()
      
      for (sc in CFG$scenarios) {
        for (cp in CFG$censor_levels) {
          all_panels[[length(all_panels) + 1]] <- tryCatch(
            make_power_panel(grid, sc, n1, n2, cp, tq),
            error = function(e) {
              cat("Panel failed:", sc, n1, n2, cp, tq, "\n  ", e$message, "\n")
              NULL
            }
          )
        }
      }
      
      all_panels <- Filter(Negate(is.null), all_panels)
      
      if (length(all_panels) == 0) {
        cat("No panels for n1=", n1, " n2=", n2, " tau_q=", tq, "\n")
        next
      }
      
      fig_title <- sprintf(
        "Empirical power | n1 = %d, n2 = %d | tau_q = %.2f | alpha = %.2f",
        n1, n2, tq, CFG$alpha
      )
      
      fig2 <- patchwork::wrap_plots(
        all_panels,
        ncol = length(CFG$censor_levels),
        nrow = length(CFG$scenarios)
      ) +
        patchwork::plot_annotation(
          title = fig_title,
          theme = theme(plot.title = element_text(face = "bold", size = 12))
        )
      
      fname <- file.path(
        CFG$out_dir,
        sprintf("Figure2_power_grid_n1_%d_n2_%d_tau%.2f.png", n1, n2, tq)
      )
      
      ggsave(
        fname,
        fig2,
        width = 5 * length(CFG$censor_levels),
        height = 4.8 * length(CFG$scenarios),
        dpi = 300
      )
      
      cat("Saved:", fname, "\n")
    }
  }
}

# ================================================================
# 10. FIGURE 3 — POWER VS CENSORING FOR EACH SAMPLE-SIZE SETTING
# ================================================================

save_figure3 <- function(grid) {
  for (i in seq_len(nrow(CFG$n_configs))) {
    n1 <- CFG$n_configs$n1[i]
    n2 <- CFG$n_configs$n2[i]
    
    df_n <- grid[grid$n1 == n1 & grid$n2 == n2, ]
    
    all_panels <- list()
    
    for (sc in CFG$scenarios) {
      for (tq in CFG$tau_levels) {
        df_sub <- df_n[df_n$Scenario == sc & df_n$TauQ == tq, ]
        
        panel <- ggplot(
          df_sub,
          aes(
            x = 100 * CensorTarget,
            y = Power,
            colour = Method,
            group = Method
          )
        ) +
          geom_line(linewidth = 1.1) +
          geom_point(size = 2.5) +
          scale_colour_manual(values = METHOD_COLS, labels = METHOD_ORDER) +
          scale_x_continuous(breaks = 100 * CFG$censor_levels) +
          coord_cartesian(ylim = c(0, 100)) +
          labs(
            title = sc,
            subtitle = sprintf("tau_q = %.2f", tq),
            x = "Target censoring (%)",
            y = "Power (%)",
            colour = "Method"
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
    }
    
    fig3 <- patchwork::wrap_plots(
      all_panels,
      ncol = length(CFG$tau_levels),
      nrow = length(CFG$scenarios),
      guides = "collect"
    ) +
      patchwork::plot_annotation(
        title = sprintf(
          "Power vs. censoring | n1 = %d, n2 = %d | alpha = %.2f",
          n1, n2, CFG$alpha
        ),
        theme = theme(
          plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom"
        )
      ) &
      theme(legend.position = "bottom")
    
    fname <- file.path(
      CFG$out_dir,
      sprintf("Figure3_power_by_censoring_n1_%d_n2_%d.png", n1, n2)
    )
    
    ggsave(
      fname,
      fig3,
      width = 5 * length(CFG$tau_levels),
      height = 4.5 * length(CFG$scenarios),
      dpi = 300
    )
    
    cat("Saved:", fname, "\n")
  }
}

# ================================================================
# 11. RUN EVERYTHING
# ================================================================

cat("\n=== CONFIG ===\n")
cat(sprintf(
  "  Scenarios     : %s\n",
  paste(CFG$scenarios, collapse = "\n                  ")
))
cat("  Sample sizes  :\n")
print(CFG$n_configs)
cat(sprintf("  tau_levels    : %s\n", paste(CFG$tau_levels, collapse = ", ")))
cat(sprintf("  censor_levels : %s\n", paste(CFG$censor_levels, collapse = ", ")))
cat(sprintf("  B_mc / B_perm : %d / %d\n", CFG$B_mc, CFG$B_perm))
cat(sprintf("  alpha         : %.2f\n", CFG$alpha))
cat(sprintf("  ref_tau       : %.2f\n", CFG$ref_tau))
cat(sprintf("  ref_censor    : %.2f\n", CFG$ref_censor))
cat(sprintf("  Output dir    : %s\n", CFG$out_dir))
cat("==============\n\n")

grid <- run_grid()

cat("\n--- Results at reference censoring and tau ---\n")
print(
  grid[
    grid$CensorTarget == CFG$ref_censor &
      grid$TauQ == CFG$ref_tau,
    c(
      "Scenario", "Method", "Power", "Lower", "Upper",
      "Valid", "Invalid", "n1", "n2"
    )
  ]
)

save_figure1()
save_figure2(grid)
save_figure3(grid)

cat("\n=== DONE ===\n")
cat("All output saved to:", CFG$out_dir, "\n")
cat(" ", basename(CFG$csv_file), "\n")
cat("  Figure1_survival_curves.png\n")
cat("  Figure2_power_grid_n1_60_n2_40_tau0.80.png\n")
cat("  Figure2_power_grid_n1_60_n2_40_tau0.90.png\n")
cat("  Figure2_power_grid_n1_80_n2_120_tau0.80.png\n")
cat("  Figure2_power_grid_n1_80_n2_120_tau0.90.png\n")
cat("  Figure3_power_by_censoring_n1_60_n2_40.png\n")
cat("  Figure3_power_by_censoring_n1_80_n2_120.png\n")

end_time <- Sys.time()
time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("Time taken: %.2f hours\n", time_taken))
