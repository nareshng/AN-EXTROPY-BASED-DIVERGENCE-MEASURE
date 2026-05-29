# =============================================================================
# Censored Real Data Analysis
# =============================================================================
#
# Cases:
#   1. Veteran lung cancer: Standard treatment vs Test treatment
#   2. Lung cancer: Male vs Female
#   3. GBSG2 breast cancer: No hormonal therapy vs Hormonal therapy
#
#
# Internal event coding:
#   status = 1 means event/death/failure
#   status = 0 means censored
#
# =============================================================================


# =============================================================================
# 0. Packages
# =============================================================================

required_packages <- c(
  "survival", "survminer", "dplyr", "ggplot2",
  "tibble", "readr", "knitr", "TH.data"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tibble)
library(readr)
library(knitr)
library(TH.data)


# =============================================================================
# 1. Output folders
# =============================================================================

out_dir <- "censored_real_data_three_cases_outputs"
fig_dir <- file.path(out_dir, "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nOutputs will be saved in:\n")
cat(normalizePath(out_dir, mustWork = FALSE), "\n\n")


# =============================================================================
# 2. Utility functions
# =============================================================================

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}


# =============================================================================
# 3. Kaplan-Meier step-function utilities
# =============================================================================

fit_km_one_group <- function(time, status) {
  survival::survfit(survival::Surv(time, status) ~ 1)
}


eval_km_step <- function(km_fit, t) {
  if (length(km_fit$time) == 0) {
    return(rep(1, length(t)))
  }
  
  stats::approx(
    x = c(0, km_fit$time),
    y = c(1, km_fit$surv),
    xout = t,
    method = "constant",
    f = 0,
    rule = 2
  )$y
}


choose_tau_from_event_quantile <- function(time1, status1,
                                           time2, status2,
                                           tau_q = 0.80) {
  event_times1 <- time1[status1 == 1]
  event_times2 <- time2[status2 == 1]
  
  if (length(event_times1) < 2 || length(event_times2) < 2) {
    stop("Too few events in one of the groups to choose tau.")
  }
  
  tau <- max(
    as.numeric(stats::quantile(event_times1, probs = tau_q, na.rm = TRUE)),
    as.numeric(stats::quantile(event_times2, probs = tau_q, na.rm = TRUE))
  )
  
  if (!is.finite(tau) || tau <= 0) {
    stop("Invalid tau selected. Check event coding and event times.")
  }
  
  tau
}


make_step_grid <- function(km1, km2, tau) {
  grid <- sort(unique(c(0, km1$time, km2$time, tau)))
  grid <- grid[grid >= 0 & grid <= tau]
  grid <- sort(unique(grid))
  
  if (length(grid) < 2) {
    stop("Integration grid has fewer than two points.")
  }
  
  if (tail(grid, 1) < tau) {
    grid <- sort(unique(c(grid, tau)))
  }
  
  grid
}


# =============================================================================
# 4. Exact computation of D_tau^KM
# =============================================================================

compute_Dtau_KM_core <- function(time1, status1,
                                 time2, status2,
                                 tau) {
  km1 <- fit_km_one_group(time1, status1)
  km2 <- fit_km_one_group(time2, status2)
  
  grid <- make_step_grid(km1, km2, tau)
  
  left_times <- grid[-length(grid)]
  right_times <- grid[-1]
  widths <- right_times - left_times
  
  S1_left <- eval_km_step(km1, left_times)
  S2_left <- eval_km_step(km2, left_times)
  
  diff_left <- S1_left - S2_left
  contribution <- diff_left^2 * widths
  
  Dtau <- sum(contribution)
  
  step_df <- tibble::tibble(
    left_time = left_times,
    right_time = right_times,
    width = widths,
    S_group1 = S1_left,
    S_group2 = S2_left,
    diff = diff_left,
    diff_squared = diff_left^2,
    contribution = contribution,
    cumulative_Dtau = cumsum(contribution)
  )
  
  list(
    Dtau = Dtau,
    tau = tau,
    km1 = km1,
    km2 = km2,
    grid = grid,
    step_df = step_df
  )
}


# =============================================================================
# 5. Greenwood plug-in variance estimator
# =============================================================================

greenwood_component <- function(km_event_fit,
                                step_df,
                                survival_col,
                                tau) {
  km_tab <- tibble::tibble(
    time = km_event_fit$time,
    n_risk = km_event_fit$n.risk,
    n_event = km_event_fit$n.event
  ) %>%
    dplyr::filter(n_event > 0, time <= tau) %>%
    dplyr::mutate(
      denom = n_risk * (n_risk - n_event),
      greenwood_increment = dplyr::if_else(
        denom > 0,
        n_event / denom,
        NA_real_
      )
    ) %>%
    dplyr::filter(is.finite(greenwood_increment))
  
  if (nrow(km_tab) == 0) {
    return(0)
  }
  
  vals <- numeric(nrow(km_tab))
  
  for (i in seq_len(nrow(km_tab))) {
    t0 <- km_tab$time[i]
    
    temp <- step_df %>%
      dplyr::filter(left_time >= t0, left_time < tau)
    
    if (nrow(temp) == 0) {
      vals[i] <- 0
    } else {
      S_col <- temp[[survival_col]]
      vals[i] <- sum(temp$diff * S_col * temp$width)
    }
  }
  
  sum(km_tab$greenwood_increment * vals^2)
}


compute_greenwood_CI <- function(time1, status1,
                                 time2, status2,
                                 tau,
                                 alpha = 0.05) {
  n1 <- length(time1)
  n2 <- length(time2)
  n <- n1 + n2
  
  core <- compute_Dtau_KM_core(time1, status1, time2, status2, tau)
  
  sigma_F_sq <- greenwood_component(
    km_event_fit = core$km1,
    step_df = core$step_df,
    survival_col = "S_group1",
    tau = tau
  )
  
  sigma_G_sq <- greenwood_component(
    km_event_fit = core$km2,
    step_df = core$step_df,
    survival_col = "S_group2",
    tau = tau
  )
  
  sigma_tau_sq <- (4 * n / n1) * sigma_F_sq +
    (4 * n / n2) * sigma_G_sq
  
  sigma_tau <- sqrt(max(sigma_tau_sq, 0))
  z <- stats::qnorm(1 - alpha / 2)
  
  lower <- core$Dtau - z * sigma_tau / sqrt(n)
  upper <- core$Dtau + z * sigma_tau / sqrt(n)
  
  list(
    Dtau = core$Dtau,
    tau = tau,
    sigma_tau_sq = sigma_tau_sq,
    sigma_tau = sigma_tau,
    lower = max(lower, 0),
    upper = upper,
    step_df = core$step_df
  )
}


# =============================================================================
# 6. Analysis function for one dataset
# =============================================================================

analyze_censored_dataset <- function(dat,
                                     dataset_name,
                                     time_col,
                                     status_col,
                                     group_col,
                                     group1,
                                     group2,
                                     group1_label = as.character(group1),
                                     group2_label = as.character(group2),
                                     tau = NULL,
                                     tau_q = 0.80,
                                     B = 1000,
                                     alpha = 0.05,
                                     seed = 123,
                                     save_plots = TRUE) {
  
  set.seed(seed)
  
  dat_clean <- dat %>%
    dplyr::select(
      time = dplyr::all_of(time_col),
      status = dplyr::all_of(status_col),
      group = dplyr::all_of(group_col)
    ) %>%
    dplyr::filter(!is.na(time), !is.na(status), !is.na(group)) %>%
    dplyr::filter(group %in% c(group1, group2)) %>%
    dplyr::mutate(
      group = factor(
        group,
        levels = c(group1, group2),
        labels = c(group1_label, group2_label)
      ),
      status = as.integer(status),
      time = as.numeric(time)
    ) %>%
    dplyr::filter(time >= 0)
  
  d1 <- dat_clean %>% dplyr::filter(group == group1_label)
  d2 <- dat_clean %>% dplyr::filter(group == group2_label)
  
  if (nrow(d1) < 5 || nrow(d2) < 5) {
    stop(paste(dataset_name, ": one group has too few observations."))
  }
  
  if (sum(d1$status == 1) < 2 || sum(d2$status == 1) < 2) {
    print(table(dat_clean$group, dat_clean$status))
    stop(paste(dataset_name, ": one group has too few events. Check event coding."))
  }
  
  if (is.null(tau)) {
    tau <- choose_tau_from_event_quantile(
      d1$time, d1$status,
      d2$time, d2$status,
      tau_q = tau_q
    )
  }
  
  core <- compute_Dtau_KM_core(
    d1$time, d1$status,
    d2$time, d2$status,
    tau = tau
  )
  
  greenwood <- compute_greenwood_CI(
    d1$time, d1$status,
    d2$time, d2$status,
    tau = tau,
    alpha = alpha
  )
  
  logrank_fit <- survival::survdiff(
    survival::Surv(time, status) ~ group,
    data = dat_clean
  )
  
  logrank_p <- 1 - stats::pchisq(
    logrank_fit$chisq,
    df = length(logrank_fit$n) - 1
  )
  
  # ---------------------------------------------------------------------------
  # Bootstrap CI
  # ---------------------------------------------------------------------------
  
  boot_vals <- numeric(B)
  
  for (b in seq_len(B)) {
    b1 <- d1[sample(seq_len(nrow(d1)), size = nrow(d1), replace = TRUE), ]
    b2 <- d2[sample(seq_len(nrow(d2)), size = nrow(d2), replace = TRUE), ]
    
    bcore <- compute_Dtau_KM_core(
      b1$time, b1$status,
      b2$time, b2$status,
      tau = tau
    )
    
    boot_vals[b] <- bcore$Dtau
  }
  
  boot_ci <- stats::quantile(
    boot_vals,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  )
  
  boot_se <- stats::sd(boot_vals, na.rm = TRUE)
  
  result_row <- tibble::tibble(
    Dataset = dataset_name,
    Comparison = paste0(group1_label, " vs ", group2_label),
    Group_1 = group1_label,
    Group_2 = group2_label,
    n1 = nrow(d1),
    n2 = nrow(d2),
    Events_1 = sum(d1$status == 1),
    Events_2 = sum(d2$status == 1),
    Censoring_1_percent = 100 * mean(d1$status == 0),
    Censoring_2_percent = 100 * mean(d2$status == 0),
    tau = tau,
    Dtau_KM = core$Dtau,
    Bootstrap_SE = boot_se,
    Bootstrap_CI_lower = boot_ci[1],
    Bootstrap_CI_upper = boot_ci[2],
    Greenwood_CI_lower = greenwood$lower,
    Greenwood_CI_upper = greenwood$upper,
    Logrank_p_value = logrank_p
  )
  
  step_df <- core$step_df %>%
    dplyr::mutate(
      Dataset = dataset_name,
      Comparison = paste0(group1_label, " vs ", group2_label)
    )
  
  plot_base_name <- clean_filename(
    paste(dataset_name, group1_label, "vs", group2_label)
  )
  
  km_plot_file <- file.path(fig_dir, paste0(plot_base_name, "_KM_curve.pdf"))
  diff_plot_file <- file.path(fig_dir, paste0(plot_base_name, "_squared_difference.pdf"))
  cum_plot_file <- file.path(fig_dir, paste0(plot_base_name, "_cumulative_divergence.pdf"))
  
  if (save_plots) {
    
    # -------------------------------------------------------------------------
    # Kaplan-Meier plot WITHOUT number-at-risk table
    # -------------------------------------------------------------------------
    
    km_fit_grouped <- survival::survfit(
      survival::Surv(time, status) ~ group,
      data = dat_clean
    )
    
    km_plot <- survminer::ggsurvplot(
      km_fit_grouped,
      data = dat_clean,
      conf.int = TRUE,
      risk.table = FALSE,      # removes bottom number-at-risk part
      pval = TRUE,
      censor = TRUE,
      xlab = "Time",
      ylab = "Survival probability",
      title = paste0(dataset_name, ": Kaplan-Meier curves"),
      legend.title = "Group",
      legend.labs = c(group1_label, group2_label),
      ggtheme = ggplot2::theme_bw()
    )
    
    grDevices::pdf(km_plot_file, width = 7.2, height = 4.8)
    print(km_plot$plot)
    grDevices::dev.off()
    
    # -------------------------------------------------------------------------
    # Pointwise squared survival-difference plot
    # -------------------------------------------------------------------------
    
    diff_plot <- ggplot2::ggplot(
      step_df,
      ggplot2::aes(x = left_time, y = diff_squared)
    ) +
      ggplot2::geom_step(linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = tau, linetype = "dashed") +
      ggplot2::labs(
        title = paste0(dataset_name, ": squared survival difference"),
        subtitle = bquote(
          tau == .(round(tau, 2)) ~ "," ~
            widehat(D)[tau]^KM == .(round(core$Dtau, 4))
        ),
        x = "Time",
        y = expression((widehat(bar(F))[n[1]](t) - widehat(bar(G))[n[2]](t))^2)
      ) +
      ggplot2::theme_bw()
    
    ggplot2::ggsave(
      filename = diff_plot_file,
      plot = diff_plot,
      width = 7.2,
      height = 4.8,
      device = "pdf"
    )
    
    # -------------------------------------------------------------------------
    # Cumulative divergence plot with proper tau notation
    # -----------------------------------------------------------ƒ--------------
    
    cum_plot <- ggplot2::ggplot(
      step_df,
      ggplot2::aes(x = right_time, y = cumulative_Dtau)
    ) +
      ggplot2::geom_step(linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = tau, linetype = "dashed") +
      ggplot2::labs(
        title = paste0(dataset_name, ": cumulative divergence"),
        subtitle = bquote(
          tau == .(round(tau, 2)) ~ "," ~
            widehat(D)[tau]^KM == .(round(core$Dtau, 4))
        ),
        x = "Time",
        y = expression(
          integral((widehat(bar(F))[n[1]](u) - widehat(bar(G))[n[2]](u))^2 * du, 0, t)
        )
      ) +
      ggplot2::theme_bw()
    
    ggplot2::ggsave(
      filename = cum_plot_file,
      plot = cum_plot,
      width = 7.2,
      height = 4.8,
      device = "pdf"
    )
  }
  
  list(
    result = result_row,
    step_df = step_df,
    boot_vals = boot_vals,
    data = dat_clean,
    tau = tau,
    km_plot_file = km_plot_file,
    diff_plot_file = diff_plot_file,
    cum_plot_file = cum_plot_file
  )
}


# =============================================================================
# 7. Prepare datasets with correct event coding
# =============================================================================

# -----------------------------------------------------------------------------
# 7.1 Veteran lung cancer dataset
# survival::veteran:
#   status: 0 = censored, 1 = event/death
#   trt: 1 = standard treatment, 2 = test treatment
# -----------------------------------------------------------------------------

data(veteran, package = "survival")

veteran_dat <- veteran %>%
  dplyr::transmute(
    time = time,
    status = status,
    trt = trt
  )


# -----------------------------------------------------------------------------
# 7.2 Lung cancer dataset
# survival::lung:
#   status: 1 = censored, 2 = event/death
#   sex: 1 = male, 2 = female
# -----------------------------------------------------------------------------

data(lung, package = "survival")

lung_dat <- lung %>%
  dplyr::transmute(
    time = time,
    status = ifelse(status == 2, 1, 0),
    sex = sex
  )


# -----------------------------------------------------------------------------
# 7.3 GBSG2 breast cancer dataset
# TH.data::GBSG2:
#   time = recurrence-free survival time
#   cens = event indicator, 1 = event, 0 = censored
#   horTh = hormonal therapy yes/no
# -----------------------------------------------------------------------------

data("GBSG2", package = "TH.data")

gbsg2_dat <- GBSG2 %>%
  dplyr::transmute(
    time = time,
    status = cens,
    horTh = horTh
  )


# =============================================================================
# 8. Check event coding
# =============================================================================

cat("\nEvent coding checks:\n")

cat("\nVeteran: treatment by status\n")
print(table(veteran_dat$trt, veteran_dat$status))

cat("\nLung: sex by status\n")
print(table(lung_dat$sex, lung_dat$status))

cat("\nGBSG2: hormonal therapy by status\n")
print(table(gbsg2_dat$horTh, gbsg2_dat$status))


# =============================================================================
# 9. Run selected analyses
# =============================================================================

# Use 1000 while testing. Use 2000 or 5000 for final paper results.
B_boot <- 1000

res_veteran <- analyze_censored_dataset(
  dat = veteran_dat,
  dataset_name = "Veteran lung cancer",
  time_col = "time",
  status_col = "status",
  group_col = "trt",
  group1 = 1,
  group2 = 2,
  group1_label = "Standard treatment",
  group2_label = "Test treatment",
  tau_q = 0.80,
  B = B_boot,
  seed = 1001
)

res_lung <- analyze_censored_dataset(
  dat = lung_dat,
  dataset_name = "Lung cancer",
  time_col = "time",
  status_col = "status",
  group_col = "sex",
  group1 = 1,
  group2 = 2,
  group1_label = "Male",
  group2_label = "Female",
  tau_q = 0.80,
  B = B_boot,
  seed = 1002
)

res_gbsg2 <- analyze_censored_dataset(
  dat = gbsg2_dat,
  dataset_name = "GBSG2 breast cancer",
  time_col = "time",
  status_col = "status",
  group_col = "horTh",
  group1 = "no",
  group2 = "yes",
  group1_label = "No hormonal therapy",
  group2_label = "Hormonal therapy",
  tau_q = 0.80,
  B = B_boot,
  seed = 1003
)


# =============================================================================
# 10. Combine results
# =============================================================================

all_results_raw <- dplyr::bind_rows(
  res_veteran$result,
  res_lung$result,
  res_gbsg2$result
)

all_step_dfs <- dplyr::bind_rows(
  res_veteran$step_df,
  res_lung$step_df,
  res_gbsg2$step_df
)


# =============================================================================
# 11. Paper-ready result table
# =============================================================================

paper_results <- all_results_raw %>%
  dplyr::mutate(
    Censoring_1_percent = round(Censoring_1_percent, 1),
    Censoring_2_percent = round(Censoring_2_percent, 1),
    tau = round(tau, 2),
    Dtau_KM = round(Dtau_KM, 4),
    Bootstrap_SE = round(Bootstrap_SE, 4),
    Bootstrap_CI = paste0(
      "(",
      round(Bootstrap_CI_lower, 4),
      ", ",
      round(Bootstrap_CI_upper, 4),
      ")"
    ),
    Greenwood_CI = paste0(
      "(",
      round(Greenwood_CI_lower, 4),
      ", ",
      round(Greenwood_CI_upper, 4),
      ")"
    ),
    Logrank_p_value = ifelse(
      Logrank_p_value < 0.001,
      "<0.001",
      as.character(signif(Logrank_p_value, 3))
    )
  ) %>%
  dplyr::select(
    Dataset,
    Comparison,
    n1,
    n2,
    Events_1,
    Events_2,
    Censoring_1_percent,
    Censoring_2_percent,
    tau,
    Dtau_KM,
    Bootstrap_CI,
    Greenwood_CI,
    Logrank_p_value
  )

cat("\nPaper-ready table:\n")
print(paper_results)


# =============================================================================
# 12. Save tables
# =============================================================================

readr::write_csv(
  all_results_raw,
  file.path(out_dir, "three_cases_censored_real_data_results_raw.csv")
)

readr::write_csv(
  paper_results,
  file.path(out_dir, "three_cases_censored_real_data_paper_table.csv")
)

readr::write_csv(
  all_step_dfs,
  file.path(out_dir, "three_cases_stepwise_divergence_values.csv")
)


# =============================================================================
# 13. Combined cumulative-divergence figure only
# =============================================================================

combined_cum_plot <- ggplot2::ggplot(
  all_step_dfs,
  ggplot2::aes(x = right_time, y = cumulative_Dtau)
) +
  ggplot2::geom_step(linewidth = 0.8) +
  ggplot2::facet_wrap(~ Dataset, scales = "free") +
  ggplot2::labs(
    title = expression("Cumulative Kaplan--Meier divergence"),
    x = "Time",
    y = expression(
      integral((widehat(bar(F))[n[1]](u) - widehat(bar(G))[n[2]](u))^2 * du, 0, t)
    )
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = file.path(fig_dir, "three_cases_combined_cumulative_divergence.pdf"),
  plot = combined_cum_plot,
  width = 10,
  height = 5.8,
  device = "pdf"
)


# =============================================================================
# 14. Print saved files
# =============================================================================

cat("\n============================================================\n")
cat("Analysis complete.\n")
cat("Tables saved in:\n")
cat(normalizePath(out_dir, mustWork = FALSE), "\n\n")

cat("Figures saved in:\n")
cat(normalizePath(fig_dir, mustWork = FALSE), "\n")
cat("============================================================\n\n")

cat("Saved files:\n")
print(list.files(out_dir, recursive = TRUE, full.names = TRUE))


