## =============================================================================
## Section 6 -- Figures 1-3 and the combined cumulative-divergence figure, drawn
## exactly as in the code that produced the submitted manuscript:
##   * Kaplan-Meier panel  : survminer::ggsurvplot(conf.int, pval, theme_bw)
##   * cumulative divergence: ggplot2::geom_step + dashed line at tau
## Exactly two PDFs per data set are written, named for the manuscript:
##   Veteran_KM_curve.pdf       Veteran_cumulative_divergence.pdf
##   Lung_cancer_KM_curve.pdf   Lung_cancer_cumulative_divergence.pdf
##   GBSG2_KM_curve.pdf         GBSG2_cumulative_divergence.pdf
## File names, page sizes, titles, subtitles and axis expressions follow the
## original script; only the numbers change, because the estimator now uses the
## corrected Greenwood scaling and the tau rule stated in Section 6.
##
## Requires ggplot2; survminer is used for the Kaplan-Meier panel when present.
## =============================================================================

## Input: a free-text label. Output: the original script's file-name slug.
fig_clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

## Input: the step table from rd_cumulative(), plus dataset and comparison names.
## Output: the step data frame in the column layout the original figures used.
fig_step_df <- function(steps, dataset, comparison) {
  data.frame(
    left_time = steps$left_time, right_time = steps$right_time,
    width = steps$width, S1_left = steps$S_group1, S2_left = steps$S_group2,
    diff_squared = steps$squared_difference, contribution = steps$contribution,
    cumulative_Dtau = steps$cumulative_D_tau,
    Dataset = dataset, Comparison = comparison, stringsAsFactors = FALSE
  )
}

## Input: the analysis row, the two samples, the step table, group labels, the
##        output directory, the file stem and the geom for the cumulative panel.
## Output: the two file paths invisibly. Exactly two PDFs are written per data
##         set: <stem>_KM_curve.pdf and <stem>_cumulative_divergence.pdf.
fig_dataset_figures <- function(row, t1, s1, t2, s2, step_df, labels, out_dir, stem,
                                cumulative_geom = c("step", "line")) {
  cumulative_geom <- match.arg(cumulative_geom)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for the manuscript figure style.", call. = FALSE)
  }
  tau <- row$tau
  subtitle <- bquote(tau == .(round(tau, 2)) ~ "," ~
                       widehat(D)[tau]^KM == .(round(row$D_tau_KM, 4)))
  km_file <- file.path(out_dir, paste0(stem, "_KM_curve.pdf"))
  cum_file <- file.path(out_dir, paste0(stem, "_cumulative_divergence.pdf"))

  ## ---- panel (a): Kaplan-Meier curves ------------------------------------------
  dat <- data.frame(time = c(t1, t2), status = c(s1, s2),
                    group = rep(labels, c(length(t1), length(t2))),
                    stringsAsFactors = FALSE)
  km_fit <- survival::survfit(survival::Surv(time, status) ~ group, data = dat)
  if (requireNamespace("survminer", quietly = TRUE)) {
    km_plot <- survminer::ggsurvplot(
      km_fit, data = dat, conf.int = TRUE, risk.table = FALSE, pval = TRUE,
      censor = TRUE, xlab = "Time", ylab = "Survival probability",
      title = paste0(row$Dataset, ": Kaplan--Meier curves"),
      legend.title = "Group", legend.labs = labels,
      ggtheme = ggplot2::theme_bw()
    )$plot
  } else {
    message("survminer is not installed: drawing the Kaplan-Meier panel with ggplot2 only.")
    km_plot <- fig_km_ggplot(km_fit, labels, row)
  }
  grDevices::pdf(km_file, width = 7.2, height = 4.8)
  print(km_plot)
  grDevices::dev.off()

  ## ---- panel (b): cumulative divergence ----------------------------------------
  cum_layer <- if (cumulative_geom == "step") {
    ggplot2::geom_step(linewidth = 0.8)
  } else {
    ggplot2::geom_line(linewidth = 0.8)
  }
  cum_plot <- ggplot2::ggplot(step_df, ggplot2::aes(x = right_time, y = cumulative_Dtau)) +
    cum_layer +
    ggplot2::geom_vline(xintercept = tau, linetype = "dashed") +
    ggplot2::labs(
      title = paste0(row$Dataset, ": cumulative divergence"),
      subtitle = subtitle, x = "Time",
      y = expression(integral((widehat(bar(F))[n[1]](u) - widehat(bar(G))[n[2]](u))^2 * du, 0, t))
    ) +
    ggplot2::theme_bw()
  ggplot2::ggsave(cum_file, cum_plot, width = 7.2, height = 4.8, device = "pdf")

  invisible(c(km = km_file, cumulative = cum_file))
}

## Input: a two-group survfit, the labels and the analysis row.
## Output: a ggplot Kaplan-Meier panel used only when survminer is unavailable.
fig_km_ggplot <- function(fit, labels, row) {
  strata <- rep(seq_along(fit$strata), fit$strata)
  curve <- do.call(rbind, lapply(seq_along(fit$strata), function(g) {
    keep <- strata == g
    data.frame(time = c(0, fit$time[keep]), surv = c(1, fit$surv[keep]),
               lower = c(1, fit$lower[keep]), upper = c(1, fit$upper[keep]),
               n.censor = c(0, fit$n.censor[keep]),
               Group = factor(labels[g], levels = labels), stringsAsFactors = FALSE)
  }))
  band <- do.call(rbind, lapply(split(curve, curve$Group), function(d) {
    d <- d[order(d$time), ]; n <- nrow(d)
    data.frame(time = c(rbind(d$time, c(d$time[-1L], d$time[n]))),
               lower = rep(d$lower, each = 2L), upper = rep(d$upper, each = 2L),
               Group = d$Group[1L], stringsAsFactors = FALSE)
  }))
  band <- band[is.finite(band$lower) & is.finite(band$upper), ]
  ggplot2::ggplot(curve, ggplot2::aes(x = time, y = surv, colour = Group)) +
    ggplot2::geom_ribbon(data = band,
                         ggplot2::aes(x = time, ymin = lower, ymax = upper,
                                      fill = Group, colour = NULL),
                         alpha = 0.2, inherit.aes = FALSE) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::geom_point(data = curve[curve$n.censor > 0, ], shape = 3, size = 1.6,
                        show.legend = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::annotate("text", x = 0, y = 0.08, hjust = 0,
                      label = sprintf("p = %.4g", row$Logrank_p)) +
    ggplot2::labs(title = paste0(row$Dataset, ": Kaplan--Meier curves"),
                  x = "Time", y = "Survival probability") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "top")
}
