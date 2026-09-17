#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
run_root <- NULL
for (arg in args) {
  if (startsWith(arg, "--run-root=")) {
    run_root <- sub("^--run-root=", "", arg)
  } else if (arg %in% c("--help", "-h")) {
    cat("Usage: Rscript make_manuscript_tables.R --run-root=PATH\n")
    quit(save = "no", status = 0L)
  } else {
    stop("Unknown option: ", arg, call. = FALSE)
  }
}
if (is.null(run_root) || !dir.exists(run_root)) {
  stop("--run-root must identify a completed master-run directory.",
       call. = FALSE)
}
run_root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)
paper_dir <- file.path(run_root, "paper_outputs")
dir.create(paper_dir, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path) || dir.exists(path) || file.info(path)$size <= 0) {
    stop("Required output is missing or empty: ", path, call. = FALSE)
  }
  path
}

available_modes <- c("full", "quick")[dir.exists(file.path(
  run_root, "uncensored", c("full", "quick")
))]
if (length(available_modes) != 1L) {
  stop("Could not identify exactly one uncensored run mode under ", run_root,
       call. = FALSE)
}
mode <- available_modes[[1L]]

## Tables 1-6 already have validated manuscript-formatted CSVs.
for (table in 1:6) {
  source_file <- need(file.path(
    run_root, "uncensored", mode, paste0("table_", table),
    paste0("Table_", table, "_manuscript.csv")
  ))
  file.copy(source_file,
            file.path(paper_dir, paste0("Table_", table, "_manuscript.csv")),
            overwrite = TRUE)
}

read_table <- function(name) {
  utils::read.csv(
    need(file.path(run_root, "censored", "paper_tables", name)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

## Tables 7-8: retain only columns printed in the manuscript.
point_sizes <- c("(20,20)", "(50,50)", "(100,100)", "(200,200)")
for (item in list(
  list(number = 7L, file = "Table7_RelMSE_Exponential.csv"),
  list(number = 8L, file = "Table8_RelMSE_Weibull.csv")
)) {
  x <- read_table(item$file)
  keep <- c("Setting", "Censoring", paste("D_KM", point_sizes))
  if (!all(keep %in% names(x))) {
    stop("Unexpected columns in ", item$file, call. = FALSE)
  }
  utils::write.csv(x[keep], file.path(
    paper_dir, paste0("Table_", item$number, "_manuscript.csv")
  ), row.names = FALSE)
}

## Tables 9-10: retain CP and AL and remove targets, MCSE and diagnostics.
ci_sizes <- c("(30,40)", "(70,50)", "(100,100)")
ci_measures <- c("Normal_CP", "Normal_AL", "Bootstrap_CP", "Bootstrap_AL")
for (item in list(
  list(number = 9L, file = "Table9_CI_Exponential.csv"),
  list(number = 10L, file = "Table10_CI_Weibull.csv")
)) {
  x <- read_table(item$file)
  keep <- c("Setting", "Censoring", unlist(lapply(
    ci_sizes, function(size) paste(ci_measures, size)
  ), use.names = FALSE))
  if (!all(keep %in% names(x))) {
    stop("Unexpected columns in ", item$file, call. = FALSE)
  }
  utils::write.csv(x[keep], file.path(
    paper_dir, paste0("Table_", item$number, "_manuscript.csv")
  ), row.names = FALSE)
}

## Table 11: reproduce the printed columns and nonnegative displayed lower bound.
x <- utils::read.csv(
  need(file.path(run_root, "real_data", "Table11_real_data_full.csv")),
  check.names = FALSE, stringsAsFactors = FALSE
)
table11 <- data.frame(
  Dataset = x$Dataset,
  Comparison = x$Comparison,
  n1 = x$n1,
  n2 = x$n2,
  `Cens. 1 (%)` = round(x$Censoring_1_pct, 1),
  `Cens. 2 (%)` = round(x$Censoring_2_pct, 1),
  tau = round(x$tau, 1),
  `D_tau_KM` = round(x$D_tau_KM, 4),
  SE = round(x$SE_Greenwood, 4),
  `95% Normal CI` = sprintf("(%.4f, %.4f)",
                            pmax(0, x$Normal_lo), x$Normal_hi),
  `95% Bootstrap CI` = sprintf("(%.4f, %.4f)",
                               x$Percentile_lo, x$Percentile_hi),
  check.names = FALSE, stringsAsFactors = FALSE
)
utils::write.csv(table11, file.path(paper_dir, "Table_11_manuscript.csv"),
                 row.names = FALSE)

## Tables 12-14: three-decimal matrices, matching the printed manuscript.
for (index in 1:3) {
  x <- utils::read.csv(
    need(file.path(run_root, "image_data",
                   paste0("Divergence_matrix_index", index, ".csv"))),
    row.names = 1L, check.names = FALSE
  )
  utils::write.csv(round(x, 3), file.path(
    paper_dir, paste0("Table_", 11L + index, "_manuscript.csv")
  ))
}

## Copy and number the manuscript figure panels.
figure_map <- c(
  Figure_1a_KM.pdf = "Veteran_KM_curve.pdf",
  Figure_1b_cumulative.pdf = "Veteran_cumulative_divergence.pdf",
  Figure_2a_KM.pdf = "Lung_cancer_KM_curve.pdf",
  Figure_2b_cumulative.pdf = "Lung_cancer_cumulative_divergence.pdf",
  Figure_3a_KM.pdf = "GBSG2_KM_curve.pdf",
  Figure_3b_cumulative.pdf = "GBSG2_cumulative_divergence.pdf"
)
for (destination in names(figure_map)) {
  source_file <- need(file.path(run_root, "real_data", figure_map[[destination]]))
  file.copy(source_file, file.path(paper_dir, destination), overwrite = TRUE)
}
file.copy(
  need(file.path(run_root, "image_data", "Figure4_MRI_images.pdf")),
  file.path(paper_dir, "Figure_4.pdf"), overwrite = TRUE
)

outputs <- list.files(paper_dir, full.names = TRUE)
manifest <- data.frame(
  file = basename(outputs), bytes = unname(file.info(outputs)$size),
  md5 = unname(tools::md5sum(outputs)), stringsAsFactors = FALSE
)
utils::write.csv(manifest,
                 file.path(paper_dir, "PAPER_OUTPUT_MANIFEST.csv"),
                 row.names = FALSE)
cat("Manuscript outputs assembled in ", paper_dir, "\n", sep = "")
