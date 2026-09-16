#!/usr/bin/env Rscript
## =============================================================================
## Master runner for Tables 1-6.
##
##   Rscript run_uncensored_tables.R --output-dir=results --cores=4
##   Rscript run_uncensored_tables.R --output-dir=check --cores=4 --quick
##
## --quick is passed to every script; --cores only to the two confidence-interval
## scripts, since Sections 5.1 and 5.2 are single-stream by design.
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)
known <- c("--output-dir=", "--cores=", "--quick", "--quiet")
unknown <- args[!startsWith(args, "--output-dir=") & !startsWith(args, "--cores=") &
                  !args %in% c("--quick", "--quiet")]
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       "\nAccepted: ", paste(known, collapse = " "), call. = FALSE)
}

full <- commandArgs(trailingOnly = FALSE)
file_arg <- full[startsWith(full, "--file=")]
root <- if (length(file_arg)) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]),
                             fixed = TRUE), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

out_arg <- args[startsWith(args, "--output-dir=")]
output_dir <- if (length(out_arg)) sub("^--output-dir=", "", out_arg[[1L]]) else "results"
if (!grepl("^([A-Za-z]:)?[/\\\\]", output_dir, perl = TRUE)) {
  output_dir <- file.path(root, output_dir)
}
if (!dir.exists(output_dir) &&
    !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

ci_extra <- args[startsWith(args, "--cores=") | args %in% c("--quick", "--quiet")]
other_extra <- args[args %in% c("--quick")]

jobs <- list(
  list(table = "1", script = "Sim_RelativeMSE_comparison_Exponential_dist.R", ci = FALSE),
  list(table = "2", script = "Sim_RelativeMSE_comparison_Weibull_dist.R", ci = FALSE),
  list(table = "3", script = "Sim_MSE_Estimators_Exponential_dist.R", ci = FALSE),
  list(table = "4", script = "Sim_MSE_Estimators_Weibull_dist.R", ci = FALSE),
  list(table = "5", script = "Sim_Confidence_Intervals_Exponential_dist.R", ci = TRUE),
  list(table = "6", script = "Sim_Confidence_Intervals_Weibull_dist.R", ci = TRUE)
)

missing <- vapply(jobs, function(j) !file.exists(file.path(root, j$script)), logical(1))
if (any(missing)) {
  stop("Missing script(s): ",
       paste(vapply(jobs[missing], `[[`, character(1), "script"), collapse = ", "),
       call. = FALSE)
}

log_dir <- file.path(output_dir, "logs")
if (!dir.exists(log_dir) &&
    !dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create log directory: ", log_dir, call. = FALSE)
}

for (job in jobs) {
  cat("Table ", job$table, ": ", job$script, "\n", sep = "")
  started <- Sys.time()
  child <- c(paste0("--output-dir=", output_dir),
             if (job$ci) ci_extra else other_extra)
  log_file <- file.path(log_dir, sub("[.]R$", ".log", job$script))
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--no-save", "--no-restore", "--no-site-file", "--no-environ",
             shQuote(file.path(root, job$script)),
             vapply(child, shQuote, character(1L))),
    stdout = log_file, stderr = paste0(log_file, ".err")
  )
  if (!identical(as.integer(status), 0L)) {
    stop("Table ", job$table, " failed with exit status ", status,
         ". See ", log_file, ".err", call. = FALSE)
  }
  cat("  completed in ",
      sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))),
      " minutes\n", sep = "")
}

expected <- c("Table1_RelMSE_Exponential.csv", "Table2_RelMSE_Weibull.csv",
              "Table3_MSE_Exponential.csv", "Table4_MSE_Weibull.csv",
              "Table5_CI_Exponential.csv", "Table6_CI_Weibull.csv")
paths <- file.path(output_dir, expected)
info <- file.info(paths)
ok <- file.exists(paths) & !dir.exists(paths) & is.finite(info$size) & info$size > 0
if (!all(ok)) {
  stop("Run finished but output(s) are missing or empty: ",
       paste(expected[!ok], collapse = ", "), call. = FALSE)
}
cat("Tables 1-6 written to ", output_dir, "\n", sep = "")
