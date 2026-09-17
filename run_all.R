#!/usr/bin/env Rscript

## Master workflow for all manuscript tables and figures.
## This script never installs packages and never overwrites an existing run.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(mode = "quick", cores = 1L, seed_offset = 0L,
            output_dir = NULL)
for (arg in args) {
  if (startsWith(arg, "--mode=")) {
    opt$mode <- sub("^--mode=", "", arg)
  } else if (startsWith(arg, "--cores=")) {
    opt$cores <- suppressWarnings(as.integer(sub("^--cores=", "", arg)))
  } else if (startsWith(arg, "--seed-offset=")) {
    opt$seed_offset <- suppressWarnings(as.integer(
      sub("^--seed-offset=", "", arg)
    ))
  } else if (startsWith(arg, "--output-dir=")) {
    opt$output_dir <- sub("^--output-dir=", "", arg)
  } else if (arg %in% c("--help", "-h")) {
    cat(
      "Usage: Rscript --vanilla run_all.R ",
      "--mode=quick|full --cores=N --seed-offset=N ",
      "[--output-dir=PATH]\n",
      sep = ""
    )
    quit(save = "no", status = 0L)
  } else {
    stop("Unknown option: ", arg, call. = FALSE)
  }
}

if (!opt$mode %in% c("quick", "full")) {
  stop("--mode must be quick or full.", call. = FALSE)
}
if (length(opt$cores) != 1L || is.na(opt$cores) || opt$cores < 1L) {
  stop("--cores must be a positive integer.", call. = FALSE)
}
if (length(opt$seed_offset) != 1L || is.na(opt$seed_offset) ||
    opt$seed_offset < 0L) {
  stop("--seed-offset must be a nonnegative integer.", call. = FALSE)
}
if (.Platform$OS.type == "windows" && opt$cores > 1L) {
  warning("Windows uses one core for this workflow.", call. = FALSE)
  opt$cores <- 1L
}

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[startsWith(all_args, "--file=")]
if (length(script_arg) != 1L) {
  stop("Invoke this file with Rscript.", call. = FALSE)
}
root <- dirname(normalizePath(
  gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE),
  winslash = "/", mustWork = TRUE
))

stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
run_root <- if (is.null(opt$output_dir)) {
  file.path(root, "results", paste0(opt$mode, "_", stamp))
} else if (grepl("^([A-Za-z]:)?[/\\\\]", opt$output_dir, perl = TRUE)) {
  opt$output_dir
} else {
  file.path(root, opt$output_dir)
}
run_root <- normalizePath(run_root, winslash = "/", mustWork = FALSE)
if (dir.exists(run_root) && length(list.files(run_root, all.files = TRUE,
                                              no.. = TRUE))) {
  stop("Output directory is not empty: ", run_root, call. = FALSE)
}
if (!dir.exists(run_root) &&
    !dir.create(run_root, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create output directory: ", run_root, call. = FALSE)
}

rscript <- file.path(R.home("bin"), "Rscript")
log_dir <- file.path(run_root, "logs")
dir.create(log_dir, showWarnings = FALSE)

run_step <- function(label, script, step_args = character()) {
  if (!file.exists(script)) {
    stop("Missing script for ", label, ": ", script, call. = FALSE)
  }
  log_stem <- file.path(log_dir, gsub("[^A-Za-z0-9]+", "_", label))
  stdout_file <- paste0(log_stem, ".stdout.log")
  stderr_file <- paste0(log_stem, ".stderr.log")
  command_args <- c("--vanilla", shQuote(script), vapply(step_args, shQuote,
                                                         character(1L)))
  started <- Sys.time()
  status <- system2(rscript, command_args, stdout = stdout_file,
                    stderr = stderr_file)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  record <- data.frame(
    step = label, status = as.integer(status), elapsed_seconds = elapsed,
    stdout_log = stdout_file, stderr_log = stderr_file,
    stringsAsFactors = FALSE
  )
  if (!identical(as.integer(status), 0L)) {
    stop(label, " failed. See ", stdout_file, " and ", stderr_file,
         call. = FALSE)
  }
  record
}

records <- list()
records[[length(records) + 1L]] <- run_step(
  "capture_environment",
  file.path(root, "capture_environment.R"),
  c("--strict", paste0("--output-dir=", file.path(run_root, "environment")))
)

records[[length(records) + 1L]] <- run_step(
  "uncensored_tables_1_6",
  file.path(root, "Uncensored Simulations", "run_tables_1_6.R"),
  c(paste0("--mode=", opt$mode),
    paste0("--seed-offset=", opt$seed_offset),
    paste0("--output-dir=", file.path(run_root, "uncensored")))
)

records[[length(records) + 1L]] <- run_step(
  "censored_tables_7_10",
  file.path(root, "Censored Simulations", "run_tables_7_10.R"),
  c(paste0("--mode=", opt$mode), paste0("--cores=", opt$cores), "--tidy",
    paste0("--seed-offset=", opt$seed_offset),
    paste0("--output-dir=", file.path(run_root, "censored")))
)

real_boot <- if (opt$mode == "full") 5000L else 199L
real_seed <- 2026L + opt$seed_offset
records[[length(records) + 1L]] <- run_step(
  "real_data_table_11_figures_1_3",
  file.path(root, "Real Data Analysis", "Censored_Real_Data_Analysis.R"),
  c(paste0("--bootstrap-reps=", real_boot), paste0("--seed=", real_seed),
    "--gbsg2-source=TH.data",
    paste0("--output-dir=", file.path(run_root, "real_data")))
)

image_manifest <- file.path(
  root, "Real Data Analysis", "Images", "image_manifest.csv"
)
if (!file.exists(image_manifest)) {
  stop("Complete and add the required image manifest: ", image_manifest,
       call. = FALSE)
}
records[[length(records) + 1L]] <- run_step(
  "image_tables_12_14_figure_4",
  file.path(root, "Real Data Analysis", "Image_based_Real_Data_Analysis.R"),
  c(paste0("--manifest=", image_manifest),
    paste0("--output-dir=", file.path(run_root, "image_data")))
)

records[[length(records) + 1L]] <- run_step(
  "assemble_manuscript_outputs",
  file.path(root, "make_manuscript_tables.R"),
  paste0("--run-root=", run_root)
)

run_log <- do.call(rbind, records)
utils::write.csv(run_log, file.path(run_root, "MASTER_RUN_LOG.csv"),
                 row.names = FALSE)

writeLines(c(
  paste("Mode:", opt$mode),
  paste("Cores:", opt$cores),
  paste("Seed offset:", opt$seed_offset),
  paste("Real-data base seed:", real_seed),
  paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(run_root, "MASTER_RUN_SETTINGS.txt"), useBytes = TRUE)

files <- list.files(run_root, recursive = TRUE, full.names = TRUE,
                    all.files = TRUE, no.. = TRUE)
files <- files[file.exists(files) & !dir.exists(files)]
manifest <- data.frame(
  file = substring(files, nchar(run_root) + 2L),
  bytes = unname(file.info(files)$size),
  md5 = unname(tools::md5sum(files)),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(run_root, "MASTER_FILE_MANIFEST.csv"),
                 row.names = FALSE)

cat("Workflow completed successfully.\n")
cat("Run directory: ", run_root, "\n", sep = "")
