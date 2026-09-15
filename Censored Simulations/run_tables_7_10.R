#!/usr/bin/env Rscript

## Master runner for the corrected right-censoring simulations in Tables 7--10.
## Uses base R only. Run from any directory with Rscript run_tables_7_10.R.

usage <- paste0(
  "Usage: Rscript run_tables_7_10.R [options]\n",
  "  --mode=quick|full   quick software check or full paper run (default: full)\n",
  "  --output-dir=PATH   output directory (default: results/<mode>_<timestamp>)\n",
  "  --cores=N           parallel CI cells on non-Windows systems (default: 1)\n",
  "  --force             allow reuse of a nonempty compatible output directory\n",
  "  --tidy              put the four table CSVs in <output-dir>/paper_tables/ and\n",
  "                      move the supporting files into <output-dir>/details/\n",
  "  --help              print this message\n\n",
  "Full mode uses B=2000 Monte Carlo replications, R_boot=1000 bootstrap\n",
  "resamples, point seed 5401, and CI seed 5402. Quick mode uses B=50 and\n",
  "R_boot=199 and is not suitable for manuscript results.\n"
)

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[startsWith(full_args, "--file=")]
if (length(script_arg) != 1L) {
  stop("This master runner must be invoked with Rscript.", call. = FALSE)
}
root <- dirname(normalizePath(
  gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE),  ## FIX: ~+~
  winslash = "/",
  mustWork = TRUE
))

opt <- list(mode = "full", output_dir = NULL, cores = 1L, force = FALSE, tidy = FALSE)
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (identical(arg, "--help")) {
    cat(usage)
    quit(save = "no", status = 0L)
  } else if (startsWith(arg, "--mode=")) {
    opt$mode <- tolower(sub("^--mode=", "", arg))
  } else if (startsWith(arg, "--output-dir=")) {
    opt$output_dir <- sub("^--output-dir=", "", arg)
  } else if (startsWith(arg, "--cores=")) {
    cores_text <- sub("^--cores=", "", arg)
    opt$cores <- if (grepl("^[1-9][0-9]*$", cores_text)) {
      suppressWarnings(as.integer(cores_text))
    } else {
      NA_integer_
    }
  } else if (identical(arg, "--force")) {
    opt$force <- TRUE
  } else if (identical(arg, "--tidy")) {
    opt$tidy <- TRUE
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
if (.Platform$OS.type == "windows" && opt$cores > 1L) {
  warning("Windows does not support forked workers; using --cores=1.",
          call. = FALSE)
  opt$cores <- 1L
}
if (getRversion() < "4.1.0") {
  stop("R >= 4.1.0 is required; found ", getRversion(), ".",
       call. = FALSE)
}

settings <- if (opt$mode == "quick") {
  list(B = 50L, R_boot = 199L)
} else {
  list(B = 2000L, R_boot = 1000L)
}
point_seed <- 5401L
ci_seed <- 5402L

if (is.null(opt$output_dir) || !nzchar(opt$output_dir)) {
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  opt$output_dir <- file.path(root, "results", paste0(opt$mode, "_", stamp))
} else {
  opt$output_dir <- path.expand(opt$output_dir)
  if (!grepl("^([A-Za-z]:)?[/\\\\]", opt$output_dir, perl = TRUE)) {
    opt$output_dir <- file.path(root, opt$output_dir)
  }
}
output_dir <- normalizePath(
  opt$output_dir,
  winslash = "/",
  mustWork = FALSE
)
path_key <- function(path) {
  if (.Platform$OS.type == "windows") tolower(path) else path
}
output_key <- path_key(output_dir)
root_key <- path_key(root)
output_contains_source <- identical(output_key, root_key) ||
  startsWith(paste0(root_key, "/"), paste0(output_key, "/"))
if (output_contains_source || identical(dirname(output_dir), output_dir)) {
  stop(
    "--output-dir must not be the source directory, a filesystem root, or ",
    "an ancestor of the source directory.",
    call. = FALSE
  )
}

if (dir.exists(output_dir)) {
  existing <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  if (length(existing) && !opt$force) {
    stop(
      "The output directory is nonempty. Use a new directory or add --force: ",
      output_dir,
      call. = FALSE
    )
  }
}
if (!dir.exists(output_dir) &&
    !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}
obsolete_point_tables <- setdiff(
  list.files(output_dir, pattern = "^Table[78]_.*[.]csv$"),
  c("Table7_RelMSE_Exponential.csv", "Table8_RelMSE_Weibull.csv")
)
if (length(obsolete_point_tables)) {
  stop(
    "Unexpected obsolete Table 7--8 file(s) detected: ",
    paste(obsolete_point_tables, collapse = ", "),
    ". Use a fresh output directory.",
    call. = FALSE
  )
}
log_dir <- file.path(output_dir, "logs")
if (!dir.exists(log_dir) &&
    !dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create log directory: ", log_dir, call. = FALSE)
}

required_sources <- file.path(root, c(
  "Point_estimation_right_censoring.R",
  "Confidence_Intervals_right_censoring.R",
  "Aggregate_right_censoring_results.R",
  "atomic_io_functions.R",
  "censored_simulation_helpers.R",
  "km_functions.R"
))
if (any(!file.exists(required_sources))) {
  stop(
    "Missing required source file(s): ",
    paste(basename(required_sources[!file.exists(required_sources)]), collapse = ", "),
    call. = FALSE
  )
}
source(file.path(root, "atomic_io_functions.R"))

## FIX: take the sample-size design from cens_sample_sizes() instead of
## hard-coding it here, so editing that one function cannot break the contracts.
design_env <- new.env(parent = globalenv())
sys.source(file.path(root, "km_functions.R"), envir = design_env)
sys.source(file.path(root, "censored_simulation_helpers.R"), envir = design_env)
size_tags <- function(analysis) {
  vapply(design_env$cens_sample_sizes(analysis),
         function(nn) sprintf("(%d,%d)", nn[1L], nn[2L]), character(1L))
}
n_cells <- function(analysis) 36L * length(size_tags(analysis))

run_driver <- function(label, script, child_args, log_name) {
  stdout_file <- tempfile("table710_stdout_", tmpdir = log_dir, fileext = ".log")
  stderr_file <- tempfile("table710_stderr_", tmpdir = log_dir, fileext = ".log")
  on.exit(unlink(c(stdout_file, stderr_file)), add = TRUE)

  command <- file.path(R.home("bin"), "Rscript")
  args <- c(
    "--no-save",
    "--no-restore",
    "--no-site-file",
    "--no-environ",
    shQuote(script),
    vapply(child_args, shQuote, character(1L))
  )
  cat(label, "\n", sep = "")
  started <- Sys.time()
  status <- system2(
    command,
    args = args,
    stdout = stdout_file,
    stderr = stderr_file
  )
  output <- c(
    if (file.exists(stdout_file)) readLines(stdout_file, warn = FALSE) else character(),
    "----- STDERR -----",
    if (file.exists(stderr_file)) readLines(stderr_file, warn = FALSE) else character()
  )
  writeLines(output, file.path(log_dir, log_name), useBytes = TRUE)
  if (!identical(as.integer(status), 0L)) {
    stop(
      label,
      " failed with exit status ",
      status,
      ". See ",
      file.path(log_dir, log_name),
      call. = FALSE
    )
  }
  cat(
    "  completed in ",
    sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))),
    " minutes\n",
    sep = ""
  )
  invisible(TRUE)
}

common <- c(
  paste0("--mode=", opt$mode),
  paste0("--B=", settings$B),
  paste0("--output-dir=", output_dir)
)
run_driver(
  "Generating Tables 7--8 ...",
  file.path(root, "Point_estimation_right_censoring.R"),
  c(common, paste0("--seed=", point_seed)),
  "tables_7_8.log"
)
run_driver(
  "Generating Tables 9--10 ...",
  file.path(root, "Confidence_Intervals_right_censoring.R"),
  c(
    common,
    paste0("--R-boot=", settings$R_boot),
    paste0("--seed=", ci_seed),
    paste0("--cores=", opt$cores)
  ),
  "tables_9_10.log"
)

expected <- c(
  "Table7_RelMSE_Exponential.csv",
  "Table8_RelMSE_Weibull.csv",
  "Table9_CI_Exponential.csv",
  "Table10_CI_Weibull.csv",
  "Section54_point_estimation_long.csv",
  "Section54_point_estimation_replications.csv",
  "Section54_CI_long.csv",
  "Section54_CI_replications.csv"
)
expected_paths <- file.path(output_dir, expected)
output_info <- file.info(expected_paths)
valid_output <- file.exists(expected_paths) & !dir.exists(expected_paths) &
  is.finite(output_info$size) & output_info$size > 0
if (!all(valid_output)) {
  stop(
    "Run completed but expected output(s) are missing or empty: ",
    paste(expected[!valid_output], collapse = ", "),
    call. = FALSE
  )
}

read_output_csv <- function(name) {
  utils::read.csv(
    file.path(output_dir, name),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

main_names <- expected[seq_len(4L)]
tables <- lapply(main_names, read_output_csv)
point_long <- read_output_csv("Section54_point_estimation_long.csv")
ci_long <- read_output_csv("Section54_CI_long.csv")

point_sizes <- size_tags("point")
ci_sizes <- size_tags("ci")
if (!identical(
      point_sizes,
      c("(20,20)", "(50,50)", "(100,100)", "(200,200)")
    )) {
  stop("Tables 7--8 sample-size design was changed unexpectedly.", call. = FALSE)
}
if (!identical(ci_sizes, c("(30,40)", "(70,50)", "(100,100)"))) {
  stop("Tables 9--10 sample-size design was changed unexpectedly.", call. = FALSE)
}
point_columns <- c(
  "Setting", "Censoring", "D_tau", "DCC_tau",
  unlist(lapply(point_sizes, function(nn) {
    paste(c("D_KM", "DCC_KM"), nn)
  }), use.names = FALSE)
)
ci_columns <- c(
  "Setting", "Censoring", "D_tau",
  unlist(lapply(ci_sizes, function(nn) {
    paste(
      c(
        "Normal_CP", "Normal_MCSE", "Normal_AL", "Bootstrap_CP",
        "Bootstrap_MCSE", "Bootstrap_AL", "Valid_%"
      ),
      nn
    )
  }), use.names = FALSE)
)

## FIX: report WHICH part of the contract failed, instead of one generic message.
validate_table <- function(x, columns, label) {
  problems <- character(0)
  if (nrow(x) != 18L) {
    problems <- c(problems, sprintf("expected 18 rows, found %d", nrow(x)))
  }
  if (!identical(names(x), columns)) {
    missing_cols <- setdiff(columns, names(x))
    extra_cols <- setdiff(names(x), columns)
    if (length(missing_cols)) {
      problems <- c(problems, paste0(
        "column(s) expected but absent: ", paste(missing_cols, collapse = ", ")
      ))
    }
    if (length(extra_cols)) {
      problems <- c(problems, paste0(
        "column(s) present but not expected: ", paste(extra_cols, collapse = ", ")
      ))
    }
    if (!length(missing_cols) && !length(extra_cols)) {
      problems <- c(problems, "columns are in a different order than expected")
    }
  }
  if (anyDuplicated(x[c("Setting", "Censoring")]) > 0L) {
    problems <- c(problems, "duplicated Setting/Censoring keys")
  }
  numeric_columns <- intersect(setdiff(columns, c("Setting", "Censoring")), names(x))
  bad_numeric <- numeric_columns[vapply(x[numeric_columns], function(value) {
    !is.numeric(value) || any(!is.finite(value))
  }, logical(1L))]
  if (length(bad_numeric)) {
    problems <- c(problems, paste0(
      "non-numeric or non-finite value(s) in: ",
      paste(bad_numeric, collapse = ", ")
    ))
  }
  if (length(problems)) {
    stop(label, " failed its publication-table contract:\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  }
}
validate_table(tables[[1L]], point_columns, "Table 7")
validate_table(tables[[2L]], point_columns, "Table 8")
validate_table(tables[[3L]], ci_columns, "Table 9")
validate_table(tables[[4L]], ci_columns, "Table 10")

point_required <- c(
  "B_requested", "failure_rate", "D_tau", "DCC_tau",
  "RelMSE_D", "RelMSE_DCC"
)
ci_required <- c(
  "B_requested", "R_boot", "failure_rate", "bootstrap_failure_rate"
)
if (!all(point_required %in% names(point_long)) ||
    nrow(point_long) != n_cells("point") ||
    any(!is.finite(point_long$B_requested)) ||
    any(point_long$B_requested != settings$B) ||
    any(!is.finite(point_long$failure_rate)) ||
    any(!vapply(point_long[c(
      "D_tau", "DCC_tau", "RelMSE_D", "RelMSE_DCC"
    )], function(value) is.numeric(value) && all(is.finite(value)), logical(1L))) ||
    any(point_long$failure_rate > 0)) {
  stop("Tables 7--8 summary contract failed.", call. = FALSE)
}
unexpected_point_tables <- setdiff(
  list.files(output_dir, pattern = "^Table[78]_.*[.]csv$"),
  expected[seq_len(2L)]
)
if (length(unexpected_point_tables)) {
  stop(
    "Unexpected obsolete Table 7--8 file(s) detected: ",
    paste(unexpected_point_tables, collapse = ", "),
    ". Use a fresh output directory.",
    call. = FALSE
  )
}
if (!all(ci_required %in% names(ci_long)) ||
    nrow(ci_long) != n_cells("ci") ||
    any(!is.finite(ci_long$B_requested)) ||
    any(!is.finite(ci_long$R_boot)) ||
    any(ci_long$B_requested != settings$B) ||
    any(ci_long$R_boot != settings$R_boot) ||
    any(!is.finite(ci_long$failure_rate)) ||
    any(!is.finite(ci_long$bootstrap_failure_rate)) ||
    any(ci_long$failure_rate > 0) ||
    any(ci_long$bootstrap_failure_rate > 0)) {
  stop("Tables 9--10 summary contract failed.", call. = FALSE)
}

main_paths <- file.path(output_dir, main_names)
manifest <- data.frame(
  table = 7:10,
  file = main_names,
  rows = vapply(tables, nrow, integer(1L)),
  columns = vapply(tables, ncol, integer(1L)),
  bytes = unname(file.info(main_paths)$size),
  md5 = unname(tools::md5sum(main_paths)),
  stringsAsFactors = FALSE
)
atomic_write_csv(
  manifest,
  file.path(output_dir, "Tables7_10_manifest.csv")
)

settings_text <- c(
  paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Mode:", opt$mode),
  paste("Monte Carlo replications B:", settings$B),
  paste("Bootstrap resamples R_boot:", settings$R_boot),
  paste("Point-estimation seed:", point_seed),
  paste("Confidence-interval seed:", ci_seed),
  paste("Requested CI cores:", opt$cores),
  paste("R version:", R.version.string),
  paste("RNG: Mersenne-Twister; Inversion; Rejection")
)
writeLines(
  settings_text,
  file.path(output_dir, "Tables7_10_run_settings.txt"),
  useBytes = TRUE
)
atomic_save_rds(
  list(
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    mode = opt$mode,
    B = settings$B,
    R_boot = settings$R_boot,
    point_seed = point_seed,
    ci_seed = ci_seed,
    cores = opt$cores,
    manifest = manifest,
    R_version = R.version.string,
    platform = R.version$platform
  ),
  file.path(output_dir, "Tables7_10_run_metadata.rds")
)

## FIX/OPTION: --tidy separates the four manuscript tables from the supporting
## output. Nothing is deleted; supporting files are only moved into details/.
if (isTRUE(opt$tidy)) {
  tables_dir <- file.path(output_dir, "paper_tables")
  details_dir <- file.path(output_dir, "details")
  dir.create(tables_dir, showWarnings = FALSE)
  dir.create(details_dir, showWarnings = FALSE)
  keep_at_top <- c("paper_tables", "details", "logs")
  for (nm in setdiff(list.files(output_dir, all.files = TRUE, no.. = TRUE), keep_at_top)) {
    file.rename(file.path(output_dir, nm), file.path(details_dir, nm))
  }
  copied <- file.copy(file.path(details_dir, main_names),
                      file.path(tables_dir, main_names), overwrite = TRUE)
  if (!all(copied)) stop("Could not assemble paper_tables/.", call. = FALSE)
  cat("Manuscript tables: ", tables_dir, "\n", sep = "")
  cat("Supporting output:  ", details_dir, "\n", sep = "")
}

cat("Tables 7--10 were generated successfully.\n")
cat("Output directory: ", output_dir, "\n", sep = "")
