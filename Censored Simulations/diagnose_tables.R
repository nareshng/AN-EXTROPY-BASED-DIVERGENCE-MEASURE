## Standalone diagnostic: run from the package folder with the failing output dir.
##   Rscript diagnose_tables.R results/quick_20260915T180654Z
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Give the output directory of the failed run.", call. = FALSE)
out <- args[1]
env <- new.env(parent = globalenv())
sys.source("km_functions.R", envir = env)
sys.source("censored_simulation_helpers.R", envir = env)
tag <- function(a) vapply(env$cens_sample_sizes(a), function(nn) sprintf("(%d,%d)", nn[1], nn[2]), character(1))
cat("cens_sample_sizes(\"point\"): ", paste(tag("point"), collapse = " "), "\n", sep = "")
cat("cens_sample_sizes(\"ci\")   : ", paste(tag("ci"), collapse = " "), "\n\n", sep = "")
runner <- paste(readLines("run_tables_7_10.R", warn = FALSE), collapse = "\n")
for (nm in c("point_sizes", "ci_sizes")) {
  hit <- regmatches(runner, regexpr(paste0("\n", nm, " <- [^\n]*"), runner))
  cat("run_tables_7_10.R ", trimws(hit), "\n", sep = "")
}
for (nm in c("point_long", "ci_long")) {
  hit <- regmatches(runner, regexpr(paste0("nrow\\(", nm, "\\) != [^ \n|]*"), runner))
  cat("run_tables_7_10.R  ", if (length(hit)) hit else "(row-count check not found)", "\n", sep = "")
}
cat("\n")
point_sizes <- tag("point")
ci_sizes <- tag("ci")
if (!identical(
      point_sizes,
      c("(20,20)", "(50,50)", "(100,100)", "(200,200)")
    )) {
  stop("Tables 7--8 sample-size design is incorrect.", call. = FALSE)
}
if (!identical(ci_sizes, c("(30,40)", "(70,50)", "(100,100)"))) {
  stop("Tables 9--10 sample-size design is incorrect.", call. = FALSE)
}
point_columns <- c(
  "Setting", "Censoring", "D_tau", "DCC_tau",
  unlist(lapply(point_sizes, function(nn) {
    c(paste("D_KM", nn), paste("DCC_KM", nn))
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
point_files <- c(
  "Table7_RelMSE_Exponential.csv",
  "Table8_RelMSE_Weibull.csv"
)
all_files <- c(
  point_files,
  "Table9_CI_Exponential.csv",
  "Table10_CI_Weibull.csv"
)
problems <- character(0)

for (f in all_files) {
  p <- file.path(out, f)
  if (!file.exists(p)) {
    problems <- c(problems, paste(f, "is missing"))
    next
  }
  x <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  num <- names(x)[-(1:2)]
  bad <- num[vapply(x[num], function(v) !is.numeric(v) || any(!is.finite(v)), logical(1))]
  cat(f, ": ", nrow(x), " rows, ", ncol(x), " cols; size blocks: ",
      paste(unique(sub("^[^(]*", "", num)), collapse = " "), "\n", sep = "")
  if (nrow(x) != 18L) {
    problems <- c(problems, sprintf("%s has %d rows instead of 18", f, nrow(x)))
  }
  if (length(bad)) {
    problems <- c(problems, paste0(
      f, " has invalid numeric column(s): ", paste(bad, collapse = ", ")
    ))
  }
  expected_columns <- if (f %in% point_files) point_columns else ci_columns
  if (!identical(names(x), expected_columns)) {
    problems <- c(problems, paste(f, "has an incorrect column schema"))
  }
}
unexpected_point_tables <- setdiff(
  list.files(out, pattern = "^Table[78]_.*[.]csv$"),
  point_files
)
if (length(unexpected_point_tables)) {
  problems <- c(problems, paste0(
    "obsolete Table 7--8 file(s): ",
    paste(unexpected_point_tables, collapse = ", ")
  ))
}
for (f in c("Section54_point_estimation_long.csv", "Section54_CI_long.csv")) {
  p <- file.path(out, f)
  if (file.exists(p)) cat(f, ": ", nrow(utils::read.csv(p)), " rows\n", sep = "")
}
if (length(problems)) {
  stop(
    "Table diagnostics failed:\n  - ",
    paste(problems, collapse = "\n  - "),
    call. = FALSE
  )
}
cat("All four publication-table contracts passed.\n")
