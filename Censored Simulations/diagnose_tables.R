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
for (f in c("Table7_RMSE_Exponential.csv", "Table9_CI_Exponential.csv", "Table10_CI_Weibull.csv")) {
  p <- file.path(out, f)
  if (!file.exists(p)) { cat(f, ": MISSING\n"); next }
  x <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  num <- names(x)[-(1:2)]
  bad <- num[vapply(x[num], function(v) !is.numeric(v) || any(!is.finite(v)), logical(1))]
  cat(f, ": ", nrow(x), " rows, ", ncol(x), " cols; size blocks: ",
      paste(unique(sub("^[^(]*", "", num)), collapse = " "), "\n", sep = "")
  if (length(bad)) cat("   non-finite/non-numeric columns: ", paste(bad, collapse = ", "), "\n", sep = "")
}
for (f in c("Section54_point_estimation_long.csv", "Section54_CI_long.csv")) {
  p <- file.path(out, f)
  if (file.exists(p)) cat(f, ": ", nrow(utils::read.csv(p)), " rows\n", sep = "")
}
