## =============================================================================
## Section 5.4.2 -- confidence intervals for fixed D_tau under right censoring
##
## Primary interval:
##   Greenwood normal interval based on the direct variance from D_km().
##
## Diagnostic comparator:
##   Ordinary within-group pairs-bootstrap percentile interval, with the same
##   fixed tau in every bootstrap resample. This interval is retained to document
##   its finite-sample behaviour; the code does not claim uniform validity near
##   the degenerate null F=G or when follow-up through tau is sparse.
##
## Exploratory intervals (long output only):
##   reflected/basic bootstrap and plug-in bias-centred normal intervals.
##
## Full run:
##   Rscript Confidence_Intervals_right_censoring.R --mode=full \
##     --output-dir=results/censored_ci
## Quick reproducibility check:
##   Rscript Confidence_Intervals_right_censoring.R --mode=quick \
##     --output-dir=results/censored_ci_quick
## Optional overrides: --B=2000 --R-boot=1000 --seed=123 --cores=1
## Positional B and R_boot are accepted for backward compatibility.
## =============================================================================

driver_args <- commandArgs(trailingOnly = FALSE)
driver_file <- gsub("~+~", " ", sub("^--file=", "", driver_args[startsWith(driver_args, "--file=")]), fixed = TRUE)  ## FIX: Rscript encodes spaces in --file= as ~+~
driver_dir <- if (length(driver_file)) {
  dirname(normalizePath(driver_file[1], mustWork = TRUE))
} else {
  getwd()
}
source(file.path(driver_dir, "atomic_io_functions.R"))
source(file.path(driver_dir, "km_functions.R"))
source(file.path(driver_dir, "censored_simulation_helpers.R"))

opt <- cens_parse_cli(commandArgs(trailingOnly = TRUE), "ci")
B <- as.integer(opt$B)
R_boot <- as.integer(opt$R_boot)
base_seed <- as.integer(opt$seed)
n_cores <- as.integer(opt$cores)
out_dir <- opt$output_dir
alpha <- 0.05
z <- stats::qnorm(1 - alpha / 2)

## The two larger balanced settings directly address the absence of visible
## convergence at n1=n2=100 in the previous revision.
sample_sizes <- cens_sample_sizes("ci")
tasks <- cens_build_tasks(sample_sizes)
source_md5 <- cens_source_md5("ci", driver_dir)

if (opt$list_tasks) {
  print(cens_task_manifest(tasks), row.names = FALSE)
  quit(save = "no", status = 0L)
}
if (opt$aggregate_only && !is.null(opt$task_ids)) {
  stop("--aggregate-only cannot be combined with --task-ids.", call. = FALSE)
}
cens_prepare_run(out_dir)

atomic_write_csv(
  cens_task_manifest(tasks),
  file.path(out_dir, "Section54_CI_task_manifest.csv")
)

if (opt$aggregate_only) {
  cens_aggregate_checkpoints(
    "ci", tasks, sample_sizes, opt$mode, B, base_seed, out_dir,
    source_md5, R_boot, alpha
  )
  message("Rebuilt Tables 9--10 from validated per-cell checkpoints.")
  quit(save = "no", status = 0L)
}

selected_tasks <- cens_select_tasks(tasks, opt$task_ids)

message(sprintf(
  "CI run: mode=%s, B=%d, R_boot=%d, cells=%d, cores=%d, output=%s",
  opt$mode,
  B,
  R_boot,
  length(selected_tasks),
  n_cores,
  normalizePath(out_dir, mustWork = TRUE)
))
t0 <- Sys.time()
cens_execute_tasks(
  "ci", selected_tasks, opt$mode, B, base_seed, out_dir,
  opt$resume, source_md5, R_boot, alpha, z, n_cores
)

if (is.null(opt$task_ids)) {
  cens_aggregate_checkpoints(
    "ci", tasks, sample_sizes, opt$mode, B, base_seed, out_dir,
    source_md5, R_boot, alpha
  )
} else {
  message(
    "Selected cells were checkpointed. Run --aggregate-only after all tasks exist."
  )
}
message(sprintf(
  "Completed in %.2f minutes. Outputs written to %s",
  as.numeric(difftime(Sys.time(), t0, units = "mins")),
  normalizePath(out_dir, mustWork = TRUE)
))
