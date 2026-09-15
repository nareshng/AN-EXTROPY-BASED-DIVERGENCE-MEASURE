## =============================================================================
## Section 5.4.1 -- point estimation of D_tau and D_CC,tau under right censoring
##
## The target tau is a fixed population quantity: the smaller of the two 80th
## lifetime percentiles. Tables 7--8 report relative mean squared error for
## each estimator against its own target:
##
##   RelMSE_D   = mean{(D_hat - D_tau)^2} / D_tau^2,
##   RelMSE_DCC = mean{(DCC_hat - D_CC,tau)^2} / D_CC,tau^2.
##
## DCC_hat uses linearly interpolated Kaplan--Meier curves and the deterministic
## one-global-crossing aggregation documented in km_functions.R. The simulated
## hazard differences change sign at most once. Since D_tau and D_CC,tau are
## different targets, their relative-MSE values are not direct efficiency
## comparisons between the two divergence measures.
##
## Tables 7--8 use (20,20), (50,50), (100,100), and (200,200). Tables 9--10
## use only (30,40), (70,50), and (100,100).
##
## Full run:
##   Rscript Point_estimation_right_censoring.R --mode=full \
##     --output-dir=results/censored_point
## Quick reproducibility check:
##   Rscript Point_estimation_right_censoring.R --mode=quick \
##     --output-dir=results/censored_point_quick
## Optional overrides: --B=2000 --seed=123
## A single positional integer is accepted as a backward-compatible B override.
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

opt <- cens_parse_cli(commandArgs(trailingOnly = TRUE), "point")
B <- as.integer(opt$B)
base_seed <- as.integer(opt$seed)
out_dir <- opt$output_dir

sample_sizes <- cens_sample_sizes("point")
tasks <- cens_build_tasks(sample_sizes)
source_md5 <- cens_source_md5("point", driver_dir)

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
  file.path(out_dir, "Section54_point_task_manifest.csv")
)

if (opt$aggregate_only) {
  cens_aggregate_checkpoints(
    "point", tasks, sample_sizes, opt$mode, B, base_seed, out_dir,
    source_md5
  )
  message("Rebuilt Tables 7--8 from validated per-cell checkpoints.")
  quit(save = "no", status = 0L)
}

selected_tasks <- cens_select_tasks(tasks, opt$task_ids)

message(sprintf("Point-estimation run: mode=%s, B=%d, cells=%d, output=%s",
                opt$mode, B, length(selected_tasks),
                normalizePath(out_dir, mustWork = TRUE)))
t0 <- Sys.time()
cens_execute_tasks(
  "point", selected_tasks, opt$mode, B, base_seed, out_dir,
  opt$resume, source_md5
)

if (is.null(opt$task_ids)) {
  cens_aggregate_checkpoints(
    "point", tasks, sample_sizes, opt$mode, B, base_seed, out_dir,
    source_md5
  )
} else {
  message(
    "Selected cells were checkpointed. Run --aggregate-only after all tasks exist."
  )
}
message(sprintf("Completed in %.2f minutes. Outputs written to %s",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                normalizePath(out_dir, mustWork = TRUE)))
