#!/usr/bin/env Rscript

## Build Tables 7--10 
## Examples:
##   Rscript Aggregate_right_censoring_results.R --analysis=point \
##     --mode=full --B=2000 --seed=5401 --output-dir=results/full/section_5_4
##   Rscript Aggregate_right_censoring_results.R --analysis=ci \
##     --mode=full --B=2000 --R-boot=1000 --seed=5402 \
##     --output-dir=results/full/section_5_4

args <- commandArgs(trailingOnly = TRUE)
analysis_arg <- args[startsWith(args, "--analysis=")]
help_requested <- any(args == "--help")
if (help_requested || length(analysis_arg) != 1L) {
  cat(paste0(
    "Usage: Rscript Aggregate_right_censoring_results.R ",
    "--analysis=point|ci [matching run options]\n",
    "Required matching options normally include --mode, --B, --seed and, ",
    "for CI, --R-boot. All cell checkpoints must be present.\n"
  ))
  quit(save = "no", status = if (help_requested) 0L else 2L)
}
analysis <- sub("^--analysis=", "", analysis_arg[[1L]])
if (!analysis %in% c("point", "ci")) {
  stop("--analysis must be point or ci.", call. = FALSE)
}
args <- args[!startsWith(args, "--analysis=")]

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[startsWith(full_args, "--file=")]
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE),  ## FIX: ~+~
                        mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
source(file.path(script_dir, "atomic_io_functions.R"))
source(file.path(script_dir, "km_functions.R"))
source(file.path(script_dir, "censored_simulation_helpers.R"))

opt <- cens_parse_cli(args, analysis)
if (!is.null(opt$task_ids)) {
  stop("Aggregation requires the complete task grid; omit --task-ids.",
       call. = FALSE)
}
sample_sizes <- cens_sample_sizes(analysis)
tasks <- cens_build_tasks(sample_sizes)
if (opt$list_tasks) {
  print(cens_task_manifest(tasks), row.names = FALSE)
  quit(save = "no", status = 0L)
}

cens_prepare_run(opt$output_dir)
source_md5 <- cens_source_md5(analysis, script_dir)
if (analysis == "point") {
  cens_aggregate_checkpoints(
    analysis, tasks, sample_sizes, opt$mode, opt$B, opt$seed,
    opt$output_dir, source_md5
  )
} else {
  cens_aggregate_checkpoints(
    analysis, tasks, sample_sizes, opt$mode, opt$B, opt$seed,
    opt$output_dir, source_md5, opt$R_boot, 0.05
  )
}
message("Aggregation completed from validated cell checkpoints only.")
