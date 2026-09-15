## Shared infrastructure for the censored point-estimation and CI simulations.
## This file assumes that km_functions.R has already been sourced.

## Input: "point" or "ci". Output: the command-line help text, invisibly.
cens_print_usage <- function(analysis) {
  analysis <- match.arg(analysis, c("point", "ci"))
  if (analysis == "point") {
    txt <- paste0(
      "Usage: Rscript Point_estimation_right_censoring.R [options]\n",
      "  --mode=quick|full       quick: B=50; full: B=2000 (default)\n",
      "  --B=N                   Monte Carlo replications (2 <= N <= 1000003)\n",
      "  --seed=N                positive integer base seed (default 123)\n",
      "  --task-ids=LIST         exact cell IDs, for example 1,4,9-12\n",
      "  --list-tasks            print the deterministic cell map and exit\n",
      "  --aggregate-only        rebuild outputs from saved cell checkpoints\n",
      "  --resume                reuse compatible cell checkpoints (default)\n",
      "  --no-resume             recompute and atomically replace selected cells\n",
      "  --output-dir=PATH       output directory\n",
      "  --help                  print this message\n"
    )
  } else {
    txt <- paste0(
      "Usage: Rscript Confidence_Intervals_right_censoring.R [options]\n",
      "  --mode=quick|full       quick: B=50, R_boot=199; ",
      "full: 2000, 1000 (default)\n",
      "  --B=N                   Monte Carlo replications (2 <= N <= 1000003)\n",
      "  --R-boot=N              bootstrap resamples per replication (N >= 19)\n",
      "  --seed=N                positive integer base seed (default 123)\n",
      "  --cores=N               forked workers on non-Windows systems (default 1)\n",
      "  --task-ids=LIST         exact cell IDs, for example 1,4,9-12\n",
      "  --list-tasks            print the deterministic cell map and exit\n",
      "  --aggregate-only        rebuild outputs from saved cell checkpoints\n",
      "  --resume                reuse compatible cell checkpoints (default)\n",
      "  --no-resume             recompute and atomically replace selected cells\n",
      "  --output-dir=PATH       output directory\n",
      "  --help                  print this message\n"
    )
  }
  cat(txt)
  invisible(txt)
}

## Input: command-line character vector and analysis name. Output: validated options list.
cens_parse_cli <- function(args, analysis) {
  analysis <- match.arg(analysis, c("point", "ci"))
  if (any(args == "--help")) {
    cens_print_usage(analysis)
    quit(save = "no", status = 0)
  }

  opt <- list(
    mode = "full",
    B = NULL,
    R_boot = NULL,
    seed = 123L,
    cores = 1L,
    output_dir = NULL,
    task_ids = NULL,
    list_tasks = FALSE,
    aggregate_only = FALSE,
    resume = TRUE
  )
  positional <- character(0)
  for (arg in args) {
    if (startsWith(arg, "--mode=")) {
      opt$mode <- substring(arg, 8L)
    } else if (startsWith(arg, "--B=")) {
      opt$B <- suppressWarnings(as.integer(substring(arg, 5L)))
    } else if (analysis == "ci" && startsWith(arg, "--R-boot=")) {
      opt$R_boot <- suppressWarnings(as.integer(substring(arg, 10L)))
    } else if (startsWith(arg, "--seed=")) {
      opt$seed <- suppressWarnings(as.integer(substring(arg, 8L)))
    } else if (analysis == "ci" && startsWith(arg, "--cores=")) {
      opt$cores <- suppressWarnings(as.integer(substring(arg, 9L)))
    } else if (startsWith(arg, "--task-ids=")) {
      opt$task_ids <- cens_parse_task_ids(substring(arg, 12L))
    } else if (identical(arg, "--list-tasks")) {
      opt$list_tasks <- TRUE
    } else if (identical(arg, "--aggregate-only")) {
      opt$aggregate_only <- TRUE
    } else if (identical(arg, "--resume")) {
      opt$resume <- TRUE
    } else if (identical(arg, "--no-resume")) {
      opt$resume <- FALSE
    } else if (startsWith(arg, "--output-dir=")) {
      opt$output_dir <- substring(arg, 14L)
    } else if (startsWith(arg, "--")) {
      stop(sprintf("Unknown option: %s", arg), call. = FALSE)
    } else {
      positional <- c(positional, arg)
    }
  }

  if (!opt$mode %in% c("quick", "full")) {
    stop("--mode must be quick or full.", call. = FALSE)
  }
  if (analysis == "point" &&
      (length(positional) > 1L ||
       (length(positional) && !is.null(opt$B)))) {
    stop("Supply B either positionally or through --B, not both.", call. = FALSE)
  }
  if (analysis == "ci" && length(positional) > 2L) {
    stop("At most two positional arguments are allowed.", call. = FALSE)
  }
  if (length(positional) >= 1L) {
    if (analysis == "ci" && !is.null(opt$B)) {
      stop("Supply B either positionally or through --B, not both.", call. = FALSE)
    }
    opt$B <- suppressWarnings(as.integer(positional[1]))
  }
  if (length(positional) == 2L) {
    if (!is.null(opt$R_boot)) {
      msg <- "Supply R_boot either positionally or through --R-boot, not both."
      stop(msg, call. = FALSE)
    }
    opt$R_boot <- suppressWarnings(as.integer(positional[2]))
  }

  if (is.null(opt$B)) {
    opt$B <- if (opt$mode == "quick") 50L else 2000L
  }
  if (analysis == "ci" && is.null(opt$R_boot)) {
    opt$R_boot <- if (opt$mode == "quick") 199L else 1000L
  }
  if (is.null(opt$output_dir) || !nzchar(opt$output_dir)) {
    suffix <- switch(
      analysis,
      point = if (opt$mode == "quick") "censored_point_quick" else "censored_point",
      ci = if (opt$mode == "quick") "censored_ci_quick" else "censored_ci"
    )
    opt$output_dir <- file.path(getwd(), "results", suffix)
  }

  ## The cell-specific seed spacing is 1,000,003. Keeping B at or below that
  ## value guarantees distinct seeds across every task/replication pair.
  assert_scalar_numeric(
    opt$B,
    "B",
    lower = 2,
    upper = 1000003,
    integer = TRUE
  )
  if (analysis == "ci") {
    assert_scalar_numeric(opt$R_boot, "R_boot", lower = 19, integer = TRUE)
  }
  assert_scalar_numeric(
    opt$seed,
    "seed",
    lower = 1,
    upper = 2147483646,
    integer = TRUE
  )
  if (analysis == "ci") {
    assert_scalar_numeric(opt$cores, "cores", lower = 1, integer = TRUE)
  }
  opt
}

## Input: comma-separated positive integers or inclusive ranges.
## Output: sorted unique task IDs, or NULL for an empty specification.
cens_parse_task_ids <- function(value) {
  if (length(value) != 1L || is.na(value)) {
    stop("--task-ids must be one comma-separated value.", call. = FALSE)
  }
  value <- trimws(value)
  if (!nzchar(value)) return(NULL)
  tokens <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (any(!nzchar(tokens))) {
    stop("--task-ids contains an empty item.", call. = FALSE)
  }
  expand_one <- function(token) {
    if (grepl("^[0-9]+$", token)) return(as.integer(token))
    if (!grepl("^[0-9]+-[0-9]+$", token)) {
      stop("Invalid --task-ids item: ", token, call. = FALSE)
    }
    bounds <- as.integer(strsplit(token, "-", fixed = TRUE)[[1L]])
    if (bounds[[1L]] > bounds[[2L]]) {
      stop("Descending --task-ids range is not allowed: ", token,
           call. = FALSE)
    }
    seq.int(bounds[[1L]], bounds[[2L]])
  }
  ids <- unlist(lapply(tokens, expand_one), use.names = FALSE)
  if (!length(ids) || anyNA(ids) || any(ids < 1L)) {
    stop("--task-ids must contain positive integers.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop("--task-ids contains duplicate cell IDs.", call. = FALSE)
  }
  sort(as.integer(ids))
}

## Input: output-directory path. Output: normalized writable path, invisibly.
cens_prepare_run <- function(output_dir) {
  if (!dir.exists(output_dir) &&
      !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop(sprintf("Could not create output directory: %s", output_dir), call. = FALSE)
  }
  if (file.access(output_dir, 2L) != 0L) {
    stop("Output directory is not writable.", call. = FALSE)
  }
  if (!isTRUE(verify_km_functions(verbose = TRUE))) {
    stop("KM implementation self-check failed; simulation was not started.",
         call. = FALSE)
  }
  do.call(RNGkind, as.list(cens_rng_signature()))
  invisible(normalizePath(output_dir, mustWork = TRUE))
}

## Input: none. Output: named RNG algorithms used by every censored simulation.
cens_rng_signature <- function() {
  c(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}

## Input: none. Output: shared censoring levels and named lifetime configurations.
cens_scenario_design <- function() {
  list(
    censoring_levels = c(0.10, 0.30, 0.50),
    exp_configs = list(
      E1 = c(0.5, 1.0),
      E2 = c(1.0, 2.0),
      E3 = c(2.0, 5.0),
      E4 = c(1.0, 0.5),
      E5 = c(5.0, 2.0),
      E6 = c(10.0, 5.0)
    ),
    weibull_configs = list(
      W1 = c(0.5, 1.0, 1.0, 1.0),
      W2 = c(1.0, 1.0, 2.0, 1.0),
      W3 = c(2.0, 1.0, 1.0, 1.0),
      W4 = c(1.5, 1.0, 3.0, 1.5),
      W5 = c(0.7, 2.0, 1.5, 1.0),
      W6 = c(3.0, 1.0, 1.2, 2.0)
    )
  )
}

## Input: "point" or "ci". Output: canonical manuscript sample-size pairs.
cens_sample_sizes <- function(analysis) {
  analysis <- match.arg(analysis, c("point", "ci"))
  if (analysis == "point") {
    list(c(20L, 20L), c(50L, 50L), c(100L, 100L), c(200L, 200L))
  } else {
    ## Tables 9-10 sample-size pairs (manuscript layout: the first three).
    ## Add c(200L, 200L), c(500L, 500L) here for the convergence rows.
    list(
      c(30L, 40L), c(70L, 50L), c(100L, 100L)
    )
  }
}

## Input: seed, cell index, replication index. Output: deterministic valid R seed.
cens_safe_seed <- function(base_seed, cell_id, replication) {
  modulus <- 2147483646
  as.integer(
    (as.double(base_seed) + 1000003 * cell_id + replication) %% modulus + 1
  )
}

## Input: numeric vector. Output: mean over finite values, or NA.
cens_safe_mean <- function(x) {
  if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
}

## Input: numeric vector. Output: SD over finite values when at least two exist.
cens_safe_sd <- function(x) {
  if (sum(is.finite(x)) >= 2L) stats::sd(x[is.finite(x)]) else NA_real_
}

## Input: numeric vector and probability. Output: type-8 finite-value quantile, or NA.
cens_safe_quantile <- function(x, p) {
  if (!any(is.finite(x))) return(NA_real_)
  as.numeric(
    stats::quantile(x[is.finite(x)], p, names = FALSE, type = 8)
  )
}

## Input: logical flags and inclusion mask. Output: conditional TRUE rate, or NA.
cens_safe_flag_rate <- function(flag, keep) {
  if (any(keep)) mean(flag[keep] %in% TRUE) else NA_real_
}

## Input: estimates, scalar truth, mask, RelMSE flag. Output: named accuracy vector.
cens_accuracy <- function(est, truth, keep, include_relmse = FALSE) {
  use <- keep & is.finite(est)
  x <- est[use]
  if (!length(x)) {
    ans <- c(
      n = 0,
      mean = NA,
      bias = NA,
      sd = NA,
      mcse_mean = NA,
      mse = NA,
      rmse = NA
    )
    if (include_relmse) ans <- c(ans, relmse = NA)
    return(ans)
  }
  mse <- mean((x - truth)^2)
  sx <- if (length(x) >= 2L) stats::sd(x) else NA_real_
  ans <- c(
    n = length(x),
    mean = mean(x),
    bias = mean(x) - truth,
    sd = sx,
    mcse_mean = if (is.finite(sx)) sx / sqrt(length(x)) else NA_real_,
    mse = mse,
    rmse = sqrt(mse)
  )
  if (include_relmse) {
    relmse <- if (truth^2 > .Machine$double.eps) mse / truth^2 else NA_real_
    ans <- c(ans, relmse = relmse)
  }
  ans
}

## Input: interval endpoints, truth, attempt mask. Output: named coverage diagnostics.
cens_interval_metrics <- function(lo, hi, truth, attempted) {
  valid <- attempted & is.finite(lo) & is.finite(hi) & lo <= hi
  covered <- valid & lo <= truth & truth <= hi
  n_valid <- sum(valid)
  conditional_p <- if (n_valid) mean(covered[valid]) else NA_real_
  all_attempt_p <- sum(covered) / length(attempted)
  c(
    n_valid = n_valid,
    availability = mean(valid),
    cp = 100 * conditional_p,
    cp_mcse = if (n_valid && is.finite(conditional_p)) {
      100 * sqrt(conditional_p * (1 - conditional_p) / n_valid)
    } else {
      NA_real_
    },
    cp_all_attempts = 100 * all_attempt_p,
    cp_all_attempts_mcse = 100 * sqrt(
      all_attempt_p * (1 - all_attempt_p) / length(attempted)
    ),
    average_length = if (n_valid) mean(hi[valid] - lo[valid]) else NA_real_
  )
}

## Input: sample sizes and shared design lists. Output: ordered simulation task list.
cens_build_tasks <- function(sample_sizes, design = cens_scenario_design()) {
  tasks <- list()
  cell_id <- 0L
  for (family in c("exp", "weibull")) {
    configs <- if (family == "exp") {
      design$exp_configs
    } else {
      design$weibull_configs
    }
    for (cfg_name in names(configs)) {
      for (cens in design$censoring_levels) {
        for (nn in sample_sizes) {
          cell_id <- cell_id + 1L
          tasks[[cell_id]] <- list(
            cell_id = cell_id,
            family = family,
            cfg_name = cfg_name,
            cfg = configs[[cfg_name]],
            cens = cens,
            n1 = nn[1],
            n2 = nn[2]
          )
        }
      }
    }
  }
  tasks
}

## Input: ordered task list. Output: one row per deterministic simulation cell.
cens_task_manifest <- function(tasks) {
  rows <- lapply(tasks, function(task) {
    values <- paste(
      formatC(task$cfg, digits = 17L, format = "g"),
      collapse = ","
    )
    parameterization <- if (task$family == "exp") {
      "rate1,rate2"
    } else {
      "shape1,scale1,shape2,scale2"
    }
    data.frame(
      task_id = task$cell_id,
      family = task$family,
      config = task$cfg_name,
      target_censoring = task$cens,
      n1 = task$n1,
      n2 = task$n2,
      parameterization = parameterization,
      parameter_values = values,
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  ans
}

## Inputs: ordered task list and optional IDs. Output: exact selected task list.
cens_select_tasks <- function(tasks, task_ids = NULL) {
  if (is.null(task_ids)) return(tasks)
  available <- vapply(tasks, `[[`, integer(1), "cell_id")
  unknown <- setdiff(task_ids, available)
  if (length(unknown)) {
    stop(
      "Unknown task ID(s): ", paste(unknown, collapse = ", "),
      ". Use --list-tasks to inspect the deterministic map.",
      call. = FALSE
    )
  }
  tasks[match(task_ids, available)]
}

## Inputs: analysis name and driver directory. Output: named source MD5 vector.
cens_source_md5 <- function(analysis, driver_dir) {
  analysis <- match.arg(analysis, c("point", "ci"))
  driver_name <- if (analysis == "point") {
    "Point_estimation_right_censoring.R"
  } else {
    "Confidence_Intervals_right_censoring.R"
  }
  paths <- file.path(
    driver_dir,
    c(
      driver_name,
      "atomic_io_functions.R",
      "censored_simulation_helpers.R",
      "km_functions.R"
    )
  )
  if (any(!file.exists(paths))) {
    stop("Cannot fingerprint all censored-simulation sources.", call. = FALSE)
  }
  stats::setNames(unname(tools::md5sum(paths)), basename(paths))
}

## Inputs: analysis name, output directory, and task ID.
## Output: deterministic per-cell RDS checkpoint path.
cens_checkpoint_path <- function(analysis, output_dir, task_id) {
  analysis <- match.arg(analysis, c("point", "ci"))
  checkpoint_dir <- file.path(output_dir, "intermediate_results", analysis)
  if (!dir.exists(checkpoint_dir) &&
      !dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create checkpoint directory: ", checkpoint_dir,
         call. = FALSE)
  }
  file.path(checkpoint_dir, sprintf("%s_cell_%03d.rds", analysis, task_id))
}

## Input: one task. Output: stable task fields used to validate a checkpoint.
cens_task_signature <- function(task) {
  list(
    cell_id = as.integer(task$cell_id),
    family = task$family,
    cfg_name = task$cfg_name,
    cfg = as.numeric(task$cfg),
    cens = as.numeric(task$cens),
    n1 = as.integer(task$n1),
    n2 = as.integer(task$n2)
  )
}

## Inputs: analysis/task/run settings and source hashes.
## Output: metadata that must match before a cell checkpoint can be resumed.
cens_checkpoint_metadata <- function(analysis, task, mode, B, base_seed,
                                     source_md5, R_boot = NA_integer_,
                                     alpha = NA_real_) {
  list(
    schema = "extropy-censored-cell-v1",
    analysis = analysis,
    task = cens_task_signature(task),
    mode = mode,
    B = as.integer(B),
    R_boot = as.integer(R_boot),
    base_seed = as.integer(base_seed),
    alpha = as.numeric(alpha),
    RNG_kind = unname(cens_rng_signature()),
    R_version = R.version.string,
    platform = R.version$platform,
    source_md5 = source_md5
  )
}

## Inputs: checkpoint path, task, analysis, and optional expected metadata.
## Output: validated checkpoint object; stops on corruption or incompatibility.
cens_read_checkpoint <- function(path, task, analysis, expected = NULL) {
  object <- tryCatch(readRDS(path), error = function(e) e)
  required_metadata <- c(
    "schema", "analysis", "task", "mode", "B", "R_boot", "base_seed",
    "alpha", "RNG_kind", "R_version", "platform", "source_md5"
  )
  if (inherits(object, "error") || !is.list(object) ||
      !identical(object$schema, "extropy-censored-cell-v1") ||
      !identical(object$analysis, analysis) ||
      !identical(object$task, cens_task_signature(task)) ||
      !is.list(object$metadata) ||
      !all(required_metadata %in% names(object$metadata)) ||
      !is.list(object$result) ||
      !is.data.frame(object$result$summary) ||
      nrow(object$result$summary) != 1L ||
      !is.data.frame(object$result$replications)) {
    stop("Invalid or corrupted cell checkpoint: ", path, call. = FALSE)
  }
  stable_fields <- setdiff(names(expected), c("R_version", "platform"))
  if (!is.null(expected) &&
      !identical(object$metadata[stable_fields], expected[stable_fields])) {
    stop(
      "Incompatible checkpoint: ", path,
      ". Use a new output directory or --no-resume.",
      call. = FALSE
    )
  }
  B <- object$metadata$B
  reps <- object$result$replications
  expected_seeds <- vapply(
    seq_len(B),
    function(index) {
      cens_safe_seed(object$metadata$base_seed, task$cell_id, index)
    },
    integer(1)
  )
  if (nrow(reps) != B || !"replication" %in% names(reps) ||
      !"seed" %in% names(reps) ||
      !identical(as.integer(reps$replication), seq_len(B)) ||
      !identical(as.integer(reps$seed), expected_seeds)) {
    stop("Checkpoint replication contract failed: ", path, call. = FALSE)
  }
  object
}

## Inputs: analysis/task/settings, checkpoint policy, and source hashes.
## Output: computed or exactly compatible saved cell result.
cens_run_checkpointed_cell <- function(analysis, task, mode, B, base_seed,
                                       output_dir, resume, source_md5,
                                       R_boot = NA_integer_,
                                       alpha = NA_real_, z = NA_real_) {
  path <- cens_checkpoint_path(analysis, output_dir, task$cell_id)
  metadata <- cens_checkpoint_metadata(
    analysis, task, mode, B, base_seed, source_md5, R_boot, alpha
  )
  if (resume && file.exists(path)) {
    saved <- cens_read_checkpoint(path, task, analysis, metadata)
    message("Reused compatible checkpoint: ", basename(path))
    return(saved$result)
  }

  result <- if (analysis == "point") {
    cens_run_point_cell(task, B, base_seed)
  } else {
    cens_run_ci_cell(
      task, B, R_boot, base_seed, alpha, z, report = FALSE
    )
  }
  checkpoint <- list(
    schema = "extropy-censored-cell-v1",
    analysis = analysis,
    task = cens_task_signature(task),
    metadata = metadata,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    result = result
  )
  atomic_save_rds(checkpoint, path)
  verified <- cens_read_checkpoint(path, task, analysis, metadata)
  verified$result
}

## Inputs: analysis, tasks, run settings, selection, and worker count.
## Output: selected cell results after each has been atomically checkpointed.
cens_execute_tasks <- function(analysis, tasks, mode, B, base_seed, output_dir,
                               resume, source_md5, R_boot = NA_integer_,
                               alpha = NA_real_, z = NA_real_, cores = 1L) {
  analysis <- match.arg(analysis, c("point", "ci"))
  worker <- function(task) {
    ans <- cens_run_checkpointed_cell(
      analysis, task, mode, B, base_seed, output_dir, resume, source_md5,
      R_boot, alpha, z
    )
    if (analysis == "point") {
      cens_report_point_cell(ans)
    } else {
      cens_report_ci_cell(ans)
    }
    ans
  }
  if (analysis == "ci" && cores > 1L && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(
      tasks,
      worker,
      mc.cores = cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    if (analysis == "ci" && cores > 1L) {
      warning("--cores>1 is unavailable on Windows; using one core.",
              call. = FALSE)
    }
    results <- lapply(tasks, worker)
  }
  if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
    stop(
      "At least one worker failed; completed cells remain safely checkpointed.",
      call. = FALSE
    )
  }
  results
}

## Input: scenario, group sizes and variance flag. Output: observed pairs and KM fit.
cens_draw_observed <- function(sc, n1, n2, want_var) {
  generated <- sc$gen(n1, n2)
  t1 <- pmin(generated$X, generated$C)
  s1 <- as.integer(generated$X <= generated$C)
  t2 <- pmin(generated$Y, generated$E)
  s2 <- as.integer(generated$Y <= generated$E)
  list(
    ok = TRUE,
    r = D_km(t1, s1, t2, s2, sc$tau, want_var = want_var),
    t1 = t1,
    s1 = s1,
    t2 = t2,
    s2 = s2,
    censor1 = mean(s1 == 0L),
    censor2 = mean(s2 == 0L)
  )
}

## Input: list of cell results. Output: row-bound summary and replication frames.
cens_bind_results <- function(results) {
  summary <- do.call(rbind, lapply(results, `[[`, "summary"))
  replications <- do.call(rbind, lapply(results, `[[`, "replications"))
  rownames(summary) <- NULL
  rownames(replications) <- NULL
  list(summary = summary, replications = replications)
}

## Input: long results frame. Output: publication-table keys and initialized frame.
cens_wide_base <- function(df) {
  keys <- unique(df[, c("config", "setting", "censoring")])
  wide <- data.frame(
    Setting = keys$setting,
    Censoring = sprintf("%.0f%%", 100 * keys$censoring),
    stringsAsFactors = FALSE
  )
  wide$D_tau <- round(vapply(seq_len(nrow(keys)), function(i) {
    keep <- df$config == keys$config[i] & df$censoring == keys$censoring[i]
    df$D_tau[keep][1]
  }, numeric(1)), 6)
  list(keys = keys, wide = wide)
}

## Input: long frame, key frame, row/sample identifiers and column. Output: scalar value.
cens_pick_cell <- function(df, keys, i, nn, column) {
  keep <- df$config == keys$config[i] &
    df$censoring == keys$censoring[i] &
    df$n1 == nn[1] &
    df$n2 == nn[2]
  x <- df[[column]][keep]
  if (length(x) == 1L) x else NA_real_
}

## Input: one task plus Monte Carlo size/base seed. Output: summary and replication frames.
cens_run_point_cell <- function(task, B, base_seed) {
  sc <- make_scenario(task$family, task$cfg, task$cens)
  est <- rep(NA_real_, B)
  censor1 <- censor2 <- rep(NA_real_, B)
  y1 <- y2 <- rep(NA_real_, B)
  zero1 <- zero2 <- terminal1 <- terminal2 <- rep(NA, B)
  support_regular <- failed <- rep(FALSE, B)
  reason <- error_message <- rep("", B)
  rep_seed <- integer(B)

  for (b in seq_len(B)) {
    rep_seed[b] <- cens_safe_seed(base_seed, task$cell_id, b)
    set.seed(rep_seed[b])
    ans <- tryCatch(
      cens_draw_observed(sc, task$n1, task$n2, want_var = FALSE),
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )

    if (!isTRUE(ans$ok)) {
      failed[b] <- TRUE
      error_message[b] <- ans$error
    } else {
      r <- ans$r
      est[b] <- r$D
      censor1[b] <- ans$censor1
      censor2[b] <- ans$censor2
      y1[b] <- r$Y1tau
      y2[b] <- r$Y2tau
      zero1[b] <- r$zero_risk_1
      zero2[b] <- r$zero_risk_2
      terminal1[b] <- r$terminal_before_tau_1
      terminal2[b] <- r$terminal_before_tau_2
      support_regular[b] <- r$support_regular
      reason[b] <- r$nonregular_reason
    }
  }

  reps <- data.frame(
    family = task$family, config = task$cfg_name, setting = sc$label,
    target_censoring = task$cens, n1 = task$n1, n2 = task$n2,
    tau = sc$tau, D_tau = sc$D, replication = seq_len(B), seed = rep_seed,
    D_hat = est,
    realized_censoring_1 = censor1, realized_censoring_2 = censor2,
    at_risk_tau_1 = y1, at_risk_tau_2 = y2,
    zero_risk_tau_1 = zero1, zero_risk_tau_2 = zero2,
    terminal_before_tau_1 = terminal1, terminal_before_tau_2 = terminal2,
    support_regular = support_regular, failed = failed,
    nonregular_reason = reason, error_message = error_message,
    stringsAsFactors = FALSE
  )

  successful <- !failed & is.finite(est)
  regular_success <- successful & support_regular
  all_acc <- cens_accuracy(est, sc$D, successful, include_relmse = TRUE)
  regular_acc <- cens_accuracy(
    est,
    sc$D,
    regular_success,
    include_relmse = TRUE
  )

  summary <- data.frame(
    family = task$family, config = task$cfg_name, setting = sc$label,
    censoring = task$cens,
    theoretical_censoring_1 = sc$cens1,
    theoretical_censoring_2 = sc$cens2,
    realized_censoring_mean_1 = cens_safe_mean(censor1),
    realized_censoring_mean_2 = cens_safe_mean(censor2),
    realized_censoring_mcse_1 = cens_safe_sd(censor1) /
      sqrt(sum(is.finite(censor1))),
    realized_censoring_mcse_2 = cens_safe_sd(censor2) /
      sqrt(sum(is.finite(censor2))),
    n1 = task$n1, n2 = task$n2, B_requested = B,
    n_success = unname(all_acc["n"]), n_regular = unname(regular_acc["n"]),
    failure_rate = mean(failed),
    support_nonregular_rate = cens_safe_flag_rate(!support_regular, successful),
    support_nonregular_rate_all_attempts = mean(successful & !support_regular),
    regular_analysis_rate = mean(regular_success),
    tau = sc$tau, D_tau = sc$D,
    expected_at_risk_tau_1 = task$n1 * sc$P_obs_gt_tau_1,
    expected_at_risk_tau_2 = task$n2 * sc$P_obs_gt_tau_2,
    mean_at_risk_tau_1 = cens_safe_mean(y1),
    mean_at_risk_tau_2 = cens_safe_mean(y2),
    median_at_risk_tau_1 = cens_safe_quantile(y1, 0.5),
    median_at_risk_tau_2 = cens_safe_quantile(y2, 0.5),
    q05_at_risk_tau_1 = cens_safe_quantile(y1, 0.05),
    q05_at_risk_tau_2 = cens_safe_quantile(y2, 0.05),
    zero_risk_rate_1 = cens_safe_flag_rate(zero1, successful),
    zero_risk_rate_2 = cens_safe_flag_rate(zero2, successful),
    zero_risk_rate_1_all_attempts = mean(zero1 %in% TRUE),
    zero_risk_rate_2_all_attempts = mean(zero2 %in% TRUE),
    terminal_rate_1 = cens_safe_flag_rate(terminal1, successful),
    terminal_rate_2 = cens_safe_flag_rate(terminal2, successful),
    terminal_rate_1_all_attempts = mean(terminal1 %in% TRUE),
    terminal_rate_2_all_attempts = mean(terminal2 %in% TRUE),
    mean_est = unname(all_acc["mean"]), bias = unname(all_acc["bias"]),
    sd = unname(all_acc["sd"]), MCSE_mean = unname(all_acc["mcse_mean"]),
    MSE = unname(all_acc["mse"]), RMSE = unname(all_acc["rmse"]),
    RelMSE = unname(all_acc["relmse"]),
    mean_est_regular = unname(regular_acc["mean"]),
    bias_regular = unname(regular_acc["bias"]),
    sd_regular = unname(regular_acc["sd"]),
    MCSE_mean_regular = unname(regular_acc["mcse_mean"]),
    MSE_regular = unname(regular_acc["mse"]),
    RMSE_regular = unname(regular_acc["rmse"]),
    RelMSE_regular = unname(regular_acc["relmse"]),
    stringsAsFactors = FALSE
  )
  list(summary = summary, replications = reps)
}

## Input: one point-estimation cell result. Output: progress message, invisibly.
cens_report_point_cell <- function(ans) {
  s <- ans$summary
  status_fmt <- paste0(
    "%-7s %s %2.0f%% n=(%d,%d): RMSE=%.5g, regular=%.1f%%, ",
    "zero-risk=(%.1f%%, %.1f%%)"
  )
  message(sprintf(
    status_fmt,
    s$family,
    s$config,
    100 * s$censoring,
    s$n1,
    s$n2,
    s$RMSE,
    100 * s$regular_analysis_rate,
    100 * s$zero_risk_rate_1,
    100 * s$zero_risk_rate_2
  ))
  invisible(ans)
}

## Input: long results, value column and sample sizes. Output: publication table.
cens_make_point_wide <- function(df, value, sample_sizes) {
  initialized <- cens_wide_base(df)
  keys <- initialized$keys
  wide <- initialized$wide
  for (nn in sample_sizes) {
    tag <- sprintf("(%d,%d)", nn[1], nn[2])
    wide[[tag]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(cens_pick_cell(df, keys, i, nn, value), 6)
    }, numeric(1))
  }
  wide
}

## Input: CI task/run settings. Output: summary and replication frames.
cens_run_ci_cell <- function(task, B, R_boot, base_seed, alpha, z,
                             report = FALSE) {
  sc <- make_scenario(task$family, task$cfg, task$cens)
  Dhat <- se <- bias_hat <- rep(NA_real_, B)
  censor1 <- censor2 <- y1 <- y2 <- rep(NA_real_, B)
  zero1 <- zero2 <- terminal1 <- terminal2 <- rep(NA, B)
  support_regular <- inference_regular <- failed <- bootstrap_failed <- rep(FALSE, B)
  reason <- error_message <- bootstrap_error <- rep("", B)
  rep_seed <- integer(B)

  normal_lo <- normal_hi <- percentile_lo <- percentile_hi <- rep(NA_real_, B)
  basic_lo <- basic_hi <- biascorr_lo <- biascorr_hi <- rep(NA_real_, B)
  boot_sd <- boot_mean <- boot_zero1 <- boot_zero2 <- rep(NA_real_, B)
  boot_terminal1 <- boot_terminal2 <- boot_nonregular <- boot_regular_n <-
    rep(NA_real_, B)

  for (b in seq_len(B)) {
    rep_seed[b] <- cens_safe_seed(base_seed, task$cell_id, b)
    set.seed(rep_seed[b])
    ans <- tryCatch(
      cens_draw_observed(sc, task$n1, task$n2, want_var = TRUE),
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )

    if (!isTRUE(ans$ok)) {
      failed[b] <- TRUE
      error_message[b] <- ans$error
      next
    }

    r <- ans$r
    Dhat[b] <- r$D
    se[b] <- r$se
    bias_hat[b] <- r$bias_hat
    censor1[b] <- ans$censor1
    censor2[b] <- ans$censor2
    y1[b] <- r$Y1tau
    y2[b] <- r$Y2tau
    zero1[b] <- r$zero_risk_1
    zero2[b] <- r$zero_risk_2
    terminal1[b] <- r$terminal_before_tau_1
    terminal2[b] <- r$terminal_before_tau_2
    support_regular[b] <- r$support_regular
    inference_regular[b] <- r$inference_regular
    reason[b] <- r$nonregular_reason

    if (!r$inference_regular) next

    normal_lo[b] <- r$D - z * r$se
    normal_hi[b] <- r$D + z * r$se
    ## Exploratory bias-centred interval; it is deliberately not truncated at 0.
    biascorr_lo[b] <- r$D - r$bias_hat - z * r$se
    biascorr_hi[b] <- r$D - r$bias_hat + z * r$se

    bt <- tryCatch(
      boot_engine(
        ans$t1,
        ans$s1,
        ans$t2,
        ans$s2,
        sc$tau,
        r$D,
        R = R_boot,
        alpha = alpha
      ),
      error = function(e) e
    )
    if (inherits(bt, "error")) {
      bootstrap_failed[b] <- TRUE
      bootstrap_error[b] <- conditionMessage(bt)
      next
    }
    percentile_lo[b] <- unname(bt["perc_lo"])
    percentile_hi[b] <- unname(bt["perc_hi"])
    basic_lo[b] <- unname(bt["basic_lo"])
    basic_hi[b] <- unname(bt["basic_hi"])
    boot_sd[b] <- unname(bt["boot_sd"])
    boot_mean[b] <- unname(bt["boot_mean"])
    boot_zero1[b] <- unname(bt["boot_zero_risk_rate_1"])
    boot_zero2[b] <- unname(bt["boot_zero_risk_rate_2"])
    boot_terminal1[b] <- unname(bt["boot_terminal_rate_1"])
    boot_terminal2[b] <- unname(bt["boot_terminal_rate_2"])
    boot_nonregular[b] <- unname(bt["boot_nonregular_rate"])
    boot_regular_n[b] <- unname(bt["boot_regular_n"])
  }

  successful <- !failed & is.finite(Dhat)
  regular_success <- successful & inference_regular
  bootstrap_attempted <- regular_success & !bootstrap_failed
  acc_all <- cens_accuracy(Dhat, sc$D, successful)
  acc_regular <- cens_accuracy(Dhat, sc$D, regular_success)
  normal <- cens_interval_metrics(normal_lo, normal_hi, sc$D, regular_success)
  percentile <- cens_interval_metrics(
    percentile_lo,
    percentile_hi,
    sc$D,
    bootstrap_attempted
  )
  basic <- cens_interval_metrics(basic_lo, basic_hi, sc$D, bootstrap_attempted)
  biascorr <- cens_interval_metrics(
    biascorr_lo,
    biascorr_hi,
    sc$D,
    regular_success
  )

  reps <- data.frame(
    family = task$family, config = task$cfg_name, setting = sc$label,
    target_censoring = task$cens, n1 = task$n1, n2 = task$n2,
    tau = sc$tau, D_tau = sc$D, replication = seq_len(B), seed = rep_seed,
    D_hat = Dhat, SE_Greenwood = se, bias_hat = bias_hat,
    realized_censoring_1 = censor1, realized_censoring_2 = censor2,
    at_risk_tau_1 = y1, at_risk_tau_2 = y2,
    zero_risk_tau_1 = zero1, zero_risk_tau_2 = zero2,
    terminal_before_tau_1 = terminal1, terminal_before_tau_2 = terminal2,
    support_regular = support_regular, inference_regular = inference_regular,
    Normal_lo = normal_lo, Normal_hi = normal_hi,
    Percentile_lo = percentile_lo, Percentile_hi = percentile_hi,
    Basic_lo = basic_lo, Basic_hi = basic_hi,
    BiasCorr_lo = biascorr_lo, BiasCorr_hi = biascorr_hi,
    bootstrap_SD = boot_sd, bootstrap_mean = boot_mean,
    bootstrap_zero_risk_rate_1 = boot_zero1,
    bootstrap_zero_risk_rate_2 = boot_zero2,
    bootstrap_terminal_rate_1 = boot_terminal1,
    bootstrap_terminal_rate_2 = boot_terminal2,
    bootstrap_nonregular_rate = boot_nonregular,
    bootstrap_regular_n = boot_regular_n,
    failed = failed, bootstrap_failed = bootstrap_failed,
    nonregular_reason = reason, error_message = error_message,
    bootstrap_error = bootstrap_error,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    family = task$family, config = task$cfg_name, setting = sc$label,
    censoring = task$cens,
    theoretical_censoring_1 = sc$cens1,
    theoretical_censoring_2 = sc$cens2,
    realized_censoring_mean_1 = cens_safe_mean(censor1),
    realized_censoring_mean_2 = cens_safe_mean(censor2),
    realized_censoring_mcse_1 = cens_safe_sd(censor1) /
      sqrt(sum(is.finite(censor1))),
    realized_censoring_mcse_2 = cens_safe_sd(censor2) /
      sqrt(sum(is.finite(censor2))),
    n1 = task$n1, n2 = task$n2, B_requested = B, R_boot = R_boot,
    n_success = unname(acc_all["n"]),
    n_inference_regular = unname(acc_regular["n"]),
    failure_rate = mean(failed),
    support_nonregular_rate = cens_safe_flag_rate(!support_regular, successful),
    support_nonregular_rate_all_attempts = mean(successful & !support_regular),
    variance_nonpositive_rate = cens_safe_flag_rate(
      support_regular & !inference_regular,
      successful
    ),
    variance_nonpositive_rate_all_attempts = mean(
      successful & support_regular & !inference_regular
    ),
    inference_regular_rate = mean(regular_success),
    bootstrap_failure_rate = cens_safe_flag_rate(
      bootstrap_failed,
      regular_success
    ),
    bootstrap_failure_rate_all_attempts = mean(regular_success & bootstrap_failed),
    tau = sc$tau, D_tau = sc$D,
    expected_at_risk_tau_1 = task$n1 * sc$P_obs_gt_tau_1,
    expected_at_risk_tau_2 = task$n2 * sc$P_obs_gt_tau_2,
    mean_at_risk_tau_1 = cens_safe_mean(y1),
    mean_at_risk_tau_2 = cens_safe_mean(y2),
    median_at_risk_tau_1 = cens_safe_quantile(y1, 0.5),
    median_at_risk_tau_2 = cens_safe_quantile(y2, 0.5),
    q05_at_risk_tau_1 = cens_safe_quantile(y1, 0.05),
    q05_at_risk_tau_2 = cens_safe_quantile(y2, 0.05),
    zero_risk_rate_1 = cens_safe_flag_rate(zero1, successful),
    zero_risk_rate_2 = cens_safe_flag_rate(zero2, successful),
    zero_risk_rate_1_all_attempts = mean(zero1 %in% TRUE),
    zero_risk_rate_2_all_attempts = mean(zero2 %in% TRUE),
    terminal_rate_1 = cens_safe_flag_rate(terminal1, successful),
    terminal_rate_2 = cens_safe_flag_rate(terminal2, successful),
    terminal_rate_1_all_attempts = mean(terminal1 %in% TRUE),
    terminal_rate_2_all_attempts = mean(terminal2 %in% TRUE),
    mean_est = unname(acc_all["mean"]), bias = unname(acc_all["bias"]),
    empirical_SD = unname(acc_all["sd"]),
    MCSE_mean = unname(acc_all["mcse_mean"]),
    MSE = unname(acc_all["mse"]), RMSE = unname(acc_all["rmse"]),
    mean_est_regular = unname(acc_regular["mean"]),
    bias_regular = unname(acc_regular["bias"]),
    empirical_SD_regular = unname(acc_regular["sd"]),
    MCSE_mean_regular = unname(acc_regular["mcse_mean"]),
    RMSE_regular = unname(acc_regular["rmse"]),
    mean_SE_Greenwood = cens_safe_mean(se[regular_success]),
    SD_to_meanSE_ratio = unname(acc_regular["sd"]) /
      cens_safe_mean(se[regular_success]),
    mean_bootstrap_SD = cens_safe_mean(boot_sd),
    mean_bias_hat = cens_safe_mean(bias_hat[regular_success]),
    mean_bootstrap_zero_risk_rate_1 = cens_safe_mean(boot_zero1),
    mean_bootstrap_zero_risk_rate_2 = cens_safe_mean(boot_zero2),
    mean_bootstrap_terminal_rate_1 = cens_safe_mean(boot_terminal1),
    mean_bootstrap_terminal_rate_2 = cens_safe_mean(boot_terminal2),
    mean_bootstrap_nonregular_rate = cens_safe_mean(boot_nonregular),
    Normal_n_valid = unname(normal["n_valid"]),
    Normal_availability = unname(normal["availability"]),
    Normal_CP = unname(normal["cp"]),
    Normal_CP_MCSE = unname(normal["cp_mcse"]),
    Normal_CP_all_attempts = unname(normal["cp_all_attempts"]),
    Normal_CP_all_attempts_MCSE = unname(normal["cp_all_attempts_mcse"]),
    Normal_AL = unname(normal["average_length"]),
    Bootstrap_n_valid = unname(percentile["n_valid"]),
    Bootstrap_availability = unname(percentile["availability"]),
    Bootstrap_CP = unname(percentile["cp"]),
    Bootstrap_CP_MCSE = unname(percentile["cp_mcse"]),
    Bootstrap_CP_all_attempts = unname(percentile["cp_all_attempts"]),
    Bootstrap_CP_all_attempts_MCSE = unname(percentile["cp_all_attempts_mcse"]),
    Bootstrap_AL = unname(percentile["average_length"]),
    Basic_CP = unname(basic["cp"]),
    Basic_CP_MCSE = unname(basic["cp_mcse"]),
    Basic_AL = unname(basic["average_length"]),
    BiasCorr_CP = unname(biascorr["cp"]),
    BiasCorr_CP_MCSE = unname(biascorr["cp_mcse"]),
    BiasCorr_AL = unname(biascorr["average_length"]),
    stringsAsFactors = FALSE
  )
  ans <- list(summary = summary, replications = reps)
  if (report) cens_report_ci_cell(ans)
  ans
}

## Input: one CI cell result. Output: progress message, invisibly.
cens_report_ci_cell <- function(ans) {
  s <- ans$summary
  status_fmt <- paste0(
    "%-7s %s %2.0f%% n=(%d,%d): all-attempt Normal %.1f%% (MCSE %.2f), ",
    "all-attempt Bootstrap %.1f%% (MCSE %.2f), regular %.1f%%"
  )
  message(sprintf(
    status_fmt,
    s$family,
    s$config,
    100 * s$censoring,
    s$n1,
    s$n2,
    s$Normal_CP_all_attempts,
    s$Normal_CP_all_attempts_MCSE,
    s$Bootstrap_CP_all_attempts,
    s$Bootstrap_CP_all_attempts_MCSE,
    100 * s$inference_regular_rate
  ))
  invisible(ans)
}

## Input: long CI results and sample sizes. Output: publication-layout CI table.
cens_make_ci_wide <- function(df, sample_sizes) {
  initialized <- cens_wide_base(df)
  keys <- initialized$keys
  wide <- initialized$wide
  pick <- function(i, nn, column) {
    cens_pick_cell(df, keys, i, nn, column)
  }
  for (nn in sample_sizes) {
    tag <- sprintf("(%d,%d)", nn[1], nn[2])
    wide[[paste0("Normal_CP ", tag)]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(pick(i, nn, "Normal_CP_all_attempts"), 1)
    }, numeric(1))
    wide[[paste0("Normal_MCSE ", tag)]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(pick(i, nn, "Normal_CP_all_attempts_MCSE"), 2)
    }, numeric(1))
    wide[[paste0("Normal_AL ", tag)]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(pick(i, nn, "Normal_AL"), 6)
    }, numeric(1))
    wide[[paste0("Bootstrap_CP ", tag)]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(pick(i, nn, "Bootstrap_CP_all_attempts"), 1)
    }, numeric(1))
    wide[[paste0("Bootstrap_MCSE ", tag)]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(pick(i, nn, "Bootstrap_CP_all_attempts_MCSE"), 2)
    }, numeric(1))
    wide[[paste0("Bootstrap_AL ", tag)]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(pick(i, nn, "Bootstrap_AL"), 6)
    }, numeric(1))
    wide[[paste0("Valid_% ", tag)]] <- vapply(seq_len(nrow(keys)), function(i) {
      round(100 * pick(i, nn, "inference_regular_rate"), 1)
    }, numeric(1))
  }
  wide
}

## Inputs: analysis name, task list, run settings, output path, and source hashes.
## Output: all validated checkpoint results and a checksum/provenance manifest.
cens_collect_checkpoints <- function(analysis, tasks, mode, B, base_seed,
                                     output_dir, source_md5,
                                     R_boot = NA_integer_,
                                     alpha = NA_real_) {
  analysis <- match.arg(analysis, c("point", "ci"))
  results <- vector("list", length(tasks))
  objects <- vector("list", length(tasks))
  paths <- character(length(tasks))
  problems <- character()

  for (i in seq_along(tasks)) {
    task <- tasks[[i]]
    path <- cens_checkpoint_path(analysis, output_dir, task$cell_id)
    paths[[i]] <- path
    if (!file.exists(path)) {
      problems <- c(problems, sprintf("task %d: missing", task$cell_id))
      next
    }
    expected <- cens_checkpoint_metadata(
      analysis, task, mode, B, base_seed, source_md5, R_boot, alpha
    )
    object <- tryCatch(
      cens_read_checkpoint(path, task, analysis, expected),
      error = function(e) e
    )
    if (inherits(object, "error")) {
      problems <- c(
        problems,
        sprintf("task %d: %s", task$cell_id, conditionMessage(object))
      )
      next
    }
    objects[[i]] <- object
    results[[i]] <- object$result
  }

  if (length(problems)) {
    stop(
      "Cannot aggregate censored ", analysis, " checkpoints:\n  - ",
      paste(problems, collapse = "\n  - "),
      call. = FALSE
    )
  }

  provenance <- lapply(objects, function(object) {
    object$metadata[c("R_version", "platform")]
  })
  if (!all(vapply(provenance, identical, logical(1), provenance[[1L]]))) {
    stop("Checkpoints were computed under mixed R/platform provenance.",
         call. = FALSE)
  }

  task_frame <- cens_task_manifest(tasks)
  info <- file.info(paths)
  checkpoint_manifest <- cbind(
    task_frame,
    data.frame(
      checkpoint = file.path(
        "intermediate_results",
        analysis,
        basename(paths)
      ),
      bytes = unname(info$size),
      md5 = unname(tools::md5sum(paths)),
      created_utc = vapply(objects, `[[`, character(1), "created_utc"),
      R_version = vapply(
        objects,
        function(object) object$metadata$R_version,
        character(1)
      ),
      platform = vapply(
        objects,
        function(object) object$metadata$platform,
        character(1)
      ),
      stringsAsFactors = FALSE
    )
  )
  list(results = results, manifest = checkpoint_manifest)
}

## Inputs: point-cell results, sample sizes, output path, and manifests.
## Output: the canonical Tables 7--8 files and diagnostics, written atomically.
cens_write_point_outputs <- function(results, sample_sizes, output_dir,
                                     task_manifest, checkpoint_manifest) {
  bound <- cens_bind_results(results)
  res <- bound$summary
  reps <- bound$replications
  atomic_write_csv(
    res,
    file.path(output_dir, "Section54_point_estimation_long.csv")
  )
  atomic_write_csv(
    reps,
    file.path(output_dir, "Section54_point_estimation_replications.csv")
  )
  atomic_write_csv(
    cens_make_point_wide(res[res$family == "exp", ], "RelMSE", sample_sizes),
    file.path(output_dir, "Table7_RelMSE_Exponential.csv")
  )
  atomic_write_csv(
    cens_make_point_wide(
      res[res$family == "weibull", ], "RelMSE", sample_sizes
    ),
    file.path(output_dir, "Table8_RelMSE_Weibull.csv")
  )
  atomic_write_csv(
    cens_make_point_wide(res[res$family == "exp", ], "RMSE", sample_sizes),
    file.path(output_dir, "Table7_RMSE_Exponential.csv")
  )
  atomic_write_csv(
    cens_make_point_wide(
      res[res$family == "weibull", ], "RMSE", sample_sizes
    ),
    file.path(output_dir, "Table8_RMSE_Weibull.csv")
  )
  atomic_write_csv(
    task_manifest,
    file.path(output_dir, "Section54_point_task_manifest.csv")
  )
  atomic_write_csv(
    checkpoint_manifest,
    file.path(output_dir, "Section54_point_checkpoint_manifest.csv")
  )
  if (any(res$failure_rate > 0)) {
    warning(
      "Some Monte Carlo replications failed; inspect the replication output.",
      call. = FALSE
    )
  }
  invisible(bound)
}

## Inputs: CI-cell results, sample sizes, output path, and manifests.
## Output: the canonical Tables 9--10 files and diagnostics, written atomically.
cens_write_ci_outputs <- function(results, sample_sizes, output_dir,
                                  task_manifest, checkpoint_manifest) {
  bound <- cens_bind_results(results)
  res <- bound$summary
  reps <- bound$replications
  atomic_write_csv(res, file.path(output_dir, "Section54_CI_long.csv"))
  atomic_write_csv(
    reps,
    file.path(output_dir, "Section54_CI_replications.csv")
  )
  atomic_write_csv(
    cens_make_ci_wide(res[res$family == "exp", ], sample_sizes),
    file.path(output_dir, "Table9_CI_Exponential.csv")
  )
  atomic_write_csv(
    cens_make_ci_wide(res[res$family == "weibull", ], sample_sizes),
    file.path(output_dir, "Table10_CI_Weibull.csv")
  )
  atomic_write_csv(
    task_manifest,
    file.path(output_dir, "Section54_CI_task_manifest.csv")
  )
  atomic_write_csv(
    checkpoint_manifest,
    file.path(output_dir, "Section54_CI_checkpoint_manifest.csv")
  )
  if (any(res$failure_rate > 0 | res$bootstrap_failure_rate > 0,
          na.rm = TRUE)) {
    warning(
      "Some replications failed; inspect Section54_CI_replications.csv.",
      call. = FALSE
    )
  }
  if (any(res$inference_regular_rate < 0.95, na.rm = TRUE)) {
    warning(
      paste(
        "At least one cell has fewer than 95% inference-regular samples;",
        "report this limitation."
      ),
      call. = FALSE
    )
  }
  invisible(bound)
}

## Inputs: analysis/task design, exact run settings, output path, and sources.
## Output: canonical long/wide results rebuilt solely from complete checkpoints.
cens_aggregate_checkpoints <- function(analysis, tasks, sample_sizes, mode, B,
                                       base_seed, output_dir, source_md5,
                                       R_boot = NA_integer_,
                                       alpha = NA_real_) {
  collected <- cens_collect_checkpoints(
    analysis, tasks, mode, B, base_seed, output_dir, source_md5,
    R_boot, alpha
  )
  task_manifest <- cens_task_manifest(tasks)
  if (analysis == "point") {
    cens_write_point_outputs(
      collected$results,
      sample_sizes,
      output_dir,
      task_manifest,
      collected$manifest
    )
  } else {
    cens_write_ci_outputs(
      collected$results,
      sample_sizes,
      output_dir,
      task_manifest,
      collected$manifest
    )
  }
}
