## =============================================================================
## Shared Section 5.3 simulation engine.
##
## Every Monte Carlo replication draws its samples after set.seed() with a seed
## derived from (base seed, cell id, replication index).  The result of a cell
## therefore does not depend on how many cells were run before it, nor on
## whether cells were run in parallel, which makes the tables reproducible cell
## by cell and safe to distribute over cores.
##
## Requires: base R; parallel (shipped with R) when cores > 1.
## =============================================================================

## Inputs: base seed, cell index, replication index.
## Output: a distinct integer seed; the spacing of 1,000,003 keeps the streams
## of different cells from overlapping for any realistic number of replications.
unc_safe_seed <- function(base_seed, cell_id, replication) {
  as.integer(
    (as.double(base_seed) + 1000003 * cell_id + replication) %% 2147483646 + 1
  )
}

## Input: none. Output: the RNG algorithms every Section 5.3 run must use.
unc_rng_signature <- function() {
  c(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}

## Inputs: parameter list, sample-size pairs, a sampler and a true-value
## function, plus Monte Carlo settings.
## Output: one row per (parameter setting, sample-size pair).
unc_run_coverage <- function(params_list, sizes, draw, true_value,
                             iterations = 2000L, B_boot = 499L,
                             grid_size = 400L, tail_mult = 8,
                             alpha = 0.05, seed = 2026L, cores = 1L,
                             verbose = TRUE, checkpoint_dir = NULL,
                             cells_wanted = NULL, resume = TRUE,
                             with_kernel = FALSE) {
  stopifnot(is.list(params_list), is.list(sizes), is.function(draw),
            is.function(true_value))
  iterations <- as.integer(iterations)
  cores <- max(1L, as.integer(cores))
  do.call(RNGkind, as.list(unc_rng_signature()))

  cells <- list()
  for (p in seq_along(params_list)) {
    for (k in seq_along(sizes)) {
      cells[[length(cells) + 1L]] <- list(
        cell_id = length(cells) + 1L,
        params = params_list[[p]],
        n1 = as.integer(sizes[[k]][1]),
        n2 = as.integer(sizes[[k]][2])
      )
    }
  }

  methods <- c("JEL", "NA_", "Emp", "Ker")

  run_cell <- function(cell) {
    theta <- true_value(cell$params)
    covered <- valid <- stats::setNames(rep(0L, 4L), methods)
    lengths <- lapply(methods, function(m) rep(NA_real_, iterations))
    names(lengths) <- methods

    for (h in seq_len(iterations)) {
      set.seed(unc_safe_seed(seed, cell$cell_id, h))
      d <- draw(cell$params, cell$n1, cell$n2)
      res <- tryCatch(
        all_CIs(X = d$X, Y = d$Y, theta = theta, B_boot = B_boot,
                alpha = alpha, grid_size = grid_size, tail_mult = tail_mult,
                with_kernel = with_kernel),
        error = function(e) list(
          cov = stats::setNames(rep(NA, 4L), methods),
          al = stats::setNames(rep(NA_real_, 4L), methods)
        )
      )
      for (m in methods) {
        if (!is.na(res$cov[m])) {
          covered[m] <- covered[m] + as.integer(isTRUE(as.logical(res$cov[m])))
          valid[m] <- valid[m] + 1L
        }
        lengths[[m]][h] <- unname(res$al[m])
      }
    }

    cp <- vapply(methods, function(m) safe_cp(covered[m], valid[m]), numeric(1))
    al <- vapply(methods, function(m) safe_mean(lengths[[m]]), numeric(1))

    data.frame(
      cell_id = cell$cell_id,
      t(stats::setNames(cell$params, paste0("par", seq_along(cell$params)))),
      true_D = theta, n1 = cell$n1, n2 = cell$n2, iterations = iterations,
      JEL_CP = unname(cp["JEL"]), JEL_AL = unname(al["JEL"]),
      NA_CP = unname(cp["NA_"]), NA_AL = unname(al["NA_"]),
      Emp_CP = unname(cp["Emp"]), Emp_AL = unname(al["Emp"]),
      Ker_CP = unname(cp["Ker"]), Ker_AL = unname(al["Ker"]),
      Valid_JEL = unname(valid["JEL"]), Valid_NA = unname(valid["NA_"]),
      Valid_Emp = unname(valid["Emp"]), Valid_Ker = unname(valid["Ker"]),
      row.names = NULL
    )
  }

  ## A checkpoint is reused only when everything that could change the numbers
  ## matches: the cell itself, the replication count, the bootstrap size, the
  ## quadrature settings, alpha, the base seed and the RNG kinds.
  signature <- function(cell) {
    list(cell_id = cell$cell_id, params = cell$params, n1 = cell$n1, n2 = cell$n2,
         iterations = iterations, B_boot = B_boot, grid_size = grid_size,
         tail_mult = tail_mult, alpha = alpha, seed = seed,
         with_kernel = isTRUE(with_kernel),
         rng = unname(unc_rng_signature()))
  }
  cell_path <- function(cell) {
    file.path(checkpoint_dir, sprintf("cell_%03d.rds", cell$cell_id))
  }
  if (!is.null(checkpoint_dir) && !dir.exists(checkpoint_dir) &&
      !dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create checkpoint directory: ", checkpoint_dir, call. = FALSE)
  }

  load_cell <- function(cell) {
    if (is.null(checkpoint_dir)) return(NULL)
    p <- cell_path(cell)
    if (!file.exists(p)) return(NULL)
    obj <- tryCatch(readRDS(p), error = function(e) NULL)
    if (is.null(obj) || !identical(obj$signature, signature(cell))) return(NULL)
    obj$row
  }
  save_cell <- function(cell, row) {
    if (is.null(checkpoint_dir)) return(invisible(NULL))
    p <- cell_path(cell)
    tmp <- tempfile(pattern = basename(p), tmpdir = dirname(p), fileext = ".rds")
    saveRDS(list(signature = signature(cell), row = row), tmp, version = 3L)
    if (!isTRUE(suppressWarnings(file.rename(tmp, p)))) {
      unlink(tmp)
      stop("Could not write checkpoint: ", p, call. = FALSE)
    }
    invisible(p)
  }

  announce <- function(cell) {
    cached <- if (resume) load_cell(cell) else NULL
    if (!is.null(cached)) {
      if (verbose) message(sprintf("cell %2d/%2d: reused checkpoint",
                                   cell$cell_id, length(cells)))
      return(cached)
    }
    if (verbose) {
      message(sprintf("cell %2d/%2d: params (%s), n = (%d, %d)",
                      cell$cell_id, length(cells),
                      paste(format(cell$params, trim = TRUE), collapse = ", "),
                      cell$n1, cell$n2))
    }
    row <- run_cell(cell)
    save_cell(cell, row)
    row
  }

  selected <- if (is.null(cells_wanted)) cells else {
    bad <- setdiff(cells_wanted, seq_along(cells))
    if (length(bad)) {
      stop("--cells refers to cell(s) outside 1..", length(cells), ": ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
    cells[sort(unique(cells_wanted))]
  }

  rows <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(selected, announce, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(selected, announce)
  }

  failed <- !vapply(rows, is.data.frame, logical(1))
  if (any(failed)) {
    stop("Simulation cell(s) failed: ",
         paste(vapply(selected[failed], `[[`, numeric(1), "cell_id"),
               collapse = ", "), call. = FALSE)
  }

  ## When cells were run in pieces, assemble the table from every checkpoint the
  ## run directory holds, and refuse to write a partial table.
  if (!is.null(checkpoint_dir) && !is.null(cells_wanted)) {
    rows <- lapply(cells, load_cell)
    if (any(vapply(rows, is.null, logical(1)))) {
      missing <- which(vapply(rows, is.null, logical(1)))
      message("Cells still to run: ", paste(missing, collapse = ", "))
      return(NULL)
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$cell_id), ]
}

## Inputs: the command line, plus defaults. Output: a validated settings list.
unc_parse_cli <- function(args, defaults) {
  opt <- defaults
  for (a in args) {
    if (startsWith(a, "--iterations=")) opt$iterations <- as.integer(sub("^--iterations=", "", a))
    else if (startsWith(a, "--B-boot=")) opt$B_boot <- as.integer(sub("^--B-boot=", "", a))
    else if (startsWith(a, "--seed=")) opt$seed <- as.integer(sub("^--seed=", "", a))
    else if (startsWith(a, "--cores=")) opt$cores <- as.integer(sub("^--cores=", "", a))
    else if (startsWith(a, "--output-dir=")) opt$output_dir <- sub("^--output-dir=", "", a)
    else if (startsWith(a, "--cells=")) opt$cells <- unc_parse_cells(sub("^--cells=", "", a))
    else if (identical(a, "--no-resume")) opt$resume <- FALSE
    else if (identical(a, "--with-kernel")) opt$with_kernel <- TRUE
    else if (identical(a, "--quick")) { opt$iterations <- 50L; opt$B_boot <- 99L }
    else if (identical(a, "--quiet")) opt$verbose <- FALSE
    else stop("Unknown argument: ", a, call. = FALSE)
  }
  for (nm in c("iterations", "B_boot", "seed", "cores")) {
    if (!is.finite(opt[[nm]]) || opt[[nm]] < 0) {
      stop("--", gsub("_", "-", nm), " must be a nonnegative integer.", call. = FALSE)
    }
  }
  if (opt$iterations < 1L) stop("--iterations must be at least 1.", call. = FALSE)
  opt$cores <- max(1L, opt$cores)
  opt
}

## Inputs: a results frame and a destination path. Output: the path, invisibly.
unc_write_csv <- function(object, path) {
  dir <- dirname(path)
  if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory: ", dir, call. = FALSE)
  }
  tmp <- tempfile(pattern = basename(path), tmpdir = dir, fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(object, tmp, row.names = FALSE)
  if (!isTRUE(suppressWarnings(file.rename(tmp, path)))) {
    stop("Could not write: ", path, call. = FALSE)
  }
  invisible(path)
}

## Input: a cell specification such as "1,4,9-12". Output: an integer vector.
unc_parse_cells <- function(spec) {
  pieces <- strsplit(spec, ",", fixed = TRUE)[[1L]]
  pieces <- trimws(pieces[nzchar(trimws(pieces))])
  if (!length(pieces)) stop("--cells is empty.", call. = FALSE)
  ids <- unlist(lapply(pieces, function(p) {
    if (grepl("^[0-9]+$", p)) return(as.integer(p))
    m <- regmatches(p, regexec("^([0-9]+)-([0-9]+)$", p))[[1L]]
    if (length(m) != 3L) stop("Cannot parse --cells entry: ", p, call. = FALSE)
    a <- as.integer(m[2L]); b <- as.integer(m[3L])
    if (a > b) stop("Reversed range in --cells: ", p, call. = FALSE)
    seq.int(a, b)
  }), use.names = FALSE)
  sort(unique(as.integer(ids)))
}
