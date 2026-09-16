# =============================================================================
# Regression test: two same-seed quick runs must be byte-identical
# =============================================================================

test_location <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Run this test with Rscript.")
  }

  test_file <- sub("^--file=", "", file_arg[[1L]])
  decoded_file <- gsub("~+~", " ", test_file, fixed = TRUE)
  if (!file.exists(test_file) && file.exists(decoded_file)) {
    test_file <- decoded_file
  }

  dirname(normalizePath(
    test_file,
    winslash = "/",
    mustWork = TRUE
  ))
}


parse_tables <- function(args) {
  table_arg <- grep("^--tables=", args, value = TRUE)
  if (length(table_arg) == 0L) {
    return(1:6)
  }
  if (length(table_arg) > 1L) {
    stop("Use --tables only once.")
  }

  tables <- suppressWarnings(as.integer(strsplit(
    sub("^--tables=", "", table_arg),
    ",",
    fixed = TRUE
  )[[1L]]))

  if (length(tables) == 0L || anyNA(tables) || any(!tables %in% 1:6)) {
    stop("--tables must be a comma-separated subset of 1,2,3,4,5,6.")
  }

  unique(tables)
}


run_quick <- function(label, root, tables) {
  runner <- file.path(test_location(), "run_tables_1_6.R")
  output_root <- file.path(root, label)
  args <- c(
    "--vanilla",
    shQuote(runner),
    "--mode=quick",
    sprintf("--tables=%s", paste(tables, collapse = ",")),
    sprintf("--output-dir=%s", shQuote(output_root))
  )

  log <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = args,
    stdout = TRUE,
    stderr = TRUE,
    wait = TRUE
  )
  status <- attr(log, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (status != 0L) {
    stop(paste(c(sprintf("Quick run %s failed:", label), log), collapse = "\n"))
  }

  file.path(output_root, "quick")
}


main <- function() {
  tables <- parse_tables(commandArgs(trailingOnly = TRUE))
  validator <- new.env(parent = baseenv())
  sys.source(file.path(test_location(), "validate_tables_1_6.R"), envir = validator)
  plan <- validator$tables_1_6_plan(test_location())
  plan <- plan[match(tables, plan$Table), , drop = FALSE]

  root <- tempfile("tables_1_6_determinism_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  first <- run_quick("first", root, tables)
  second <- run_quick("second", root, tables)

  relative_files <- unlist(lapply(seq_len(nrow(plan)), function(i) {
    table_dir <- sprintf("table_%d", plan$Table[[i]])
    c(
      file.path(table_dir, plan$Output[[i]]),
      file.path(table_dir, sprintf("Table_%d_manuscript.csv", plan$Table[[i]]))
    )
  }))

  first_files <- file.path(first, relative_files)
  second_files <- file.path(second, relative_files)
  first_hash <- unname(tools::md5sum(first_files))
  second_hash <- unname(tools::md5sum(second_files))
  different <- relative_files[is.na(first_hash) | is.na(second_hash) | first_hash != second_hash]

  if (length(different) > 0L) {
    stop(paste(
      c("Same-seed quick runs were not identical:", paste0("- ", different)),
      collapse = "\n"
    ))
  }

  cat(sprintf(
    "Determinism test passed for Table(s) %s (%d files compared).\n",
    paste(tables, collapse = ", "),
    length(relative_files)
  ))
  invisible(TRUE)
}


if (sys.nframe() == 0L) {
  main()
}
