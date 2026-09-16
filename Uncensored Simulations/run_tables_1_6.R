# =============================================================================
# Reproducible runner for manuscript Tables 1--6
# =============================================================================

runner_location <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Run this file with Rscript so its location can be determined.")
  }

  # Some Rscript builds encode spaces in --file= as the literal token "~+~".
  # Decode only when the raw path does not exist and the decoded path does.
  runner_file <- sub("^--file=", "", file_arg[[1L]])
  decoded_file <- gsub("~+~", " ", runner_file, fixed = TRUE)
  if (!file.exists(runner_file) && file.exists(decoded_file)) {
    runner_file <- decoded_file
  }

  dirname(normalizePath(
    runner_file,
    winslash = "/",
    mustWork = TRUE
  ))
}


print_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript run_tables_1_6.R [options]",
    "",
    "Options:",
    "  --mode=full|quick       full: 2,000 replications; quick: 20 (default: full)",
    "  --tables=1,2,...,6      run only the listed tables (default: all)",
    "  --output-dir=PATH       output root (default: results/tables_1_6)",
    "  --overwrite             permit replacement of an existing table CSV",
    "  --help                  show this message",
    "",
    "Each table runs in a clean R process. Quick results are smoke tests and are",
    "never manuscript results.",
    sep = "\n"
  ))
}


parse_runner_args <- function(args) {
  options <- list(
    mode = "full",
    tables = 1:6,
    output_dir = NULL,
    overwrite = FALSE,
    help = FALSE
  )

  for (arg in args) {
    if (identical(arg, "--help")) {
      options$help <- TRUE
    } else if (identical(arg, "--overwrite")) {
      options$overwrite <- TRUE
    } else if (grepl("^--mode=", arg)) {
      options$mode <- sub("^--mode=", "", arg)
    } else if (grepl("^--tables=", arg)) {
      value <- sub("^--tables=", "", arg)
      pieces <- strsplit(value, ",", fixed = TRUE)[[1L]]
      options$tables <- suppressWarnings(as.integer(pieces))
    } else if (grepl("^--output-dir=", arg)) {
      options$output_dir <- sub("^--output-dir=", "", arg)
    } else {
      stop("Unknown argument: ", arg, "\nUse --help for usage.")
    }
  }

  if (!options$mode %in% c("full", "quick")) {
    stop("--mode must be 'full' or 'quick'.")
  }

  if (length(options$tables) == 0L || anyNA(options$tables) ||
      any(!options$tables %in% 1:6)) {
    stop("--tables must be a comma-separated subset of 1,2,3,4,5,6.")
  }

  options$tables <- unique(options$tables)
  options
}


absolute_path <- function(path, relative_to = getwd()) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(file.path(relative_to, path), winslash = "/", mustWork = FALSE)
}


write_run_metadata <- function(output_dir, plan, mode, replications) {
  source_paths <- file.path(runner_location(), plan$Script)
  source_hash <- unname(tools::md5sum(source_paths))
  output_relative <- file.path(
    sprintf("table_%d", plan$Table),
    plan$Output
  )
  output_paths <- file.path(output_dir, output_relative)
  output_hash <- rep(NA_character_, length(output_paths))
  exists <- file.exists(output_paths)
  output_hash[exists] <- unname(tools::md5sum(output_paths[exists]))
  manuscript_relative <- file.path(
    sprintf("table_%d", plan$Table),
    sprintf("Table_%d_manuscript.csv", plan$Table)
  )
  manuscript_paths <- file.path(output_dir, manuscript_relative)
  manuscript_hash <- rep(NA_character_, length(manuscript_paths))
  manuscript_exists <- file.exists(manuscript_paths)
  manuscript_hash[manuscript_exists] <- unname(tools::md5sum(manuscript_paths[manuscript_exists]))
  quadrature_relative <- ifelse(
    !is.na(plan$Quadrature_output) & identical(mode, "full"),
    file.path(sprintf("table_%d", plan$Table), plan$Quadrature_output),
    NA_character_
  )
  quadrature_paths <- file.path(output_dir, quadrature_relative)
  quadrature_hash <- rep(NA_character_, length(quadrature_paths))
  quadrature_exists <- !is.na(quadrature_relative) & file.exists(quadrature_paths)
  quadrature_hash[quadrature_exists] <- unname(tools::md5sum(
    quadrature_paths[quadrature_exists]
  ))

  manifest <- data.frame(
    Table = plan$Table,
    Description = plan$Description,
    Mode = mode,
    Replications = replications,
    Seed = plan$Seed,
    Quadrature_nodes = plan$N_quad,
    R_version = as.character(getRversion()),
    emplik_version = if (requireNamespace("emplik", quietly = TRUE)) {
      as.character(utils::packageVersion("emplik"))
    } else {
      NA_character_
    },
    Script = plan$Script,
    Source_MD5 = source_hash,
    Output = output_relative,
    Output_MD5 = output_hash,
    Manuscript_Output = manuscript_relative,
    Manuscript_Output_MD5 = manuscript_hash,
    Quadrature_Check = quadrature_relative,
    Quadrature_Check_MD5 = quadrature_hash,
    stringsAsFactors = FALSE
  )

  manifest_path <- file.path(output_dir, "run_manifest.csv")
  if (file.exists(manifest_path)) {
    existing <- utils::read.csv(manifest_path, check.names = FALSE)
    if (!identical(names(existing), names(manifest))) {
      stop(
        "The existing run_manifest.csv uses an incompatible schema. ",
        "Move it or choose a new --output-dir."
      )
    }
    existing <- existing[!existing$Table %in% manifest$Table, , drop = FALSE]
    manifest <- rbind(existing, manifest)
    manifest <- manifest[order(manifest$Table), , drop = FALSE]
    rownames(manifest) <- NULL
  }
  utils::write.csv(manifest, manifest_path, row.names = FALSE)

  session <- c(
    sprintf("Run completed (UTC): %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    sprintf("Mode: %s", mode),
    sprintf("Replications per scenario: %d", replications),
    capture.output(sessionInfo())
  )
  writeLines(session, file.path(output_dir, "sessionInfo.txt"), useBytes = TRUE)

  invisible(manifest)
}


format_manuscript_table <- function(table_number, raw_csv, table_dir) {
  raw <- utils::read.csv(raw_csv, check.names = FALSE)
  output <- file.path(table_dir, sprintf("Table_%d_manuscript.csv", table_number))

  format_fixed <- function(values, digits) {
    formatted <- rep(NA_character_, length(values))
    finite <- is.finite(values)
    formatted[finite] <- formatC(
      values[finite],
      format = "f",
      digits = digits
    )
    formatted
  }

  if (table_number == 1L) {
    parameter_columns <- c("lambda1", "lambda2")
  } else if (table_number == 2L) {
    parameter_columns <- c("shape1", "scale1", "shape2", "scale2")
  } else if (table_number == 3L) {
    parameter_columns <- c("lambda1", "lambda2")
  } else if (table_number == 4L) {
    parameter_columns <- c("shape1", "scale1", "shape2", "scale2")
  } else if (table_number == 5L) {
    parameter_columns <- c("lambda1", "lambda2")
  } else if (table_number == 6L) {
    parameter_columns <- c("shape1", "scale1", "shape2", "scale2")
  } else {
    stop("Unknown table number: ", table_number)
  }

  make_wide <- function(data, parameters, labels_to_columns, label_name, digits) {
    parameter_grid <- unique(data[, parameters, drop = FALSE])
    sample_grid <- unique(data[, c("n1", "n2"), drop = FALSE])
    sample_names <- sprintf("n_%d_%d", sample_grid$n1, sample_grid$n2)
    rows <- vector("list", nrow(parameter_grid) * length(labels_to_columns))
    row_index <- 1L

    for (parameter_index in seq_len(nrow(parameter_grid))) {
      keep <- rep(TRUE, nrow(data))
      for (parameter in parameters) {
        keep <- keep & data[[parameter]] == parameter_grid[[parameter]][[parameter_index]]
      }
      block <- data[keep, , drop = FALSE]

      for (label in names(labels_to_columns)) {
        values <- rep(NA_real_, nrow(sample_grid))
        for (sample_index in seq_len(nrow(sample_grid))) {
          value_index <- which(
            block$n1 == sample_grid$n1[[sample_index]] &
              block$n2 == sample_grid$n2[[sample_index]]
          )
          if (length(value_index) != 1L) {
            stop(sprintf(
              "Expected exactly one row for sample pair (%d,%d).",
              sample_grid$n1[[sample_index]],
              sample_grid$n2[[sample_index]]
            ))
          }
          values[[sample_index]] <- block[[labels_to_columns[[label]]]][[value_index]]
        }

        row <- data.frame(
          parameter_grid[parameter_index, , drop = FALSE],
          label,
          as.list(stats::setNames(format_fixed(values, digits), sample_names)),
          check.names = FALSE
        )
        names(row)[length(parameters) + 1L] <- label_name
        rows[[row_index]] <- row
        row_index <- row_index + 1L
      }
    }

    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
  }

  if (table_number <= 2L) {
    measure_map <- c(
      D_Kernel = "RelMSE_D",
      D_CC = "RelMSE_DCC",
      KL = "RelMSE_KL"
    )
    manuscript <- make_wide(
      raw,
      parameter_columns,
      measure_map,
      "Measure",
      4L
    )
  } else if (table_number <= 4L) {
    estimator_map <- c(
      Kernel = "MSE_Kernel",
      Empirical = "MSE_Emp",
      U_statistic = "MSE_Ustat"
    )
    manuscript <- make_wide(
      raw,
      parameter_columns,
      estimator_map,
      "Estimator",
      4L
    )
  } else {
    manuscript <- data.frame(
      raw[, c(parameter_columns, "n1", "n2"), drop = FALSE],
      JEL_CP_pct = format_fixed(raw$JEL_CP, 1L),
      JEL_AL = format_fixed(raw$JEL_AL, 3L),
      NA_CP_pct = format_fixed(raw$NA_CP, 1L),
      NA_AL = format_fixed(raw$NA_AL, 3L),
      Emp_CP_pct = format_fixed(raw$Emp_CP, 1L),
      Emp_AL = format_fixed(raw$Emp_AL, 3L),
      check.names = FALSE
    )
  }

  rownames(manuscript) <- NULL
  utils::write.csv(manuscript, output, row.names = FALSE, na = "")
  invisible(output)
}


run_one_table <- function(row, output_dir, mode, replications, overwrite) {
  table_number <- row$Table[[1L]]
  table_dir <- file.path(output_dir, sprintf("table_%d", table_number))
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  expected_csv <- file.path(table_dir, row$Output[[1L]])
  manuscript_csv <- file.path(table_dir, sprintf("Table_%d_manuscript.csv", table_number))
  log_path <- file.path(table_dir, "run.log")
  quadrature_csv <- if (identical(mode, "full") && !is.na(row$Quadrature_output[[1L]])) {
    file.path(table_dir, row$Quadrature_output[[1L]])
  } else {
    character()
  }
  protected <- c(expected_csv, manuscript_csv, log_path, quadrature_csv)
  existing <- protected[file.exists(protected)]
  if (length(existing) > 0L && !overwrite) {
    stop(
      sprintf(
        "Table %d artifact(s) already exist:\n%s\nUse --overwrite or choose another --output-dir.",
        table_number,
        paste0("- ", existing, collapse = "\n")
      )
    )
  }

  source_path <- file.path(runner_location(), row$Script[[1L]])
  rscript <- file.path(R.home("bin"), "Rscript")

  env <- c(
    sprintf("EXTROPY_MC_REPS=%d", replications),
    sprintf("EXTROPY_SEED=%d", row$Seed[[1L]]),
    "EXTROPY_OUTPUT_DIR=.",
    sprintf("EXTROPY_VERBOSE=%s", if (identical(mode, "full")) "true" else "false"),
    "TZ=UTC"
  )
  if (is.finite(row$N_quad[[1L]])) {
    env <- c(env, sprintf("EXTROPY_N_QUAD=%d", row$N_quad[[1L]]))
  }
  if (table_number <= 2L) {
    env <- c(
      env,
      sprintf(
        "EXTROPY_QUAD_CHECK=%s",
        if (identical(mode, "full")) "true" else "false"
      )
    )
  }

  cat(sprintf("Table %d: %s\n", table_number, row$Description[[1L]]))
  cat(sprintf("  script: %s\n", row$Script[[1L]]))
  cat(sprintf("  output: %s\n", expected_csv))

  old_dir <- setwd(table_dir)
  on.exit(setwd(old_dir), add = TRUE)

  status <- system2(
    command = rscript,
    args = c("--vanilla", shQuote(source_path)),
    stdout = "run.log",
    stderr = "run.log",
    env = env,
    wait = TRUE
  )

  setwd(old_dir)
  on.exit(NULL, add = FALSE)

  if (length(status) != 1L || !is.finite(status) || status != 0L) {
    stop(sprintf(
      "Table %d failed with exit status %d. See %s",
      table_number,
      status,
      log_path
    ))
  }

  if (!file.exists(expected_csv)) {
    stop(sprintf(
      "Table %d finished but did not create its expected CSV: %s",
      table_number,
      expected_csv
    ))
  }
  if (length(quadrature_csv) == 1L && !file.exists(quadrature_csv)) {
    stop(sprintf(
      "Table %d finished but did not create its quadrature check: %s",
      table_number,
      quadrature_csv
    ))
  }

  invisible(expected_csv)
}


main <- function() {
  options <- parse_runner_args(commandArgs(trailingOnly = TRUE))
  if (options$help) {
    print_usage()
    return(invisible(TRUE))
  }

  base_dir <- runner_location()
  validation_env <- new.env(parent = baseenv())
  sys.source(file.path(base_dir, "validate_tables_1_6.R"), envir = validation_env)
  plan <- validation_env$validate_source_files(base_dir, quiet = TRUE)

  selected_plan <- plan[plan$Table %in% options$tables, , drop = FALSE]
  selected_plan <- selected_plan[match(options$tables, selected_plan$Table), , drop = FALSE]

  if (getRversion() < "4.1.0") {
    stop("This workflow requires R >= 4.1.0.")
  }
  if (getRversion() != "4.6.1") {
    warning(
      "The archival reference environment uses R 4.6.1; this run uses R ",
      as.character(getRversion()),
      ". The exact version will be recorded in the manifest."
    )
  }

  if (any(selected_plan$Table >= 5L) && !requireNamespace("emplik", quietly = TRUE)) {
    stop(
      "Package 'emplik' 1.3-3 is required for Tables 5--6. ",
      "See README.md for the pinned installation command."
    )
  }
  if (any(selected_plan$Table >= 5L) &&
      utils::packageVersion("emplik") != "1.3.3") {
    stop(
      "Tables 5--6 require the pinned emplik 1.3-3 release; installed version is ",
      as.character(utils::packageVersion("emplik")),
      ". See README.md for installation instructions."
    )
  }

  replications <- if (identical(options$mode, "full")) 2000L else 20L
  output_root <- if (is.null(options$output_dir)) {
    file.path(base_dir, "results", "tables_1_6")
  } else {
    absolute_path(options$output_dir)
  }
  output_dir <- file.path(output_root, options$mode)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("Mode: %s (%d replications per scenario)\n", options$mode, replications))
  cat(sprintf("Tables: %s\n", paste(options$tables, collapse = ", ")))
  cat(sprintf("Output directory: %s\n\n", output_dir))

  completed <- logical(nrow(selected_plan))
  error_messages <- rep(NA_character_, nrow(selected_plan))

  for (i in seq_len(nrow(selected_plan))) {
    outcome <- tryCatch(
      {
        run_one_table(
          selected_plan[i, , drop = FALSE],
          output_dir,
          options$mode,
          replications,
          options$overwrite
        )
        validation_env$validate_table_output(
          selected_plan$Table[[i]],
          file.path(
            output_dir,
            sprintf("table_%d", selected_plan$Table[[i]]),
            selected_plan$Output[[i]]
          ),
          selected_plan$Expected_rows[[i]],
          replications
        )
        table_dir <- file.path(
          output_dir,
          sprintf("table_%d", selected_plan$Table[[i]])
        )
        manuscript_csv <- format_manuscript_table(
          selected_plan$Table[[i]],
          file.path(table_dir, selected_plan$Output[[i]]),
          table_dir
        )
        validation_env$validate_manuscript_output(
          selected_plan$Table[[i]],
          manuscript_csv
        )
        if (identical(options$mode, "full") &&
            !is.na(selected_plan$Quadrature_output[[i]])) {
          validation_env$validate_quadrature_output(
            selected_plan$Table[[i]],
            file.path(table_dir, selected_plan$Quadrature_output[[i]]),
            selected_plan$Expected_rows[[i]]
          )
        }
        TRUE
      },
      error = function(error) {
        error_messages[[i]] <<- conditionMessage(error)
        FALSE
      }
    )
    completed[[i]] <- outcome
  }

  write_run_metadata(output_dir, selected_plan, options$mode, replications)

  if (any(!completed)) {
    failures <- sprintf(
      "- Table %d: %s",
      selected_plan$Table[!completed],
      error_messages[!completed]
    )
    stop(paste(c("One or more table runs failed:", failures), collapse = "\n"))
  }

  cat("\nAll selected tables completed and passed output validation.\n")
  cat(sprintf("Manifest: %s\n", file.path(output_dir, "run_manifest.csv")))
  invisible(TRUE)
}


if (sys.nframe() == 0L) {
  main()
}
