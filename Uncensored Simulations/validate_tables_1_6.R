# =============================================================================
# Static and output validation for manuscript Tables 1--6
# =============================================================================

script_location <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }

  script_file <- sub("^--file=", "", file_arg[[1L]])
  decoded_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  if (!file.exists(script_file) && file.exists(decoded_file)) {
    script_file <- decoded_file
  }

  dirname(normalizePath(
    script_file,
    winslash = "/",
    mustWork = TRUE
  ))
}


tables_1_6_plan <- function(base_dir = script_location()) {
  data.frame(
    Table = 1:6,
    Description = c(
      "Kernel-divergence comparison, exponential",
      "Kernel-divergence comparison, Weibull",
      "Estimator comparison, exponential",
      "Estimator comparison, Weibull",
      "Confidence intervals, exponential",
      "Confidence intervals, Weibull"
    ),
    Script = c(
      "Sim_RelativeMSE_comparison_Exponential_dist.R",
      "Sim_RelativeMSE_comparison_Weibull_dist.R",
      "Sim_MSE_Estimators_Exponential_dist.R",
      "Sim_MSE_Estimators_Weibull_dist.R",
      "Sim_Confidence_Intervals_Exponential_dist.R",
      "Sim_Confidence_Intervals_Weibull_dist.R"
    ),
    Output = c(
      "MSE_and_Relative_MSE_Section_5_1.csv",
      "MSE_and_Relative_MSE_Section_5_1_Weibull.csv",
      "MSE_and_Relative_MSE_results.csv",
      "MSE_and_Relative_MSE_results_Weibull.csv",
      "Coverage_and_Average_Length_Section_5_3_Exponential_Efficient_Corrected.csv",
      "Coverage_and_Average_Length_Section_5_3_Weibull.csv"
    ),
    Seed = c(2024L, 2024L, 2026L, 2024L, 2026L, 2026L),
    N_quad = c(200L, 200L, 80L, 80L, NA_integer_, NA_integer_),
    Quadrature_output = c(
      "Quadrature_Check_Section_5_1.csv",
      "Quadrature_Check_Section_5_1_Weibull.csv",
      NA_character_, NA_character_, NA_character_, NA_character_
    ),
    Expected_rows = c(24L, 24L, 42L, 42L, 20L, 24L),
    stringsAsFactors = FALSE
  )
}


required_columns <- function(table_number) {
  switch(
    as.character(table_number),
    `1` = c(
      "lambda1", "lambda2", "n1", "n2", "true_D", "true_DCC", "true_KL",
      "MSE_D", "MSE_DCC", "MSE_KL", "RelMSE_D", "RelMSE_DCC",
      "RelMSE_KL", "NA_D", "NA_DCC", "NA_KL"
    ),
    `2` = c(
      "shape1", "scale1", "shape2", "scale2", "n1", "n2", "true_D",
      "true_DCC", "true_KL", "MSE_D", "MSE_DCC", "MSE_KL", "RelMSE_D",
      "RelMSE_DCC", "RelMSE_KL", "NA_D", "NA_DCC", "NA_KL"
    ),
    `3` = c(
      "lambda1", "lambda2", "true_D", "n1", "n2", "MSE_Kernel",
      "MSE_Emp", "MSE_Ustat", "RelMSE_Kernel", "RelMSE_Emp",
      "RelMSE_Ustat", "Kernel_NA_Count"
    ),
    `4` = c(
      "shape1", "scale1", "shape2", "scale2", "true_D", "n1", "n2",
      "MSE_Kernel", "MSE_Emp", "MSE_Ustat", "RelMSE_Kernel",
      "RelMSE_Emp", "RelMSE_Ustat", "Kernel_NA_Count"
    ),
    `5` = c(
      "lambda1", "lambda2", "true_D", "n1", "n2", "JEL_CP", "JEL_AL",
      "NA_CP", "NA_AL", "Emp_CP", "Emp_AL", "Valid_JEL", "Valid_NA",
      "Valid_Emp", "Failed_JEL", "Failed_NA", "Failed_Emp", "Covered_JEL",
      "Covered_NA", "Covered_Emp"
    ),
    `6` = c(
      "shape1", "scale1", "shape2", "scale2", "true_D", "n1", "n2",
      "JEL_CP", "JEL_AL", "NA_CP", "NA_AL", "Emp_CP", "Emp_AL",
      "Valid_JEL", "Valid_NA", "Valid_Emp", "Failed_JEL", "Failed_NA",
      "Failed_Emp", "Covered_JEL", "Covered_NA", "Covered_Emp"
    ),
    stop("Unknown table number: ", table_number)
  )
}


expected_scenarios <- function(table_number) {
  sample_grid <- if (table_number %in% 1:2) {
    data.frame(n1 = c(10, 50, 80, 200), n2 = c(10, 40, 100, 200))
  } else if (table_number %in% 3:4) {
    data.frame(
      n1 = c(5, 10, 20, 20, 40, 30, 50),
      n2 = c(10, 10, 10, 20, 30, 40, 50)
    )
  } else if (table_number == 5L) {
    data.frame(n1 = c(20, 30, 70, 100), n2 = c(10, 40, 50, 100))
  } else if (table_number == 6L) {
    data.frame(n1 = c(10, 30, 70, 100), n2 = c(10, 40, 50, 100))
  } else {
    stop("Unknown table number: ", table_number)
  }

  parameter_grid <- if (table_number == 1L) {
    data.frame(
      lambda1 = c(0.1, 0.5, 2, 3, 7, 10),
      lambda2 = c(0.2, 1, 1, 5, 5, 7)
    )
  } else if (table_number %in% c(2L, 6L)) {
    data.frame(
      shape1 = c(0.5, 1, 2, 1.5, 0.7, 3),
      scale1 = c(1, 1, 1, 1, 2, 1),
      shape2 = c(1, 2, 1, 3, 1.5, 1.2),
      scale2 = c(1, 1, 1, 1.5, 1, 2)
    )
  } else if (table_number == 3L) {
    data.frame(
      lambda1 = c(0.5, 0.1, 0.1, 1, 2, 10),
      lambda2 = c(0.1, 0.5, 1, 0.5, 5, 5)
    )
  } else if (table_number == 4L) {
    data.frame(
      shape1 = c(0.5, 1, 1.5, 2, 0.7, 3),
      scale1 = c(1, 1, 1, 1, 2, 1),
      shape2 = c(1, 2, 3, 2, 1.5, 1.2),
      scale2 = c(1, 1, 1, 2, 1, 2)
    )
  } else if (table_number == 5L) {
    data.frame(
      lambda1 = c(0.2, 0.1, 1, 1, 10),
      lambda2 = c(0.1, 0.5, 0.5, 2, 1)
    )
  } else {
    stop("Unknown table number: ", table_number)
  }

  rows <- lapply(seq_len(nrow(parameter_grid)), function(i) {
    cbind(
      parameter_grid[rep(i, nrow(sample_grid)), , drop = FALSE],
      sample_grid,
      row.names = NULL
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


row_keys <- function(data) {
  if (nrow(data) == 0L) {
    return(character())
  }
  do.call(paste, c(lapply(data, as.character), sep = "\r"))
}


assert_exact_scenarios <- function(table_number, result) {
  expected <- expected_scenarios(table_number)
  actual <- result[, names(expected), drop = FALSE]
  if (!identical(sort(row_keys(actual)), sort(row_keys(expected)))) {
    stop(sprintf(
      "Table %d output does not contain the exact manuscript scenario grid.",
      table_number
    ))
  }
  invisible(TRUE)
}


validate_source_files <- function(base_dir = script_location(), quiet = FALSE) {
  plan <- tables_1_6_plan(base_dir)
  failures <- character()

  required_env <- c(
    "EXTROPY_MC_REPS",
    "EXTROPY_SEED",
    "EXTROPY_OUTPUT_DIR",
    "EXTROPY_VERBOSE"
  )

  forbidden_portability <- c(
    "(?m)^[[:space:]]*setwd\\s*\\(",
    "(?m)^[[:space:]]*install\\.packages\\s*\\(",
    "/Users/",
    "[A-Za-z]:[/\\\\]Users[/\\\\]"
  )

  for (i in seq_len(nrow(plan))) {
    table_number <- plan$Table[[i]]
    path <- file.path(base_dir, plan$Script[[i]])

    if (!file.exists(path)) {
      failures <- c(failures, sprintf("Table %d script is missing: %s", table_number, path))
      next
    }

    parsed <- tryCatch(parse(file = path), error = identity)
    if (inherits(parsed, "error")) {
      failures <- c(
        failures,
        sprintf("Table %d script does not parse: %s", table_number, conditionMessage(parsed))
      )
      next
    }

    text <- paste(readLines(path, warn = FALSE), collapse = "\n")

    missing_env <- required_env[!vapply(required_env, grepl, logical(1L), x = text, fixed = TRUE)]
    if (length(missing_env) > 0L) {
      failures <- c(
        failures,
        sprintf(
          "Table %d script is missing runner interface variable(s): %s",
          table_number,
          paste(missing_env, collapse = ", ")
        )
      )
    }

    for (pattern in forbidden_portability) {
      if (grepl(pattern, text, perl = TRUE)) {
        failures <- c(
          failures,
          sprintf("Table %d script contains non-portable pattern: %s", table_number, pattern)
        )
      }
    }

    if (!grepl(plan$Output[[i]], text, fixed = TRUE)) {
      failures <- c(
        failures,
        sprintf("Table %d script does not name its expected CSV: %s", table_number, plan$Output[[i]])
      )
    }

    if (table_number >= 5L) {
      kernel_ci_tokens <- c(
        "kernel_bootstrap_ci", "make_kernel_object", "B_boot", "Ker_CP",
        "Ker_AL", "Valid_Ker"
      )
      present <- kernel_ci_tokens[vapply(kernel_ci_tokens, grepl, logical(1L), x = text, fixed = TRUE)]
      if (length(present) > 0L) {
        failures <- c(
          failures,
          sprintf(
            "Table %d still contains kernel-CI token(s): %s",
            table_number,
            paste(present, collapse = ", ")
          )
        )
      }

      compact <- gsub("[[:space:]]+", "", text)
      jel_signatures <- c(
        "cx<-(n/(n-2))*((n2-1)*(2/n1)-1)",
        "cy<-(n/(n-2))*((n1-1)*(2/n2)-1)",
        "V_JEL_X<-n*sc$D_Ustat-(n-1)*loo$D_Ustat_loo_X",
        "V_JEL_Y<-n*sc$D_Ustat-(n-1)*loo$D_Ustat_loo_Y",
        "emplik::el.test(w,mu=0)"
      )
      missing_jel <- jel_signatures[
        !vapply(jel_signatures, grepl, logical(1L), x = compact, fixed = TRUE)
      ]
      if (length(missing_jel) > 0L) {
        failures <- c(
          failures,
          sprintf(
            "Table %d no longer contains the retained JEL implementation signature(s): %s",
            table_number,
            paste(missing_jel, collapse = "; ")
          )
        )
      }
    }
  }

  if (length(failures) > 0L) {
    stop(paste(c("Tables 1--6 source validation failed:", paste0("- ", failures)), collapse = "\n"))
  }

  if (!quiet) {
    cat(sprintf("Source validation passed for %d scripts.\n", nrow(plan)))
  }

  invisible(plan)
}


validate_table_output <- function(
    table_number,
    csv_path,
    expected_rows = NULL,
    expected_replications = NULL
) {
  if (!file.exists(csv_path)) {
    stop(sprintf("Table %d output is missing: %s", table_number, csv_path))
  }

  result <- utils::read.csv(csv_path, check.names = FALSE)
  required <- required_columns(table_number)
  missing <- setdiff(required, names(result))

  if (length(missing) > 0L) {
    stop(sprintf(
      "Table %d output is missing column(s): %s",
      table_number,
      paste(missing, collapse = ", ")
    ))
  }

  if (!is.null(expected_rows) && nrow(result) != expected_rows) {
    stop(sprintf(
      "Table %d output has %d rows; expected %d.",
      table_number,
      nrow(result),
      expected_rows
    ))
  }

  assert_exact_scenarios(table_number, result)

  key_columns <- intersect(
    c("lambda1", "lambda2", "shape1", "scale1", "shape2", "scale2", "n1", "n2"),
    names(result)
  )
  if (anyDuplicated(result[, key_columns, drop = FALSE])) {
    stop(sprintf("Table %d output contains duplicate scenario rows.", table_number))
  }

  if (any(!is.finite(result$n1)) || any(!is.finite(result$n2)) ||
      any(result$n1 <= 0) || any(result$n2 <= 0)) {
    stop(sprintf("Table %d output contains invalid sample sizes.", table_number))
  }

  if (table_number <= 4L) {
    metric_columns <- grep("^(MSE|RelMSE)", names(result), value = TRUE)
    for (column in metric_columns) {
      values <- result[[column]]
      if (any(!is.finite(values)) || any(values < 0)) {
        stop(sprintf(
          "Table %d output contains a nonfinite or negative value in %s.",
          table_number,
          column
        ))
      }
    }

    failure_columns <- intersect(
      c("NA_D", "NA_DCC", "NA_KL", "Kernel_NA_Count"),
      names(result)
    )
    for (column in failure_columns) {
      values <- result[[column]]
      if (any(!is.finite(values)) || any(values != round(values)) || any(values != 0)) {
        stop(sprintf(
          "Table %d publication output has a nonzero or invalid failure count in %s.",
          table_number,
          column
        ))
      }
    }
  } else {
    if (any(grepl("Ker|Kernel", names(result), ignore.case = TRUE))) {
      stop(sprintf("Table %d output still contains a kernel-CI column.", table_number))
    }

    cp_columns <- grep("_CP$", names(result), value = TRUE)
    al_columns <- grep("_AL$", names(result), value = TRUE)

    for (column in cp_columns) {
      values <- result[[column]]
      if (any(!is.finite(values)) || any(values < 0 | values > 100)) {
        stop(sprintf("Table %d output has coverage outside [0,100] in %s.", table_number, column))
      }
    }

    for (column in al_columns) {
      method <- sub("_AL$", "", column)
      valid <- result[[paste0("Valid_", method)]]
      values <- result[[column]]
      if (any(values[is.finite(values)] < 0) || any(valid > 0 & !is.finite(values))) {
        stop(sprintf(
          "Table %d output contains an invalid interval length in %s.",
          table_number,
          column
        ))
      }
    }

    methods <- c("JEL", "NA", "Emp")
    totals <- lapply(methods, function(method) {
      valid <- result[[paste0("Valid_", method)]]
      failed <- result[[paste0("Failed_", method)]]
      covered <- result[[paste0("Covered_", method)]]
      counts <- c(valid, failed, covered)
      if (any(!is.finite(counts)) || any(counts < 0) ||
          any(abs(counts - round(counts)) > 1e-12)) {
        stop(sprintf(
          "Table %d contains invalid valid/failed counts for %s.",
          table_number,
          method
        ))
      }
      if (any(covered > valid)) {
        stop(sprintf(
          "Table %d has Covered_%s greater than Valid_%s.",
          table_number,
          method,
          method
        ))
      }
      valid + failed
    })
    totals <- do.call(cbind, totals)

    if (any(totals != totals[[1L]])) {
      stop(sprintf(
        "Table %d valid plus failed counts do not use a common planned denominator.",
        table_number
      ))
    }
    if (!is.null(expected_replications) && any(totals != expected_replications)) {
      stop(sprintf(
        "Table %d valid plus failed counts do not equal %d planned replications.",
        table_number,
        expected_replications
      ))
    }

    planned <- totals[, 1L]
    for (method in methods) {
      expected_cp <- 100 * result[[paste0("Covered_", method)]] / planned
      actual_cp <- result[[paste0(method, "_CP")]]
      if (any(abs(actual_cp - expected_cp) > 1e-10)) {
        stop(sprintf(
          "Table %d %s coverage does not equal 100 * Covered / planned replications.",
          table_number,
          method
        ))
      }
    }
  }

  invisible(result)
}


validate_manuscript_output <- function(table_number, csv_path) {
  if (!file.exists(csv_path)) {
    stop(sprintf("Table %d manuscript CSV is missing: %s", table_number, csv_path))
  }

  result <- utils::read.csv(csv_path, check.names = FALSE)

  if (table_number <= 2L) {
    parameter_columns <- if (table_number == 1L) {
      c("lambda1", "lambda2")
    } else {
      c("shape1", "scale1", "shape2", "scale2")
    }
    expected_names <- c(
      parameter_columns,
      "Measure", "n_10_10", "n_50_40", "n_80_100", "n_200_200"
    )
    expected_rows <- 18L
    rounded_columns <- grep("^n_", names(result), value = TRUE)
    digits <- 4L
  } else if (table_number <= 4L) {
    parameter_columns <- if (table_number == 3L) {
      c("lambda1", "lambda2")
    } else {
      c("shape1", "scale1", "shape2", "scale2")
    }
    expected_names <- c(
      parameter_columns,
      "Estimator", "n_5_10", "n_10_10", "n_20_10", "n_20_20",
      "n_40_30", "n_30_40", "n_50_50"
    )
    expected_rows <- 18L
    rounded_columns <- grep("^n_", names(result), value = TRUE)
    digits <- 4L
  } else {
    parameter_columns <- if (table_number == 5L) {
      c("lambda1", "lambda2")
    } else {
      c("shape1", "scale1", "shape2", "scale2")
    }
    expected_names <- c(
      parameter_columns,
      "n1", "n2", "JEL_CP_pct", "JEL_AL", "NA_CP_pct", "NA_AL",
      "Emp_CP_pct", "Emp_AL"
    )
    expected_rows <- if (table_number == 5L) 20L else 24L
    rounded_columns <- c("JEL_AL", "NA_AL", "Emp_AL")
    digits <- 3L
  }

  if (!identical(names(result), expected_names)) {
    stop(sprintf(
      "Table %d manuscript CSV columns differ from the exact contract. Expected: %s",
      table_number,
      paste(expected_names, collapse = ", ")
    ))
  }

  if (nrow(result) != expected_rows) {
    stop(sprintf(
      "Table %d manuscript CSV has %d rows; expected %d.",
      table_number,
      nrow(result),
      expected_rows
    ))
  }

  if (table_number <= 4L) {
    scenario_grid <- unique(expected_scenarios(table_number)[, parameter_columns, drop = FALSE])
    labels <- if (table_number <= 2L) {
      c("D_Kernel", "D_CC", "KL")
    } else {
      c("Kernel", "Empirical", "U_statistic")
    }
    label_name <- if (table_number <= 2L) "Measure" else "Estimator"
    expected_keys <- do.call(rbind, lapply(seq_len(nrow(scenario_grid)), function(i) {
      data.frame(
        scenario_grid[rep(i, length(labels)), , drop = FALSE],
        stats::setNames(list(labels), label_name),
        check.names = FALSE,
        row.names = NULL
      )
    }))
    actual_keys <- result[, c(parameter_columns, label_name), drop = FALSE]
    if (!identical(row_keys(actual_keys), row_keys(expected_keys))) {
      stop(sprintf(
        "Table %d manuscript CSV has an incorrect parameter/row-label grid or order.",
        table_number
      ))
    }
  } else {
    expected <- expected_scenarios(table_number)
    actual <- result[, names(expected), drop = FALSE]
    if (!identical(row_keys(actual), row_keys(expected))) {
      stop(sprintf(
        "Table %d manuscript CSV has an incorrect scenario grid or order.",
        table_number
      ))
    }
  }

  for (column in rounded_columns) {
    values <- result[[column]]
    if (any(!is.finite(values)) ||
        any(abs(values - round(values, digits)) > 1e-12)) {
      stop(sprintf(
        "Table %d manuscript column %s is not rounded to %d decimals.",
        table_number,
        column,
        digits
      ))
    }
  }

  if (table_number >= 5L) {
    cp_columns <- c("JEL_CP_pct", "NA_CP_pct", "Emp_CP_pct")
    for (column in cp_columns) {
      values <- result[[column]]
      if (any(!is.finite(values)) || any(values < 0 | values > 100) ||
          any(abs(values - round(values, 1L)) > 1e-12)) {
        stop(sprintf(
          "Table %d manuscript column %s must be a percentage rounded to one decimal.",
          table_number,
          column
        ))
      }
    }

    if (any(grepl("Ker|Kernel", names(result), ignore.case = TRUE))) {
      stop(sprintf("Table %d manuscript CSV contains a kernel-CI column.", table_number))
    }
  }

  invisible(result)
}


validate_quadrature_output <- function(table_number, csv_path, expected_rows) {
  if (!table_number %in% 1:2) {
    stop("Quadrature validation is defined only for Tables 1--2.")
  }
  if (!file.exists(csv_path)) {
    stop(sprintf("Table %d quadrature check is missing: %s", table_number, csv_path))
  }

  result <- utils::read.csv(csv_path, check.names = FALSE)
  scenario <- expected_scenarios(table_number)
  difference_columns <- c("Diff_RelMSE_D", "Diff_RelMSE_DCC", "Diff_RelMSE_KL")
  expected_names <- c(names(scenario), "n_quad_coarse", "n_quad_fine", difference_columns)

  if (!identical(names(result), expected_names)) {
    stop(sprintf("Table %d quadrature-check columns differ from the exact contract.", table_number))
  }
  if (nrow(result) != expected_rows) {
    stop(sprintf(
      "Table %d quadrature check has %d rows; expected %d.",
      table_number,
      nrow(result),
      expected_rows
    ))
  }
  if (!identical(row_keys(result[, names(scenario), drop = FALSE]), row_keys(scenario))) {
    stop(sprintf("Table %d quadrature check uses an incorrect scenario grid.", table_number))
  }
  if (any(!is.finite(result$n_quad_coarse)) || any(!is.finite(result$n_quad_fine)) ||
      any(result$n_quad_fine <= result$n_quad_coarse)) {
    stop(sprintf("Table %d quadrature check has invalid node counts.", table_number))
  }
  for (column in difference_columns) {
    if (any(!is.finite(result[[column]])) || any(result[[column]] < 0)) {
      stop(sprintf(
        "Table %d quadrature check has an invalid value in %s.",
        table_number,
        column
      ))
    }
  }

  invisible(result)
}


validate_output_tree <- function(output_dir, base_dir = script_location(), quiet = FALSE) {
  plan <- tables_1_6_plan(base_dir)
  manifest_path <- file.path(output_dir, "run_manifest.csv")
  manifest <- if (file.exists(manifest_path)) {
    utils::read.csv(manifest_path, check.names = FALSE)
  } else {
    NULL
  }

  for (i in seq_len(nrow(plan))) {
    csv_path <- file.path(
      output_dir,
      sprintf("table_%d", plan$Table[[i]]),
      plan$Output[[i]]
    )
    expected_replications <- NULL
    if (!is.null(manifest) && all(c("Table", "Replications") %in% names(manifest))) {
      matched <- manifest$Replications[manifest$Table == plan$Table[[i]]]
      if (length(matched) == 1L) {
        expected_replications <- matched
      }
    }
    validate_table_output(
      plan$Table[[i]],
      csv_path,
      plan$Expected_rows[[i]],
      expected_replications
    )
    validate_manuscript_output(
      plan$Table[[i]],
      file.path(
        output_dir,
        sprintf("table_%d", plan$Table[[i]]),
        sprintf("Table_%d_manuscript.csv", plan$Table[[i]])
      )
    )
    if (!is.null(manifest) && all(c("Table", "Mode") %in% names(manifest))) {
      mode <- manifest$Mode[manifest$Table == plan$Table[[i]]]
      if (length(mode) == 1L && identical(mode, "full") &&
          !is.na(plan$Quadrature_output[[i]])) {
        validate_quadrature_output(
          plan$Table[[i]],
          file.path(
            output_dir,
            sprintf("table_%d", plan$Table[[i]]),
            plan$Quadrature_output[[i]]
          ),
          plan$Expected_rows[[i]]
        )
      }
    }
  }

  if (!quiet) {
    cat("Output validation passed for Tables 1--6.\n")
  }

  invisible(TRUE)
}


if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_arg <- grep("^--output-dir=", args, value = TRUE)

  validate_source_files()

  if (length(output_arg) == 1L) {
    validate_output_tree(sub("^--output-dir=", "", output_arg[[1L]]))
  } else if (length(output_arg) > 1L) {
    stop("Use --output-dir only once.")
  }
}
