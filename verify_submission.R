#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
run_root <- NULL
for (arg in args) {
  if (startsWith(arg, "--run-root=")) {
    run_root <- sub("^--run-root=", "", arg)
  } else if (arg %in% c("--help", "-h")) {
    cat("Usage: Rscript verify_submission.R --run-root=PATH\n")
    quit(save = "no", status = 0L)
  } else {
    stop("Unknown option: ", arg, call. = FALSE)
  }
}

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[startsWith(all_args, "--file=")]
if (length(script_arg) != 1L) {
  stop("Invoke this file with Rscript.", call. = FALSE)
}
root <- dirname(normalizePath(
  gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE),
  winslash = "/", mustWork = TRUE
))
if (is.null(run_root) || !dir.exists(run_root)) {
  stop("--run-root must identify a completed full-run directory.",
       call. = FALSE)
}
run_root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)

checks <- data.frame(check = character(), passed = logical(),
                     details = character(), stringsAsFactors = FALSE)
add <- function(name, passed, details = "") {
  checks[nrow(checks) + 1L, ] <<- list(name, isTRUE(passed), details)
}

required_root <- c(
  "README.txt", "renv.lock", "DATA_PROVENANCE.txt",
  "CLEAN_MACHINE_TEST_RECORD.txt", "MANUSCRIPT_CODE_CROSSWALK.csv",
  "REPRODUCIBILITY_CHECKLIST_STATUS.csv", "MONTE_CARLO_STABILITY_RECORD.txt"
)
for (file in required_root) {
  path <- file.path(root, file)
  add(paste("required file", file), file.exists(path) && !dir.exists(path), path)
}

manifest <- file.path(root, "Real Data Analysis", "Images", "image_manifest.csv")
add("completed image manifest", file.exists(manifest), manifest)

if (file.exists(manifest)) {
  image_record <- tryCatch(
    utils::read.csv(manifest, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(error) error
  )
  if (inherits(image_record, "error")) {
    add("readable image manifest", FALSE, conditionMessage(image_record))
  } else {
    required_manifest_columns <- c(
      "file", "analysis_group", "analysis_index", "source_dataset_class",
      "source_record_id", "source_url", "selection_rule", "retrieval_date",
      "rights_holder", "license", "redistribution_permitted", "original_md5",
      "original_width", "original_height"
    )
    add("image manifest schema",
        all(required_manifest_columns %in% names(image_record)),
        paste(setdiff(required_manifest_columns, names(image_record)),
              collapse = ", "))
    if (all(c("file", "original_md5") %in% names(image_record))) {
      image_paths <- file.path(dirname(manifest), image_record$file)
      images_exist <- file.exists(image_paths) & !dir.exists(image_paths)
      add("all manifested images exist", nrow(image_record) == 9L &&
            all(images_exist) && !anyDuplicated(image_record$file),
          paste(image_record$file[!images_exist], collapse = ", "))
      if (all(images_exist)) {
        observed_md5 <- unname(tools::md5sum(image_paths))
        add("image checksums match manifest",
            identical(tolower(observed_md5),
                      tolower(as.character(image_record$original_md5))), "")
      }
    }
  }
}

text_files <- c(file.path(root, "README.txt"),
                file.path(root, "DATA_PROVENANCE.txt"),
                file.path(root, "CLEAN_MACHINE_TEST_RECORD.txt"),
                file.path(root, "MONTE_CARLO_STABILITY_RECORD.txt"), manifest)
text_files <- text_files[file.exists(text_files)]
if (length(text_files)) {
  text <- unlist(lapply(text_files, readLines, warn = FALSE), use.names = FALSE)
  placeholders <- grep("AUTHOR ACTION REQUIRED|AUTHOR_MUST_COMPLETE|TO_COMPLETE",
                       text, value = TRUE)
  add("no unresolved author placeholders", !length(placeholders),
      paste(head(placeholders, 5L), collapse = " | "))
}

paper_dir <- file.path(run_root, "paper_outputs")
tables <- file.path(paper_dir, paste0("Table_", 1:14, "_manuscript.csv"))
figures <- file.path(paper_dir, c(
  "Figure_1a_KM.pdf", "Figure_1b_cumulative.pdf",
  "Figure_2a_KM.pdf", "Figure_2b_cumulative.pdf",
  "Figure_3a_KM.pdf", "Figure_3b_cumulative.pdf", "Figure_4.pdf"
))
for (path in c(tables, figures)) {
  good <- file.exists(path) && !dir.exists(path) && file.info(path)$size > 0
  add(paste("paper output", basename(path)), good, path)
}

environment_files <- file.path(run_root, "environment",
                               c("sessionInfo.txt", "environment.txt",
                                 "required_package_versions.csv"))
for (path in environment_files) {
  add(paste("environment record", basename(path)),
      file.exists(path) && file.info(path)$size > 0, path)
}

full_results <- c(
  file.path(run_root, "uncensored", "full", "run_manifest.csv"),
  file.path(run_root, "censored", "details", "Tables7_10_manifest.csv"),
  file.path(run_root, "MASTER_RUN_LOG.csv"),
  file.path(run_root, "MASTER_RUN_SETTINGS.txt"),
  file.path(run_root, "MASTER_FILE_MANIFEST.csv")
)
for (path in full_results) {
  add(paste("full-run record", basename(path)),
      file.exists(path) && file.info(path)$size > 0, path)
}

r_files <- list.files(root, pattern = "[.]R$", recursive = TRUE,
                      full.names = TRUE)
r_text <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
add("no automatic package installation",
    !any(grepl("^[[:space:]]*install[.]packages[(]", r_text)), "")
add("no global-workspace deletion",
    !any(grepl("rm[[:space:]]*[(][[:space:]]*list[[:space:]]*=[[:space:]]*ls",
               r_text)), "")

report <- file.path(run_root, "SUBMISSION_VERIFICATION.csv")
utils::write.csv(checks, report, row.names = FALSE)
failed <- checks[!checks$passed, , drop = FALSE]
if (nrow(failed)) {
  print(failed, row.names = FALSE)
  stop(nrow(failed), " submission verification check(s) failed. See ", report,
       call. = FALSE)
}
cat("All submission verification checks passed. Report: ", report, "\n",
    sep = "")
