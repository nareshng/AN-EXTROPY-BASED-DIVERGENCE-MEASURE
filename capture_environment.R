#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- "environment"
strict <- FALSE
for (arg in args) {
  if (startsWith(arg, "--output-dir=")) {
    out_dir <- sub("^--output-dir=", "", arg)
  } else if (identical(arg, "--strict")) {
    strict <- TRUE
  } else if (arg %in% c("--help", "-h")) {
    cat("Usage: Rscript capture_environment.R [--strict] ",
        "[--output-dir=PATH]\n", sep = "")
    quit(save = "no", status = 0L)
  } else {
    stop("Unknown option: ", arg, call. = FALSE)
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
required <- c(
  emplik = "1.3-3",
  survival = "3.5-8",
  `TH.data` = "1.1-2",
  jpeg = "0.1-10",
  ggplot2 = "3.4.4",
  survminer = "0.4.9"
)

installed <- vapply(names(required), requireNamespace, logical(1L),
                    quietly = TRUE)
if (any(!installed)) {
  stop("Missing required package(s): ",
       paste(names(required)[!installed], collapse = ", "), call. = FALSE)
}

actual <- vapply(names(required), function(package) {
  as.character(utils::packageVersion(package))
}, character(1L))
version_ok <- actual == required
version_ok[names(required) == "emplik"] <-
  actual[names(required) == "emplik"] %in% c("1.3-3", "1.3.3")

versions <- data.frame(
  package = names(required), required = unname(required), actual = actual,
  matches = version_ok, stringsAsFactors = FALSE
)
utils::write.csv(versions, file.path(out_dir, "required_package_versions.csv"),
                 row.names = FALSE)

if (strict && getRversion() != "4.6.1") {
  stop("Strict reproduction requires R 4.6.1; found ", getRversion(),
       call. = FALSE)
}
if (strict && any(!version_ok)) {
  bad <- versions[!versions$matches, , drop = FALSE]
  stop("Package-version mismatch: ",
       paste(paste0(bad$package, " required ", bad$required,
                    ", found ", bad$actual), collapse = "; "),
       call. = FALSE)
}

for (package in names(required)) {
  suppressPackageStartupMessages(
    library(package, character.only = TRUE)
  )
}

session_file <- file.path(out_dir, "sessionInfo.txt")
sink(session_file)
print(sessionInfo())
sink()

environment_lines <- c(
  paste("Captured UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R:", R.version.string),
  paste("Platform:", R.version$platform),
  paste("OS:", paste(names(Sys.info()), Sys.info(), sep = "=", collapse = "; ")),
  paste("Locale:", paste(Sys.getlocale(), collapse = "; ")),
  paste("Working directory:", normalizePath(getwd(), winslash = "/"))
)
writeLines(environment_lines, file.path(out_dir, "environment.txt"),
           useBytes = TRUE)

cat("Environment recorded in ", normalizePath(out_dir, winslash = "/"),
    "\n", sep = "")

