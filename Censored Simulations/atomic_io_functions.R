## Cross-platform atomic file writers used by censored simulations.
## Sourcing this file only defines functions and has no side effects.

## Inputs: temporary path, destination path, and human-readable file label.
## Output: destination path invisibly after atomic replacement.
atomic_replace_file <- function(tmp, path, label) {
  if (isTRUE(suppressWarnings(file.rename(tmp, path)))) {
    return(invisible(path))
  }
  if (!file.exists(path)) {
    unlink(tmp)
    stop("Could not move temporary ", label, " to: ", path, call. = FALSE)
  }

  backup <- tempfile(
    pattern = paste0(basename(path), "_backup_"),
    tmpdir = dirname(path)
  )
  if (!isTRUE(suppressWarnings(file.rename(path, backup)))) {
    unlink(tmp)
    stop(
      "Could not preserve the existing ", label, " before replacement: ",
      path,
      call. = FALSE
    )
  }

  if (!isTRUE(suppressWarnings(file.rename(tmp, path)))) {
    restored <- isTRUE(suppressWarnings(file.rename(backup, path)))
    unlink(tmp)
    stop(
      "Could not move temporary ", label, " to: ", path,
      if (restored) {
        "; the previous file was restored."
      } else {
        paste0("; the previous file remains recoverable at ", backup, ".")
      },
      call. = FALSE
    )
  }

  unlink(backup)
  invisible(path)
}

## Inputs: R object and destination path.
## Output: destination path invisibly after atomically saving an RDS file.
atomic_save_rds <- function(object, path) {
  tmp <- tempfile(
    pattern = paste0(basename(path), "_"),
    tmpdir = dirname(path),
    fileext = ".rds"
  )
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3L)
  atomic_replace_file(tmp, path, "checkpoint")
}

## Inputs: table-like R object and destination path.
## Output: destination path invisibly after atomically writing a CSV file.
atomic_write_csv <- function(object, path) {
  tmp <- tempfile(
    pattern = paste0(basename(path), "_"),
    tmpdir = dirname(path),
    fileext = ".csv"
  )
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(
    object,
    tmp,
    row.names = FALSE,
    na = "NA",
    quote = TRUE
  )
  atomic_replace_file(tmp, path, "output")
}
