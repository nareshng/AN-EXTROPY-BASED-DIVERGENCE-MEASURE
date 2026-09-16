#!/usr/bin/env Rscript
## =============================================================================
## Section 6.1 -- Tables 12-14 and Figure 4: divergence between MRI images.
##
## For each image index, the script forms the three pairwise divergences among
## the NT, BT and MT images. Images are converted to grayscale in [0, 1] and,
## by default, resized to the explicit common size of 227 x 255 pixels
## (width x height) using deterministic bilinear interpolation.
##
## The primary default is the empirical estimator in equation (2.5). The
## nonnegative plug-in integral is also saved as a labelled diagnostic; the two
## quantities are never silently interchanged.
## =============================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
analysis_script_path <- if (length(script_arg)) {
  normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE),
                winslash = "/", mustWork = TRUE)
} else {
  NA_character_
}
script_dir <- if (!is.na(analysis_script_path)) {
  dirname(analysis_script_path)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

function_candidates <- file.path(script_dir, c("image_functions.R", "image_functions(2).R"))
function_path <- function_candidates[file.exists(function_candidates)][1L]
if (is.na(function_path)) {
  stop("Could not find image_functions.R beside this script.", call. = FALSE)
}
source(function_path)

usage <- c(
  "Usage:",
  "  Rscript Image_based_Real_Data_Analysis.R [options]",
  "  Rscript Image_based_Real_Data_Analysis.R --self-test",
  "",
  "Options:",
  "  --image-dir=PATH       folder containing <group><index>.<ext>",
  "  --indices=LIST         comma-separated indices (default: complete triples)",
  "  --extension=EXT        jpg, jpeg or png (default: jpg)",
  "  --estimator=NAME       empirical (equation 2.5) or plugin (default: empirical)",
  "  --resize=SPEC          WIDTHxHEIGHT or none (default: 227x255)",
  "  --resize-method=NAME   bilinear or nearest (default: bilinear)",
  "  --manifest=PATH        optional CSV with one row per input file and a file column",
  "  --digits=N             digits printed to the console (default: 3)",
  "  --output-dir=PATH      output folder (default: image_outputs beside this script)",
  "  --no-figure            skip Figure 4",
  "  --self-test            run deterministic checks and exit",
  "  --help"
)

GROUPS <- c("NT", "BT", "MT")
default_image_dir <- local({
  candidates <- file.path(script_dir, c("Images", "images"))
  hit <- candidates[dir.exists(candidates)]
  if (length(hit)) hit[1L] else candidates[1L]
})
CFG <- list(
  image_dir = default_image_dir,
  indices = NULL,
  extension = "jpg",
  estimator = "empirical",
  resize = "227x255",
  resize_method = "bilinear",
  manifest = NULL,
  digits = 3L,
  output_dir = file.path(script_dir, "image_outputs"),
  figure = TRUE,
  self_test = FALSE
)

for (arg in commandArgs(trailingOnly = TRUE)) {
  if (identical(arg, "--help")) {
    cat(paste(usage, collapse = "\n"), "\n")
    quit(save = "no", status = 0L)
  } else if (identical(arg, "--self-test")) {
    CFG$self_test <- TRUE
  } else if (identical(arg, "--no-figure")) {
    CFG$figure <- FALSE
  } else if (startsWith(arg, "--image-dir=")) {
    CFG$image_dir <- sub("^--image-dir=", "", arg)
  } else if (startsWith(arg, "--indices=")) {
    CFG$indices <- as.integer(strsplit(sub("^--indices=", "", arg), ",", fixed = TRUE)[[1L]])
  } else if (startsWith(arg, "--extension=")) {
    CFG$extension <- tolower(sub("^--extension=", "", arg))
  } else if (startsWith(arg, "--estimator=")) {
    CFG$estimator <- tolower(sub("^--estimator=", "", arg))
  } else if (startsWith(arg, "--resize=")) {
    CFG$resize <- tolower(sub("^--resize=", "", arg))
  } else if (startsWith(arg, "--resize-method=")) {
    CFG$resize_method <- tolower(sub("^--resize-method=", "", arg))
  } else if (startsWith(arg, "--manifest=")) {
    CFG$manifest <- sub("^--manifest=", "", arg)
  } else if (startsWith(arg, "--digits=")) {
    CFG$digits <- as.integer(sub("^--digits=", "", arg))
  } else if (startsWith(arg, "--output-dir=")) {
    CFG$output_dir <- sub("^--output-dir=", "", arg)
  } else {
    stop("Unknown option: ", arg, "\n\n", paste(usage, collapse = "\n"), call. = FALSE)
  }
}

if (CFG$self_test) {
  ok <- img_self_test()
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}
if (!CFG$extension %in% c("jpg", "jpeg", "png")) {
  stop("--extension must be jpg, jpeg or png.", call. = FALSE)
}
if (!CFG$estimator %in% c("empirical", "plugin")) {
  stop("--estimator must be empirical or plugin.", call. = FALSE)
}
if (!CFG$resize_method %in% c("bilinear", "nearest")) {
  stop("--resize-method must be bilinear or nearest.", call. = FALSE)
}
if (length(CFG$digits) != 1L || is.na(CFG$digits) || CFG$digits < 0L) {
  stop("--digits must be a nonnegative integer.", call. = FALSE)
}
if (!dir.exists(CFG$image_dir)) {
  stop("Image directory not found: ", CFG$image_dir,
       "\nLooked for Images/ and images/ beside the script; ",
       "give the location with --image-dir=PATH.", call. = FALSE)
}

target_size <- NULL
if (!identical(CFG$resize, "none")) {
  resize_match <- regexec("^([1-9][0-9]*)x([1-9][0-9]*)$", CFG$resize)
  resize_parts <- regmatches(CFG$resize, resize_match)[[1L]]
  if (length(resize_parts) != 3L) {
    stop("--resize must be WIDTHxHEIGHT (for example, 227x255) or none.", call. = FALSE)
  }
  target_size <- c(width = as.integer(resize_parts[2L]), height = as.integer(resize_parts[3L]))
}

path_of <- function(group, index) {
  file.path(CFG$image_dir, paste0(group, index, ".", CFG$extension))
}
if (is.null(CFG$indices)) {
  pattern <- paste0("^(NT|BT|MT)[0-9]+\\.", CFG$extension, "$")
  files <- list.files(CFG$image_dir, pattern = pattern, ignore.case = TRUE)
  found <- as.integer(unique(sub("^(NT|BT|MT)([0-9]+)\\..*$", "\\2", toupper(files))))
  CFG$indices <- sort(found[vapply(found, function(i) {
    all(file.exists(vapply(GROUPS, path_of, character(1), index = i)))
  }, logical(1))])
} else {
  if (anyNA(CFG$indices) || any(CFG$indices < 1L)) {
    stop("--indices must contain positive integers.", call. = FALSE)
  }
  CFG$indices <- unique(CFG$indices)
}
if (!length(CFG$indices)) {
  stop("No complete NT/BT/MT triple was found in ", CFG$image_dir, call. = FALSE)
}
missing <- unlist(lapply(CFG$indices, function(i) {
  paths <- vapply(GROUPS, path_of, character(1), index = i)
  paths[!file.exists(paths)]
}))
if (length(missing)) {
  stop("Missing image file(s): ", paste(basename(missing), collapse = ", "), call. = FALSE)
}

out_dir <- CFG$output_dir
if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create the output directory: ", out_dir, call. = FALSE)
}

## ---- read and preprocess once ------------------------------------------------
images <- list()
values <- list()
meta <- list()
for (g in GROUPS) for (i in CFG$indices) {
  id <- paste0(g, i)
  path <- path_of(g, i)
  raw <- img_read(path)
  analysed <- if (is.null(target_size)) {
    raw$matrix
  } else {
    img_resize(raw$matrix, target_size[["width"]], target_size[["height"]],
               method = CFG$resize_method)
  }
  images[[id]] <- analysed
  values[[id]] <- as.numeric(analysed)
  info <- file.info(path)
  meta[[length(meta) + 1L]] <- data.frame(
    image = id,
    group = g,
    index = i,
    file = basename(path),
    md5 = unname(tools::md5sum(path)),
    file_size_bytes = unname(info$size),
    original_width = raw$width,
    original_height = raw$height,
    analysed_width = ncol(analysed),
    analysed_height = nrow(analysed),
    resized = raw$width != ncol(analysed) || raw$height != nrow(analysed),
    original_pixels = length(raw$values),
    analysed_pixels = length(values[[id]]),
    original_distinct_intensities = raw$n_distinct,
    analysed_distinct_intensities = length(unique(values[[id]])),
    mean_intensity = mean(values[[id]]),
    stringsAsFactors = FALSE
  )
}
meta <- do.call(rbind, meta)

## Attach optional provenance fields without changing the input order.
if (is.null(CFG$manifest)) {
  candidate_manifest <- file.path(CFG$image_dir, "image_manifest.csv")
  if (file.exists(candidate_manifest)) CFG$manifest <- candidate_manifest
}
if (!is.null(CFG$manifest)) {
  if (!file.exists(CFG$manifest)) stop("Manifest not found: ", CFG$manifest, call. = FALSE)
  manifest <- utils::read.csv(CFG$manifest, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"file" %in% names(manifest)) stop("The manifest must contain a 'file' column.", call. = FALSE)
  if (anyDuplicated(manifest$file)) stop("The manifest contains duplicate file names.", call. = FALSE)
  row_id <- match(meta$file, manifest$file)
  if (anyNA(row_id)) {
    stop("Manifest entries are missing for: ", paste(meta$file[is.na(row_id)], collapse = ", "),
         call. = FALSE)
  }
  manifest_rows <- manifest[row_id, setdiff(names(manifest), "file"), drop = FALSE]
  if ("analysis_group" %in% names(manifest_rows)) {
    supplied_group <- as.character(manifest_rows$analysis_group)
    if (anyNA(supplied_group) || any(supplied_group != meta$group)) {
      stop("Manifest analysis_group values do not agree with the filenames.", call. = FALSE)
    }
  }
  if ("analysis_index" %in% names(manifest_rows)) {
    supplied_index <- suppressWarnings(as.integer(manifest_rows$analysis_index))
    if (anyNA(supplied_index) || any(supplied_index != meta$index)) {
      stop("Manifest analysis_index values do not agree with the filenames.", call. = FALSE)
    }
  }
  meta <- cbind(meta, manifest_rows)
} else {
  warning("No image_manifest.csv was supplied; byte-level checksums are recorded, but dataset label provenance cannot be verified.",
          call. = FALSE)
}

## ---- estimate and write Tables 12-14 -----------------------------------------
long <- list()
for (i in CFG$indices) {
  ids <- paste0(GROUPS, i)
  matrix_empirical <- img_matrix(values, ids, img_divergence_empirical)
  matrix_plugin <- img_matrix(values, ids, img_divergence_plugin)
  matrix_primary <- if (CFG$estimator == "empirical") matrix_empirical else matrix_plugin

  utils::write.csv(round(matrix_primary, 6),
                   file.path(out_dir, sprintf("Divergence_matrix_index%d.csv", i)))
  utils::write.csv(round(matrix_plugin, 6),
                   file.path(out_dir, sprintf("Divergence_matrix_index%d_plugin.csv", i)))

  rows <- img_long(matrix_primary, i)
  rows$estimator <- CFG$estimator
  rows$D_hat_Emp <- img_long(matrix_empirical, i)$D_hat
  rows$D_hat_plugin <- img_long(matrix_plugin, i)$D_hat
  rows <- rows[, c("index", "group_1", "group_2", "estimator", "D_hat",
                   "D_hat_Emp", "D_hat_plugin")]
  long[[length(long) + 1L]] <- rows

  if (any(matrix_primary[upper.tri(matrix_primary)] < 0)) {
    message("Index ", i, ": equation (2.5) produced a negative finite-sample estimate.")
  }
  cat(sprintf("\nIndex %d (images %s; primary estimator: %s):\n",
              i, paste(ids, collapse = ", "), CFG$estimator))
  print(round(matrix_primary, CFG$digits))
}
long <- do.call(rbind, long)
utils::write.csv(long, file.path(out_dir, "Image_divergences_long.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(out_dir, "Image_inventory.csv"), row.names = FALSE)

## ---- Figure 4: rows are NT/BT/MT; columns are image indices ------------------
if (CFG$figure) {
  grDevices::pdf(file.path(out_dir, "Figure4_MRI_images.pdf"),
                 width = 3 * length(CFG$indices), height = 3 * length(GROUPS))
  graphics::par(mfrow = c(length(GROUPS), length(CFG$indices)), mar = c(1, 1, 2.5, 1))
  for (g in GROUPS) for (i in CFG$indices) {
    id <- paste0(g, i)
    matrix_image <- images[[id]]
    graphics::image(t(matrix_image[nrow(matrix_image):1L, ]),
                    col = grDevices::gray.colors(256, 0, 1), axes = FALSE,
                    main = id, useRaster = TRUE, asp = 1)
    graphics::box()
  }
  grDevices::dev.off()
}

estimator_description <- if (CFG$estimator == "empirical") {
  "equation (2.5): plug-in integral minus mean(x)/n1 minus mean(y)/n2"
} else {
  "exact plug-in integral of (Fbar_n1 - Gbar_n2)^2"
}
resize_description <- if (is.null(target_size)) {
  "none; native image dimensions retained"
} else {
  paste0(target_size[["width"]], "x", target_size[["height"]], " (width x height)")
}
manifest_description <- if (is.null(CFG$manifest)) "not supplied" else normalizePath(CFG$manifest)
analysis_md5 <- if (is.na(analysis_script_path)) "not available" else unname(tools::md5sum(analysis_script_path))
settings <- c(
  paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Image directory:", normalizePath(CFG$image_dir)),
  paste("Indices analysed:", paste(CFG$indices, collapse = ", ")),
  paste("Groups:", paste(GROUPS, collapse = ", ")),
  paste("Primary estimator:", estimator_description),
  paste("Resize:", resize_description),
  paste("Resize method:", if (is.null(target_size)) "not applicable" else CFG$resize_method),
  "Grayscale conversion: unweighted mean of RGB channels; intensities scaled to [0, 1]",
  paste("Provenance manifest:", manifest_description),
  paste("Analysis script MD5:", analysis_md5),
  paste("Function script MD5:", unname(tools::md5sum(function_path))),
  paste("R version:", R.version.string)
)
writeLines(settings, file.path(out_dir, "Section6_1_run_settings.txt"))
utils::capture.output(utils::sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("\nOutputs written to ", normalizePath(out_dir), "\n", sep = "")
