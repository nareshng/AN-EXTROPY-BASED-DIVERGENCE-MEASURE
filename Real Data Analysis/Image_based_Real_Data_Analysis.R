#!/usr/bin/env Rscript
## =============================================================================
## Section 6.1 -- Tables 12-14 and Figure 4: divergence between MRI images.
##
## For every image index the script forms the three pairwise divergences between
## the normal (NT), benign (BT) and malignant (MT) images and writes a 3 x 3
## matrix.  Pixel intensities are converted to grayscale in [0, 1] and the
## estimator is evaluated exactly on the pooled distinct intensities, which is
## the tie-safe form needed for 8-bit image data.
##
## Usage (run from this folder):
##   Rscript Image_based_Real_Data_Analysis.R
##   Rscript Image_based_Real_Data_Analysis.R --image-dir=Images --indices=1,2,4
##   Rscript Image_based_Real_Data_Analysis.R --self-test
##
## Options
##   --image-dir=PATH     directory holding <group><index>.<ext> (default Images/)
##   --indices=LIST       image indices to analyse (default: every complete triple)
##   --extension=EXT      jpg, jpeg or png (default jpg)
##   --digits=N           digits in the printed matrices (default 3)
##   --output-dir=PATH    output directory (default image_outputs)
##   --no-figure          skip Figure 4
##   --self-test          run the deterministic estimator checks and exit
##   --help
##
## Requires: jpeg (for .jpg) or png (for .png).  No other package is needed.
## =============================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE),
                        winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(script_dir, "image_functions.R"))

GROUPS <- c("NT", "BT", "MT")
CFG <- list(image_dir = file.path(script_dir, "Images"), indices = NULL, extension = "jpg",
            digits = 3L, output_dir = "image_outputs", figure = TRUE, self_test = FALSE)
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (identical(arg, "--help")) {
    cat(paste(readLines(file.path(script_dir, "Image_based_Real_Data_Analysis.R"), n = 28)[-1L],
              collapse = "\n"), "\n")
    quit(save = "no", status = 0L)
  }
  else if (identical(arg, "--self-test")) CFG$self_test <- TRUE
  else if (identical(arg, "--no-figure")) CFG$figure <- FALSE
  else if (startsWith(arg, "--image-dir=")) CFG$image_dir <- sub("^--image-dir=", "", arg)
  else if (startsWith(arg, "--indices=")) CFG$indices <- as.integer(strsplit(sub("^--indices=", "", arg), ",")[[1L]])
  else if (startsWith(arg, "--extension=")) CFG$extension <- tolower(sub("^--extension=", "", arg))
  else if (startsWith(arg, "--digits=")) CFG$digits <- as.integer(sub("^--digits=", "", arg))
  else if (startsWith(arg, "--output-dir=")) CFG$output_dir <- sub("^--output-dir=", "", arg)
  else stop("Unknown option: ", arg, call. = FALSE)
}
if (CFG$self_test) {
  ok <- img_self_test()
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}
if (!CFG$extension %in% c("jpg", "jpeg", "png")) {
  stop("--extension must be jpg, jpeg or png.", call. = FALSE)
}
if (!dir.exists(CFG$image_dir)) {
  stop("Image directory not found: ", CFG$image_dir,
       "\nGive its location with --image-dir=PATH.", call. = FALSE)
}

## ---- which indices are available ---------------------------------------------
path_of <- function(group, index) {
  file.path(CFG$image_dir, paste0(group, index, ".", CFG$extension))
}
if (is.null(CFG$indices)) {
  files <- list.files(CFG$image_dir, pattern = paste0("^(NT|BT|MT)[0-9]+\\.", CFG$extension, "$"))
  found <- as.integer(unique(sub("^(NT|BT|MT)([0-9]+)\\..*$", "\\2", files)))
  CFG$indices <- sort(found[vapply(found, function(i)
    all(file.exists(vapply(GROUPS, path_of, character(1), index = i))), logical(1))])
}
if (!length(CFG$indices)) {
  stop("No complete NT/BT/MT triple was found in ", CFG$image_dir, call. = FALSE)
}
missing <- unlist(lapply(CFG$indices, function(i) {
  p <- vapply(GROUPS, path_of, character(1), index = i); p[!file.exists(p)]
}))
if (length(missing)) stop("Missing image file(s): ", paste(basename(missing), collapse = ", "), call. = FALSE)

out_dir <- CFG$output_dir
if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create the output directory: ", out_dir, call. = FALSE)
}

## ---- read, estimate, write ----------------------------------------------------
images <- list(); meta <- list(); long <- list()
for (i in CFG$indices) {
  values <- list()
  for (g in GROUPS) {
    id <- paste0(g, i)
    im <- img_read(path_of(g, i))
    values[[id]] <- im$values
    images[[id]] <- im$matrix
    meta[[length(meta) + 1L]] <- data.frame(
      image = id, file = basename(path_of(g, i)), width = im$width, height = im$height,
      pixels = length(im$values), distinct_intensities = im$n_distinct,
      mean_intensity = mean(im$values), stringsAsFactors = FALSE)
  }
  ids <- paste0(GROUPS, i)
  m_plugin <- img_matrix(values, ids, img_divergence)
  m_25 <- img_matrix(values, ids, img_divergence_25)
  utils::write.csv(round(m_plugin, 6), file.path(out_dir, sprintf("Divergence_matrix_index%d.csv", i)))
  rows <- img_long(m_plugin, i)
  rows$D_hat_eq2_5 <- img_long(m_25, i)$D_hat
  long[[length(long) + 1L]] <- rows
  cat(sprintf("\nIndex %d (images %s):\n", i, paste(ids, collapse = ", ")))
  print(round(m_plugin, CFG$digits))
}
long <- do.call(rbind, long)
meta <- do.call(rbind, meta)
utils::write.csv(long, file.path(out_dir, "Image_divergences_long.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(out_dir, "Image_inventory.csv"), row.names = FALSE)

## ---- Figure 4: the image panel ------------------------------------------------
if (CFG$figure) {
  grDevices::pdf(file.path(out_dir, "Figure4_MRI_images.pdf"), width = 9,
                 height = 3 * length(CFG$indices))
  graphics::par(mfrow = c(length(CFG$indices), 3L), mar = c(1, 1, 2.5, 1))
  for (i in CFG$indices) for (g in GROUPS) {
    m <- images[[paste0(g, i)]]
    graphics::image(t(m[nrow(m):1L, ]), col = grDevices::gray.colors(256, 0, 1),
                    axes = FALSE, main = paste0(g, i), useRaster = TRUE)
    graphics::box()
  }
  grDevices::dev.off()
}

writeLines(c(
  paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Image directory:", normalizePath(CFG$image_dir)),
  paste("Indices analysed:", paste(CFG$indices, collapse = ", ")),
  paste("Estimator: exact integral of (Fbar - Gbar)^2 on the pooled distinct intensities"),
  paste("R version:", R.version.string)
), file.path(out_dir, "Section6_1_run_settings.txt"))
utils::capture.output(utils::sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("\nOutputs written to ", normalizePath(out_dir), "\n", sep = "")
