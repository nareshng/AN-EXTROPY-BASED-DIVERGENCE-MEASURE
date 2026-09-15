## =============================================================================
## Section 6.1 -- shared functions for the image-based analysis (Tables 12-14,
## Figure 4).
##
## Pixel intensities are treated as the two samples.  Image data have many ties
## (an 8-bit image takes at most 256 distinct values), so the estimator is
## evaluated exactly as the integral of the squared difference of the two
## empirical survival functions over the pooled distinct intensities; a
## rank formula with an arbitrary tie-breaking rule does not evaluate it.
##
## Requires: jpeg for .jpg/.jpeg input, png for .png input.
## =============================================================================

## Input: a numeric matrix or array read from an image file.
## Output: a grayscale matrix in [0, 1] (channel mean for colour images).
img_to_grayscale <- function(img) {
  if (!is.numeric(img)) stop("The image must be numeric.", call. = FALSE)
  if (length(dim(img)) == 3L) {
    channels <- dim(img)[3L]
    gray <- if (channels >= 3L) {
      (img[, , 1L] + img[, , 2L] + img[, , 3L]) / 3   # unweighted channel mean
    } else {
      img[, , 1L]
    }
  } else if (length(dim(img)) == 2L) {
    gray <- img
  } else {
    stop("Unsupported image dimensions.", call. = FALSE)
  }
  if (max(gray) > 1) gray <- gray / 255              # 8-bit input
  gray[] <- pmin(pmax(gray, 0), 1)
  gray
}

## Input: a file path. Output: list(values, matrix, width, height, n_distinct).
img_read <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  img <- switch(ext,
    jpg = , jpeg = {
      if (!requireNamespace("jpeg", quietly = TRUE)) {
        stop("Package 'jpeg' is required to read ", basename(path), ".", call. = FALSE)
      }
      jpeg::readJPEG(path)
    },
    png = {
      if (!requireNamespace("png", quietly = TRUE)) {
        stop("Package 'png' is required to read ", basename(path), ".", call. = FALSE)
      }
      png::readPNG(path)
    },
    stop("Unsupported image format: ", ext, call. = FALSE)
  )
  gray <- img_to_grayscale(img)
  values <- as.numeric(gray)
  if (!length(values) || anyNA(values)) stop("No usable pixels in ", basename(path), ".", call. = FALSE)
  list(values = values, matrix = gray, width = ncol(gray), height = nrow(gray),
       n_distinct = length(unique(values)))
}

## Input: two nonnegative samples.
## Output: the exact integral of {Fbar_n1(u) - Gbar_n2(u)}^2 du over [0, Inf).
## This is the V-statistic form of D and is tie-safe, symmetric and nonnegative.
img_divergence <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (!length(x) || !length(y) || anyNA(x) || anyNA(y)) {
    stop("Both samples must be nonempty and free of missing values.", call. = FALSE)
  }
  if (any(x < 0) || any(y < 0)) stop("Intensities must be nonnegative.", call. = FALSE)
  x <- sort(x); y <- sort(y)
  knots <- sort(unique(c(0, x, y)))
  if (length(knots) == 1L) return(0)
  left <- knots[-length(knots)]
  width <- diff(knots)
  Sx <- 1 - findInterval(left, x) / length(x)        # right-continuous survival
  Sy <- 1 - findInterval(left, y) / length(y)
  sum(width * (Sx - Sy)^2)
}

## Input: two nonnegative samples.
## Output: the off-diagonal estimator of equation (2.5), computed exactly.
## D_Emp = D_plugin - mean(x)/n1 - mean(y)/n2.
img_divergence_25 <- function(x, y) {
  img_divergence(x, y) - mean(x) / length(x) - mean(y) / length(y)
}

## Input: a named list of intensity vectors and the group order.
## Output: the symmetric matrix of pairwise divergences.
img_matrix <- function(values, ids, estimator = img_divergence) {
  k <- length(ids)
  out <- matrix(0, k, k, dimnames = list(ids, ids))
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k) {
    out[i, j] <- out[j, i] <- estimator(values[[ids[i]]], values[[ids[j]]])
  }
  out
}

## Input: a divergence matrix and its index. Output: one row per unordered pair.
img_long <- function(m, index) {
  ids <- rownames(m)
  pairs <- t(utils::combn(length(ids), 2L))
  data.frame(index = index, group_1 = ids[pairs[, 1L]], group_2 = ids[pairs[, 2L]],
             D_hat = m[pairs], stringsAsFactors = FALSE)
}

## Input: none. Output: TRUE invisibly; deterministic checks of the estimator.
img_self_test <- function() {
  checks <- logical(0)
  x <- c(0.1, 0.4, 0.4, 0.9); y <- c(0.2, 0.4, 0.7, 0.7)
  ## Closed form on a small example: integrate the step functions by hand.
  knots <- sort(unique(c(0, x, y)))
  left <- knots[-length(knots)]
  Sx <- 1 - findInterval(left, sort(x)) / 4
  Sy <- 1 - findInterval(left, sort(y)) / 4
  checks["matches explicit step integration"] <-
    isTRUE(all.equal(img_divergence(x, y), sum(diff(knots) * (Sx - Sy)^2)))
  checks["identical samples give zero"] <- isTRUE(all.equal(img_divergence(x, x), 0))
  checks["symmetric"] <- isTRUE(all.equal(img_divergence(x, y), img_divergence(y, x)))
  checks["nonnegative under heavy ties"] <- {
    set.seed(3)
    a <- sample(0:255, 5000, TRUE) / 255
    b <- sample(0:255, 5000, TRUE, prob = 1 + (0:255) / 255) / 255
    img_divergence(a, b) >= 0
  }
  ## Against the U-statistic representation D = E|min| terms, on a small sample.
  set.seed(7); a <- round(runif(60), 2); b <- round(runif(50), 2)
  v_form <- mean(outer(a, a, pmin)) + mean(outer(b, b, pmin)) - 2 * mean(outer(a, b, pmin))
  checks["equals the V-statistic form"] <- isTRUE(all.equal(img_divergence(a, b), v_form))
  checks["equation (2.5) differs by the diagonal terms"] <- isTRUE(all.equal(
    img_divergence_25(a, b), v_form - mean(a) / length(a) - mean(b) / length(b)))
  ## Scale equivariance: D(c x, c y) = c D(x, y).
  checks["scale equivariant"] <- isTRUE(all.equal(img_divergence(2 * a, 2 * b),
                                                  2 * img_divergence(a, b)))
  ## Grayscale conversion: 8-bit input is rescaled, colour is averaged.
  arr <- array(c(matrix(255, 2, 2), matrix(0, 2, 2), matrix(255, 2, 2)), dim = c(2, 2, 3))
  checks["grayscale conversion"] <- isTRUE(all.equal(as.numeric(img_to_grayscale(arr)),
                                                     rep(2 / 3, 4)))
  message("Image estimator self-checks: ", sum(checks), "/", length(checks), " passed.")
  if (!all(checks)) message("Failed: ", paste(names(checks)[!checks], collapse = ", "))
  invisible(all(checks))
}
