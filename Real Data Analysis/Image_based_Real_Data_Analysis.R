# ================================================================
# IMAGE-BASED EMPIRICAL DIVERGENCE TABLES
#
# Output:
#   image_divergence_final_tables.csv
#
# Images expected:
#   NT1.jpg, BT1.jpg, MT1.jpg
#   NT2.jpg, BT2.jpg, MT2.jpg
#   NT4.jpg, BT4.jpg, MT4.jpg
#   NT5.jpg, BT5.jpg, MT5.jpg
# ================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(png)
  library(jpeg)
})

# ================================================================
# 1. CONFIGURATION
# ================================================================

CFG <- list(
  image_dir = "/Users/gargn2/Downloads",
  out_csv   = "image_divergence_final_tables.csv",
  digits    = 3
)

GROUPS  <- c("NT", "BT", "MT")
INDICES <- c(1, 2, 4, 5)

# ================================================================
# 2. IMAGE LOADING
# ================================================================

load_grayscale_image <- function(path) {
  
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  
  ext <- tolower(tools::file_ext(path))
  
  img <- switch(
    ext,
    "png"  = readPNG(path),
    "jpg"  = readJPEG(path),
    "jpeg" = readJPEG(path),
    stop("Unsupported image format: ", ext)
  )
  
  # Convert RGB/RGBA image to grayscale
  if (length(dim(img)) == 3) {
    img <- apply(img[, , 1:3, drop = FALSE], c(1, 2), mean)
  }
  
  x <- as.vector(img)
  x <- x[is.finite(x)]
  
  # Ensure pixel intensities are in [0, 1]
  x <- pmin(pmax(x, 0), 1)
  
  return(x)
}

make_image_path <- function(group, index, image_dir) {
  file.path(image_dir, paste0(group, index, ".jpg"))
}

load_image_set <- function(groups, indices, image_dir) {
  
  images <- list()
  
  for (idx in indices) {
    for (grp in groups) {
      
      id   <- paste0(grp, idx)
      path <- make_image_path(grp, idx, image_dir)
      
      images[[id]] <- load_grayscale_image(path)
    }
  }
  
  return(images)
}

# ================================================================
# 3. EMPIRICAL ESTIMATOR OF D
#
# Based on equation:
# D_hat_Emp =
#   (2/n1) sum_j X_(j) {S_j/n2 - j(1/n1 + 1/n2)}
# + (2/n2) sum_j Y_(j) {R_j/n1 - j(1/n1 + 1/n2)}
#
# Ties are common in image data, so ties.method = "max" is used
# to match the right-continuous empirical CDF convention.
# ================================================================

compute_D_emp <- function(x, y) {
  
  x <- sort(x[is.finite(x)])
  y <- sort(y[is.finite(y)])
  
  n1 <- length(x)
  n2 <- length(y)
  
  if (n1 < 2 || n2 < 2) {
    stop("Both samples must contain at least two observations.")
  }
  
  pooled <- c(x, y)
  ranks  <- rank(pooled, ties.method = "max")
  
  S <- ranks[seq_len(n1)]
  R <- ranks[(n1 + 1):(n1 + n2)]
  
  # Reorder ranks according to sorted x and sorted y
  S <- S[order(pooled[seq_len(n1)])]
  R <- R[order(pooled[(n1 + 1):(n1 + n2)])]
  
  j1 <- seq_len(n1)
  j2 <- seq_len(n2)
  
  term_x <- (2 / n1) * sum(
    x * (S / n2 - j1 * (1 / n1 + 1 / n2))
  )
  
  term_y <- (2 / n2) * sum(
    y * (R / n1 - j2 * (1 / n1 + 1 / n2))
  )
  
  D <- term_x + term_y
  
  # Numerical/tie effects can create tiny negative values
  max(D, 0)
}

# ================================================================
# 4. CREATE 3 x 3 DIVERGENCE MATRIX FOR ONE INDEX
# ================================================================

make_divergence_matrix <- function(images, index, groups = GROUPS, digits = 3) {
  
  ids <- paste0(groups, index)
  
  mat <- matrix(
    0,
    nrow = length(ids),
    ncol = length(ids),
    dimnames = list(ids, ids)
  )
  
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      
      if (i == j) {
        mat[i, j] <- 0
      } else {
        mat[i, j] <- compute_D_emp(images[[ids[i]]], images[[ids[j]]])
      }
    }
  }
  
  round(mat, digits)
}

# ================================================================
# 5. CONVERT MATRIX TO LONG CSV FORMAT
# ================================================================

matrix_to_long <- function(mat, index) {
  
  out <- expand.grid(
    RowImage = rownames(mat),
    ColImage = colnames(mat),
    stringsAsFactors = FALSE
  )
  
  out$Index      <- index
  out$Comparison <- paste0(out$RowImage, "_vs_", out$ColImage)
  out$D_emp      <- as.vector(mat)
  
  out <- out[, c("Index", "RowImage", "ColImage", "Comparison", "D_emp")]
  
  return(out)
}

# ================================================================
# 6. RUN ANALYSIS
# ================================================================

cat("\n=== Loading images ===\n")

images <- load_image_set(
  groups    = GROUPS,
  indices   = INDICES,
  image_dir = CFG$image_dir
)

cat("Loaded images:\n")
print(names(images))

cat("\n=== Computing divergence matrices ===\n")

all_results <- list()

for (idx in INDICES) {
  
  mat <- make_divergence_matrix(
    images = images,
    index  = idx,
    groups = GROUPS,
    digits = CFG$digits
  )
  
  cat("\nIndex:", idx, "\n")
  print(mat)
  
  all_results[[as.character(idx)]] <- matrix_to_long(mat, idx)
}

final_results <- do.call(rbind, all_results)
rownames(final_results) <- NULL

write.csv(final_results, CFG$out_csv, row.names = FALSE)

cat("\nSaved final CSV:", CFG$out_csv, "\n")
cat("\n=== DONE ===\n")
