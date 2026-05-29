# ================================================================
# IMAGE-BASED DIVERGENCE ANALYSIS
#
# Measures:
#   1. D_emp  : empirical estimator of D
#   2. D_grid : grid-based integral of survival-curve difference
#   3. D_CC   : Cox--Czanner-type divergence
#   4. KL     : kernel plug-in Kullback--Leibler divergence
#
# Data:
#   Grayscale image pixel intensities scaled to [0, 1]
# ================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(png)
  library(jpeg)
  library(pracma)
})

# ================================================================
# 1. CONFIGURATION
# ================================================================

CFG <- list(
  image_dir   = "/Users/gargn2/Downloads",
  grid_points = 512,
  eps         = 1e-10
)

# Image files to compare
IMAGE_FILES <- c(
  NT5 = "NT5.jpg",
  BT5 = "BT5.jpg",
  MT5 = "MT5.jpg",
  V1  = "T1.jpg",
  V2  = "T2.jpg",
  V3  = "T3.jpg",
  V4  = "T4.jpg",
  V5  = "T5.jpeg",
  V6  = "V3.jpg"
)

# Pairwise comparisons
COMPARISONS <- list(
  c("V1", "V1"),
  c("V1", "V2"),
  c("V1", "V3"),
  c("V1", "V4"),
  c("V2", "V3")
)

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
  
  # Ensure intensities lie in [0, 1]
  x <- pmin(pmax(x, 0), 1)
  
  return(x)
}

load_image_set <- function(image_dir, image_files) {
  
  out <- lapply(image_files, function(fname) {
    load_grayscale_image(file.path(image_dir, fname))
  })
  
  names(out) <- names(image_files)
  return(out)
}

# ================================================================
# 3. BASIC ESTIMATION UTILITIES
# ================================================================

empirical_survival <- function(x, grid) {
  vapply(grid, function(t) mean(x > t), numeric(1))
}

kernel_density_on_grid <- function(x, grid, eps = 1e-10) {
  
  dens <- density(
    x,
    from = min(grid),
    to   = max(grid),
    n    = length(grid),
    na.rm = TRUE
  )
  
  fx <- approx(dens$x, dens$y, xout = grid, rule = 2)$y
  pmax(fx, eps)
}

# ================================================================
# 4. DIVERGENCE ESTIMATORS
# ================================================================

# ------------------------------------------------
# 4.1 Grid-based estimator of D
#     D = int_0^1 {S1(x) - S2(x)}^2 dx
# ------------------------------------------------

compute_D_grid <- function(x, y, grid_points = 512) {
  
  grid <- seq(0, 1, length.out = grid_points)
  
  Sx <- empirical_survival(x, grid)
  Sy <- empirical_survival(y, grid)
  
  trapz(grid, (Sx - Sy)^2)
}

# ------------------------------------------------
# 4.2 Weighted version Dw
#     Dw = int_0^1 w(x){S1(x) - S2(x)}^2 dx
# ------------------------------------------------

compute_Dw_grid <- function(x, y, weight_fun, grid_points = 512) {
  
  grid <- seq(0, 1, length.out = grid_points)
  
  Sx <- empirical_survival(x, grid)
  Sy <- empirical_survival(y, grid)
  wx <- weight_fun(grid)
  
  if (length(wx) != length(grid)) {
    stop("weight_fun must return a vector of the same length as grid.")
  }
  
  trapz(grid, wx * (Sx - Sy)^2)
}

# ------------------------------------------------
# 4.3 Empirical estimator of D from equation (2.5)
# ------------------------------------------------

compute_D_emp <- function(x, y) {
  
  x <- sort(x[is.finite(x)])
  y <- sort(y[is.finite(y)])
  
  n1 <- length(x)
  n2 <- length(y)
  
  if (n1 < 2 || n2 < 2) {
    stop("Both samples must contain at least two observations.")
  }
  
  # Ranks in pooled sample.
  # For continuous data ties are unlikely. For pixel data ties are common,
  # so ties.method = 'max' gives the right-continuous empirical CDF convention.
  pooled <- c(x, y)
  
  S <- rank(pooled, ties.method = "max")[seq_len(n1)]
  R <- rank(pooled, ties.method = "max")[(n1 + 1):(n1 + n2)]
  
  # Sort ranks according to sorted x and sorted y
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
  
  D_emp <- term_x + term_y
  
  # Numerical/tie effects can produce tiny negative values
  max(D_emp, 0)
}

# ------------------------------------------------
# 4.4 Cox--Czanner-type divergence
#     D_CC = int |S1(x) f2(x) - S2(x) f1(x)| dx
# ------------------------------------------------

compute_DCC <- function(x, y, grid_points = 512, eps = 1e-10) {
  
  grid <- seq(0, 1, length.out = grid_points)
  
  Sx <- empirical_survival(x, grid)
  Sy <- empirical_survival(y, grid)
  
  fx <- kernel_density_on_grid(x, grid, eps)
  fy <- kernel_density_on_grid(y, grid, eps)
  
  trapz(grid, abs(Sx * fy - Sy * fx))
}

# ------------------------------------------------
# 4.5 KL divergence
#     KL(f || g) = int f log(f/g)
# ------------------------------------------------

compute_KL <- function(x, y, grid_points = 512, eps = 1e-10) {
  
  grid <- seq(0, 1, length.out = grid_points)
  
  fx <- kernel_density_on_grid(x, grid, eps)
  fy <- kernel_density_on_grid(y, grid, eps)
  
  trapz(grid, fx * log(fx / fy))
}

# ------------------------------------------------
# 4.6 Symmetric KL / Jeffreys divergence
# ------------------------------------------------

compute_KL_sym <- function(x, y, grid_points = 512, eps = 1e-10) {
  compute_KL(x, y, grid_points, eps) +
    compute_KL(y, x, grid_points, eps)
}

# ================================================================
# 5. PAIRWISE COMPARISON WRAPPER
# ================================================================

compare_pair <- function(images, id1, id2, cfg = CFG) {
  
  x <- images[[id1]]
  y <- images[[id2]]
  
  data.frame(
    Image1 = id1,
    Image2 = id2,
    D_emp  = compute_D_emp(x, y),
    D_grid = compute_D_grid(x, y, cfg$grid_points),
    D_CC   = compute_DCC(x, y, cfg$grid_points, cfg$eps),
    KL     = compute_KL(x, y, cfg$grid_points, cfg$eps),
    KL_sym = compute_KL_sym(x, y, cfg$grid_points, cfg$eps),
    stringsAsFactors = FALSE
  )
}

run_comparisons <- function(images, comparisons, cfg = CFG) {
  
  out <- lapply(comparisons, function(pair) {
    compare_pair(images, pair[1], pair[2], cfg)
  })
  
  do.call(rbind, out)
}

# ================================================================
# 6. RUN ANALYSIS
# ================================================================

cat("\n=== Loading images ===\n")

images <- load_image_set(CFG$image_dir, IMAGE_FILES)

cat("Loaded images:\n")
print(names(images))

cat("\n=== Running pairwise comparisons ===\n")

results <- run_comparisons(images, COMPARISONS, CFG)

print(results)

# Optional: save results
write.csv(results, "image_divergence_results.csv", row.names = FALSE)

cat("\nSaved: image_divergence_results.csv\n")
cat("\n=== DONE ===\n")