# This function take both 1D (numeric vector) and multidimensional (matrix) data;
# The algorithm chooses a random data point in "treat" and picks up the most similar data point from "ctrl" (the one with the minimal Euclidean distance);
# A "ctrl" data point cannot be chosen twice (sampling w/o replacement);
# A data point in "treat" can be selected twice, if the query size exceeds the size of "treat"; 
# The input data are expected to be transformed to normal (sqrt, log, inverse...) outside of the function;
# The normalized columns are standardized (zero mean, unit variance) inside of the function;


standardize_matrix_cols <- function(mat, col_mean, col_sd) {
  stopifnot("matrix" %in% class(mat))
  stopifnot(length(col_mean) == ncol(mat))
  stopifnot(length(col_sd) == ncol(mat))
  if (ncol(mat) == 1) {
    scaled <- (mat - col_mean) / col_sd
  } else {
    shifted <- t(apply(mat, 1, `-`, col_mean))
    scaled <- t(apply(shifted, 1, `/`, col_sd))
  }
  return(scaled)
}

find_uniq_min_with_mask <- function(vec, mask) {
  # mask is expected to contain TRUE for the values to use
  min <- min(vec[mask])
  idx <- which(vec == min & mask)
  if (length(idx) > 1) {
    idx <- sample(idx, size = 1)
  }
  return(idx)
}

calculate_distance <- function(x_vec, y_mat) {
  if (ncol(y_mat) == 1) {
    dist <- abs(as.numeric(y_mat) - x_vec)
  } else {
    dist <- apply(y_mat, 1, `-`, x_vec) %>% t() %>% `^`(2) %>% rowSums() %>% sqrt()
  }
  return(dist)
}

find_matched_control <- function(treat, ctrl, size = NULL, replace = FALSE) {
  # Force input to matrix (in case of 1D data):
  x <- as.matrix(treat)
  y <- as.matrix(ctrl)
  # Sanity checks:
  stopifnot(ncol(x) == ncol(y))
  n1 <- nrow(x)
  n2 <- nrow(y)
  if (is.null(size)) {
    size <- n1
  }
  stopifnot(size < n2)
  if (size > n1) {
    replace <- TRUE
  }
  # Combine the treat and control data points:
  xy <- rbind(x, y)
  idx_x <- seq(1, n1)
  idx_y <- seq(1, n2) + n1
  # Standardize:
  col_mean <- colMeans(xy)
  col_sd <- apply(xy, 2, sd)
  xy_norm <- standardize_matrix_cols(xy, col_mean, col_sd)
  x_norm <- xy_norm[idx_x, ] %>% as.matrix()
  y_norm <- xy_norm[idx_y, ] %>% as.matrix()
  # Define the order of queries in x:
  idx_query <- sample(idx_x, size = size, replace = replace)
  # Sample from the control dataset:
  idx_out <- vector("integer", size)
  unused <- rep(TRUE, times = n2)
  for (i in seq_along(idx_query)) {
    idx <- idx_query[[i]]
    x_vec <- x_norm[idx, ]
    vals <- calculate_distance(x_vec, y_norm)
    chosen <- find_uniq_min_with_mask(vals, unused) # index of the nearest data point in y
    idx_out[[i]] <- chosen
    unused[[chosen]] <- FALSE
  }
  return(idx_out)
}

