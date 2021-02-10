merge_and_normalize_GRanges <- function(x, norm = 1000000, precision = 6, min_support = 1) {
  stranded_input <- unlist(lapply(x, function(gr) { all(strand(gr) %in% c("+", "-")) }))
  if (any(!stranded_input)) {
    stop("Input GRanges are unstranded!")
  }
  # x is expected to be a list of stranded GRanges with "score" column
  cov_fw <- lapply(x, function(x) { coverage(x[strand(x) == "+"], weight = "score") })
  cov_fw <- Reduce(`+`, cov_fw)
  cov_rev <- lapply(x, function(x) { coverage(x[strand(x) == "-"], weight = "score") })
  cov_rev <- Reduce(`+`, cov_rev)
  if (min_support > 1) {
    sup_fw <- lapply(x, function(x) { coverage(x[strand(x) == "+"]) })
    sup_fw <- Reduce(`+`, sup_fw)
    cov_fw[sup_fw > 0 & sup_fw < min_support] <- 0
    sup_rev <- lapply(x, function(x) { coverage(x[strand(x) == "-"]) })
    sup_rev <- Reduce(`+`, sup_rev)
    cov_rev[sup_rev > 0 & sup_rev < min_support] <- 0
  }
  cov_fw <- bindAsGRanges(score = cov_fw)
  strand(cov_fw) <- "+"
  cov_rev <- bindAsGRanges(score = cov_rev)
  strand(cov_rev) <- "-"
  out <- sort(c(cov_fw, cov_rev))
  if (is.numeric(norm) && length(norm) == 1) {
    auc <- sum(width(out) * score(out))
    norm_factor <- auc / norm
    score(out) <- score(out) / norm_factor
  }
  score(out) <- round(score(out), precision)
  out <- out[score(out) > 0]
  return(out)
}