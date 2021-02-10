convert_GRanges_to_coverage <- function(gr, merge.strands = FALSE, flip.strands = FALSE, normalize = FALSE, norm_to = 1000000L, norm_by = "input_reads", precision = 5) {
  if (isTRUE(merge.strands)) {
    cov_all <- coverage(gr)
    out <- bindAsGRanges(score = cov_all)
  } else {
    cov_fw <- coverage(gr[strand(gr) == "+"])
    cov_rev <- coverage(gr[strand(gr) == "-"])
    out_fw <- bindAsGRanges(score = cov_fw)
    strand(out_fw) <- "+"
    out_rev <- bindAsGRanges(score = cov_rev)
    strand(out_rev) <- "-"
    out <- c(out_fw, out_rev)
    if (isTRUE(flip.strands)) {
      strand(out) <- ifelse(strand(out) == "+", "-", "+")
    }
  }
  out <- out[score(out) > 0]
  out <- sortSeqlevels(out)
  out <- sort(out)
  if (isTRUE(normalize)) {
    if (norm_by == "input_reads") {
      norm_factor <- length(gr) / norm_to
    } else if (norm_by == "total_area") {
      norm_factor <- sum(width(out) * score(out)) / norm_to
    }
    score(out) <- round(score(out) / norm_factor, precision)
  }
  return(out)
}