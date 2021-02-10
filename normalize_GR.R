# GRanges is expected to have scores
normalize_GR <- function(gr, norm = 1000000, norm_factor = NULL, precision = 6, by = NULL) {
  gr <- gr[mcols(gr)$score > 0]
  if (is.null(norm_factor)) {
    if (!is.null(by) && class(by) == "GRanges") {
      rel_reads <- gr[gr %over% by]
    } else {
      rel_reads <- gr
    }
    auc <- sum(mcols(rel_reads)$score * width(rel_reads))
    norm_factor <- auc / norm
  } else if (!is.numeric(norm_factor) || length(norm_factor) != 1) {
    stop("Invalid norm_factor!")
  }
  mcols(gr)$score <- round(mcols(gr)$score / norm_factor, precision)
  return(gr)
}
