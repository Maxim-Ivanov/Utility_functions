RNAseq_BAM_to_normalized_bedGraph <- function(bamfile, mode = "PE", stranded = TRUE, switch_strand = TRUE, normalize = TRUE, norm_to = 1e06, skip_over = NULL) {
  stopifnot(mode %in% c("PE", "SE"))
  if (mode == "PE") {
    ga <- readGAlignmentPairs(bamfile)
  } else {
    ga <- readGAlignments(bamfile)
  }
  if (!is.null(skip_over) && class(skip_over) == "GRanges") {
    if (!isTRUE(stranded)) {
      strand(skip_over) <- "*"
    } else if (isTRUE(switch_strand)) {
      strand(skip_over) <- ifelse(strand(skip_over) == "+", "-", "+")
    }
    bad <- ga %over% skip_over
    message(sum(bad), " (", round(mean(bad) * 100, 1), "%) reads were skipped due to overlap with unwanted intervals;")
    ga <- ga[!bad]
  }
  cov <- convert_GAlignments_to_coverage(ga, merge.strands = !stranded, flip.strands = switch_strand, normalize = normalize, norm_to = norm_to)
  return(cov)
}
