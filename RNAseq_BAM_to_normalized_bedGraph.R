RNAseq_BAM_to_normalized_bedGraph <- function(bamfile, mode = "PE", stranded = TRUE, switch_strand = TRUE, normalize = TRUE, norm_to = 1e06) {
  stopifnot(mode %in% c("PE", "SE"))
  if (mode == "PE") {
    ga <- readGAlignmentPairs(bamfile)
  } else {
    ga <- readGAlignments(bamfile)
  }
  cov <- convert_GAlignments_to_coverage(ga, merge.strands = !stranded, flip.strands = switch_strand, normalize = normalize, norm_to = norm_to)
  return(cov)
}
