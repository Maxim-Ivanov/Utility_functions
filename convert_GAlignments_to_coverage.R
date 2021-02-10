convert_GAlignments_to_coverage <- function(ga, split = TRUE, mode = "whole_read", merge.strands = FALSE, 
                                         flip.strands = FALSE, normalize = FALSE, norm_to = 1000000L, precision = 6) {
  require(GenomicAlignments)
  if (mode %in% c("start", "end")) {
    gr <- resize(granges(ga), width = 1, fix = mode)
  } else if (mode == "whole_read") {
    if (isTRUE(split)) {
      gr <- unlist(grglist(ga), use.names = FALSE)
    } else {
      gr <- granges(ga)
    }
  }
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
  out <- out[order(seqnames(out), start(out))]
  if (isTRUE(normalize)) {
    norm_factor <- length(ga) / norm_to
    score(out) <- round(score(out) / norm_factor, precision)
  }
  return(out)
}
