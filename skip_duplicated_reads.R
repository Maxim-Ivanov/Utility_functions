skip_duplicated_reads <- function(bam) {
  nms <- names(bam)
  if (is.null(nms)) {
    stop("GAlignments object must have read names!")
  }
  message(length(bam), " input BAM records;")
  idx <- which(duplicated(nms) | duplicated(nms, fromLast = TRUE))
  if (length(idx) > 0) {
    message("Skipped ", length(idx), " BAM records with duplicated names (", length(unique(nms[idx])), " chimeric reads);")
    bam <- bam[-idx]
  }
  return(bam)
}
