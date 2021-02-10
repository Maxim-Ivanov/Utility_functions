skip_duplicated_reads <- function(bam) {
  nms <- names(bam)
  idx <- which(duplicated(nms) | duplicated(nms, fromLast = TRUE))
  return(bam[-idx])
}