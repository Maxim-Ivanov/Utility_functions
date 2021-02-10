# The function getOverlappingScores() integrates Bedgraph/Bigwig signal over each genomic interval of interest. It calls findOverlaps() function from GenomicRanges; 
# Input 1: intervals of interest (as GRanges object);
# Input 2: signal (as GrangesList; each GRanges corresponds to a separate Bedgraph of Bigwig file);
# If signal is stranded, then "+" and "-" values are expected to be pooled within a single GRanges object. By default, findOverlaps() takes the strand info into account when looking for overlaps;
# Input 3: number of characters to strip from the beginning and the end of each filename;
# Output: Granges object with additional metadata columns containing the integrated Bedgraph signal;
# Works correctly even when stranded intervals are combined with unstranded data, or vice versa; 

get_overlapping_scores <- function(intervals, signal_grl, value = "mcols", round = 6) {
  require(GenomicRanges)
  results <- vector("list", length(signal_grl))
  old_names <- names(intervals)
  names(intervals) <- 1:length(intervals)
  int_plus <- intervals[strand(intervals)=="+"]
  int_minus <- intervals[strand(intervals)=="-"]
  int_star <- intervals[strand(intervals)=="*"]
  intervals_split <- list(int_plus, int_minus, int_star)
  strands <- list("+", "-", c("+", "-"))
  for (i in seq_along(signal_grl)) {
    signal <- signal_grl[[i]]
    message(names(signal_grl)[[i]]); flush.console()
    res <- vector("list", 3)
    for (j in 1:3) {
      curr_strands <- c(strands[[j]], "*")
      curr_signal <- signal[strand(signal) %in% curr_strands]
      curr_intervals <- intervals_split[[j]]
      curr_cov <- coverage(curr_signal, weight="score")
      cov_by_intervals <- curr_cov[curr_intervals]
      numlist <- as(cov_by_intervals, "NumericList")
      integrated <- lapply(numlist, FUN=sum)
      names(integrated) <- names(curr_intervals)
      res[[j]] <- integrated
    }
    res <- Reduce(c, res)
    res <- res[order(as.numeric(names(res)))]
    res <- as.numeric(res)
    if (is.numeric(round) && length(round) == 1) {
      res <- round(res, round)
    }
    results[[i]] <- res
  }
  results <- t(do.call(rbind, results))
  colnames(results) <- names(signal_grl)
  if (value == "mcols") {
    mcols(intervals) <- cbind(mcols(intervals), as.data.frame(results))
    names(intervals) <- old_names
    return(intervals)
  } else if (value == "count_matrix") {
    rownames(results) <- old_names
    return(results)
  }
}

