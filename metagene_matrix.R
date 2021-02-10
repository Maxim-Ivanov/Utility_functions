# The intervals are expected to be stranded;
# The signal can be either stranded (e.g. NET-Seq) or unstranded (e.g. ChIP-Seq);
# In the antisense mode (antisenseMode=TRUE) the intervals are flipped to the opposite strand to investigate the antisense transcription;
#
# If at least one interval has a different length, then intervals can be either: a) scaled (extended or shrinked) to a certain common length (scaling = TRUE), or b) remain unscaled (scaling = FALSE; default).
# In the latter case (scaling = FALSE), longer intervals are trimmed and shorter intervals are filled with NA values. The intervals can be anchored either at the start (default), or at the end (anchor = c("start", "end"));
# In both cases it is required to determine the length of the output matrix. The matrix.length argument can be either an arbitrary number, or "max" (width of the longest interval), or "min" (width of the shortest interval), or the mean interval length (rounded to the nearest integer), otherwise the median interval length is used by default;
#
# na.as.zeros parameter allows to substitute NA values with zeros. By default it is FALSE, thus the output matrix may contain NA values (given that scaling = FALSE). This has to be taken into account when applying statistical functions to the output matrix;
#
# skip.zeros parameter allows to remove intervals with zero signal (not a single tag along the whole interval);
#
# skip.outliers parameter allows to remove intervals with the highest signal (e.g. the top 0.5% intervals are considered as potential outliers with skip.outliers=0.995);
#
# equal.weights allows to normalize the integral signal over each interval to 1 (i.e. intervals with low coverage and with high coverage contribute equally to the final metagene profile). This may be useful when the signal is expected to consist of sharp discrete peaks (e.g. TSS-Seq), and their relative positions within intervals are of interest (whereas their height is considered not important);

### ADD CHECK: if matrix.length is different from interval sizes, then scaling must automatically turn on. 

findOutOfBounds <- function(gr) {
  lookup <- as.list(seqlengths(seqinfo(gr)))
  lengths <- lookup[as.character(seqnames(gr))]
  out <- (start(gr) < 0) | (end(gr) > lengths)
  return(out)
}

calcMatrixLength <- function(intervals, matrix.length) {
  max_width <- max(width(intervals))
  min_width <- min(width(intervals))
  if (is.numeric(matrix.length)) {
    out <- matrix.length
  } else if (max_width == min_width) {
    out <- max_width
  } else if (matrix.length == "max") {
    out <- max_width
  } else if (matrix.length == "min") {
    out <- min_width
  } else if (matrix.length == "mean") {
    out <- round(mean(width(intervals)))
  } else if (matrix.length == "median") {
    out <- round(median(width(intervals)))
  } else {
    message("Not sure how to determine the number of bins. Using median interval length by default..."); flush.console()
    out <- round(median(width(intervals)))
  }
  message("Matrix length = ", out); flush.console()
  return(out)
}

expandOrShrink <- function(x, mlen) {
  if (length(x) == mlen) {
    return(as.numeric(x))
  } else {
    runLength(x) <- runLength(x) * mlen
    return(colMeans(matrix(x, ncol = mlen)))###
  }
}

trimOrFill <- function(x, mlen, anchor, na.as.zeros) {
  if (length(x) == mlen) {
    out <- x
  } else if (length(x) > mlen) {
    if (anchor == "start") {
      out <- x[1:mlen]
    } else if (anchor == "end") {
      out <- x[(length(x)-mlen+1):length(x)]
    }
  } else {
    if (isTRUE(na.as.zeros)) { vals <- 0 } else { vals <- NA }
    to_add <- Rle(values = vals, lengths = mlen - length(x))
    if (anchor == "start") {
      out <- c(x, to_add)
    } else if (anchor == "end") {
      out <- c(to_add, x)
    }
  }
  return(as.numeric(out))
}


calcMatrix <- function(signal, intervals, strand = NULL, max_w, min_w, scaling, mlen, anchor, na.as.zeros) {
  if (!is.null(strand)) {
    intervals <- intervals[strand(intervals) == strand]
    signal <- signal[strand(signal) == strand]
  }
  rlelist <- round(coverage(signal, weight = "score"), 8) # to avoid small negative and positive values instead of zeros
  rlelist_split <- rlelist[intervals]
  rlelist_split <- revElements(rlelist_split, strand(intervals) == "-")
  if ((max_w != min_w) || (max_w != mlen)) {
    if (isTRUE(scaling)) {
      if (length(rlelist_split) > 0) { message("Expanding and shrinking (this can be slow)..."); flush.console() }
      numlist <- lapply(rlelist_split, expandOrShrink, mlen = mlen)
    } else {
      if (length(rlelist_split) > 0) { message("Trimming and filling..."); flush.console() }
      numlist <- lapply(rlelist_split, trimOrFill, mlen = mlen, anchor = anchor, na.as.zeros = na.as.zeros)
    }
  } else {
    numlist <- as(rlelist_split, "NumericList")
  }
  mat <- do.call(rbind, numlist)
  return(mat)
}

average_bin <- function(x, shrink_method) {
  if (all(is.na(x))) {
    return(NA)
  } else if (shrink_method == "mean") {
    return(sum(x, na.rm = TRUE) / length(x))
  } else if (shrink_method == "median") {
    return(median(x, na.rm = TRUE))
  } else if (shrink_method == "mode") {
    freq <- table(x)
    return(as.numeric(names(freq)[which.max(freq)]))# the most frequently observed value; can be very slow
    #return(as.numeric(names(sort(-table(x)))[1]))
  }
}

shrink_row <- function(x, mask, shrink_method) {
  return(as.numeric(tapply(x, mask, average_bin, shrink_method = shrink_method, simplify = FALSE)))
}

shrink_to_bins <- function(mat, binsize, shrink_method) {
  mask <- rep(seq(1, ncol(mat) / binsize), each = binsize)
  out <- t(apply(mat, 1, shrink_row, mask = mask, shrink_method = shrink_method))
  return(out)
}

metagene_matrix <- function(signal, intervals, scaling = FALSE, matrix.length = NA, anchor = "start", na.as.zeros = FALSE, 
                           skip.zeros = TRUE, skip.outliers = 0.995, skip.top.obs = FALSE, n.top.obs = 3, 
                           equal.weights = FALSE, antisenseMode = FALSE, shrink = FALSE, binsize = 5, shrink_method = "mean") {
  library(GenomicRanges)
  out_of_bound <- findOutOfBounds(intervals)
  if (sum(out_of_bound) > 0) {
    message(sum(out_of_bound), " intervals were out-of-bounds;"); flush.console()
    intervals <- intervals[!out_of_bound]
  }
  names(intervals) <- 1:length(intervals) # enumerate intervals
  mlen <- calcMatrixLength(intervals, matrix.length)
  max_w <- max(width(intervals)); min_w <- min(width(intervals))
  if (any(strand(intervals)=="*")) {
    message(round(sum(strand(intervals)=="*") / length(intervals) * 100, 1), "% intervals are unstranded!"); flush.console()
  }
  if (any(strand(signal)=="*")) {
    message(round(sum(strand(signal)=="*") / length(signal) * 100, 1), "% signal values are unstranded!"); flush.console()
  }
  if (isTRUE(antisenseMode)) {
    stranded <- strand(signal) %in% c("+", "-")
    strand(signal)[stranded] <- ifelse(strand(signal)[stranded]=="+", "-", "+")
    message("Antisense mode: signal was flipped to the opposite strand;"); flush.console()
  }
  interval_strands <- list("+", "-", "*")
  signal_strands <- list(c("+", "*"), c("-", "*"), c("+", "-", "*"))
  matlist <- vector("list", 3)
  for (i in 1:3) {
    curr_signal <- signal[strand(signal) %in% signal_strands[[i]]]
    curr_intervals <- intervals[strand(intervals) %in% interval_strands[[i]]]
    if (length(curr_intervals) > 0) {
      curr_mat <- calcMatrix(signal = curr_signal, intervals = curr_intervals, max_w = max_w, min_w = min_w, 
                             scaling = scaling, mlen = mlen, anchor = anchor, na.as.zeros = na.as.zeros)
      rownames(curr_mat) <- names(curr_intervals) # to trace back the original interval (because matrix rows are permuted at this step)
      matlist[[i]] <- curr_mat
    }
  }
  mat <- do.call(rbind, matlist)
  mat <- mat[order(as.numeric(rownames(mat))), ] # restore the original order of intervals
  if (isTRUE(shrink) && is.numeric(binsize) && length(binsize) == 1 && shrink_method %in% c("mean", "median", "mode")) {
    if (ncol(mat) %% binsize == 0) {
      message("Averaging signal within ", binsize, " bp bins;")
      mat <- shrink_to_bins(mat, binsize, shrink_method) # when it is required to average signal within N bp bins but scaling == FALSE
      message("Now matrix length is ", ncol(mat), "!")
    } else {
      message("Check the binsize parameter!")
    }
  }
  gene_cov <- rowSums(mat, na.rm = TRUE)
  if (is.numeric(skip.outliers) & length(skip.outliers) == 1) {
    q <- quantile(gene_cov, skip.outliers)
    outliers <- gene_cov > q
    mat <- mat[!outliers, ]
    message("Skipped ", sum(outliers), " potential outliers;"); flush.console()
    gene_cov <- rowSums(mat, na.rm = TRUE) # recalculate gene_cov
  }
  if (isTRUE(skip.top.obs) && is.numeric(n.top.obs) && length(n.top.obs) == 1) {
    mat_sorted <- apply(mat, 2, sort, decreasing = TRUE, na.last = TRUE) # sort each column independently from other columns
    max_vals <- mat_sorted[n.top.obs, ] # find N'th top value
    for (j in 1:ncol(mat)) {
      to_skip <- mat[, j] >= max_vals[j]
      mat[to_skip, j] <- NA # in each column, change N top observations to NA
    }
    message("Skipped ", n.top.obs, " top observations in each bin;")
    gene_cov <- rowSums(mat, na.rm = TRUE)
  }
  if (isTRUE(skip.zeros)) {
    zeros <- gene_cov == 0
    if (sum(!zeros) < 100) {
      stop("Too little intervals with non-zero signal! Consider changing skip.zeros to FALSE!")
    }
    if (sum(zeros) > 0) {
      mat <- mat[!zeros, ]
      message("Skipped ", sum(zeros), " intervals with zero signal;"); flush.console()
    }
  }
  if (isTRUE(equal.weights)) {
    mat <- t(apply(mat, 1, function(x) { x / sum(x, na.rm=TRUE) }))
    message("All intervals were assigned equal weights;"); flush.console()
  }
  return(mat)
}

