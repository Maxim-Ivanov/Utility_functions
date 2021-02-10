# sort(grl) sorts only ranges within list elements;
# This function reorders list elements by their range;

sort_grl <- function(grl) {
  orig_names <- names(grl)
  names(grl) <- 1:length(grl)
  nms <- names(sort(unlist(range(grl))))
  idx <- match(nms, names(grl))
  out <- grl[idx]
  names(out) <- orig_names[idx]
  return(out)
}
