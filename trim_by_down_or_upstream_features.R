trim_by_down_or_upstream_features <- function(windows, features, mode = "down", offset = 10, ignore.strand = FALSE, check_parent = TRUE, trim_to_zero_width = TRUE) {
  stopifnot(mode %in% c("down", "up"))
  old_names <- names(windows)
  if (isTRUE(check_parent)) {
    stopifnot(!is.null(names(windows)) && !is.null(names(features))) # both windows and features are expected to have names
  }
  mcols(windows)$orig_order <- 1:length(windows) # enumerate windows
  if (offset > 0) {
    if (mode == "down") {
      features <- trim(resize(features, width = width(features) + offset, fix = "end")) # extend features upstream by the offset value
    } else if (mode == "up") {
      features <- trim(resize(features, width = width(features) + offset, fix = "start"))
    }
  }
  over_any <- overlapsAny(windows, features, ignore.strand = ignore.strand)
  out1 <- windows[!over_any] # exclude and save windows which do not overlap any feature
  win_rem <- windows[over_any]
  hits <- findOverlaps(win_rem, features, ignore.strand = ignore.strand)
  win_par <- win_rem[queryHits(hits)]
  feat_par <- features[subjectHits(hits)]
  if (isTRUE(check_parent)) {
    over_parent <- names(win_par) == names(feat_par) # names should match between a window and its parental feature to exclude guaranteed overlaps
    only_parent <- as.logical(tapply(over_parent, mcols(win_par)$orig_order, all)) # find windows which overlap with the parent only
    out2 <- win_rem[only_parent] # exclude and save windows which overlap with the parent only
    win_rem <- win_rem[!only_parent]
    win_par <- win_par[!over_parent] # continue with windows which overlap at least one non-parental feature
    feat_par <- feat_par[!over_parent]
  } else {
    out2 <- win_rem[NULL]
  }
  over_type <- parallelOverlapType(win_par, feat_par) # detect the type of overlap
  if (mode == "down") {
    #bad <- as.logical(tapply(over_type, list(mcols(win_par)$orig_order), function(x) { any(x == "up") | any(x == "contains") }))
    bad <- as.logical(BiocGenerics::tapply(over_type, list(S4Vectors::mcols(win_par)$orig_order), function(x) { any(x %in% c("up", "contains", "exact", "no_down")) }))
  } else if (mode == "up") {
    #bad <- as.logical(tapply(over_type, list(mcols(win_par)$orig_order), function(x) { any(x == "down") | any(x == "contains") }))
    bad <- as.logical(BiocGenerics::tapply(over_type, list(S4Vectors::mcols(win_par)$orig_order), function(x) { any(x %in% c("down", "contains", "exact", "no_up")) }))
  }
  out3 <- win_rem[bad] # exclude and save windows which overlap any feature in undesired orientation
  if (isTRUE(trim_to_zero_width)) {
    out3 <- resize(out3, width = 0) # trim such windows to zero width
  }
  out4 <- win_rem[!bad] # continue with windows which overlap at least one downstream or upstream feature
  if (mode == "down") {
    out4_zw <- resize(out4, width = 0, fix = "start") # resize both windows and features to 1 bp width around start
    feat_zw <- resize(features, width = 0, fix = "start")
    if (isTRUE(ignore.strand)) { strand(feat_zw) <- "*" }
    feat_near <- features[precede(out4_zw, feat_zw)] # extract the nearest downstream feature start
  } else if (mode == "up") {
    out4_zw <- resize(out4, width = 0, fix = "end")
    feat_zw <- resize(features, width = 0, fix = "end")
    if (isTRUE(ignore.strand)) { strand(feat_zw) <- "*" }
    feat_near <- features[follow(out4_zw, feat_zw)] # extract the nearest upstream feature start
  }
  out4_trimmed <- pgap(out4_zw, feat_near, ignore.strand = ignore.strand) # compute truncated windows as gap between the starts
  mcols(out4_trimmed) <- mcols(out4) # restore mcols which were removed in the pgap operation
  out <- c(out1, out2, out3, out4_trimmed) # combine all processed windows
  out <- out[order(mcols(out)$orig_order)] # restore the original order
  mcols(out)$orig_order <- NULL
  names(out) <- old_names
  return(out)
}
