extract_first_or_last_element <- function(nested_list, mode = "first") {
  unlisted <- unlist(nested_list)
  idx <- rep(1:length(nested_list), times = elementNROWS(nested_list))
  if (mode == "first") {
    choose <- !duplicated(idx)
  } else if (mode == "last") {
    choose <- !duplicated(idx, fromLast = TRUE)
  }
  return(unlisted[choose])
}

filter_ga_by_terminal_subalignments <- function(ga, abs_threshold = 30, rel_threshold = 0.05) {
  grl <- grglist(ga)
  spl <- elementNROWS(grl) >= 2
  grl_spl <- grl[spl]
  unal <- psetdiff(unlist(range(grl_spl)), grl_spl)
  w_aln <- width(grl_spl)
  w_unal <- width(unal)
  first_aln <- extract_first_or_last_element(w_aln)
  last_aln <- extract_first_or_last_element(w_aln, mode = "last")
  first_unal <- extract_first_or_last_element(w_unal)
  last_unal <- extract_first_or_last_element(w_unal, mode = "last")
  if (!is.null(abs_threshold)) {
    good_abs <- first_aln >= abs_threshold & last_aln >= abs_threshold
    good_unspl <- width(unlist(grl[!spl])) >= abs_threshold
  } else {
    good_abs <- rep(TRUE, sum(spl))
    good_unspl <- rep(TRUE, sum(!spl))
  }
  if (!is.null(rel_threshold)) {
    good_rel <- (first_aln / first_unal >= rel_threshold) & (last_aln / last_unal >= rel_threshold)
  } else (
    good_rel <- rep(TRUE, length(grl_spl))
  )
  good_spl <- good_abs & good_rel
  final_good <- rep(FALSE, length(ga))
  final_good[!spl][good_unspl] <- TRUE
  final_good[spl][good_spl] <- TRUE
  message(length(ga), " -> ", sum(final_good))
  return(ga[final_good])
}
