remove_first_and_last_exons <- function(ebg, mode = "3") {
  ebg <- sort(ebg)
  ebg_unl <- unlist(ebg)
  ebg_red <- GenomicRanges::reduce(ebg)
  ebg_red_unl <- unlist(ebg_red)
  enum_fw <- unlist(sapply(lengths(ebg_red), function(x) { seq(1, x) }))
  enum_rev <- unlist(sapply(lengths(ebg_red), function(x) { seq(x, 1) }))
  first_red <- ebg_red_unl[as.logical((strand(ebg_red_unl) == "+" & enum_fw == 1) | (strand(ebg_red_unl) == "-" & enum_rev == 1))]
  h1 <- findOverlaps(ebg_unl, first_red)
  same_gene_1 <- names(ebg_unl[queryHits(h1)]) == names(first_red[subjectHits(h1)])
  idx_1 <- unique(queryHits(h1[same_gene_1]))
  first_exons <- ebg_unl[idx_1]
  exons_wo_first <- ebg_unl[-idx_1]
  last_red <- ebg_red_unl[as.logical((strand(ebg_red_unl) == "+" & enum_rev == 1) | (strand(ebg_red_unl) == "-" & enum_fw == 1))]
  h2 <- findOverlaps(ebg_unl, last_red)
  same_gene_2 <- names(ebg_unl[queryHits(h2)]) == names(last_red[subjectHits(h2)])
  idx_2 <- unique(queryHits(h2[same_gene_2]))
  last_exons <- ebg_unl[idx_2]
  exons_wo_last <- ebg_unl[-idx_2]
  internal_exons <- ebg_unl[-unique(c(idx_1, idx_2))]
  if (mode == "1") {
    out <- internal_exons
  } else if (mode == "3") {
    out <- list(exons_wo_first, exons_wo_last, internal_exons)
  } else if (mode == "5") {
    out <- list(exons_wo_first, exons_wo_last, internal_exons, first_exons, last_exons)
  }
  return(out)
}

find_splice_sites <- function(ebg) {
  res <- remove_first_and_last_exons(ebg)
  exons_wo_first <- res[[1]]
  exons_wo_last <- res[[2]]
  ss_5p <- suppressWarnings(trim(reduce(resize(resize(exons_wo_last, 1, "end"), 3, "center"))))
  ss_3p <- suppressWarnings(trim(reduce(resize(flank(exons_wo_first, 1), 3, "center"))))
  return(list(ss_5p, ss_3p))
}
