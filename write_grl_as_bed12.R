write_grl_as_bed12 <- function(grl, filename, name_in_mcols = FALSE) {
  gr <- range(grl) %>% unlist(use.names = FALSE)
  grl_unl <- unlist(grl, use.names = FALSE)
  idx_fw <- lapply(elementNROWS(grl), function(x) { seq(1, x) }) %>% unlist() %>% unname()
  idx_rev <- lapply(elementNROWS(grl), function(x) { seq(x, 1) }) %>% unlist() %>% unname()
  idx <- ifelse(strand(grl_unl) == "+", idx_fw, idx_rev)
  mc <- grl_unl[idx == 1] %>% mcols()
  chrom <- seqnames(gr) %>% as.character()
  chromStart <- start(gr)
  chromEnd <- end(gr)
  if (isTRUE(name_in_mcols)) {
    if (is.null(mc$name)) {
      name <- "."
    } else {
      name <- mc$name
    }
  } else {
    if (is.null(names(grl))) {
      name <- "."
    } else {
      name <- names(grl)
    }
  }
  if (is.null(mc$score)) {
    score <- "."
  } else {
    score <- mc$score
  }
  strand <- strand(gr) %>% as.character()
  if (is.null(mc$thick)) {
    thickStart <- chromStart
    thickEnd <- chromEnd
  } else {
    thickStart <- mc$thick %>% start()
    thickEnd <- mc$thick %>% end()
  }
  if (is.null(mc$itemRgb)) {
    itemRgb <- "."
  } else {
    itemRgb <- mc$itemRgb
  }
  blockCount <- elementNROWS(grl)
  blockSizes <- width(grl) %>% lapply(str_c, collapse = ",") %>% unlist()
  blockStarts <- start(grl) %>% `-`(start(gr)) %>% lapply(str_c, collapse = ",") %>% unlist()
  tbl <- tibble(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)
  tbl <- arrange(tbl, chrom, chromStart, chromEnd)
  write_tsv(tbl, filename, col_names = FALSE)
}
