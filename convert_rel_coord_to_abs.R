# gr = stranded GRanges with absolute (genomic) coordinates of all exons within a transcript;
# ir = IRanges with relative (transcriptomic) coordinates of any intervals of interest within given transcript;
# convert_rel_coord_to_abs() converts transcriptomic coordinates into genomic ones;
# The strandedness of gr is taken into account (i.e. relative coordinates are interpreted in 5'->3' direction);
# Usage example:
# grl <- exonsBy(txdb, by = "tx", use.names = TRUE)
# irl <- ranges(tx_intervals_grl)
# out_grl <- mapply(convert_rel_coord_to_abs, grl, irl, SIMPLIFY = FALSE) %>% as("CompressedGRangesList")

convert_positions <- function(pos_ir, gr, cumw, str_gr) {
  # Find indexes of parent ranges in gr:
  idx <- lapply(pos_ir, function(x) { idx <- which(x - cumw > 0); return(idx[length(idx)]) }) %>% unlist()
  # Offset from the parent range border:
  offset <- pos_ir - cumw[idx]
  # Convert relative position to absolute position:
  if (str_gr == "+") {
    abs_pos <- start(gr)[idx] + offset - 1
  } else {
    abs_pos <- end(gr)[idx] - offset + 1
  }
  return(abs_pos)
}

convert_rel_coord_to_abs <- function(gr, ir) {
  str_gr <- unique(strand(gr))
  chr_gr <- unique(seqnames(gr))
  # Check if all ranges in gr are on the same strand of the same chromosome:
  stopifnot(length(str_gr) == 1 && str_gr %in% c("+", "-") && length(chr_gr) == 1)
  if (length(gr) > length(GenomicRanges::reduce(gr))) {
    warning("Overlapping ranges detected in gr!")
    gr <- GenomicRanges::reduce(gr)
  }
  # Ensure that ranges in gr are properly sorted:
  if (str_gr == "+") {
    gr <- sort(gr)
  } else {
    gr <- sort(gr, decreasing = TRUE)
  }
  # Convert relative coordinated to absolute coordinates:
  cumw <- cumsum(width(gr)) - width(gr) # cumulative width
  abs_start <- convert_positions(start(ir), gr, cumw, str_gr)
  abs_end <- convert_positions(end(ir), gr, cumw, str_gr)
  # Construct output GRanges:
  out <- GRanges(seqnames = chr_gr, IRanges(start = abs_start, end = abs_end), strand = str_gr, seqinfo = seqinfo(gr))
  names(out) <- names(ir)
  mcols(out) <- mcols(ir)
  # Split ranges in out which overlap with gaps in gr (otherwise they will extend over introns):
  gaps <- GenomicRanges::setdiff(range(gr), gr)
  if (length(gaps) > 0 && any(out %over% gaps)) {
    p1 <- out[out %outside% gaps]
    p2 <- out[out %over% gaps]
    hits <- findOverlaps(p2, gaps)
    # It may happen so that one range covers more than one gap, thus unique() and split():
    p2_par <- p2[unique(queryHits(hits))]
    gaps_par <- gaps[subjectHits(hits)] %>% split(queryHits(hits))
    p2_spl <- GenomicRanges::psetdiff(p2_par, gaps_par)
    # Restore mcols which was lost at psetdiff() call:
    mc <- mcols(p2_par)
    idx <- rep(1:nrow(mc), times = elementNROWS(p2_spl))
    p3 <- unlist(p2_spl)
    mcols(p3) <- mc[idx, ]
    names(p3) <- names(p2_par)[idx]
    out <- c(p1, p3)
    if (str_gr == "+") {
      out <- sort(out)
    } else {
      out <- sort(out, decreasing = TRUE)
    }
  }
  return(out)
}
