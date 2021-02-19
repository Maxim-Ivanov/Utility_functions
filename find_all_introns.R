# Input: exons grouped by gene;
# Traditional solution: psetdiff(unlist(range(ebg)), ebg)
# This solution does not work good for genes containing overlapping exons...
# The function below works correct even for genes with overlapping exons.
# The output grl is parallel to the input grl (can contain elements with zero length for intronless genes); 

find_all_introns <- function(ebg, verbose = FALSE) {
  old_names <- names(ebg)
  names(ebg) <- 1:length(ebg)
  # Divide the input ebg into genes with and without overlapping exons:
  over <- lengths(GenomicRanges::reduce(ebg)) != lengths(disjoin(ebg))
  grl0 <- ebg[!over]
  grl <- ebg[over]
  if (verbose == TRUE) {
    message(length(grl), " (", round(mean(over) * 100, 2), "%) input genes have overlapping exons;")
  }
  # Process grl0 in the traditional way:
  out0 <- psetdiff(unlist(range(grl0)), grl0)
  # Unlist grl:
  gene_idx <- rep(seq_along(grl), lengths(grl)) # index of original genes
  gr <- granges(unlist(grl, use.names = FALSE))
  abs_idx <- seq_along(gr) # absolute index
  # Unlist also the reduced version of grl:
  grl_red <- GenomicRanges::reduce(grl)
  subject_idx <- rep(seq_along(grl_red), lengths(grl_red))
  gr_red <- unlist(grl_red, use.names = FALSE)
  # Find all overlaps:
  hits <- findOverlaps(gr, gr_red)
  # Filter to hits within the same gene:
  same_gene <- gene_idx[queryHits(hits)] == subject_idx[subjectHits(hits)]
  hits <- hits[same_gene]
  over_idx <- subjectHits(hits) # overlap index
  # Split the absolute index by the overlap index:
  abs_by_over <- split(abs_idx, over_idx)
  # Shorten the gene index to fit abs_by_over:
  gene_by_over <- split(gene_idx, over_idx)
  gene_by_over <- unlist(vapply(gene_by_over, `[[`, integer(1), 1)) # extract first elements
  # Compute outer product from absolute indexes of each two adjacent exons:
  abs_wo_last <- abs_by_over[-length(abs_by_over)]
  abs_wo_first <- abs_by_over[-1]
  res0 <- Map(function(x, y) { outer(x, y, paste, sep = ",") }, abs_wo_last, abs_wo_first)
  # Create parallel vector of gene indexes:
  gene_wo_last <- gene_by_over[-length(gene_by_over)]
  gene_wo_first <- gene_by_over[-1]
  res1 <- Map(function(x, y) { if (x == y) return(x)}, gene_wo_last, gene_wo_first) # NULL values where x != y
  # Remove elements which correspond to "introns" between adjacent genes:
  remove <- vapply(res1, is.null, logical(1))
  res <- res0[!remove]
  r1 <- unlist(res)
  # Expand gene indexes to fit the unlisted outer products:
  r2 <- rep(unlist(res1), lengths(res))
  # Convert pairwise combinations to integer vectors:
  spl <- strsplit(r1, ",")
  idx1 <- lapply(spl, `[[`, 1) # absolute indexes of first exons in each pair
  idx2 <- lapply(spl, `[[`, 2) # ... of second exons ...
  # Extract first and second exons from the input gr:
  exons1 <- gr[as.integer(idx1)]
  exons2 <- gr[as.integer(idx2)]
  # Calculate gaps between exons:
  introns <- pgap(exons1, exons2)
  # Get names of original genes:
  gene_names <- names(grl)[as.integer(r2)]
  # Convert introns to grl:
  out1 <- splitAsList(introns, gene_names)
  # Remove duplicated introns within each gene:
  out1 <- unique(out1)
  # There can be also genes which have overlapping exons but do not produce any introns:
  unused <- setdiff(names(grl), names(out1))
  out2 <- split(GRanges(), unused)
  # Restore the original order of genes:
  out <- c(out0, out1, out2)
  out <- out[order(as.integer(names(out)))]
  names(out) <- old_names
  return(out)
}
