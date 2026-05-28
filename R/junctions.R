#' Read junctions overlapping a region for a sample subset
#'
#' Filters in-memory rowData first (cheap), then performs ONE HDF5 read
#' of the resulting (junction x sample) submatrix. `start >= w_start &
#' end <= w_end` is fully-contained (matching how the SE was built);
#' switch to overlap if you want junctions whose anchors sit outside
#' the window too.
#'
#' Empty data.tables are returned with the correct schema so callers
#' can `nrow() == 0` rather than NULL-check.
#'
#' @keywords internal
read_region_junctions <- function(ctx, chr, w_start, w_end, samples,
                                  min_reads) {
  empty <- data.table::data.table(
    jid          = character(0),
    sample       = character(0),
    count        = integer(0),
    start        = integer(0),
    end          = integer(0),
    strand       = character(0),
    annot        = integer(0),
    strand_annot = character(0)
  )
  if (is.null(ctx$sj_se)) return(empty)
  
  # rowData has columns named chr/start/end that collide with the
  # function args of the same name. Rename locals so data.table's
  # NSE finds the column on the LHS and the variable on the RHS,
  # avoiding the fragile `..chr` prefix (which is documented for `j`
  # but not reliably supported in `i` across data.table versions).
  rd <- ctx$sj_row_meta
  q_chr   <- chr
  q_start <- w_start
  q_end   <- w_end
  hits <- rd[chr == q_chr & start >= q_start & end <= q_end]
  if (nrow(hits) == 0) return(empty)
  
  sample_cols <- ctx$sj_col_for_sample[samples]
  sample_cols <- sample_cols[!is.na(sample_cols)]
  if (length(sample_cols) == 0) return(empty)
  
  m <- as.matrix(
    SummarizedExperiment::assay(ctx$sj_se, "counts")[hits$row_idx,
                                                     sample_cols,
                                                     drop = FALSE]
  )
  rownames(m) <- hits$jid
  colnames(m) <- names(sample_cols)
  
  m_dt <- data.table::as.data.table(m, keep.rownames = "jid")
  long <- data.table::melt(m_dt, id.vars = "jid",
                           variable.name = "sample",
                           value.name    = "count")
  long[, sample := as.character(sample)]
  long <- long[count >= min_reads]
  if (nrow(long) == 0) return(empty)
  
  long <- merge(long,
                hits[, .(jid, start, end, strand, annot)],
                by = "jid")
  
  long[, strand_annot := paste0(strand, "/",
                                ifelse(annot == 1L, "annot", "novel"))]
  long
}

#' Strand x annotation color palette for the junction layer
#'
#' Saturated = annotated, pale = novel. Blue = +, red = -, grey = `*`.
#'
#' @keywords internal
junction_palette <- function() {
  c(
    "+/annot" = "royalblue4",
    "+/novel" = "royalblue1",
    "-/annot" = "tomato4",
    "-/novel" = "tomato1",
    "*/annot" = "seagreen4",
    "*/novel" = "seagreen1"
  )
}

#' Attach y-positions and visual attrs to junctions
#'
#' Junctions are placed in the thin band between each sample's
#' baseline (`local_idx`) and the shifted wiggle zero line
#' (`local_idx + junc_band`). Within that band, junctions are
#' stacked into a few sub-rows ordered by read count; extras wrap
#' onto earlier sub-rows (lowest-count first, so the dominant
#' junctions stay on their own line).
#'
#' @keywords internal
attach_junction_positions <- function(junctions, unique_samples,
                                      junc_band = 0.25) {
  if (nrow(junctions) == 0) return(junctions)
  
  junc <- merge(junctions, unique_samples, by = "sample")
  data.table::setorder(junc, combined_facet, local_idx, -count)
  junc[, sub_idx := seq_len(.N) - 1L, by = .(combined_facet, local_idx)]
  
  sub_spacing <- 0.045
  n_visible   <- max(1L, as.integer(floor((junc_band - 0.08) / sub_spacing)))
  junc[, sub_idx := sub_idx %% n_visible]
  junc[, junc_y  := local_idx + 0.04 + sub_idx * sub_spacing]
  
  #junc[, junc_lw      := pmin(2.2, 0.4 + log10(count + 1) * 0.7)]
  junc[, junc_lw := 0.4]
  junc[, junc_tooltip := paste0(
    "<b>Junction:</b> ", jid,
    "<br><b>Sample:</b> ",  sample,
    "<br><b>Reads:</b> ",   count,
    "<br><b>Strand:</b> ",  strand,
    "<br><b>Status:</b> ",  ifelse(annot == 1L, "annotated", "novel")
  )]
  junc
}