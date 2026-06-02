#' Read junctions overlapping a region for a sample subset
#'
#' Filters in-memory rowData first (cheap), then performs a dual HDF5 read
#' of both the PSI and raw count (junction x sample) submatrices. `start >= w_start &
#' end <= w_end` is fully-contained (matching how the SE was built).
#'
#' Empty data.tables are returned with the correct schema so callers
#' can `nrow() == 0` rather than NULL-check.
#'
#' @keywords internal
read_region_junctions <- function(ctx, chr, w_start, w_end, samples,
                                  min_psi5 = 0, min_psi3 = 0) {
  empty <- data.table::data.table(
    jid          = character(0),
    sample       = character(0),
    count        = numeric(0),
    raw_count    = integer(0),
    psi5         = numeric(0),
    psi3         = numeric(0),
    start        = integer(0),
    end          = integer(0),
    strand       = character(0),
    annot        = integer(0),
    strand_annot = character(0),
    cluster_5    = character(0),
    cluster_3    = character(0)
  )
  if (is.null(ctx$sj_se)) return(empty)
  
  rd <- ctx$sj_row_meta
  q_chr   <- chr
  q_start <- w_start
  q_end   <- w_end
  
  # CHANGE: Allow partially overlapping junctions by checking interval intersection
  hits <- rd[chr == q_chr & start <= q_end & end >= q_start]
  if (nrow(hits) == 0) return(empty)
  
  sample_cols <- ctx$sj_col_for_sample[samples]
  sample_cols <- sample_cols[!is.na(sample_cols)]
  if (length(sample_cols) == 0) return(empty)
  
  # Fetch BOTH fractional splice usage matrix values
  m_psi5 <- as.matrix(SummarizedExperiment::assay(ctx$sj_se, "psi5")[hits$row_idx, sample_cols, drop = FALSE])
  m_psi3 <- as.matrix(SummarizedExperiment::assay(ctx$sj_se, "psi3")[hits$row_idx, sample_cols, drop = FALSE])
  rownames(m_psi5) <- hits$jid; colnames(m_psi5) <- names(sample_cols)
  rownames(m_psi3) <- hits$jid; colnames(m_psi3) <- names(sample_cols)
  
  dt_psi5 <- data.table::as.data.table(m_psi5, keep.rownames = "jid")
  long_psi5 <- data.table::melt(dt_psi5, id.vars = "jid", variable.name = "sample", value.name = "psi5_val")
  
  dt_psi3 <- data.table::as.data.table(m_psi3, keep.rownames = "jid")
  long_psi3 <- data.table::melt(dt_psi3, id.vars = "jid", variable.name = "sample", value.name = "psi3_val")
  
  long_psi <- merge(long_psi5, long_psi3, by = c("jid", "sample"))
  long_psi[, sample := as.character(sample)]
  
  if ("counts" %in% SummarizedExperiment::assayNames(ctx$sj_se)) {
    m_cts <- as.matrix(SummarizedExperiment::assay(ctx$sj_se, "counts")[hits$row_idx, sample_cols, drop = FALSE])
    rownames(m_cts) <- hits$jid; colnames(m_cts) <- names(sample_cols)
    m_cts_dt <- data.table::as.data.table(m_cts, keep.rownames = "jid")
    long_cts <- data.table::melt(m_cts_dt, id.vars = "jid", variable.name = "sample", value.name = "raw_count")
    long_cts[, sample := as.character(sample)]
    long <- merge(long_psi, long_cts, by = c("jid", "sample"))
  } else {
    long <- long_psi
    long[, raw_count := NA_integer_]
  }
  
  # Convert 0-10000 scaled integers to UI percentages (0-100%)
  long[, psi5 := psi5_val / 100]
  long[, psi3 := psi3_val / 100]
  
  # Enforce strict AND filtration boundaries AND strip absolute zeros
  if (!all(is.na(long$raw_count))) {
    long <- long[(psi5 >= min_psi5 & psi3 >= min_psi3) & raw_count > 0]
  } else {
    long <- long[psi5 >= min_psi5 & psi3 >= min_psi3]
  }
  
  if (nrow(long) == 0) return(empty)
  
  # Set line thickness priority ordering to whichever value is higher
  long[, count := pmax(psi5, psi3)]
  
  extra_metadata <- intersect(c("SYMBOL", "cluster_5", "cluster_3"), colnames(hits))
  long <- merge(long, hits[, c("jid", "start", "end", "strand", "annot", extra_metadata), with = FALSE], by = "jid")
  long[, strand_annot := paste0(strand, "/", ifelse(annot == 1L, "annot", "novel"))]
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

#' Attach y-positions and visual attrs to junctions with Interval Packing
#'
#' @keywords internal
attach_junction_positions <- function(junctions, unique_samples,
                                      junc_band = 0.65) {
  if (nrow(junctions) == 0) return(junctions)
  
  junc <- merge(junctions, unique_samples, by = "sample")
  
  # Helper function to compute optimal layout tracks using greedy interval packing
  pack_lanes <- function(start_vec, end_vec) {
    if (length(start_vec) == 0) return(integer(0))
    
    # Sort left-to-right by start coordinate
    ord <- order(start_vec)
    s_sorted <- start_vec[ord]
    e_sorted <- end_vec[ord]
    
    # Track the rightmost coordinate mapped to each lane
    lane_ends <- numeric(0)
    assigned_lanes <- integer(length(start_vec))
    
    # 25bp visual buffer prevents adjacent lines from colliding
    buffer <- 25L 
    
    for (i in seq_along(s_sorted)) {
      placed <- FALSE
      # Find the first available lane that ends before this junction starts
      for (l in seq_along(lane_ends)) {
        if (s_sorted[i] > (lane_ends[l] + buffer)) {
          lane_ends[l] <- e_sorted[i]
          assigned_lanes[ord[i]] <- l - 1L
          placed <- TRUE
          break
        }
      }
      # If all current lanes are blocked by overlapping junctions, open a new lane
      if (!placed) {
        lane_ends <- c(lane_ends, e_sorted[i])
        assigned_lanes[ord[i]] <- length(lane_ends) - 1L
      }
    }
    return(assigned_lanes)
  }
  
  # Calculate packed lanes independently per sample track
  junc[, sub_idx := pack_lanes(start, end), by = .(combined_facet, local_idx)]
  
  # --- Inverted vertical layout --------
  sub_spacing <- 0.060  # Snugged up line spacing slightly for better density
  n_visible   <- max(1L, as.integer(floor((junc_band - 0.08) / sub_spacing)))
  junc[, sub_idx := sub_idx %% n_visible]
  
  # Instead of adding to floor, subtract downward from the wiggle track baseline
  junc[, junc_y  := local_idx + junc_band - 0.07 - (sub_idx * sub_spacing)]
  
  junc[, junc_lw := 0.4]
  
  # ---- Vectorized Rich Tooltip Generation  ----
  junc[, junc_tooltip := paste0(
    "<b>Junction:</b> ", jid,
    "<br><b>Sample:</b> ",  sample,
    "<br><b>PSI5 (5' Donor Focus):</b> ", round(psi5, 2), "%",
    "<br><b>PSI3 (3' Acceptor Focus):</b> ", round(psi3, 2), "%"
  )]
  
  if ("raw_count" %in% colnames(junc)) {
    junc[, junc_tooltip := paste0(junc_tooltip, "<br><b>Raw Count:</b> ", 
                                  data.table::fifelse(is.na(raw_count), "N/A", as.character(raw_count)))]
  }
  
  if ("cluster_5" %in% colnames(junc)) {
    junc[, junc_tooltip := paste0(junc_tooltip, "<br><b>5' Cluster (Donor):</b> ", cluster_5)]
  } 
  if ("cluster_3" %in% colnames(junc)) {
    junc[, junc_tooltip := paste0(junc_tooltip, "<br><b>3' Cluster (Acceptor):</b> ", cluster_3)]
  } 
  
  junc[, junc_tooltip := paste0(
    junc_tooltip,
    "<br><b>Strand:</b> ",  strand,
    "<br><b>Status:</b> ",  ifelse(annot == 1L, "annotated", "novel")
  )]
  
  junc
}