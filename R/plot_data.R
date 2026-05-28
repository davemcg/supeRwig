#' Subset annotation to features overlapping the window
#'
#' @keywords internal
subset_region_annotation <- function(anno_dt, chr, w_start, w_end) {
  region_anno <- anno_dt[seqnames == chr &
                           start <= w_end &
                           end   >= w_start &
                           type %in% c("transcript", "exon")]
  has_tx <- nrow(region_anno) > 0
  
  if (has_tx) {
    region_anno[, tx_label := paste0(gene_name, " - ", transcript_id)]
    region_anno[, tx_idx   := as.numeric(as.factor(tx_label))]
    tx_base  <- region_anno[type == "transcript"]
    tx_exons <- region_anno[type == "exon"]
    exon_gr  <- GenomicRanges::makeGRangesFromDataFrame(tx_exons)
    exon_hi  <- data.table::as.data.table(GenomicRanges::reduce(exon_gr))
  } else {
    exon_hi  <- data.table::data.table(start = numeric(0), end = numeric(0))
    tx_base  <- data.table::data.table()
    tx_exons <- data.table::data.table()
  }
  
  list(
    region_anno     = region_anno,
    tx_base         = tx_base,
    tx_exons        = tx_exons,
    exon_highlights = exon_hi,
    has_transcripts = has_tx
  )
}

#' Join coverage with metadata, compute per-sample y-offsets and tooltips
#'
#' @keywords internal
build_plot_data <- function(dt_full, meta_cur, facet_cols, overlap_factor,
                            summary_type, junc_band = 0) {
  if (is.null(facet_cols) || length(facet_cols) == 0) {
    meta_cur$dummy_facet <- "All Samples"
    facet_cols <- "dummy_facet"
  }
  
  tissue_map <- unique(meta_cur)
  tissue_map <- tissue_map[!duplicated(sample_accession)]
  facet_header <- paste(facet_cols, collapse = " - ")
  tissue_map[, combined_facet := paste0(
    facet_header, "\n",
    do.call(paste, c(.SD, sep = " - "))
  ), .SDcols = facet_cols]
  
  meta_cols <- setdiff(colnames(meta_cur),
                       c("sample_accession", "dummy_facet", "combined_facet"))
  
  tissue_map[, static_tooltip := paste0("<b>Sample:</b> ", sample_accession)]
  for (col in meta_cols) {
    tissue_map[, static_tooltip := paste0(static_tooltip,
                                          "<br><b>", col, ":</b> ", get(col))]
  }
  
  dt_full <- merge(dt_full, tissue_map,
                   by.x = "sample", by.y = "sample_accession", all.x = TRUE)
  
  unique_samples <- unique(dt_full[, c("sample", "combined_facet"),
                                   with = FALSE])
  data.table::setorderv(unique_samples, c("combined_facet", "sample"))
  unique_samples[, local_idx := seq_len(.N), by = combined_facet]
  
  pd <- merge(dt_full, unique_samples, by = c("sample", "combined_facet"))
  
  pd[, log_val := log2(value + 1)]
  max_log <- max(pd$log_val, na.rm = TRUE)
  if (max_log == 0 || is.na(max_log)) max_log <- 1
  
  # Shift wiggle baseline up by junc_band so the strip [local_idx,
  # local_idx + junc_band] is reserved for the junction layer.
  pd[, offset_y := (log_val / max_log) * overlap_factor +
       local_idx + junc_band]
  pd[, plot_x   := (binned_pos + bin_end) / 2]
  data.table::setorderv(pd, c("combined_facet", "local_idx", "plot_x"))
  pd[, tooltip_text := paste0(static_tooltip,
                              "<br><b>", summary_type, ":</b> ",
                              round(value, 2))]
  
  list(plot_data = pd, unique_samples = unique_samples)
}

#' Build per-region tooltip text for uploaded BED highlights
#'
#' Defensive about which BED columns are present: only `seqnames`,
#' `start`, `end` are required; `name`, `score`, `strand` and any
#' extras are appended if present. Adds a `bed_tooltip` column.
#'
#' @keywords internal
build_bed_tooltips <- function(bed_dt) {
  if (nrow(bed_dt) == 0) {
    bed_dt[, bed_tooltip := character(0)]
    return(bed_dt)
  }
  width_kb <- (bed_dt$end - bed_dt$start + 1) / 1000
  tip <- paste0(
    "<b>", bed_dt$seqnames, ":</b> ",
    format(bed_dt$start, big.mark = ","), "-",
    format(bed_dt$end,   big.mark = ","),
    " (", formatC(width_kb, format = "f", digits = 2), " kb)"
  )
  if ("name"   %in% colnames(bed_dt))
    tip <- paste0(tip, "<br><b>Name:</b> ",   bed_dt$name)
  if ("score"  %in% colnames(bed_dt))
    tip <- paste0(tip, "<br><b>Score:</b> ",  bed_dt$score)
  if ("strand" %in% colnames(bed_dt))
    tip <- paste0(tip, "<br><b>Strand:</b> ", as.character(bed_dt$strand))
  bed_dt[, bed_tooltip := tip]
  bed_dt
}