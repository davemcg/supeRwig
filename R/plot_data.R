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

#' Decide which metadata columns can safely go into the tooltip
#'
#' Drops any column whose median character length exceeds
#' `max_median_chars`. This is a safety net for SRA-style fields like
#' `sample_attribute` and `study_abstract` that grow to multiple KB
#' per sample -- embedding those in every interactive SVG element
#' causes ggiraph to produce enormous SVGs and the browser to hang.
#'
#' Only structural columns are excluded by name: `sample_accession`
#' is already shown explicitly in the tooltip, and `dummy_facet` /
#' `combined_facet` are internal bookkeeping. Everything else is
#' judged on length alone.
#'
#' @keywords internal
.tooltip_safe_cols <- function(meta_cur, max_median_chars = 120) {
  structural <- c("sample_accession", "dummy_facet", "combined_facet")
  candidate  <- setdiff(colnames(meta_cur), structural)
  too_long   <- vapply(candidate, function(c) {
    v <- meta_cur[[c]]
    if (!is.character(v) && !is.factor(v)) return(FALSE)
    # type = "bytes" never validates multibyte sequences, so it cannot
    # throw "invalid multibyte string" on dirty SRA fields. Byte length
    # is a safe upper bound for character length, which is what we want
    # here anyway (the goal is to keep tooltip SVG small).
    stats::median(nchar(as.character(v), type = "bytes"),
                  na.rm = TRUE) > max_median_chars
  }, logical(1))
  candidate[!too_long]
}

#' Truncate long values for safe embedding in HTML tooltips
#' @keywords internal
.trunc_value <- function(x, max_chars = 120) {
  s <- as.character(x)
  n <- nchar(s, type = "bytes")
  long <- !is.na(s) & n > max_chars
  s[long] <- paste0(substr(s[long], 1, max_chars - 1), "\u2026")
  s
}

#' Join coverage with metadata, compute per-sample y-offsets and tooltips
#'
#' Tooltip text is generated from a vetted subset of metadata columns
#' (see `.tooltip_safe_cols`) with each value truncated to
#' `tooltip_max_chars`. Without this, free-text SRA fields like
#' `sample_attribute` and `study_abstract` blow up the SVG by 100x+ and
#' the browser hangs during render.
#'
#' @keywords internal
build_plot_data <- function(dt_full, meta_cur, facet_cols, overlap_factor,
                            summary_type, junc_band = 0,
                            tooltip_max_chars = 120) {
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

  # Only include short, plottable columns in the tooltip. Long
  # free-text fields are silently dropped (you can still see them by
  # inspecting the SE directly).
  tooltip_cols <- .tooltip_safe_cols(meta_cur, max_median_chars = tooltip_max_chars)
  tooltip_cols <- setdiff(tooltip_cols, c("dummy_facet", "combined_facet"))

  tissue_map[, static_tooltip := paste0("<b>Sample:</b> ", sample_accession)]
  for (col in tooltip_cols) {
    tissue_map[, static_tooltip := paste0(
      static_tooltip,
      "<br><b>", col, ":</b> ",
      .trunc_value(get(col), tooltip_max_chars)
    )]
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