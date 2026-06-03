#' Subset annotation to features overlapping the window
#' @keywords internal
subset_region_annotation <- function(anno_dt, chr, w_start, w_end) {
  region_anno <- anno_dt[seqnames == chr & start <= w_end & end >= w_start & type %in% c("transcript", "exon")]
  has_tx <- nrow(region_anno) > 0
  
  if (has_tx) {
    region_anno[, `:=`(tx_label = paste0(gene_name, " - ", transcript_id), 
                       tx_idx = as.numeric(as.factor(paste0(gene_name, " - ", transcript_id))))]
    
    # Find transcript IDs that explicitly contain the principal tag
    principal_ids <- character(0)
    if ("tag" %in% colnames(region_anno)) {
      principal_ids <- unique(region_anno[grepl("GENCODE_Primary", tag), transcript_id])
    } else if ("appris" %in% colnames(region_anno)) {
      principal_ids <- unique(region_anno[grepl("appris_principal_1", appris), transcript_id])
    } 
    
    # Broadcast flag to both the transcript lines and exon blocks
    region_anno[, is_principal := transcript_id %in% principal_ids]
    # ---------------------------------
    
    exon_hi <- data.table::as.data.table(GenomicRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(region_anno[type == "exon"])))
  } else {
    exon_hi <- tx_base <- tx_exons <- data.table::data.table()
  }
  
  list(
    region_anno = region_anno, 
    tx_base = region_anno[type == "transcript"], 
    tx_exons = region_anno[type == "exon"], 
    exon_highlights = exon_hi, 
    has_transcripts = has_tx
  )
}

#' @keywords internal
.tooltip_safe_cols <- function(meta_cur, max_median_chars = 120) {
  candidate <- setdiff(colnames(meta_cur), c("sample_accession", "dummy_facet", "combined_facet"))
  too_long <- vapply(candidate, function(c) {
    v <- meta_cur[[c]]
    if (!is.character(v) && !is.factor(v)) return(FALSE)
    stats::median(nchar(as.character(v), type = "bytes"), na.rm = TRUE) > max_median_chars
  }, logical(1))
  candidate[!too_long]
}

#' @keywords internal
build_plot_data <- function(dt_full, meta_cur, facet_cols, overlap_factor, junc_band = 0, tooltip_max_chars = 120) {
  if (length(facet_cols) == 0) { meta_cur$dummy_facet <- "All Samples"; facet_cols <- "dummy_facet" }
  
  tissue_map <- unique(meta_cur)[!duplicated(sample_accession)]
  tissue_map[, combined_facet := do.call(paste, c(.SD, sep = " - ")),
             .SDcols = facet_cols]
  
  # Vectorized tooltip construction
  tips <- paste0("<b>Sample:</b> ", tissue_map$sample_accession)
  for (col in .tooltip_safe_cols(meta_cur, tooltip_max_chars)) {
    s <- as.character(tissue_map[[col]])
    long <- !is.na(s) & nchar(s, type = "bytes") > tooltip_max_chars
    s[long] <- paste0(substr(s[long], 1, tooltip_max_chars - 1), "\u2026")
    tips <- paste0(tips, "<br><b>", col, ":</b> ", s)
  }
  tissue_map[, static_tooltip := tips]
  
  pd <- merge(dt_full, tissue_map, by.x = "sample", by.y = "sample_accession", all.x = TRUE)
  unique_samples <- unique(pd[, .(sample, combined_facet)])[order(combined_facet, sample)][, local_idx := seq_len(.N), by = combined_facet]
  
  pd <- merge(pd, unique_samples, by = c("sample", "combined_facet"))
  pd[, `:=`(log_val = log2(value + 1), plot_x = (binned_pos + bin_end) / 2)]
  max_log <- max(pd$log_val, na.rm = TRUE); if (max_log == 0 || is.na(max_log)) max_log <- 1
  
  pd[, offset_y := (log_val / max_log) * overlap_factor + local_idx + junc_band]
  pd[, tooltip_text := paste0(static_tooltip, ": ", round(value, 2))]
  data.table::setorderv(pd, c("combined_facet", "local_idx", "plot_x"))
  
  list(plot_data = pd, unique_samples = unique_samples)
}

#' @keywords internal
build_bed_tooltips <- function(bed_dt) {
  if (nrow(bed_dt) == 0) return(bed_dt[, bed_tooltip := character(0)])
  tip <- sprintf("<b>%s:</b> %s-%s (%.2f kb)", bed_dt$seqnames, format(bed_dt$start, big.mark = ","), format(bed_dt$end, big.mark = ","), (bed_dt$end - bed_dt$start + 1) / 1000)
  if ("name" %in% colnames(bed_dt)) tip <- paste0(tip, "<br><b>Name:</b> ", bed_dt$name)
  if ("score" %in% colnames(bed_dt)) tip <- paste0(tip, "<br><b>Score:</b> ", bed_dt$score)
  if ("strand" %in% colnames(bed_dt)) tip <- paste0(tip, "<br><b>Strand:</b> ", bed_dt$strand)
  bed_dt[, bed_tooltip := tip]
}