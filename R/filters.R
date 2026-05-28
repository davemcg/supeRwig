#' Filter samples by log2(CPM+1) of a target gene
#'
#' Runs before the BigWig read so we skip I/O for samples that won't
#' meet the cutoff anyway. Returns the input unchanged if no gene is
#' selected or the gene isn't in the SE.
#'
#' @keywords internal
cpm_filter_samples <- function(samples, gene_name, ctx, min_log_cpm) {
  if (is.null(gene_name) || !nzchar(gene_name)) return(samples)
  
  se_rows <- ctx$gene_to_se_row[[gene_name]]
  if (is.null(se_rows)) {
    warning("Gene '", gene_name, "' not in SE; skipping CPM filter.")
    return(samples)
  }
  
  cpm_vec <- if (length(se_rows) == 1) {
    as.numeric(ctx$cpm_assay[se_rows, ])
  } else {
    colSums(as.matrix(ctx$cpm_assay[se_rows, , drop = FALSE]))
  }
  names(cpm_vec) <- names(ctx$se_col_for_sample)
  
  log_cpm <- log2(cpm_vec[samples] + 1)
  keep    <- !is.na(log_cpm) & log_cpm >= min_log_cpm
  
  cat(sprintf("CPM filter (log2(CPM+1) >= %g for %s): %d/%d kept\n",
              min_log_cpm, gene_name, sum(keep), length(keep)))
  
  samples[keep]
}

#' Cap samples per (study x facet-group) combination
#'
#' @keywords internal
downsample_by_study <- function(meta, facet_cols, n) {
  if (n <= 0) return(meta)
  ds_by <- intersect(c("study_accession", facet_cols), colnames(meta))
  if (length(ds_by) == 0) return(meta)
  meta <- meta[order(sample_accession)]
  meta[, utils::head(.SD, n), by = ds_by]
}

#' Escape characters that would break ggiraph tooltips
#'
#' @keywords internal
sanitize_metadata <- function(meta) {
  char_cols <- names(meta)[vapply(meta, is.character, logical(1))]
  for (col in char_cols) {
    safe <- meta[[col]]
    # Force the encoding marker to UTF-8 before round-tripping. Without
    # this, strings marked "unknown" get interpreted in the current
    # locale by iconv and bad bytes survive sub="".
    Encoding(safe) <- "UTF-8"
    safe <- iconv(safe, from = "UTF-8", to = "UTF-8", sub = "")
    # If iconv still returns NA for a non-NA input (rare, version-dependent),
    # fall back to empty string so downstream nchar/gsub don't trip.
    bad <- is.na(safe) & !is.na(meta[[col]])
    if (any(bad)) safe[bad] <- ""
    safe <- gsub("'",  "&#39;",  safe, useBytes = TRUE)
    safe <- gsub("\"", "&quot;", safe, useBytes = TRUE)
    Encoding(safe) <- "UTF-8"   # gsub(useBytes=TRUE) drops the mark
    data.table::set(meta, j = col, value = safe)
  }
  meta
}