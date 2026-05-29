#' Read one BigWig as a binned summary
#'
#' @keywords internal
.bw_summary_one <- function(path, sname, gr, n_bins, type) {
  tryCatch({
    bwf  <- rtracklayer::BigWigFile(path)
    bins <- rtracklayer::summary(bwf, gr, size = n_bins, type = type)[[1]]
    scores <- as.numeric(bins$score)
    
    list(
      sample     = rep.int(sname, length(scores)),
      binned_pos = GenomicRanges::start(bins),
      bin_end    = GenomicRanges::end(bins),
      value      = scores
    )
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s",
                    basename(path), conditionMessage(e)))
    NULL
  })
}

#' Read a genomic region across many BigWigs in parallel
#'
#' @keywords internal
read_region_bigwigs <- function(samples, bw_file_map, chr, start, end,
                                n_bins, bp_backend) {
  target_files  <- bw_file_map[samples]
  missing_files <- samples[is.na(target_files)]
  if (length(missing_files) > 0) {
    warning(sprintf("No BigWig file for %d samples (e.g. %s)",
                    length(missing_files),
                    paste(utils::head(missing_files, 3), collapse = ", ")))
  }
  target_files <- target_files[!is.na(target_files)]
  if (length(target_files) == 0)
    stop("No BigWig files matched the post-filter samples.")
  
  cat(sprintf("Reading %d BigWigs at %d bins...\n",
              length(target_files), n_bins))
  t0 <- Sys.time()
  
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = start, end = end)
  )
  
  results <- BiocParallel::bpmapply(
    FUN       = .bw_summary_one,
    path      = unname(target_files),
    sname     = names(target_files),
    MoreArgs  = list(gr     = gr,
                     n_bins = n_bins,
                     type   = 'max'),
    SIMPLIFY  = FALSE,
    USE.NAMES = FALSE,
    BPPARAM   = bp_backend
  )
  
  dt <- data.table::rbindlist(results)
  if (nrow(dt) == 0)
    stop("No coverage data could be read for this region.")
  dt[is.na(value) | is.nan(value), value := 0]
  
  cat(sprintf("Read complete in %.2f s\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  dt
}