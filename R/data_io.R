#' Assemble the app context (loaded data + worker pool)
#' @keywords internal
build_context <- function(se_dir, bigwig_base, anno_fst, bigwig_ext, n_workers, sj_se_dir = NULL) {
  se_bits <- load_se(se_dir)
  c(se_bits, load_annotation(anno_fst), load_junction_se(sj_se_dir), list(
    bw_file_map = build_bigwig_map(se_bits$eiad_meta$sample_accession, bigwig_base, bigwig_ext),
    bp_backend  = register_bp_backend(n_workers)
  ))
}

#' @keywords internal
register_bp_backend <- function(n_workers) {
  is_mac_or_win <- .Platform$OS.type == "windows" || Sys.info()["sysname"] == "Darwin"
  
  bp <- if (is_mac_or_win) {
    BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
  } else {
    BiocParallel::MulticoreParam(workers = n_workers)
  }
  
  BiocParallel::register(bp)
  bp
}

#' @keywords internal
load_se <- function(se_dir) {
  cat("Loading HDF5SummarizedExperiment from ", se_dir, "...\n")
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = se_dir)
  if (!"cpm" %in% SummarizedExperiment::assayNames(se)) stop("Assay 'cpm' not found.")
  
  rd <- SummarizedExperiment::rowData(se)
  gene_names <- if ("gene_name" %in% colnames(rd)) as.character(rd$gene_name) else rownames(se)
  
  eiad_meta <- unique(data.table::as.data.table(SummarizedExperiment::colData(se)))
  if (!"sample_accession" %in% colnames(eiad_meta)) eiad_meta[, sample_accession := colnames(se)]
  
  list(
    se = se, cpm_assay = SummarizedExperiment::assay(se, "cpm"),
    gene_to_se_row = split(seq_len(nrow(se)), gene_names),
    eiad_meta = eiad_meta, se_col_for_sample = setNames(seq_len(ncol(se)), eiad_meta$sample_accession)
  )
}

#' @keywords internal
load_junction_se <- function(sj_se_dir) {
  if (is.null(sj_se_dir)) return(list(sj_se = NULL, sj_row_meta = NULL, sj_col_for_sample = NULL))
  cat("Loading junction HDF5SummarizedExperiment from ", sj_se_dir, "...\n")
  sj_se <- HDF5Array::loadHDF5SummarizedExperiment(dir = sj_se_dir)
  
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(sj_se))
  if (!all(c("chr", "start", "end", "strand", "annot") %in% colnames(rd))) stop("Missing columns.")
  rd[, row_idx := .I]
  if (!"jid" %in% colnames(rd)) rd[, jid := rownames(sj_se)]
  data.table::setkey(rd, chr, start, end)
  
  list(sj_se = sj_se, sj_row_meta = rd, sj_col_for_sample = setNames(seq_len(ncol(sj_se)), colnames(sj_se)))
}

#' @keywords internal
load_annotation <- function(anno_fst) {
  cat("Loading annotation from ", anno_fst, "...\n")
  anno_dt <- fst::read_fst(anno_fst, as.data.table = TRUE)
  if ("type" %in% colnames(anno_dt)) anno_dt <- anno_dt[type %in% c("gene", "transcript", "exon")]
  if (!"transcript_id" %in% colnames(anno_dt)) anno_dt[, transcript_id := NA_character_]
  data.table::setkey(anno_dt, seqnames, start, end)
  list(anno_dt = anno_dt, unique_genes = sort(unique(anno_dt[type == "gene"]$gene_name)))
}

#' @keywords internal
build_bigwig_map <- function(sample_ids, bigwig_base, bigwig_ext) {
  is_url <- grepl("^https?://", bigwig_base, ignore.case = TRUE)
  if (!grepl("/$", bigwig_base)) bigwig_base <- paste0(bigwig_base, "/")
  paths <- setNames(paste0(bigwig_base, sample_ids, bigwig_ext), sample_ids)
  
  if (!is_url && any(missing <- !file.exists(paths))) {
    warning(sprintf("%d BigWig files missing.", sum(missing)))
  }
  paths
}

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
    Encoding(safe) <- "UTF-8"
    safe <- iconv(safe, "UTF-8", "UTF-8", sub = "")
    safe[is.na(safe) & !is.na(meta[[col]])] <- ""
    safe <- gsub("\"", "&quot;", gsub("'", "&#39;", safe, useBytes = TRUE), useBytes = TRUE)
    Encoding(safe) <- "UTF-8"
    data.table::set(meta, j = col, value = safe)
  }
  meta
}
