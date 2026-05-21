#' Assemble the app context (loaded data + worker pool)
#'
#' Bundles everything the UI and server need so we can pass a single
#' object through the call graph instead of carrying ~8 globals.
#'
#' @keywords internal
build_context <- function(se_dir, bigwig_base, anno_fst, bigwig_ext,
                          n_workers) {
  bp_backend <- register_bp_backend(n_workers)
  se_bits    <- load_se(se_dir)
  anno_bits  <- load_annotation(anno_fst)
  # Build BigWig URLs/paths from the SE's sample_accession column.
  # Must run AFTER load_se because we need se_bits$eiad_meta.
  bw_map     <- build_bigwig_map(se_bits$eiad_meta$sample_accession,
                                 bigwig_base, bigwig_ext)
  
  c(se_bits, anno_bits, list(
    bw_file_map = bw_map,
    bp_backend  = bp_backend
  ))
}

#' @keywords internal
register_bp_backend <- function(n_workers) {
  bp <- if (.Platform$OS.type == "windows") {
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
  cat(sprintf("  %d rows x %d samples\n", nrow(se), ncol(se)))
  
  if (!"cpm" %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay 'cpm' not found. Available assays: ",
         paste(SummarizedExperiment::assayNames(se), collapse = ", "))
  }
  cpm_assay <- SummarizedExperiment::assay(se, "cpm")
  
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  gene_names <- if ("gene_name" %in% colnames(rd)) {
    as.character(rd$gene_name)
  } else {
    rownames(se)
  }
  gene_to_se_row <- split(seq_len(nrow(se)), gene_names)
  
  cd <- data.table::as.data.table(
    as.data.frame(SummarizedExperiment::colData(se))
  )
  if (!"sample_accession" %in% colnames(cd)) {
    cd[, sample_accession := colnames(se)]
  }
  eiad_meta <- unique(cd)
  data.table::setDT(eiad_meta)
  se_col_for_sample <- setNames(seq_len(ncol(se)), eiad_meta$sample_accession)
  
  list(
    se                = se,
    cpm_assay         = cpm_assay,
    gene_to_se_row    = gene_to_se_row,
    eiad_meta         = eiad_meta,
    se_col_for_sample = se_col_for_sample
  )
}

#' @keywords internal
load_annotation <- function(anno_fst) {
  cat("Loading annotation from ", anno_fst, "...\n")
  anno_dt <- fst::read_fst(anno_fst, as.data.table = TRUE)
  if ("type" %in% colnames(anno_dt)) {
    anno_dt <- anno_dt[type %in% c("gene", "transcript", "exon")]
  }
  for (req in c("type", "seqnames", "start", "end", "strand", "gene_name")) {
    if (!req %in% colnames(anno_dt))
      stop("Annotation fst is missing required column: ", req)
  }
  if (!"transcript_id" %in% colnames(anno_dt))
    anno_dt[, transcript_id := NA_character_]
  
  list(
    anno_dt      = anno_dt,
    unique_genes = sort(unique(anno_dt[type == "gene"]$gene_name))
  )
}

#' Compose per-sample BigWig URLs (or local paths) from the SE
#'
#' BigWig file locations are constructed as
#' `<bigwig_base><sample_accession><bigwig_ext>` rather than discovered
#' by directory scan. This lets the app read directly from HTTPS,
#' S3-style endpoints, or any other URL scheme rtracklayer's
#' BigWigFile understands, without requiring the files to be
#' co-located on the Shiny host.
#'
#' For local prefixes only, the function checks file existence and
#' warns if any are missing. URL prefixes are not pre-checked because
#' that would mean one HEAD request per sample at startup.
#'
#' @keywords internal
build_bigwig_map <- function(sample_ids, bigwig_base, bigwig_ext) {
  is_url <- grepl("^https?://", bigwig_base, ignore.case = TRUE)
  # Be lenient about a forgotten trailing separator.
  if (!grepl("/$", bigwig_base)) bigwig_base <- paste0(bigwig_base, "/")
  
  paths <- paste0(bigwig_base, sample_ids, bigwig_ext)
  names(paths) <- sample_ids
  
  cat(sprintf("Built BigWig map for %d samples (%s%s)\n",
              length(paths),
              if (is_url) "remote: " else "local: ",
              bigwig_base))
  
  if (!is_url) {
    missing <- !file.exists(paths)
    if (any(missing)) {
      ex <- utils::head(paths[missing], 3)
      warning(sprintf(
        "%d of %d BigWig files do not exist on disk. Examples: %s",
        sum(missing), length(paths), paste(ex, collapse = ", ")
      ))
    }
  }
  
  paths
}