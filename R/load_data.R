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
  bp <- if (.Platform$OS.type == "windows") BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
  else BiocParallel::MulticoreParam(workers = n_workers)
  BiocParallel::register(bp); bp
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