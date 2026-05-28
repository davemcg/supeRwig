#' Launch supeRwig
#' 
#' @param bigwig_base Prefix prepended to each sample accession to
#'   build the BigWig URL or path. Trailing slash is added if absent.
#'   Use a local path like `"/data/bigwig/"` or a URL like
#'   `"hpc.nih.gov/~mcgaugheyd/eyeIntegration/2026/bigwig/"`. Files are assumed to be
#'   named `<bigwig_base><sample_accession><bigwig_ext>`.
#' @param bigwig_ext  BigWig file extension (default ".bw").
#' @param n_workers   Number of parallel workers for BigWig reads.
#' @param sj_se_dir Optional path to an HDF5-backed SummarizedExperiment
#'   of junction counts (as written by the SJ pipeline). When provided,
#'   a per-sample junction track is drawn below each wiggle. When NULL,
#'   the junction UI panel is hidden and the app behaves as before.
#' @return A [shiny::shinyApp()] object.
#' @export
supeRwig <- function(se_dir, bigwig_base, anno_fst,
                     sj_se_dir = NULL,
                     bigwig_ext = ".bw",
                     n_workers = max(1L, parallel::detectCores() - 1L)) {
  ctx <- build_context(
    se_dir      = se_dir,
    bigwig_base = bigwig_base,
    anno_fst    = anno_fst,
    sj_se_dir   = sj_se_dir,
    bigwig_ext  = bigwig_ext,
    n_workers   = n_workers
  )
  shiny::shinyApp(
    ui     = build_ui(ctx),
    server = build_server(ctx)
  )
}