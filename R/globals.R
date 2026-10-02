# Pull data.table fully into the package namespace. Without this, the
# [.data.table S3 method isn't reliably dispatched from inside package
# functions and `dt[col == val]` silently falls through to [.data.frame,
# which evaluates `col` as a variable and errors with "object 'col'
# not found". This is the canonical fix recommended in
# vignette("datatable-importing", package = "data.table").
#' @import data.table
NULL

# Declare data.table column references that look like undefined variables
# to R CMD check. Without this, the package generates dozens of
# "no visible binding for global variable" notes.

utils::globalVariables(c(
  # existing
  "type", "gene_name", "transcript_id", "seqnames", "start", "end", "strand",
  "sample", "sample_accession", "value", "binned_pos", "bin_end",
  "log_val", "tx_idx", "tx_label", "local_idx",
  "combined_facet", "static_tooltip", "line_color", "offset_y", "plot_x",
  "tooltip_text", "dummy_facet", "draw_start", "draw_end", "label_x",
  ".SD", ".N", ".I", "is_principal", "fill_val",
  # junction-layer additions
  "jid", "annot", "count", "strand_annot",
  "sub_idx", "junc_y", "junc_lw", "junc_tooltip", "row_idx",
  # mean-averaging additions
  "n_samples", "n_samples_total", "col_val"
))


#' Parse a UCSC-style region string ("chr1:93,992,834-94,121,148")
#' @keywords internal
parse_ucsc_region <- function(s) {
  if (is.null(s) || !nzchar(s)) return(NULL)
  s <- gsub("[,[:space:]]", "", s)
  m <- regmatches(s, regexec("^([^:]+):(\\d+)-(\\d+)$", s))[[1]]
  if (length(m) != 4L) return(NULL)
  start <- suppressWarnings(as.integer(m[3]))
  end   <- suppressWarnings(as.integer(m[4]))
  if (is.na(start) || is.na(end) || start >= end) return(NULL)
  list(chr = m[2], start = start, end = end)
}

#' Format chr/start/end back into UCSC-style with thousands separators
#' @keywords internal
format_ucsc_region <- function(chr, start, end) {
  sprintf("%s:%s-%s", chr,
          format(start, big.mark = ",", scientific = FALSE),
          format(end,   big.mark = ",", scientific = FALSE))
}


#' Shared plot theme
#' @keywords internal
theme_panel_only <- function() {
  cowplot::theme_minimal_vgrid() +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(t = 5, r = 8, b = 5, l = 8)
    )
}

#' Categorical palette for arbitrary numbers of levels
#' @keywords internal
cat_palette <- function(levels) {
  n <- length(levels)
  if (n == 0) return(character(0))
  pal <- c(pals::cols25()[-c(6,7,13,14)], 
           pals::polychrome()[-c(1,2,20)],
           pals::glasbey(),
           pals::okabe())[seq_len(n)]
  setNames(unname(pal), levels)
}

# Single-line null-coalescing operator to clean up default parameter handling
`%||%` <- function(a, b) if (is.null(a)) b else a