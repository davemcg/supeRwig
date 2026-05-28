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
  ".SD", ".N", ".I",
  # junction-layer additions
  "jid", "annot", "count", "strand_annot",
  "sub_idx", "junc_y", "junc_lw", "junc_tooltip", "row_idx"
))