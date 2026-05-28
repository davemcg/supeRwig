#' Build the main wiggle ggplot
#'
#' @param plot_data        Per-bin coverage with offsets and tooltips.
#' @param exon_highlights  Reduced exon ranges to shade behind the traces.
#' @param bed_highlights   Optional user-uploaded BED ranges (subset to
#'                         the current window). Empty data.table = no
#'                         BED layer drawn.
#' @param bed_color  Fill color and opacity for BED bands.
#' @param junctions Optional data.table with columns `sample, start,
#'   end, strand, annot, count, strand_annot, junc_y, junc_lw,
#'   junc_tooltip, jid` as produced by `attach_junction_positions`.
#'   NULL or zero-row = no junction layer drawn.
#' @keywords internal
build_main_plot <- function(plot_data, exon_highlights,
                            junctions = NULL,
                            bed_highlights = data.table::data.table(
                              start = numeric(0), end = numeric(0)
                            ),
                            bed_color = "#B22222",
                            chr, w_start, w_end,
                            overlap_factor, color_var, summary_type) {
  lc <- resolve_line_colors(plot_data, color_var)
  plot_data <- lc$plot_data
  
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_rect(
      data = exon_highlights,
      ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
                   ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey85", alpha = 0.5
    ) +
    ggiraph::geom_rect_interactive(
      data = bed_highlights,
      ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
                   ymin = -Inf, ymax = Inf,
                   tooltip = bed_tooltip,
                   data_id = bed_tooltip),
      inherit.aes = FALSE,
      fill  = bed_color %||% "#B22222",
      alpha = 1
    ) +
    ggiraph::geom_step_interactive(
      ggplot2::aes(x = plot_x, y = offset_y, group = sample,
                   tooltip = tooltip_text, data_id = sample,
                   color = line_color),
      direction = "mid",
      linewidth = 0.4
    ) +
    ggplot2::scale_y_continuous(
      breaks = NULL,
      expand = ggplot2::expansion(add = c(0.1, overlap_factor - 1))
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::coord_cartesian(xlim = c(w_start, w_end)) +
    theme_panel_only() +
    ggplot2::labs(
      title = sprintf("Region: %s:%s-%s", chr,
                      format(w_start, big.mark = ","),
                      format(w_end,   big.mark = ",")),
      x = "Genomic Position"
    ) +
    ggforce::facet_col(ggplot2::vars(combined_facet),
                       scales = "free_y", space = "free", shrink = TRUE)
  
  # Wiggle color scale -- colors are applied but the legend is always
  # suppressed (color_var still drives the line colors)
  p <- p + ggplot2::scale_color_manual(values = lc$colors, guide = "none")
  
  # ---- Optional junction layer ---------------------------------------------
  # Open a fresh color scale via ggnewscale so the junction
  # strand/annot palette doesn't collide with the wiggle's color_var.
  # Straight horizontal segments, color = strand_annot, thickness =
  # log10(count). Each junction is its own interactive element.
  if (!is.null(junctions) && nrow(junctions) > 0) {
    pal <- junction_palette()
    p <- p +
      ggnewscale::new_scale_color() +
      ggiraph::geom_segment_interactive(
        data = junctions,
        ggplot2::aes(x = start, xend = end,
                     y = junc_y, yend = junc_y,
                     color    = strand_annot,
                     linewidth = junc_lw,
                     tooltip  = junc_tooltip,
                     data_id  = jid),
        inherit.aes = FALSE,
        lineend = "round"
      ) +
      ggplot2::scale_color_manual(
        values = pal,
        breaks = names(pal),       # stable legend order
        name   = "Junction (strand/annot)",
        drop   = TRUE
      ) +
      ggplot2::scale_linewidth_identity() +
      ggplot2::guides(color = ggplot2::guide_legend(
        nrow = 1,
        override.aes = list(linewidth = 2)
      )) +
      ggplot2::theme(legend.position = "bottom")
  }
  
  p
}
# Null-coalescing operator. Useful for the rare case where bed_color or
# bed_alpha arrives as NULL (e.g. if the package is invoked from a
# script without the BED UI controls populated).
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @keywords internal
resolve_line_colors <- function(plot_data, color_var) {
  if (!is.null(color_var) && color_var %in% colnames(plot_data)) {
    vals <- as.character(plot_data[[color_var]])
    vals[is.na(vals)] <- "NA"
    plot_data[, line_color := vals]
    list(plot_data = plot_data,
         colors    = cat_palette(sort(unique(plot_data$line_color))))
  } else {
    plot_data[, line_color := "_default"]
    list(plot_data = plot_data, colors = c("_default" = "grey15"))
  }
}