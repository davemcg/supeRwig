#' Build the main wiggle ggplot
#'
#' @param plot_data        Per-bin coverage with offsets and tooltips.
#' @param exon_highlights  Reduced exon ranges to shade behind the traces.
#' @param bed_highlights   Optional user-uploaded BED ranges (subset to
#'                         the current window). Empty data.table = no
#'                         BED layer drawn.
#' @param bed_color  Fill color and opacity for BED bands.
#' @keywords internal
build_main_plot <- function(plot_data, exon_highlights,
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
    # User BED layer -- interactive so hovering a band shows its
    # coordinates / name / score. Empty bed_highlights yields zero
    # geoms, so this is free when nothing is uploaded.
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
  
  if (!is.null(color_var)) {
    p + ggplot2::scale_color_manual(values = lc$colors, name = color_var) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title    = ggplot2::element_text(size = 9),
        legend.text     = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(0.4, "cm")
      )
  } else {
    p + ggplot2::scale_color_manual(values = lc$colors, guide = "none") +
      ggplot2::theme(legend.position = "none")
  }
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