#' Build the transcript/exon minimap
#'
#' @keywords internal
build_minimap <- function(region_anno, tx_base, tx_exons, has_transcripts,
                          chr, w_start, w_end) {
  if (!has_transcripts) {
    return(list(
      plot          = empty_minimap(w_start, w_end),
      num_tx        = 1,
      tx_hover_info = data.table::data.table()
    ))
  }
  
  tx_base[, draw_start := ifelse(strand == "-", end + 0.5, start - 0.5)]
  tx_base[, draw_end   := ifelse(strand == "-", start - 0.5, end + 0.5)]
  tx_base[, label_x    := pmax(w_start, pmin(start, end))]
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tx_base,
      ggplot2::aes(x = draw_start, xend = draw_end,
                   y = tx_idx, yend = tx_idx),
      color = "black",
      arrow = ggplot2::arrow(length = ggplot2::unit(0.08, "inches"),
                             type = "closed")
    ) +
    ggplot2::geom_rect(
      data = tx_exons,
      ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
                   ymin = tx_idx - 0.25, ymax = tx_idx + 0.25),
      fill = "black"
    ) +
    ggplot2::geom_text(
      data = tx_base,
      ggplot2::aes(x = label_x, y = tx_idx + 0.4, label = tx_label),
      hjust = 0, vjust = 0, size = 2.8, color = "grey25"
    ) +
    ggplot2::scale_y_continuous(
      breaks = NULL,
      expand = ggplot2::expansion(add = c(0.4, 0.9))
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::coord_cartesian(xlim = c(w_start, w_end), clip = "off") +
    theme_panel_only() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)
  
  list(
    plot          = p,
    num_tx        = length(unique(region_anno$tx_idx)),
    tx_hover_info = tx_base[, .(tx_label, seqnames, start, end, strand,
                                tx_idx)]
  )
}

#' @keywords internal
empty_minimap <- function(w_start, w_end) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = (w_start + w_end) / 2, y = 1,
                      label = "No transcripts/exons in region") +
    ggplot2::coord_cartesian(xlim = c(w_start, w_end)) +
    theme_panel_only() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)
}

#' Compute pixel heights for the main plot and minimap
#'
#' @keywords internal
compute_plot_dimensions <- function(n_samples, n_facets, n_tx,
                                    has_color, minimap_override,
                                    show_junctions = FALSE) {
  legend_px  <- if (has_color)      60 else 0
  # Junctions live inside the per-sample band; this small bump is just
  # enough to keep them legible when many samples stack up.
  junc_px    <- if (show_junctions) n_samples * 20 else 0
  main_px    <- (n_samples * 15) + (n_facets * 35) + 100 +
    legend_px + junc_px
  
  minimap_px <- if (!is.na(minimap_override) && minimap_override > 0) {
    as.integer(minimap_override)
  } else {
    min(300L, (n_tx * 18) + 50)
  }
  
  list(main_px = main_px, minimap_px = minimap_px)
}
