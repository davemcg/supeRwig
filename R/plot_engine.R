#' Generate Interactive Combined Region Plot (Coverage & Minimap) and Extract Data
#'
#' @param ctx Application context list returned by `build_context()`.
#' @param chr Chromosome name (e.g., "chr1").
#' @param start Genomic start coordinate.
#' @param end Genomic end coordinate.
#' @param facet_cols Character vector of metadata columns to group/facet by.
#' @param color_by Character string specifying metadata column for line colors.
#' @param target_gene Optional gene symbol used for CPM filtering.
#' @param max_samples Maximum samples per facet group (0 for unlimited).
#' @param min_expr Minimum log2(CPM+1) threshold for sample filtering.
#' @param bin_size Bin resolution in base pairs (0 for auto-calculation).
#' @param show_junctions Logical, whether to render splice junction tracks.
#' @param min_psi5 Minimum 5' PSI threshold for junctions (0 to 1).
#' @param min_psi3 Minimum 3' PSI threshold for junctions (0 to 1).
#' @param overlap_factor Height scaling factor for overlapping wiggle tracks.
#' @param output_file Output HTML filename for the combined coverage and minimap plot.
#' @param return_data Logical; if TRUE, includes underlying data.tables in return list.
#'
#' @return A named list containing `$plot` (combined girafe widget)
#'   and optionally `$data` if `return_data = TRUE`.
#' @export
plot_wiggle <- function(ctx,
                        chr,
                        start,
                        end,
                        facet_cols     = NULL,
                        color_by       = NULL,
                        target_gene    = NULL,
                        max_samples    = 3,
                        min_expr       = 5,
                        bin_size       = 0,
                        show_junctions = TRUE,
                        min_psi5       = 1,
                        min_psi3       = 1,
                        overlap_factor = 1.2,
                        output_file    = NULL,
                        return_data    = FALSE) {

  # 1. Filter Metadata & CPM Filter
  meta <- sanitize_metadata(downsample_by_study(ctx$eiad_meta, facet_cols, max_samples))
  samples <- cpm_filter_samples(meta$sample_accession, target_gene, ctx, min_expr)
  if (length(samples) == 0) stop("No samples passed expression cutoff.")
  meta <- meta[sample_accession %in% samples]

  # 2. Read BigWigs
  bp_wide <- end - start
  b_size  <- if (bin_size > 0) bin_size else as.integer(max(1, bp_wide / 750))
  n_bins  <- as.integer(max(1, ceiling(bp_wide / b_size)))
  bw_dt   <- read_region_bigwigs(samples, ctx$bw_file_map, chr, start, end, n_bins, ctx$bp_backend)

  # 3. Read Splice Junctions
  junc_band <- if (show_junctions && !is.null(ctx$sj_se)) 0.45 else 0
  raw_junc  <- if (junc_band > 0) {
    read_region_junctions(ctx, chr, start, end, samples, min_psi5 = min_psi5, min_psi3 = min_psi3)
  } else NULL

  # 4. Format Plotting Datasets
  pd_res   <- build_plot_data(bw_dt, meta, facet_cols, overlap_factor, junc_band = junc_band)
  anno_res <- subset_region_annotation(ctx$anno_dt, chr, start, end)

  junc_data <- if (!is.null(raw_junc) && nrow(raw_junc) > 0) {
    attach_junction_positions(raw_junc, pd_res$unique_samples, junc_band = 0.45)
  } else NULL

  # 5. Build ggplot Objects
  main_p <- build_main_plot(
    plot_data       = pd_res$plot_data,
    exon_highlights = anno_res$exon_highlights,
    junctions       = junc_data,
    chr             = chr,
    w_start         = start,
    w_end           = end,
    overlap_factor  = overlap_factor,
    color_var       = color_by,
    bed_highlights  = data.table::data.table(
      start = numeric(0), end = numeric(0))
  )

  mm_res <- build_minimap(
    anno_res$region_anno, anno_res$tx_base, anno_res$tx_exons,
    anno_res$has_transcripts, chr, start, end
  )

  # 6. Calculate Independent Pixel Dimensions
  dims <- compute_plot_dimensions(
    n_samples        = nrow(pd_res$unique_samples),
    n_facets         = length(unique(pd_res$unique_samples$combined_facet)),
    n_tx             = mm_res$num_tx,
    has_color        = !is.null(color_by),
    minimap_override = NA,
    show_junctions   = (junc_band > 0)
  )

  # 7. Combine Coverage and Minimap Track Objects
  combined_p <- patchwork::wrap_plots(
    main_p,
    mm_res$plot,
    ncol = 1,
    heights = c(dims$main_px, dims$minimap_px)
  )

  # 8. Render Unified Widget
  total_height_px <- dims$main_px + dims$minimap_px
  combined_widget <- ggiraph::girafe(
    ggobj      = combined_p,
    width_svg  = 14,
    height_svg = total_height_px / 72,
    options    = list(
      ggiraph::opts_sizing(rescale = TRUE, width = 1),
      ggiraph::opts_toolbar(saveaspng = FALSE, hidden = c("lasso_select", "lasso_deselect")),
      ggiraph::opts_selection(type = "none"),
      ggiraph::opts_tooltip(
        css = paste0("background-color: rgba(255,255,255,0.95); ",
                     "color: black; padding: 10px; border-radius: 5px; ",
                     "box-shadow: 2px 2px 5px rgba(0,0,0,0.2); ",
                     "font-family: Arial, sans-serif;"),
        use_fill = FALSE
      ),
      ggiraph::opts_hover(css = "stroke-width: 4px; stroke: #FF6700;")
    )
  )

  # 9. Save Single HTML File (if path requested)
  if (!is.null(output_file)) {
    htmlwidgets::saveWidget(combined_widget, file = output_file, selfcontained = TRUE)
    message(sprintf("Wiggle plot saved to: %s", output_file))
  }

  # 10. Assemble Output Bundle
  out <- list(
    plot = combined_widget
  )

  if (isTRUE(return_data)) {
    out$data <- list(
      coverage_data   = pd_res$plot_data,
      unique_samples  = pd_res$unique_samples,
      annotation      = anno_res,
      junctions       = junc_data,
      sample_metadata = meta,
      coordinates     = list(chr = chr, start = start, end = end)
    )
  }

  out
}


#' Build the main wiggle ggplot
#' @keywords internal
build_main_plot <- function(plot_data,
                            exon_highlights,
                            junctions      = NULL,
                            bed_highlights = data.table::data.table(start = numeric(0), end = numeric(0)),
                            bed_tooltip    = character(0),
                            bed_color      = "#B22222",
                            chr            = NULL,
                            w_start        = NULL,
                            w_end          = NULL,
                            overlap_factor = 1.2,
                            color_var      = NULL) {
  lc <- resolve_line_colors(plot_data, color_var)
  plot_data <- lc$plot_data

  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_rect(data = exon_highlights, ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey85", alpha = 0.5) +
    ggiraph::geom_rect_interactive(data = bed_highlights, ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf, tooltip = bed_tooltip, data_id = bed_tooltip), inherit.aes = FALSE, fill = bed_color %||% "#B22222", alpha = 1) +
    ggiraph::geom_step_interactive(ggplot2::aes(x = plot_x, y = offset_y, group = sample, tooltip = tooltip_text, data_id = sample, color = line_color), direction = "mid", linewidth = 0.4) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(add = c(0.1, overlap_factor - 1))) +
    ggplot2::scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::scale_color_manual(values = lc$colors, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(w_start, w_end)) +
    theme_panel_only() +
    ggplot2::labs(title = sprintf("Region: %s:%s-%s", chr, format(w_start, big.mark = ","), format(w_end, big.mark = ",")), x = "Genomic Position") +
    ggforce::facet_col(ggplot2::vars(combined_facet), scales = "free_y", space = "free", shrink = TRUE)

  # ---- Optional junction layer ---------------------------------------------
  if (!is.null(junctions) && nrow(junctions) > 0) {
    pal <- junction_palette()
    p <- p +
      ggnewscale::new_scale_color() +
      ggiraph::geom_segment_interactive(data = junctions, ggplot2::aes(x = start, xend = end, y = junc_y, yend = junc_y, color = strand_annot, linewidth = junc_lw, tooltip = junc_tooltip, data_id = jid), inherit.aes = FALSE, lineend = "round") +
      ggplot2::scale_color_manual(values = pal, breaks = names(pal), name = "Junction (strand/annot)", drop = TRUE) +
      ggplot2::scale_linewidth_identity() +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 1, override.aes = list(linewidth = 2))) +
      ggplot2::theme(legend.position = "bottom")
  }

  p
}

#' Build the reactive graph that feeds the plot outputs
#' @keywords internal
build_plot_reactive <- function(input, ctx, rv, bed_data, timings_rv) {

  # ---- Gated snapshot ----------------------------------------------------
  gated_r <- shiny::eventReactive(rv$trigger, {
    shiny::req(rv$trigger > 0, rv$chr, rv$start, rv$end)
    if (rv$start >= rv$end)
      stop("Start position must be less than End position.")
    list(
      chr             = rv$chr,
      start           = as.integer(rv$start),
      end             = as.integer(rv$end),
      facet_group     = input$facet_group,
      mean_average    = isTRUE(input$mean_average),
      groupings       = input$groupings,
      dynamic_filters = setNames(
        lapply(input$groupings,
               function(g) input[[paste0("dynamic_filter_", g)]]),
        input$groupings
      ),
      max_samples     = input$max_samples,
      target_gene     = input$target_gene,
      min_expr        = input$min_expr,
      bin_size        = input$bin_size,
      show_junctions  = !is.null(ctx$sj_se) && isTRUE(input$show_junctions),
      min_psi5        = as.numeric(input$min_psi5 %||% 0.0),
      min_psi3        = as.numeric(input$min_psi3 %||% 0.0)
    )
  })

  # ---- Gated nodes -------------------------------------------------------
  filtered_meta_r <- shiny::reactive({
    g <- gated_r()
    meta <- data.table::copy(ctx$eiad_meta)
    for (name in g$groupings) {
      sel <- g$dynamic_filters[[name]]
      if (length(sel) > 0) {
        meta <- if ("NA" %in% sel) {
          meta[is.na(get(name)) | get(name) %in% sel[sel != "NA"]]
        } else {
          meta[get(name) %in% sel]
        }
      }
    }
    sanitize_metadata(downsample_by_study(meta, g$facet_group, g$max_samples))
  })

  cpm_samples_r <- shiny::reactive({
    g       <- gated_r()
    meta    <- filtered_meta_r()
    samples <- cpm_filter_samples(meta$sample_accession,
                                  g$target_gene, ctx, g$min_expr)
    if (length(samples) == 0) stop("No samples passed expression cutoff.")
    list(meta = meta[sample_accession %in% samples], samples = samples)
  })

  bigwig_r <- shiny::reactive({
    g  <- gated_r()
    cs <- cpm_samples_r()
    bp_wide <- g$end - g$start
    b_size  <- if (g$bin_size > 0) g$bin_size
    else as.integer(max(1, bp_wide / 750))
    n_bins  <- as.integer(max(1, ceiling(bp_wide / b_size)))
    read_region_bigwigs(cs$samples, ctx$bw_file_map,
                        g$chr, g$start, g$end,
                        n_bins, ctx$bp_backend)
  }) |> shiny::bindCache(
    cpm_samples_r()$samples,
    gated_r()$chr, gated_r()$start, gated_r()$end,
    gated_r()$bin_size
  )

  junctions_raw_r <- shiny::reactive({
    g <- gated_r()
    if (!g$show_junctions) return(NULL)
    cs <- cpm_samples_r()
    read_region_junctions(ctx, g$chr, g$start, g$end,
                          samples = cs$samples,
                          min_psi5 = g$min_psi5,
                          min_psi3 = g$min_psi3)
  }) |> shiny::bindCache(
    gated_r()$show_junctions,
    cpm_samples_r()$samples,
    gated_r()$chr, gated_r()$start, gated_r()$end,
    gated_r()$min_psi5,
    gated_r()$min_psi3
  )

  annotation_r <- shiny::reactive({
    g <- gated_r()
    subset_region_annotation(ctx$anno_dt, g$chr, g$start, g$end)
  })

  # ---- Live nodes --------------------------------------------------------
  bed_in_region_r <- shiny::reactive({
    g       <- gated_r()
    bed_sub <- bed_data()
    if (is.null(bed_sub) || nrow(bed_sub) == 0) {
      return(data.table::data.table(start = numeric(0), end = numeric(0),
                                    bed_tooltip = character(0)))
    }
    build_bed_tooltips(
      bed_sub[seqnames == g$chr & start <= g$end & end >= g$start]
    )
  })

  plot_data_r <- shiny::reactive({
    g  <- gated_r()
    cs <- cpm_samples_r()
    junc_band <- if (g$show_junctions) 0.45 else 0

    t_bw <- Sys.time()
    bw   <- bigwig_r()
    timings_rv$bigwig <- as.numeric(Sys.time() - t_bw)

    t_pd <- Sys.time()
    out_sample <- build_plot_data(bw, cs$meta, g$facet_group,
                                  input$overlap_factor, junc_band = junc_band,
                                  mean_average = FALSE)

    out_wiggle <- if (g$mean_average) {
      build_plot_data(bw, cs$meta, g$facet_group,
                      input$overlap_factor, junc_band = junc_band,
                      mean_average = TRUE)
    } else {
      out_sample
    }

    timings_rv$plot_data <- as.numeric(Sys.time() - t_pd)
    list(sample_pd = out_sample, wiggle_pd = out_wiggle)
  })

  junctions_positioned_r <- shiny::reactive({
    g <- gated_r()
    if (!g$show_junctions) return(NULL)
    raw <- junctions_raw_r()
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    cs <- cpm_samples_r()
    pd <- plot_data_r()
    attach_junction_positions(raw, pd$wiggle_pd$unique_samples,
                              junc_band = 0.45,
                              mean_average = g$mean_average,
                              meta = cs$meta,
                              facet_cols = g$facet_group)
  })

  minimap_r <- shiny::reactive({
    g <- gated_r()
    a <- annotation_r()
    build_minimap(a$region_anno, a$tx_base, a$tx_exons, a$has_transcripts,
                  g$chr, g$start, g$end)
  })

  dimensions_r <- shiny::reactive({
    pd <- plot_data_r()$wiggle_pd
    mm <- minimap_r()
    compute_plot_dimensions(
      nrow(pd$unique_samples),
      length(unique(pd$unique_samples$combined_facet)),
      mm$num_tx,
      isTRUE(nzchar(input$color_by)),
      input$minimap_height,
      gated_r()$show_junctions
    )
  })

  heatmap_dimensions_r <- shiny::reactive({
    pd <- plot_data_r()$sample_pd
    compute_heatmap_dimensions(
      nrow(pd$unique_samples),
      length(unique(pd$unique_samples$combined_facet))
    )
  })

  # Final bundle consumed by register_outputs
  shiny::reactive({
    g         <- gated_r()
    pd        <- plot_data_r()
    cs        <- cpm_samples_r()
    bw        <- bigwig_r()
    a         <- annotation_r()
    mm        <- minimap_r()
    dims      <- dimensions_r()
    color_var <- if (isTRUE(nzchar(input$color_by))) input$color_by else NULL

    t_mp <- Sys.time()

    # 1. Build Wiggle ggplot (uses wiggle_pd which respects mean_average)
    main_plot <- build_main_plot(
      plot_data       = pd$wiggle_pd$plot_data,
      exon_highlights = a$exon_highlights,
      junctions       = junctions_positioned_r(),
      bed_highlights  = bed_in_region_r(),
      bed_color       = input$bed_color,
      chr             = g$chr,
      w_start         = g$start,
      w_end           = g$end,
      overlap_factor  = input$overlap_factor,
      color_var       = color_var
    )

    # 2. Build heatmap ggplot (always uses sample_pd for individual sample rows)
    heatmap_plot <- build_heatmap_ggplot(
      plot_data  = pd$sample_pd$plot_data,
      scale_rows = isTRUE(input$scale_rows_hm),
      chr        = g$chr,
      w_start    = g$start,
      w_end      = g$end
    )

    timings_rv$main_plot <- as.numeric(Sys.time() - t_mp)

    g_strand <- "+"
    if (!is.null(g$target_gene) && !is.null(ctx$anno_dt)) {
      df_sub <- ctx$anno_dt[type == "gene" & gene_name == g$target_gene, strand]
      if (length(df_sub) > 0) g_strand <- as.character(df_sub[1])
    }
    gene_oriented_label <- if (g_strand == "-") {
      paste0(g$target_gene, " (3' ← 5')")
    } else {
      paste0(g$target_gene, " (5' → 3')")
    }

    list(
      plot            = main_plot,
      heatmap_plot    = heatmap_plot,
      facet_cols      = g$facet_group,
      target_gene     = gene_oriented_label,
      auto_height     = dims$main_px,
      heatmap_height  = heatmap_dimensions_r()$main_px,
      plot_minimap    = mm$plot,
      minimap_px      = dims$minimap_px,
      tx_hover_info   = mm$tx_hover_info,
      has_transcripts = a$has_transcripts,
      gated_params    = g,
      plot_data       = pd$wiggle_pd$plot_data,
      exon_highlights = a$exon_highlights,
      junctions       = junctions_positioned_r(),
      bed_highlights  = bed_in_region_r(),
      color_var       = color_var
    )
  })
}


#' @keywords internal
resolve_line_colors <- function(plot_data, color_var) {
  if (!is.null(color_var) && color_var %in% colnames(plot_data)) {
    plot_data[, line_color := data.table::fcoalesce(as.character(get(color_var)), "NA")]
    list(plot_data = plot_data, colors = cat_palette(sort(unique(plot_data$line_color))))
  } else {
    plot_data[, line_color := "_default"]
    list(plot_data = plot_data, colors = c("_default" = "grey15"))
  }
}


#' Build the transcript/exon minimap
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
                   y = tx_idx, yend = tx_idx,
                   color = is_principal),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.08, "inches"),
                             type = "closed")
    ) +
    ggplot2::geom_rect(
      data = tx_exons,
      ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
                   ymin = tx_idx - 0.25, ymax = tx_idx + 0.25,
                   color = is_principal, fill = is_principal)
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
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "gray25"), guide = "none") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "firebrick1", "FALSE" = "grey25"), guide = "none")

  meta_cols <- intersect(
    c("tag", "transcript_type"),
    colnames(tx_base)
  )
  target_cols <- c("tx_label", "seqnames", "start", "end", "strand", "tx_idx", meta_cols)

  list(
    plot          = p,
    num_tx        = length(unique(region_anno$tx_idx)),
    tx_hover_info = tx_base[, target_cols, with = FALSE]
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
#' @keywords internal
compute_plot_dimensions <- function(n_samples, n_facets, n_tx,
                                    has_color, minimap_override,
                                    show_junctions = FALSE) {
  legend_px  <- if (has_color)      60 else 0
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


#' Pixel height for the raster heatmap
#' @keywords internal
compute_heatmap_dimensions <- function(n_samples, n_facets) {
  strip_px  <- n_facets * 15    # top facet strips (one per group)
  rows_px   <- n_samples * 6    # compact raster bands
  chrome_px <- 100              # title + bottom legend + axis
  list(main_px = max(300L, as.integer(rows_px + strip_px + chrome_px)))
}


#' Subset annotation to features overlapping the window
#' @keywords internal
subset_region_annotation <- function(anno_dt, chr, w_start, w_end) {
  region_anno <- anno_dt[seqnames == chr & start <= w_end & end >= w_start & type %in% c("transcript", "exon")]
  has_tx <- nrow(region_anno) > 0

  if (has_tx) {
    region_anno[, `:=`(tx_label = paste0(gene_name, " - ", transcript_id),
                       tx_idx = as.numeric(as.factor(paste0(gene_name, " - ", transcript_id))))]

    principal_ids <- character(0)
    if ("tag" %in% colnames(region_anno)) {
      principal_ids <- unique(region_anno[grepl("GENCODE_Primary|CCDS", tag), transcript_id])
    } else if ("appris" %in% colnames(region_anno)) {
      principal_ids <- unique(region_anno[grepl("appris_principal_1", appris), transcript_id])
    }

    region_anno[, is_principal := transcript_id %in% principal_ids]

    exon_hi <- data.table::as.data.table(GenomicRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(region_anno[type == "exon"])))
  } else {
    exon_hi <- tx_base <- tx_exons <- data.table::data.table()
  }

  list(
    region_anno = region_anno,
    tx_base = region_anno[type == "transcript"],
    tx_exons = region_anno[type == "exon"],
    exon_highlights = exon_hi,
    has_transcripts = has_tx
  )
}

#' @keywords internal
.tooltip_safe_cols <- function(meta_cur, max_median_chars = 120) {
  candidate <- setdiff(colnames(meta_cur), c("sample_accession", "dummy_facet", "combined_facet"))
  too_long <- vapply(candidate, function(c) {
    v <- meta_cur[[c]]
    if (!is.character(v) && !is.factor(v)) return(FALSE)
    stats::median(nchar(as.character(v), type = "bytes"), na.rm = TRUE) > max_median_chars
  }, logical(1))
  candidate[!too_long]
}

#' @keywords internal
build_plot_data <- function(dt_full, meta_cur, facet_cols, overlap_factor, junc_band = 0, tooltip_max_chars = 120, mean_average = FALSE) {
  if (length(facet_cols) == 0) { meta_cur$dummy_facet <- "All Samples"; facet_cols <- "dummy_facet" }

  meta_work <- data.table::copy(meta_cur)
  tissue_map <- unique(meta_work)[!duplicated(sample_accession)]
  tissue_map[, combined_facet := do.call(paste, c(.SD, sep = " - ")),
             .SDcols = facet_cols]

  if (isTRUE(mean_average)) {
    pd_raw <- merge(dt_full, tissue_map, by.x = "sample", by.y = "sample_accession", all.x = TRUE)

    pd <- pd_raw[, .(
      value     = mean(value, na.rm = TRUE),
      n_samples = data.table::uniqueN(sample)
    ), by = .(combined_facet, binned_pos, bin_end)]

    facet_meta <- tissue_map[, .(
      n_samples_total = .N
    ), by = combined_facet]

    safe_cols <- .tooltip_safe_cols(tissue_map, tooltip_max_chars)
    for (col in safe_cols) {
      if (col %in% colnames(tissue_map)) {
        val_dt <- tissue_map[, .(
          col_val = if (data.table::uniqueN(get(col)) == 1) as.character(get(col)[1]) else paste0(data.table::uniqueN(get(col)), " values")
        ), by = combined_facet]
        data.table::setnames(val_dt, "col_val", col)
        facet_meta <- merge(facet_meta, val_dt, by = "combined_facet", all.x = TRUE)
      }
    }

    tips <- paste0("<b>Group / Facet:</b> ", facet_meta$combined_facet,
                   "<br><b>Samples:</b> ", facet_meta$n_samples_total)
    for (col in setdiff(safe_cols, c("combined_facet", "dummy_facet"))) {
      s <- as.character(facet_meta[[col]])
      long <- !is.na(s) & nchar(s, type = "bytes") > tooltip_max_chars
      s[long] <- paste0(substr(s[long], 1, tooltip_max_chars - 1), "\u2026")
      tips <- paste0(tips, "<br><b>", col, ":</b> ", s)
    }
    facet_meta[, static_tooltip := tips]
    facet_meta[, sample := combined_facet]

    pd <- merge(pd, facet_meta, by = "combined_facet", all.x = TRUE)

    unique_samples <- unique(pd[, .(sample, combined_facet)])[order(combined_facet, sample)][, local_idx := seq_len(.N), by = combined_facet]
    pd <- merge(pd, unique_samples, by = c("sample", "combined_facet"))

    pd[, `:=`(log_val = log2(value + 1), plot_x = (binned_pos + bin_end) / 2)]
    max_log <- max(pd$log_val, na.rm = TRUE); if (max_log == 0 || is.na(max_log)) max_log <- 1

    pd[, offset_y := (log_val / max_log) * overlap_factor + local_idx + junc_band]
    pd[, tooltip_text := paste0(static_tooltip, "<br><b>Mean Coverage:</b> ", round(value, 2))]
    data.table::setorderv(pd, c("combined_facet", "local_idx", "plot_x"))

    list(plot_data = pd, unique_samples = unique_samples)
  } else {
    tips <- paste0("<b>Sample:</b> ", tissue_map$sample_accession)
    for (col in .tooltip_safe_cols(tissue_map, tooltip_max_chars)) {
      s <- as.character(tissue_map[[col]])
      long <- !is.na(s) & nchar(s, type = "bytes") > tooltip_max_chars
      s[long] <- paste0(substr(s[long], 1, tooltip_max_chars - 1), "\u2026")
      tips <- paste0(tips, "<br><b>", col, ":</b> ", s)
    }
    tissue_map[, static_tooltip := tips]

    pd <- merge(dt_full, tissue_map, by.x = "sample", by.y = "sample_accession", all.x = TRUE)
    unique_samples <- unique(pd[, .(sample, combined_facet)])[order(combined_facet, sample)][, local_idx := seq_len(.N), by = combined_facet]

    pd <- merge(pd, unique_samples, by = c("sample", "combined_facet"))
    pd[, `:=`(log_val = log2(value + 1), plot_x = (binned_pos + bin_end) / 2)]
    max_log <- max(pd$log_val, na.rm = TRUE); if (max_log == 0 || is.na(max_log)) max_log <- 1

    pd[, offset_y := (log_val / max_log) * overlap_factor + local_idx + junc_band]
    pd[, tooltip_text := paste0(static_tooltip, ": ", round(value, 2))]
    data.table::setorderv(pd, c("combined_facet", "local_idx", "plot_x"))

    list(plot_data = pd, unique_samples = unique_samples)
  }
}

#' @keywords internal
build_bed_tooltips <- function(bed_dt) {
  if (nrow(bed_dt) == 0) return(bed_dt[, bed_tooltip := character(0)])
  tip <- sprintf("<b>%s:</b> %s-%s (%.2f kb)", bed_dt$seqnames, format(bed_dt$start, big.mark = ","), format(bed_dt$end, big.mark = ","), (bed_dt$end - bed_dt$start + 1) / 1000)
  if ("name" %in% colnames(bed_dt)) tip <- paste0(tip, "<br><b>Name:</b> ", bed_dt$name)
  if ("score" %in% colnames(bed_dt)) tip <- paste0(tip, "<br><b>Score:</b> ", bed_dt$score)
  if ("strand" %in% colnames(bed_dt)) tip <- paste0(tip, "<br><b>Strand:</b> ", bed_dt$strand)
  bed_dt[, bed_tooltip := tip]
}


#' Strand x annotation color palette for the junction layer
#' @keywords internal
junction_palette <- function() {
  c(
    "+/annot" = "royalblue4",
    "+/novel" = "royalblue1",
    "-/annot" = "tomato4",
    "-/novel" = "tomato1",
    "*/annot" = "seagreen4",
    "*/novel" = "seagreen1"
  )
}

#' @keywords internal
magma_palette <- function(n = 256) {
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    viridisLite::magma(n)
  } else if (requireNamespace("viridis", quietly = TRUE)) {
    viridis::magma(n)
  } else {
    grDevices::hcl.colors(n, palette = "Magma")
  }
}

#' Attach y-positions and visual attrs to junctions with Interval Packing
#' @keywords internal
attach_junction_positions <- function(junctions, unique_samples,
                                      junc_band = 0.65,
                                      mean_average = FALSE,
                                      meta = NULL,
                                      facet_cols = NULL) {
  if (nrow(junctions) == 0) return(junctions)

  if (isTRUE(mean_average) && !is.null(meta)) {
    s_map <- unique(meta)[!duplicated(sample_accession)]
    if (length(facet_cols) == 0) { s_map$dummy_facet <- "All Samples"; facet_cols <- "dummy_facet" }
    s_map[, combined_facet := do.call(paste, c(.SD, sep = " - ")), .SDcols = facet_cols]

    junc_mapped <- merge(junctions, s_map[, .(sample_accession, combined_facet)],
                         by.x = "sample", by.y = "sample_accession")
    if (nrow(junc_mapped) == 0) return(junctions[0])

    extra_cols <- intersect(c("cluster_5", "cluster_3"), colnames(junc_mapped))
    by_cols <- c("combined_facet", "jid", "start", "end", "strand", "annot", "strand_annot", extra_cols)

    junc_avg <- junc_mapped[, .(
      psi5      = mean(psi5, na.rm = TRUE),
      psi3      = mean(psi3, na.rm = TRUE),
      raw_count = if (all(is.na(raw_count))) NA_real_ else mean(raw_count, na.rm = TRUE),
      n_samples = .N
    ), by = by_cols]

    junc_avg[, count := pmax(psi5, psi3)]
    junc_avg[, sample := combined_facet]

    junc <- merge(junc_avg, unique_samples, by = c("sample", "combined_facet"))
  } else {
    junc <- merge(junctions, unique_samples, by = "sample")
  }

  if (nrow(junc) == 0) return(junc)

  pack_lanes <- function(start_vec, end_vec) {
    if (length(start_vec) == 0) return(integer(0))

    ord <- order(start_vec)
    s_sorted <- start_vec[ord]
    e_sorted <- end_vec[ord]

    lane_ends <- numeric(0)
    assigned_lanes <- integer(length(start_vec))
    buffer <- 200L

    for (i in seq_along(s_sorted)) {
      placed <- FALSE
      for (l in seq_along(lane_ends)) {
        if (s_sorted[i] > (lane_ends[l] + buffer)) {
          lane_ends[l] <- e_sorted[i]
          assigned_lanes[ord[i]] <- l - 1L
          placed <- TRUE
          break
        }
      }
      if (!placed) {
        lane_ends <- c(lane_ends, e_sorted[i])
        assigned_lanes[ord[i]] <- length(lane_ends) - 1L
      }
    }
    return(assigned_lanes)
  }

  junc[, sub_idx := pack_lanes(start, end), by = .(combined_facet, local_idx)]

  sub_spacing <- 0.060
  n_visible   <- max(1L, as.integer(floor((junc_band - 0.08) / sub_spacing)))
  junc[, sub_idx := sub_idx %% n_visible]

  junc[, junc_y  := local_idx + junc_band - 0.07 - (sub_idx * sub_spacing)]
  junc[, junc_lw := 0.8]

  if (isTRUE(mean_average)) {
    junc[, junc_tooltip := paste0(
      "<b>Junction:</b> ", jid,
      "<br><b>Group / Facet:</b> ", combined_facet,
      "<br><b>Mean PSI5 (5' Donor Focus):</b> ", round(psi5, 2), "%",
      "<br><b>Mean PSI3 (3' Acceptor Focus):</b> ", round(psi3, 2), "%",
      "<br><b>Samples with Junction:</b> ", n_samples
    )]
    if ("raw_count" %in% colnames(junc) && !all(is.na(junc$raw_count))) {
      junc[, junc_tooltip := paste0(junc_tooltip, "<br><b>Mean Raw Count:</b> ", round(raw_count, 1))]
    }
  } else {
    junc[, junc_tooltip := paste0(
      "<b>Junction:</b> ", jid,
      "<br><b>Sample:</b> ",  sample,
      "<br><b>PSI5 (5' Donor Focus):</b> ", round(psi5, 2), "%",
      "<br><b>PSI3 (3' Acceptor Focus):</b> ", round(psi3, 2), "%"
    )]
    if ("raw_count" %in% colnames(junc)) {
      junc[, junc_tooltip := paste0(junc_tooltip, "<br><b>Raw Count:</b> ",
                                    data.table::fifelse(is.na(raw_count), "N/A", as.character(raw_count)))]
    }
  }

  if ("cluster_5" %in% colnames(junc)) {
    junc[, junc_tooltip := paste0(junc_tooltip, "<br><b>5' Cluster (Donor):</b> ", cluster_5)]
  }
  if ("cluster_3" %in% colnames(junc)) {
    junc[, junc_tooltip := paste0(junc_tooltip, "<br><b>3' Cluster (Acceptor):</b> ", cluster_3)]
  }

  junc[, junc_tooltip := paste0(
    junc_tooltip,
    "<br><b>Strand:</b> ",  strand,
    "<br><b>Status:</b> ",  ifelse(annot == 1L, "annotated", "novel")
  )]

  junc
}


#' Generate Heatmap Visualization of Read Coverage with Aligned Transcripts
#'
#' @export
plot_region_heatmap <- function(ctx,
                                chr,
                                start,
                                end,
                                facet_cols       = NULL,
                                target_gene      = NULL,
                                max_samples      = 3,
                                min_expr         = 5,
                                bin_size         = 0,
                                scale_rows       = TRUE,
                                palette          = "magma",
                                show_transcripts = TRUE,
                                show_row_names   = FALSE,
                                show_column_axis = TRUE,
                                title            = NULL,
                                return_data      = FALSE) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required for heatmap visualization.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required for heatmap color mapping.")
  }

  # 1. Filter Metadata & CPM Filter
  meta <- sanitize_metadata(downsample_by_study(ctx$eiad_meta, facet_cols, max_samples))
  samples <- cpm_filter_samples(meta$sample_accession, target_gene, ctx, min_expr)
  if (length(samples) == 0) stop("No samples passed expression cutoff.")
  meta <- meta[sample_accession %in% samples]

  # 2. Read BigWigs
  bp_wide <- end - start
  b_size  <- if (bin_size > 0) bin_size else as.integer(max(1, bp_wide / 750))
  n_bins  <- as.integer(max(1, ceiling(bp_wide / b_size)))
  bw_dt   <- read_region_bigwigs(samples, ctx$bw_file_map, chr, start, end, n_bins, ctx$bp_backend)

  # 3. Handle Faceting & Metadata Alignment
  if (is.null(facet_cols) || length(facet_cols) == 0) {
    meta$dummy_facet <- "All Samples"
    facet_cols <- "dummy_facet"
  }

  meta_unique <- unique(meta)[!duplicated(sample_accession)]
  meta_unique[, combined_facet := do.call(paste, c(.SD, sep = " - ")), .SDcols = facet_cols]

  dt_merged <- merge(bw_dt, meta_unique, by.x = "sample", by.y = "sample_accession", all.x = TRUE)

  # 4. Reshape to Wide Matrix (Always Sample-Level)
  mat_dt <- data.table::dcast(
    dt_merged,
    sample ~ binned_pos,
    value.var = "value",
    fun.aggregate = mean,
    fill = 0
  )
  mat_samples <- mat_dt$sample
  mat <- as.matrix(mat_dt[, -1, with = FALSE])
  rownames(mat) <- mat_samples
  sample_meta <- meta_unique[match(rownames(mat), sample_accession)]
  row_split   <- sample_meta$combined_facet

  # Reorder columns numerically by genomic start position
  col_positions <- as.numeric(colnames(mat))
  col_ord <- order(col_positions)
  mat <- mat[, col_ord, drop = FALSE]

  # 5. min max row scaling
  if (isTRUE(scale_rows)) {
    row_mins <- apply(mat, 1, min, na.rm = TRUE)
    row_maxs <- apply(mat, 1, max, na.rm = TRUE)
    row_range <- row_maxs - row_mins
    row_range[row_range == 0 | is.na(row_range)] <- 1
    mat <- (mat - row_mins) / row_range
  }

  # 6. Color Mapping
  magma_colors <- magma_palette(100)

  min_val <- min(mat, na.rm = TRUE)
  max_val <- max(mat, na.rm = TRUE)
  if (min_val == max_val) max_val <- min_val + 1
  breaks <- seq(min_val, max_val, length.out = 100)
  col_fun <- circlize::colorRamp2(breaks, magma_colors)

  # 8. Subset Region Annotation & Build Aligned Transcripts Annotation
  anno_res <- subset_region_annotation(ctx$anno_dt, chr, start, end)

  top_anno <- if (isTRUE(show_transcripts)) {
    build_aligned_transcript_annotation(anno_res, start, end)
  } else NULL

  bottom_anno <- if (isTRUE(show_column_axis)) {
    build_coordinate_axis_annotation(start, end)
  } else NULL

  legend_label <- "Min Max"
  plot_title <- title %||% sprintf("Region: %s:%s-%s", chr, format(start, big.mark = ","), format(end, big.mark = ","))

  # 9. Build ComplexHeatmap Object
  ht <- ComplexHeatmap::Heatmap(
    matrix               = mat,
    name                 = legend_label,
    col                  = col_fun,
    cluster_rows         = FALSE,
    cluster_columns      = FALSE,
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    row_split            = row_split,
    row_title_rot        = 0,
    row_title_gp         = grid::gpar(fontsize = 9, fontface = "bold"),
    row_names_gp         = grid::gpar(fontsize = 8),
    top_annotation       = top_anno,
    bottom_annotation    = bottom_anno,
    column_title         = plot_title,
    column_title_gp      = grid::gpar(fontsize = 11, fontface = "bold"),
    use_raster           = TRUE,
    raster_resize_mat    = FALSE,
    raster_quality       = 5
  )

  if (isTRUE(return_data)) {
    return(list(
      heatmap         = ht,
      matrix          = mat,
      sample_metadata = sample_meta,
      annotation      = anno_res,
      coordinates     = list(chr = chr, start = start, end = end)
    ))
  }

  ht
}

#' Build the coverage heatmap as a ggplot (raster), aligned to the wiggle/minimap
#' @keywords internal
build_heatmap_ggplot <- function(plot_data, scale_rows, chr, w_start, w_end, title = NULL) {
  pd <- data.table::copy(plot_data)

  if (isTRUE(scale_rows)) {
    pd[, fill_val := {
      lo  <- min(value, na.rm = TRUE)
      hi  <- max(value, na.rm = TRUE)
      rng <- hi - lo
      if (!is.finite(rng) || rng == 0) rep(0, .N) else (value - lo) / rng
    }, by = sample]
    legend_label <- "Min Max"
  } else {
    pd[, fill_val := value]
    legend_label <- "Coverage"
  }

  pd[, combined_facet := factor(
    combined_facet,
    levels = sort(unique(as.character(combined_facet)))
  )]

  sample_rows <- unique(pd[, .(sample, combined_facet, local_idx, static_tooltip)])

  plot_title <- title %||% sprintf(
    "Region: %s:%s-%s", chr,
    format(w_start, big.mark = ","), format(w_end, big.mark = ",")
  )

  ggplot2::ggplot(pd, ggplot2::aes(x = plot_x, y = local_idx)) +
    ggplot2::geom_raster(ggplot2::aes(fill = fill_val), interpolate = FALSE) +
    ggiraph::geom_rect_interactive(
      data = sample_rows,
      ggplot2::aes(
        xmin = w_start, xmax = w_end,
        ymin = local_idx - 0.5, ymax = local_idx + 0.5,
        tooltip = static_tooltip, data_id = sample
      ),
      inherit.aes = FALSE,
      fill = "white", alpha = 0.001
    ) +
    ggplot2::scale_fill_gradientn(colours = magma_palette(256), name = legend_label) +
    ggplot2::scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(0)) +
    ggplot2::coord_cartesian(xlim = c(w_start, w_end)) +
    theme_panel_only() +
    ggplot2::labs(title = plot_title, x = "Genomic Position", y = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barwidth  = grid::unit(6, "lines"),
        barheight = grid::unit(0.4, "lines"),
        title.vjust = 0.85
      )
    ) +
    ggforce::facet_col(ggplot2::vars(combined_facet), scales = "free_y", space = "free", shrink = TRUE) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.spacing.y = grid::unit(0, "pt"),
      strip.text      = ggplot2::element_text(
        hjust = 0.5,
        margin = ggplot2::margin(t = 1, r = 0, b = 1, l = 0, unit = "pt"))
    )
}

#' Build Aligned Transcript Column Annotation for ComplexHeatmap
#' @keywords internal
build_aligned_transcript_annotation <- function(anno_res, w_start, w_end) {
  if (!anno_res$has_transcripts) {
    empty_fn <- ComplexHeatmap::AnnotationFunction(
      fun = function(index, k, n) {
        grid::grid.text("No transcripts in region", x = 0.5, y = 0.5,
                        gp = grid::gpar(fontsize = 9, col = "grey50"))
      },
      height = grid::unit(0.8, "cm")
    )
    return(ComplexHeatmap::HeatmapAnnotation(
      transcripts = empty_fn,
      show_annotation_name = FALSE
    ))
  }

  tx_base   <- anno_res$tx_base
  tx_exons  <- anno_res$tx_exons
  num_tx    <- length(unique(anno_res$region_anno$tx_idx))
  height_cm <- min(6.0, max(1.2, num_tx * 0.45 + 0.3))

  anno_fn <- ComplexHeatmap::AnnotationFunction(
    fun = function(index, k, n) {
      w_span <- w_end - w_start
      if (w_span <= 0) return(NULL)

      for (i in seq_len(num_tx)) {
        y_center <- 1 - (i - 0.5) / num_tx
        box_h    <- 0.45 / num_tx

        curr_tx <- tx_base[tx_base$tx_idx == i, ]
        if (nrow(curr_tx) > 0) {
          t_start  <- curr_tx$start[1]
          t_end    <- curr_tx$end[1]
          t_strand <- curr_tx$strand[1]
          t_label  <- curr_tx$tx_label[1]
          is_princ <- isTRUE(curr_tx$is_principal[1])

          x1_npc <- pmax(0, pmin(1, (t_start - w_start) / w_span))
          x2_npc <- pmax(0, pmin(1, (t_end - w_start) / w_span))

          line_col <- if (is_princ) "firebrick3" else "grey35"

          grid::grid.lines(
            x = grid::unit(c(x1_npc, x2_npc), "npc"),
            y = grid::unit(c(y_center, y_center), "npc"),
            gp = grid::gpar(col = line_col, lwd = 1.5)
          )

          if ((x2_npc - x1_npc) > 0.05) {
            n_arrows  <- max(2, min(5, floor((x2_npc - x1_npc) * 8)))
            arrow_pos <- seq(x1_npc + 0.02, x2_npc - 0.02, length.out = n_arrows)
            dx <- if (t_strand == "-") -0.003 else 0.003
            for (ap in arrow_pos) {
              grid::grid.lines(
                x = grid::unit(c(ap - dx, ap + dx), "npc"),
                y = grid::unit(c(y_center, y_center), "npc"),
                arrow = grid::arrow(length = grid::unit(0.04, "inches"), type = "closed"),
                gp = grid::gpar(col = line_col, fill = line_col)
              )
            }
          }

          grid::grid.text(
            label = t_label,
            x = grid::unit(pmax(0.005, x1_npc), "npc"),
            y = grid::unit(y_center + box_h * 0.55, "npc"),
            just = c("left", "bottom"),
            gp = grid::gpar(fontsize = 7, col = "grey20", fontface = if (is_princ) "bold" else "plain")
          )
        }

        curr_exons <- tx_exons[tx_exons$tx_idx == i, ]
        if (nrow(curr_exons) > 0) {
          for (j in seq_len(nrow(curr_exons))) {
            ex_s <- curr_exons$start[j]
            ex_e <- curr_exons$end[j]
            is_princ <- isTRUE(curr_exons$is_principal[j])

            ex_x1 <- pmax(0, pmin(1, (ex_s - w_start) / w_span))
            ex_x2 <- pmax(0, pmin(1, (ex_e - w_start) / w_span))
            ex_w  <- ex_x2 - ex_x1

            if (ex_w > 0) {
              fill_c   <- if (is_princ) "firebrick1" else "grey35"
              border_c <- if (is_princ) "firebrick4" else "grey15"

              grid::grid.rect(
                x = grid::unit(ex_x1, "npc"),
                y = grid::unit(y_center, "npc"),
                width = grid::unit(ex_w, "npc"),
                height = grid::unit(box_h, "npc"),
                just = c("left", "center"),
                gp = grid::gpar(fill = fill_c, col = border_c, lwd = 0.8)
              )
            }
          }
        }
      }
    },
    height = grid::unit(height_cm, "cm")
  )

  ComplexHeatmap::HeatmapAnnotation(
    transcripts = anno_fn,
    show_annotation_name = FALSE
  )
}

#' Build Genomic Coordinate Axis Bottom Annotation for ComplexHeatmap
#' @keywords internal
build_coordinate_axis_annotation <- function(w_start, w_end) {
  anno_axis <- ComplexHeatmap::AnnotationFunction(
    fun = function(index, k, n) {
      ticks  <- stats::quantile(c(w_start, w_end), probs = seq(0, 1, length.out = 6))
      ticks  <- unique(round(ticks))
      at_npc <- (ticks - w_start) / (w_end - w_start)
      labels <- format(ticks, big.mark = ",", scientific = FALSE)

      grid::grid.xaxis(
        at = at_npc,
        label = labels,
        gp = grid::gpar(fontsize = 8)
      )
    },
    height = grid::unit(0.7, "cm")
  )

  ComplexHeatmap::HeatmapAnnotation(
    genomic_position = anno_axis,
    show_annotation_name = FALSE
  )
}
