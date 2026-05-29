#' Build the reactive graph that feeds the plot outputs
#'
#' Returns a single reactive of the bundle that `register_outputs` expects.
#' The graph is split into "gated" nodes (re-evaluated only when the user
#' clicks Generate Plot or hits Reset Zoom) and "live" nodes (re-evaluated
#' whenever a cosmetic input changes, including the brush-driven view).
#'
#' @keywords internal
build_plot_reactive <- function(input, ctx, rv, bed_data) {
  
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
      min_junc_reads  = as.integer(input$min_junc_reads %||% 5L)
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
    else as.integer(max(1, bp_wide / 1000)) # bigger values (e.g.2000) increases the resolution of wiggle plot (and increases plotting time)
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
                          samples = cs$samples, min_reads = g$min_junc_reads)
  }) |> shiny::bindCache(
    gated_r()$show_junctions,
    cpm_samples_r()$samples,
    gated_r()$chr, gated_r()$start, gated_r()$end,
    gated_r()$min_junc_reads
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
    t0 <- Sys.time()
    g  <- gated_r()
    cs <- cpm_samples_r()
    junc_band <- if (g$show_junctions) 0.4 else 0
    out <- build_plot_data(bigwig_r(), cs$meta, g$facet_group,
                           input$overlap_factor, 
                           junc_band = junc_band)
    cat("build_plot_data:", round(as.numeric(Sys.time() - t0), 2), "s\n")
    out
  })
  
  junctions_positioned_r <- shiny::reactive({
    g <- gated_r()
    if (!g$show_junctions) return(NULL)
    raw <- junctions_raw_r()
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    attach_junction_positions(raw, plot_data_r()$unique_samples,
                              junc_band = 0.4)
  })
  
  minimap_r <- shiny::reactive({
    g <- gated_r()
    a <- annotation_r()
    build_minimap(a$region_anno, a$tx_base, a$tx_exons, a$has_transcripts,
                  g$chr, g$start, g$end)
  })
  
  dimensions_r <- shiny::reactive({
    pd <- plot_data_r()
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
  
  # Final bundle consumed by register_outputs
  shiny::reactive({
    g    <- gated_r()
    pd   <- plot_data_r()
    a    <- annotation_r()
    mm   <- minimap_r()
    dims <- dimensions_r()
    color_var <- if (isTRUE(nzchar(input$color_by))) input$color_by else NULL
    list(
      plot = build_main_plot(
        pd$plot_data, a$exon_highlights, junctions_positioned_r(),
        bed_in_region_r(), input$bed_color,
        g$chr, g$start, g$end,
        input$overlap_factor, color_var
      ),
      facet_cols   = g$facet_group,
      target_gene     = g$target_gene,
      auto_height     = dims$main_px,
      auto_height     = dims$main_px,
      plot_minimap    = mm$plot,
      minimap_px      = dims$minimap_px,
      tx_hover_info   = mm$tx_hover_info,
      has_transcripts = a$has_transcripts
    )
  })
}