#' Build the app server function
#'
#' @keywords internal
build_server <- function(ctx) {
  function(input, output, session) {
    shiny::updateSelectizeInput(
      session, "target_gene",
      choices  = ctx$unique_genes, server = TRUE,
      selected = if ("ABCA4" %in% ctx$unique_genes) "ABCA4"
      else ctx$unique_genes[1]
    )
    
    rv <- shiny::reactiveValues(
      chr = NULL, start = NULL, end = NULL, trigger = 0,
      home_chr = NULL, home_start = NULL, home_end = NULL
    )
    
    register_navigation_observers(input, session, ctx, rv)
    register_dynamic_filters(input, output, ctx)
    
    bed_data <- bed_data_reactive(input)
    
    # Auto-refresh when a new BED is uploaded, but only if a plot
    # already exists -- otherwise we'd trigger before chr/start/end
    # are set and stop() out.
    shiny::observeEvent(input$bed_file, {
      if (!is.null(rv$chr)) rv$trigger <- rv$trigger + 1
    })
    
    plot_reactive <- build_plot_reactive(input, ctx, rv, bed_data)
    register_outputs(input, output, session, plot_reactive)
  }
}

#' Parse the uploaded BED file
#'
#' Returns NULL when no file is uploaded. Parse errors surface as a
#' Shiny notification rather than crashing the app.
#'
#' @keywords internal
bed_data_reactive <- function(input) {
  shiny::reactive({
    if (is.null(input$bed_file)) return(NULL)
    tryCatch({
      bed_gr <- rtracklayer::import.bed(input$bed_file$datapath)
      data.table::as.data.table(as.data.frame(bed_gr))
    }, error = function(e) {
      shiny::showNotification(
        paste("Failed to parse BED:", conditionMessage(e)),
        type = "error", duration = 8
      )
      NULL
    })
  })
}

#' @keywords internal
register_navigation_observers <- function(input, session, ctx, rv) {
  shiny::observeEvent(input$target_gene, {
    shiny::req(input$target_gene)
    gene_rows <- ctx$anno_dt[type == "gene" & gene_name == input$target_gene]
    if (nrow(gene_rows) > 0) {
      t_chr   <- as.character(gene_rows$seqnames)[1]
      t_start <- min(gene_rows$start)
      t_end   <- max(gene_rows$end)
      shiny::updateTextInput(session,    "target_chr",   value = t_chr)
      shiny::updateNumericInput(session, "target_start", value = t_start)
      shiny::updateNumericInput(session, "target_end",   value = t_end)
      rv$home_chr   <- t_chr
      rv$home_start <- t_start
      rv$home_end   <- t_end
    }
  })
  
  shiny::observeEvent(input$plot_btn, {
    rv$chr     <- input$target_chr
    rv$start   <- as.integer(input$target_start)
    rv$end     <- as.integer(input$target_end)
    rv$trigger <- rv$trigger + 1
  })
  
  shiny::observeEvent(input$reset_btn, {
    shiny::req(rv$home_chr, rv$home_start, rv$home_end)
    rv$chr   <- rv$home_chr
    rv$start <- rv$home_start
    rv$end   <- rv$home_end
    shiny::updateTextInput(session,    "target_chr",   value = rv$chr)
    shiny::updateNumericInput(session, "target_start", value = rv$start)
    shiny::updateNumericInput(session, "target_end",   value = rv$end)
    rv$trigger <- rv$trigger + 1
  })
  
  shiny::observeEvent(input$region_brush, {
    brush <- input$region_brush
    rv$start <- as.integer(brush$xmin)
    rv$end   <- as.integer(brush$xmax)
    shiny::updateNumericInput(session, "target_start", value = rv$start)
    shiny::updateNumericInput(session, "target_end",   value = rv$end)
    session$resetBrush("region_brush")
    rv$trigger <- rv$trigger + 1
  })
}

#' @keywords internal
register_dynamic_filters <- function(input, output, ctx) {
  output$dynamic_group_filters_ui <- shiny::renderUI({
    shiny::req(length(input$groupings) > 0)
    inputs <- lapply(input$groupings, function(group_col) {
      values <- ctx$eiad_meta[[group_col]]
      char_values <- as.character(values)
      char_values[is.na(values)] <- "NA"
      choices <- sort(unique(char_values))
      shiny::selectizeInput(
        inputId  = paste0("dynamic_filter_", group_col),
        label    = paste("Filter by", group_col, ":"),
        choices  = choices,
        multiple = TRUE
      )
    })
    shiny::tagList(inputs)
  })
}

#' @keywords internal
filtered_meta_reactive <- function(input, ctx) {
  shiny::reactive({
    meta <- data.table::copy(ctx$eiad_meta)
    for (g in input$groupings) {
      sel <- input[[paste0("dynamic_filter_", g)]]
      if (is.null(sel) || length(sel) == 0) next
      if ("NA" %in% sel) {
        other <- sel[sel != "NA"]
        meta <- meta[is.na(get(g)) | get(g) %in% other]
      } else {
        meta <- meta[get(g) %in% sel]
      }
    }
    meta
  })
}

#' @keywords internal
build_plot_reactive <- function(input, ctx, rv, bed_data) {
  filtered_meta <- filtered_meta_reactive(input, ctx)
  
  shiny::eventReactive(rv$trigger, {
    shiny::req(rv$trigger > 0, rv$chr, rv$start, rv$end)
    t_chr   <- rv$chr
    t_start <- rv$start
    t_end   <- rv$end
    if (t_start >= t_end)
      stop("Start position must be less than End position.")
    
    # ---- Metadata pipeline
    meta_cur <- filtered_meta()
    meta_cur <- downsample_by_study(meta_cur, input$facet_group,
                                    input$max_samples)
    meta_cur <- sanitize_metadata(meta_cur)
    
    # ---- Annotation pipeline
    anno_bits <- subset_region_annotation(ctx$anno_dt,
                                          t_chr, t_start, t_end)
    
    # ---- Optional user-uploaded BED highlights
    bed_dt <- bed_data()
    bed_in_region <- if (is.null(bed_dt) || nrow(bed_dt) == 0) {
      data.table::data.table(start = numeric(0), end = numeric(0),
                             bed_tooltip = character(0))
    } else {
      sub <- bed_dt[seqnames == t_chr & start <= t_end & end >= t_start]
      build_bed_tooltips(sub)
    }
    
    # ---- Sample selection (CPM filter before BigWig I/O)
    target_samples <- cpm_filter_samples(
      meta_cur$sample_accession,
      input$target_gene, ctx, input$min_expr
    )
    if (length(target_samples) == 0)
      stop("No samples passed the log2(CPM+1) cutoff for ",
           input$target_gene, ".")
    meta_cur <- meta_cur[sample_accession %in% target_samples]
    
    # ---- Bin geometry
    bp_wide <- t_end - t_start
    b_size  <- if (input$bin_size > 0) input$bin_size
    else as.integer(max(1, bp_wide / 2000))
    n_bins  <- as.integer(max(1, ceiling(bp_wide / b_size)))
    
    # ---- BigWig read
    dt_full <- read_region_bigwigs(
      target_samples, ctx$bw_file_map,
      t_chr, t_start, t_end, n_bins,
      input$summary_type, ctx$bp_backend
    )
    
    # ---- Junction track? Reserve a thin per-sample band so junctions
    # don't bleed into the previous sample's wiggle.
    show_junctions <- !is.null(ctx$sj_se) && isTRUE(input$show_junctions)
    junc_band <- if (show_junctions) 0.4 else 0
    
    # ---- Plot data assembly (junc_band shifts the wiggle baseline up)
    pd_bits <- build_plot_data(dt_full, meta_cur, input$facet_group,
                               input$overlap_factor, input$summary_type,
                               junc_band = junc_band)
    
    junctions <- if (show_junctions) {
      raw <- read_region_junctions(
        ctx, t_chr, t_start, t_end,
        samples        = target_samples,
        min_reads      = as.integer(input$min_junc_reads)
      )
      attach_junction_positions(raw, pd_bits$unique_samples,
                                junc_band = junc_band)
    } else NULL
    
    color_var <- if (isTRUE(nzchar(input$color_by))) input$color_by else NULL
    
    # ---- Plots
    main_plot <- build_main_plot(
      pd_bits$plot_data,
      exon_highlights = anno_bits$exon_highlights,
      junctions       = junctions,
      bed_highlights  = bed_in_region,
      bed_color       = input$bed_color,
      chr             = t_chr,
      w_start         = t_start,
      w_end           = t_end,
      overlap_factor  = input$overlap_factor,
      color_var       = color_var,
      summary_type    = input$summary_type
    )
    
    mm <- build_minimap(
      anno_bits$region_anno, anno_bits$tx_base, anno_bits$tx_exons,
      anno_bits$has_transcripts, t_chr, t_start, t_end
    )
    
    dims <- compute_plot_dimensions(
      n_samples = nrow(pd_bits$unique_samples),
      n_facets  = length(unique(pd_bits$unique_samples$combined_facet)),
      n_tx      = mm$num_tx,
      has_color = !is.null(color_var),
      minimap_override = input$minimap_height,
      show_junctions   = show_junctions
    )
    
    list(
      plot            = main_plot,
      auto_height     = dims$main_px,
      plot_minimap    = mm$plot,
      minimap_px      = dims$minimap_px,
      tx_hover_info   = mm$tx_hover_info,
      has_transcripts = anno_bits$has_transcripts
    )
  })
}

#' @keywords internal
register_outputs <- function(input, output, session, plot_reactive) {
  output$minimap_container <- shiny::renderUI({
    shiny::req(plot_reactive())
    h <- plot_reactive()$minimap_px
    shiny::div(
      style = "position: relative;",
      shiny::plotOutput(
        "minimap",
        height = paste0(h, "px"),
        width  = "100%",
        brush  = shiny::brushOpts(id = "region_brush",
                                  direction = "x", resetOnNew = TRUE),
        hover  = shiny::hoverOpts(id = "minimap_hover",
                                  delay = 100, delayType = "debounce",
                                  nullOutside = TRUE)
      ),
      shiny::uiOutput("hover_info")
    )
  })
  
  output$minimap <- shiny::renderPlot(
    { shiny::req(plot_reactive()); plot_reactive()$plot_minimap },
    width  = function() {
      w <- session$clientData$output_bw_plot_width
      if (is.null(w) || w == 0) return(NULL)
      w
    },
    height = function() {
      h <- plot_reactive()$minimap_px
      if (is.null(h)) 150 else h
    },
    res = 72
  )
  
  output$hover_info <- shiny::renderUI({
    hover <- input$minimap_hover
    if (is.null(hover)) return(NULL)
    shiny::req(plot_reactive())
    info <- plot_reactive()$tx_hover_info
    if (is.null(info) || nrow(info) == 0) return(NULL)
    hits <- info[
      abs(tx_idx - hover$y) < 0.35 &
        pmin(start, end) <= hover$x &
        pmax(start, end) >= hover$x
    ]
    if (nrow(hits) == 0) return(NULL)
    hits <- hits[order(abs(tx_idx - hover$y))][1]
    
    style <- sprintf(
      paste0("position: absolute; left: %.0fpx; top: %.0fpx; ",
             "background: rgba(255,255,255,0.97); color: #222; ",
             "padding: 8px 10px; border-radius: 5px; ",
             "box-shadow: 2px 2px 5px rgba(0,0,0,0.2); ",
             "font-family: Arial, sans-serif; font-size: 0.82rem; ",
             "pointer-events: none; z-index: 100; max-width: 320px;"),
      hover$coords_css$x + 12, hover$coords_css$y + 12
    )
    shiny::div(
      style = style,
      shiny::HTML(paste0(
        "<b>Transcript:</b> ", hits$tx_label, "<br>",
        "<b>Coordinates:</b> ", hits$seqnames, ":",
        format(hits$start, big.mark = ","), "-",
        format(hits$end,   big.mark = ","), "<br>",
        "<b>Strand:</b> ", hits$strand
      ))
    )
  })
  
  output$bw_plot <- ggiraph::renderGirafe({
    shiny::req(plot_reactive())
    res <- plot_reactive()
    target_px <- if (!is.na(input$plot_height) && input$plot_height > 0) {
      input$plot_height
    } else {
      res$auto_height
    }
    ggiraph::girafe(
      ggobj      = res$plot,
      width_svg  = 16,
      height_svg = target_px / 72,
      options = list(
        ggiraph::opts_sizing(rescale = TRUE, width = 1),
        ggiraph::opts_toolbar(saveaspng = FALSE,
                              hidden = c("lasso_select", "lasso_deselect")),
        ggiraph::opts_selection(type = "none"),
        ggiraph::opts_tooltip(
          css = paste0("background-color: rgba(255,255,255,0.95); ",
                       "color: black; padding: 10px; border-radius: 5px; ",
                       "box-shadow: 2px 2px 5px rgba(0,0,0,0.2); ",
                       "font-family: Arial, sans-serif;"),
          use_fill = FALSE
        ),
        ggiraph::opts_hover(css = "stroke-width: 2px; stroke: #3A5836;")
      )
    )
  })
}