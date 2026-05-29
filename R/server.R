#' Build the app server function
#' @keywords internal
build_server <- function(ctx) {
  function(input, output, session) {
    shiny::updateSelectizeInput(
      session, "target_gene",
      choices  = ctx$unique_genes,
      server   = TRUE,
      selected = if ("ABCA4" %in% ctx$unique_genes) "ABCA4" else ctx$unique_genes[1]
    )
    rv <- shiny::reactiveValues(
      chr = NULL, start = NULL, end = NULL, trigger = 0,
      home_chr = NULL, home_start = NULL, home_end = NULL,
      view_start = NULL, view_end = NULL
    )
    
    register_navigation_observers(input, session, ctx, rv)
    register_dynamic_filters(input, output, ctx)
    bed_data <- bed_data_reactive(input)
    
    register_outputs(
      input, output, session,
      build_plot_reactive(input, ctx, rv, bed_data)
    )
  }
}

#' @keywords internal
bed_data_reactive <- function(input) {
  shiny::reactive({
    if (is.null(input$bed_file)) return(NULL)
    tryCatch(
      data.table::as.data.table(
        as.data.frame(rtracklayer::import.bed(input$bed_file$datapath))
      ),
      error = function(e) {
        shiny::showNotification(
          paste("Failed to parse BED:", conditionMessage(e)),
          type = "error", duration = 8
        )
        NULL
      }
    )
  })
}

#' @keywords internal
register_navigation_observers <- function(input, session, ctx, rv) {
  shiny::observeEvent(input$target_gene, {
    shiny::req(input$target_gene)
    g_rows <- ctx$anno_dt[type == "gene" & gene_name == input$target_gene]
    if (nrow(g_rows) > 0) {
      rv$home_chr   <- as.character(g_rows$seqnames)[1]
      rv$home_start <- min(g_rows$start)
      rv$home_end   <- max(g_rows$end)
      shiny::updateTextInput(session,    "target_chr",   value = rv$home_chr)
      shiny::updateNumericInput(session, "target_start", value = rv$home_start)
      shiny::updateNumericInput(session, "target_end",   value = rv$home_end)
    }
  })
  
  # Read region: the slab of BigWig that's actually been pulled. Only
  # changes (and re-reads) on Generate Plot or Reset Zoom. View region
  # resets to match the read region on each pull.
  set_read_region <- function(c, s, e) {
    rv$chr        <- c
    rv$start      <- as.integer(s)
    rv$end        <- as.integer(e)
    rv$view_start <- rv$start
    rv$view_end   <- rv$end
    rv$trigger    <- rv$trigger + 1
  }
  
  shiny::observeEvent(input$plot_btn, {
    set_read_region(input$target_chr, input$target_start, input$target_end)
  })
  
  shiny::observeEvent(input$reset_btn, {
    shiny::req(rv$home_chr)
    shiny::updateTextInput(session,    "target_chr",   value = rv$home_chr)
    shiny::updateNumericInput(session, "target_start", value = rv$home_start)
    shiny::updateNumericInput(session, "target_end",   value = rv$home_end)
    set_read_region(rv$home_chr, rv$home_start, rv$home_end)
  })
  
  # Brushing the minimap only moves the view window — no re-read.
  # Inputs reflect the new view so Generate Plot would re-pull at
  # higher resolution for the zoomed-in coords if the user wants that.
  shiny::observeEvent(input$region_brush, {
    rv$view_start <- as.integer(input$region_brush$xmin)
    rv$view_end   <- as.integer(input$region_brush$xmax)
    shiny::updateNumericInput(session, "target_start", value = rv$view_start)
    shiny::updateNumericInput(session, "target_end",   value = rv$view_end)
    session$resetBrush("region_brush")
  })
}

#' @keywords internal
register_dynamic_filters <- function(input, output, ctx) {
  output$dynamic_group_filters_ui <- shiny::renderUI({
    shiny::req(length(input$groupings) > 0)
    shiny::tagList(lapply(input$groupings, function(g) {
      choices <- sort(unique(
        ifelse(is.na(ctx$eiad_meta[[g]]), "NA", as.character(ctx$eiad_meta[[g]]))
      ))
      shiny::selectizeInput(
        inputId  = paste0("dynamic_filter_", g),
        label    = paste("Filter by", g, ":"),
        choices  = choices,
        multiple = TRUE
      )
    }))
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
    t0 <- Sys.time()
    g <- ggiraph::girafe(
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
    cat("girafe:", round(as.numeric(Sys.time() - t0), 2), "s\n")
    g
  })
}