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

#' @keywords internal
register_navigation_observers <- function(input, session, ctx, rv) {
  shiny::observeEvent(input$target_gene, {
    shiny::req(input$target_gene)
    g_rows <- ctx$anno_dt[type == "gene" & gene_name == input$target_gene]
    if (nrow(g_rows) > 0) {
      rv$home_chr   <- as.character(g_rows$seqnames)[1]
      rv$home_start <- min(g_rows$start)
      rv$home_end   <- max(g_rows$end)
      shiny::updateTextInput(
        session, "target_region",
        value = format_ucsc_region(rv$home_chr, rv$home_start, rv$home_end)
      )
    }
  })
  
  set_read_region <- function(c, s, e) {
    rv$chr     <- c
    rv$start   <- as.integer(s)
    rv$end     <- as.integer(e)
    rv$trigger <- rv$trigger + 1
  }
  
  shiny::observeEvent(input$plot_btn, {
    parsed <- parse_ucsc_region(input$target_region)
    if (is.null(parsed)) {
      shiny::showNotification(
        "Region must look like chr1:93,992,834-94,121,148 (start < end).",
        type = "error", duration = 6
      )
      return()
    }
    set_read_region(parsed$chr, parsed$start, parsed$end)
  })
  
  shiny::observeEvent(input$reset_btn, {
    shiny::req(rv$home_chr)
    shiny::updateTextInput(
      session, "target_region",
      value = format_ucsc_region(rv$home_chr, rv$home_start, rv$home_end)
    )
    set_read_region(rv$home_chr, rv$home_start, rv$home_end)
  })
  
  shiny::observeEvent(input$region_brush, {
    s <- as.integer(input$region_brush$xmin)
    e <- as.integer(input$region_brush$xmax)
    shiny::updateTextInput(
      session, "target_region",
      value = format_ucsc_region(rv$chr, s, e)
    )
    set_read_region(rv$chr, s, e)
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
register_outputs <- function(input, output, session, plot_reactive, timings_rv) {
  output$minimap_container <- shiny::renderUI({
    shiny::req(plot_reactive())
    h <- plot_reactive()$minimap_px
    shiny::div(
      style = "position: relative;",
      shiny::plotOutput(
        "minimap",
        height = paste0(h, "px"),
        width  = "100%",
        brush  = shiny::brushOpts(id = "region_brush", direction = "x", resetOnNew = TRUE),
        hover  = shiny::hoverOpts(id = "minimap_hover", delay = 100, delayType = "debounce", nullOutside = TRUE)
      ),
      shiny::uiOutput("hover_info")
    )
  })
  
  output$bw_facet_label <- shiny::renderText({
    shiny::req(plot_reactive())
    fc <- plot_reactive()$facet_cols
    if (length(fc) == 0) "" else paste(fc, collapse = " - ")
  })
  
  output$bw_gene_label <- shiny::renderText({
    shiny::req(plot_reactive())
    g <- plot_reactive()$target_gene
    if (is.null(g) || !nzchar(g)) "" else g
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
    
    # --- ADDED: DYNAMIC BOUNDARY DETECTION ---
    # Detect total height allocated to the minimap
    total_h <- plot_reactive()$minimap_px
    if (is.null(total_h) || total_h == 0) total_h <- 150 # Safe fallback height
    
    # Calculate remaining headroom below the mouse pointer
    space_below <- total_h - hover$coords_css$y
    
    # Default state: draw tooltip offset slightly down
    y_pos <- hover$coords_css$y + 12
    transform_css <- ""
    
    # If a multi-line tooltip has less than 120px layout clearance, flip it UP
    if (space_below < 120) {
      y_pos <- hover$coords_css$y - 12
      transform_css <- "transform: translateY(-100%); "
    }
    
    # Rebuilt style string via safe concatenation to prevent sprintf string injection bugs
    style <- paste0(
      "position: absolute; ",
      "left: ", hover$coords_css$x + 12, "px; ",
      "top: ", y_pos, "px; ",
      transform_css,
      "background: rgba(255,255,255,0.97); color: #222; ",
      "padding: 8px 10px; border-radius: 5px; ",
      "box-shadow: 2px 2px 5px rgba(0,0,0,0.2); ",
      "font-family: Arial, sans-serif; font-size: 0.82rem; ",
      "pointer-events: none; z-index: 100; max-width: 320px;"
    )
    # ----------------------------------------
    
    # Extract rich metadata targets
    type_col <- intersect(c("transcript_type", "transcript_biotype", "biotype"), colnames(hits))
    tx_type  <- if (length(type_col) > 0) as.character(hits[[type_col[1]]]) else "N/A"
    if (is.na(tx_type) || tx_type == "") tx_type <- "N/A"
    
    tx_tag <- if ("tag" %in% colnames(hits)) as.character(hits$tag) else "N/A"
    if (is.na(tx_tag) || tx_tag == "") tx_tag <- "N/A"
    
    shiny::div(
      style = style,
      shiny::HTML(paste0(
        "<b>Transcript:</b> ", hits$tx_label, "<br>",
        "<b>Coordinates:</b> ", hits$seqnames, ":",
        format(hits$start, big.mark = ","), "-",
        format(hits$end,   big.mark = ","), "<br>",
        "<b>Strand:</b> ", hits$strand, "<br>",
        "<b>Type:</b> ", tx_type, "<br>",
        "<b>Tag:</b> ", tx_tag
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
    timings_rv$girafe <- as.numeric(Sys.time() - t0)
    g
  })
  
  output$timing_info <- shiny::renderUI({
    fmt <- function(x) if (is.na(x)) "—" else sprintf("%.2f s", x)
    shiny::HTML(paste0(
      "<div style='font-family: monospace; text-align: left; white-space: pre; font-size: 0.8rem; line-height: 1.4;'>",
      "BigWig read     : ", fmt(timings_rv$bigwig), "<br>",
      "build_plot_data : ", fmt(timings_rv$plot_data), "<br>",
      "build_main_plot : ", fmt(timings_rv$main_plot), "<br>",
      "girafe (server) : ", fmt(timings_rv$girafe), "<br>",
      "------------------------<br>",
      "Total (server)  : ", fmt(sum(c(timings_rv$bigwig, timings_rv$plot_data,
                                      timings_rv$main_plot, timings_rv$girafe), na.rm = TRUE)),
      "</div>"
    ))
  })
  
  # ---- Data Download Handler ----
  output$download_data <- shiny::downloadHandler(
    filename = function() {
      shiny::req(plot_reactive())
      g <- plot_reactive()$gated_params
      paste0("supeRwig_data_", g$chr, "_", g$start, "_", g$end, ".rds")
    },
    content = function(file) {
      shiny::req(plot_reactive())
      res <- plot_reactive()
      export_bundle <- list(
        plot_data       = res$plot_data,
        exon_highlights = res$exon_highlights,
        junctions       = res$junctions,
        bed_highlights  = res$bed_highlights,
        color_var       = res$color_var,
        chr             = res$gated_params$chr,
        w_start         = res$gated_params$start,
        w_end           = res$gated_params$end,
        overlap_factor  = input$overlap_factor,
        bed_color       = input$bed_color
      )
      saveRDS(export_bundle, file)
    }
  )
  
  # ---- Script Download Handler ----
  output$download_script <- shiny::downloadHandler(
    filename = function() {
      shiny::req(plot_reactive())
      g <- plot_reactive()$gated_params
      paste0("supeRwig_plot_", g$chr, "_", g$start, "_", g$end, ".R")
    },
    content = function(file) {
      shiny::req(plot_reactive())
      res <- plot_reactive()
      g   <- res$gated_params
      data_file_name <- paste0("supeRwig_data_", g$chr, "_", g$start, "_", g$end, ".rds")
      
      script_lines <- c(
        "# ==========================================================================",
        "# supeRwig: Local Plot Customization Script",
        paste0("# Region: ", g$chr, ":", g$start, "-", g$end),
        "# ==========================================================================",
        "",
        "library(ggplot2)",
        "library(data.table)",
        "if (!requireNamespace('ggforce', quietly = TRUE)) install.packages('ggforce')",
        "if (!requireNamespace('ggnewscale', quietly = TRUE)) install.packages('ggnewscale')",
        "",
        "# 1. Load Data Bundle",
        paste0("data_path <- '", data_file_name, "'"),
        "if (!file.exists(data_path)) {",
        "  stop(paste('Could not find data file:', data_path, '\\nEnsure it is in your current working directory.'))",
        "}",
        "bundle <- readRDS(data_path)",
        "",
        "# 2. Re-map Plot Attributes",
        "plot_data       <- bundle$plot_data",
        "exon_highlights <- bundle$exon_highlights",
        "junctions       <- bundle$junctions",
        "bed_highlights  <- bundle$bed_highlights",
        "color_var       <- bundle$color_var",
        "chr             <- bundle$chr",
        "w_start         <- bundle$w_start",
        "w_end           <- bundle$w_end",
        "overlap_factor  <- bundle$overlap_factor",
        "bed_color       <- bundle$bed_color",
        "",
        "# 3. Handle Color Profiles",
        "if (!is.null(color_var) && color_var %in% colnames(plot_data)) {",
        "  plot_data[, line_color := data.table::fcoalesce(as.character(get(color_var)), 'NA')]",
        "  unique_groups  <- sort(unique(plot_data$line_color))",
        "  line_colors    <- scales::hue_pal()(length(unique_groups))",
        "  names(line_colors) <- unique_groups",
        "} else {",
        "  plot_data[, line_color := '_default']",
        "  line_colors    <- c('_default' = 'grey15')",
        "}",
        "",
        "# Junction Palette Definition",
        "junc_palette <- c(",
        "  '+/annot' = 'royalblue4', '+/novel' = 'royalblue1',",
        "  '-/annot' = 'tomato4',    '-/novel' = 'tomato1',",
        "  '*/annot' = 'seagreen4',  '*/novel' = 'seagreen1'",
        ")",
        "",
        "# 4. Build ggplot Construction",
        "p <- ggplot(plot_data) +",
        "  geom_rect(data = exon_highlights,",
        "            aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf),",
        "            inherit.aes = FALSE, fill = 'grey85', alpha = 0.5) +",
        "  geom_rect(data = bed_highlights,",
        "            aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf),",
        "            inherit.aes = FALSE, fill = bed_color, alpha = 0.3) +",
        "  geom_step(aes(x = plot_x, y = offset_y, group = sample, color = line_color),",
        "            direction = 'mid', linewidth = 0.4) +",
        "  scale_y_continuous(breaks = NULL, expand = expansion(add = c(0.1, overlap_factor - 1))) +",
        "  scale_x_continuous(labels = function(x) format(x, big.mark = ',', scientific = FALSE)) +",
        "  scale_color_manual(values = line_colors, name = color_var) +",
        "  coord_cartesian(xlim = c(w_start, w_end)) +",
        "  labs(",
        "    title = sprintf('Region: %s:%s-%s', chr, format(w_start, big.mark = ','), format(w_end, big.mark = ',')),",
        "    x = 'Genomic Position', y = 'Coverage Profiles'",
        "  ) +",
        "  theme_minimal() +",
        "  theme(",
        "    panel.spacing.y  = unit(0.2, 'lines'),",
        "    strip.text.y     = element_text(angle = 0, hjust = 0),",
        "    panel.grid.major = element_blank(),",
        "    panel.grid.minor = element_blank(),",
        "    legend.position  = 'right'",
        "  ) +",
        "  ggforce::facet_col(vars(combined_facet), scales = 'free_y', space = 'free', shrink = TRUE)",
        "",
        "# 5. Superimpose Junction Tracks (if present)",
        "if (!is.null(junctions) && nrow(junctions) > 0) {",
        "  p <- p +",
        "    ggnewscale::new_scale_color() +",
        "    geom_segment(data = junctions,",
        "                 aes(x = start, xend = end, y = junc_y, yend = junc_y, color = strand_annot, linewidth = junc_lw),",
        "                 inherit.aes = FALSE, lineend = 'round') +",
        "    scale_color_manual(values = junc_palette, name = 'Junction Status') +",
        "    scale_linewidth_identity() +",
        "    guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2))) +",
        "    theme(legend.position = 'bottom')",
        "}",
        "",
        "print(p)",
        "NULL"
      )
      writeLines(script_lines, file)
    }
  )
}

#' Filter samples by log2(CPM+1) of a target gene
#'
#' Runs before the BigWig read so we skip I/O for samples that won't
#' meet the cutoff anyway. Returns the input unchanged if no gene is
#' selected or the gene isn't in the SE.
#'
#' @keywords internal
cpm_filter_samples <- function(samples, gene_name, ctx, min_log_cpm) {
  if (is.null(gene_name) || !nzchar(gene_name)) return(samples)
  
  se_rows <- ctx$gene_to_se_row[[gene_name]]
  if (is.null(se_rows)) {
    warning("Gene '", gene_name, "' not in SE; skipping CPM filter.")
    return(samples)
  }
  
  cpm_vec <- if (length(se_rows) == 1) {
    as.numeric(ctx$cpm_assay[se_rows, ])
  } else {
    colSums(as.matrix(ctx$cpm_assay[se_rows, , drop = FALSE]))
  }
  names(cpm_vec) <- names(ctx$se_col_for_sample)
  
  log_cpm <- log2(cpm_vec[samples] + 1)
  keep    <- !is.na(log_cpm) & log_cpm >= min_log_cpm
  
  cat(sprintf("CPM filter (log2(CPM+1) >= %g for %s): %d/%d kept\n",
              min_log_cpm, gene_name, sum(keep), length(keep)))
  
  samples[keep]
}

#' Cap samples per (study x facet-group) combination
#'
#' @keywords internal
downsample_by_study <- function(meta, facet_cols, n) {
  if (n <= 0) return(meta)
  ds_by <- intersect(c("study_accession", facet_cols), colnames(meta))
  if (length(ds_by) == 0) return(meta)
  meta <- meta[order(sample_accession)]
  meta[, utils::head(.SD, n), by = ds_by]
}

#' Escape characters that would break ggiraph tooltips
#'
#' @keywords internal
sanitize_metadata <- function(meta) {
  char_cols <- names(meta)[vapply(meta, is.character, logical(1))]
  for (col in char_cols) {
    safe <- meta[[col]]
    Encoding(safe) <- "UTF-8"
    safe <- iconv(safe, "UTF-8", "UTF-8", sub = "")
    safe[is.na(safe) & !is.na(meta[[col]])] <- ""
    safe <- gsub("\"", "&quot;", gsub("'", "&#39;", safe, useBytes = TRUE), useBytes = TRUE)
    Encoding(safe) <- "UTF-8"
    data.table::set(meta, j = col, value = safe)
  }
  meta
}