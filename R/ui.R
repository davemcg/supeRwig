#' Build the app UI
#'
#' @keywords internal
build_ui <- function(ctx) {
  bslib::page_sidebar(
    title   = "supeRwig",
    theme   = bslib::bs_theme(version = 5, primary = "#3A5836"),
    sidebar = build_sidebar(ctx),
    build_plot_viewer_panel()
  )
}

#' @keywords internal
build_sidebar <- function(ctx) {
  meta_cols <- colnames(ctx$eiad_meta)
  
  # Junction panel only appears when a junction SE was supplied.
  # Junction panel updated for dual simultaneous filtering
  junction_panel <- if (!is.null(ctx$sj_se)) {
    bslib::accordion_panel(
      "Junction Track",
      shiny::checkboxInput("show_junctions",
                           "Show per-sample junction track",
                           value = FALSE),
      shiny::numericInput("min_psi5",
                          "Min PSI5 % Cutoff (5' Donor):",
                          value = 1, min = 0.01, max = 100, step = 0.1),
      shiny::numericInput("min_psi3",
                          "Min PSI3 % Cutoff (3' Acceptor):",
                          value = 1, min = 0.01, max = 100, step = 0.1),
      shiny::helpText(
        "Junctions appear as straight horizontal lines. Traces must pass ",
        "both cutoffs to be displayed."
      )
    )
  } else NULL
  
  bslib::sidebar(
    width = 350,
    shiny::h4("Target & Parameters"),
    shiny::selectizeInput("target_gene", "Search by Gene:",
                          choices = NULL, multiple = FALSE),
    
    shiny::textInput("target_region", "Search by Region:",
                     value = "chr1:93,992,834-94,121,148",
                     placeholder = "chr1:93,992,834-94,121,148"),
    shiny::hr(class = "my-2"),
    shiny::selectizeInput(
      "facet_group", "Facet / Group By:",
      choices  = meta_cols,
      selected = intersect(c("Tissue", "Sub_Tissue", "Source",
                             "Age", "Perturbation"),
                           meta_cols),
      multiple = TRUE
    ),
    shiny::selectizeInput(
      "color_by", "Color Wiggle Lines By (optional):",
      choices  = c("None" = "", meta_cols),
      selected = if ("study_accession" %in% meta_cols) "study_accession"
      else "",
      multiple = FALSE
    ),
    do.call(bslib::accordion, c(
      list(open = FALSE),
      list(
        bslib::accordion_panel(
          "Data Filters",
          shiny::selectizeInput("groupings", "Metadata to Filter By:",
                                choices = meta_cols, multiple = TRUE),
          shiny::uiOutput("dynamic_group_filters_ui"),
          shiny::hr(class = "my-2"),
          shiny::numericInput("max_samples",
                              "Max Samples per Study (0 = All):",
                              value = 3, min = 0)
        ),
        bslib::accordion_panel(
          "Plot Settings",
          shiny::numericInput("min_expr", "Min log2(CPM+1) for gene:",
                              value = 5, min = 0, step = 0.5),
          shiny::numericInput("bin_size",
                              "Bin Size (base pairs) [0 = Auto]:",
                              value = 0, min = 0),
          shiny::numericInput("overlap_factor", "Overlap Factor:",
                              value = 1.2, step = 0.1),
          shiny::numericInput("plot_height",
                              "Plot Height (pixels) [0 = Auto]:",
                              value = 0, min = 0),
          shiny::numericInput("minimap_height",
                              "Minimap Height (pixels) [0 = Auto]:",
                              value = 0, min = 0)
        ),
        bslib::accordion_panel(
          "BED Highlights",
          shiny::fileInput(
            "bed_file", "Upload BED file:",
            accept = c(".bed", ".bed.gz", ".txt", ".tsv")
          ),
          colourpicker::colourInput(
            "bed_color", "Highlight color:",
            value = "#B22222", showColour = "background"
          ),
          shiny::helpText(
            "Regions in the uploaded BED that overlap the current ",
            "window are drawn as vertical bands behind the wiggle ",
            "traces. New uploads refresh automatically; color and ",
            "opacity apply on next plot generation."
          )
        )
      ),
      # junction_panel goes next if present
      if (!is.null(junction_panel)) list(junction_panel) else list(),
      # --- EXPORT ---
      list(
        bslib::accordion_panel(
          "Export Options",
          shiny::p("Download data and code to reproduce and customize this plot locally:"),
          shiny::downloadButton("download_data", "Download Plot Data (.rds)", class = "btn-outline-secondary w-100 mb-2"),
          shiny::downloadButton("download_script", "Download R Script (.R)", class = "btn-outline-secondary w-100")
        )
      )
    )),
    shiny::actionButton("plot_btn", "Generate Plot",
                        class = "btn-primary w-100")
  )
}

#' @keywords internal
build_plot_viewer_panel <- function() {
  shiny::tagList(
    bslib::card(
      fill  = FALSE,
      class = "p-0 mb-2",
      bslib::card_header(
        class = paste("bg-light py-1 px-2",
                      "d-flex justify-content-between align-items-center"),
        shiny::span(paste0("Minimap (drag = zoom X-Axis  ·  ",
                           "hover transcripts for details)")),
        shiny::actionButton(
          "reset_btn", "Reset Zoom",
          class = "btn-sm btn-outline-secondary",
          style = "padding: 0.1rem 0.5rem; font-size: 0.8rem;"
        )
      ),
      bslib::card_body(fill = FALSE, padding = 0,
                       shiny::uiOutput("minimap_container"))
    ),
    bslib::card(
      full_screen = TRUE,
      class = "p-0",
      bslib::card_header(
        class = "bg-light py-1 px-2",
        shiny::div(
          class = "d-flex align-items-center w-100",
          style = "gap: 0.5rem;",
          shiny::div(
            style = "flex: 1; min-width: 0; text-align: left; display: flex; align-items: center; gap: 0.35rem;",
            shiny::span(style = "white-space: nowrap; overflow: hidden; text-overflow: ellipsis;", "Coverage (cpm)"),
            bslib::tooltip(
              shiny::span(style = "cursor: pointer; font-size: 0.9rem; color: #6c757d;", "ℹ"),
              shiny::uiOutput("timing_info"),
              placement = "right"
            )
          ),
          shiny::div(
            style = paste("flex: 2; min-width: 0; text-align: center;",
                          "white-space: nowrap; overflow: hidden;",
                          "text-overflow: ellipsis;"),
            shiny::textOutput("bw_facet_label", inline = TRUE)
          ),
          shiny::div(
            style = paste("flex: 1; min-width: 0; text-align: right;",
                          "white-space: nowrap; overflow: hidden;",
                          "text-overflow: ellipsis;"),
            shiny::textOutput("bw_gene_label", inline = TRUE)
          )
        )
      ),
      bslib::card_body(
        fill = TRUE, padding = 0,
        ggiraph::girafeOutput("bw_plot", width = "100%", height = "100%")
      )
    )
  )
}