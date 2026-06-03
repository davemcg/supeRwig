#' Build the app server function
#' @keywords internal
build_server <- function(ctx) {
  function(input, output, session) {
    # Initialize the primary typeahead gene search field
    shiny::updateSelectizeInput(
      session, "target_gene",
      choices  = ctx$unique_genes,
      server   = TRUE,
      selected = if ("ABCA4" %in% ctx$unique_genes) "ABCA4" else ctx$unique_genes[1]
    )
    
    # Coordinates zoom and genomic window navigation state
    rv <- shiny::reactiveValues(
      chr = NULL, start = NULL, end = NULL, trigger = 0,
      home_chr = NULL, home_start = NULL, home_end = NULL
    )
    
    # Internal benchmarking metrics across backend operations
    timings_rv <- shiny::reactiveValues(
      bigwig = NA_real_, plot_data = NA_real_,
      main_plot = NA_real_, girafe = NA_real_
    )
    
    # Delegate behavior layout tasks to utility routines
    register_navigation_observers(input, session, ctx, rv)
    register_dynamic_filters(input, output, ctx)
    bed_data <- bed_data_reactive(input)
    
    # Fire up rendering graphs and structural export targets
    register_outputs(
      input, output, session,
      build_plot_reactive(input, ctx, rv, bed_data, timings_rv),
      timings_rv
    )
  }
}