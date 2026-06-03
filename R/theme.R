#' Shared plot theme
#'
#' @keywords internal
theme_panel_only <- function() {
  cowplot::theme_minimal_vgrid() +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(t = 5, r = 8, b = 5, l = 8)
    )
}

#' Categorical palette for arbitrary numbers of levels
#'
#' Concatenates several `pals` palettes to handle large categorical
#' variables (study_accession etc.) without recycling.
#'
#' @keywords internal
cat_palette <- function(levels) {
  n <- length(levels)
  if (n == 0) return(character(0))
  pal <- c(pals::cols25()[-c(6,7,13,14)], 
           pals::polychrome()[-c(1,2,20)],
           pals::glasbey(),
           pals::okabe())[seq_len(n)]
  setNames(unname(pal), levels)
}