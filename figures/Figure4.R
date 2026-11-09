draw_Figure4 <- function(workbook, out, config) {
  require_packages(c("ggplot2", "svglite"))
  d <- read_source(workbook, "contrasts")
  d <- plot_numbers(d, c("Delta_IE_per_1000", "q"))
  d$mark <- saved_stars(d$q)
  # Missing q is left missing. Never borrow IHD-T2D q for another contrast.
  p <- ggplot2::ggplot(d, ggplot2::aes(contrast, indicator_id, fill = Delta_IE_per_1000)) +
    ggplot2::geom_tile(colour = "white", linewidth = .4) +
    ggplot2::geom_text(ggplot2::aes(label = mark), size = 3.5, family = "Arial") +
    ggplot2::scale_fill_gradient2(low = config$colours$negative, mid = "white", high = config$colours$positive,
      midpoint = 0, name = expression(Delta*IE)) +
    ggplot2::labs(x = NULL, y = NULL) + paper_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_panel(p, out, "Figure4_SixContrasts", height = 150)
}
