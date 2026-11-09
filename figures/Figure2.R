draw_Figure2 <- function(workbook, out, config) {
  require_packages(c("ggplot2", "svglite"))
  d <- read_source(workbook, "phewas")
  d <- plot_numbers(d, c("beta", "p", "OR", "CI_low", "CI_high", "descriptive_BH_q"))
  d$position <- seq_len(nrow(d)); d$signed_logP <- sign(d$beta) * -log10(d$p)
  categories <- unique(d$category); d$category <- factor(d$category, levels = categories)
  p <- ggplot2::ggplot(d, ggplot2::aes(position, signed_logP, colour = category)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = .2) +
    ggplot2::geom_point(size = .75, alpha = .65, na.rm = TRUE) +
    ggplot2::labs(x = NULL, y = expression(sign(beta) %*% -log[10](P))) + paper_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "right")
  save_panel(p, out, "Figure2_PheWAS")
  selected <- d[is.finite(d$descriptive_BH_q) & d$descriptive_BH_q < .05, , drop = FALSE]
  forest <- ggplot2::ggplot(selected, ggplot2::aes(OR, reorder(description, OR))) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey65", linetype = 2, linewidth = .25) +
    ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = reorder(description, OR)), linewidth = .35) +
    ggplot2::geom_point(ggplot2::aes(colour = beta > 0), size = 1.8) +
    ggplot2::scale_colour_manual(values = c(config$colours$negative, config$colours$positive), guide = "none") +
    ggplot2::labs(x = "Odds ratio (95% CI)", y = NULL) + paper_theme()
  save_panel(forest, out, "Figure2_Associations", height = 100)
}
