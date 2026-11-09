draw_Figure3 <- function(workbook, out, config) {
  require_packages(c("ggplot2", "svglite"))
  d <- read_source(workbook, "models")
  d <- d[d$sample_set == "common_sample", , drop = FALSE]
  d <- plot_numbers(d, c("OR", "CI_low", "CI_high", "q_BH15"))
  p <- ggplot2::ggplot(d, ggplot2::aes(OR, disease, colour = model)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey65", linewidth = .25, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .15,
      position = ggplot2::position_dodge(.65), linewidth = .35) +
    ggplot2::geom_point(position = ggplot2::position_dodge(.65), size = 1.7) +
    ggplot2::labs(x = "Odds ratio (95% CI)", y = NULL, colour = NULL) + paper_theme()
  save_panel(p, out, "Figure3_Adjustment", height = 100)
  cox <- read_source(workbook, "cox")
  cox <- plot_numbers(cox, c("effect_value", "raw_HR", "raw_CI_low", "raw_CI_high"))
  cox <- cox[cox$trait %in% config$disease_codes, , drop = FALSE]
  # The reviewed table deliberately leaves pre-review HR/CI columns empty.
  # Use its explicitly approved effect_value and the matching saved raw CI.
  stopifnot(nrow(cox)==5, all(cox$numeric_status=="ESTIMABLE"),
    all(abs(cox$effect_value-cox$raw_HR)<1e-12), all(is.finite(cox$raw_CI_low)), all(is.finite(cox$raw_CI_high)))
  p <- ggplot2::ggplot(cox, ggplot2::aes(effect_value, trait)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey65", linewidth = .25, linetype = 2) +
    ggplot2::geom_segment(ggplot2::aes(x = raw_CI_low, xend = raw_CI_high, yend = trait), linewidth = .35) +
    ggplot2::geom_point(size = 1.8, colour = config$colours$exposure) +
    ggplot2::labs(x = "Hazard ratio (95% CI)", y = NULL) + paper_theme()
  save_panel(p, out, "Figure3_Incident", height = 90)
}
