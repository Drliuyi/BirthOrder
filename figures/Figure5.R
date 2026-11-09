# Figure 5b quantitative panels. Author artwork/path-diagram composition stays separate.
joint_panel <- function(domains, forest, model, config) {
  ids <- c("SMOKING", "BODY_COMPOSITION", "CLINICAL_LIPIDS", "GLYCAEMIA", "NMR_ALL", "CANDIDATE_PROTEINS")
  z <- domains[domains$model_id == model, , drop = FALSE]
  z <- z[order(match(z$domain_id, ids)), , drop = FALSE]
  rows <- do.call(rbind, lapply(c("IHD", "T2D", "Delta"), function(e) data.frame(
    outcome = e, row = rev(seq_len(nrow(z))), domain = z$domain_id,
    effect = as.numeric(z[[paste0(e,"_component_point_per_1000")]]),
    q = as.numeric(z[[paste0(e,"_q_BH")]]))))
  rows$outcome <- factor(rows$outcome, levels = c("IHD", "T2D", "Delta"))
  rows$label <- sprintf("%.2f", rows$effect); rows$mark <- saved_stars(rows$q)
  f <- forest[forest$model_id == model, , drop = FALSE]
  name <- c(S44_BASE11 = "11 Risk Factors", S44_NMR_PLUS40 = "+29 Metabolites", S44_PROT_PLUS15 = "+4 Proteins")[[model]]
  colour <- config$colours$population[[name]]
  side <- unique(rows[c("row", "domain")])
  heat <- ggplot2::ggplot(rows, ggplot2::aes(as.numeric(outcome), row, fill = effect)) +
    ggplot2::geom_tile(width = 1, height = 1, colour = "white", linewidth = .5) +
    ggplot2::geom_text(ggplot2::aes(label = label), family = "Arial", size = 2.7, colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = mark), nudge_x = .29, nudge_y = .22, family = "Arial", size = 3.2, colour = "black") +
    ggplot2::geom_rect(data = side, ggplot2::aes(xmin = .29, xmax = .47, ymin = row-.5, ymax = row+.5),
      inherit.aes = FALSE, fill = config$colours$domain[side$domain]) +
    ggplot2::annotate("rect", xmin = .5, xmax = 3.5, ymin = nrow(z)+.52, ymax = nrow(z)+.70, fill = colour) +
    ggplot2::scale_fill_gradient2(low = config$colours$negative, mid = "white", high = config$colours$positive,
      midpoint = 0, limits = c(-1.5, 1.5), name = "IE (per 1,000)", oob = scales::squish) +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = c("IHD", "T2D", "ΔIE"), position = "top", limits = c(.28, 3.5), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(.5,nrow(z)+.75), expand = c(0,0)) +
    ggplot2::labs(x = NULL, y = NULL) + paper_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), plot.margin = ggplot2::margin(0,0,0,0))
  pm <- as.numeric(f$PM)
  bar <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = 1, fill = "#DDE1E4") +
    ggplot2::annotate("rect", xmin = 0, xmax = pm, ymin = 0, ymax = 1, fill = colour) +
    ggplot2::annotate("text", x = 50, y = .5, label = sprintf("%.1f%%",pm), family = "Arial", fontface = "bold", size = 2.8) +
    ggplot2::scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) + ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0))
  point <- as.numeric(f$point_B0_per_1000); low <- as.numeric(f$CI_type7_low_per_1000); high <- as.numeric(f$CI_type7_high_per_1000)
  risk <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = 0, xend = 0, y = .1, yend = .8, linewidth = .25, linetype = 2, colour = "grey65") +
    ggplot2::annotate("segment", x = low, xend = high, y = .48, yend = .48, linewidth = .4, colour = colour) +
    ggplot2::annotate("point", x = point, y = .48, size = 1.9, colour = colour) +
    ggplot2::annotate("text", x = 2.5, y = .02,
      label = sprintf("%.3f (95%% CI %.3f–%.3f)",point,low,high), family = "Arial", size = 2.45) +
    ggplot2::scale_x_continuous(limits = c(0,5), breaks = c(0,5)) +
    ggplot2::scale_y_continuous(limits = c(-.15,1)) + ggplot2::labs(x = "Joint ΔIE", y = NULL) + paper_theme(7) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  patchwork::wrap_plots(heat, bar, risk, ncol = 1, heights = c(nrow(z),.43,1.6))
}

draw_Figure5 <- function(workbook, out, config) {
  require_packages(c("ggplot2", "patchwork", "svglite", "scales"))
  domains <- read_source(workbook, "heatmap"); forest <- read_source(workbook, "forest")
  # q comes from saved, model-by-effect families, never a new global adjustment.
  marks <- sum(vapply(c("IHD","T2D","Delta"),function(e)
    sum(suppressWarnings(as.numeric(domains[[paste0(e,"_q_BH")]]))<.05,na.rm=TRUE),integer(1)))
  stopifnot(marks == 14)
  models <- c("S44_BASE11", "S44_NMR_PLUS40", "S44_PROT_PLUS15")
  panels <- lapply(models, function(m) joint_panel(domains, forest, m, config))
  p <- patchwork::wrap_plots(panels, nrow = 1, guides = "collect") & ggplot2::theme(legend.position = "right")
  save_panel(p, out, "Figure5_JointMediation", height = 92)
}
