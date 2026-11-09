# Plotting only. This directory never reads RDS or fits models.
read_source <- function(path, sheet) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required", call. = FALSE)
  as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text", trim_ws = FALSE, .name_repair = "minimal"))
}
plot_numbers <- function(d, columns) {
  for (name in columns) d[[name]] <- suppressWarnings(as.numeric(d[[name]]))
  d
}
paper_theme <- function(size = 8) {
  ggplot2::theme_classic(base_family = "Arial", base_size = size) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(linewidth = .25), axis.ticks = ggplot2::element_line(linewidth = .25),
      plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(), legend.key = ggplot2::element_blank())
}
save_panel <- function(plot, out, name, width = 180, height = 110) {
  for (ext in c("png", "pdf", "svg")) {
    path <- file.path(out, paste0(name, ".", ext))
    if (file.exists(path)) stop("Preserve existing figure: ", path, call. = FALSE)
    pdf_device <- if (Sys.info()[["sysname"]] == "Darwin")
      function(filename, width, height, bg = "white", ...) grDevices::quartz(file = filename, type = "pdf", width = width, height = height, bg = bg, family = "Arial") else grDevices::cairo_pdf
    device <- if (ext == "pdf") pdf_device else if (ext == "svg") svglite::svglite else "png"
    ggplot2::ggsave(path, plot, width = width, height = height, units = "mm", dpi = 600, device = device, bg = "white")
    if (!file.exists(path) || file.info(path)$size == 0) stop("Figure export failed: ", path, call. = FALSE)
  }
}
saved_stars <- function(q) ifelse(!is.finite(q), "", ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", ""))))
