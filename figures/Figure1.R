# Figure 1 contains author-composed illustrations, not an automatically fitted graph.
draw_Figure1 <- function(workbook, out, config) {
  asset <- Sys.getenv("BIRTH_FIGURE1_ARTWORK", "")
  if (!nzchar(asset) || !file.exists(asset)) stop("Supply BIRTH_FIGURE1_ARTWORK: the approved original editable artwork. Source numbers are in Figure1.xlsx; no concept image is invented.", call. = FALSE)
  target <- file.path(out, paste0("Figure1.", tools::file_ext(asset)))
  stopifnot(!file.exists(target), file.copy(asset, target, overwrite = FALSE))
}
