# Figure 1 uses a separately supplied illustration.
draw_Figure1 <- function(workbook, out, config) {
  asset <- Sys.getenv("BIRTH_FIGURE1_ARTWORK", "")
  if (!nzchar(asset) || !file.exists(asset)) stop("Set BIRTH_FIGURE1_ARTWORK to the source illustration.", call. = FALSE)
  target <- file.path(out, paste0("Figure1.", tools::file_ext(asset)))
  stopifnot(!file.exists(target), file.copy(asset, target, overwrite = FALSE))
}
