assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  assert(!length(missing), paste("Install required packages before running:", paste(missing, collapse = ", ")))
}
read_table <- function(path) {
  assert(file.exists(path), paste("Missing input:", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}
write_table <- function(data, out, name) {
  path <- file.path(out, name)
  assert(!file.exists(path), paste("Preserve existing output:", path))
  utils::write.csv(data, path, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  invisible(path)
}
read_frozen <- function(config, figure, name) read_table(file.path(config$paths$frozen, figure, name))
join_by_id <- function(left, right, id = "eid") {
  assert(id %in% names(left) && id %in% names(right), paste("Missing ID:", id))
  for (d in list(left, right)) assert(!anyNA(d[[id]]) && !anyDuplicated(as.character(d[[id]])), "Missing or duplicate IDs")
  columns <- setdiff(names(right), id)
  assert(!length(intersect(columns, names(left))), "Join would silently replace an existing column")
  cbind(left, right[match(as.character(left[[id]]), as.character(right[[id]])), columns, drop = FALSE])
}
new_output <- function(path, config, code_root) {
  norm <- function(x) tolower(normalizePath(x, winslash = "/", mustWork = FALSE))
  within <- function(x, y) identical(x, y) || startsWith(paste0(x, "/"), paste0(y, "/"))
  p <- norm(path)
  for (protected in c(config$paths$phenotype, config$paths$frozen, config$paths$helpers, code_root)) {
    q <- norm(protected)
    assert(!within(p, q) && !within(q, p), paste("Output overlaps protected input/code:", protected))
  }
  assert(!file.exists(path), paste("Choose a new output directory:", path))
  assert(dir.create(path, recursive = TRUE), "Could not create output directory")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
numeric_value <- function(x) suppressWarnings(as.numeric(as.character(x)))
bh_fixed <- function(p, family_size) {
  ok <- is.finite(p)
  assert(all(p[ok] >= 0 & p[ok] <= 1) && sum(ok) <= family_size, "Invalid correction family")
  q <- rep(NA_real_, length(p))
  if (any(ok)) q[ok] <- stats::p.adjust(p[ok], method = "BH", n = family_size)
  q
}
stars <- function(q) ifelse(!is.finite(q), "", ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", ""))))
