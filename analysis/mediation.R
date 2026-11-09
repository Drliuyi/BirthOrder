# Effects are 11-year absolute risk differences; never infer b as IE/a.
mediation_effects <- function(R00, R10, R11) {
  r <- c(R00, R10, R11)
  assert(length(R00) == length(R10) && length(R10) == length(R11) && all(is.finite(r) & r >= 0 & r <= 1), "Invalid 11-year risks")
  TE <- R11 - R00; DE <- R10 - R00; IE <- R11 - R10
  data.frame(TE = 1000 * TE, DE = 1000 * DE, IE = 1000 * IE, PM = ifelse(TE == 0, NA_real_, 100 * IE / TE))
}
select_mediators <- function(master, check_counts = TRUE) {
  columns <- c("indicator_id", "IHD_q_BH_global605_exploratory", "T2D_q_BH_global605_exploratory",
    "paired_Delta_IE_per_1000", "paired_q_BH_within_contrast121_exploratory", "IHD_PM_percent", "T2D_PM_percent")
  assert(nrow(master) == 121 && !anyDuplicated(master$indicator_id) && all(columns %in% names(master)), "Complete current 121-indicator master required")
  x <- lapply(master[setdiff(columns, "indicator_id")], numeric_value)
  master$pass_IE <- x$IHD_q_BH_global605_exploratory < .05 | x$T2D_q_BH_global605_exploratory < .05
  master$pass_direction <- x$paired_Delta_IE_per_1000 > 0
  master$pass_pair <- x$paired_q_BH_within_contrast121_exploratory < .05
  master$pass_three <- master$pass_IE & master$pass_direction & master$pass_pair
  master$pass_PM <- x$IHD_PM_percent > 5 | x$T2D_PM_percent > 5
  master$selected <- master$pass_three & master$pass_PM
  if (check_counts) assert(sum(master$pass_three, na.rm = TRUE) == 64 && sum(master$selected, na.rm = TRUE) == 44,
    "Frozen 121 -> 64 -> 44 selection differs; review input version/fields, never change rules")
  master
}
paired_effects <- function(first, second, value = "IE", id = "bootstrap_id") {
  assert(!anyDuplicated(first[[id]]) && !anyDuplicated(second[[id]]), "Duplicate Bootstrap IDs")
  keep <- intersect(first[[id]][is.finite(first[[value]])], second[[id]][is.finite(second[[value]])])
  assert(all(keep %in% 1:1000), "Invalid saved Bootstrap IDs")
  a <- first[[value]][match(keep, first[[id]])]; b <- second[[value]][match(keep, second[[id]])]
  data.frame(bootstrap_id = keep, first = a, second = b, difference = a - b)
}
read_paper_sheet <- function(config, sheet) {
  require_packages("readxl")
  assert(file.exists(config$workbook), paste("Missing frozen supplementary workbook:", config$workbook))
  assert(sheet %in% readxl::excel_sheets(config$workbook), paste("Missing full source sheet:", sheet))
  as.data.frame(readxl::read_excel(config$workbook, sheet = sheet, skip = 3, col_types = "text", trim_ws = FALSE, .name_repair = "minimal"))
}
