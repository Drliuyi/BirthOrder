# Table 1: current available baseline records. No artificial one-person padding.
describe_groups <- function(d, variables, group = "birth_group") {
  assert(group %in% names(d) && all(variables$column %in% names(d)), "Missing descriptive-table variables")
  groups <- c("Overall", levels(d[[group]])); result <- list()
  for (g in groups) {
    z <- if (g == "Overall") d else d[!is.na(d[[group]]) & d[[group]] == g, , drop = FALSE]
    for (j in seq_len(nrow(variables))) {
      v <- variables[j, ]; x <- z[[v$column]]; known <- !is.na(x)
      row <- data.frame(variable = v$column, label = v$label, group = g, group_N = nrow(z), valid_N = sum(known), missing_N = sum(!known),
        mean = NA_real_, SD = NA_real_, median = NA_real_, Q1 = NA_real_, Q3 = NA_real_, n = NA_integer_, percent = NA_real_)
      if (v$type == "continuous") {
        x <- numeric_value(x); x <- x[is.finite(x)]
        row$valid_N <- length(x); row$missing_N <- nrow(z) - length(x)
        if (length(x)) { row$mean <- mean(x); row$SD <- sd(x); row$median <- median(x); row$Q1 <- unname(quantile(x, .25)); row$Q3 <- unname(quantile(x, .75)) }
        row$level <- ""; result[[length(result) + 1L]] <- row
      } else {
        assert(v$type == "categorical", "Unknown Table 1 variable type")
        x <- as.character(x); x[is.na(x) | !nzchar(x)] <- "Missing"
        for (level in sort(unique(x))) {
          count <- row; count$level <- level; count$n <- sum(x == level)
          count$percent <- if (nrow(z)) 100 * count$n / nrow(z) else NA_real_
          result[[length(result) + 1L]] <- count
        }
      }
    }
  }
  do.call(rbind, result)
}
analyse_Table1 <- function(config, code_root) {
  d <- read_baseline(config)
  assert(file.exists(config$paths$withdrawals), "Current withdrawal list is required for Table 1")
  ids <- utils::read.csv(config$paths$withdrawals, header = FALSE, colClasses = "character")[[1]]
  d <- d[!as.character(d$eid) %in% ids, , drop = FALSE]
  complete <- complete.cases(d[c(config$primary_covariates, "S")])
  d <- d[complete & is.finite(d$age) & is.finite(d$S), , drop = FALSE]
  variables <- read_table(file.path(code_root, "reference", "table1.csv"))
  list(baseline = describe_groups(d, variables))
}
