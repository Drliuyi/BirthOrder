# Figure 3: sequential adjustment, full-follow-up Cox and probability gradients.
analyse_Figure3 <- function(config, code_root) {
  assert(file.exists(config$paths$incident), "Authoritative full-follow-up incident risk sets are required before starting any fits")
  risks <- readRDS(config$paths$incident)
  assert(is.data.frame(risks) && all(c("eid", "outcome", "entry", "exit", "event", "eligible") %in% names(risks)), "Invalid incident input schema; no models started")
  d <- read_baseline(config)
  assert(file.exists(config$paths$withdrawals), "Current withdrawal list is required for the updated adjustment roster")
  withdrawn <- utils::read.csv(config$paths$withdrawals, header = FALSE, colClasses = "character")[[1]]
  withdrawn <- setdiff(withdrawn, c("eid", "participant_id"))
  assert(all(grepl("^[0-9]+$", withdrawn)), "Invalid withdrawal ID schema")
  d <- d[!as.character(d$eid) %in% withdrawn, , drop = FALSE]
  d <- join_by_id(d, read_diseases(config))
  d <- adjustment_covariates(d, config)
  result <- adjustment_models(d, config$disease_codes)
  result$cox <- incident_cox(d, risks, config)
  result
}
