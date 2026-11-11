# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ev_diabetes_prior_specs <- function() data.frame(column = c("fod_gp2_dm", "fod_gp2_t2dm", "fod_gp3_dm", "fod_gp3_t2dm", 
    "fod_gp_dm", "fod_gp_t2dm", "fod_srd_dm", "fod_srd_t2dm", "fod_ref_t2dm"), source = c(rep("local_GP_mapping_positive_only", 
    6L), rep("SRD_year_resolution_positive_only", 2L), "UKB130708_first_occurrence_E11"), precision = c(rep("day_or_placeholder", 
    6L), rep("year", 2L), "reported_date_mixed_source"), stringsAsFactors = FALSE)

ev_diabetes_prior_contract_example <- function() list(verified = FALSE, evidence = "", source_columns = character(), 
    source_path = "", baseline_comparison_verified = FALSE)

ev_diabetes_prior_date <- function(x, label) {
    if (inherits(x, "Date")) {
        if (any(!is.na(x) & !is.finite(as.numeric(x)))) 
            stop(label, ": nonfinite date")
        return(x)
    }
    if (inherits(x, "POSIXt")) 
        return(as.Date(x))
    z <- trimws(as.character(x))
    missing <- is.na(x) | toupper(z) %in% c("", "NA", "NAN", "NULL")
    good <- missing | grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z)
    if (any(!good)) 
        stop(label, ": unverified non-ISO date precision")
    out <- rep(as.Date(NA), length(z))
    out[!missing] <- as.Date(z[!missing], format = "%Y-%m-%d")
    if (any(!missing & is.na(out))) 
        stop(label, ": invalid date")
    out
}

ev_diabetes_prior_sources <- function(all, eids, entry, birth_year = NULL, contract = ev_diabetes_prior_contract_example(), 
    specs = ev_diabetes_prior_specs()) {
    all <- as.data.frame(all, stringsAsFactors = FALSE)
    if (!"eid" %in% names(all) || anyNA(all$eid) || anyDuplicated(all$eid) || anyDuplicated(names(all))) 
        stop("Supplementary DM source requires unique nonmissing eid/column names")
    if (anyNA(eids) || anyDuplicated(eids)) 
        stop("Requested supplementary DM eid invalid")
    n <- length(eids)
    if (length(entry) != n) 
        stop("Supplementary DM baseline date length mismatch")
    entry <- ev_diabetes_prior_date(entry, "baseline entry")
    if (!all(c("column", "source", "precision") %in% names(specs)) || anyNA(specs$column) || anyDuplicated(specs$column) || 
        !nrow(specs)) 
        stop("Invalid reviewed supplementary DM specs")
    approved <- ev_diabetes_prior_specs()
    if (any(!specs$column %in% approved$column)) 
        stop("Unreviewed DM source column; do not invent mapping")
    if (!identical(as.character(specs$precision), approved$precision[match(specs$column, approved$column)])) 
        stop("Supplementary DM source precision differs from reviewed construction")
    evidence <- if (is.null(contract$evidence)) 
        ""
    else as.character(contract$evidence)
    present <- specs$column %in% names(all)
    gated <- isTRUE(contract$verified) && isTRUE(contract$baseline_comparison_verified) && length(evidence) == 
        1L && !is.na(evidence) && nzchar(trimws(evidence)) && all(present) && setequal(as.character(contract$source_columns), 
        as.character(specs$column))
    if (!gated) 
        return(list(verified = FALSE, prior = rep(FALSE, n), timing_uncertain = rep(TRUE, n), exclude = rep(TRUE, 
            n), audit = data.frame(column = specs$column, present = present), scope = "BLOCKED supplementary DM provenance; no missing date was interpreted as negative"))
    ix <- match(as.character(eids), as.character(all$eid))
    linked <- !is.na(ix)
    if (is.null(birth_year) && "birth_year" %in% names(all)) 
        birth_year <- all$birth_year[ix]
    if (is.null(birth_year)) 
        birth_year <- rep(NA_real_, n)
    if (length(birth_year) != n) 
        stop("Supplementary DM birth-year length mismatch")
    yy <- suppressWarnings(as.numeric(as.character(birth_year)))
    valid_yob <- !is.na(yy) & is.finite(yy) & yy == floor(yy) & yy >= 1900 & yy <= 2100
    earliest_birth <- rep(as.Date(NA), n)
    earliest_birth[valid_yob] <- as.Date(paste0(as.integer(yy[valid_yob]), "-01-01"))
    pseudo <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
    prior <- rep(FALSE, n)
    uncertain <- !linked | is.na(entry)
    audit <- list()
    source_flags <- list()
    for (j in seq_len(nrow(specs))) {
        col <- specs$column[j]
        d <- ev_diabetes_prior_date(all[[col]][ix], col)
        exists <- !is.na(d)
        special <- exists & d %in% pseudo
        prebirth <- exists & !is.na(earliest_birth) & d < earliest_birth
        uncertain_j <- special | prebirth
        regular <- exists & !special & !prebirth & !is.na(entry)
        if (specs$precision[j] == "year") {
            yr <- as.integer(format(d, "%Y"))
            byr <- as.integer(format(entry, "%Y"))
            prior_j <- regular & !is.na(yr) & !is.na(byr) & yr < byr
            uncertain_j <- uncertain_j | (regular & !is.na(yr) & !is.na(byr) & yr == byr)
            after_j <- regular & !is.na(yr) & !is.na(byr) & yr > byr
        }
        else {
            prior_j <- regular & d <= entry
            after_j <- regular & d > entry
        }
        prior_j[is.na(prior_j)] <- FALSE
        uncertain_j[is.na(uncertain_j)] <- TRUE
        after_j[is.na(after_j)] <- FALSE
        prior <- prior | prior_j
        uncertain <- uncertain | uncertain_j
        source_flags[[col]] <- list(prior = prior_j, timing_uncertain = uncertain_j, date_recorded = exists, after_baseline_record = after_j)
        audit[[j]] <- data.frame(column = col, source = specs$source[j], precision = specs$precision[j], n_requested = n, 
            n_linked = sum(linked), n_date_recorded = sum(exists), n_prior = sum(prior_j), n_timing_uncertain = sum(uncertain_j), 
            n_pseudo_date = sum(special), n_prebirth_year = sum(prebirth), n_after_baseline_record = sum(after_j), 
            n_date_missing = sum(!exists), missing_date_is_negative = FALSE, stringsAsFactors = FALSE)
    }
    list(verified = TRUE, eid = eids, prior = prior, timing_uncertain = uncertain, exclude = prior | uncertain, 
        source_flags = source_flags, audit = do.call(rbind, audit), scope = paste("Positive-only supplementary exclusions from locally mapped GP,", 
            "year-resolution SRD and UKB E11 first-occurrence records. Missing dates", "do not establish absence or source coverage; no diabetes type inferred;", 
            "SRD source may include retrospective later-visit reports of earlier onset.", "Not all-source incidence or the original Eastwood incident algorithm."))
}

ev_diabetes_prior_tests <- function() {
    a <- data.frame(eid = paste0("synthetic_", 1:8), birth_year = 1950, fod_gp2_dm = as.Date(c("2005-01-01", NA, 
        "1902-02-02", NA, NA, NA, NA, NA)), fod_gp3_dm = as.Date(c(NA, "2015-01-01", "2007-01-01", NA, NA, NA, 
        NA, NA)), fod_srd_dm = as.Date(c(NA, NA, NA, "2008-07-01", "2007-07-01", "2009-07-01", "1940-07-01", NA)), 
        fod_ref_t2dm = as.Date(c(NA, NA, NA, NA, NA, NA, NA, "2037-07-07")))
    s <- ev_diabetes_prior_specs()
    s <- s[s$column %in% setdiff(names(a), c("eid", "birth_year")), ]
    c <- list(verified = TRUE, evidence = "synthetic only", source_columns = s$column, baseline_comparison_verified = TRUE)
    original <- a
    entry <- rep(as.Date("2008-05-01"), 8L)
    z <- ev_diabetes_prior_sources(a, a$eid, entry, contract = c, specs = s)
    stopifnot(z$verified, identical(a, original), identical(z$prior, c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, 
        FALSE)), identical(z$timing_uncertain, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)), !z$exclude[2], 
        !z$exclude[6], all(!z$audit$missing_date_is_negative))
    k <- ev_diabetes_prior_sources(a, c(a$eid[1], "synthetic_absent"), entry[1:2], contract = c, specs = s)
    stopifnot(k$prior[1], !k$prior[2], k$timing_uncertain[2], k$exclude[2])
    stopifnot(!ev_diabetes_prior_sources(a, a$eid, entry, specs = s)$verified)
    invisible(TRUE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ev_diabetes_prior_specs <- function() data.frame(column = c("fod_gp2_dm", "fod_gp2_t2dm", "fod_gp3_dm", 
            "fod_gp3_t2dm", "fod_gp_dm", "fod_gp_t2dm", "fod_srd_dm", "fod_srd_t2dm", "fod_ref_t2dm"), source = c(rep("local_GP_mapping_positive_only", 
            6L), rep("SRD_year_resolution_positive_only", 2L), "UKB130708_first_occurrence_E11"), precision = c(rep("day_or_placeholder", 
            6L), rep("year", 2L), "reported_date_mixed_source"), stringsAsFactors = FALSE)
        ev_diabetes_prior_contract_example <- function() list(verified = FALSE, evidence = "", source_columns = character(), 
            source_path = "", baseline_comparison_verified = FALSE)
        ev_diabetes_prior_date <- function(x, label) {
            if (inherits(x, "Date")) {
                if (any(!is.na(x) & !is.finite(as.numeric(x)))) stop(label, ": nonfinite date")
                return(x)
            }
            if (inherits(x, "POSIXt")) return(as.Date(x))
            z <- trimws(as.character(x))
            missing <- is.na(x) | toupper(z) %in% c("", "NA", "NAN", "NULL")
            good <- missing | grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z)
            if (any(!good)) stop(label, ": unverified non-ISO date precision")
            out <- rep(as.Date(NA), length(z))
            out[!missing] <- as.Date(z[!missing], format = "%Y-%m-%d")
            if (any(!missing & is.na(out))) stop(label, ": invalid date")
            out
        }
        ev_diabetes_prior_sources <- function(all, eids, entry, birth_year = NULL, contract = ev_diabetes_prior_contract_example(), 
            specs = ev_diabetes_prior_specs()) {
            all <- as.data.frame(all, stringsAsFactors = FALSE)
            if (!"eid" %in% names(all) || anyNA(all$eid) || anyDuplicated(all$eid) || anyDuplicated(names(all))) stop("Supplementary DM source requires unique nonmissing eid/column names")
            if (anyNA(eids) || anyDuplicated(eids)) stop("Requested supplementary DM eid invalid")
            n <- length(eids)
            if (length(entry) != n) stop("Supplementary DM baseline date length mismatch")
            entry <- ev_diabetes_prior_date(entry, "baseline entry")
            if (!all(c("column", "source", "precision") %in% names(specs)) || anyNA(specs$column) || anyDuplicated(specs$column) || 
                !nrow(specs)) stop("Invalid reviewed supplementary DM specs")
            approved <- ev_diabetes_prior_specs()
            if (any(!specs$column %in% approved$column)) stop("Unreviewed DM source column; do not invent mapping")
            if (!identical(as.character(specs$precision), approved$precision[match(specs$column, approved$column)])) stop("Supplementary DM source precision differs from reviewed construction")
            evidence <- if (is.null(contract$evidence)) "" else as.character(contract$evidence)
            present <- specs$column %in% names(all)
            gated <- isTRUE(contract$verified) && isTRUE(contract$baseline_comparison_verified) && length(evidence) == 
                1L && !is.na(evidence) && nzchar(trimws(evidence)) && all(present) && setequal(as.character(contract$source_columns), 
                as.character(specs$column))
            if (!gated) return(list(verified = FALSE, prior = rep(FALSE, n), timing_uncertain = rep(TRUE, n), exclude = rep(TRUE, 
                n), audit = data.frame(column = specs$column, present = present), scope = "BLOCKED supplementary DM provenance; no missing date was interpreted as negative"))
            ix <- match(as.character(eids), as.character(all$eid))
            linked <- !is.na(ix)
            if (is.null(birth_year) && "birth_year" %in% names(all)) birth_year <- all$birth_year[ix]
            if (is.null(birth_year)) birth_year <- rep(NA_real_, n)
            if (length(birth_year) != n) stop("Supplementary DM birth-year length mismatch")
            yy <- suppressWarnings(as.numeric(as.character(birth_year)))
            valid_yob <- !is.na(yy) & is.finite(yy) & yy == floor(yy) & yy >= 1900 & yy <= 2100
            earliest_birth <- rep(as.Date(NA), n)
            earliest_birth[valid_yob] <- as.Date(paste0(as.integer(yy[valid_yob]), "-01-01"))
            pseudo <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
            prior <- rep(FALSE, n)
            uncertain <- !linked | is.na(entry)
            audit <- list()
            source_flags <- list()
            for (j in seq_len(nrow(specs))) {
                col <- specs$column[j]
                d <- ev_diabetes_prior_date(all[[col]][ix], col)
                exists <- !is.na(d)
                special <- exists & d %in% pseudo
                prebirth <- exists & !is.na(earliest_birth) & d < earliest_birth
                uncertain_j <- special | prebirth
                regular <- exists & !special & !prebirth & !is.na(entry)
                if (specs$precision[j] == "year") {
                  yr <- as.integer(format(d, "%Y"))
                  byr <- as.integer(format(entry, "%Y"))
                  prior_j <- regular & !is.na(yr) & !is.na(byr) & yr < byr
                  uncertain_j <- uncertain_j | (regular & !is.na(yr) & !is.na(byr) & yr == byr)
                  after_j <- regular & !is.na(yr) & !is.na(byr) & yr > byr
                } else {
                  prior_j <- regular & d <= entry
                  after_j <- regular & d > entry
                }
                prior_j[is.na(prior_j)] <- FALSE
                uncertain_j[is.na(uncertain_j)] <- TRUE
                after_j[is.na(after_j)] <- FALSE
                prior <- prior | prior_j
                uncertain <- uncertain | uncertain_j
                source_flags[[col]] <- list(prior = prior_j, timing_uncertain = uncertain_j, date_recorded = exists, 
                  after_baseline_record = after_j)
                audit[[j]] <- data.frame(column = col, source = specs$source[j], precision = specs$precision[j], 
                  n_requested = n, n_linked = sum(linked), n_date_recorded = sum(exists), n_prior = sum(prior_j), 
                  n_timing_uncertain = sum(uncertain_j), n_pseudo_date = sum(special), n_prebirth_year = sum(prebirth), 
                  n_after_baseline_record = sum(after_j), n_date_missing = sum(!exists), missing_date_is_negative = FALSE, 
                  stringsAsFactors = FALSE)
            }
            list(verified = TRUE, eid = eids, prior = prior, timing_uncertain = uncertain, exclude = prior | uncertain, 
                source_flags = source_flags, audit = do.call(rbind, audit), scope = paste("Positive-only supplementary exclusions from locally mapped GP,", 
                  "year-resolution SRD and UKB E11 first-occurrence records. Missing dates", "do not establish absence or source coverage; no diabetes type inferred;", 
                  "SRD source may include retrospective later-visit reports of earlier onset.", "Not all-source incidence or the original Eastwood incident algorithm."))
        }
        ev_diabetes_prior_tests <- function() {
            a <- data.frame(eid = paste0("synthetic_", 1:8), birth_year = 1950, fod_gp2_dm = as.Date(c("2005-01-01", 
                NA, "1902-02-02", NA, NA, NA, NA, NA)), fod_gp3_dm = as.Date(c(NA, "2015-01-01", "2007-01-01", 
                NA, NA, NA, NA, NA)), fod_srd_dm = as.Date(c(NA, NA, NA, "2008-07-01", "2007-07-01", "2009-07-01", 
                "1940-07-01", NA)), fod_ref_t2dm = as.Date(c(NA, NA, NA, NA, NA, NA, NA, "2037-07-07")))
            s <- ev_diabetes_prior_specs()
            s <- s[s$column %in% setdiff(names(a), c("eid", "birth_year")), ]
            c <- list(verified = TRUE, evidence = "synthetic only", source_columns = s$column, baseline_comparison_verified = TRUE)
            original <- a
            entry <- rep(as.Date("2008-05-01"), 8L)
            z <- ev_diabetes_prior_sources(a, a$eid, entry, contract = c, specs = s)
            stopifnot(z$verified, identical(a, original), identical(z$prior, c(TRUE, FALSE, TRUE, FALSE, TRUE, 
                FALSE, FALSE, FALSE)), identical(z$timing_uncertain, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, 
                TRUE, TRUE)), !z$exclude[2], !z$exclude[6], all(!z$audit$missing_date_is_negative))
            k <- ev_diabetes_prior_sources(a, c(a$eid[1], "synthetic_absent"), entry[1:2], contract = c, specs = s)
            stopifnot(k$prior[1], !k$prior[2], k$timing_uncertain[2], k$exclude[2])
            stopifnot(!ev_diabetes_prior_sources(a, a$eid, entry, specs = s)$verified)
            invisible(TRUE)
        }
    }), envir = context)
}

