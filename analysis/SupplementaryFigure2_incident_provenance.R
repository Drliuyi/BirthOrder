# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ev_prov_missing <- function(x) {
    z <- as.character(x)
    is.na(x) | toupper(trimws(z)) %in% c("", "NA", "NAN", "NULL")
}

ev_prov_unique_eid <- function(x, label) {
    if (!is.data.frame(x) || !"eid" %in% names(x) || anyNA(x$eid) || anyDuplicated(as.character(x$eid))) 
        stop(label, ": missing/duplicate eid")
    invisible(TRUE)
}

ev_prov_date <- function(x, label) {
    if (inherits(x, "Date")) 
        return(as.Date(x))
    z <- as.character(x)
    missing <- ev_prov_missing(x)
    if (any(!missing & !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z))) 
        stop(label, ": unverified date precision/encoding")
    out <- rep(as.Date(NA), length(z))
    out[!missing] <- as.Date(z[!missing])
    if (any(!missing & is.na(out))) 
        stop(label, ": invalid date")
    out
}

ev_prov_pseudodate <- function(x) {
    as.character(x) %in% c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07")
}

ev_prov_registry_dates <- function(full_raw, eids) {
    ev_prov_unique_eid(full_raw, "P0 full raw")
    required <- c("p40000_i0", "p40000_i1", "p191")
    if (!all(required %in% names(full_raw))) 
        stop("P0 lacks uncapped death/loss fields: ", paste(setdiff(required, names(full_raw)), collapse = ","))
    ix <- match(as.character(eids), as.character(full_raw$eid))
    if (anyNA(ix)) 
        stop("Some target eids absent from complete raw census")
    a <- ev_prov_date(full_raw$p40000_i0[ix], "death instance0")
    b <- ev_prov_date(full_raw$p40000_i1[ix], "death instance1")
    if (any(!is.na(a) & !is.na(b) & a != b)) 
        stop("Conflicting uncapped death instances: reconcile before modeling")
    a[is.na(a)] <- b[is.na(a)]
    data.frame(eid = as.character(eids), date_death_uncapped = a, date_lost_registry = ev_prov_date(full_raw$p191[ix], 
        "recorded loss"), registry_source_row_present = TRUE, stringsAsFactors = FALSE)
}

ev_prov_uniform_coverage <- function(full_raw, code_roster, date_roster) {
    ev_prov_unique_eid(full_raw, "P0 full raw")
    ev_prov_unique_eid(code_roster, "hospital code roster")
    ev_prov_unique_eid(date_roster, "hospital date roster")
    if (!identical(as.character(code_roster$eid), as.character(date_roster$eid))) 
        stop("Hospital code/date rosters differ")
    if (!setequal(as.character(full_raw$eid), as.character(code_roster$eid))) 
        stop("Full raw and hospital linkage-population rosters differ")
    out <- data.frame(eid = as.character(full_raw$eid), eligible = TRUE, hospital_coverage_group = "UKB_institutional_uniform", 
        stringsAsFactors = FALSE)
    attr(out, "linkage_model") <- "Institutional full-cohort linkage, no person-level success certificate"
    attr(out, "cutoff") <- "2022-05-31"
    attr(out, "release_limitation") <- "Original UKB release ID unrecorded; identify local snapshot by P0 SHA256"
    out
}

ev_prov_opcs_coronary_codes <- function() {
    suffix <- list(K40 = c(1, 2, 3, 4, 8, 9), K41 = c(1, 2, 3, 4, 8, 9), K42 = c(1, 2, 3, 4, 8, 9), K43 = c(1, 
        2, 3, 4, 8, 9), K44 = c(1, 2, 8, 9), K45 = c(1, 2, 3, 4, 5, 6, 8, 9), K46 = c(1, 2, 3, 4, 5, 8, 9), K49 = c(1, 
        2, 3, 4, 8, 9), K50 = c(1, 2, 3, 4, 8, 9), K75 = c(1, 2, 3, 4, 8, 9))
    unname(unlist(Map(function(prefix, suffixes) paste0(prefix, suffixes), names(suffix), suffix), use.names = FALSE))
}

ev_prov_icd9_prior <- function(indir, eids, entry) {
    if (!exists("ev_inc_read_pair", mode = "function") || !exists("ev_inc_scan", mode = "function")) 
        stop("Source incident.R before provenance helpers")
    pair <- ev_inc_read_pair(file.path(indir, "Rdata", "icd9.code0.rds"), file.path(indir, "Rdata", "icd9.date0.rds"), 
        eids)
    scan <- ev_inc_scan(pair, "^41[0-4][0-9A-Z]*$", entry)
    flag <- as.integer(scan$prior)
    flag[!scan$known] <- NA_integer_
    data.frame(eid = as.character(eids), prior_ihd_icd9 = flag, prior_ihd_icd9_undated = scan$undated, stringsAsFactors = FALSE)
}

ev_prov_firstocc_prior <- function(firstocc, eids, entry) {
    ev_prov_unique_eid(firstocc, "UKB first occurrence")
    ix <- match(as.character(eids), as.character(firstocc$eid))
    prefixes <- paste0("p", c(131296, 131298, 131300, 131302, 131304, 131306))
    cols <- vapply(prefixes, function(prefix) {
        found <- grep(paste0("^", prefix, "(_i0)?$"), names(firstocc), value = TRUE)
        if (length(found) != 1L) 
            stop("First-occurrence field absent/ambiguous: ", prefix)
        found
    }, character(1))
    prior <- pseudo <- rep(FALSE, length(eids))
    known <- !is.na(ix) & !is.na(entry)
    for (col in cols) {
        date <- ev_prov_date(firstocc[[col]][ix], col)
        bad_date <- ev_prov_pseudodate(date)
        pseudo <- pseudo | bad_date
        date[bad_date] <- as.Date(NA)
        prior <- prior | (!is.na(date) & !is.na(entry) & date <= entry)
    }
    flag <- as.integer(prior)
    flag[!known | (pseudo & !prior)] <- NA_integer_
    data.frame(eid = as.character(eids), prior_ihd_ukb_firstocc = flag, firstocc_ihd_date_unknown = pseudo, stringsAsFactors = FALSE)
}

ev_prov_gp_positive_prior <- function(gp, eids, entry) {
    ev_prov_unique_eid(gp, "partial GP history")
    if (!"fod_gp_cvd_cad" %in% names(gp)) 
        stop("Audited GP CAD date absent")
    ix <- match(as.character(eids), as.character(gp$eid))
    date <- ev_prov_date(gp$fod_gp_cvd_cad[ix], "GP CAD first-record date")
    pseudo <- ev_prov_pseudodate(date)
    date[pseudo] <- as.Date(NA)
    flag <- as.integer(!is.na(date) & !is.na(entry) & date <= entry)
    flag[is.na(entry) | pseudo] <- NA_integer_
    data.frame(eid = as.character(eids), prior_ihd_gp_available_positive = flag, gp_history_source_row_present = !is.na(ix), 
        gp_ihd_date_unknown = pseudo, stringsAsFactors = FALSE)
}

ev_prov_touchscreen_prior <- function(baseline_ts, eids) {
    ev_prov_unique_eid(baseline_ts, "baseline touchscreen cardiovascular history")
    if (!"p6150_i0" %in% names(baseline_ts)) 
        stop("Require intact RAW p6150_i0; transformed cvd_dr lost none-vs-unknown")
    ix <- match(as.character(eids), as.character(baseline_ts$eid))
    raw <- baseline_ts$p6150_i0[ix]
    flag <- rep(NA_integer_, length(eids))
    available <- !is.na(ix) & !ev_prov_missing(raw)
    for (i in which(available)) {
        vals <- trimws(strsplit(as.character(raw[i]), "|", fixed = TRUE)[[1]])
        if (any(!vals %in% c("1", "2", "3", "4", "-7", "-3"))) 
            stop("Unparsed baseline field6150 response")
        if (any(vals %in% c("-7", "-3")) && length(vals) != 1L) 
            stop("Exclusive none/PNA field6150 code combined with other responses")
        if (identical(vals, "-3")) 
            next
        flag[i] <- as.integer(any(vals %in% c("1", "2")))
    }
    data.frame(eid = as.character(eids), prior_ihd_touchscreen = flag, touchscreen_ihd_unknown = is.na(flag), stringsAsFactors = FALSE)
}

ev_prov_operation_prior <- function(baseline_ops, eids) {
    ev_prov_unique_eid(baseline_ops, "baseline self-reported operations")
    cols <- paste0("p20004_i0_a", 0:31)
    required <- c("p136_i0", cols)
    if (!all(required %in% names(baseline_ops))) 
        stop("Require intact baseline operation count and all32 instance0 code positions")
    ix <- match(as.character(eids), as.character(baseline_ops$eid))
    n <- suppressWarnings(as.numeric(as.character(baseline_ops$p136_i0[ix])))
    count_valid <- !is.na(n) & n >= 0 & n <= 32 & n == floor(n)
    flag <- rep(NA_integer_, length(eids))
    for (i in which(!is.na(ix))) {
        raw <- unlist(baseline_ops[ix[i], cols, drop = FALSE], use.names = FALSE)
        raw <- raw[!ev_prov_missing(raw)]
        codes <- suppressWarnings(as.numeric(as.character(raw)))
        if (anyNA(codes)) 
            stop("Non-numeric baseline operation code")
        if (any(codes %in% c(1070, 1095, 1523))) {
            flag[i] <- 1L
            next
        }
        uncertain <- any(codes %in% c(1069, 99999) | codes < 0)
        if (count_valid[i] && n[i] == length(codes) && !uncertain) 
            flag[i] <- 0L
    }
    data.frame(eid = as.character(eids), prior_ihd_baseline_operation = flag, baseline_operation_ihd_unknown = is.na(flag), 
        stringsAsFactors = FALSE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ev_prov_missing <- function(x) {
            z <- as.character(x)
            is.na(x) | toupper(trimws(z)) %in% c("", "NA", "NAN", "NULL")
        }
        ev_prov_unique_eid <- function(x, label) {
            if (!is.data.frame(x) || !"eid" %in% names(x) || anyNA(x$eid) || anyDuplicated(as.character(x$eid))) stop(label, 
                ": missing/duplicate eid")
            invisible(TRUE)
        }
        ev_prov_date <- function(x, label) {
            if (inherits(x, "Date")) return(as.Date(x))
            z <- as.character(x)
            missing <- ev_prov_missing(x)
            if (any(!missing & !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z))) stop(label, ": unverified date precision/encoding")
            out <- rep(as.Date(NA), length(z))
            out[!missing] <- as.Date(z[!missing])
            if (any(!missing & is.na(out))) stop(label, ": invalid date")
            out
        }
        ev_prov_pseudodate <- function(x) {
            as.character(x) %in% c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07")
        }
        ev_prov_registry_dates <- function(full_raw, eids) {
            ev_prov_unique_eid(full_raw, "P0 full raw")
            required <- c("p40000_i0", "p40000_i1", "p191")
            if (!all(required %in% names(full_raw))) stop("P0 lacks uncapped death/loss fields: ", paste(setdiff(required, 
                names(full_raw)), collapse = ","))
            ix <- match(as.character(eids), as.character(full_raw$eid))
            if (anyNA(ix)) stop("Some target eids absent from complete raw census")
            a <- ev_prov_date(full_raw$p40000_i0[ix], "death instance0")
            b <- ev_prov_date(full_raw$p40000_i1[ix], "death instance1")
            if (any(!is.na(a) & !is.na(b) & a != b)) stop("Conflicting uncapped death instances: reconcile before modeling")
            a[is.na(a)] <- b[is.na(a)]
            data.frame(eid = as.character(eids), date_death_uncapped = a, date_lost_registry = ev_prov_date(full_raw$p191[ix], 
                "recorded loss"), registry_source_row_present = TRUE, stringsAsFactors = FALSE)
        }
        ev_prov_uniform_coverage <- function(full_raw, code_roster, date_roster) {
            ev_prov_unique_eid(full_raw, "P0 full raw")
            ev_prov_unique_eid(code_roster, "hospital code roster")
            ev_prov_unique_eid(date_roster, "hospital date roster")
            if (!identical(as.character(code_roster$eid), as.character(date_roster$eid))) stop("Hospital code/date rosters differ")
            if (!setequal(as.character(full_raw$eid), as.character(code_roster$eid))) stop("Full raw and hospital linkage-population rosters differ")
            out <- data.frame(eid = as.character(full_raw$eid), eligible = TRUE, hospital_coverage_group = "UKB_institutional_uniform", 
                stringsAsFactors = FALSE)
            attr(out, "linkage_model") <- "Institutional full-cohort linkage, no person-level success certificate"
            attr(out, "cutoff") <- "2022-05-31"
            attr(out, "release_limitation") <- "Original UKB release ID unrecorded; identify local snapshot by P0 SHA256"
            out
        }
        ev_prov_opcs_coronary_codes <- function() {
            suffix <- list(K40 = c(1, 2, 3, 4, 8, 9), K41 = c(1, 2, 3, 4, 8, 9), K42 = c(1, 2, 3, 4, 8, 9), K43 = c(1, 
                2, 3, 4, 8, 9), K44 = c(1, 2, 8, 9), K45 = c(1, 2, 3, 4, 5, 6, 8, 9), K46 = c(1, 2, 3, 4, 5, 8, 
                9), K49 = c(1, 2, 3, 4, 8, 9), K50 = c(1, 2, 3, 4, 8, 9), K75 = c(1, 2, 3, 4, 8, 9))
            unname(unlist(Map(function(prefix, suffixes) paste0(prefix, suffixes), names(suffix), suffix), use.names = FALSE))
        }
        ev_prov_icd9_prior <- function(indir, eids, entry) {
            if (!exists("ev_inc_read_pair", mode = "function") || !exists("ev_inc_scan", mode = "function")) stop("Source incident.R before provenance helpers")
            pair <- ev_inc_read_pair(file.path(indir, "Rdata", "icd9.code0.rds"), file.path(indir, "Rdata", "icd9.date0.rds"), 
                eids)
            scan <- ev_inc_scan(pair, "^41[0-4][0-9A-Z]*$", entry)
            flag <- as.integer(scan$prior)
            flag[!scan$known] <- NA_integer_
            data.frame(eid = as.character(eids), prior_ihd_icd9 = flag, prior_ihd_icd9_undated = scan$undated, 
                stringsAsFactors = FALSE)
        }
        ev_prov_firstocc_prior <- function(firstocc, eids, entry) {
            ev_prov_unique_eid(firstocc, "UKB first occurrence")
            ix <- match(as.character(eids), as.character(firstocc$eid))
            prefixes <- paste0("p", c(131296, 131298, 131300, 131302, 131304, 131306))
            cols <- vapply(prefixes, function(prefix) {
                found <- grep(paste0("^", prefix, "(_i0)?$"), names(firstocc), value = TRUE)
                if (length(found) != 1L) stop("First-occurrence field absent/ambiguous: ", prefix)
                found
            }, character(1))
            prior <- pseudo <- rep(FALSE, length(eids))
            known <- !is.na(ix) & !is.na(entry)
            for (col in cols) {
                date <- ev_prov_date(firstocc[[col]][ix], col)
                bad_date <- ev_prov_pseudodate(date)
                pseudo <- pseudo | bad_date
                date[bad_date] <- as.Date(NA)
                prior <- prior | (!is.na(date) & !is.na(entry) & date <= entry)
            }
            flag <- as.integer(prior)
            flag[!known | (pseudo & !prior)] <- NA_integer_
            data.frame(eid = as.character(eids), prior_ihd_ukb_firstocc = flag, firstocc_ihd_date_unknown = pseudo, 
                stringsAsFactors = FALSE)
        }
        ev_prov_gp_positive_prior <- function(gp, eids, entry) {
            ev_prov_unique_eid(gp, "partial GP history")
            if (!"fod_gp_cvd_cad" %in% names(gp)) stop("Audited GP CAD date absent")
            ix <- match(as.character(eids), as.character(gp$eid))
            date <- ev_prov_date(gp$fod_gp_cvd_cad[ix], "GP CAD first-record date")
            pseudo <- ev_prov_pseudodate(date)
            date[pseudo] <- as.Date(NA)
            flag <- as.integer(!is.na(date) & !is.na(entry) & date <= entry)
            flag[is.na(entry) | pseudo] <- NA_integer_
            data.frame(eid = as.character(eids), prior_ihd_gp_available_positive = flag, gp_history_source_row_present = !is.na(ix), 
                gp_ihd_date_unknown = pseudo, stringsAsFactors = FALSE)
        }
        ev_prov_touchscreen_prior <- function(baseline_ts, eids) {
            ev_prov_unique_eid(baseline_ts, "baseline touchscreen cardiovascular history")
            if (!"p6150_i0" %in% names(baseline_ts)) stop("Require intact RAW p6150_i0; transformed cvd_dr lost none-vs-unknown")
            ix <- match(as.character(eids), as.character(baseline_ts$eid))
            raw <- baseline_ts$p6150_i0[ix]
            flag <- rep(NA_integer_, length(eids))
            available <- !is.na(ix) & !ev_prov_missing(raw)
            for (i in which(available)) {
                vals <- trimws(strsplit(as.character(raw[i]), "|", fixed = TRUE)[[1]])
                if (any(!vals %in% c("1", "2", "3", "4", "-7", "-3"))) stop("Unparsed baseline field6150 response")
                if (any(vals %in% c("-7", "-3")) && length(vals) != 1L) stop("Exclusive none/PNA field6150 code combined with other responses")
                if (identical(vals, "-3")) next
                flag[i] <- as.integer(any(vals %in% c("1", "2")))
            }
            data.frame(eid = as.character(eids), prior_ihd_touchscreen = flag, touchscreen_ihd_unknown = is.na(flag), 
                stringsAsFactors = FALSE)
        }
        ev_prov_operation_prior <- function(baseline_ops, eids) {
            ev_prov_unique_eid(baseline_ops, "baseline self-reported operations")
            cols <- paste0("p20004_i0_a", 0:31)
            required <- c("p136_i0", cols)
            if (!all(required %in% names(baseline_ops))) stop("Require intact baseline operation count and all32 instance0 code positions")
            ix <- match(as.character(eids), as.character(baseline_ops$eid))
            n <- suppressWarnings(as.numeric(as.character(baseline_ops$p136_i0[ix])))
            count_valid <- !is.na(n) & n >= 0 & n <= 32 & n == floor(n)
            flag <- rep(NA_integer_, length(eids))
            for (i in which(!is.na(ix))) {
                raw <- unlist(baseline_ops[ix[i], cols, drop = FALSE], use.names = FALSE)
                raw <- raw[!ev_prov_missing(raw)]
                codes <- suppressWarnings(as.numeric(as.character(raw)))
                if (anyNA(codes)) stop("Non-numeric baseline operation code")
                if (any(codes %in% c(1070, 1095, 1523))) {
                  flag[i] <- 1L
                  next
                }
                uncertain <- any(codes %in% c(1069, 99999) | codes < 0)
                if (count_valid[i] && n[i] == length(codes) && !uncertain) flag[i] <- 0L
            }
            data.frame(eid = as.character(eids), prior_ihd_baseline_operation = flag, baseline_operation_ihd_unknown = is.na(flag), 
                stringsAsFactors = FALSE)
        }
    }), envir = context)
}

