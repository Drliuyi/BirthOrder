# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ev_operation_prior <- function(raw, eids) {
    raw <- ev_ids(raw, "baseline operation source")
    cols <- grep("^p20004_i0_a[0-9]+$", names(raw), value = TRUE)
    if (!length(cols) || !"p136_i0" %in% names(raw)) 
        stop("Baseline operation/count contract missing")
    ix <- match(eids, raw$eid)
    d <- raw[ix, , drop = FALSE]
    n <- length(eids)
    count <- ev_num(d$p136_i0)
    used <- integer(n)
    positive <- uncertain <- rep(FALSE, n)
    for (v in cols) {
        z <- ev_num(d[[v]])
        used <- used + as.integer(!is.na(z))
        positive <- positive | (!is.na(z) & z %in% c(1070, 1095, 1523))
        uncertain <- uncertain | (!is.na(z) & z %in% c(1069, 99999, -1, -3))
    }
    known <- !is.na(ix) & !is.na(count) & count >= 0 & count == floor(count) & count == used
    flag <- as.integer(positive)
    flag[(!known | uncertain) & !positive] <- NA_integer_
    data.frame(eid = eids, prior_ihd_operation_report = flag, operation_source_unknown = is.na(flag), operation_count_mismatch = !known, 
        stringsAsFactors = FALSE)
}

ev_operation_tests <- function() {
    x <- data.frame(eid = as.character(1:6), p136_i0 = c(0, 1, 1, 1, NA, 2), p20004_i0_a0 = c(NA, 1070, 1095, 1069, 
        NA, 1000), p20004_i0_a1 = NA)
    z <- ev_operation_prior(x, x$eid)
    stopifnot(identical(z$prior_ihd_operation_report, c(0L, 1L, 1L, NA_integer_, NA_integer_, NA_integer_)))
    TRUE
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ev_operation_prior <- function(raw, eids) {
            raw <- ev_ids(raw, "baseline operation source")
            cols <- grep("^p20004_i0_a[0-9]+$", names(raw), value = TRUE)
            if (!length(cols) || !"p136_i0" %in% names(raw)) stop("Baseline operation/count contract missing")
            ix <- match(eids, raw$eid)
            d <- raw[ix, , drop = FALSE]
            n <- length(eids)
            count <- ev_num(d$p136_i0)
            used <- integer(n)
            positive <- uncertain <- rep(FALSE, n)
            for (v in cols) {
                z <- ev_num(d[[v]])
                used <- used + as.integer(!is.na(z))
                positive <- positive | (!is.na(z) & z %in% c(1070, 1095, 1523))
                uncertain <- uncertain | (!is.na(z) & z %in% c(1069, 99999, -1, -3))
            }
            known <- !is.na(ix) & !is.na(count) & count >= 0 & count == floor(count) & count == used
            flag <- as.integer(positive)
            flag[(!known | uncertain) & !positive] <- NA_integer_
            data.frame(eid = eids, prior_ihd_operation_report = flag, operation_source_unknown = is.na(flag), operation_count_mismatch = !known, 
                stringsAsFactors = FALSE)
        }
        ev_operation_tests <- function() {
            x <- data.frame(eid = as.character(1:6), p136_i0 = c(0, 1, 1, 1, NA, 2), p20004_i0_a0 = c(NA, 1070, 
                1095, 1069, NA, 1000), p20004_i0_a1 = NA)
            z <- ev_operation_prior(x, x$eid)
            stopifnot(identical(z$prior_ihd_operation_report, c(0L, 1L, 1L, NA_integer_, NA_integer_, NA_integer_)))
            TRUE
        }
    }), envir = context)
}

