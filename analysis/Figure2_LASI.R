# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

S8_ENDPOINTS <- c("diabetes", "heart")

S8_EXPOSURES <- c("relative_rank", "absolute_rank", "firstlast")

s8_cap <- function(d, cohort) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    required <- c("id", "family_size", "birth_rank", "exposure_eligible")
    if (length(setdiff(required, names(d)))) 
        stop("Missing S8 input fields")
    if (anyNA(d$id) || any(!nzchar(as.character(d$id))) || anyDuplicated(d$id)) 
        stop("Invalid participant IDs")
    if (!is.logical(d$exposure_eligible) || anyNA(d$exposure_eligible)) 
        stop("Eligibility must be nonmissing logical")
    S <- ra_numeric(d$family_size)
    R <- ra_numeric(d$birth_rank)
    valid_s <- is.finite(S) & S >= 2 & S == floor(S)
    valid_r <- is.finite(R) & R >= 1 & R <= S & R == floor(R)
    if (cohort == "SHARE") {
        if (!"last" %in% names(d)) 
            stop("SHARE requires reviewed reported first/last")
        valid <- valid_s & d$last %in% 0:1
    }
    else valid <- valid_s & valid_r
    valid[is.na(valid)] <- FALSE
    d$eligible_before_s8 <- d$exposure_eligible & valid
    d$exposure_eligible <- d$eligible_before_s8 & S <= 8
    d$family_size <- S
    d$birth_rank <- R
    if (cohort != "SHARE") {
        rr <- ifelse(valid, (R - 1)/(S - 1), NA_real_)
        if ("relative_rank" %in% names(d)) {
            supplied <- ra_numeric(d$relative_rank)
            bad <- d$eligible_before_s8 & is.finite(supplied) & abs(supplied - rr) > 1e-10
            if (any(bad)) 
                stop("Stored relative position disagrees with R/S; no automatic repair")
        }
        d$relative_rank <- ifelse(d$exposure_eligible, rr, NA_real_)
        d$last <- ifelse(d$exposure_eligible & R == 1, 0L, ifelse(d$exposure_eligible & R == S, 1L, NA_integer_))
        d$absolute_f <- factor(R, levels = sort(unique(R[d$exposure_eligible])))
        if ("1" %in% levels(d$absolute_f)) 
            d$absolute_f <- stats::relevel(d$absolute_f, "1")
    }
    else {
        d$relative_rank <- NA_real_
        d$absolute_f <- factor(rep(NA_character_, nrow(d)))
        d$last[!d$exposure_eligible] <- NA_integer_
    }
    d$family_f <- factor(S)
    d
}

s8_sample <- function(d, cohort, endpoint, exposure, sample = "primary") {
    if (!endpoint %in% S8_ENDPOINTS || !exposure %in% S8_EXPOSURES) 
        stop("Out-of-scope endpoint/design")
    if (cohort == "SHARE" && exposure != "firstlast") 
        stop("SHARE has no numeric rank; do not invent one")
    y <- d[[endpoint]]
    if (any(!is.na(y) & !y %in% 0:1)) 
        stop("Unknown disease cannot become a control")
    direct <- cohort == "CLHLS" && sample != "all_answers_sensitivity"
    v <- switch(exposure, relative_rank = "relative_rank", absolute_rank = "absolute_f", firstlast = "last")
    ok <- rep(TRUE, nrow(d))
    flow <- list()
    add <- function(stage, condition, note = "") {
        condition[is.na(condition)] <- FALSE
        before <- sum(ok)
        ok <<- ok & condition
        flow[[length(flow) + 1L]] <<- data.frame(cohort, endpoint, exposure, sample, stage, n = sum(ok), excluded = before - 
            sum(ok), note, stringsAsFactors = FALSE)
    }
    add("prepared", rep(TRUE, nrow(d)))
    add("previous_valid_exposure", d$eligible_before_s8, "Frozen original QC; WLS fixed selected family member, no redraw")
    add("family_size_2_to_8", d$exposure_eligible, if (cohort %in% c("HRS", "KLoSA", "SHARE")) 
        "2<=S_proxy<=8; living-sibling proxy, not verified total birth-family S"
    else "2<=reported total children incl respondent<=8")
    if (direct) 
        add("rank_direct_response", d$direct_response %in% 1L, "CLHLS same direct-response primary as previous version")
    if (cohort == "CLHLS" && !direct) 
        add("known_proxy_status", !is.na(d$proxy_f))
    if (exposure == "firstlast") 
        add("first_or_last", !is.na(d$last), "Only this supplemental design excludes middle positions")
    add("known_endpoint", !is.na(y), "Missing or unknown is not disease-free")
    needed <- c("id", v, "family_size", ra_covariates(cohort, "M1", direct))
    if (length(setdiff(needed, names(d)))) 
        stop("Missing locked covariates for ", cohort)
    add("M0_M1_common_complete_cases", complete.cases(d[, needed, drop = FALSE]) & is.finite(d$age0), "Identical endpoint/design-specific M1-complete participants in M0 and M1")
    z <- d[ok, , drop = FALSE]
    z$y <- z[[endpoint]]
    for (k in names(z)) if (is.factor(z[[k]])) {
        z[[k]] <- droplevels(z[[k]])
        if (nlevels(z[[k]]) > 1L) 
            contrasts(z[[k]]) <- contr.treatment(nlevels(z[[k]]))
    }
    if ("1" %in% levels(z$absolute_f)) 
        z$absolute_f <- relevel(z$absolute_f, "1")
    if (nrow(z) && any(z$family_size < 2 | z$family_size > 8)) 
        stop("S8 sample invariant failed")
    list(data = z, flow = do.call(rbind, flow), direct = direct, exposure_var = v, structural_notes = "S8-only restriction; frozen covariates/cutpoints")
}

s8_adjust <- function(g) {
    required <- c("cohort", "endpoint", "model", "exposure", "sample", "status", "p_value")
    if (length(setdiff(required, names(g)))) 
        stop("Missing global-test columns")
    if (any(!g$endpoint %in% S8_ENDPOINTS)) 
        stop("Unplanned endpoint")
    g$p_holm_cohort2 <- g$p_holm_crosscohort_supplementary <- NA_real_
    g$failed_slots_counted_as_p1 <- TRUE
    g$test_family <- paste(g$cohort, g$model, g$exposure, g$sample, sep = "|")
    for (k in unique(g$test_family)) {
        ix <- which(g$test_family == k)
        if (length(ix) != 2L || !setequal(g$endpoint[ix], S8_ENDPOINTS)) 
            stop("Expected exactly two endpoint slots per cohort/design/model/sample")
        good <- g$status[ix] == "ok" & is.finite(g$p_value[ix])
        p <- rep(1, 2)
        p[good] <- g$p_value[ix][good]
        a <- p.adjust(p, "holm")
        g$p_holm_cohort2[ix[good]] <- a[good]
    }
    cross <- paste(g$model, g$exposure, g$sample, sep = "|")
    for (k in unique(cross)) {
        ix <- which(cross == k)
        if (anyDuplicated(paste(g$cohort[ix], g$endpoint[ix]))) 
            stop("Duplicated cross-cohort slots")
        good <- g$status[ix] == "ok" & is.finite(g$p_value[ix])
        p <- rep(1, length(ix))
        p[good] <- g$p_value[ix][good]
        a <- p.adjust(p, "holm")
        g$p_holm_crosscohort_supplementary[ix[good]] <- a[good]
    }
    g$analysis_role <- ifelse(g$sample != "primary" | g$model != "M1", "sensitivity", ifelse(g$exposure == "relative_rank", 
        "primary_relative_position", "supplemental_design"))
    g
}

s8_bind <- function(x) {
    if (!length(x)) 
        return(data.frame())
    cols <- unique(unlist(lapply(x, names)))
    do.call(rbind, lapply(x, function(d) {
        for (k in setdiff(cols, names(d))) d[[k]] <- NA
        d[, cols, drop = FALSE]
    }))
}

s8_release <- function(r, cells = NULL) {
    for (k in c("pure_factor_levels", "unused_factor_levels", "family_min", "family_max", "birth_rank_max")) r[[k]] <- NULL
    r$family_restriction <- "2<=S<=8; proxy cohorts explicitly labelled S_proxy"
    sensitive <- intersect(c("n", "cases", "noncases", "beta", "OR", "CI_lower", "CI_upper", "p_value", "p_holm_cohort2", 
        "p_holm_crosscohort_supplementary", "statistic", "full_penalized_loglik", "restricted_penalized_loglik"), 
        names(r))
    for (i in seq_len(nrow(r))) {
        counts <- unlist(r[i, intersect(c("n", "cases", "noncases"), names(r)), drop = FALSE], use.names = FALSE)
        if (any(is.finite(counts) & counts <= 5)) {
            r[i, sensitive] <- NA
            r$status[i] <- "SUPPRESSED_SMALL_TOTAL"
        }
    }
    if (!is.null(cells) && "rank" %in% names(r) && any(r$exposure == "absolute_rank")) {
        risky <- as.character(cells$group[cells$n <= 5 | cells$cases <= 5 | cells$noncases <= 5])
        reference_risky <- "1" %in% risky
        if (length(risky) == 1L && nrow(cells) > 1L) {
            other <- which(!as.character(cells$group) %in% risky)
            if (length(other)) 
                risky <- c(risky, as.character(cells$group[other[which.min(cells$n[other])]]))
        }
        ix <- which(r$exposure == "absolute_rank" & (reference_risky | as.character(r$rank) %in% risky))
        if (length(ix)) {
            fields <- intersect(c("beta", "OR", "CI_lower", "CI_upper", "p_value"), names(r))
            r[ix, fields] <- NA
            r$status[ix] <- "SUPPRESSED_SMALL_CELL"
        }
    }
    r
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        S8_ENDPOINTS <- c("diabetes", "heart")
        S8_EXPOSURES <- c("relative_rank", "absolute_rank", "firstlast")
        s8_cap <- function(d, cohort) {
            d <- as.data.frame(d, stringsAsFactors = FALSE)
            required <- c("id", "family_size", "birth_rank", "exposure_eligible")
            if (length(setdiff(required, names(d)))) stop("Missing S8 input fields")
            if (anyNA(d$id) || any(!nzchar(as.character(d$id))) || anyDuplicated(d$id)) stop("Invalid participant IDs")
            if (!is.logical(d$exposure_eligible) || anyNA(d$exposure_eligible)) stop("Eligibility must be nonmissing logical")
            S <- ra_numeric(d$family_size)
            R <- ra_numeric(d$birth_rank)
            valid_s <- is.finite(S) & S >= 2 & S == floor(S)
            valid_r <- is.finite(R) & R >= 1 & R <= S & R == floor(R)
            if (cohort == "SHARE") {
                if (!"last" %in% names(d)) stop("SHARE requires reviewed reported first/last")
                valid <- valid_s & d$last %in% 0:1
            } else valid <- valid_s & valid_r
            valid[is.na(valid)] <- FALSE
            d$eligible_before_s8 <- d$exposure_eligible & valid
            d$exposure_eligible <- d$eligible_before_s8 & S <= 8
            d$family_size <- S
            d$birth_rank <- R
            if (cohort != "SHARE") {
                rr <- ifelse(valid, (R - 1)/(S - 1), NA_real_)
                if ("relative_rank" %in% names(d)) {
                  supplied <- ra_numeric(d$relative_rank)
                  bad <- d$eligible_before_s8 & is.finite(supplied) & abs(supplied - rr) > 1e-10
                  if (any(bad)) stop("Stored relative position disagrees with R/S; no automatic repair")
                }
                d$relative_rank <- ifelse(d$exposure_eligible, rr, NA_real_)
                d$last <- ifelse(d$exposure_eligible & R == 1, 0L, ifelse(d$exposure_eligible & R == S, 1L, NA_integer_))
                d$absolute_f <- factor(R, levels = sort(unique(R[d$exposure_eligible])))
                if ("1" %in% levels(d$absolute_f)) d$absolute_f <- stats::relevel(d$absolute_f, "1")
            } else {
                d$relative_rank <- NA_real_
                d$absolute_f <- factor(rep(NA_character_, nrow(d)))
                d$last[!d$exposure_eligible] <- NA_integer_
            }
            d$family_f <- factor(S)
            d
        }
        s8_sample <- function(d, cohort, endpoint, exposure, sample = "primary") {
            if (!endpoint %in% S8_ENDPOINTS || !exposure %in% S8_EXPOSURES) stop("Out-of-scope endpoint/design")
            if (cohort == "SHARE" && exposure != "firstlast") stop("SHARE has no numeric rank; do not invent one")
            y <- d[[endpoint]]
            if (any(!is.na(y) & !y %in% 0:1)) stop("Unknown disease cannot become a control")
            direct <- cohort == "CLHLS" && sample != "all_answers_sensitivity"
            v <- switch(exposure, relative_rank = "relative_rank", absolute_rank = "absolute_f", firstlast = "last")
            ok <- rep(TRUE, nrow(d))
            flow <- list()
            add <- function(stage, condition, note = "") {
                condition[is.na(condition)] <- FALSE
                before <- sum(ok)
                ok <<- ok & condition
                flow[[length(flow) + 1L]] <<- data.frame(cohort, endpoint, exposure, sample, stage, n = sum(ok), 
                  excluded = before - sum(ok), note, stringsAsFactors = FALSE)
            }
            add("prepared", rep(TRUE, nrow(d)))
            add("previous_valid_exposure", d$eligible_before_s8, "Frozen original QC; WLS fixed selected family member, no redraw")
            add("family_size_2_to_8", d$exposure_eligible, if (cohort %in% c("HRS", "KLoSA", "SHARE")) "2<=S_proxy<=8; living-sibling proxy, not verified total birth-family S" else "2<=reported total children incl respondent<=8")
            if (direct) add("rank_direct_response", d$direct_response %in% 1L, "CLHLS same direct-response primary as previous version")
            if (cohort == "CLHLS" && !direct) add("known_proxy_status", !is.na(d$proxy_f))
            if (exposure == "firstlast") add("first_or_last", !is.na(d$last), "Only this supplemental design excludes middle positions")
            add("known_endpoint", !is.na(y), "Missing or unknown is not disease-free")
            needed <- c("id", v, "family_size", ra_covariates(cohort, "M1", direct))
            if (length(setdiff(needed, names(d)))) stop("Missing locked covariates for ", cohort)
            add("M0_M1_common_complete_cases", complete.cases(d[, needed, drop = FALSE]) & is.finite(d$age0), "Identical endpoint/design-specific M1-complete participants in M0 and M1")
            z <- d[ok, , drop = FALSE]
            z$y <- z[[endpoint]]
            for (k in names(z)) if (is.factor(z[[k]])) {
                z[[k]] <- droplevels(z[[k]])
                if (nlevels(z[[k]]) > 1L) contrasts(z[[k]]) <- contr.treatment(nlevels(z[[k]]))
            }
            if ("1" %in% levels(z$absolute_f)) z$absolute_f <- relevel(z$absolute_f, "1")
            if (nrow(z) && any(z$family_size < 2 | z$family_size > 8)) stop("S8 sample invariant failed")
            list(data = z, flow = do.call(rbind, flow), direct = direct, exposure_var = v, structural_notes = "S8-only restriction; frozen covariates/cutpoints")
        }
        s8_adjust <- function(g) {
            required <- c("cohort", "endpoint", "model", "exposure", "sample", "status", "p_value")
            if (length(setdiff(required, names(g)))) stop("Missing global-test columns")
            if (any(!g$endpoint %in% S8_ENDPOINTS)) stop("Unplanned endpoint")
            g$p_holm_cohort2 <- g$p_holm_crosscohort_supplementary <- NA_real_
            g$failed_slots_counted_as_p1 <- TRUE
            g$test_family <- paste(g$cohort, g$model, g$exposure, g$sample, sep = "|")
            for (k in unique(g$test_family)) {
                ix <- which(g$test_family == k)
                if (length(ix) != 2L || !setequal(g$endpoint[ix], S8_ENDPOINTS)) stop("Expected exactly two endpoint slots per cohort/design/model/sample")
                good <- g$status[ix] == "ok" & is.finite(g$p_value[ix])
                p <- rep(1, 2)
                p[good] <- g$p_value[ix][good]
                a <- p.adjust(p, "holm")
                g$p_holm_cohort2[ix[good]] <- a[good]
            }
            cross <- paste(g$model, g$exposure, g$sample, sep = "|")
            for (k in unique(cross)) {
                ix <- which(cross == k)
                if (anyDuplicated(paste(g$cohort[ix], g$endpoint[ix]))) stop("Duplicated cross-cohort slots")
                good <- g$status[ix] == "ok" & is.finite(g$p_value[ix])
                p <- rep(1, length(ix))
                p[good] <- g$p_value[ix][good]
                a <- p.adjust(p, "holm")
                g$p_holm_crosscohort_supplementary[ix[good]] <- a[good]
            }
            g$analysis_role <- ifelse(g$sample != "primary" | g$model != "M1", "sensitivity", ifelse(g$exposure == 
                "relative_rank", "primary_relative_position", "supplemental_design"))
            g
        }
        s8_bind <- function(x) {
            if (!length(x)) return(data.frame())
            cols <- unique(unlist(lapply(x, names)))
            do.call(rbind, lapply(x, function(d) {
                for (k in setdiff(cols, names(d))) d[[k]] <- NA
                d[, cols, drop = FALSE]
            }))
        }
        s8_release <- function(r, cells = NULL) {
            for (k in c("pure_factor_levels", "unused_factor_levels", "family_min", "family_max", "birth_rank_max")) r[[k]] <- NULL
            r$family_restriction <- "2<=S<=8; proxy cohorts explicitly labelled S_proxy"
            sensitive <- intersect(c("n", "cases", "noncases", "beta", "OR", "CI_lower", "CI_upper", "p_value", 
                "p_holm_cohort2", "p_holm_crosscohort_supplementary", "statistic", "full_penalized_loglik", "restricted_penalized_loglik"), 
                names(r))
            for (i in seq_len(nrow(r))) {
                counts <- unlist(r[i, intersect(c("n", "cases", "noncases"), names(r)), drop = FALSE], use.names = FALSE)
                if (any(is.finite(counts) & counts <= 5)) {
                  r[i, sensitive] <- NA
                  r$status[i] <- "SUPPRESSED_SMALL_TOTAL"
                }
            }
            if (!is.null(cells) && "rank" %in% names(r) && any(r$exposure == "absolute_rank")) {
                risky <- as.character(cells$group[cells$n <= 5 | cells$cases <= 5 | cells$noncases <= 5])
                reference_risky <- "1" %in% risky
                if (length(risky) == 1L && nrow(cells) > 1L) {
                  other <- which(!as.character(cells$group) %in% risky)
                  if (length(other)) risky <- c(risky, as.character(cells$group[other[which.min(cells$n[other])]]))
                }
                ix <- which(r$exposure == "absolute_rank" & (reference_risky | as.character(r$rank) %in% risky))
                if (length(ix)) {
                  fields <- intersect(c("beta", "OR", "CI_lower", "CI_upper", "p_value"), names(r))
                  r[ix, fields] <- NA
                  r$status[ix] <- "SUPPRESSED_SMALL_CELL"
                }
            }
            r
        }
    }), envir = context)
}

