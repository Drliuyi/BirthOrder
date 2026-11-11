# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

as_num <- function(x) {
    if (is.numeric(x) || is.logical(x)) 
        return(as.numeric(x))
    suppressWarnings(as.numeric(as.character(x)))
}

clean_count <- function(x) {
    z <- as_num(x)
    z[!is.finite(z) | z < 0 | z != floor(z)] <- NA_real_
    z
}

make_birth_order <- function(dat) {
    req <- c("cnt_older", "cnt_brother", "cnt_sister", "adopted")
    miss <- setdiff(req, names(dat))
    if (length(miss)) 
        stop("Missing birth-order columns in all.rds: ", paste(miss, collapse = ", "))
    dat %>% mutate(cnt_older = clean_count(cnt_older), cnt_brother = clean_count(cnt_brother), cnt_sister = clean_count(cnt_sister), 
        adopted_num = if_else(as_num(adopted) %in% c(0, 1), as_num(adopted), NA_real_), sibling_n = cnt_brother + 
            cnt_sister, sibship_size = sibling_n + 1, birth_rank = cnt_older + 1, impossible_rank = !is.na(cnt_older) & 
            !is.na(sibling_n) & cnt_older > sibling_n, eligible_border = !is.na(cnt_older) & !is.na(cnt_brother) & 
            !is.na(cnt_sister) & !impossible_rank & (is.na(adopted_num) | adopted_num != 1), border12 = case_when(!eligible_border | 
            sibship_size != 2 ~ NA_real_, birth_rank == 1 ~ 0, birth_rank == 2 ~ 1, TRUE ~ NA_real_), border1n = case_when(!eligible_border | 
            sibship_size < 2 ~ NA_real_, birth_rank == 1 ~ 0, birth_rank == sibship_size ~ 1, TRUE ~ NA_real_), 
        border_rank = if_else(eligible_border & sibship_size >= 2, birth_rank, NA_real_), border_relative = if_else(eligible_border & 
            sibship_size >= 2, (birth_rank - 1)/(sibship_size - 1), NA_real_))
}

linear_exposures <- c("border12", "border1n", "border_rank", "border_relative")

linear_covariates <- c("age", "sex", "ethnicity", "center", "assessment_year")

linear_units <- c(border12 = "second-born vs first-born; two-child full sibships", border1n = "last-born vs first-born; full sibships >=2; middle excluded", 
    border_rank = "per one-position later birth rank; singleton excluded; middle included", border_relative = "relative position 0 to 1; singleton excluded; middle included")

linear_formulas <- function(exposure) {
    if (length(exposure) != 1L || !exposure %in% linear_exposures) 
        stop("Unknown branch")
    stats::reformulate(c(exposure, linear_covariates, if (exposure != "border12") "factor(S)"), response = "y")
}

ethnicity_top <- function(x) {
    z <- suppressWarnings(as.numeric(as.character(x)))
    labels <- c("White", "Mixed", "Asian or Asian British", "Black or Black British", "Chinese", "Other ethnic group")
    groups <- list(c(1, 1001:1003), c(2, 2001:2004), c(3, 3001:3004), c(4, 4001:4003), 5, 6)
    out <- rep(NA_character_, length(z))
    for (i in seq_along(groups)) out[z %in% groups[[i]]] <- labels[i]
    factor(out, levels = labels)
}

linear_covars <- function(dat) {
    required <- c("eid", "age", "sex", "ethnicity", "center", "date_attend", "sibship_size", linear_exposures)
    absent <- setdiff(required, names(dat))
    if (length(absent)) 
        stop("Missing required baseline adapter fields: ", paste(absent, collapse = ", "))
    d <- as.data.frame(dat)
    if (anyNA(d$eid) || any(!nzchar(as.character(d$eid))) || anyDuplicated(as.character(d$eid))) 
        stop("Duplicate/missing participant IDs")
    d$age <- as_num(d$age)
    d$age[!is.finite(d$age)] <- NA_real_
    s <- as_num(d$sex)
    s[!s %in% c(0, 1)] <- NA_real_
    d$sex <- factor(s, levels = c(0, 1))
    d$ethnicity <- ethnicity_top(d$ethnicity)
    center <- as.character(d$center)
    center[is.na(d$center) | !nzchar(center) | center %in% c("-1", "-3")] <- NA_character_
    d$center <- factor(center)
    date <- if (inherits(d$date_attend, "Date")) 
        d$date_attend
    else suppressWarnings(as.Date(as.character(d$date_attend), format = "%Y-%m-%d"))
    d$assessment_year <- factor(format(date, "%Y"))
    d$S <- as_num(d$sibship_size)
    d$S[!is.finite(d$S) | d$S < 1 | d$S != floor(d$S)] <- NA_real_
    d
}

linear_sample <- function(dat, y, exposure, sex_rule = "both") {
    if (!exposure %in% linear_exposures || !sex_rule %in% c("both", "female", "male")) 
        stop("Invalid branch/sex rule")
    if (is.character(y) && length(y) == 1L) {
        if (!is.na(y) && y %in% names(dat)) 
            y <- dat[[y]]
        else if (!is.na(y) && !y %in% c("0", "1")) 
            stop("Explicit outcome column missing: ", y)
    }
    if (length(y) != nrow(dat)) 
        stop("Outcome vector length does not match cohort")
    if (anyNA(dat$eid) || any(!nzchar(as.character(dat$eid))) || anyDuplicated(as.character(dat$eid))) 
        stop("Duplicate/missing participant IDs")
    if (any(!is.na(y) & !y %in% c(0, 1))) 
        stop("Outcome must retain binary 0/1/NA; invalid states are not controls")
    d <- as.data.frame(dat)
    d$y <- as_num(y)
    f <- linear_formulas(exposure)
    required <- unique(all.vars(f))
    absent <- setdiff(c("eid", required), names(d))
    if (length(absent)) 
        stop("Required model field missing: ", paste(absent, collapse = ", "))
    unknown_n <- sum(is.na(d$y))
    eligible <- !is.na(d[[exposure]]) & is.finite(d[[exposure]])
    if (sex_rule != "both") 
        eligible <- eligible & !is.na(d$sex) & (as.character(d$sex) == if (sex_rule == "female") 
            "0"
        else "1")
    d <- d[eligible, , drop = FALSE]
    keep <- stats::complete.cases(d[, required, drop = FALSE])
    for (k in required) if (is.numeric(d[[k]])) 
        keep <- keep & is.finite(d[[k]])
    d <- droplevels(d[keep, , drop = FALSE])
    terms <- attr(stats::terms(f), "term.labels")
    remove <- character()
    if (sex_rule != "both") 
        remove <- c(remove, "sex")
    if (nrow(d)) 
        for (term in setdiff(terms, exposure)) {
            column <- if (term == "factor(S)") 
                "S"
            else term
            if (length(unique(d[[column]])) <= 1L) 
                remove <- unique(c(remove, term))
        }
    f <- stats::reformulate(setdiff(terms, remove), response = "y")
    list(data = d, formula = f, removed_terms = remove, counts = list(N = nrow(d), cases = sum(d$y == 1), noncases = sum(d$y == 
        0), unknown_outcome_before_cc = unknown_n), note = paste("Removed constant/sex-restricted terms:", paste(remove, 
        collapse = ";")))
}

adjust_branch <- function(results, M, complete = TRUE) {
    if (length(M) != 1L || !is.finite(M) || M < 0 || M != floor(M)) 
        stop("Invalid fixed planned M")
    z <- results
    z$M <- as.integer(M)
    z$q <- z$bonf <- NA_real_
    z$correction_status <- if (complete) 
        "FINAL_FIXED_BRANCH_M"
    else "PARTIAL_NO_FINAL_Q"
    valid <- is.finite(z$p) & z$p >= 0 & z$p <= 1
    if (sum(valid) > M) 
        stop("Successful P count exceeds registered planned denominator")
    if (complete && any(valid) && M > 0) {
        z$q[valid] <- stats::p.adjust(z$p[valid], method = "BH", n = M)
        z$bonf[valid] <- stats::p.adjust(z$p[valid], method = "bonferroni", n = M)
    }
    z
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        suppressPackageStartupMessages(library(dplyr))
        as_num <- function(x) {
            if (is.numeric(x) || is.logical(x)) return(as.numeric(x))
            suppressWarnings(as.numeric(as.character(x)))
        }
        clean_count <- function(x) {
            z <- as_num(x)
            z[!is.finite(z) | z < 0 | z != floor(z)] <- NA_real_
            z
        }
        make_birth_order <- function(dat) {
            req <- c("cnt_older", "cnt_brother", "cnt_sister", "adopted")
            miss <- setdiff(req, names(dat))
            if (length(miss)) stop("Missing birth-order columns in all.rds: ", paste(miss, collapse = ", "))
            dat %>% mutate(cnt_older = clean_count(cnt_older), cnt_brother = clean_count(cnt_brother), cnt_sister = clean_count(cnt_sister), 
                adopted_num = if_else(as_num(adopted) %in% c(0, 1), as_num(adopted), NA_real_), sibling_n = cnt_brother + 
                  cnt_sister, sibship_size = sibling_n + 1, birth_rank = cnt_older + 1, impossible_rank = !is.na(cnt_older) & 
                  !is.na(sibling_n) & cnt_older > sibling_n, eligible_border = !is.na(cnt_older) & !is.na(cnt_brother) & 
                  !is.na(cnt_sister) & !impossible_rank & (is.na(adopted_num) | adopted_num != 1), border12 = case_when(!eligible_border | 
                  sibship_size != 2 ~ NA_real_, birth_rank == 1 ~ 0, birth_rank == 2 ~ 1, TRUE ~ NA_real_), border1n = case_when(!eligible_border | 
                  sibship_size < 2 ~ NA_real_, birth_rank == 1 ~ 0, birth_rank == sibship_size ~ 1, TRUE ~ NA_real_), 
                border_rank = if_else(eligible_border & sibship_size >= 2, birth_rank, NA_real_), border_relative = if_else(eligible_border & 
                  sibship_size >= 2, (birth_rank - 1)/(sibship_size - 1), NA_real_))
        }
        linear_exposures <- c("border12", "border1n", "border_rank", "border_relative")
        linear_covariates <- c("age", "sex", "ethnicity", "center", "assessment_year")
        linear_units <- c(border12 = "second-born vs first-born; two-child full sibships", border1n = "last-born vs first-born; full sibships >=2; middle excluded", 
            border_rank = "per one-position later birth rank; singleton excluded; middle included", border_relative = "relative position 0 to 1; singleton excluded; middle included")
        linear_formulas <- function(exposure) {
            if (length(exposure) != 1L || !exposure %in% linear_exposures) stop("Unknown branch")
            stats::reformulate(c(exposure, linear_covariates, if (exposure != "border12") "factor(S)"), response = "y")
        }
        ethnicity_top <- function(x) {
            z <- suppressWarnings(as.numeric(as.character(x)))
            labels <- c("White", "Mixed", "Asian or Asian British", "Black or Black British", "Chinese", "Other ethnic group")
            groups <- list(c(1, 1001:1003), c(2, 2001:2004), c(3, 3001:3004), c(4, 4001:4003), 5, 6)
            out <- rep(NA_character_, length(z))
            for (i in seq_along(groups)) out[z %in% groups[[i]]] <- labels[i]
            factor(out, levels = labels)
        }
        linear_covars <- function(dat) {
            required <- c("eid", "age", "sex", "ethnicity", "center", "date_attend", "sibship_size", linear_exposures)
            absent <- setdiff(required, names(dat))
            if (length(absent)) stop("Missing required baseline adapter fields: ", paste(absent, collapse = ", "))
            d <- as.data.frame(dat)
            if (anyNA(d$eid) || any(!nzchar(as.character(d$eid))) || anyDuplicated(as.character(d$eid))) stop("Duplicate/missing participant IDs")
            d$age <- as_num(d$age)
            d$age[!is.finite(d$age)] <- NA_real_
            s <- as_num(d$sex)
            s[!s %in% c(0, 1)] <- NA_real_
            d$sex <- factor(s, levels = c(0, 1))
            d$ethnicity <- ethnicity_top(d$ethnicity)
            center <- as.character(d$center)
            center[is.na(d$center) | !nzchar(center) | center %in% c("-1", "-3")] <- NA_character_
            d$center <- factor(center)
            date <- if (inherits(d$date_attend, "Date")) d$date_attend else suppressWarnings(as.Date(as.character(d$date_attend), 
                format = "%Y-%m-%d"))
            d$assessment_year <- factor(format(date, "%Y"))
            d$S <- as_num(d$sibship_size)
            d$S[!is.finite(d$S) | d$S < 1 | d$S != floor(d$S)] <- NA_real_
            d
        }
        linear_sample <- function(dat, y, exposure, sex_rule = "both") {
            if (!exposure %in% linear_exposures || !sex_rule %in% c("both", "female", "male")) stop("Invalid branch/sex rule")
            if (is.character(y) && length(y) == 1L) {
                if (!is.na(y) && y %in% names(dat)) y <- dat[[y]] else if (!is.na(y) && !y %in% c("0", "1")) stop("Explicit outcome column missing: ", 
                  y)
            }
            if (length(y) != nrow(dat)) stop("Outcome vector length does not match cohort")
            if (anyNA(dat$eid) || any(!nzchar(as.character(dat$eid))) || anyDuplicated(as.character(dat$eid))) stop("Duplicate/missing participant IDs")
            if (any(!is.na(y) & !y %in% c(0, 1))) stop("Outcome must retain binary 0/1/NA; invalid states are not controls")
            d <- as.data.frame(dat)
            d$y <- as_num(y)
            f <- linear_formulas(exposure)
            required <- unique(all.vars(f))
            absent <- setdiff(c("eid", required), names(d))
            if (length(absent)) stop("Required model field missing: ", paste(absent, collapse = ", "))
            unknown_n <- sum(is.na(d$y))
            eligible <- !is.na(d[[exposure]]) & is.finite(d[[exposure]])
            if (sex_rule != "both") eligible <- eligible & !is.na(d$sex) & (as.character(d$sex) == if (sex_rule == 
                "female") "0" else "1")
            d <- d[eligible, , drop = FALSE]
            keep <- stats::complete.cases(d[, required, drop = FALSE])
            for (k in required) if (is.numeric(d[[k]])) keep <- keep & is.finite(d[[k]])
            d <- droplevels(d[keep, , drop = FALSE])
            terms <- attr(stats::terms(f), "term.labels")
            remove <- character()
            if (sex_rule != "both") remove <- c(remove, "sex")
            if (nrow(d)) for (term in setdiff(terms, exposure)) {
                column <- if (term == "factor(S)") "S" else term
                if (length(unique(d[[column]])) <= 1L) remove <- unique(c(remove, term))
            }
            f <- stats::reformulate(setdiff(terms, remove), response = "y")
            list(data = d, formula = f, removed_terms = remove, counts = list(N = nrow(d), cases = sum(d$y == 1), 
                noncases = sum(d$y == 0), unknown_outcome_before_cc = unknown_n), note = paste("Removed constant/sex-restricted terms:", 
                paste(remove, collapse = ";")))
        }
        adjust_branch <- function(results, M, complete = TRUE) {
            if (length(M) != 1L || !is.finite(M) || M < 0 || M != floor(M)) stop("Invalid fixed planned M")
            z <- results
            z$M <- as.integer(M)
            z$q <- z$bonf <- NA_real_
            z$correction_status <- if (complete) "FINAL_FIXED_BRANCH_M" else "PARTIAL_NO_FINAL_Q"
            valid <- is.finite(z$p) & z$p >= 0 & z$p <= 1
            if (sum(valid) > M) stop("Successful P count exceeds registered planned denominator")
            if (complete && any(valid) && M > 0) {
                z$q[valid] <- stats::p.adjust(z$p[valid], method = "BH", n = M)
                z$bonf[valid] <- stats::p.adjust(z$p[valid], method = "bonferroni", n = M)
            }
            z
        }
    }), envir = context)
}

