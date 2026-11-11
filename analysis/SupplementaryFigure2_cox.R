# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ip_cox_assert <- function(ok, message) {
    if (!isTRUE(ok)) 
        stop(message, call. = FALSE)
    invisible(TRUE)
}

ip_cox_capture <- function(expr) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = identity)
    list(value = value, warnings = unique(warnings))
}

ip_cox_dates <- function(x) {
    if (inherits(x, "Date")) 
        return(x)
    if (inherits(x, "POSIXt")) 
        return(as.Date(x, tz = "UTC"))
    z <- trimws(as.character(x))
    out <- suppressWarnings(as.Date(z, format = "%Y-%m-%d"))
    valid <- !is.na(z) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z)
    out[!valid] <- as.Date(NA)
    out
}

ip_cox_bad_date <- function(raw, parsed) {
    !is.na(raw) & nzchar(trimws(as.character(raw))) & is.na(parsed)
}

ip_cox_calendar2 <- function(date) {
    if (!length(date)) 
        return(as.Date(character()))
    target <- paste0(as.integer(format(date, "%Y")) + 2L, format(date, "-%m-%d"))
    out <- suppressWarnings(as.Date(target, format = "%Y-%m-%d"))
    leap_day <- !is.na(date) & format(date, "%m-%d") == "02-29"
    if (any(leap_day)) {
        out[leap_day] <- as.Date(paste0(as.integer(format(date[leap_day], "%Y")) + 2L, "-02-28"))
    }
    out[is.na(date)] <- as.Date(NA)
    out
}

ip_riskset <- function(d, event, baseline_status, variant, admin = "2022-05-31") {
    ip_cox_assert(is.data.frame(d), "RISKSET_REQUIRES_DATA_FRAME")
    ip_cox_assert(length(variant) == 1L && variant %in% c("primary", "lag2", "S_le_8"), "UNKNOWN_INCIDENT_VARIANT")
    required <- c("eid", "date_attend", "S", "date_death_uncapped", "date_lost_registry")
    ip_cox_assert(all(required %in% names(d)), paste("RISKSET_FIELDS_MISSING", paste(setdiff(required, names(d)), 
        collapse = ",")))
    n <- nrow(d)
    ip_cox_assert(length(event) == n && inherits(event, "Date"), "EVENT_MUST_BE_A_LENGTH_N_DATE_VECTOR")
    ip_cox_assert(length(baseline_status) == n && all(is.na(baseline_status) | baseline_status %in% c(0, 1)), "BASELINE_STATUS_MUST_RETAIN_0_1_NA")
    ip_cox_assert(!anyNA(d$eid) && all(nzchar(as.character(d$eid))) && !anyDuplicated(as.character(d$eid)), "DUPLICATE_OR_MISSING_EID")
    admin_date <- ip_cox_dates(admin)
    ip_cox_assert(length(admin_date) == 1L && !is.na(admin_date), "ADMIN_CENSOR_MUST_BE_ONE_VALID_DATE")
    attend <- ip_cox_dates(d$date_attend)
    death <- ip_cox_dates(d$date_death_uncapped)
    loss <- ip_cox_dates(d$date_lost_registry)
    entry <- if (variant == "lag2") 
        ip_cox_calendar2(attend)
    else attend
    censor <- pmin(rep(as.numeric(admin_date), n), ifelse(is.na(death), Inf, as.numeric(death)), ifelse(is.na(loss), 
        Inf, as.numeric(loss)))
    reason <- rep("INCLUDED", n)
    flag <- function(condition, label) {
        condition[is.na(condition)] <- FALSE
        reason[reason == "INCLUDED" & condition] <<- label
    }
    flag(is.na(attend), "MISSING_OR_INVALID_BASELINE_DATE")
    flag(is.na(baseline_status), "UNKNOWN_SOURCE_OR_DISEASE_STATUS")
    flag(baseline_status == 1 | (!is.na(event) & event <= attend), "BASELINE_DISEASE_EXCLUDED")
    flag(ip_cox_bad_date(d$date_death_uncapped, death) | ip_cox_bad_date(d$date_lost_registry, loss), "INVALID_CENSOR_DATE")
    flag((!is.na(death) & death < attend) | (!is.na(loss) & loss < attend) | admin_date < attend, "CENSOR_BEFORE_BASELINE")
    if (variant == "S_le_8") {
        ip_cox_assert(is.numeric(d$S), "S_MUST_BE_NUMERIC")
        flag(!is.finite(d$S) | d$S > 8, "OUTSIDE_S_RESTRICTION")
    }
    if (variant == "lag2") 
        flag(!is.na(event) & event <= entry, "EVENT_ON_OR_BEFORE_LANDMARK")
    flag(!is.finite(censor) | is.na(entry) | censor <= as.numeric(entry), "NO_POSITIVE_FOLLOWUP")
    observed <- !is.na(event) & as.numeric(event) > as.numeric(entry) & as.numeric(event) <= censor
    observed[is.na(observed)] <- FALSE
    exit <- ifelse(observed, as.numeric(event), censor)
    keep <- reason == "INCLUDED"
    dd <- d[keep, , drop = FALSE]
    dd$.entry_time <- (as.numeric(entry[keep]) - as.numeric(attend[keep]))/365.25
    dd$.exit_time <- (exit[keep] - as.numeric(attend[keep]))/365.25
    dd$.event <- as.integer(observed[keep])
    ip_cox_assert(all(is.finite(dd$.entry_time)) && all(is.finite(dd$.exit_time)) && all(dd$.exit_time > dd$.entry_time), 
        "INVALID_TIME_ON_STUDY_RISKSET")
    labels <- c("MISSING_OR_INVALID_BASELINE_DATE", "UNKNOWN_SOURCE_OR_DISEASE_STATUS", "BASELINE_DISEASE_EXCLUDED", 
        "INVALID_CENSOR_DATE", "CENSOR_BEFORE_BASELINE", "OUTSIDE_S_RESTRICTION", "EVENT_ON_OR_BEFORE_LANDMARK", 
        "NO_POSITIVE_FOLLOWUP", "INCLUDED")
    counts <- vapply(labels, function(x) sum(reason == x), integer(1))
    flow <- data.frame(variant = variant, step = c("INPUT", labels), count = c(n, counts), stringsAsFactors = FALSE)
    tally <- function(x) sum(x, na.rm = TRUE)
    metrics <- c(input = n, included = sum(keep), excluded = sum(!keep), unknown_status = tally(is.na(baseline_status)), 
        baseline_status_positive = tally(baseline_status == 1), event_on_or_before_baseline = tally(!is.na(event) & 
            event <= attend), event_after_death = tally(!is.na(event) & !is.na(death) & event > death), event_after_loss = tally(!is.na(event) & 
            !is.na(loss) & event > loss), event_after_admin = tally(!is.na(event) & event > admin_date), event_on_censor_included = tally(keep & 
            observed & as.numeric(event) == censor), event_on_death_included = tally(keep & observed & !is.na(death) & 
            event == death), event_on_loss_included = tally(keep & observed & !is.na(loss) & event == loss), event_on_admin_included = tally(keep & 
            observed & event == admin_date), death_is_earliest_censor = tally(keep & !is.na(death) & as.numeric(death) == 
            censor), loss_is_earliest_censor = tally(keep & !is.na(loss) & as.numeric(loss) == censor), admin_is_earliest_censor = tally(keep & 
            as.numeric(admin_date) == censor), landmark_feb29_rolled_to_feb28 = if (variant == "lag2") tally(!is.na(attend) & 
            format(attend, "%m-%d") == "02-29") else 0, events = sum(dd$.event), non_events = sum(dd$.event == 
            0L), person_years = sum(dd$.exit_time - dd$.entry_time))
    audit <- data.frame(variant = variant, metric = names(metrics), value = unname(metrics), stringsAsFactors = FALSE)
    list(data = dd, flow = flow, audit = audit)
}

ip_cox_blank_result <- function(exposure, trait, variant) {
    data.frame(exposure = exposure, trait = trait, variant = variant, test_id = paste(exposure, trait, variant, 
        sep = "__"), N_riskset = NA_integer_, N = NA_integer_, events = NA_integer_, non_events = NA_integer_, 
        excluded_exposure = NA_integer_, excluded_complete_case = NA_integer_, person_years = NA_real_, min_events = 20L, 
        min_non_events = 20L, fit_N = NA_integer_, fit_events = NA_integer_, sample_count_verified = FALSE, event_count_verified = FALSE, 
        raw_beta = NA_real_, raw_se = NA_real_, raw_HR = NA_real_, raw_CI_low = NA_real_, raw_CI_high = NA_real_, 
        raw_p = NA_real_, beta = NA_real_, se = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, 
        p = NA_real_, p_raw = NA_real_, rank = NA_integer_, columns = NA_integer_, iterations = NA_integer_, ph_global_p = NA_real_, 
        ph_exposure_p = NA_real_, ph_min_p = NA_real_, ph_alpha = 0.05, ph_status = "NOT_RUN", status = "NOT_ESTIMABLE", 
        error = "", warning = "", fit_warning = "", ph_warning = "", fit_reliable = FALSE, release_allowed = FALSE, 
        inference = "PENDING", formula = "", requested_formula = "", removed_terms = "", sample_hash = "", sample_hash_algorithm = "sha256", 
        time_scale = "time_on_study_left_truncated_for_landmark", age_handling = "reported_baseline_age_linear_covariate", 
        ties = "efron", stringsAsFactors = FALSE)
}

ip_cox_ph_rows <- function(exposure, trait, variant, status = "NOT_RUN", error = "") {
    data.frame(exposure = exposure, trait = trait, variant = variant, term = c(exposure, "GLOBAL"), chisq = NA_real_, 
        df = NA_real_, p = NA_real_, N = NA_integer_, events = NA_integer_, status = status, error = error, stringsAsFactors = FALSE)
}

ip_cox_coefficient_rows <- function(exposure, trait, variant, terms = exposure) {
    data.frame(exposure = exposure, trait = trait, variant = variant, term = terms, is_target = terms == exposure, 
        raw_beta = NA_real_, raw_se = NA_real_, raw_variance = NA_real_, raw_z = NA_real_, raw_HR = NA_real_, raw_CI_low = NA_real_, 
        raw_CI_high = NA_real_, raw_p = NA_real_, finite_beta = FALSE, positive_finite_variance = FALSE, status = "NOT_FITTED", 
        error = "", warning = "", stringsAsFactors = FALSE)
}

ip_fit <- function(d, exposure, trait, variant, reuse = NULL) {
    r <- ip_cox_blank_result(exposure, trait, variant)
    ph <- ip_cox_ph_rows(exposure, trait, variant)
    coefficients <- ip_cox_coefficient_rows(exposure, trait, variant)
    ph_object <- NULL
    influence <- data.frame()
    r$estimand <- "cause_specific_HR_first_hospital_record"
    r$numeric_status <- "NOT_FITTED"
    r$ph_target_flag <- r$ph_global_flag <- r$ph_covariate_flag <- NA
    r$ph_covariate_terms <- ""
    tryCatch({
        ip_cox_assert(requireNamespace("survival", quietly = TRUE), "SURVIVAL_PACKAGE_REQUIRED")
        ip_cox_assert(length(exposure) == 1L && exposure %in% c("border12", "border1n", "border_rank", "border_relative"), 
            "UNKNOWN_EXPOSURE")
        ip_cox_assert(length(variant) == 1L && variant %in% c("primary", "lag2", "S_le_8"), "UNKNOWN_INCIDENT_VARIANT")
        covars <- c("age", "sex", "ethnicity", "center", "assessment_year", if (exposure != "border12") "factor(S)")
        response <- "survival::Surv(.entry_time, .exit_time, .event)"
        f <- stats::reformulate(c(exposure, covars), response = response)
        r$requested_formula <- r$formula <- paste(deparse(f), collapse = " ")
        required <- unique(c("eid", all.vars(f)))
        ip_cox_assert(is.data.frame(d) && all(required %in% names(d)), paste("MODEL_FIELDS_MISSING", paste(setdiff(required, 
            names(d)), collapse = ",")))
        ip_cox_assert(!anyNA(d$eid) && all(nzchar(as.character(d$eid))) && !anyDuplicated(as.character(d$eid)), 
            "DUPLICATE_OR_MISSING_EID")
        ip_cox_assert(is.numeric(d[[exposure]]) && is.numeric(d$age), "EXPOSURE_AND_BASELINE_AGE_MUST_BE_NUMERIC")
        ip_cox_assert(all(vapply(d[c("sex", "ethnicity", "center", "assessment_year")], is.factor, logical(1))), 
            "CATEGORICAL_COVARIATES_MUST_BE_FACTORS")
        ip_cox_assert(all(is.na(d$.event) | d$.event %in% c(0, 1)), "INVALID_EVENT_INDICATOR")
        r$N_riskset <- nrow(d)
        eligible <- !is.na(d[[exposure]]) & is.finite(d[[exposure]])
        r$excluded_exposure <- sum(!eligible)
        dd <- d[eligible, , drop = FALSE]
        complete <- stats::complete.cases(dd[, required, drop = FALSE])
        for (k in required) if (is.numeric(dd[[k]])) 
            complete <- complete & is.finite(dd[[k]])
        r$excluded_complete_case <- sum(!complete)
        dd <- droplevels(dd[complete, , drop = FALSE])
        r$N <- nrow(dd)
        r$events <- sum(dd$.event)
        r$non_events <- r$N - r$events
        r$person_years <- sum(dd$.exit_time - dd$.entry_time)
        if (exists("bi_hash", mode = "function")) {
            r$sample_hash <- bi_hash(sort(as.character(dd$eid), method = "radix"))
        }
        else {
            ip_cox_assert(requireNamespace("digest", quietly = TRUE), "DIGEST_PACKAGE_REQUIRED")
            r$sample_hash <- digest::digest(sort(as.character(dd$eid), method = "radix"), algo = "sha256")
        }
        removed <- character()
        if (nrow(dd)) 
            for (term in covars) {
                column <- if (term == "factor(S)") 
                  "S"
                else term
                if (length(unique(dd[[column]])) <= 1L) 
                  removed <- c(removed, term)
            }
        r$removed_terms <- paste(removed, collapse = ";")
        f <- stats::reformulate(c(exposure, setdiff(covars, removed)), response = response)
        r$formula <- paste(deparse(f), collapse = " ")
        if (r$events < r$min_events || r$non_events < r$min_non_events) {
            r$status <- "LOW_EVENTS"
            stop(sprintf("INSUFFICIENT_EVENT_OR_NON_EVENT_SUPPORT: events=%d (minimum=%d); non_events=%d (minimum=%d)", 
                r$events, r$min_events, r$non_events, r$min_non_events), call. = FALSE)
        }
        ip_cox_assert(all(dd$.exit_time > dd$.entry_time) && all(dd$.entry_time >= 0), "INVALID_FOLLOWUP_INTERVAL")
        ip_cox_assert(length(unique(dd[[exposure]])) > 1L, "CONSTANT_EXPOSURE")
        mm <- stats::model.matrix(f, data = dd)
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
        coefficients <- ip_cox_coefficient_rows(exposure, trait, variant, colnames(mm))
        r$columns <- ncol(mm)
        ip_cox_assert(all(is.finite(mm)), "COX_NONFINITE_DESIGN")
        r$rank <- qr(mm)$rank
        ip_cox_assert(r$rank == r$columns, "COX_RANK_DEFICIENT_DESIGN")
        if (!is.null(reuse)) {
            q <- reuse$result
            ip_cox_assert(q$N == r$N && q$events == r$events && q$sample_hash == r$sample_hash && q$formula == 
                r$formula, "REUSE_SAMPLE_OR_FORMULA_MISMATCH")
            ip_cox_assert(identical(reuse$coefficients$term, colnames(mm)), "REUSE_COEFFICIENT_ORDER_MISMATCH")
            captured <- ip_cox_capture(survival::coxph(f, data = dd, ties = "efron", init = reuse$coefficients$raw_beta, 
                singular.ok = FALSE, x = TRUE, y = TRUE, model = TRUE, na.action = stats::na.fail, control = survival::coxph.control(iter.max = 0)))
            ip_cox_assert(!length(captured$warnings), "FIXED_COEFFICIENT_DIAGNOSTIC_WARNING")
            captured$warnings <- if (nzchar(q$fit_warning)) 
                q$fit_warning
            else character()
            r$action <- "REUSED_ESTIMATES_FIXED_COEFFICIENT_DIAGNOSTICS_NO_OPTIMISATION"
        }
        else {
            captured <- ip_cox_capture(survival::coxph(f, data = dd, ties = "efron", singular.ok = FALSE, x = TRUE, 
                y = TRUE, model = TRUE, na.action = stats::na.fail, control = survival::coxph.control(iter.max = 50)))
            r$action <- "NEW_ORDINARY_COX_FIT"
        }
        r$fit_warning <- paste(captured$warnings, collapse = "; ")
        r$warning <- r$fit_warning
        if (inherits(captured$value, "error")) 
            stop(conditionMessage(captured$value), call. = FALSE)
        fit <- captured$value
        r$fit_N <- if (length(fit$n) == 1L) 
            as.integer(fit$n)
        else NA_integer_
        r$fit_events <- if (length(fit$nevent) == 1L) 
            as.integer(fit$nevent)
        else NA_integer_
        r$sample_count_verified <- isTRUE(length(fit$n) == 1L && is.finite(fit$n) && fit$n == r$N)
        r$event_count_verified <- isTRUE(length(fit$nevent) == 1L && is.finite(fit$nevent) && fit$nevent == r$events)
        b <- stats::coef(fit)
        variance <- diag(stats::vcov(fit))
        se <- suppressWarnings(sqrt(variance))
        coefficients <- ip_cox_coefficient_rows(exposure, trait, variant, names(b))
        coefficients$raw_beta <- unname(b)
        coefficients$raw_se <- unname(se)
        coefficients$raw_variance <- unname(variance)
        coefficients$raw_z <- unname(b/se)
        coefficients$raw_HR <- unname(exp(b))
        coefficients$raw_CI_low <- unname(exp(b - stats::qnorm(0.975) * se))
        coefficients$raw_CI_high <- unname(exp(b + stats::qnorm(0.975) * se))
        coefficients$raw_p <- unname(2 * stats::pnorm(-abs(b/se)))
        coefficients$finite_beta <- is.finite(b)
        coefficients$positive_finite_variance <- is.finite(variance) & variance > 0
        coefficients$status <- "FITTED_RAW_DIAGNOSTIC"
        coefficients$warning <- r$fit_warning
        if (length(fit$iter)) 
            r$iterations <- as.integer(max(fit$iter))
        target <- match(exposure, coefficients$term)
        ip_cox_assert(!is.na(target), "COX_TARGET_COEFFICIENT_MISSING")
        for (k in c("raw_beta", "raw_se", "raw_HR", "raw_CI_low", "raw_CI_high", "raw_p")) r[[k]] <- coefficients[[k]][target]
        if (!is.null(reuse)) {
            ip_cox_assert(isTRUE(all.equal(unname(b), reuse$coefficients$raw_beta, tolerance = 1e-10)), "REUSE_BETA_CHANGED")
            ip_cox_assert(isTRUE(all.equal(unname(se), reuse$coefficients$raw_se, tolerance = 1e-06)), "REUSE_INFORMATION_CHANGED")
            for (k in c("raw_beta", "raw_se", "raw_HR", "raw_CI_low", "raw_CI_high", "raw_p")) r[[k]] <- reuse$result[[k]]
            r$iterations <- reuse$result$iterations
            r$reuse_exact <- TRUE
        }
        else r$reuse_exact <- FALSE
        issues <- character()
        if (!r$sample_count_verified) 
            issues <- c(issues, sprintf("COX_SAMPLE_COUNT_MISMATCH: complete_N=%d; fit_n=%s", r$N, r$fit_N))
        if (!r$event_count_verified) 
            issues <- c(issues, sprintf("COX_EVENT_COUNT_MISMATCH: complete_events=%d; fit_nevent=%s", r$events, 
                r$fit_events))
        if (length(captured$warnings)) 
            issues <- c(issues, "COX_FIT_WARNING_REQUIRES_REVIEW")
        if (!all(coefficients$finite_beta & coefficients$positive_finite_variance)) 
            issues <- c(issues, "COX_NONFINITE_OR_UNIDENTIFIED_COEFFICIENT")
        if (!all(is.finite(c(r$raw_beta, r$raw_se, r$raw_HR, r$raw_CI_low, r$raw_CI_high, r$raw_p))) || r$raw_HR <= 
            0 || r$raw_CI_low <= 0 || r$raw_se <= 0) 
            issues <- c(issues, "COX_TARGET_EFFECT_OR_INTERVAL_UNRELIABLE")
        r$fit_reliable <- !length(issues)
        r$numeric_status <- if (r$fit_reliable) 
            "ESTIMABLE"
        else "FIT_UNRELIABLE"
        if (!r$fit_reliable) {
            r$status <- "FIT_UNRELIABLE"
            r$error <- paste(issues, collapse = "; ")
            r$ph_status <- "NOT_INTERPRETED_UNRELIABLE_FIT"
            ph <- ip_cox_ph_rows(exposure, trait, variant, r$ph_status, "Resolve numerical estimation before interpreting PH diagnostics")
        }
        else {
            checked <- ip_cox_capture(survival::cox.zph(fit, transform = "km", terms = TRUE, global = TRUE))
            r$ph_warning <- paste(checked$warnings, collapse = "; ")
            r$warning <- paste(c(captured$warnings, checked$warnings), collapse = "; ")
            if (inherits(checked$value, "error")) {
                r$ph_status <- "PH_DIAGNOSTIC_FAILED"
                r$error <- conditionMessage(checked$value)
                ph <- ip_cox_ph_rows(exposure, trait, variant, r$ph_status, r$error)
            }
            else {
                pt <- as.data.frame(checked$value$table)
                valid <- all(c("chisq", "df", "p") %in% names(pt)) && exposure %in% rownames(pt) && "GLOBAL" %in% 
                  rownames(pt) && all(is.finite(as.matrix(pt))) && all(pt$p >= 0 & pt$p <= 1)
                if (!valid) {
                  r$ph_status <- "PH_DIAGNOSTIC_FAILED"
                  r$error <- "PH_DIAGNOSTIC_NONFINITE_OR_INCOMPLETE"
                  ph <- ip_cox_ph_rows(exposure, trait, variant, r$ph_status, r$error)
                }
                else {
                  ph_object <- checked$value
                  r$ph_global_p <- pt["GLOBAL", "p"]
                  r$ph_exposure_p <- pt[exposure, "p"]
                  r$ph_min_p <- min(pt$p)
                  nuisance <- setdiff(rownames(pt), c(exposure, "GLOBAL"))
                  r$ph_target_flag <- r$ph_exposure_p < r$ph_alpha
                  r$ph_global_flag <- r$ph_global_p < r$ph_alpha
                  flagged <- nuisance[pt[nuisance, "p"] < r$ph_alpha]
                  r$ph_covariate_flag <- length(flagged) > 0L
                  r$ph_covariate_terms <- paste(flagged, collapse = ";")
                  r$ph_status <- sc_ph_classify(r$ph_target_flag, r$ph_global_flag, r$ph_covariate_flag, length(checked$warnings) > 
                    0L)
                  ph <- data.frame(exposure = exposure, trait = trait, variant = variant, term = rownames(pt), 
                    chisq = pt$chisq, df = pt$df, p = pt$p, N = r$N, events = r$events, status = r$ph_status, error = "", 
                    stringsAsFactors = FALSE)
                }
            }
            r$status <- "ESTIMABLE_PENDING_REVIEW"
            r$inference <- "Cause-specific HR; PH diagnostics are review flags, not automatic failure"
            influence_check <- ip_cox_capture({
                v <- stats::vcov(fit)
                corr <- stats::cov2cor(v)
                db <- as.matrix(stats::residuals(fit, type = "dfbetas"))
                ip_cox_assert(ncol(db) == length(b), "DFBETAS_COLUMN_COUNT_MISMATCH")
                colnames(db) <- names(b)
                zz <- abs(db[, exposure])
                other <- setdiff(colnames(corr), exposure)
                data.frame(exposure = exposure, trait = trait, variant = variant, max_abs_target_dfbetas = max(zz), 
                  q99_abs_target_dfbetas = unname(quantile(zz, 0.99)), max_abs_target_nuisance_correlation = max(abs(corr[exposure, 
                    other])), most_correlated_term = other[which.max(abs(corr[exposure, other]))], sex_coefficient_correlation = if ("sex1" %in% 
                    other) 
                    corr[exposure, "sex1"]
                  else NA_real_, diagnostic = "Influence summaries are descriptive; no deletion threshold")
            })
            if (inherits(influence_check$value, "error")) {
                r$influence_error <- conditionMessage(influence_check$value)
            }
            else influence <- influence_check$value
        }
    }, error = function(e) {
        r$error <<- conditionMessage(e)
        if (r$status != "LOW_EVENTS") 
            r$status <<- "NOT_ESTIMABLE"
        r$fit_reliable <<- FALSE
        r$release_allowed <<- FALSE
    })
    ph$N <- r$N
    ph$events <- r$events
    if (r$ph_status == "NOT_RUN") 
        ph$error <- r$error
    coefficients$status <- if (all(coefficients$status == "NOT_FITTED")) 
        "NOT_FITTED"
    else r$status
    coefficients$error <- r$error
    coefficients$warning <- r$fit_warning
    list(result = r, ph = ph, coefficients = coefficients, ph_object = ph_object, influence = influence)
}

sc_ph_classify <- function(target, global, nuisance, warned = FALSE) {
    if (warned) 
        return("PH_WARNING_REVIEW")
    if (target) 
        return("EXPOSURE_PH_REVIEW")
    if (global) 
        return("GLOBAL_PH_REVIEW")
    if (nuisance) 
        return("COVARIATE_PH_REVIEW")
    "NO_PH_SIGNAL"
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ip_cox_assert <- function(ok, message) {
            if (!isTRUE(ok)) stop(message, call. = FALSE)
            invisible(TRUE)
        }
        ip_cox_capture <- function(expr) {
            warnings <- character()
            value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }), error = identity)
            list(value = value, warnings = unique(warnings))
        }
        ip_cox_dates <- function(x) {
            if (inherits(x, "Date")) return(x)
            if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))
            z <- trimws(as.character(x))
            out <- suppressWarnings(as.Date(z, format = "%Y-%m-%d"))
            valid <- !is.na(z) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z)
            out[!valid] <- as.Date(NA)
            out
        }
        ip_cox_bad_date <- function(raw, parsed) {
            !is.na(raw) & nzchar(trimws(as.character(raw))) & is.na(parsed)
        }
        ip_cox_calendar2 <- function(date) {
            if (!length(date)) return(as.Date(character()))
            target <- paste0(as.integer(format(date, "%Y")) + 2L, format(date, "-%m-%d"))
            out <- suppressWarnings(as.Date(target, format = "%Y-%m-%d"))
            leap_day <- !is.na(date) & format(date, "%m-%d") == "02-29"
            if (any(leap_day)) {
                out[leap_day] <- as.Date(paste0(as.integer(format(date[leap_day], "%Y")) + 2L, "-02-28"))
            }
            out[is.na(date)] <- as.Date(NA)
            out
        }
        ip_riskset <- function(d, event, baseline_status, variant, admin = "2022-05-31") {
            ip_cox_assert(is.data.frame(d), "RISKSET_REQUIRES_DATA_FRAME")
            ip_cox_assert(length(variant) == 1L && variant %in% c("primary", "lag2", "S_le_8"), "UNKNOWN_INCIDENT_VARIANT")
            required <- c("eid", "date_attend", "S", "date_death_uncapped", "date_lost_registry")
            ip_cox_assert(all(required %in% names(d)), paste("RISKSET_FIELDS_MISSING", paste(setdiff(required, 
                names(d)), collapse = ",")))
            n <- nrow(d)
            ip_cox_assert(length(event) == n && inherits(event, "Date"), "EVENT_MUST_BE_A_LENGTH_N_DATE_VECTOR")
            ip_cox_assert(length(baseline_status) == n && all(is.na(baseline_status) | baseline_status %in% c(0, 
                1)), "BASELINE_STATUS_MUST_RETAIN_0_1_NA")
            ip_cox_assert(!anyNA(d$eid) && all(nzchar(as.character(d$eid))) && !anyDuplicated(as.character(d$eid)), 
                "DUPLICATE_OR_MISSING_EID")
            admin_date <- ip_cox_dates(admin)
            ip_cox_assert(length(admin_date) == 1L && !is.na(admin_date), "ADMIN_CENSOR_MUST_BE_ONE_VALID_DATE")
            attend <- ip_cox_dates(d$date_attend)
            death <- ip_cox_dates(d$date_death_uncapped)
            loss <- ip_cox_dates(d$date_lost_registry)
            entry <- if (variant == "lag2") ip_cox_calendar2(attend) else attend
            censor <- pmin(rep(as.numeric(admin_date), n), ifelse(is.na(death), Inf, as.numeric(death)), ifelse(is.na(loss), 
                Inf, as.numeric(loss)))
            reason <- rep("INCLUDED", n)
            flag <- function(condition, label) {
                condition[is.na(condition)] <- FALSE
                reason[reason == "INCLUDED" & condition] <<- label
            }
            flag(is.na(attend), "MISSING_OR_INVALID_BASELINE_DATE")
            flag(is.na(baseline_status), "UNKNOWN_SOURCE_OR_DISEASE_STATUS")
            flag(baseline_status == 1 | (!is.na(event) & event <= attend), "BASELINE_DISEASE_EXCLUDED")
            flag(ip_cox_bad_date(d$date_death_uncapped, death) | ip_cox_bad_date(d$date_lost_registry, loss), "INVALID_CENSOR_DATE")
            flag((!is.na(death) & death < attend) | (!is.na(loss) & loss < attend) | admin_date < attend, "CENSOR_BEFORE_BASELINE")
            if (variant == "S_le_8") {
                ip_cox_assert(is.numeric(d$S), "S_MUST_BE_NUMERIC")
                flag(!is.finite(d$S) | d$S > 8, "OUTSIDE_S_RESTRICTION")
            }
            if (variant == "lag2") flag(!is.na(event) & event <= entry, "EVENT_ON_OR_BEFORE_LANDMARK")
            flag(!is.finite(censor) | is.na(entry) | censor <= as.numeric(entry), "NO_POSITIVE_FOLLOWUP")
            observed <- !is.na(event) & as.numeric(event) > as.numeric(entry) & as.numeric(event) <= censor
            observed[is.na(observed)] <- FALSE
            exit <- ifelse(observed, as.numeric(event), censor)
            keep <- reason == "INCLUDED"
            dd <- d[keep, , drop = FALSE]
            dd$.entry_time <- (as.numeric(entry[keep]) - as.numeric(attend[keep]))/365.25
            dd$.exit_time <- (exit[keep] - as.numeric(attend[keep]))/365.25
            dd$.event <- as.integer(observed[keep])
            ip_cox_assert(all(is.finite(dd$.entry_time)) && all(is.finite(dd$.exit_time)) && all(dd$.exit_time > 
                dd$.entry_time), "INVALID_TIME_ON_STUDY_RISKSET")
            labels <- c("MISSING_OR_INVALID_BASELINE_DATE", "UNKNOWN_SOURCE_OR_DISEASE_STATUS", "BASELINE_DISEASE_EXCLUDED", 
                "INVALID_CENSOR_DATE", "CENSOR_BEFORE_BASELINE", "OUTSIDE_S_RESTRICTION", "EVENT_ON_OR_BEFORE_LANDMARK", 
                "NO_POSITIVE_FOLLOWUP", "INCLUDED")
            counts <- vapply(labels, function(x) sum(reason == x), integer(1))
            flow <- data.frame(variant = variant, step = c("INPUT", labels), count = c(n, counts), stringsAsFactors = FALSE)
            tally <- function(x) sum(x, na.rm = TRUE)
            metrics <- c(input = n, included = sum(keep), excluded = sum(!keep), unknown_status = tally(is.na(baseline_status)), 
                baseline_status_positive = tally(baseline_status == 1), event_on_or_before_baseline = tally(!is.na(event) & 
                  event <= attend), event_after_death = tally(!is.na(event) & !is.na(death) & event > death), event_after_loss = tally(!is.na(event) & 
                  !is.na(loss) & event > loss), event_after_admin = tally(!is.na(event) & event > admin_date), 
                event_on_censor_included = tally(keep & observed & as.numeric(event) == censor), event_on_death_included = tally(keep & 
                  observed & !is.na(death) & event == death), event_on_loss_included = tally(keep & observed & 
                  !is.na(loss) & event == loss), event_on_admin_included = tally(keep & observed & event == admin_date), 
                death_is_earliest_censor = tally(keep & !is.na(death) & as.numeric(death) == censor), loss_is_earliest_censor = tally(keep & 
                  !is.na(loss) & as.numeric(loss) == censor), admin_is_earliest_censor = tally(keep & as.numeric(admin_date) == 
                  censor), landmark_feb29_rolled_to_feb28 = if (variant == "lag2") tally(!is.na(attend) & format(attend, 
                  "%m-%d") == "02-29") else 0, events = sum(dd$.event), non_events = sum(dd$.event == 0L), person_years = sum(dd$.exit_time - 
                  dd$.entry_time))
            audit <- data.frame(variant = variant, metric = names(metrics), value = unname(metrics), stringsAsFactors = FALSE)
            list(data = dd, flow = flow, audit = audit)
        }
        ip_cox_blank_result <- function(exposure, trait, variant) {
            data.frame(exposure = exposure, trait = trait, variant = variant, test_id = paste(exposure, trait, 
                variant, sep = "__"), N_riskset = NA_integer_, N = NA_integer_, events = NA_integer_, non_events = NA_integer_, 
                excluded_exposure = NA_integer_, excluded_complete_case = NA_integer_, person_years = NA_real_, 
                min_events = 20L, min_non_events = 20L, fit_N = NA_integer_, fit_events = NA_integer_, sample_count_verified = FALSE, 
                event_count_verified = FALSE, raw_beta = NA_real_, raw_se = NA_real_, raw_HR = NA_real_, raw_CI_low = NA_real_, 
                raw_CI_high = NA_real_, raw_p = NA_real_, beta = NA_real_, se = NA_real_, HR = NA_real_, CI_low = NA_real_, 
                CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, rank = NA_integer_, columns = NA_integer_, 
                iterations = NA_integer_, ph_global_p = NA_real_, ph_exposure_p = NA_real_, ph_min_p = NA_real_, 
                ph_alpha = 0.05, ph_status = "NOT_RUN", status = "NOT_ESTIMABLE", error = "", warning = "", fit_warning = "", 
                ph_warning = "", fit_reliable = FALSE, release_allowed = FALSE, inference = "PENDING", formula = "", 
                requested_formula = "", removed_terms = "", sample_hash = "", sample_hash_algorithm = "sha256", 
                time_scale = "time_on_study_left_truncated_for_landmark", age_handling = "reported_baseline_age_linear_covariate", 
                ties = "efron", stringsAsFactors = FALSE)
        }
        ip_cox_ph_rows <- function(exposure, trait, variant, status = "NOT_RUN", error = "") {
            data.frame(exposure = exposure, trait = trait, variant = variant, term = c(exposure, "GLOBAL"), chisq = NA_real_, 
                df = NA_real_, p = NA_real_, N = NA_integer_, events = NA_integer_, status = status, error = error, 
                stringsAsFactors = FALSE)
        }
        ip_cox_coefficient_rows <- function(exposure, trait, variant, terms = exposure) {
            data.frame(exposure = exposure, trait = trait, variant = variant, term = terms, is_target = terms == 
                exposure, raw_beta = NA_real_, raw_se = NA_real_, raw_variance = NA_real_, raw_z = NA_real_, raw_HR = NA_real_, 
                raw_CI_low = NA_real_, raw_CI_high = NA_real_, raw_p = NA_real_, finite_beta = FALSE, positive_finite_variance = FALSE, 
                status = "NOT_FITTED", error = "", warning = "", stringsAsFactors = FALSE)
        }
        ip_fit <- function(d, exposure, trait, variant, reuse = NULL) {
            r <- ip_cox_blank_result(exposure, trait, variant)
            ph <- ip_cox_ph_rows(exposure, trait, variant)
            coefficients <- ip_cox_coefficient_rows(exposure, trait, variant)
            ph_object <- NULL
            influence <- data.frame()
            r$estimand <- "cause_specific_HR_first_hospital_record"
            r$numeric_status <- "NOT_FITTED"
            r$ph_target_flag <- r$ph_global_flag <- r$ph_covariate_flag <- NA
            r$ph_covariate_terms <- ""
            tryCatch({
                ip_cox_assert(requireNamespace("survival", quietly = TRUE), "SURVIVAL_PACKAGE_REQUIRED")
                ip_cox_assert(length(exposure) == 1L && exposure %in% c("border12", "border1n", "border_rank", 
                  "border_relative"), "UNKNOWN_EXPOSURE")
                ip_cox_assert(length(variant) == 1L && variant %in% c("primary", "lag2", "S_le_8"), "UNKNOWN_INCIDENT_VARIANT")
                covars <- c("age", "sex", "ethnicity", "center", "assessment_year", if (exposure != "border12") "factor(S)")
                response <- "survival::Surv(.entry_time, .exit_time, .event)"
                f <- stats::reformulate(c(exposure, covars), response = response)
                r$requested_formula <- r$formula <- paste(deparse(f), collapse = " ")
                required <- unique(c("eid", all.vars(f)))
                ip_cox_assert(is.data.frame(d) && all(required %in% names(d)), paste("MODEL_FIELDS_MISSING", paste(setdiff(required, 
                  names(d)), collapse = ",")))
                ip_cox_assert(!anyNA(d$eid) && all(nzchar(as.character(d$eid))) && !anyDuplicated(as.character(d$eid)), 
                  "DUPLICATE_OR_MISSING_EID")
                ip_cox_assert(is.numeric(d[[exposure]]) && is.numeric(d$age), "EXPOSURE_AND_BASELINE_AGE_MUST_BE_NUMERIC")
                ip_cox_assert(all(vapply(d[c("sex", "ethnicity", "center", "assessment_year")], is.factor, logical(1))), 
                  "CATEGORICAL_COVARIATES_MUST_BE_FACTORS")
                ip_cox_assert(all(is.na(d$.event) | d$.event %in% c(0, 1)), "INVALID_EVENT_INDICATOR")
                r$N_riskset <- nrow(d)
                eligible <- !is.na(d[[exposure]]) & is.finite(d[[exposure]])
                r$excluded_exposure <- sum(!eligible)
                dd <- d[eligible, , drop = FALSE]
                complete <- stats::complete.cases(dd[, required, drop = FALSE])
                for (k in required) if (is.numeric(dd[[k]])) complete <- complete & is.finite(dd[[k]])
                r$excluded_complete_case <- sum(!complete)
                dd <- droplevels(dd[complete, , drop = FALSE])
                r$N <- nrow(dd)
                r$events <- sum(dd$.event)
                r$non_events <- r$N - r$events
                r$person_years <- sum(dd$.exit_time - dd$.entry_time)
                if (exists("bi_hash", mode = "function")) {
                  r$sample_hash <- bi_hash(sort(as.character(dd$eid), method = "radix"))
                } else {
                  ip_cox_assert(requireNamespace("digest", quietly = TRUE), "DIGEST_PACKAGE_REQUIRED")
                  r$sample_hash <- digest::digest(sort(as.character(dd$eid), method = "radix"), algo = "sha256")
                }
                removed <- character()
                if (nrow(dd)) for (term in covars) {
                  column <- if (term == "factor(S)") "S" else term
                  if (length(unique(dd[[column]])) <= 1L) removed <- c(removed, term)
                }
                r$removed_terms <- paste(removed, collapse = ";")
                f <- stats::reformulate(c(exposure, setdiff(covars, removed)), response = response)
                r$formula <- paste(deparse(f), collapse = " ")
                if (r$events < r$min_events || r$non_events < r$min_non_events) {
                  r$status <- "LOW_EVENTS"
                  stop(sprintf("INSUFFICIENT_EVENT_OR_NON_EVENT_SUPPORT: events=%d (minimum=%d); non_events=%d (minimum=%d)", 
                    r$events, r$min_events, r$non_events, r$min_non_events), call. = FALSE)
                }
                ip_cox_assert(all(dd$.exit_time > dd$.entry_time) && all(dd$.entry_time >= 0), "INVALID_FOLLOWUP_INTERVAL")
                ip_cox_assert(length(unique(dd[[exposure]])) > 1L, "CONSTANT_EXPOSURE")
                mm <- stats::model.matrix(f, data = dd)
                mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
                coefficients <- ip_cox_coefficient_rows(exposure, trait, variant, colnames(mm))
                r$columns <- ncol(mm)
                ip_cox_assert(all(is.finite(mm)), "COX_NONFINITE_DESIGN")
                r$rank <- qr(mm)$rank
                ip_cox_assert(r$rank == r$columns, "COX_RANK_DEFICIENT_DESIGN")
                if (!is.null(reuse)) {
                  q <- reuse$result
                  ip_cox_assert(q$N == r$N && q$events == r$events && q$sample_hash == r$sample_hash && q$formula == 
                    r$formula, "REUSE_SAMPLE_OR_FORMULA_MISMATCH")
                  ip_cox_assert(identical(reuse$coefficients$term, colnames(mm)), "REUSE_COEFFICIENT_ORDER_MISMATCH")
                  captured <- ip_cox_capture(survival::coxph(f, data = dd, ties = "efron", init = reuse$coefficients$raw_beta, 
                    singular.ok = FALSE, x = TRUE, y = TRUE, model = TRUE, na.action = stats::na.fail, control = survival::coxph.control(iter.max = 0)))
                  ip_cox_assert(!length(captured$warnings), "FIXED_COEFFICIENT_DIAGNOSTIC_WARNING")
                  captured$warnings <- if (nzchar(q$fit_warning)) q$fit_warning else character()
                  r$action <- "REUSED_ESTIMATES_FIXED_COEFFICIENT_DIAGNOSTICS_NO_OPTIMISATION"
                } else {
                  captured <- ip_cox_capture(survival::coxph(f, data = dd, ties = "efron", singular.ok = FALSE, 
                    x = TRUE, y = TRUE, model = TRUE, na.action = stats::na.fail, control = survival::coxph.control(iter.max = 50)))
                  r$action <- "NEW_ORDINARY_COX_FIT"
                }
                r$fit_warning <- paste(captured$warnings, collapse = "; ")
                r$warning <- r$fit_warning
                if (inherits(captured$value, "error")) stop(conditionMessage(captured$value), call. = FALSE)
                fit <- captured$value
                r$fit_N <- if (length(fit$n) == 1L) as.integer(fit$n) else NA_integer_
                r$fit_events <- if (length(fit$nevent) == 1L) as.integer(fit$nevent) else NA_integer_
                r$sample_count_verified <- isTRUE(length(fit$n) == 1L && is.finite(fit$n) && fit$n == r$N)
                r$event_count_verified <- isTRUE(length(fit$nevent) == 1L && is.finite(fit$nevent) && fit$nevent == 
                  r$events)
                b <- stats::coef(fit)
                variance <- diag(stats::vcov(fit))
                se <- suppressWarnings(sqrt(variance))
                coefficients <- ip_cox_coefficient_rows(exposure, trait, variant, names(b))
                coefficients$raw_beta <- unname(b)
                coefficients$raw_se <- unname(se)
                coefficients$raw_variance <- unname(variance)
                coefficients$raw_z <- unname(b/se)
                coefficients$raw_HR <- unname(exp(b))
                coefficients$raw_CI_low <- unname(exp(b - stats::qnorm(0.975) * se))
                coefficients$raw_CI_high <- unname(exp(b + stats::qnorm(0.975) * se))
                coefficients$raw_p <- unname(2 * stats::pnorm(-abs(b/se)))
                coefficients$finite_beta <- is.finite(b)
                coefficients$positive_finite_variance <- is.finite(variance) & variance > 0
                coefficients$status <- "FITTED_RAW_DIAGNOSTIC"
                coefficients$warning <- r$fit_warning
                if (length(fit$iter)) r$iterations <- as.integer(max(fit$iter))
                target <- match(exposure, coefficients$term)
                ip_cox_assert(!is.na(target), "COX_TARGET_COEFFICIENT_MISSING")
                for (k in c("raw_beta", "raw_se", "raw_HR", "raw_CI_low", "raw_CI_high", "raw_p")) r[[k]] <- coefficients[[k]][target]
                if (!is.null(reuse)) {
                  ip_cox_assert(isTRUE(all.equal(unname(b), reuse$coefficients$raw_beta, tolerance = 1e-10)), "REUSE_BETA_CHANGED")
                  ip_cox_assert(isTRUE(all.equal(unname(se), reuse$coefficients$raw_se, tolerance = 1e-06)), "REUSE_INFORMATION_CHANGED")
                  for (k in c("raw_beta", "raw_se", "raw_HR", "raw_CI_low", "raw_CI_high", "raw_p")) r[[k]] <- reuse$result[[k]]
                  r$iterations <- reuse$result$iterations
                  r$reuse_exact <- TRUE
                } else r$reuse_exact <- FALSE
                issues <- character()
                if (!r$sample_count_verified) issues <- c(issues, sprintf("COX_SAMPLE_COUNT_MISMATCH: complete_N=%d; fit_n=%s", 
                  r$N, r$fit_N))
                if (!r$event_count_verified) issues <- c(issues, sprintf("COX_EVENT_COUNT_MISMATCH: complete_events=%d; fit_nevent=%s", 
                  r$events, r$fit_events))
                if (length(captured$warnings)) issues <- c(issues, "COX_FIT_WARNING_REQUIRES_REVIEW")
                if (!all(coefficients$finite_beta & coefficients$positive_finite_variance)) issues <- c(issues, 
                  "COX_NONFINITE_OR_UNIDENTIFIED_COEFFICIENT")
                if (!all(is.finite(c(r$raw_beta, r$raw_se, r$raw_HR, r$raw_CI_low, r$raw_CI_high, r$raw_p))) || 
                  r$raw_HR <= 0 || r$raw_CI_low <= 0 || r$raw_se <= 0) issues <- c(issues, "COX_TARGET_EFFECT_OR_INTERVAL_UNRELIABLE")
                r$fit_reliable <- !length(issues)
                r$numeric_status <- if (r$fit_reliable) "ESTIMABLE" else "FIT_UNRELIABLE"
                if (!r$fit_reliable) {
                  r$status <- "FIT_UNRELIABLE"
                  r$error <- paste(issues, collapse = "; ")
                  r$ph_status <- "NOT_INTERPRETED_UNRELIABLE_FIT"
                  ph <- ip_cox_ph_rows(exposure, trait, variant, r$ph_status, "Resolve numerical estimation before interpreting PH diagnostics")
                } else {
                  checked <- ip_cox_capture(survival::cox.zph(fit, transform = "km", terms = TRUE, global = TRUE))
                  r$ph_warning <- paste(checked$warnings, collapse = "; ")
                  r$warning <- paste(c(captured$warnings, checked$warnings), collapse = "; ")
                  if (inherits(checked$value, "error")) {
                    r$ph_status <- "PH_DIAGNOSTIC_FAILED"
                    r$error <- conditionMessage(checked$value)
                    ph <- ip_cox_ph_rows(exposure, trait, variant, r$ph_status, r$error)
                  } else {
                    pt <- as.data.frame(checked$value$table)
                    valid <- all(c("chisq", "df", "p") %in% names(pt)) && exposure %in% rownames(pt) && "GLOBAL" %in% 
                      rownames(pt) && all(is.finite(as.matrix(pt))) && all(pt$p >= 0 & pt$p <= 1)
                    if (!valid) {
                      r$ph_status <- "PH_DIAGNOSTIC_FAILED"
                      r$error <- "PH_DIAGNOSTIC_NONFINITE_OR_INCOMPLETE"
                      ph <- ip_cox_ph_rows(exposure, trait, variant, r$ph_status, r$error)
                    } else {
                      ph_object <- checked$value
                      r$ph_global_p <- pt["GLOBAL", "p"]
                      r$ph_exposure_p <- pt[exposure, "p"]
                      r$ph_min_p <- min(pt$p)
                      nuisance <- setdiff(rownames(pt), c(exposure, "GLOBAL"))
                      r$ph_target_flag <- r$ph_exposure_p < r$ph_alpha
                      r$ph_global_flag <- r$ph_global_p < r$ph_alpha
                      flagged <- nuisance[pt[nuisance, "p"] < r$ph_alpha]
                      r$ph_covariate_flag <- length(flagged) > 0L
                      r$ph_covariate_terms <- paste(flagged, collapse = ";")
                      r$ph_status <- sc_ph_classify(r$ph_target_flag, r$ph_global_flag, r$ph_covariate_flag, length(checked$warnings) > 
                        0L)
                      ph <- data.frame(exposure = exposure, trait = trait, variant = variant, term = rownames(pt), 
                        chisq = pt$chisq, df = pt$df, p = pt$p, N = r$N, events = r$events, status = r$ph_status, 
                        error = "", stringsAsFactors = FALSE)
                    }
                  }
                  r$status <- "ESTIMABLE_PENDING_REVIEW"
                  r$inference <- "Cause-specific HR; PH diagnostics are review flags, not automatic failure"
                  influence_check <- ip_cox_capture({
                    v <- stats::vcov(fit)
                    corr <- stats::cov2cor(v)
                    db <- as.matrix(stats::residuals(fit, type = "dfbetas"))
                    ip_cox_assert(ncol(db) == length(b), "DFBETAS_COLUMN_COUNT_MISMATCH")
                    colnames(db) <- names(b)
                    zz <- abs(db[, exposure])
                    other <- setdiff(colnames(corr), exposure)
                    data.frame(exposure = exposure, trait = trait, variant = variant, max_abs_target_dfbetas = max(zz), 
                      q99_abs_target_dfbetas = unname(quantile(zz, 0.99)), max_abs_target_nuisance_correlation = max(abs(corr[exposure, 
                        other])), most_correlated_term = other[which.max(abs(corr[exposure, other]))], sex_coefficient_correlation = if ("sex1" %in% 
                        other) corr[exposure, "sex1"] else NA_real_, diagnostic = "Influence summaries are descriptive; no deletion threshold")
                  })
                  if (inherits(influence_check$value, "error")) {
                    r$influence_error <- conditionMessage(influence_check$value)
                  } else influence <- influence_check$value
                }
            }, error = function(e) {
                r$error <<- conditionMessage(e)
                if (r$status != "LOW_EVENTS") r$status <<- "NOT_ESTIMABLE"
                r$fit_reliable <<- FALSE
                r$release_allowed <<- FALSE
            })
            ph$N <- r$N
            ph$events <- r$events
            if (r$ph_status == "NOT_RUN") ph$error <- r$error
            coefficients$status <- if (all(coefficients$status == "NOT_FITTED")) "NOT_FITTED" else r$status
            coefficients$error <- r$error
            coefficients$warning <- r$fit_warning
            list(result = r, ph = ph, coefficients = coefficients, ph_object = ph_object, influence = influence)
        }
        sc_ph_classify <- function(target, global, nuisance, warned = FALSE) {
            if (warned) return("PH_WARNING_REVIEW")
            if (target) return("EXPOSURE_PH_REVIEW")
            if (global) return("GLOBAL_PH_REVIEW")
            if (nuisance) return("COVARIATE_PH_REVIEW")
            "NO_PH_SIGNAL"
        }
    }), envir = context)
}

