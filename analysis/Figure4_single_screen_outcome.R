# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ff_y_stop <- function(code, ...) stop(paste0(code, if (length(list(...))) ": " else "", paste0(..., collapse = "")), 
    call. = FALSE)

ff_y_assert <- function(ok, code, ...) {
    if (!isTRUE(ok)) 
        ff_y_stop(code, ...)
    invisible(TRUE)
}

ff_y_capture <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = unique(warnings))
}

ff_y_covariates <- c("age", "sex", "ethnicity", "center", "assessment_year", "S")

ff_y_runtime <- function(require = FALSE) {
    packages <- c("survival", "riskRegression", "prodlim", "digest")
    installed <- vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
    versions <- stats::setNames(lapply(packages, function(p) {
        if (installed[[p]]) 
            as.character(utils::packageVersion(p))
        else NA_character_
    }), packages)
    out <- list(installed = as.list(installed), versions = versions, fit_interface = "TWO_SEPARATE_survival::coxph_EFRON_ORIGINAL_SELECTION_WEIGHTS", 
        predict_interface = "survival::basehaz/survfit/predict.coxph + EXACT_REUSED_cp_b_aj_PRODUCT_LIMIT", nonunit_selection_weights = "PRESERVED_NATIVE_SURVIVAL_USER_AUTHORIZED_20261027", 
        robust_cluster_fit = "PRESERVED_NATIVE_SURVIVAL_WHEN_CLUSTER_IS_BOUND", reference_source = "https://github.com/tagteam/riskRegression/blob/master/R/predictCox.R", 
        reference_checked_utc = "2026-10-27")
    out$AJ_backend <- if (exists("ff_aj_status", mode = "function")) 
        ff_aj_status()
    else list(ready = FALSE, state = "FROZEN_R_REFERENCE_ONLY_SYNTHETIC_CONTEXT")
    if (isTRUE(getOption("ff.y.require_compiled_aj", FALSE))) 
        ff_y_assert(isTRUE(out$AJ_backend$ready), "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
    if (require) {
        needed <- c("survival", "digest")
        ff_y_assert(all(installed[needed]), "FF_Y_DEPENDENCY_BLOCKED", paste(needed[!installed[needed]], collapse = ", "))
        ff_y_assert(identical(digest::digest(body(ff_y_aj_original), algo = "sha256"), ff_y_kernel_provenance()$exact_body_sha256), 
            "FF_Y_REUSED_AJ_BODY_HASH_CHANGED")
    }
    out
}

ff_y_validate_endpoint <- function(data, tau_report, time_col = "time", status_col = "status") {
    ff_y_assert(is.data.frame(data) && nrow(data) > 1L, "FF_Y_INPUT_DATA_FRAME_REQUIRED")
    ff_y_assert(length(time_col) == 1L && length(status_col) == 1L && all(c(time_col, status_col) %in% names(data)), 
        "FF_Y_FROZEN_ENDPOINT_COLUMNS_MISSING")
    time <- data[[time_col]]
    status <- data[[status_col]]
    ff_y_assert(is.numeric(time) && length(time) == nrow(data) && all(is.finite(time)) && all(time > 0), "FF_Y_INVALID_FROZEN_FULLFOLLOWUP_TIME")
    ff_y_assert(is.numeric(status) && length(status) == nrow(data) && !anyNA(status) && all(status %in% 0:2), "FF_Y_INVALID_FROZEN_COMPETING_STATUS")
    ff_y_assert(is.numeric(tau_report) && length(tau_report) == 1L && is.finite(tau_report) && tau_report > 0, 
        "FF_Y_TAU_REPORT_NOT_FROZEN_POSITIVE_TIME")
    ff_y_assert(any(time >= tau_report), "FF_Y_TAU_CENSORING_SUPPORT_ABSENT")
    ff_y_assert(any(status == 1L), "FF_Y_NO_DISEASE_EVENTS")
    ff_y_assert(any(status == 2L), "FF_Y_NO_COMPETING_DEATH_EVENTS_TWO_COX_UNIDENTIFIED")
    list(time = time, status = status, n = nrow(data), events_full_followup = sum(status == 1L), deaths_full_followup = sum(status == 
        2L), nondeath_censors_full_followup = sum(status == 0L), max_time = max(time), events_after_tau = sum(status == 
        1L & time > tau_report), deaths_after_tau = sum(status == 2L & time > tau_report), observed_at_risk_at_tau = sum(time >= 
        tau_report), tau_report = tau_report, fit_followup_cap = Inf, time_col = time_col, status_col = status_col)
}

ff_y_validate_weights <- function(data, weight_col = NULL, cluster_col = NULL) {
    if (is.null(weight_col)) {
        present <- intersect(c(".w", "selection_weight", "selection_w", "ipw"), names(data))
        ff_y_assert(!length(present), "FF_Y_EXISTING_WEIGHTS_REQUIRE_EXPLICIT_BINDING", paste(present, collapse = ","))
        w <- rep(1, nrow(data))
    }
    else {
        ff_y_assert(is.character(weight_col) && length(weight_col) == 1L && weight_col %in% names(data), "FF_Y_WEIGHT_COLUMN_MISSING")
        w <- data[[weight_col]]
        ff_y_assert(is.numeric(w) && length(w) == nrow(data) && all(is.finite(w)) && all(w > 0), "FF_Y_INVALID_SELECTION_WEIGHT")
    }
    if (!is.null(cluster_col)) {
        ff_y_assert(is.character(cluster_col) && length(cluster_col) == 1L && cluster_col %in% names(data) && !anyNA(data[[cluster_col]]), 
            "FF_Y_CLUSTER_COLUMN_MISSING_OR_INVALID")
    }
    w
}

ff_y_design <- function(data, covariates, mediator_terms = character(), reference = NULL) {
    ff_y_assert(is.data.frame(data), "FF_Y_INPUT_DATA_FRAME_REQUIRED")
    ff_y_assert(is.character(covariates) && !anyNA(covariates) && !anyDuplicated(covariates), "FF_Y_INVALID_COVARIATE_LIST")
    covariates <- covariates[covariates != "B"]
    ff_y_assert(setequal(covariates, ff_y_covariates), "FF_Y_COVARIATES_NOT_FROZEN_LINEAR_AGE_FACTOR_C")
    ff_y_assert(is.character(mediator_terms) && !anyNA(mediator_terms) && !anyDuplicated(mediator_terms) && !any(mediator_terms %in% 
        c("B", ff_y_covariates)), "FF_Y_INVALID_FROZEN_MEDIATOR_BASIS_NAMES")
    needed <- c("B", covariates, mediator_terms)
    ff_y_assert(all(needed %in% names(data)), "FF_Y_FROZEN_PREDICTOR_COLUMNS_MISSING", paste(setdiff(needed, names(data)), 
        collapse = ","))
    d <- as.data.frame(data[, needed, drop = FALSE])
    for (nm in c("B", "age")) ff_y_assert(is.numeric(d[[nm]]) && all(is.finite(d[[nm]])), "FF_Y_CONTINUOUS_LINEAR_PREDICTOR_REQUIRED", 
        nm)
    ff_y_assert(all(d$B >= 0 & d$B <= 1), "FF_Y_B_OUTSIDE_CONTINUOUS_ZERO_ONE")
    if (is.null(reference)) {
        factor_cols <- union(ff_y_categorical_covariates, mediator_terms[vapply(d[mediator_terms], function(x) is.factor(x) || 
            is.character(x), logical(1L))])
        levels <- stats::setNames(lapply(factor_cols, function(nm) {
            x <- d[[nm]]
            ff_y_assert(!anyNA(x), "FF_Y_MISSING_FACTOR_VALUE", nm)
            lev <- if (is.factor(x)) 
                base::levels(x)
            else sort(unique(as.character(x)))
            ff_y_assert(length(lev) >= 2L && !anyDuplicated(lev), "FF_Y_FACTOR_LEVEL_SUPPORT_ABSENT", nm)
            as.character(lev)
        }), factor_cols)
    }
    else {
        ff_y_assert(identical(needed, reference$needed) && identical(mediator_terms, reference$mediator_terms), 
            "FF_Y_PREDICTOR_SPECIFICATION_DRIFT")
        factor_cols <- reference$factor_cols
        levels <- reference$levels
    }
    for (nm in factor_cols) {
        ff_y_assert(!anyNA(d[[nm]]) && all(as.character(d[[nm]]) %in% levels[[nm]]), "FF_Y_NEW_OR_MISSING_FACTOR_LEVEL", 
            nm)
        d[[nm]] <- factor(as.character(d[[nm]]), levels = levels[[nm]], ordered = FALSE)
    }
    numeric_cols <- setdiff(needed, factor_cols)
    for (nm in numeric_cols) ff_y_assert(is.numeric(d[[nm]]) && all(is.finite(d[[nm]])), "FF_Y_NONFINITE_OR_NONNUMERIC_FROZEN_BASIS", 
        nm)
    f <- stats::reformulate(sprintf("`%s`", needed))
    mf <- stats::model.frame(f, data = d, na.action = stats::na.fail, drop.unused.levels = FALSE)
    contrast <- stats::setNames(lapply(levels, function(lev) stats::contr.treatment(lev, base = 1L)), factor_cols)
    X <- stats::model.matrix(f, data = mf, contrasts.arg = contrast)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    ff_y_assert(nrow(X) == nrow(data) && all(is.finite(X)), "FF_Y_DESIGN_ROW_LOSS_OR_NONFINITE")
    if (is.null(reference)) {
        ff_y_assert(ncol(X) > 0L && qr(X, tol = 1e-11)$rank == ncol(X), "FF_Y_DESIGN_RANK_DEFICIENT_NO_CATEGORY_DELETION")
        reference <- list(needed = needed, mediator_terms = mediator_terms, factor_cols = factor_cols, levels = levels, 
            columns = colnames(X), safe_names = sprintf(".ff_x%03d", seq_len(ncol(X))), contrasts = "contr.treatment baseline=first frozen level", 
            age = "continuous_linear_baseline_age", mediator_transform = "NONE_CALLER_SUPPLIED_V5_FROZEN_BASIS")
    }
    else ff_y_assert(identical(colnames(X), reference$columns), "FF_Y_PREDICTION_DESIGN_DRIFT")
    work <- as.data.frame(X)
    names(work) <- reference$safe_names
    list(X = X, work = work, reference = reference)
}

ff_y_ph <- function(fit, scientific_names = NULL, safe_names = NULL) {
    got <- tryCatch(ff_y_capture(survival::cox.zph(fit, transform = "km")), error = identity)
    if (inherits(got, "error")) 
        return(list(status = "PH_DIAGNOSTIC_UNAVAILABLE", reason = conditionMessage(got)))
    tab <- as.data.frame(got$value$table)
    tab$term <- rownames(tab)
    rownames(tab) <- NULL
    if (length(scientific_names) && length(scientific_names) == length(safe_names)) 
        tab$term[match(safe_names, tab$term, nomatch = 0L)] <- scientific_names[match(safe_names, tab$term, nomatch = 0L) > 
            0L]
    list(status = "DESCRIPTIVE_PH_NO_AUTOMATIC_P_VALUE_REJECTION", table = tab, warnings = got$warnings)
}

ff_y_kernel_provenance <- function() list(source = "scripts/birth/border_competing_paths_figures_v1/f/3.distribution.R::cp_b_aj", 
    source_file_sha256 = birth_private_setting("PRIVATE_HASH_1"), exact_body_sha256 = birth_private_setting("PRIVATE_HASH_2"), 
    reuse = "EXACT_PREVIOUSLY_VALIDATED_BODY_COMPILED_NOT_NEW_ESTIMATOR")

ff_y_risk <- function(risk, tau_report) {
    ff_y_assert(is.data.frame(risk), "FF_Y_RISK_DATA_FRAME_REQUIRED")
    if (all(c("time", "status") %in% names(risk))) {
        rr <- data.frame(time = risk$time, status = risk$status)
    }
    else {
        ff_y_assert(all(c("exit", "event", "death_event") %in% names(risk)), "FF_Y_RISK_COLUMNS_MISSING")
        ff_y_assert(!anyNA(risk$event) && !anyNA(risk$death_event) && all(risk$event %in% 0:1) && all(risk$death_event %in% 
            0:1) && all(risk$event + risk$death_event <= 1), "FF_Y_INVALID_COMPETING_EVENT_BINDING")
        rr <- data.frame(time = risk$exit, status = as.integer(risk$event + 2L * risk$death_event))
    }
    if ("entry" %in% names(risk)) 
        ff_y_assert(all(is.finite(risk$entry)) && all(risk$entry == 0), "FF_Y_NONZERO_ENTRY_UNSUPPORTED_NOT_DISCARDED")
    ff_y_validate_endpoint(rr, tau_report)
}

ff_y_matrix <- function(X, reference = NULL) {
    ff_y_assert(is.matrix(X) && is.numeric(X) && nrow(X) >= 1L && all(is.finite(X)) && !is.null(colnames(X)) && 
        !anyNA(colnames(X)) && !anyDuplicated(colnames(X)), "FF_Y_INVALID_FROZEN_DESIGN_MATRIX")
    if ("(Intercept)" %in% colnames(X)) {
        ff_y_assert(all(X[, "(Intercept)"] == 1), "FF_Y_INVALID_INTERCEPT_COLUMN")
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    ff_y_assert(all(c("B", "age") %in% colnames(X)) && all(X[, "B"] >= 0 & X[, "B"] <= 1) && !any(grepl("ns\\(age|bs\\(age|rcs\\(age|rcspline.*age|attained", 
        colnames(X))), "FF_Y_FROZEN_CONTINUOUS_B_LINEAR_AGE_MATRIX_REQUIRED")
    if (is.null(reference)) {
        ff_y_assert(nrow(X) > 1L && qr(X, tol = 1e-11)$rank == ncol(X), "FF_Y_DESIGN_RANK_DEFICIENT_NO_CATEGORY_DELETION")
        reference <- list(columns = colnames(X), safe_names = sprintf(".ff_x%03d", seq_len(ncol(X))), source = "EXACT_CALLER_SUPPLIED_FROZEN_V5_BASIS_NO_REPARAMETERIZATION")
    }
    else ff_y_assert(identical(colnames(X), reference$columns), "FF_Y_PREDICTION_DESIGN_DRIFT")
    work <- as.data.frame(X)
    names(work) <- reference$safe_names
    list(X = X, work = work, reference = reference)
}

ff_y_fit_native_one <- function(work, safe_names, cause, cluster_bound = FALSE) {
    response <- if (cause == 1L) 
        "survival::Surv(.ff_time, .ff_disease)"
    else "survival::Surv(.ff_time, .ff_death)"
    f <- stats::reformulate(safe_names, response, env = environment())
    if (cluster_bound) 
        ff_y_capture(survival::coxph(f, data = work, weights = .ff_weight, cluster = .ff_cluster, robust = TRUE, 
            ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE, 
            control = survival::coxph.control(iter.max = 50L)))
    else ff_y_capture(survival::coxph(f, data = work, weights = .ff_weight, robust = FALSE, ties = "efron", singular.ok = FALSE, 
        na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE, control = survival::coxph.control(iter.max = 50L)))
}

ff_y_native_baseline <- function(fit, safe_names) {
    zero <- as.data.frame(stats::setNames(as.list(rep(0, length(safe_names))), safe_names))
    bh <- ff_y_capture(survival::basehaz(fit, centered = FALSE))
    sf <- ff_y_capture(survival::survfit(fit, newdata = zero, se.fit = FALSE, conf.type = "none", stype = 2, ctype = 2))
    ff_y_assert(!"strata" %in% names(bh$value) && identical(as.numeric(bh$value$time), as.numeric(sf$value$time)), 
        "FF_Y_NATIVE_BASELINE_GRID_OR_STRATA_CHANGED")
    delta <- max(abs(bh$value$hazard - as.numeric(sf$value$cumhaz)), 0)
    ff_y_assert(is.finite(delta) && delta <= 1e-10 && all(is.finite(bh$value$hazard)) && all(diff(c(0, bh$value$hazard)) >= 
        0), "FF_Y_NATIVE_EFRON_BASEHAZ_SURVFIT_VALIDATION_FAILED")
    list(time = bh$value$time, cumulative = bh$value$hazard, max_abs_basehaz_vs_survfit = delta, tolerance = 1e-10, 
        warnings = unique(c(bh$warnings, sf$warnings)))
}

ff_fit_y_matrix <- function(X, risk, w, cluster = NULL, tau_report, ties = "efron", diagnostics = FALSE) {
    ff_y_assert(identical(ties, "efron"), "FF_Y_TIES_DIFFER_FROM_FROZEN_INCIDENT_EFRON")
    runtime <- ff_y_runtime(require = TRUE)
    dd <- ff_y_matrix(X)
    endpoint <- ff_y_risk(risk, tau_report)
    ff_y_assert(nrow(dd$X) == endpoint$n && is.numeric(w) && length(w) == endpoint$n && all(is.finite(w)) && all(w > 
        0), "FF_Y_MATRIX_RISK_OR_ORIGINAL_WEIGHT_BINDING_INVALID")
    if (!is.null(cluster)) 
        ff_y_assert(length(cluster) == endpoint$n && !anyNA(cluster), "FF_Y_CLUSTER_BINDING_INVALID")
    work <- dd$work
    work$.ff_time <- endpoint$time
    work$.ff_disease <- as.integer(endpoint$status == 1L)
    work$.ff_death <- as.integer(endpoint$status == 2L)
    work$.ff_weight <- w
    if (!is.null(cluster)) 
        work$.ff_cluster <- as.character(cluster)
    fits <- baselines <- vector("list", 2L)
    warnings <- character()
    for (cause in c(1L, 2L)) {
        got <- ff_y_fit_native_one(work, dd$reference$safe_names, cause, !is.null(cluster))
        m <- got$value
        warnings <- c(warnings, got$warnings)
        ff_y_assert(!any(grepl("infinite|converg|singular|iteration", got$warnings, ignore.case = TRUE)) && inherits(m, 
            "coxph") && !inherits(m, "coxphms") && identical(m$method, "efron") && all(is.finite(stats::coef(m))) && 
            all(is.finite(stats::vcov(m))), "FF_Y_NONCONVERGED_OR_UNRELIABLE_NATIVE_EFRON_COX", paste(got$warnings, 
            collapse = " | "))
        mf <- stats::model.frame(m)
        ff_y_assert(as.integer(m$n) == endpoint$n && identical(as.numeric(m$y[, "time"]), as.numeric(endpoint$time)) && 
            identical(as.numeric(m$y[, "status"]), as.numeric(endpoint$status == cause)) && identical(as.numeric(stats::model.weights(mf)), 
            as.numeric(w)) && isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), 
            "FF_Y_FROZEN_TIME_EVENT_SAMPLE_WEIGHT_OR_DESIGN_CHANGED")
        if (!is.null(cluster)) 
            ff_y_assert(!is.null(m$naive.var) && identical(as.character(stats::model.extract(mf, "cluster")), as.character(cluster)), 
                "FF_Y_ROBUST_CLUSTER_BINDING_CHANGED")
        fits[[cause]] <- m
        baselines[[cause]] <- ff_y_native_baseline(m, dd$reference$safe_names)
    }
    names(fits) <- names(baselines) <- c("disease", "death")
    grid <- sort(unique(unlist(lapply(baselines, `[[`, "time"))))
    cumulative <- vapply(baselines, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
    hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
        , drop = FALSE])
    ff_y_assert(all(is.finite(hazard)) && all(hazard >= 0), "FF_Y_INVALID_WEIGHTED_EFRON_HAZARD_INCREMENT")
    keep <- rowSums(hazard) > 0
    engine <- list(time = grid[keep], hazard = hazard[keep, , drop = FALSE], coefficients = vapply(fits, stats::coef, 
        numeric(ncol(dd$X))), provenance = ff_y_kernel_provenance(), baseline_public_validation = baselines)
    coefficient_table <- data.frame(term = dd$reference$columns, internal_term = dd$reference$safe_names, disease_beta = engine$coefficients[, 
        1L], death_beta = engine$coefficients[, 2L], disease_se = sqrt(diag(stats::vcov(fits[[1L]]))), death_se = sqrt(diag(stats::vcov(fits[[2L]]))), 
        row.names = NULL)
    out <- list(fit = fits, engine = engine, matrix_reference = dd$reference, coefficients = coefficient_table, 
        endpoint = endpoint, tau_report = tau_report, ties = ties, runtime = runtime, binding = list(weight = w, 
            cluster = if (is.null(cluster)) NULL else as.character(cluster)), warnings = unique(c(warnings, unlist(lapply(baselines, 
            `[[`, "warnings")))), method = "FULLFOLLOWUP_WEIGHTED_SEPARATE_COX_EFRON_NATIVE_HAZARDS_EXACT_REUSED_AJ", 
        weight = list(min = min(w), max = max(w), sum = sum(w), ESS = sum(w)^2/sum(w^2), interpretation = "ORIGINAL_V5_SELECTION_WEIGHTS_NO_RESCALE_NO_OUTCOME_IPCW"), 
        representative_validation_required = TRUE)
    if (isTRUE(diagnostics)) 
        out$PH <- lapply(fits, ff_y_ph, scientific_names = dd$reference$columns, safe_names = dd$reference$safe_names)
    class(out) <- c("ff_y_matrix_fit", "ff_y_fit")
    probe <- unique(as.integer(round(seq(1, nrow(dd$X), length.out = min(16L, nrow(dd$X))))))
    nd <- dd$work[probe, , drop = FALSE]
    lp_native <- vapply(fits, function(m) as.numeric(stats::predict(m, newdata = nd, type = "lp", reference = "zero")), 
        numeric(length(probe)))
    lp_matrix <- dd$X[probe, , drop = FALSE] %*% engine$coefficients
    out$native_validation <- list(lp_max_abs = max(abs(lp_native - lp_matrix)), lp_tolerance = 1e-10, basehaz_survfit_max_abs = max(vapply(baselines, 
        `[[`, numeric(1L), "max_abs_basehaz_vs_survfit")), basehaz_survfit_tolerance = 1e-10, passed = FALSE)
    ff_y_assert(is.finite(out$native_validation$lp_max_abs) && out$native_validation$lp_max_abs <= out$native_validation$lp_tolerance, 
        "FF_Y_NATIVE_LP_MATRIX_BINDING_VALIDATION_FAILED")
    ff_predict_y_matrix(out, dd$X[probe, , drop = FALSE])
    out$native_validation$passed <- TRUE
    out
}

ff_fit_y <- function(data, covariates, mediator_terms, tau_report, diagnostics = FALSE, ties = "efron", time_col = "time", 
    status_col = "status", weight_col = NULL, cluster_col = NULL) {
    endpoint <- ff_y_validate_endpoint(data, tau_report, time_col, status_col)
    weights <- ff_y_validate_weights(data, weight_col, cluster_col)
    design <- ff_y_design(data, covariates, mediator_terms)
    out <- ff_fit_y_matrix(design$X, data.frame(time = endpoint$time, status = endpoint$status), weights, if (is.null(cluster_col)) 
        NULL
    else data[[cluster_col]], tau_report, ties, diagnostics)
    out$reference <- design$reference
    out$covariates <- covariates
    out$mediator_terms <- mediator_terms
    out$weight$column <- weight_col
    out
}

ff_y_probability_check <- function(p) {
    ff_y_assert(is.matrix(p) && identical(colnames(p), c("disease", "death", "survival")) && all(is.finite(p)) && 
        all(p >= 0 & p <= 1) && max(abs(rowSums(p) - 1), 0) < 1e-08, "FF_Y_INVALID_CIF_PROBABILITY_OR_CLOSURE_NO_CLIPPING")
    invisible(TRUE)
}

ff_y_predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
    64)) {
    ff_y_assert(inherits(object, "ff_y_fit"), "FF_Y_FITTED_OBJECT_REQUIRED")
    ff_y_assert(length(tau_report) == 1L && is.finite(tau_report) && tau_report >= 0 && tau_report <= object$endpoint$max_time, 
        "FF_Y_PREDICTION_TIME_OUTSIDE_SUPPORT")
    ff_y_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1L && length(max_prediction_mb) == 
        1L && is.finite(max_prediction_mb) && max_prediction_mb > 0, "FF_Y_INVALID_PREDICTION_MEMORY_BUDGET")
    dd <- ff_y_matrix(X, object$matrix_reference)
    actual_chunk <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
    p <- matrix(NA_real_, nrow(dd$X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
    warnings <- character()
    h <- object$engine$hazard[object$engine$time <= tau_report, , drop = FALSE]
    for (start in seq.int(1L, nrow(dd$X), by = actual_chunk)) {
        ii <- seq.int(start, min(nrow(dd$X), start + actual_chunk - 1L))
        nd <- dd$work[ii, , drop = FALSE]
        lp <- vapply(object$fit, function(m) {
            got <- ff_y_capture(stats::predict(m, newdata = nd, type = "lp", reference = "zero"))
            warnings <<- c(warnings, got$warnings)
            as.numeric(got$value)
        }, numeric(length(ii)))
        lp <- matrix(lp, nrow = length(ii), ncol = 2L)
        rr <- exp(lp)
        ff_y_assert(all(is.finite(rr)) && all(rr > 0), "FF_Y_NONFINITE_NATIVE_RELATIVE_RISK_NO_CLIPPING")
        if (exists("ff_aj_status", mode = "function") && isTRUE(ff_aj_status()$ready)) {
            p[ii, ] <- ff_aj_backend(h, rr)
        }
        else {
            ff_y_assert(!isTRUE(getOption("ff.y.require_compiled_aj", FALSE)), "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
            p[ii, ] <- ff_y_aj(h, cbind(rr, rr), rep(1L, nrow(h)))
        }
        ff_y_probability_check(p[ii, , drop = FALSE])
    }
    ff_y_probability_check(p)
    attr(p, "warnings") <- unique(warnings)
    attr(p, "chunk_size") <- actual_chunk
    attr(p, "product.limit") <- TRUE
    attr(p, "risk_range") <- apply(p, 2L, range)
    p
}

ff_predict_y_matrix <- function(object, X, tau_report = object$tau_report, chunk_size = getOption("ff.y.chunk_size", 
    4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 64)) {
    ff_y_assert(identical(as.numeric(tau_report), as.numeric(object$tau_report)), "FF_Y_REPORT_TIME_DRIFT_NO_TIMEPOINT_SEARCH")
    ff_y_predict_at(object, X, tau_report, chunk_size, max_prediction_mb)
}

ff_predict_y <- function(object, newdata, tau_report = object$tau_report, chunk_size = getOption("ff.y.chunk_size", 
    4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 64)) {
    ff_y_assert(!is.null(object$reference), "FF_Y_RAW_PREDICT_REQUIRES_RAW_FIT_REFERENCE")
    dd <- ff_y_design(newdata, object$covariates, object$mediator_terms, object$reference)
    ff_predict_y_matrix(object, dd$X, tau_report, chunk_size, max_prediction_mb)
}

ff_y_effect_closure <- function(risks, tolerance = 1e-12) {
    ff_y_assert(is.numeric(risks) && all(c("R00", "R10", "R11") %in% names(risks)) && all(is.finite(risks)) && 
        all(risks >= 0 & risks <= 1), "FF_Y_INVALID_INTEGRATED_SCENARIO_RISKS")
    TE <- unname(risks["R11"] - risks["R00"])
    IDE <- unname(risks["R10"] - risks["R00"])
    IIE <- unname(risks["R11"] - risks["R10"])
    difference <- abs(TE - IDE - IIE)
    ff_y_assert(difference <= tolerance, "FF_Y_INTEGRATED_EFFECT_CLOSURE_FAILED")
    list(TE = TE, IDE = IDE, IIE = IIE, PM = IIE/TE, max_abs = difference, tolerance = tolerance, tolerance_basis = "DOUBLE_PRECISION_ADDITION_1E_12_ABSOLUTE_PROBABILITY", 
        passed = TRUE)
}

ff_validate_y_matrix <- function(object, X, risk, w, cluster = NULL, require_rr = TRUE) {
    coefficient_tolerance <- cif_tolerance <- 1e-09
    dd <- ff_y_matrix(X, object$matrix_reference)
    endpoint <- ff_y_risk(risk, object$tau_report)
    ff_y_assert(identical(as.numeric(endpoint$time), as.numeric(object$endpoint$time)) && identical(endpoint$status, 
        object$endpoint$status) && identical(as.numeric(w), as.numeric(object$binding$weight)) && identical(if (is.null(cluster)) 
        NULL
    else as.character(cluster), object$binding$cluster), "FF_Y_VALIDATION_INPUT_BINDING_DRIFT")
    checks <- list()
    add <- function(name, difference, tolerance, reference) {
        checks[[length(checks) + 1L]] <<- data.frame(check = name, max_abs_difference = difference, tolerance = tolerance, 
            passed = is.finite(difference) && difference <= tolerance, reference = reference)
    }
    reference <- ff_fit_y_matrix(dd$X, risk, w, cluster, object$tau_report)
    add("native_weighted_coefficient", max(abs(object$engine$coefficients - reference$engine$coefficients)), coefficient_tolerance, 
        "INDEPENDENT_PUBLIC_SURVIVAL_COXPH_EFRON_SAME_ORIGINAL_WEIGHTS")
    for (j in 1:2) {
        a <- object$engine$baseline_public_validation[[j]]
        b <- reference$engine$baseline_public_validation[[j]]
        ff_y_assert(identical(a$time, b$time), "FF_Y_VALIDATION_NATIVE_GRID_CHANGED")
        add(paste0("native_weighted_cumulative_hazard_cause", j), max(abs(a$cumulative - b$cumulative)), 1e-10, 
            "INDEPENDENT_PUBLIC_SURVIVAL_BASEHAZ_VS_PUBLIC_SURVFIT_EFRON")
    }
    probe <- unique(as.integer(round(seq(1, nrow(dd$X), length.out = min(16L, nrow(dd$X))))))
    grid <- unique(c(0, object$tau_report * c(0.25, 0.5, 0.75, 1)))
    curve <- lapply(grid, function(t) ff_y_predict_at(object, dd$X[probe, , drop = FALSE], t))
    ranges <- do.call(rbind, lapply(seq_along(grid), function(k) data.frame(time = grid[k], state = colnames(curve[[k]]), 
        min_probability = apply(curve[[k]], 2L, min), max_probability = apply(curve[[k]], 2L, max))))
    add("three_state_raw_probability_range", max(0, -ranges$min_probability, ranges$max_probability - 1), 0, "STRICT_ZERO_ONE_NO_CLIPPING_NO_EPSILON")
    add("disease_death_survival_sum", max(vapply(curve, function(p) max(abs(rowSums(p) - 1)), numeric(1L))), 1e-09, 
        "THREE_STATE_PROBABILITY_CLOSURE_NO_CLIPPING")
    monotonic_violation <- 0
    for (k in 2:length(curve)) {
        delta <- curve[[k]] - curve[[k - 1L]]
        monotonic_violation <- max(monotonic_violation, -delta[, "disease"], -delta[, "death"], delta[, "survival"])
    }
    add("three_state_time_monotonicity", monotonic_violation, 1e-12, "PRESPECIFIED_REPORT_TIME_QUARTILES_NOT_EFFECT_TIMEPOINT_SELECTION")
    rr_ready <- requireNamespace("riskRegression", quietly = TRUE) && requireNamespace("prodlim", quietly = TRUE)
    rr_status <- "NOT_EXECUTED_DEPENDENCY_MISSING"
    if (rr_ready) {
        unit <- if (all(w == 1) && is.null(cluster)) 
            object
        else ff_fit_y_matrix(dd$X, risk, rep(1, endpoint$n), NULL, object$tau_report)
        work <- dd$work
        work$.ff_time <- endpoint$time
        work$.ff_status <- endpoint$status
        f <- stats::reformulate(names(dd$work), "prodlim::Hist(.ff_time, .ff_status)")
        got <- ff_y_capture(riskRegression::CSC(f, data = work, cause = 1, surv.type = "hazard", fitter = "coxph", 
            ties = "efron", model = TRUE, singular.ok = FALSE, robust = FALSE, na.action = stats::na.fail, control = survival::coxph.control(iter.max = 50L)))
        ff_y_assert(!length(got$warnings), "FF_Y_UNIT_RR_REFERENCE_FIT_WARNING", paste(got$warnings, collapse = " | "))
        public <- got$value
        for (j in 1:2) add(paste0("unit_coefficient_vs_RR_cause", j), max(abs(stats::coef(unit$fit[[j]]) - stats::coef(public$models[[paste("Cause", 
            j)]]))), coefficient_tolerance, "INDEPENDENT_UNIT_WEIGHT_PUBLIC_RISKREGRESSION_CSC_EFRON")
        for (t in grid[-1L]) {
            p <- ff_y_predict_at(unit, dd$X[probe, , drop = FALSE], t)
            for (j in 1:3) {
                type <- if (j == 3L) 
                  "survival"
                else "absRisk"
                got <- ff_y_capture(stats::predict(public, newdata = dd$work[probe, , drop = FALSE], times = t, 
                  cause = if (j == 2L) 
                    2
                  else 1, type = type, product.limit = TRUE, se = FALSE, band = FALSE, iid = FALSE, confint = FALSE))
                ff_y_assert(!length(got$warnings), "FF_Y_UNIT_RR_REFERENCE_PREDICT_WARNING", paste(got$warnings, 
                  collapse = " | "))
                add(paste0("unit_", colnames(p)[j], "_vs_RR_t", format(t)), max(abs(p[, j] - as.numeric(got$value[[type]]))), 
                  cif_tolerance, "PUBLIC_RISKREGRESSION_PREDICT_CAUSESPECIFICCOX_PRODUCT_LIMIT_TRUE")
            }
        }
        rr_status <- "EXECUTED"
    }
    table <- do.call(rbind, checks)
    passed <- all(table$passed) && (rr_ready || !require_rr)
    out <- list(passed = passed, certified_for_representative_gate = all(table$passed) && rr_ready, native_checks_passed = all(table$passed), 
        checks = table, risk_ranges_by_time = ranges, riskRegression_unit_reference_status = rr_status, runtime = ff_y_runtime(), 
        kernel_provenance = ff_y_kernel_provenance(), tolerance_basis = paste("ABSOLUTE_DOUBLE_PRECISION_NUMERICAL_EQUIVALENCE:", 
            "coefficients/CIF 1e-9; public cumhaz/LP 1e-10; monotonicity/effect closure 1e-12;", "strict raw ranges and jump bounds; no clipping/epsilon; not statistical/P-value gates"), 
        integrated_effect_validation = "CALL_ff_y_effect_closure_ON_ACTUAL_A_INTEGRATED_R00_R10_R11", representative_scope = "PRESPECIFIED_ENGINEERING_COMBINATIONS_POINT_AND_B1_20_NOT_NEW_DRAWS")
    out$compiled_backend_acceptance <- if (exists("ff_aj_status", mode = "function")) 
        ff_aj_status()
    else list(ready = FALSE, state = "FROZEN_R_REFERENCE_ONLY_SYNTHETIC_CONTEXT")
    out$certified_for_representative_gate <- isTRUE(out$certified_for_representative_gate) && isTRUE(out$compiled_backend_acceptance$ready) && 
        isTRUE(out$compiled_backend_acceptance$acceptance$passed)
    if (isTRUE(getOption("ff.y.require_compiled_aj", FALSE))) {
        ff_y_assert(isTRUE(out$compiled_backend_acceptance$ready), "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
    }
    if (!all(table$passed) || (require_rr && !rr_ready)) {
        message <- if (!all(table$passed)) 
            paste0("FF_Y_REPRESENTATIVE_NUMERICAL_VALIDATION_FAILED: ", paste(table$check[!table$passed], collapse = ", "))
        else "FF_Y_REPRESENTATIVE_UNIT_RR_VALIDATION_DEPENDENCY_MISSING"
        stop(structure(list(message = message, call = NULL, validation = out), class = c("ff_y_validation_error", 
            "error", "condition")))
    }
    out
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ff_y_stop <- function(code, ...) stop(paste0(code, if (length(list(...))) ": " else "", paste0(..., collapse = "")), 
            call. = FALSE)
        ff_y_assert <- function(ok, code, ...) {
            if (!isTRUE(ok)) ff_y_stop(code, ...)
            invisible(TRUE)
        }
        ff_y_capture <- function(expr) {
            warnings <- character()
            value <- withCallingHandlers(expr, warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
            list(value = value, warnings = unique(warnings))
        }
        ff_y_covariates <- c("age", "sex", "ethnicity", "center", "assessment_year", "S")
        ff_y_categorical_covariates <- setdiff(ff_y_covariates, "age")
        ff_y_runtime <- function(require = FALSE) {
            packages <- c("survival", "riskRegression", "prodlim", "digest")
            installed <- vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
            versions <- stats::setNames(lapply(packages, function(p) {
                if (installed[[p]]) as.character(utils::packageVersion(p)) else NA_character_
            }), packages)
            out <- list(installed = as.list(installed), versions = versions, fit_interface = "TWO_SEPARATE_survival::coxph_EFRON_ORIGINAL_SELECTION_WEIGHTS", 
                predict_interface = "survival::basehaz/survfit/predict.coxph + EXACT_REUSED_cp_b_aj_PRODUCT_LIMIT", 
                nonunit_selection_weights = "PRESERVED_NATIVE_SURVIVAL_USER_AUTHORIZED_20261027", robust_cluster_fit = "PRESERVED_NATIVE_SURVIVAL_WHEN_CLUSTER_IS_BOUND", 
                reference_source = "https://github.com/tagteam/riskRegression/blob/master/R/predictCox.R", reference_checked_utc = "2026-10-27")
            out$AJ_backend <- if (exists("ff_aj_status", mode = "function")) ff_aj_status() else list(ready = FALSE, 
                state = "FROZEN_R_REFERENCE_ONLY_SYNTHETIC_CONTEXT")
            if (isTRUE(getOption("ff.y.require_compiled_aj", FALSE))) ff_y_assert(isTRUE(out$AJ_backend$ready), 
                "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
            if (require) {
                needed <- c("survival", "digest")
                ff_y_assert(all(installed[needed]), "FF_Y_DEPENDENCY_BLOCKED", paste(needed[!installed[needed]], 
                  collapse = ", "))
                ff_y_assert(identical(digest::digest(body(ff_y_aj_original), algo = "sha256"), ff_y_kernel_provenance()$exact_body_sha256), 
                  "FF_Y_REUSED_AJ_BODY_HASH_CHANGED")
            }
            out
        }
        ff_y_validate_endpoint <- function(data, tau_report, time_col = "time", status_col = "status") {
            ff_y_assert(is.data.frame(data) && nrow(data) > 1L, "FF_Y_INPUT_DATA_FRAME_REQUIRED")
            ff_y_assert(length(time_col) == 1L && length(status_col) == 1L && all(c(time_col, status_col) %in% 
                names(data)), "FF_Y_FROZEN_ENDPOINT_COLUMNS_MISSING")
            time <- data[[time_col]]
            status <- data[[status_col]]
            ff_y_assert(is.numeric(time) && length(time) == nrow(data) && all(is.finite(time)) && all(time > 0), 
                "FF_Y_INVALID_FROZEN_FULLFOLLOWUP_TIME")
            ff_y_assert(is.numeric(status) && length(status) == nrow(data) && !anyNA(status) && all(status %in% 
                0:2), "FF_Y_INVALID_FROZEN_COMPETING_STATUS")
            ff_y_assert(is.numeric(tau_report) && length(tau_report) == 1L && is.finite(tau_report) && tau_report > 
                0, "FF_Y_TAU_REPORT_NOT_FROZEN_POSITIVE_TIME")
            ff_y_assert(any(time >= tau_report), "FF_Y_TAU_CENSORING_SUPPORT_ABSENT")
            ff_y_assert(any(status == 1L), "FF_Y_NO_DISEASE_EVENTS")
            ff_y_assert(any(status == 2L), "FF_Y_NO_COMPETING_DEATH_EVENTS_TWO_COX_UNIDENTIFIED")
            list(time = time, status = status, n = nrow(data), events_full_followup = sum(status == 1L), deaths_full_followup = sum(status == 
                2L), nondeath_censors_full_followup = sum(status == 0L), max_time = max(time), events_after_tau = sum(status == 
                1L & time > tau_report), deaths_after_tau = sum(status == 2L & time > tau_report), observed_at_risk_at_tau = sum(time >= 
                tau_report), tau_report = tau_report, fit_followup_cap = Inf, time_col = time_col, status_col = status_col)
        }
        ff_y_validate_weights <- function(data, weight_col = NULL, cluster_col = NULL) {
            if (is.null(weight_col)) {
                present <- intersect(c(".w", "selection_weight", "selection_w", "ipw"), names(data))
                ff_y_assert(!length(present), "FF_Y_EXISTING_WEIGHTS_REQUIRE_EXPLICIT_BINDING", paste(present, 
                  collapse = ","))
                w <- rep(1, nrow(data))
            } else {
                ff_y_assert(is.character(weight_col) && length(weight_col) == 1L && weight_col %in% names(data), 
                  "FF_Y_WEIGHT_COLUMN_MISSING")
                w <- data[[weight_col]]
                ff_y_assert(is.numeric(w) && length(w) == nrow(data) && all(is.finite(w)) && all(w > 0), "FF_Y_INVALID_SELECTION_WEIGHT")
            }
            if (!is.null(cluster_col)) {
                ff_y_assert(is.character(cluster_col) && length(cluster_col) == 1L && cluster_col %in% names(data) && 
                  !anyNA(data[[cluster_col]]), "FF_Y_CLUSTER_COLUMN_MISSING_OR_INVALID")
            }
            w
        }
        ff_y_design <- function(data, covariates, mediator_terms = character(), reference = NULL) {
            ff_y_assert(is.data.frame(data), "FF_Y_INPUT_DATA_FRAME_REQUIRED")
            ff_y_assert(is.character(covariates) && !anyNA(covariates) && !anyDuplicated(covariates), "FF_Y_INVALID_COVARIATE_LIST")
            covariates <- covariates[covariates != "B"]
            ff_y_assert(setequal(covariates, ff_y_covariates), "FF_Y_COVARIATES_NOT_FROZEN_LINEAR_AGE_FACTOR_C")
            ff_y_assert(is.character(mediator_terms) && !anyNA(mediator_terms) && !anyDuplicated(mediator_terms) && 
                !any(mediator_terms %in% c("B", ff_y_covariates)), "FF_Y_INVALID_FROZEN_MEDIATOR_BASIS_NAMES")
            needed <- c("B", covariates, mediator_terms)
            ff_y_assert(all(needed %in% names(data)), "FF_Y_FROZEN_PREDICTOR_COLUMNS_MISSING", paste(setdiff(needed, 
                names(data)), collapse = ","))
            d <- as.data.frame(data[, needed, drop = FALSE])
            for (nm in c("B", "age")) ff_y_assert(is.numeric(d[[nm]]) && all(is.finite(d[[nm]])), "FF_Y_CONTINUOUS_LINEAR_PREDICTOR_REQUIRED", 
                nm)
            ff_y_assert(all(d$B >= 0 & d$B <= 1), "FF_Y_B_OUTSIDE_CONTINUOUS_ZERO_ONE")
            if (is.null(reference)) {
                factor_cols <- union(ff_y_categorical_covariates, mediator_terms[vapply(d[mediator_terms], function(x) is.factor(x) || 
                  is.character(x), logical(1L))])
                levels <- stats::setNames(lapply(factor_cols, function(nm) {
                  x <- d[[nm]]
                  ff_y_assert(!anyNA(x), "FF_Y_MISSING_FACTOR_VALUE", nm)
                  lev <- if (is.factor(x)) base::levels(x) else sort(unique(as.character(x)))
                  ff_y_assert(length(lev) >= 2L && !anyDuplicated(lev), "FF_Y_FACTOR_LEVEL_SUPPORT_ABSENT", nm)
                  as.character(lev)
                }), factor_cols)
            } else {
                ff_y_assert(identical(needed, reference$needed) && identical(mediator_terms, reference$mediator_terms), 
                  "FF_Y_PREDICTOR_SPECIFICATION_DRIFT")
                factor_cols <- reference$factor_cols
                levels <- reference$levels
            }
            for (nm in factor_cols) {
                ff_y_assert(!anyNA(d[[nm]]) && all(as.character(d[[nm]]) %in% levels[[nm]]), "FF_Y_NEW_OR_MISSING_FACTOR_LEVEL", 
                  nm)
                d[[nm]] <- factor(as.character(d[[nm]]), levels = levels[[nm]], ordered = FALSE)
            }
            numeric_cols <- setdiff(needed, factor_cols)
            for (nm in numeric_cols) ff_y_assert(is.numeric(d[[nm]]) && all(is.finite(d[[nm]])), "FF_Y_NONFINITE_OR_NONNUMERIC_FROZEN_BASIS", 
                nm)
            f <- stats::reformulate(sprintf("`%s`", needed))
            mf <- stats::model.frame(f, data = d, na.action = stats::na.fail, drop.unused.levels = FALSE)
            contrast <- stats::setNames(lapply(levels, function(lev) stats::contr.treatment(lev, base = 1L)), factor_cols)
            X <- stats::model.matrix(f, data = mf, contrasts.arg = contrast)
            X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            ff_y_assert(nrow(X) == nrow(data) && all(is.finite(X)), "FF_Y_DESIGN_ROW_LOSS_OR_NONFINITE")
            if (is.null(reference)) {
                ff_y_assert(ncol(X) > 0L && qr(X, tol = 1e-11)$rank == ncol(X), "FF_Y_DESIGN_RANK_DEFICIENT_NO_CATEGORY_DELETION")
                reference <- list(needed = needed, mediator_terms = mediator_terms, factor_cols = factor_cols, 
                  levels = levels, columns = colnames(X), safe_names = sprintf(".ff_x%03d", seq_len(ncol(X))), 
                  contrasts = "contr.treatment baseline=first frozen level", age = "continuous_linear_baseline_age", 
                  mediator_transform = "NONE_CALLER_SUPPLIED_V5_FROZEN_BASIS")
            } else ff_y_assert(identical(colnames(X), reference$columns), "FF_Y_PREDICTION_DESIGN_DRIFT")
            work <- as.data.frame(X)
            names(work) <- reference$safe_names
            list(X = X, work = work, reference = reference)
        }
        ff_y_ph <- function(fit, scientific_names = NULL, safe_names = NULL) {
            got <- tryCatch(ff_y_capture(survival::cox.zph(fit, transform = "km")), error = identity)
            if (inherits(got, "error")) return(list(status = "PH_DIAGNOSTIC_UNAVAILABLE", reason = conditionMessage(got)))
            tab <- as.data.frame(got$value$table)
            tab$term <- rownames(tab)
            rownames(tab) <- NULL
            if (length(scientific_names) && length(scientific_names) == length(safe_names)) tab$term[match(safe_names, 
                tab$term, nomatch = 0L)] <- scientific_names[match(safe_names, tab$term, nomatch = 0L) > 0L]
            list(status = "DESCRIPTIVE_PH_NO_AUTOMATIC_P_VALUE_REJECTION", table = tab, warnings = got$warnings)
        }
        ff_y_kernel_provenance <- function() list(source = "scripts/birth/border_competing_paths_figures_v1/f/3.distribution.R::cp_b_aj", 
            source_file_sha256 = birth_private_setting("PRIVATE_HASH_3"), exact_body_sha256 = birth_private_setting("PRIVATE_HASH_4"), 
            reuse = "EXACT_PREVIOUSLY_VALIDATED_BODY_COMPILED_NOT_NEW_ESTIMATOR")
        ff_y_aj_original <- local({
            cp_b_assert <- ff_y_assert
            function(hazard, rr, period) {
                n <- nrow(rr)
                s <- rep(1, n)
                f1 <- f2 <- numeric(n)
                for (j in seq_len(nrow(hazard))) {
                  k <- if (period[[j]] == 1L) 1L else 3L
                  h1 <- hazard[j, 1L] * rr[, k]
                  h2 <- hazard[j, 2L] * rr[, k + 1L]
                  cp_b_assert(all(is.finite(h1 + h2)) && all(h1 + h2 <= 1), "B_AJ_JUMP_PROBABILITY_INVALID")
                  f1 <- f1 + s * h1
                  f2 <- f2 + s * h2
                  s <- s * (1 - h1 - h2)
                }
                cp_b_assert(all(f1 >= 0 & f2 >= 0 & s >= 0) && max(abs(f1 + f2 + s - 1)) < 1e-09, "B_PIECEWISE_COMPETING_PROBABILITY_FAILURE")
                cbind(disease = f1, death = f2, survival = s)
            }
        })
        ff_y_aj <- compiler::cmpfun(ff_y_aj_original)
        ff_y_risk <- function(risk, tau_report) {
            ff_y_assert(is.data.frame(risk), "FF_Y_RISK_DATA_FRAME_REQUIRED")
            if (all(c("time", "status") %in% names(risk))) {
                rr <- data.frame(time = risk$time, status = risk$status)
            } else {
                ff_y_assert(all(c("exit", "event", "death_event") %in% names(risk)), "FF_Y_RISK_COLUMNS_MISSING")
                ff_y_assert(!anyNA(risk$event) && !anyNA(risk$death_event) && all(risk$event %in% 0:1) && all(risk$death_event %in% 
                  0:1) && all(risk$event + risk$death_event <= 1), "FF_Y_INVALID_COMPETING_EVENT_BINDING")
                rr <- data.frame(time = risk$exit, status = as.integer(risk$event + 2L * risk$death_event))
            }
            if ("entry" %in% names(risk)) ff_y_assert(all(is.finite(risk$entry)) && all(risk$entry == 0), "FF_Y_NONZERO_ENTRY_UNSUPPORTED_NOT_DISCARDED")
            ff_y_validate_endpoint(rr, tau_report)
        }
        ff_y_matrix <- function(X, reference = NULL) {
            ff_y_assert(is.matrix(X) && is.numeric(X) && nrow(X) >= 1L && all(is.finite(X)) && !is.null(colnames(X)) && 
                !anyNA(colnames(X)) && !anyDuplicated(colnames(X)), "FF_Y_INVALID_FROZEN_DESIGN_MATRIX")
            if ("(Intercept)" %in% colnames(X)) {
                ff_y_assert(all(X[, "(Intercept)"] == 1), "FF_Y_INVALID_INTERCEPT_COLUMN")
                X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            }
            ff_y_assert(all(c("B", "age") %in% colnames(X)) && all(X[, "B"] >= 0 & X[, "B"] <= 1) && !any(grepl("ns\\(age|bs\\(age|rcs\\(age|rcspline.*age|attained", 
                colnames(X))), "FF_Y_FROZEN_CONTINUOUS_B_LINEAR_AGE_MATRIX_REQUIRED")
            if (is.null(reference)) {
                ff_y_assert(nrow(X) > 1L && qr(X, tol = 1e-11)$rank == ncol(X), "FF_Y_DESIGN_RANK_DEFICIENT_NO_CATEGORY_DELETION")
                reference <- list(columns = colnames(X), safe_names = sprintf(".ff_x%03d", seq_len(ncol(X))), source = "EXACT_CALLER_SUPPLIED_FROZEN_V5_BASIS_NO_REPARAMETERIZATION")
            } else ff_y_assert(identical(colnames(X), reference$columns), "FF_Y_PREDICTION_DESIGN_DRIFT")
            work <- as.data.frame(X)
            names(work) <- reference$safe_names
            list(X = X, work = work, reference = reference)
        }
        ff_y_fit_native_one <- function(work, safe_names, cause, cluster_bound = FALSE) {
            response <- if (cause == 1L) "survival::Surv(.ff_time, .ff_disease)" else "survival::Surv(.ff_time, .ff_death)"
            f <- stats::reformulate(safe_names, response, env = environment())
            if (cluster_bound) ff_y_capture(survival::coxph(f, data = work, weights = .ff_weight, cluster = .ff_cluster, 
                robust = TRUE, ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, x = TRUE, 
                y = TRUE, control = survival::coxph.control(iter.max = 50L))) else ff_y_capture(survival::coxph(f, 
                data = work, weights = .ff_weight, robust = FALSE, ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, 
                model = TRUE, x = TRUE, y = TRUE, control = survival::coxph.control(iter.max = 50L)))
        }
        ff_y_native_baseline <- function(fit, safe_names) {
            zero <- as.data.frame(stats::setNames(as.list(rep(0, length(safe_names))), safe_names))
            bh <- ff_y_capture(survival::basehaz(fit, centered = FALSE))
            sf <- ff_y_capture(survival::survfit(fit, newdata = zero, se.fit = FALSE, conf.type = "none", stype = 2, 
                ctype = 2))
            ff_y_assert(!"strata" %in% names(bh$value) && identical(as.numeric(bh$value$time), as.numeric(sf$value$time)), 
                "FF_Y_NATIVE_BASELINE_GRID_OR_STRATA_CHANGED")
            delta <- max(abs(bh$value$hazard - as.numeric(sf$value$cumhaz)), 0)
            ff_y_assert(is.finite(delta) && delta <= 1e-10 && all(is.finite(bh$value$hazard)) && all(diff(c(0, 
                bh$value$hazard)) >= 0), "FF_Y_NATIVE_EFRON_BASEHAZ_SURVFIT_VALIDATION_FAILED")
            list(time = bh$value$time, cumulative = bh$value$hazard, max_abs_basehaz_vs_survfit = delta, tolerance = 1e-10, 
                warnings = unique(c(bh$warnings, sf$warnings)))
        }
        ff_fit_y_matrix <- function(X, risk, w, cluster = NULL, tau_report, ties = "efron", diagnostics = FALSE) {
            ff_y_assert(identical(ties, "efron"), "FF_Y_TIES_DIFFER_FROM_FROZEN_INCIDENT_EFRON")
            runtime <- ff_y_runtime(require = TRUE)
            dd <- ff_y_matrix(X)
            endpoint <- ff_y_risk(risk, tau_report)
            ff_y_assert(nrow(dd$X) == endpoint$n && is.numeric(w) && length(w) == endpoint$n && all(is.finite(w)) && 
                all(w > 0), "FF_Y_MATRIX_RISK_OR_ORIGINAL_WEIGHT_BINDING_INVALID")
            if (!is.null(cluster)) ff_y_assert(length(cluster) == endpoint$n && !anyNA(cluster), "FF_Y_CLUSTER_BINDING_INVALID")
            work <- dd$work
            work$.ff_time <- endpoint$time
            work$.ff_disease <- as.integer(endpoint$status == 1L)
            work$.ff_death <- as.integer(endpoint$status == 2L)
            work$.ff_weight <- w
            if (!is.null(cluster)) work$.ff_cluster <- as.character(cluster)
            fits <- baselines <- vector("list", 2L)
            warnings <- character()
            for (cause in c(1L, 2L)) {
                got <- ff_y_fit_native_one(work, dd$reference$safe_names, cause, !is.null(cluster))
                m <- got$value
                warnings <- c(warnings, got$warnings)
                ff_y_assert(!any(grepl("infinite|converg|singular|iteration", got$warnings, ignore.case = TRUE)) && 
                  inherits(m, "coxph") && !inherits(m, "coxphms") && identical(m$method, "efron") && all(is.finite(stats::coef(m))) && 
                  all(is.finite(stats::vcov(m))), "FF_Y_NONCONVERGED_OR_UNRELIABLE_NATIVE_EFRON_COX", paste(got$warnings, 
                  collapse = " | "))
                mf <- stats::model.frame(m)
                ff_y_assert(as.integer(m$n) == endpoint$n && identical(as.numeric(m$y[, "time"]), as.numeric(endpoint$time)) && 
                  identical(as.numeric(m$y[, "status"]), as.numeric(endpoint$status == cause)) && identical(as.numeric(stats::model.weights(mf)), 
                  as.numeric(w)) && isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), 
                  "FF_Y_FROZEN_TIME_EVENT_SAMPLE_WEIGHT_OR_DESIGN_CHANGED")
                if (!is.null(cluster)) ff_y_assert(!is.null(m$naive.var) && identical(as.character(stats::model.extract(mf, 
                  "cluster")), as.character(cluster)), "FF_Y_ROBUST_CLUSTER_BINDING_CHANGED")
                fits[[cause]] <- m
                baselines[[cause]] <- ff_y_native_baseline(m, dd$reference$safe_names)
            }
            names(fits) <- names(baselines) <- c("disease", "death")
            grid <- sort(unique(unlist(lapply(baselines, `[[`, "time"))))
            cumulative <- vapply(baselines, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
            hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
                , drop = FALSE])
            ff_y_assert(all(is.finite(hazard)) && all(hazard >= 0), "FF_Y_INVALID_WEIGHTED_EFRON_HAZARD_INCREMENT")
            keep <- rowSums(hazard) > 0
            engine <- list(time = grid[keep], hazard = hazard[keep, , drop = FALSE], coefficients = vapply(fits, 
                stats::coef, numeric(ncol(dd$X))), provenance = ff_y_kernel_provenance(), baseline_public_validation = baselines)
            coefficient_table <- data.frame(term = dd$reference$columns, internal_term = dd$reference$safe_names, 
                disease_beta = engine$coefficients[, 1L], death_beta = engine$coefficients[, 2L], disease_se = sqrt(diag(stats::vcov(fits[[1L]]))), 
                death_se = sqrt(diag(stats::vcov(fits[[2L]]))), row.names = NULL)
            out <- list(fit = fits, engine = engine, matrix_reference = dd$reference, coefficients = coefficient_table, 
                endpoint = endpoint, tau_report = tau_report, ties = ties, runtime = runtime, binding = list(weight = w, 
                  cluster = if (is.null(cluster)) NULL else as.character(cluster)), warnings = unique(c(warnings, 
                  unlist(lapply(baselines, `[[`, "warnings")))), method = "FULLFOLLOWUP_WEIGHTED_SEPARATE_COX_EFRON_NATIVE_HAZARDS_EXACT_REUSED_AJ", 
                weight = list(min = min(w), max = max(w), sum = sum(w), ESS = sum(w)^2/sum(w^2), interpretation = "ORIGINAL_V5_SELECTION_WEIGHTS_NO_RESCALE_NO_OUTCOME_IPCW"), 
                representative_validation_required = TRUE)
            if (isTRUE(diagnostics)) out$PH <- lapply(fits, ff_y_ph, scientific_names = dd$reference$columns, safe_names = dd$reference$safe_names)
            class(out) <- c("ff_y_matrix_fit", "ff_y_fit")
            probe <- unique(as.integer(round(seq(1, nrow(dd$X), length.out = min(16L, nrow(dd$X))))))
            nd <- dd$work[probe, , drop = FALSE]
            lp_native <- vapply(fits, function(m) as.numeric(stats::predict(m, newdata = nd, type = "lp", reference = "zero")), 
                numeric(length(probe)))
            lp_matrix <- dd$X[probe, , drop = FALSE] %*% engine$coefficients
            out$native_validation <- list(lp_max_abs = max(abs(lp_native - lp_matrix)), lp_tolerance = 1e-10, basehaz_survfit_max_abs = max(vapply(baselines, 
                `[[`, numeric(1L), "max_abs_basehaz_vs_survfit")), basehaz_survfit_tolerance = 1e-10, passed = FALSE)
            ff_y_assert(is.finite(out$native_validation$lp_max_abs) && out$native_validation$lp_max_abs <= out$native_validation$lp_tolerance, 
                "FF_Y_NATIVE_LP_MATRIX_BINDING_VALIDATION_FAILED")
            ff_predict_y_matrix(out, dd$X[probe, , drop = FALSE])
            out$native_validation$passed <- TRUE
            out
        }
        ff_fit_y <- function(data, covariates, mediator_terms, tau_report, diagnostics = FALSE, ties = "efron", 
            time_col = "time", status_col = "status", weight_col = NULL, cluster_col = NULL) {
            endpoint <- ff_y_validate_endpoint(data, tau_report, time_col, status_col)
            weights <- ff_y_validate_weights(data, weight_col, cluster_col)
            design <- ff_y_design(data, covariates, mediator_terms)
            out <- ff_fit_y_matrix(design$X, data.frame(time = endpoint$time, status = endpoint$status), weights, 
                if (is.null(cluster_col)) NULL else data[[cluster_col]], tau_report, ties, diagnostics)
            out$reference <- design$reference
            out$covariates <- covariates
            out$mediator_terms <- mediator_terms
            out$weight$column <- weight_col
            out
        }
        ff_y_probability_check <- function(p) {
            ff_y_assert(is.matrix(p) && identical(colnames(p), c("disease", "death", "survival")) && all(is.finite(p)) && 
                all(p >= 0 & p <= 1) && max(abs(rowSums(p) - 1), 0) < 1e-08, "FF_Y_INVALID_CIF_PROBABILITY_OR_CLOSURE_NO_CLIPPING")
            invisible(TRUE)
        }
        ff_y_predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
            64)) {
            ff_y_assert(inherits(object, "ff_y_fit"), "FF_Y_FITTED_OBJECT_REQUIRED")
            ff_y_assert(length(tau_report) == 1L && is.finite(tau_report) && tau_report >= 0 && tau_report <= object$endpoint$max_time, 
                "FF_Y_PREDICTION_TIME_OUTSIDE_SUPPORT")
            ff_y_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1L && length(max_prediction_mb) == 
                1L && is.finite(max_prediction_mb) && max_prediction_mb > 0, "FF_Y_INVALID_PREDICTION_MEMORY_BUDGET")
            dd <- ff_y_matrix(X, object$matrix_reference)
            actual_chunk <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
            p <- matrix(NA_real_, nrow(dd$X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
            warnings <- character()
            h <- object$engine$hazard[object$engine$time <= tau_report, , drop = FALSE]
            for (start in seq.int(1L, nrow(dd$X), by = actual_chunk)) {
                ii <- seq.int(start, min(nrow(dd$X), start + actual_chunk - 1L))
                nd <- dd$work[ii, , drop = FALSE]
                lp <- vapply(object$fit, function(m) {
                  got <- ff_y_capture(stats::predict(m, newdata = nd, type = "lp", reference = "zero"))
                  warnings <<- c(warnings, got$warnings)
                  as.numeric(got$value)
                }, numeric(length(ii)))
                lp <- matrix(lp, nrow = length(ii), ncol = 2L)
                rr <- exp(lp)
                ff_y_assert(all(is.finite(rr)) && all(rr > 0), "FF_Y_NONFINITE_NATIVE_RELATIVE_RISK_NO_CLIPPING")
                if (exists("ff_aj_status", mode = "function") && isTRUE(ff_aj_status()$ready)) {
                  p[ii, ] <- ff_aj_backend(h, rr)
                } else {
                  ff_y_assert(!isTRUE(getOption("ff.y.require_compiled_aj", FALSE)), "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
                  p[ii, ] <- ff_y_aj(h, cbind(rr, rr), rep(1L, nrow(h)))
                }
                ff_y_probability_check(p[ii, , drop = FALSE])
            }
            ff_y_probability_check(p)
            attr(p, "warnings") <- unique(warnings)
            attr(p, "chunk_size") <- actual_chunk
            attr(p, "product.limit") <- TRUE
            attr(p, "risk_range") <- apply(p, 2L, range)
            p
        }
        ff_predict_y_matrix <- function(object, X, tau_report = object$tau_report, chunk_size = getOption("ff.y.chunk_size", 
            4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 64)) {
            ff_y_assert(identical(as.numeric(tau_report), as.numeric(object$tau_report)), "FF_Y_REPORT_TIME_DRIFT_NO_TIMEPOINT_SEARCH")
            ff_y_predict_at(object, X, tau_report, chunk_size, max_prediction_mb)
        }
        ff_predict_y <- function(object, newdata, tau_report = object$tau_report, chunk_size = getOption("ff.y.chunk_size", 
            4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 64)) {
            ff_y_assert(!is.null(object$reference), "FF_Y_RAW_PREDICT_REQUIRES_RAW_FIT_REFERENCE")
            dd <- ff_y_design(newdata, object$covariates, object$mediator_terms, object$reference)
            ff_predict_y_matrix(object, dd$X, tau_report, chunk_size, max_prediction_mb)
        }
        ff_y_effect_closure <- function(risks, tolerance = 1e-12) {
            ff_y_assert(is.numeric(risks) && all(c("R00", "R10", "R11") %in% names(risks)) && all(is.finite(risks)) && 
                all(risks >= 0 & risks <= 1), "FF_Y_INVALID_INTEGRATED_SCENARIO_RISKS")
            TE <- unname(risks["R11"] - risks["R00"])
            IDE <- unname(risks["R10"] - risks["R00"])
            IIE <- unname(risks["R11"] - risks["R10"])
            difference <- abs(TE - IDE - IIE)
            ff_y_assert(difference <= tolerance, "FF_Y_INTEGRATED_EFFECT_CLOSURE_FAILED")
            list(TE = TE, IDE = IDE, IIE = IIE, PM = IIE/TE, max_abs = difference, tolerance = tolerance, tolerance_basis = "DOUBLE_PRECISION_ADDITION_1E_12_ABSOLUTE_PROBABILITY", 
                passed = TRUE)
        }
        ff_validate_y_matrix <- function(object, X, risk, w, cluster = NULL, require_rr = TRUE) {
            coefficient_tolerance <- cif_tolerance <- 1e-09
            dd <- ff_y_matrix(X, object$matrix_reference)
            endpoint <- ff_y_risk(risk, object$tau_report)
            ff_y_assert(identical(as.numeric(endpoint$time), as.numeric(object$endpoint$time)) && identical(endpoint$status, 
                object$endpoint$status) && identical(as.numeric(w), as.numeric(object$binding$weight)) && identical(if (is.null(cluster)) NULL else as.character(cluster), 
                object$binding$cluster), "FF_Y_VALIDATION_INPUT_BINDING_DRIFT")
            checks <- list()
            add <- function(name, difference, tolerance, reference) {
                checks[[length(checks) + 1L]] <<- data.frame(check = name, max_abs_difference = difference, tolerance = tolerance, 
                  passed = is.finite(difference) && difference <= tolerance, reference = reference)
            }
            reference <- ff_fit_y_matrix(dd$X, risk, w, cluster, object$tau_report)
            add("native_weighted_coefficient", max(abs(object$engine$coefficients - reference$engine$coefficients)), 
                coefficient_tolerance, "INDEPENDENT_PUBLIC_SURVIVAL_COXPH_EFRON_SAME_ORIGINAL_WEIGHTS")
            for (j in 1:2) {
                a <- object$engine$baseline_public_validation[[j]]
                b <- reference$engine$baseline_public_validation[[j]]
                ff_y_assert(identical(a$time, b$time), "FF_Y_VALIDATION_NATIVE_GRID_CHANGED")
                add(paste0("native_weighted_cumulative_hazard_cause", j), max(abs(a$cumulative - b$cumulative)), 
                  1e-10, "INDEPENDENT_PUBLIC_SURVIVAL_BASEHAZ_VS_PUBLIC_SURVFIT_EFRON")
            }
            probe <- unique(as.integer(round(seq(1, nrow(dd$X), length.out = min(16L, nrow(dd$X))))))
            grid <- unique(c(0, object$tau_report * c(0.25, 0.5, 0.75, 1)))
            curve <- lapply(grid, function(t) ff_y_predict_at(object, dd$X[probe, , drop = FALSE], t))
            ranges <- do.call(rbind, lapply(seq_along(grid), function(k) data.frame(time = grid[k], state = colnames(curve[[k]]), 
                min_probability = apply(curve[[k]], 2L, min), max_probability = apply(curve[[k]], 2L, max))))
            add("three_state_raw_probability_range", max(0, -ranges$min_probability, ranges$max_probability - 1), 
                0, "STRICT_ZERO_ONE_NO_CLIPPING_NO_EPSILON")
            add("disease_death_survival_sum", max(vapply(curve, function(p) max(abs(rowSums(p) - 1)), numeric(1L))), 
                1e-09, "THREE_STATE_PROBABILITY_CLOSURE_NO_CLIPPING")
            monotonic_violation <- 0
            for (k in 2:length(curve)) {
                delta <- curve[[k]] - curve[[k - 1L]]
                monotonic_violation <- max(monotonic_violation, -delta[, "disease"], -delta[, "death"], delta[, 
                  "survival"])
            }
            add("three_state_time_monotonicity", monotonic_violation, 1e-12, "PRESPECIFIED_REPORT_TIME_QUARTILES_NOT_EFFECT_TIMEPOINT_SELECTION")
            rr_ready <- requireNamespace("riskRegression", quietly = TRUE) && requireNamespace("prodlim", quietly = TRUE)
            rr_status <- "NOT_EXECUTED_DEPENDENCY_MISSING"
            if (rr_ready) {
                unit <- if (all(w == 1) && is.null(cluster)) object else ff_fit_y_matrix(dd$X, risk, rep(1, endpoint$n), 
                  NULL, object$tau_report)
                work <- dd$work
                work$.ff_time <- endpoint$time
                work$.ff_status <- endpoint$status
                f <- stats::reformulate(names(dd$work), "prodlim::Hist(.ff_time, .ff_status)")
                got <- ff_y_capture(riskRegression::CSC(f, data = work, cause = 1, surv.type = "hazard", fitter = "coxph", 
                  ties = "efron", model = TRUE, singular.ok = FALSE, robust = FALSE, na.action = stats::na.fail, 
                  control = survival::coxph.control(iter.max = 50L)))
                ff_y_assert(!length(got$warnings), "FF_Y_UNIT_RR_REFERENCE_FIT_WARNING", paste(got$warnings, collapse = " | "))
                public <- got$value
                for (j in 1:2) add(paste0("unit_coefficient_vs_RR_cause", j), max(abs(stats::coef(unit$fit[[j]]) - 
                  stats::coef(public$models[[paste("Cause", j)]]))), coefficient_tolerance, "INDEPENDENT_UNIT_WEIGHT_PUBLIC_RISKREGRESSION_CSC_EFRON")
                for (t in grid[-1L]) {
                  p <- ff_y_predict_at(unit, dd$X[probe, , drop = FALSE], t)
                  for (j in 1:3) {
                    type <- if (j == 3L) "survival" else "absRisk"
                    got <- ff_y_capture(stats::predict(public, newdata = dd$work[probe, , drop = FALSE], times = t, 
                      cause = if (j == 2L) 2 else 1, type = type, product.limit = TRUE, se = FALSE, band = FALSE, 
                      iid = FALSE, confint = FALSE))
                    ff_y_assert(!length(got$warnings), "FF_Y_UNIT_RR_REFERENCE_PREDICT_WARNING", paste(got$warnings, 
                      collapse = " | "))
                    add(paste0("unit_", colnames(p)[j], "_vs_RR_t", format(t)), max(abs(p[, j] - as.numeric(got$value[[type]]))), 
                      cif_tolerance, "PUBLIC_RISKREGRESSION_PREDICT_CAUSESPECIFICCOX_PRODUCT_LIMIT_TRUE")
                  }
                }
                rr_status <- "EXECUTED"
            }
            table <- do.call(rbind, checks)
            passed <- all(table$passed) && (rr_ready || !require_rr)
            out <- list(passed = passed, certified_for_representative_gate = all(table$passed) && rr_ready, native_checks_passed = all(table$passed), 
                checks = table, risk_ranges_by_time = ranges, riskRegression_unit_reference_status = rr_status, 
                runtime = ff_y_runtime(), kernel_provenance = ff_y_kernel_provenance(), tolerance_basis = paste("ABSOLUTE_DOUBLE_PRECISION_NUMERICAL_EQUIVALENCE:", 
                  "coefficients/CIF 1e-9; public cumhaz/LP 1e-10; monotonicity/effect closure 1e-12;", "strict raw ranges and jump bounds; no clipping/epsilon; not statistical/P-value gates"), 
                integrated_effect_validation = "CALL_ff_y_effect_closure_ON_ACTUAL_A_INTEGRATED_R00_R10_R11", representative_scope = "PRESPECIFIED_ENGINEERING_COMBINATIONS_POINT_AND_B1_20_NOT_NEW_DRAWS")
            out$compiled_backend_acceptance <- if (exists("ff_aj_status", mode = "function")) ff_aj_status() else list(ready = FALSE, 
                state = "FROZEN_R_REFERENCE_ONLY_SYNTHETIC_CONTEXT")
            out$certified_for_representative_gate <- isTRUE(out$certified_for_representative_gate) && isTRUE(out$compiled_backend_acceptance$ready) && 
                isTRUE(out$compiled_backend_acceptance$acceptance$passed)
            if (isTRUE(getOption("ff.y.require_compiled_aj", FALSE))) {
                ff_y_assert(isTRUE(out$compiled_backend_acceptance$ready), "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
            }
            if (!all(table$passed) || (require_rr && !rr_ready)) {
                message <- if (!all(table$passed)) paste0("FF_Y_REPRESENTATIVE_NUMERICAL_VALIDATION_FAILED: ", 
                  paste(table$check[!table$passed], collapse = ", ")) else "FF_Y_REPRESENTATIVE_UNIT_RR_VALIDATION_DEPENDENCY_MISSING"
                stop(structure(list(message = message, call = NULL, validation = out), class = c("ff_y_validation_error", 
                  "error", "condition")))
            }
            out
        }
    }), envir = context)
}

