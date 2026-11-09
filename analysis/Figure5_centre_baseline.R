# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ff_ce_method_spec <- function() {
    text <- paste(c("CENTER_ETHNICITY_V4_20261028", "endpoints=250,250.2,411,411.3,411.4;tau=11;full_followup_unchanged", 
        "disease_411.3=center_strata+Asian_Chinese_fixed_pool", "disease_other=EXACT_ORIGINAL_FROZEN_DESIGN", "death_all=center_strata+Asian_Chinese_fixed_pool", 
        "center_strata=exact_mutually_exclusive_frozen_treatment_dummies_plus_reference", "A+sample+draw+selection_weights+Efron+robust_cluster=UNCHANGED", 
        "every_observed_center_requires_max_followup_ge_11;zero_event_supported_baseline=0", "public_coxph+basehaz+survfit+predict_zero;original_AJ_kernel_no_clipping", 
        "one_retry_only_for_native_infinite_warning;init=original_beta;eps=1e-11;iter.max=50;toler.inf=sqrt(1e-9)", 
        "original_attempt_audit_retained;unresolved_failure_remains_failure;no_significance_selection"), collapse = "\n")
    list(id = "CENTER_ETHNICITY_V4_20261028", canonical = text, sha256 = digest::digest(text, algo = "sha256", 
        serialize = FALSE))
}

ff_ce_fail <- function(code, receipt = NULL) stop(structure(list(message = code, call = NULL, numerical_audit = receipt), 
    class = c("ff_ce_error", "error", "condition")))

ff_ce_bad_warning <- function(w) any(grepl("infinite|converg|singular|iteration", w, ignore.case = TRUE))

ff_ce_centers <- function(X, reference = NULL, centre_actual = NULL) {
    columns <- grep("^centerc[0-9]+$", colnames(X), value = TRUE)
    ff_y_assert(length(columns) > 0L && !length(setdiff(grep("^center", colnames(X), value = TRUE), columns)), 
        "CE_FROZEN_CENTER_COLUMN_MAPPING_NOT_EXACT")
    if (!is.null(reference)) 
        ff_y_assert(identical(columns, reference$columns), "CE_CENTER_COLUMN_DRIFT")
    z <- X[, columns, drop = FALSE]
    ff_y_assert(all(z %in% 0:1) && all(rowSums(z) <= 1), "CE_CENTER_DUMMIES_NOT_MUTUALLY_EXCLUSIVE")
    labels <- c("center_REFERENCE_FROZEN", columns)
    index <- ifelse(rowSums(z) == 0L, 1L, max.col(z, ties.method = "first") + 1L)
    group <- labels[index]
    if (!is.null(centre_actual)) {
        ff_y_assert(length(centre_actual) == nrow(X) && !anyNA(centre_actual), "CE_ACTUAL_CENTER_BINDING_INVALID")
        actual <- as.character(centre_actual)
        ff_y_assert(all(vapply(split(actual, group), function(v) length(unique(v)) == 1L, logical(1L))) && all(vapply(split(group, 
            actual), function(v) length(unique(v)) == 1L, logical(1L))), "CE_ACTUAL_CENTER_GROUPING_DISAGREES_WITH_FROZEN_X")
    }
    if (is.null(reference)) 
        reference <- list(columns = columns, all_design_levels = labels, observed_levels = labels[labels %in% group], 
            source = "EXACT_FROZEN_CENTER_TREATMENT_DUMMIES_NO_ROW_DELETION")
    ff_y_assert(all(group %in% reference$observed_levels), "CE_PREDICTION_CENTER_NOT_OBSERVED_IN_FIT")
    list(group = group, reference = reference, actual_center_grouping_checked = !is.null(centre_actual))
}

ff_ce_columns <- function(X, center_columns) {
    asian <- "ethnicityAsian or Asian British"
    chinese <- "ethnicityChinese"
    ff_y_assert(all(c(asian, chinese) %in% colnames(X)), "CE_EXPECTED_FROZEN_ETHNICITY_COLUMNS_ABSENT")
    ff_y_assert(all(X[, asian] %in% 0:1) && all(X[, chinese] %in% 0:1) && all(X[, asian] + X[, chinese] <= 1), 
        "CE_ETHNICITY_DUMMIES_NOT_DISJOINT")
    z <- X[, !colnames(X) %in% c(center_columns, chinese), drop = FALSE]
    z[, asian] <- X[, asian] + X[, chinese]
    colnames(z)[colnames(z) == asian] <- "ethnicityAsian_or_Chinese_pooled"
    z
}

ff_ce_support <- function(ep, centers, tau) do.call(rbind, lapply(centers$reference$observed_levels, function(g) {
    keep <- centers$group == g
    data.frame(center_design_label = g, N = sum(keep), disease_events = sum(ep$status[keep] == 1L), death_events = sum(ep$status[keep] == 
        2L), max_followup = max(ep$time[keep]), at_risk_at_tau = sum(ep$time[keep] >= tau), supported_at_tau = any(ep$time[keep] >= 
        tau), stringsAsFactors = FALSE)
}))

ff_ce_native <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL) {
    if (is.null(init)) 
        init <- rep(0, length(safe_names))
    strata <- survival::strata
    f <- stats::reformulate(c(safe_names, if (stratified) "strata(.ce_center)"), response = if (cause == 1L) 
        "survival::Surv(.ff_time,.ff_disease)"
    else "survival::Surv(.ff_time,.ff_death)", env = environment())
    if (cluster_bound) 
        ff_y_capture(survival::coxph(f, data = work, weights = .ff_weight, cluster = .ff_cluster, robust = TRUE, 
            ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE, 
            control = control, init = init))
    else ff_y_capture(survival::coxph(f, data = work, weights = .ff_weight, robust = FALSE, ties = "efron", singular.ok = FALSE, 
        na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE, control = control, init = init))
}

ff_ce_score <- function(fit) {
    detail <- survival::coxph.detail(fit, riskmat = FALSE)
    p <- length(stats::coef(fit))
    u <- if (p == 1L) 
        sum(detail$score)
    else colSums(detail$score)
    I <- if (p == 1L) 
        matrix(sum(detail$imat), 1L)
    else apply(detail$imat, c(1L, 2L), sum)
    V <- if (is.null(fit$naive.var)) 
        stats::vcov(fit)
    else fit$naive.var
    step <- as.numeric(u %*% V)
    eig <- eigen((I + t(I))/2, symmetric = TRUE, only.values = TRUE)$values
    list(score_max_abs = max(abs(u)), newton_step_max_abs = max(abs(step)), newton_step_over_SE_max_abs = max(abs(step/sqrt(diag(V)))), 
        min_information_eigenvalue = min(eig), information_condition = max(eig)/min(eig), inverse_information_vs_naive_V_max_abs = max(abs(solve(I) - 
            V)), public_recomputed_not_bitwise_native_terminal_score = TRUE)
}

ff_ce_baselines <- function(fit, safe_names, stratified, center_levels, support, cause) {
    if (!stratified) 
        return(list(GLOBAL = ff_y_native_baseline(fit, safe_names)))
    bh <- ff_y_capture(survival::basehaz(fit, centered = FALSE))
    ff_y_assert(!ff_ce_bad_warning(bh$warnings) && "strata" %in% names(bh$value), "CE_STRATIFIED_BASEHAZ_FAILED")
    labels <- as.character(bh$value$strata)
    ff_y_assert(setequal(unique(labels), center_levels), "CE_NATIVE_BASELINE_STRATA_LABEL_BINDING_FAILED")
    result <- lapply(center_levels, function(g) {
        b <- bh$value[labels == g, , drop = FALSE]
        zero <- as.data.frame(stats::setNames(as.list(rep(0, length(safe_names))), safe_names))
        zero$.ce_center <- factor(g, levels = center_levels)
        sf <- ff_y_capture(survival::survfit(fit, newdata = zero, se.fit = FALSE, conf.type = "none", stype = 2, 
            ctype = 2))
        difference <- if (identical(as.numeric(b$time), as.numeric(sf$value$time))) 
            max(abs(b$hazard - as.numeric(sf$value$cumhaz)), 0)
        else Inf
        ff_y_assert(!ff_ce_bad_warning(sf$warnings) && is.finite(difference) && difference <= 1e-10 && all(is.finite(b$hazard)) && 
            all(diff(c(0, b$hazard)) >= 0), "CE_NATIVE_BASEHAZ_SURVFIT_STRATUM_VALIDATION_FAILED", g)
        s <- support[support$center_design_label == g, , drop = FALSE]
        events <- if (cause == 1L) 
            s$disease_events
        else s$death_events
        if (events == 0L) 
            ff_y_assert(all(b$hazard == 0) && s$supported_at_tau, "CE_ZERO_EVENT_STRATUM_NOT_NATIVE_ZERO_OR_SUPPORTED")
        list(time = b$time, cumulative = b$hazard, max_abs_basehaz_vs_survfit = difference, tolerance = 1e-10, 
            warnings = unique(c(bh$warnings, sf$warnings)), zero_event_supported = events == 0L, support = list(N = s$N, 
                max_followup = s$max_followup, at_risk_at_tau = s$at_risk_at_tau, events = events))
    })
    stats::setNames(result, center_levels)
}

ff_ce_merge_baselines <- function(b) {
    grid <- sort(unique(unlist(lapply(b, `[[`, "time"))))
    ff_y_assert(length(grid) > 0L, "CE_NO_NATIVE_BASELINE_GRID")
    cumulative <- vapply(b, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
    cumulative <- matrix(cumulative, ncol = 2L)
    hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
        , drop = FALSE])
    ff_y_assert(all(is.finite(hazard)) && all(hazard >= 0), "CE_INVALID_NATIVE_HAZARD_INCREMENT")
    keep <- rowSums(hazard) > 0
    list(time = grid[keep], hazard = hazard[keep, , drop = FALSE])
}

ff_ce_fit <- function(X, risk, w, cluster = NULL, tau_report, ties = "efron", diagnostics = FALSE, endpoint_id, 
    centre_actual = NULL, audit_callback = NULL) {
    ff_y_assert(endpoint_id %in% c("250", "250.2", "411", "411.3", "411.4") && length(endpoint_id) == 1L, "CE_FIXED_FIVE_ENDPOINT_SCOPE_REQUIRED")
    ff_y_assert(identical(as.numeric(tau_report), 11) && identical(ties, "efron"), "CE_FIXED_TAU11_AND_EFRON_REQUIRED")
    runtime <- ff_y_runtime(require = TRUE)
    original <- ff_y_matrix(X)
    ep <- ff_y_risk(risk, tau_report)
    ff_y_assert(nrow(original$X) == ep$n && is.numeric(w) && length(w) == ep$n && all(is.finite(w) & w > 0), "CE_ORIGINAL_WEIGHT_SAMPLE_INVALID")
    if (!is.null(cluster)) 
        ff_y_assert(length(cluster) == ep$n && !anyNA(cluster), "CE_CLUSTER_BINDING_INVALID")
    centers <- ff_ce_centers(original$X, centre_actual = centre_actual)
    support <- ff_ce_support(ep, centers, tau_report)
    ff_y_assert(all(support$supported_at_tau), "CE_CENTER_TAU11_SUPPORT_ABSENT_NO_EXTRAPOLATION")
    changed <- ff_y_matrix(ff_ce_columns(original$X, centers$reference$columns))
    stratified <- c(disease = identical(endpoint_id, "411.3"), death = TRUE)
    designs <- list(disease = if (stratified[1L]) changed else original, death = changed)
    fits <- baselines <- audits <- vector("list", 2L)
    warnings <- character()
    lp_difference <- numeric(2L)
    names(fits) <- names(baselines) <- names(audits) <- names(designs)
    probe <- unique(as.integer(round(seq(1, ep$n, length.out = min(16L, ep$n)))))
    ctrl0 <- survival::coxph.control(iter.max = 50L)
    ctrl1 <- survival::coxph.control(eps = 1e-11, iter.max = 50L, toler.inf = ctrl0$toler.inf)
    emit <- function(rec) if (is.function(audit_callback)) 
        audit_callback(rec)
    for (j in 1:2) {
        dd <- designs[[j]]
        work <- dd$work
        work$.ce_center <- factor(centers$group, levels = centers$reference$observed_levels)
        work$.ff_time <- ep$time
        work$.ff_disease <- as.integer(ep$status == 1L)
        work$.ff_death <- as.integer(ep$status == 2L)
        work$.ff_weight <- w
        if (!is.null(cluster)) 
            work$.ff_cluster <- as.character(cluster)
        binding_hash <- digest::digest(list(X = dd$X, time = ep$time, status = ep$status, weights = w, center = centers$group, 
            cluster = cluster), algo = "sha256")
        got <- ff_ce_native(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), ctrl0)
        m <- got$value
        rec <- list(cause = names(designs)[j], binding_hash = binding_hash, original_control = ctrl0, original_warnings = got$warnings, 
            original_coefficients = stats::coef(m), original_vcov = stats::vcov(m), original_iterations = m$iter, 
            original_final_loglik = tail(m$loglik, 1L), retry_attempted = FALSE)
        emit(c(rec, list(stage = "ORIGINAL_ATTEMPT_PRESERVED")))
        finite <- all(is.finite(stats::coef(m))) && all(is.finite(stats::vcov(m)))
        if (ff_ce_bad_warning(got$warnings)) {
            eligible <- length(got$warnings) > 0L && finite && all(grepl("Loglik converged before variable.*coefficient may be infinite", 
                got$warnings))
            if (!eligible) 
                ff_ce_fail("CE_UNRESOLVED_ORIGINAL_NATIVE_WARNING_NO_ELIGIBLE_RETRY", rec)
            old <- m
            rec$retry_attempted <- TRUE
            rec$retry_control <- ctrl1
            rec$original_score_information <- ff_ce_score(old)
            got <- ff_ce_native(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), ctrl1, stats::coef(old))
            m <- got$value
            rec$retry_warnings <- got$warnings
            rec$retry_iterations <- m$iter
            rec$retry_final_loglik <- tail(m$loglik, 1L)
            rec$warm_start_initial_statistics_are_not_null_statistics <- TRUE
            rec$coefficient_max_abs_difference <- max(abs(stats::coef(m) - stats::coef(old)))
            rec$coefficient_difference_over_old_SE_max_abs <- max(abs((stats::coef(m) - stats::coef(old))/sqrt(diag(stats::vcov(old)))))
            rec$SE_max_abs_difference <- max(abs(sqrt(diag(stats::vcov(m))) - sqrt(diag(stats::vcov(old)))))
            rec$final_loglik_difference <- tail(m$loglik, 1L) - tail(old$loglik, 1L)
            rec$LP_max_abs_difference <- max(abs(as.numeric(predict(m, newdata = work, type = "lp", reference = "zero")) - 
                as.numeric(predict(old, newdata = work, type = "lp", reference = "zero"))))
            rec$retry_score_information <- ff_ce_score(m)
            emit(c(rec, list(stage = "SINGLE_NUMERICAL_RETRY_COMPLETED")))
            if (ff_ce_bad_warning(got$warnings)) 
                ff_ce_fail("CE_UNRESOLVED_NATIVE_WARNING_AFTER_SINGLE_RETRY", rec)
        }
        ff_y_assert(inherits(m, "coxph") && !inherits(m, "coxphms") && identical(m$method, "efron") && all(is.finite(stats::coef(m))) && 
            all(is.finite(stats::vcov(m))), "CE_NATIVE_NONFINITE_OR_WRONG_ESTIMATOR")
        mf <- stats::model.frame(m)
        ff_y_assert(m$n == ep$n && identical(as.numeric(m$y[, "time"]), as.numeric(ep$time)) && identical(as.numeric(m$y[, 
            "status"]), as.numeric(ep$status == j)) && identical(as.numeric(stats::model.weights(mf)), as.numeric(w)) && 
            isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), "CE_FROZEN_TIME_EVENT_WEIGHT_OR_DESIGN_CHANGED")
        if (stratified[j]) 
            ff_y_assert(identical(as.character(mf[["strata(.ce_center)"]]), centers$group), "CE_NATIVE_STRATUM_BINDING_CHANGED")
        if (!is.null(cluster)) 
            ff_y_assert(!is.null(m$naive.var) && identical(as.character(stats::model.extract(mf, "cluster")), as.character(cluster)), 
                "CE_ROBUST_CLUSTER_BINDING_CHANGED")
        baselines[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, stratified[j], centers$reference$observed_levels, 
            support, j)
        if (isTRUE(rec$retry_attempted)) {
            oldbase <- ff_ce_baselines(old, dd$reference$safe_names, stratified[j], centers$reference$observed_levels, 
                support, j)
            rec$baseline_max_abs_difference <- max(vapply(names(oldbase), function(g) {
                ff_y_assert(identical(oldbase[[g]]$time, baselines[[j]][[g]]$time), "CE_RETRY_BASELINE_GRID_CHANGED")
                max(abs(oldbase[[g]]$cumulative - baselines[[j]][[g]]$cumulative), 0)
            }, numeric(1L)))
        }
        rec$binding_hash_after <- digest::digest(list(X = dd$X, time = ep$time, status = ep$status, weights = w, 
            center = centers$group, cluster = cluster), algo = "sha256")
        ff_y_assert(identical(binding_hash, rec$binding_hash_after), "CE_RETRY_BINDING_MUTATION")
        emit(c(rec, list(stage = "ACCEPTED_NATIVE_ATTEMPT_NOT_FULL_PIPELINE_VALIDITY")))
        fits[[j]] <- m
        audits[[j]] <- rec
        warnings <- c(warnings, got$warnings)
        nd <- work[probe, , drop = FALSE]
        lp <- ff_y_capture(as.numeric(predict(m, newdata = nd, type = "lp", reference = "zero")))
        ff_y_assert(!ff_ce_bad_warning(lp$warnings), "CE_NATIVE_LP_WARNING")
        lp_difference[j] <- max(abs(lp$value - as.numeric(dd$X[probe, , drop = FALSE] %*% stats::coef(m))))
        ff_y_assert(lp_difference[j] <= 1e-10, "CE_NATIVE_LP_MATRIX_BINDING_FAILED")
    }
    engines <- stats::setNames(lapply(centers$reference$observed_levels, function(g) ff_ce_merge_baselines(list(if (stratified[1L]) baselines$disease[[g]] else baselines$disease$GLOBAL, 
        baselines$death[[g]]))), centers$reference$observed_levels)
    coefficient_table <- do.call(rbind, lapply(1:2, function(j) data.frame(cause = names(fits)[j], term = designs[[j]]$reference$columns, 
        internal_term = designs[[j]]$reference$safe_names, coefficient = unname(stats::coef(fits[[j]])), SE = sqrt(diag(stats::vcov(fits[[j]]))), 
        row.names = NULL)))
    all_baselines <- unlist(baselines, recursive = FALSE)
    out <- list(fit = fits, engine = list(by_center = engines, baseline_public_validation = baselines, provenance = ff_y_kernel_provenance()), 
        matrix_reference = original$reference, cause_references = lapply(designs, `[[`, "reference"), center_reference = centers$reference, 
        center_support = support, stratum_support = support, centre_actual_grouping_checked = centers$actual_center_grouping_checked, 
        stratified = stratified, coefficients = coefficient_table, endpoint = ep, tau_report = tau_report, endpoint_id = endpoint_id, 
        ties = ties, runtime = runtime, binding = list(weight = w, cluster = cluster), numerical_audit = audits, 
        numerical_retry = audits, method = ff_ce_method_spec()$id, method_spec = ff_ce_method_spec(), method_spec_hash = ff_ce_method_spec()$sha256, 
        warnings = unique(c(warnings, unlist(lapply(all_baselines, `[[`, "warnings")))), weight = list(min = min(w), 
            max = max(w), sum = sum(w), ESS = sum(w)^2/sum(w^2), interpretation = "UNCHANGED_ORIGINAL_SELECTION_WEIGHTS"), 
        representative_validation_required = TRUE, formal_ready = FALSE, native_validation = list(passed = FALSE, 
            lp_max_abs = max(lp_difference), basehaz_survfit_max_abs = max(vapply(all_baselines, `[[`, numeric(1L), 
                "max_abs_basehaz_vs_survfit"))))
    class(out) <- c("ff_ce_matrix_fit", "ff_y_matrix_fit", "ff_y_fit")
    if (isTRUE(diagnostics)) 
        out$PH <- lapply(1:2, function(j) ff_y_ph(fits[[j]], designs[[j]]$reference$columns, designs[[j]]$reference$safe_names))
    curve <- lapply(c(0, 0.25, 0.5, 0.75, 1) * tau_report, function(t) ff_ce_predict_at(out, original$X[probe, 
        , drop = FALSE], t))
    violation <- 0
    for (k in 2:length(curve)) {
        delta <- curve[[k]] - curve[[k - 1L]]
        violation <- max(violation, -delta[, "disease"], -delta[, "death"], delta[, "survival"])
    }
    ff_y_assert(violation <= 1e-12, "CE_CIF_MONOTONICITY_FAILED")
    out$native_validation$monotonicity_violation <- violation
    out$native_validation$three_state_closure <- max(vapply(curve, function(z) max(abs(rowSums(z) - 1)), numeric(1L)))
    out$native_validation$passed <- TRUE
    out
}

ff_ce_predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
    64)) {
    ff_y_assert(inherits(object, "ff_ce_matrix_fit"), "CE_FITTED_OBJECT_REQUIRED")
    ff_y_assert(length(tau_report) == 1L && is.finite(tau_report) && tau_report >= 0 && tau_report <= object$tau_report, 
        "CE_TIME_OUTSIDE_FIXED_TAU11_SUPPORT")
    ff_y_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1L && length(max_prediction_mb) == 
        1L && is.finite(max_prediction_mb) && max_prediction_mb > 0, "CE_PREDICTION_MEMORY_BUDGET_INVALID")
    original <- ff_y_matrix(X, object$matrix_reference)
    centers <- ff_ce_centers(original$X, object$center_reference)
    changed <- ff_ce_columns(original$X, object$center_reference$columns)
    designs <- list(disease = ff_y_matrix(if (object$stratified[1L]) changed else original$X, object$cause_references$disease), 
        death = ff_y_matrix(changed, object$cause_references$death))
    step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
    p <- matrix(NA_real_, nrow(original$X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
    warnings <- character()
    for (g in unique(centers$group)) {
        support <- object$center_support[object$center_support$center_design_label == g, , drop = FALSE]
        ff_y_assert(nrow(support) == 1L && support$max_followup >= tau_report, "CE_PREDICTION_STRATUM_UNSUPPORTED_NO_EXTRAPOLATION")
        rows <- which(centers$group == g)
        eng <- object$engine$by_center[[g]]
        h <- eng$hazard[eng$time <= tau_report, , drop = FALSE]
        for (start in seq.int(1L, length(rows), by = step)) {
            ii <- rows[seq.int(start, min(length(rows), start + step - 1L))]
            lp <- vapply(1:2, function(j) {
                nd <- designs[[j]]$work[ii, , drop = FALSE]
                nd$.ce_center <- factor(g, levels = object$center_reference$observed_levels)
                got <- ff_y_capture(as.numeric(predict(object$fit[[j]], newdata = nd, type = "lp", reference = "zero")))
                warnings <<- c(warnings, got$warnings)
                got$value
            }, numeric(length(ii)))
            rr <- exp(matrix(lp, nrow = length(ii), ncol = 2L))
            ff_y_assert(!ff_ce_bad_warning(warnings) && all(is.finite(rr)) && all(rr > 0), "CE_INVALID_NATIVE_RR_OR_WARNING_NO_CLIPPING")
            if (exists("ff_aj_status", mode = "function") && isTRUE(ff_aj_status()$ready)) 
                p[ii, ] <- ff_aj_backend(h, rr)
            else {
                ff_y_assert(!isTRUE(getOption("ff.y.require_compiled_aj", FALSE)), "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
                p[ii, ] <- ff_y_aj(h, cbind(rr, rr), rep(1L, nrow(h)))
            }
            ff_y_probability_check(p[ii, , drop = FALSE])
        }
    }
    ff_y_probability_check(p)
    attr(p, "warnings") <- unique(warnings)
    attr(p, "chunk_size") <- step
    attr(p, "product.limit") <- TRUE
    attr(p, "risk_range") <- apply(p, 2L, range)
    p
}

ff_ce_predict_y_matrix <- function(object, X, tau_report = object$tau_report, chunk_size = getOption("ff.y.chunk_size", 
    4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 64)) {
    ff_y_assert(identical(as.numeric(tau_report), as.numeric(object$tau_report)), "CE_REPORT_TIME_DRIFT")
    ff_ce_predict_at(object, X, tau_report, chunk_size, max_prediction_mb)
}

ff_ce_install <- function(endpoint_id, envir = parent.frame(), audit_callback = NULL) {
    force(endpoint_id)
    force(envir)
    force(audit_callback)
    ff_y_assert(endpoint_id %in% c("250", "250.2", "411", "411.3", "411.4"), "CE_FIXED_ENDPOINT_REQUIRED")
    if (!exists(".ff_ce_original_predict_at", envir = envir, inherits = FALSE)) {
        assign(".ff_ce_original_predict_at", get("ff_y_predict_at", envir = envir, inherits = TRUE), envir = envir)
        assign(".ff_ce_original_fit_y_matrix", get("ff_fit_y_matrix", envir = envir, inherits = TRUE), envir = envir)
    }
    previous_predict <- get(".ff_ce_original_predict_at", envir = envir, inherits = FALSE)
    fit <- function(X, risk, w, cluster = NULL, tau_report, ties = "efron", diagnostics = FALSE, centre_actual = NULL) ff_ce_fit(X, 
        risk, w, cluster, tau_report, ties, diagnostics, endpoint_id, centre_actual, audit_callback)
    predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
        64)) {
        if (inherits(object, "ff_ce_matrix_fit")) 
            ff_ce_predict_at(object, X, tau_report, chunk_size, max_prediction_mb)
        else previous_predict(object, X, tau_report, chunk_size, max_prediction_mb)
    }
    assign("ff_fit_y_matrix", fit, envir = envir)
    assign("ff_y_predict_at", predict_at, envir = envir)
    invisible(list(endpoint_id = endpoint_id, method_spec = ff_ce_method_spec(), formal_ready = FALSE, installed_functions = c("ff_fit_y_matrix", 
        "ff_y_predict_at")))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ff_ce_method_spec <- function() {
            text <- paste(c("CENTER_ETHNICITY_V4_20261028", "endpoints=250,250.2,411,411.3,411.4;tau=11;full_followup_unchanged", 
                "disease_411.3=center_strata+Asian_Chinese_fixed_pool", "disease_other=EXACT_ORIGINAL_FROZEN_DESIGN", 
                "death_all=center_strata+Asian_Chinese_fixed_pool", "center_strata=exact_mutually_exclusive_frozen_treatment_dummies_plus_reference", 
                "A+sample+draw+selection_weights+Efron+robust_cluster=UNCHANGED", "every_observed_center_requires_max_followup_ge_11;zero_event_supported_baseline=0", 
                "public_coxph+basehaz+survfit+predict_zero;original_AJ_kernel_no_clipping", "one_retry_only_for_native_infinite_warning;init=original_beta;eps=1e-11;iter.max=50;toler.inf=sqrt(1e-9)", 
                "original_attempt_audit_retained;unresolved_failure_remains_failure;no_significance_selection"), 
                collapse = "\n")
            list(id = "CENTER_ETHNICITY_V4_20261028", canonical = text, sha256 = digest::digest(text, algo = "sha256", 
                serialize = FALSE))
        }
        ff_ce_fail <- function(code, receipt = NULL) stop(structure(list(message = code, call = NULL, numerical_audit = receipt), 
            class = c("ff_ce_error", "error", "condition")))
        ff_ce_bad_warning <- function(w) any(grepl("infinite|converg|singular|iteration", w, ignore.case = TRUE))
        ff_ce_centers <- function(X, reference = NULL, centre_actual = NULL) {
            columns <- grep("^centerc[0-9]+$", colnames(X), value = TRUE)
            ff_y_assert(length(columns) > 0L && !length(setdiff(grep("^center", colnames(X), value = TRUE), columns)), 
                "CE_FROZEN_CENTER_COLUMN_MAPPING_NOT_EXACT")
            if (!is.null(reference)) ff_y_assert(identical(columns, reference$columns), "CE_CENTER_COLUMN_DRIFT")
            z <- X[, columns, drop = FALSE]
            ff_y_assert(all(z %in% 0:1) && all(rowSums(z) <= 1), "CE_CENTER_DUMMIES_NOT_MUTUALLY_EXCLUSIVE")
            labels <- c("center_REFERENCE_FROZEN", columns)
            index <- ifelse(rowSums(z) == 0L, 1L, max.col(z, ties.method = "first") + 1L)
            group <- labels[index]
            if (!is.null(centre_actual)) {
                ff_y_assert(length(centre_actual) == nrow(X) && !anyNA(centre_actual), "CE_ACTUAL_CENTER_BINDING_INVALID")
                actual <- as.character(centre_actual)
                ff_y_assert(all(vapply(split(actual, group), function(v) length(unique(v)) == 1L, logical(1L))) && 
                  all(vapply(split(group, actual), function(v) length(unique(v)) == 1L, logical(1L))), "CE_ACTUAL_CENTER_GROUPING_DISAGREES_WITH_FROZEN_X")
            }
            if (is.null(reference)) reference <- list(columns = columns, all_design_levels = labels, observed_levels = labels[labels %in% 
                group], source = "EXACT_FROZEN_CENTER_TREATMENT_DUMMIES_NO_ROW_DELETION")
            ff_y_assert(all(group %in% reference$observed_levels), "CE_PREDICTION_CENTER_NOT_OBSERVED_IN_FIT")
            list(group = group, reference = reference, actual_center_grouping_checked = !is.null(centre_actual))
        }
        ff_ce_columns <- function(X, center_columns) {
            asian <- "ethnicityAsian or Asian British"
            chinese <- "ethnicityChinese"
            ff_y_assert(all(c(asian, chinese) %in% colnames(X)), "CE_EXPECTED_FROZEN_ETHNICITY_COLUMNS_ABSENT")
            ff_y_assert(all(X[, asian] %in% 0:1) && all(X[, chinese] %in% 0:1) && all(X[, asian] + X[, chinese] <= 
                1), "CE_ETHNICITY_DUMMIES_NOT_DISJOINT")
            z <- X[, !colnames(X) %in% c(center_columns, chinese), drop = FALSE]
            z[, asian] <- X[, asian] + X[, chinese]
            colnames(z)[colnames(z) == asian] <- "ethnicityAsian_or_Chinese_pooled"
            z
        }
        ff_ce_support <- function(ep, centers, tau) do.call(rbind, lapply(centers$reference$observed_levels, function(g) {
            keep <- centers$group == g
            data.frame(center_design_label = g, N = sum(keep), disease_events = sum(ep$status[keep] == 1L), death_events = sum(ep$status[keep] == 
                2L), max_followup = max(ep$time[keep]), at_risk_at_tau = sum(ep$time[keep] >= tau), supported_at_tau = any(ep$time[keep] >= 
                tau), stringsAsFactors = FALSE)
        }))
        ff_ce_native <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL) {
            if (is.null(init)) init <- rep(0, length(safe_names))
            strata <- survival::strata
            f <- stats::reformulate(c(safe_names, if (stratified) "strata(.ce_center)"), response = if (cause == 
                1L) "survival::Surv(.ff_time,.ff_disease)" else "survival::Surv(.ff_time,.ff_death)", env = environment())
            if (cluster_bound) ff_y_capture(survival::coxph(f, data = work, weights = .ff_weight, cluster = .ff_cluster, 
                robust = TRUE, ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, x = TRUE, 
                y = TRUE, control = control, init = init)) else ff_y_capture(survival::coxph(f, data = work, weights = .ff_weight, 
                robust = FALSE, ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, 
                x = TRUE, y = TRUE, control = control, init = init))
        }
        ff_ce_score <- function(fit) {
            detail <- survival::coxph.detail(fit, riskmat = FALSE)
            p <- length(stats::coef(fit))
            u <- if (p == 1L) sum(detail$score) else colSums(detail$score)
            I <- if (p == 1L) matrix(sum(detail$imat), 1L) else apply(detail$imat, c(1L, 2L), sum)
            V <- if (is.null(fit$naive.var)) stats::vcov(fit) else fit$naive.var
            step <- as.numeric(u %*% V)
            eig <- eigen((I + t(I))/2, symmetric = TRUE, only.values = TRUE)$values
            list(score_max_abs = max(abs(u)), newton_step_max_abs = max(abs(step)), newton_step_over_SE_max_abs = max(abs(step/sqrt(diag(V)))), 
                min_information_eigenvalue = min(eig), information_condition = max(eig)/min(eig), inverse_information_vs_naive_V_max_abs = max(abs(solve(I) - 
                  V)), public_recomputed_not_bitwise_native_terminal_score = TRUE)
        }
        ff_ce_baselines <- function(fit, safe_names, stratified, center_levels, support, cause) {
            if (!stratified) return(list(GLOBAL = ff_y_native_baseline(fit, safe_names)))
            bh <- ff_y_capture(survival::basehaz(fit, centered = FALSE))
            ff_y_assert(!ff_ce_bad_warning(bh$warnings) && "strata" %in% names(bh$value), "CE_STRATIFIED_BASEHAZ_FAILED")
            labels <- as.character(bh$value$strata)
            ff_y_assert(setequal(unique(labels), center_levels), "CE_NATIVE_BASELINE_STRATA_LABEL_BINDING_FAILED")
            result <- lapply(center_levels, function(g) {
                b <- bh$value[labels == g, , drop = FALSE]
                zero <- as.data.frame(stats::setNames(as.list(rep(0, length(safe_names))), safe_names))
                zero$.ce_center <- factor(g, levels = center_levels)
                sf <- ff_y_capture(survival::survfit(fit, newdata = zero, se.fit = FALSE, conf.type = "none", stype = 2, 
                  ctype = 2))
                difference <- if (identical(as.numeric(b$time), as.numeric(sf$value$time))) max(abs(b$hazard - 
                  as.numeric(sf$value$cumhaz)), 0) else Inf
                ff_y_assert(!ff_ce_bad_warning(sf$warnings) && is.finite(difference) && difference <= 1e-10 && 
                  all(is.finite(b$hazard)) && all(diff(c(0, b$hazard)) >= 0), "CE_NATIVE_BASEHAZ_SURVFIT_STRATUM_VALIDATION_FAILED", 
                  g)
                s <- support[support$center_design_label == g, , drop = FALSE]
                events <- if (cause == 1L) s$disease_events else s$death_events
                if (events == 0L) ff_y_assert(all(b$hazard == 0) && s$supported_at_tau, "CE_ZERO_EVENT_STRATUM_NOT_NATIVE_ZERO_OR_SUPPORTED")
                list(time = b$time, cumulative = b$hazard, max_abs_basehaz_vs_survfit = difference, tolerance = 1e-10, 
                  warnings = unique(c(bh$warnings, sf$warnings)), zero_event_supported = events == 0L, support = list(N = s$N, 
                    max_followup = s$max_followup, at_risk_at_tau = s$at_risk_at_tau, events = events))
            })
            stats::setNames(result, center_levels)
        }
        ff_ce_merge_baselines <- function(b) {
            grid <- sort(unique(unlist(lapply(b, `[[`, "time"))))
            ff_y_assert(length(grid) > 0L, "CE_NO_NATIVE_BASELINE_GRID")
            cumulative <- vapply(b, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
            cumulative <- matrix(cumulative, ncol = 2L)
            hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
                , drop = FALSE])
            ff_y_assert(all(is.finite(hazard)) && all(hazard >= 0), "CE_INVALID_NATIVE_HAZARD_INCREMENT")
            keep <- rowSums(hazard) > 0
            list(time = grid[keep], hazard = hazard[keep, , drop = FALSE])
        }
        ff_ce_fit <- function(X, risk, w, cluster = NULL, tau_report, ties = "efron", diagnostics = FALSE, endpoint_id, 
            centre_actual = NULL, audit_callback = NULL) {
            ff_y_assert(endpoint_id %in% c("250", "250.2", "411", "411.3", "411.4") && length(endpoint_id) == 1L, 
                "CE_FIXED_FIVE_ENDPOINT_SCOPE_REQUIRED")
            ff_y_assert(identical(as.numeric(tau_report), 11) && identical(ties, "efron"), "CE_FIXED_TAU11_AND_EFRON_REQUIRED")
            runtime <- ff_y_runtime(require = TRUE)
            original <- ff_y_matrix(X)
            ep <- ff_y_risk(risk, tau_report)
            ff_y_assert(nrow(original$X) == ep$n && is.numeric(w) && length(w) == ep$n && all(is.finite(w) & w > 
                0), "CE_ORIGINAL_WEIGHT_SAMPLE_INVALID")
            if (!is.null(cluster)) ff_y_assert(length(cluster) == ep$n && !anyNA(cluster), "CE_CLUSTER_BINDING_INVALID")
            centers <- ff_ce_centers(original$X, centre_actual = centre_actual)
            support <- ff_ce_support(ep, centers, tau_report)
            ff_y_assert(all(support$supported_at_tau), "CE_CENTER_TAU11_SUPPORT_ABSENT_NO_EXTRAPOLATION")
            changed <- ff_y_matrix(ff_ce_columns(original$X, centers$reference$columns))
            stratified <- c(disease = identical(endpoint_id, "411.3"), death = TRUE)
            designs <- list(disease = if (stratified[1L]) changed else original, death = changed)
            fits <- baselines <- audits <- vector("list", 2L)
            warnings <- character()
            lp_difference <- numeric(2L)
            names(fits) <- names(baselines) <- names(audits) <- names(designs)
            probe <- unique(as.integer(round(seq(1, ep$n, length.out = min(16L, ep$n)))))
            ctrl0 <- survival::coxph.control(iter.max = 50L)
            ctrl1 <- survival::coxph.control(eps = 1e-11, iter.max = 50L, toler.inf = ctrl0$toler.inf)
            emit <- function(rec) if (is.function(audit_callback)) audit_callback(rec)
            for (j in 1:2) {
                dd <- designs[[j]]
                work <- dd$work
                work$.ce_center <- factor(centers$group, levels = centers$reference$observed_levels)
                work$.ff_time <- ep$time
                work$.ff_disease <- as.integer(ep$status == 1L)
                work$.ff_death <- as.integer(ep$status == 2L)
                work$.ff_weight <- w
                if (!is.null(cluster)) work$.ff_cluster <- as.character(cluster)
                binding_hash <- digest::digest(list(X = dd$X, time = ep$time, status = ep$status, weights = w, 
                  center = centers$group, cluster = cluster), algo = "sha256")
                got <- ff_ce_native(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), ctrl0)
                m <- got$value
                rec <- list(cause = names(designs)[j], binding_hash = binding_hash, original_control = ctrl0, original_warnings = got$warnings, 
                  original_coefficients = stats::coef(m), original_vcov = stats::vcov(m), original_iterations = m$iter, 
                  original_final_loglik = tail(m$loglik, 1L), retry_attempted = FALSE)
                emit(c(rec, list(stage = "ORIGINAL_ATTEMPT_PRESERVED")))
                finite <- all(is.finite(stats::coef(m))) && all(is.finite(stats::vcov(m)))
                if (ff_ce_bad_warning(got$warnings)) {
                  eligible <- length(got$warnings) > 0L && finite && all(grepl("Loglik converged before variable.*coefficient may be infinite", 
                    got$warnings))
                  if (!eligible) ff_ce_fail("CE_UNRESOLVED_ORIGINAL_NATIVE_WARNING_NO_ELIGIBLE_RETRY", rec)
                  old <- m
                  rec$retry_attempted <- TRUE
                  rec$retry_control <- ctrl1
                  rec$original_score_information <- ff_ce_score(old)
                  got <- ff_ce_native(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), ctrl1, 
                    stats::coef(old))
                  m <- got$value
                  rec$retry_warnings <- got$warnings
                  rec$retry_iterations <- m$iter
                  rec$retry_final_loglik <- tail(m$loglik, 1L)
                  rec$warm_start_initial_statistics_are_not_null_statistics <- TRUE
                  rec$coefficient_max_abs_difference <- max(abs(stats::coef(m) - stats::coef(old)))
                  rec$coefficient_difference_over_old_SE_max_abs <- max(abs((stats::coef(m) - stats::coef(old))/sqrt(diag(stats::vcov(old)))))
                  rec$SE_max_abs_difference <- max(abs(sqrt(diag(stats::vcov(m))) - sqrt(diag(stats::vcov(old)))))
                  rec$final_loglik_difference <- tail(m$loglik, 1L) - tail(old$loglik, 1L)
                  rec$LP_max_abs_difference <- max(abs(as.numeric(predict(m, newdata = work, type = "lp", reference = "zero")) - 
                    as.numeric(predict(old, newdata = work, type = "lp", reference = "zero"))))
                  rec$retry_score_information <- ff_ce_score(m)
                  emit(c(rec, list(stage = "SINGLE_NUMERICAL_RETRY_COMPLETED")))
                  if (ff_ce_bad_warning(got$warnings)) ff_ce_fail("CE_UNRESOLVED_NATIVE_WARNING_AFTER_SINGLE_RETRY", 
                    rec)
                }
                ff_y_assert(inherits(m, "coxph") && !inherits(m, "coxphms") && identical(m$method, "efron") && 
                  all(is.finite(stats::coef(m))) && all(is.finite(stats::vcov(m))), "CE_NATIVE_NONFINITE_OR_WRONG_ESTIMATOR")
                mf <- stats::model.frame(m)
                ff_y_assert(m$n == ep$n && identical(as.numeric(m$y[, "time"]), as.numeric(ep$time)) && identical(as.numeric(m$y[, 
                  "status"]), as.numeric(ep$status == j)) && identical(as.numeric(stats::model.weights(mf)), as.numeric(w)) && 
                  isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), "CE_FROZEN_TIME_EVENT_WEIGHT_OR_DESIGN_CHANGED")
                if (stratified[j]) ff_y_assert(identical(as.character(mf[["strata(.ce_center)"]]), centers$group), 
                  "CE_NATIVE_STRATUM_BINDING_CHANGED")
                if (!is.null(cluster)) ff_y_assert(!is.null(m$naive.var) && identical(as.character(stats::model.extract(mf, 
                  "cluster")), as.character(cluster)), "CE_ROBUST_CLUSTER_BINDING_CHANGED")
                baselines[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, stratified[j], centers$reference$observed_levels, 
                  support, j)
                if (isTRUE(rec$retry_attempted)) {
                  oldbase <- ff_ce_baselines(old, dd$reference$safe_names, stratified[j], centers$reference$observed_levels, 
                    support, j)
                  rec$baseline_max_abs_difference <- max(vapply(names(oldbase), function(g) {
                    ff_y_assert(identical(oldbase[[g]]$time, baselines[[j]][[g]]$time), "CE_RETRY_BASELINE_GRID_CHANGED")
                    max(abs(oldbase[[g]]$cumulative - baselines[[j]][[g]]$cumulative), 0)
                  }, numeric(1L)))
                }
                rec$binding_hash_after <- digest::digest(list(X = dd$X, time = ep$time, status = ep$status, weights = w, 
                  center = centers$group, cluster = cluster), algo = "sha256")
                ff_y_assert(identical(binding_hash, rec$binding_hash_after), "CE_RETRY_BINDING_MUTATION")
                emit(c(rec, list(stage = "ACCEPTED_NATIVE_ATTEMPT_NOT_FULL_PIPELINE_VALIDITY")))
                fits[[j]] <- m
                audits[[j]] <- rec
                warnings <- c(warnings, got$warnings)
                nd <- work[probe, , drop = FALSE]
                lp <- ff_y_capture(as.numeric(predict(m, newdata = nd, type = "lp", reference = "zero")))
                ff_y_assert(!ff_ce_bad_warning(lp$warnings), "CE_NATIVE_LP_WARNING")
                lp_difference[j] <- max(abs(lp$value - as.numeric(dd$X[probe, , drop = FALSE] %*% stats::coef(m))))
                ff_y_assert(lp_difference[j] <= 1e-10, "CE_NATIVE_LP_MATRIX_BINDING_FAILED")
            }
            engines <- stats::setNames(lapply(centers$reference$observed_levels, function(g) ff_ce_merge_baselines(list(if (stratified[1L]) baselines$disease[[g]] else baselines$disease$GLOBAL, 
                baselines$death[[g]]))), centers$reference$observed_levels)
            coefficient_table <- do.call(rbind, lapply(1:2, function(j) data.frame(cause = names(fits)[j], term = designs[[j]]$reference$columns, 
                internal_term = designs[[j]]$reference$safe_names, coefficient = unname(stats::coef(fits[[j]])), 
                SE = sqrt(diag(stats::vcov(fits[[j]]))), row.names = NULL)))
            all_baselines <- unlist(baselines, recursive = FALSE)
            out <- list(fit = fits, engine = list(by_center = engines, baseline_public_validation = baselines, 
                provenance = ff_y_kernel_provenance()), matrix_reference = original$reference, cause_references = lapply(designs, 
                `[[`, "reference"), center_reference = centers$reference, center_support = support, stratum_support = support, 
                centre_actual_grouping_checked = centers$actual_center_grouping_checked, stratified = stratified, 
                coefficients = coefficient_table, endpoint = ep, tau_report = tau_report, endpoint_id = endpoint_id, 
                ties = ties, runtime = runtime, binding = list(weight = w, cluster = cluster), numerical_audit = audits, 
                numerical_retry = audits, method = ff_ce_method_spec()$id, method_spec = ff_ce_method_spec(), method_spec_hash = ff_ce_method_spec()$sha256, 
                warnings = unique(c(warnings, unlist(lapply(all_baselines, `[[`, "warnings")))), weight = list(min = min(w), 
                  max = max(w), sum = sum(w), ESS = sum(w)^2/sum(w^2), interpretation = "UNCHANGED_ORIGINAL_SELECTION_WEIGHTS"), 
                representative_validation_required = TRUE, formal_ready = FALSE, native_validation = list(passed = FALSE, 
                  lp_max_abs = max(lp_difference), basehaz_survfit_max_abs = max(vapply(all_baselines, `[[`, numeric(1L), 
                    "max_abs_basehaz_vs_survfit"))))
            class(out) <- c("ff_ce_matrix_fit", "ff_y_matrix_fit", "ff_y_fit")
            if (isTRUE(diagnostics)) out$PH <- lapply(1:2, function(j) ff_y_ph(fits[[j]], designs[[j]]$reference$columns, 
                designs[[j]]$reference$safe_names))
            curve <- lapply(c(0, 0.25, 0.5, 0.75, 1) * tau_report, function(t) ff_ce_predict_at(out, original$X[probe, 
                , drop = FALSE], t))
            violation <- 0
            for (k in 2:length(curve)) {
                delta <- curve[[k]] - curve[[k - 1L]]
                violation <- max(violation, -delta[, "disease"], -delta[, "death"], delta[, "survival"])
            }
            ff_y_assert(violation <= 1e-12, "CE_CIF_MONOTONICITY_FAILED")
            out$native_validation$monotonicity_violation <- violation
            out$native_validation$three_state_closure <- max(vapply(curve, function(z) max(abs(rowSums(z) - 1)), 
                numeric(1L)))
            out$native_validation$passed <- TRUE
            out
        }
        ff_ce_predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
            64)) {
            ff_y_assert(inherits(object, "ff_ce_matrix_fit"), "CE_FITTED_OBJECT_REQUIRED")
            ff_y_assert(length(tau_report) == 1L && is.finite(tau_report) && tau_report >= 0 && tau_report <= object$tau_report, 
                "CE_TIME_OUTSIDE_FIXED_TAU11_SUPPORT")
            ff_y_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1L && length(max_prediction_mb) == 
                1L && is.finite(max_prediction_mb) && max_prediction_mb > 0, "CE_PREDICTION_MEMORY_BUDGET_INVALID")
            original <- ff_y_matrix(X, object$matrix_reference)
            centers <- ff_ce_centers(original$X, object$center_reference)
            changed <- ff_ce_columns(original$X, object$center_reference$columns)
            designs <- list(disease = ff_y_matrix(if (object$stratified[1L]) changed else original$X, object$cause_references$disease), 
                death = ff_y_matrix(changed, object$cause_references$death))
            step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
            p <- matrix(NA_real_, nrow(original$X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
            warnings <- character()
            for (g in unique(centers$group)) {
                support <- object$center_support[object$center_support$center_design_label == g, , drop = FALSE]
                ff_y_assert(nrow(support) == 1L && support$max_followup >= tau_report, "CE_PREDICTION_STRATUM_UNSUPPORTED_NO_EXTRAPOLATION")
                rows <- which(centers$group == g)
                eng <- object$engine$by_center[[g]]
                h <- eng$hazard[eng$time <= tau_report, , drop = FALSE]
                for (start in seq.int(1L, length(rows), by = step)) {
                  ii <- rows[seq.int(start, min(length(rows), start + step - 1L))]
                  lp <- vapply(1:2, function(j) {
                    nd <- designs[[j]]$work[ii, , drop = FALSE]
                    nd$.ce_center <- factor(g, levels = object$center_reference$observed_levels)
                    got <- ff_y_capture(as.numeric(predict(object$fit[[j]], newdata = nd, type = "lp", reference = "zero")))
                    warnings <<- c(warnings, got$warnings)
                    got$value
                  }, numeric(length(ii)))
                  rr <- exp(matrix(lp, nrow = length(ii), ncol = 2L))
                  ff_y_assert(!ff_ce_bad_warning(warnings) && all(is.finite(rr)) && all(rr > 0), "CE_INVALID_NATIVE_RR_OR_WARNING_NO_CLIPPING")
                  if (exists("ff_aj_status", mode = "function") && isTRUE(ff_aj_status()$ready)) p[ii, ] <- ff_aj_backend(h, 
                    rr) else {
                    ff_y_assert(!isTRUE(getOption("ff.y.require_compiled_aj", FALSE)), "FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY")
                    p[ii, ] <- ff_y_aj(h, cbind(rr, rr), rep(1L, nrow(h)))
                  }
                  ff_y_probability_check(p[ii, , drop = FALSE])
                }
            }
            ff_y_probability_check(p)
            attr(p, "warnings") <- unique(warnings)
            attr(p, "chunk_size") <- step
            attr(p, "product.limit") <- TRUE
            attr(p, "risk_range") <- apply(p, 2L, range)
            p
        }
        ff_ce_predict_y_matrix <- function(object, X, tau_report = object$tau_report, chunk_size = getOption("ff.y.chunk_size", 
            4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 64)) {
            ff_y_assert(identical(as.numeric(tau_report), as.numeric(object$tau_report)), "CE_REPORT_TIME_DRIFT")
            ff_ce_predict_at(object, X, tau_report, chunk_size, max_prediction_mb)
        }
        ff_ce_install <- function(endpoint_id, envir = parent.frame(), audit_callback = NULL) {
            force(endpoint_id)
            force(envir)
            force(audit_callback)
            ff_y_assert(endpoint_id %in% c("250", "250.2", "411", "411.3", "411.4"), "CE_FIXED_ENDPOINT_REQUIRED")
            if (!exists(".ff_ce_original_predict_at", envir = envir, inherits = FALSE)) {
                assign(".ff_ce_original_predict_at", get("ff_y_predict_at", envir = envir, inherits = TRUE), envir = envir)
                assign(".ff_ce_original_fit_y_matrix", get("ff_fit_y_matrix", envir = envir, inherits = TRUE), 
                  envir = envir)
            }
            previous_predict <- get(".ff_ce_original_predict_at", envir = envir, inherits = FALSE)
            fit <- function(X, risk, w, cluster = NULL, tau_report, ties = "efron", diagnostics = FALSE, centre_actual = NULL) ff_ce_fit(X, 
                risk, w, cluster, tau_report, ties, diagnostics, endpoint_id, centre_actual, audit_callback)
            predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
                64)) {
                if (inherits(object, "ff_ce_matrix_fit")) ff_ce_predict_at(object, X, tau_report, chunk_size, max_prediction_mb) else previous_predict(object, 
                  X, tau_report, chunk_size, max_prediction_mb)
            }
            assign("ff_fit_y_matrix", fit, envir = envir)
            assign("ff_y_predict_at", predict_at, envir = envir)
            invisible(list(endpoint_id = endpoint_id, method_spec = ff_ce_method_spec(), formal_ready = FALSE, 
                installed_functions = c("ff_fit_y_matrix", "ff_y_predict_at")))
        }
    }), envir = context)
}

