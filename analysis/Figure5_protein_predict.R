# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

pry_predict_assert <- function(ok, message) {
    if (!isTRUE(ok)) 
        stop(paste0("PRY_PREDICT_", message), call. = FALSE)
    invisible(TRUE)
}

pry_predict_tolerances <- function() list(schema = "SELECTED44_PROTEIN_RIDGE_PREDICTION_TOLERANCES_V1", fixed_before_comparisons = TRUE, 
    native_lp_abs = 1e-10, cumhaz_abs = 1e-12, cumhaz_rel = 1e-09, probability_abs = 1e-12, probability_rel = 1e-09, 
    scalar_R_probability_abs = 1e-12, probability_mass_abs = 1e-10, time_monotonicity_abs = 1e-12, basis = paste("native LP inherited fd_attach_fast_prediction lp_tolerance=1e-10;", 
        "scalar R probability inherited cif_tolerance=1e-12; probability/cumhaz", "abs1e-12+rel1e-9 fixed by repair specification; mass1e-10 inherited", 
        "c39_num_spec; no clipping, renormalization, NA-to-zero, or adaptive relaxation"))

pry_predict_capture <- function(expr, label) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (length(warnings)) 
        stop(structure(list(message = paste0("PRY_PREDICT_", label, ": ", paste(unique(warnings), collapse = " | ")), 
            call = NULL, warnings = unique(warnings)), class = c("pry_prediction_error", "error", "condition")))
    value
}

pry_predict_close <- function(actual, expected, atol, rtol, label) {
    pry_predict_assert(identical(dim(actual), dim(expected)) && length(actual) == length(expected) && all(is.finite(actual)) && 
        all(is.finite(expected)), paste0(label, "_SHAPE_OR_FINITE"))
    error <- abs(actual - expected)
    pry_predict_assert(all(error <= atol + rtol * pmax(abs(actual), abs(expected))), label)
    max(error, 0)
}

pry_predict_check_kernel <- function(engine) {
    pry_predict_assert(is.environment(engine) && is.function(engine$ff_aj_status) && is.function(engine$ff_aj_backend), 
        "ORIGINAL_COMPILED_KERNEL_REQUIRED")
    spec <- engine$ff_aj_status()
    pry_predict_assert(isTRUE(spec$ready) && identical(spec$computation_precision, "DOUBLE_ONLY_BINARY64") && identical(spec$method, 
        "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL") && identical(spec$clipping, FALSE) && 
        identical(spec$renormalization, FALSE), "EXACT_FROZEN_DOUBLE_EXPONENTIAL_KERNEL_REQUIRED")
    spec
}

pry_predict_kernel <- function(engine, hazard, relative, tolerances) {
    pry_predict_assert(is.matrix(hazard) && ncol(hazard) == 2L && all(is.finite(hazard) & hazard >= 0) && is.matrix(relative) && 
        nrow(relative) >= 1L && ncol(relative) == 2L && all(is.finite(relative) & relative > 0), "INVALID_HAZARD_OR_RELATIVE_NO_CLIPPING")
    p <- engine$ff_aj_backend(hazard, relative)
    pry_predict_assert(is.matrix(p) && nrow(p) == nrow(relative) && ncol(p) == 3L && all(is.finite(p) & p >= 0 & 
        p <= 1) && all(abs(rowSums(p) - 1) <= tolerances$probability_mass_abs), "INVALID_THREE_STATE_PROBABILITY")
    colnames(p) <- c("disease", "death", "survival")
    p
}

pry_predict_piece_check <- function(piece, tau) {
    fit <- piece$fit
    pry_predict_assert(inherits(fit, "coxph.penal") && inherits(fit, "coxph") && !inherits(fit, "coxphms") && identical(fit$method, 
        "efron") && is.matrix(fit$x) && ncol(fit$x) == length(stats::coef(fit)) && all(is.finite(stats::coef(fit))) && 
        all(is.finite(stats::vcov(fit))) && is.null(fit$naive.var) && !length(fit$offset) && is.null(fit$frail) && 
        !length(fit$na.action), "FINITE_NATIVE_PENAL_EFRON_NONROBUST_FIT_REQUIRED")
    pry_predict_assert(isTRUE(piece$diagnostics$passed), "FITTING_DIAGNOSTICS_MUST_PASS_FIRST")
    labels <- attr(fit$terms, "term.labels")
    pry_predict_assert(any(grepl("strata\\(\\.pry_center\\)", labels)) && all(attr(fit$terms, "order") == 1L) && 
        !length(attr(fit$terms, "offset")), "BOTH_CAUSES_EXPLICIT_CENTER_STRATA_NO_INTERACTIONS_OR_OFFSETS")
    support <- piece$center_support
    required <- c("center", "N", "max_followup", "at_risk_at_tau", "disease_events", "death_events", "supported_at_tau")
    pry_predict_assert(is.data.frame(support) && all(required %in% names(support)) && !anyDuplicated(as.character(support$center)), 
        "CENTER_SUPPORT_TABLE_REQUIRED")
    observed <- as.character(piece$observed_center_levels)
    pry_predict_assert(length(observed) > 0L && !anyNA(observed) && !anyDuplicated(observed) && setequal(observed, 
        as.character(support$center[support$N > 0])) && all(observed %in% as.character(piece$center_levels)), "ACTUAL_VS_FROZEN_CENTER_LEVEL_BINDING")
    s <- support[match(observed, as.character(support$center)), , drop = FALSE]
    pry_predict_assert(all(is.finite(s$max_followup) & s$max_followup >= tau) && all(is.finite(s$at_risk_at_tau) & 
        s$at_risk_at_tau > 0) && all(s$supported_at_tau), "CENTER_LACKS_ACTUAL_TAU11_SUPPORT_NO_TAIL_EXTRAPOLATION")
    mf <- stats::model.frame(fit)
    native_center <- as.character(mf[["strata(.pry_center)"]])
    weight <- stats::model.weights(mf)
    pry_predict_assert(is.matrix(fit$y) && ncol(fit$y) == 2L && all(c("time", "status") %in% colnames(fit$y)) && 
        length(native_center) == nrow(fit$y) && setequal(native_center, observed) && length(weight) == nrow(fit$y) && 
        all(is.finite(weight) & weight > 0), "ACTUAL_NATIVE_RIGHT_CENSORED_SAMPLE_WEIGHT_CENTER_BINDING")
    for (g in observed) {
        ii <- which(native_center == g)
        row <- support[as.character(support$center) == g, , drop = FALSE]
        pry_predict_assert(row$N == length(ii) && row$max_followup == max(fit$y[ii, "time"]) && row$at_risk_at_tau == 
            sum(fit$y[ii, "time"] >= tau) && row[[paste0(piece$cause, "_events")]] == sum(fit$y[ii, "status"] == 
            1L), "SUPPORT_AND_ZERO_EVENT_COUNTS_MUST_MATCH_NATIVE_FITTED_ROWS")
    }
    pry_predict_assert(is.matrix(piece$reference_X) && is.numeric(piece$reference_X) && all(is.finite(piece$reference_X)) && 
        nrow(piece$reference_X) == length(piece$reference_center) && !anyDuplicated(as.character(piece$reference_center)) && 
        setequal(as.character(piece$reference_center), observed), "ONE_ACTUAL_FINITE_REFERENCE_ROW_PER_OBSERVED_CENTER_REQUIRED")
    invisible(observed)
}

pry_predict_design <- function(piece, X, center, transform_X) {
    pry_predict_assert(is.matrix(X) && is.numeric(X) && nrow(X) > 0L && all(is.finite(X)) && length(center) == 
        nrow(X) && !anyNA(center), "FINITE_ORIGINAL_DESIGN_AND_ROW_CENTERS_REQUIRED")
    center <- as.character(center)
    pry_predict_assert(all(center %in% as.character(piece$observed_center_levels)), "ABSENT_OR_UNKNOWN_CENTER_IS_NOT_ZERO_EVENT_CENTER")
    nd <- transform_X(X, piece$transform, center)
    pry_predict_assert(is.data.frame(nd) && nrow(nd) == nrow(X) && identical(as.character(nd$.pry_center), center), 
        "SAVED_TRANSFORM_ROW_OR_CENTER_DRIFT")
    rownames(nd) <- sprintf("pry_row_%08d", seq_len(nrow(nd)))
    mm <- pry_predict_capture(stats::model.matrix(piece$fit, data = nd), "NATIVE_MODEL_MATRIX_WARNING")
    pry_predict_assert(is.matrix(mm) && nrow(mm) == nrow(X) && identical(colnames(mm), colnames(piece$fit$x)) && 
        ncol(mm) == length(stats::coef(piece$fit)) && all(is.finite(mm)), "NATIVE_PENAL_DESIGN_COLUMN_POSITION_BINDING_FAILED")
    lp <- as.numeric(mm %*% as.numeric(stats::coef(piece$fit)))
    pry_predict_assert(all(is.finite(lp)), "NONFINITE_TRANSFORMED_LINEAR_PREDICTOR")
    list(newdata = nd, matrix = mm, lp = lp)
}

pry_predict_native_lp <- function(piece, nd) {
    lp <- pry_predict_capture(as.numeric(stats::predict(piece$fit, newdata = nd, type = "lp", reference = "zero", 
        se.fit = FALSE)), "NATIVE_PENAL_LP_WARNING")
    pry_predict_assert(length(lp) == nrow(nd) && all(is.finite(lp)), "NATIVE_PENAL_LP_SHAPE")
    lp
}

pry_predict_native_curves <- function(piece, nd) {
    sf <- pry_predict_capture(survival::survfit(piece$fit, newdata = nd, se.fit = FALSE, conf.type = "none", stype = 2L, 
        ctype = 2L, censor = TRUE), "NATIVE_PENAL_SURVFIT_WARNING")
    time <- sf$time
    H <- sf$cumhaz
    n <- nrow(nd)
    pry_predict_assert(is.numeric(time) && length(time) > 0L && all(is.finite(time) & time >= 0) && !is.null(H) && 
        all(is.finite(H)), "MISSING_NATIVE_CUMHAZ_NOT_ZERO_FILLED")
    curves <- vector("list", n)
    if (is.matrix(H)) {
        pry_predict_assert(nrow(H) == length(time) && ncol(H) == n, "NATIVE_MATRIX_CURVE_SHAPE")
        for (i in seq_len(n)) curves[[i]] <- list(time = as.numeric(time), cumulative = as.numeric(H[, i]))
    }
    else {
        counts <- sf$strata
        if (is.null(counts) && n == 1L) 
            counts <- length(time)
        pry_predict_assert(is.numeric(H) && length(H) == length(time) && length(counts) == n && all(is.finite(counts) & 
            counts > 0 & counts == as.integer(counts)) && sum(counts) == length(time), "NATIVE_PACKED_CURVE_ROW_LENGTHS")
        if (!is.null(names(counts))) 
            pry_predict_assert(identical(names(counts), rownames(nd)), "NATIVE_PACKED_CURVE_ROW_ORDER")
        ends <- cumsum(counts)
        starts <- c(1L, head(ends, -1L) + 1L)
        for (i in seq_len(n)) {
            ii <- seq.int(starts[i], ends[i])
            curves[[i]] <- list(time = as.numeric(time[ii]), cumulative = as.numeric(H[ii]))
        }
    }
    for (curve in curves) pry_predict_assert(!is.unsorted(curve$time, strictly = TRUE) && length(curve$time) == 
        length(curve$cumulative) && all(curve$cumulative >= 0) && all(diff(curve$cumulative) >= 0), "NATIVE_CUMHAZ_NONNEGATIVE_NONDECREASING_REQUIRED")
    curves
}

pry_predict_merge_curves <- function(curves, tau) {
    pry_predict_assert(is.list(curves) && length(curves) == 2L, "TWO_CAUSE_CURVES_REQUIRED")
    grid <- sort(unique(unlist(lapply(curves, function(z) z$time[z$time <= tau]), use.names = FALSE)))
    if (!length(grid)) 
        return(list(time = numeric(), hazard = matrix(numeric(), 0L, 2L)))
    cumulative <- vapply(curves, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
    cumulative <- matrix(cumulative, ncol = 2L)
    hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
        , drop = FALSE])
    pry_predict_assert(all(is.finite(hazard) & hazard >= 0), "INVALID_NATIVE_HAZARD_INCREMENT_NO_CLIPPING")
    keep <- rowSums(hazard) > 0
    list(time = grid[keep], hazard = hazard[keep, , drop = FALSE])
}

pry_compile_endpoint <- function(pieces, engine, tau = 11, transform_X = pry_transform_X, tolerances = pry_predict_tolerances(), 
    audit_callback = NULL) {
    pry_predict_assert(is.list(pieces) && identical(names(pieces), c("disease", "death")) && identical(as.numeric(tau), 
        11) && is.function(transform_X), "PAIRED_CAUSES_FIXED_TAU11_REQUIRED")
    pry_predict_assert(identical(tolerances, pry_predict_tolerances()), "PREDECLARED_TOLERANCES_MAY_NOT_DRIFT")
    tolerance_prerecord <- list(stage = "PREDICTION_TOLERANCES_FIXED_BEFORE_COMPARISON", tolerances = tolerances, 
        native_interface = "survfit(coxph.penal,newdata=explicit,se.fit=FALSE,stype=2,ctype=2)")
    if (!is.null(audit_callback)) {
        pry_predict_assert(is.function(audit_callback), "PREDICTION_AUDIT_CALLBACK_MUST_BE_FUNCTION")
        audit_callback(tolerance_prerecord)
    }
    kernel <- pry_predict_check_kernel(engine)
    observed <- lapply(pieces, pry_predict_piece_check, tau = tau)
    pry_predict_assert(identical(observed[[1L]], observed[[2L]]) && identical(pieces$disease$model_id, pieces$death$model_id) && 
        identical(pieces$disease$b, pieces$death$b) && identical(pieces$disease$endpoint, pieces$death$endpoint), 
        "PAIRED_CAUSE_SAMPLE_MODEL_ENDPOINT_BINDING")
    frames <- lapply(pieces, function(p) stats::model.frame(p$fit))
    pry_predict_assert(identical(as.numeric(pieces$disease$fit$y[, "time"]), as.numeric(pieces$death$fit$y[, "time"])) && 
        identical(as.numeric(stats::model.weights(frames$disease)), as.numeric(stats::model.weights(frames$death))) && 
        identical(as.character(frames$disease[["strata(.pry_center)"]]), as.character(frames$death[["strata(.pry_center)"]])), 
        "TWO_CAUSES_MUST_SHARE_EXACT_ROWS_WEIGHTS_TIMES_AND_CENTER_ORDER")
    levels <- observed[[1L]]
    anchors <- baselines <- vector("list", 2L)
    names(anchors) <- names(baselines) <- names(pieces)
    classes <- lapply(pieces, function(z) class(z$fit))
    lp_error <- 0
    for (cause in names(pieces)) {
        piece <- pieces[[cause]]
        anchors[[cause]] <- stats::setNames(numeric(length(levels)), levels)
        baselines[[cause]] <- stats::setNames(vector("list", length(levels)), levels)
        for (g in levels) {
            row <- match(g, as.character(piece$reference_center))
            prepared <- pry_predict_design(piece, piece$reference_X[row, , drop = FALSE], g, transform_X)
            native <- pry_predict_native_lp(piece, prepared$newdata)
            lp_error <- max(lp_error, pry_predict_close(prepared$lp, native, tolerances$native_lp_abs, 0, "REFERENCE_NATIVE_LP_EQUIVALENCE"))
            curve <- pry_predict_native_curves(piece, prepared$newdata)[[1L]]
            support <- piece$center_support[as.character(piece$center_support$center) == g, , drop = FALSE]
            events <- support[[paste0(cause, "_events")]]
            pry_predict_assert(length(events) == 1L && is.finite(events) && events >= 0, "ACTUAL_CENTER_CAUSE_EVENT_COUNT_REQUIRED")
            if (events == 0) 
                pry_predict_assert(all(curve$cumulative == 0), "ZERO_EVENT_CENTER_MUST_BE_NATIVE_ZERO_NOT_MISSING")
            anchors[[cause]][[g]] <- prepared$lp
            baselines[[cause]][[g]] <- curve
        }
    }
    hazards <- stats::setNames(lapply(levels, function(g) pry_predict_merge_curves(lapply(baselines, `[[`, g), 
        tau)), levels)
    state <- new.env(parent = emptyenv())
    state$passed <- FALSE
    state$receipt <- NULL
    structure(list(schema = "SELECTED44_PROTEIN_RIDGE_COMPILED_ENDPOINT_V1", pieces = pieces, endpoint = pieces$disease$endpoint, 
        model_id = pieces$disease$model_id, b = pieces$disease$b, tau = tau, center_levels = levels, anchors = anchors, 
        baselines = baselines, hazards = hazards, transform_X = transform_X, tolerances = tolerances, kernel_spec = kernel, 
        tolerance_prerecord = tolerance_prerecord, original_fit_classes = classes, reference_lp_max_abs = lp_error, 
        validation_state = state, compiled_from_explicit_finite_native_newdata = TRUE, native_stype = 2L, native_ctype = 2L, 
        no_clipping = TRUE), class = "pry_ridge_compiled_endpoint")
}

pry_predict_core <- function(compiled, X, center, engine, chunk_size = 512L, max_prediction_mb = 64, tau = compiled$tau) {
    pry_predict_assert(inherits(compiled, "pry_ridge_compiled_endpoint") && is.matrix(X) && nrow(X) > 0L && length(center) == 
        nrow(X) && !anyNA(center) && all(as.character(center) %in% compiled$center_levels), "COMPILED_FIT_DESIGN_AND_KNOWN_PRESENT_CENTER_REQUIRED")
    pry_predict_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1 && chunk_size == as.integer(chunk_size) && 
        is.finite(max_prediction_mb) && max_prediction_mb > 0 && length(tau) == 1L && is.finite(tau) && tau >= 
        0 && tau <= compiled$tau, "BOUNDED_PREDICTION_REQUEST_REQUIRED")
    pry_predict_check_kernel(engine)
    largest <- max(vapply(compiled$pieces, function(p) ncol(p$fit$x), integer(1L)))
    step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * (largest + 64L))))))
    probability <- matrix(NA_real_, nrow(X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
    center <- as.character(center)
    for (g in unique(center)) {
        rows <- which(center == g)
        eng <- compiled$hazards[[g]]
        h <- eng$hazard[eng$time <= tau, , drop = FALSE]
        for (start in seq.int(1L, length(rows), by = step)) {
            ii <- rows[seq.int(start, min(length(rows), start + step - 1L))]
            lp <- vapply(names(compiled$pieces), function(cause) {
                p <- pry_predict_design(compiled$pieces[[cause]], X[ii, , drop = FALSE], center[ii], compiled$transform_X)
                p$lp - compiled$anchors[[cause]][[g]]
            }, numeric(length(ii)))
            relative <- exp(matrix(lp, nrow = length(ii), ncol = 2L))
            probability[ii, ] <- pry_predict_kernel(engine, h, relative, compiled$tolerances)
        }
    }
    attr(probability, "chunk_size") <- step
    attr(probability, "integration_method") <- compiled$kernel_spec$method
    attr(probability, "no_clipping") <- TRUE
    probability
}

pry_predict_matrix <- function(compiled, X, center, engine, chunk_size = 512L, max_prediction_mb = 64) {
    pry_predict_assert(isTRUE(compiled$validation_state$passed), "NATIVE_AND_TWO_CAUSE_VALIDATION_REQUIRED_BEFORE_PRODUCTION")
    pry_predict_core(compiled, X, center, engine, chunk_size, max_prediction_mb)
}

pry_scalar_R_probability <- function(hazard) {
    pry_predict_assert(is.matrix(hazard) && ncol(hazard) == 2L && all(is.finite(hazard) & hazard >= 0), "R_REFERENCE_HAZARD_INVALID")
    survival <- 1
    cif <- compensation <- c(0, 0)
    for (i in seq_len(nrow(hazard))) {
        a <- hazard[i, 1L]
        b <- hazard[i, 2L]
        if (a == 0 && b == 0) 
            next
        total <- a + b
        if (is.finite(total)) {
            share <- c(a, b)/total
            mass <- -expm1(-total)
            next_survival <- survival * exp(-total)
        }
        else {
            scale <- max(a, b)
            share <- (c(a, b)/scale)/sum(c(a, b)/scale)
            mass <- 1
            next_survival <- 0
        }
        increment <- survival * mass * share
        for (j in 1:2) {
            y <- increment[j] - compensation[j]
            t <- cif[j] + y
            compensation[j] <- (t - cif[j]) - y
            cif[j] <- t
        }
        survival <- next_survival
    }
    c(disease = cif[1L], death = cif[2L], survival = survival)
}

pry_validation_probe <- function(X, center, ethnicity, extra_counterfactual = list()) {
    pry_predict_assert(is.matrix(X) && "B" %in% colnames(X) && nrow(X) == length(center) && length(center) == length(ethnicity) && 
        !anyNA(center) && !anyNA(ethnicity), "VALIDATION_SCOPE_METADATA_REQUIRED")
    cells <- paste(center, ethnicity, sep = "\r")
    rows <- sort(unique(c(match(unique(cells), cells), as.integer(round(seq(1L, nrow(X), length.out = min(16L, 
        nrow(X))))))))
    original <- X[rows, , drop = FALSE]
    zero <- one <- original
    zero[, "B"] <- 0
    one[, "B"] <- 1
    make <- function(xx, label, kind) list(X = xx, center = as.character(center[rows]), ethnicity = as.character(ethnicity[rows]), 
        label = label, kind = kind)
    c(list(observed = make(original, "FIXED_OBSERVED_ROWS", "observed"), B0 = make(zero, "FIXED_B0_COUNTERFACTUAL_ROWS", 
        "counterfactual"), B1 = make(one, "FIXED_B1_COUNTERFACTUAL_ROWS", "counterfactual")), extra_counterfactual)
}

pry_validate_prediction <- function(compiled, cases, engine, expected_ethnicity, native_rows_per_call = 8L, require_counterfactual = TRUE) {
    compiled$validation_state$passed <- FALSE
    pry_predict_assert(is.list(cases) && length(cases) > 0L && is.numeric(native_rows_per_call), "EXPLICIT_BOUNDED_VALIDATION_CASES_REQUIRED")
    pry_predict_assert(length(native_rows_per_call) == 1L && native_rows_per_call >= 1L && native_rows_per_call == 
        as.integer(native_rows_per_call) && native_rows_per_call <= 16L, "BOUNDED_NATIVE_CURVE_BATCH_MAX16")
    tolerance <- compiled$tolerances
    errors <- c(native_lp = compiled$reference_lp_max_abs, native_cumhaz = 0, compiled_vs_native_CIF = 0, C_vs_scalar_R = 0, 
        probability_mass = 0, time_monotonicity = 0)
    coverage <- list()
    rows_checked <- 0L
    kinds <- character()
    for (case in cases) {
        X <- case$X
        center <- as.character(case$center)
        ethnicity <- as.character(case$ethnicity)
        pry_predict_assert(is.matrix(X) && "B" %in% colnames(X) && nrow(X) == length(center) && length(ethnicity) == 
            nrow(X) && !anyNA(ethnicity) && case$kind %in% c("observed", "counterfactual"), "FIXED_VALIDATION_ROW_METADATA_REQUIRED")
        kinds <- c(kinds, case$kind)
        coverage[[length(coverage) + 1L]] <- data.frame(center = center, ethnicity = ethnicity, B = X[, "B"], stringsAsFactors = FALSE)
        fast <- pry_predict_core(compiled, X, center, engine, chunk_size = native_rows_per_call)
        errors["probability_mass"] <- max(errors["probability_mass"], abs(rowSums(fast) - 1))
        for (g in unique(center)) {
            group_rows <- which(center == g)
            for (start in seq.int(1L, length(group_rows), by = native_rows_per_call)) {
                ii <- group_rows[seq.int(start, min(length(group_rows), start + native_rows_per_call - 1L))]
                native <- vector("list", 2L)
                names(native) <- names(compiled$pieces)
                for (cause in names(compiled$pieces)) {
                  piece <- compiled$pieces[[cause]]
                  prepared <- pry_predict_design(piece, X[ii, , drop = FALSE], center[ii], compiled$transform_X)
                  native_lp <- pry_predict_native_lp(piece, prepared$newdata)
                  errors["native_lp"] <- max(errors["native_lp"], pry_predict_close(prepared$lp, native_lp, tolerance$native_lp_abs, 
                    0, "NATIVE_LP_POSITIONAL_EQUIVALENCE"))
                  native[[cause]] <- pry_predict_native_curves(piece, prepared$newdata)
                  for (k in seq_along(ii)) {
                    curve <- native[[cause]][[k]]
                    ref <- compiled$baselines[[cause]][[g]]
                    keep <- curve$time <= compiled$tau
                    expected <- c(0, ref$cumulative)[findInterval(curve$time[keep], ref$time) + 1L] * exp(prepared$lp[k] - 
                      compiled$anchors[[cause]][[g]])
                    errors["native_cumhaz"] <- max(errors["native_cumhaz"], pry_predict_close(curve$cumulative[keep], 
                      expected, tolerance$cumhaz_abs, tolerance$cumhaz_rel, "NATIVE_ANCHORED_CUMHAZ_EQUIVALENCE"))
                  }
                }
                for (k in seq_along(ii)) {
                  merged <- pry_predict_merge_curves(lapply(native, `[[`, k), compiled$tau)
                  native_probability <- pry_predict_kernel(engine, merged$hazard, matrix(c(1, 1), 1L), tolerance)
                  scalar <- matrix(pry_scalar_R_probability(merged$hazard), nrow = 1L)
                  actual <- matrix(as.numeric(fast[ii[k], ]), nrow = 1L)
                  errors["compiled_vs_native_CIF"] <- max(errors["compiled_vs_native_CIF"], pry_predict_close(actual, 
                    native_probability, tolerance$probability_abs, tolerance$probability_rel, "NATIVE_TWO_CAUSE_CIF_EQUIVALENCE"))
                  errors["C_vs_scalar_R"] <- max(errors["C_vs_scalar_R"], pry_predict_close(native_probability, 
                    scalar, tolerance$scalar_R_probability_abs, 0, "FROZEN_C_KERNEL_VS_SCALAR_R_REFERENCE"))
                }
            }
        }
        before <- pry_predict_core(compiled, X, center, engine, chunk_size = native_rows_per_call, tau = 0)
        for (time in c(0.25, 0.5, 0.75, 1) * compiled$tau) {
            after <- pry_predict_core(compiled, X, center, engine, chunk_size = native_rows_per_call, tau = time)
            delta <- after - before
            errors["time_monotonicity"] <- max(errors["time_monotonicity"], -delta[, "disease"], -delta[, "death"], 
                delta[, "survival"])
            before <- after
        }
        rows_checked <- rows_checked + nrow(X)
    }
    scope <- do.call(rbind, coverage)
    pry_predict_assert(setequal(unique(scope$center), compiled$center_levels) && length(expected_ethnicity) > 0L && 
        setequal(unique(scope$ethnicity), as.character(expected_ethnicity)), "EVERY_PRESENT_CENTER_AND_ETHNICITY_MUST_BE_VALIDATED")
    pry_predict_assert("observed" %in% kinds && (!require_counterfactual || "counterfactual" %in% kinds), "OBSERVED_AND_COUNTERFACTUAL_VALIDATION_REQUIRED")
    pairs <- unique(scope[, c("center", "ethnicity"), drop = FALSE])
    for (k in seq_len(nrow(pairs))) pry_predict_assert(all(c(0, 1) %in% scope$B[scope$center == pairs$center[k] & 
        scope$ethnicity == pairs$ethnicity[k]]), "B0_AND_B1_REQUIRED_FOR_EVERY_VALIDATED_CENTER_ETHNICITY_CELL")
    pry_predict_assert(errors["time_monotonicity"] <= tolerance$time_monotonicity_abs && errors["probability_mass"] <= 
        tolerance$probability_mass_abs, "TIME_MONOTONICITY_OR_MASS_FAILURE")
    pry_predict_assert(identical(lapply(compiled$pieces, function(z) class(z$fit)), compiled$original_fit_classes), 
        "NATIVE_COXPH_PENAL_CLASSES_MUST_REMAIN_UNCHANGED")
    receipt <- list(schema = "SELECTED44_PROTEIN_RIDGE_NATIVE_PREDICTION_VALIDATION_V1", passed = TRUE, rows_checked = rows_checked, 
        fixed_cases = lapply(cases, function(z) list(label = z$label, kind = z$kind, rows = nrow(z$X))), centers = compiled$center_levels, 
        ethnicities = as.character(expected_ethnicity), B_interventions = c(0, 1), max_abs = as.list(errors), tolerances = tolerance, 
        native_rows_per_call = native_rows_per_call, native_interface = "survfit(coxph.penal,newdata=explicit,se.fit=FALSE,stype=2,ctype=2)", 
        baseline_reference = "ACTUAL_FINITE_ROW_PER_CENTER;exp(LP_new-LP_anchor);same_saved_transform", compiled_kernel = compiled$kernel_spec$method, 
        original_classes = compiled$original_fit_classes, zero_event_supported_centers = lapply(compiled$pieces, 
            function(p) as.character(p$center_support$center[p$center_support$N > 0 & p$center_support[[paste0(p$cause, 
                "_events")]] == 0])), no_class_deletion = TRUE, no_clipping = TRUE, no_NA_fill = TRUE, fits_performed = 0L, 
        uses_RNG = FALSE)
    compiled$validation_state$receipt <- receipt
    compiled$validation_state$passed <- TRUE
    receipt
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        pry_predict_assert <- function(ok, message) {
            if (!isTRUE(ok)) stop(paste0("PRY_PREDICT_", message), call. = FALSE)
            invisible(TRUE)
        }
        pry_predict_tolerances <- function() list(schema = "SELECTED44_PROTEIN_RIDGE_PREDICTION_TOLERANCES_V1", 
            fixed_before_comparisons = TRUE, native_lp_abs = 1e-10, cumhaz_abs = 1e-12, cumhaz_rel = 1e-09, probability_abs = 1e-12, 
            probability_rel = 1e-09, scalar_R_probability_abs = 1e-12, probability_mass_abs = 1e-10, time_monotonicity_abs = 1e-12, 
            basis = paste("native LP inherited fd_attach_fast_prediction lp_tolerance=1e-10;", "scalar R probability inherited cif_tolerance=1e-12; probability/cumhaz", 
                "abs1e-12+rel1e-9 fixed by repair specification; mass1e-10 inherited", "c39_num_spec; no clipping, renormalization, NA-to-zero, or adaptive relaxation"))
        pry_predict_capture <- function(expr, label) {
            warnings <- character()
            value <- withCallingHandlers(expr, warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
            if (length(warnings)) stop(structure(list(message = paste0("PRY_PREDICT_", label, ": ", paste(unique(warnings), 
                collapse = " | ")), call = NULL, warnings = unique(warnings)), class = c("pry_prediction_error", 
                "error", "condition")))
            value
        }
        pry_predict_close <- function(actual, expected, atol, rtol, label) {
            pry_predict_assert(identical(dim(actual), dim(expected)) && length(actual) == length(expected) && all(is.finite(actual)) && 
                all(is.finite(expected)), paste0(label, "_SHAPE_OR_FINITE"))
            error <- abs(actual - expected)
            pry_predict_assert(all(error <= atol + rtol * pmax(abs(actual), abs(expected))), label)
            max(error, 0)
        }
        pry_predict_check_kernel <- function(engine) {
            pry_predict_assert(is.environment(engine) && is.function(engine$ff_aj_status) && is.function(engine$ff_aj_backend), 
                "ORIGINAL_COMPILED_KERNEL_REQUIRED")
            spec <- engine$ff_aj_status()
            pry_predict_assert(isTRUE(spec$ready) && identical(spec$computation_precision, "DOUBLE_ONLY_BINARY64") && 
                identical(spec$method, "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL") && identical(spec$clipping, 
                FALSE) && identical(spec$renormalization, FALSE), "EXACT_FROZEN_DOUBLE_EXPONENTIAL_KERNEL_REQUIRED")
            spec
        }
        pry_predict_kernel <- function(engine, hazard, relative, tolerances) {
            pry_predict_assert(is.matrix(hazard) && ncol(hazard) == 2L && all(is.finite(hazard) & hazard >= 0) && 
                is.matrix(relative) && nrow(relative) >= 1L && ncol(relative) == 2L && all(is.finite(relative) & 
                relative > 0), "INVALID_HAZARD_OR_RELATIVE_NO_CLIPPING")
            p <- engine$ff_aj_backend(hazard, relative)
            pry_predict_assert(is.matrix(p) && nrow(p) == nrow(relative) && ncol(p) == 3L && all(is.finite(p) & 
                p >= 0 & p <= 1) && all(abs(rowSums(p) - 1) <= tolerances$probability_mass_abs), "INVALID_THREE_STATE_PROBABILITY")
            colnames(p) <- c("disease", "death", "survival")
            p
        }
        pry_predict_piece_check <- function(piece, tau) {
            fit <- piece$fit
            pry_predict_assert(inherits(fit, "coxph.penal") && inherits(fit, "coxph") && !inherits(fit, "coxphms") && 
                identical(fit$method, "efron") && is.matrix(fit$x) && ncol(fit$x) == length(stats::coef(fit)) && 
                all(is.finite(stats::coef(fit))) && all(is.finite(stats::vcov(fit))) && is.null(fit$naive.var) && 
                !length(fit$offset) && is.null(fit$frail) && !length(fit$na.action), "FINITE_NATIVE_PENAL_EFRON_NONROBUST_FIT_REQUIRED")
            pry_predict_assert(isTRUE(piece$diagnostics$passed), "FITTING_DIAGNOSTICS_MUST_PASS_FIRST")
            labels <- attr(fit$terms, "term.labels")
            pry_predict_assert(any(grepl("strata\\(\\.pry_center\\)", labels)) && all(attr(fit$terms, "order") == 
                1L) && !length(attr(fit$terms, "offset")), "BOTH_CAUSES_EXPLICIT_CENTER_STRATA_NO_INTERACTIONS_OR_OFFSETS")
            support <- piece$center_support
            required <- c("center", "N", "max_followup", "at_risk_at_tau", "disease_events", "death_events", "supported_at_tau")
            pry_predict_assert(is.data.frame(support) && all(required %in% names(support)) && !anyDuplicated(as.character(support$center)), 
                "CENTER_SUPPORT_TABLE_REQUIRED")
            observed <- as.character(piece$observed_center_levels)
            pry_predict_assert(length(observed) > 0L && !anyNA(observed) && !anyDuplicated(observed) && setequal(observed, 
                as.character(support$center[support$N > 0])) && all(observed %in% as.character(piece$center_levels)), 
                "ACTUAL_VS_FROZEN_CENTER_LEVEL_BINDING")
            s <- support[match(observed, as.character(support$center)), , drop = FALSE]
            pry_predict_assert(all(is.finite(s$max_followup) & s$max_followup >= tau) && all(is.finite(s$at_risk_at_tau) & 
                s$at_risk_at_tau > 0) && all(s$supported_at_tau), "CENTER_LACKS_ACTUAL_TAU11_SUPPORT_NO_TAIL_EXTRAPOLATION")
            mf <- stats::model.frame(fit)
            native_center <- as.character(mf[["strata(.pry_center)"]])
            weight <- stats::model.weights(mf)
            pry_predict_assert(is.matrix(fit$y) && ncol(fit$y) == 2L && all(c("time", "status") %in% colnames(fit$y)) && 
                length(native_center) == nrow(fit$y) && setequal(native_center, observed) && length(weight) == 
                nrow(fit$y) && all(is.finite(weight) & weight > 0), "ACTUAL_NATIVE_RIGHT_CENSORED_SAMPLE_WEIGHT_CENTER_BINDING")
            for (g in observed) {
                ii <- which(native_center == g)
                row <- support[as.character(support$center) == g, , drop = FALSE]
                pry_predict_assert(row$N == length(ii) && row$max_followup == max(fit$y[ii, "time"]) && row$at_risk_at_tau == 
                  sum(fit$y[ii, "time"] >= tau) && row[[paste0(piece$cause, "_events")]] == sum(fit$y[ii, "status"] == 
                  1L), "SUPPORT_AND_ZERO_EVENT_COUNTS_MUST_MATCH_NATIVE_FITTED_ROWS")
            }
            pry_predict_assert(is.matrix(piece$reference_X) && is.numeric(piece$reference_X) && all(is.finite(piece$reference_X)) && 
                nrow(piece$reference_X) == length(piece$reference_center) && !anyDuplicated(as.character(piece$reference_center)) && 
                setequal(as.character(piece$reference_center), observed), "ONE_ACTUAL_FINITE_REFERENCE_ROW_PER_OBSERVED_CENTER_REQUIRED")
            invisible(observed)
        }
        pry_predict_design <- function(piece, X, center, transform_X) {
            pry_predict_assert(is.matrix(X) && is.numeric(X) && nrow(X) > 0L && all(is.finite(X)) && length(center) == 
                nrow(X) && !anyNA(center), "FINITE_ORIGINAL_DESIGN_AND_ROW_CENTERS_REQUIRED")
            center <- as.character(center)
            pry_predict_assert(all(center %in% as.character(piece$observed_center_levels)), "ABSENT_OR_UNKNOWN_CENTER_IS_NOT_ZERO_EVENT_CENTER")
            nd <- transform_X(X, piece$transform, center)
            pry_predict_assert(is.data.frame(nd) && nrow(nd) == nrow(X) && identical(as.character(nd$.pry_center), 
                center), "SAVED_TRANSFORM_ROW_OR_CENTER_DRIFT")
            rownames(nd) <- sprintf("pry_row_%08d", seq_len(nrow(nd)))
            mm <- pry_predict_capture(stats::model.matrix(piece$fit, data = nd), "NATIVE_MODEL_MATRIX_WARNING")
            pry_predict_assert(is.matrix(mm) && nrow(mm) == nrow(X) && identical(colnames(mm), colnames(piece$fit$x)) && 
                ncol(mm) == length(stats::coef(piece$fit)) && all(is.finite(mm)), "NATIVE_PENAL_DESIGN_COLUMN_POSITION_BINDING_FAILED")
            lp <- as.numeric(mm %*% as.numeric(stats::coef(piece$fit)))
            pry_predict_assert(all(is.finite(lp)), "NONFINITE_TRANSFORMED_LINEAR_PREDICTOR")
            list(newdata = nd, matrix = mm, lp = lp)
        }
        pry_predict_native_lp <- function(piece, nd) {
            lp <- pry_predict_capture(as.numeric(stats::predict(piece$fit, newdata = nd, type = "lp", reference = "zero", 
                se.fit = FALSE)), "NATIVE_PENAL_LP_WARNING")
            pry_predict_assert(length(lp) == nrow(nd) && all(is.finite(lp)), "NATIVE_PENAL_LP_SHAPE")
            lp
        }
        pry_predict_native_curves <- function(piece, nd) {
            sf <- pry_predict_capture(survival::survfit(piece$fit, newdata = nd, se.fit = FALSE, conf.type = "none", 
                stype = 2L, ctype = 2L, censor = TRUE), "NATIVE_PENAL_SURVFIT_WARNING")
            time <- sf$time
            H <- sf$cumhaz
            n <- nrow(nd)
            pry_predict_assert(is.numeric(time) && length(time) > 0L && all(is.finite(time) & time >= 0) && !is.null(H) && 
                all(is.finite(H)), "MISSING_NATIVE_CUMHAZ_NOT_ZERO_FILLED")
            curves <- vector("list", n)
            if (is.matrix(H)) {
                pry_predict_assert(nrow(H) == length(time) && ncol(H) == n, "NATIVE_MATRIX_CURVE_SHAPE")
                for (i in seq_len(n)) curves[[i]] <- list(time = as.numeric(time), cumulative = as.numeric(H[, 
                  i]))
            } else {
                counts <- sf$strata
                if (is.null(counts) && n == 1L) counts <- length(time)
                pry_predict_assert(is.numeric(H) && length(H) == length(time) && length(counts) == n && all(is.finite(counts) & 
                  counts > 0 & counts == as.integer(counts)) && sum(counts) == length(time), "NATIVE_PACKED_CURVE_ROW_LENGTHS")
                if (!is.null(names(counts))) pry_predict_assert(identical(names(counts), rownames(nd)), "NATIVE_PACKED_CURVE_ROW_ORDER")
                ends <- cumsum(counts)
                starts <- c(1L, head(ends, -1L) + 1L)
                for (i in seq_len(n)) {
                  ii <- seq.int(starts[i], ends[i])
                  curves[[i]] <- list(time = as.numeric(time[ii]), cumulative = as.numeric(H[ii]))
                }
            }
            for (curve in curves) pry_predict_assert(!is.unsorted(curve$time, strictly = TRUE) && length(curve$time) == 
                length(curve$cumulative) && all(curve$cumulative >= 0) && all(diff(curve$cumulative) >= 0), "NATIVE_CUMHAZ_NONNEGATIVE_NONDECREASING_REQUIRED")
            curves
        }
        pry_predict_merge_curves <- function(curves, tau) {
            pry_predict_assert(is.list(curves) && length(curves) == 2L, "TWO_CAUSE_CURVES_REQUIRED")
            grid <- sort(unique(unlist(lapply(curves, function(z) z$time[z$time <= tau]), use.names = FALSE)))
            if (!length(grid)) return(list(time = numeric(), hazard = matrix(numeric(), 0L, 2L)))
            cumulative <- vapply(curves, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
            cumulative <- matrix(cumulative, ncol = 2L)
            hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
                , drop = FALSE])
            pry_predict_assert(all(is.finite(hazard) & hazard >= 0), "INVALID_NATIVE_HAZARD_INCREMENT_NO_CLIPPING")
            keep <- rowSums(hazard) > 0
            list(time = grid[keep], hazard = hazard[keep, , drop = FALSE])
        }
        pry_compile_endpoint <- function(pieces, engine, tau = 11, transform_X = pry_transform_X, tolerances = pry_predict_tolerances(), 
            audit_callback = NULL) {
            pry_predict_assert(is.list(pieces) && identical(names(pieces), c("disease", "death")) && identical(as.numeric(tau), 
                11) && is.function(transform_X), "PAIRED_CAUSES_FIXED_TAU11_REQUIRED")
            pry_predict_assert(identical(tolerances, pry_predict_tolerances()), "PREDECLARED_TOLERANCES_MAY_NOT_DRIFT")
            tolerance_prerecord <- list(stage = "PREDICTION_TOLERANCES_FIXED_BEFORE_COMPARISON", tolerances = tolerances, 
                native_interface = "survfit(coxph.penal,newdata=explicit,se.fit=FALSE,stype=2,ctype=2)")
            if (!is.null(audit_callback)) {
                pry_predict_assert(is.function(audit_callback), "PREDICTION_AUDIT_CALLBACK_MUST_BE_FUNCTION")
                audit_callback(tolerance_prerecord)
            }
            kernel <- pry_predict_check_kernel(engine)
            observed <- lapply(pieces, pry_predict_piece_check, tau = tau)
            pry_predict_assert(identical(observed[[1L]], observed[[2L]]) && identical(pieces$disease$model_id, 
                pieces$death$model_id) && identical(pieces$disease$b, pieces$death$b) && identical(pieces$disease$endpoint, 
                pieces$death$endpoint), "PAIRED_CAUSE_SAMPLE_MODEL_ENDPOINT_BINDING")
            frames <- lapply(pieces, function(p) stats::model.frame(p$fit))
            pry_predict_assert(identical(as.numeric(pieces$disease$fit$y[, "time"]), as.numeric(pieces$death$fit$y[, 
                "time"])) && identical(as.numeric(stats::model.weights(frames$disease)), as.numeric(stats::model.weights(frames$death))) && 
                identical(as.character(frames$disease[["strata(.pry_center)"]]), as.character(frames$death[["strata(.pry_center)"]])), 
                "TWO_CAUSES_MUST_SHARE_EXACT_ROWS_WEIGHTS_TIMES_AND_CENTER_ORDER")
            levels <- observed[[1L]]
            anchors <- baselines <- vector("list", 2L)
            names(anchors) <- names(baselines) <- names(pieces)
            classes <- lapply(pieces, function(z) class(z$fit))
            lp_error <- 0
            for (cause in names(pieces)) {
                piece <- pieces[[cause]]
                anchors[[cause]] <- stats::setNames(numeric(length(levels)), levels)
                baselines[[cause]] <- stats::setNames(vector("list", length(levels)), levels)
                for (g in levels) {
                  row <- match(g, as.character(piece$reference_center))
                  prepared <- pry_predict_design(piece, piece$reference_X[row, , drop = FALSE], g, transform_X)
                  native <- pry_predict_native_lp(piece, prepared$newdata)
                  lp_error <- max(lp_error, pry_predict_close(prepared$lp, native, tolerances$native_lp_abs, 0, 
                    "REFERENCE_NATIVE_LP_EQUIVALENCE"))
                  curve <- pry_predict_native_curves(piece, prepared$newdata)[[1L]]
                  support <- piece$center_support[as.character(piece$center_support$center) == g, , drop = FALSE]
                  events <- support[[paste0(cause, "_events")]]
                  pry_predict_assert(length(events) == 1L && is.finite(events) && events >= 0, "ACTUAL_CENTER_CAUSE_EVENT_COUNT_REQUIRED")
                  if (events == 0) pry_predict_assert(all(curve$cumulative == 0), "ZERO_EVENT_CENTER_MUST_BE_NATIVE_ZERO_NOT_MISSING")
                  anchors[[cause]][[g]] <- prepared$lp
                  baselines[[cause]][[g]] <- curve
                }
            }
            hazards <- stats::setNames(lapply(levels, function(g) pry_predict_merge_curves(lapply(baselines, `[[`, 
                g), tau)), levels)
            state <- new.env(parent = emptyenv())
            state$passed <- FALSE
            state$receipt <- NULL
            structure(list(schema = "SELECTED44_PROTEIN_RIDGE_COMPILED_ENDPOINT_V1", pieces = pieces, endpoint = pieces$disease$endpoint, 
                model_id = pieces$disease$model_id, b = pieces$disease$b, tau = tau, center_levels = levels, anchors = anchors, 
                baselines = baselines, hazards = hazards, transform_X = transform_X, tolerances = tolerances, kernel_spec = kernel, 
                tolerance_prerecord = tolerance_prerecord, original_fit_classes = classes, reference_lp_max_abs = lp_error, 
                validation_state = state, compiled_from_explicit_finite_native_newdata = TRUE, native_stype = 2L, 
                native_ctype = 2L, no_clipping = TRUE), class = "pry_ridge_compiled_endpoint")
        }
        pry_predict_core <- function(compiled, X, center, engine, chunk_size = 512L, max_prediction_mb = 64, tau = compiled$tau) {
            pry_predict_assert(inherits(compiled, "pry_ridge_compiled_endpoint") && is.matrix(X) && nrow(X) > 0L && 
                length(center) == nrow(X) && !anyNA(center) && all(as.character(center) %in% compiled$center_levels), 
                "COMPILED_FIT_DESIGN_AND_KNOWN_PRESENT_CENTER_REQUIRED")
            pry_predict_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1 && chunk_size == 
                as.integer(chunk_size) && is.finite(max_prediction_mb) && max_prediction_mb > 0 && length(tau) == 
                1L && is.finite(tau) && tau >= 0 && tau <= compiled$tau, "BOUNDED_PREDICTION_REQUEST_REQUIRED")
            pry_predict_check_kernel(engine)
            largest <- max(vapply(compiled$pieces, function(p) ncol(p$fit$x), integer(1L)))
            step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * (largest + 64L))))))
            probability <- matrix(NA_real_, nrow(X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
            center <- as.character(center)
            for (g in unique(center)) {
                rows <- which(center == g)
                eng <- compiled$hazards[[g]]
                h <- eng$hazard[eng$time <= tau, , drop = FALSE]
                for (start in seq.int(1L, length(rows), by = step)) {
                  ii <- rows[seq.int(start, min(length(rows), start + step - 1L))]
                  lp <- vapply(names(compiled$pieces), function(cause) {
                    p <- pry_predict_design(compiled$pieces[[cause]], X[ii, , drop = FALSE], center[ii], compiled$transform_X)
                    p$lp - compiled$anchors[[cause]][[g]]
                  }, numeric(length(ii)))
                  relative <- exp(matrix(lp, nrow = length(ii), ncol = 2L))
                  probability[ii, ] <- pry_predict_kernel(engine, h, relative, compiled$tolerances)
                }
            }
            attr(probability, "chunk_size") <- step
            attr(probability, "integration_method") <- compiled$kernel_spec$method
            attr(probability, "no_clipping") <- TRUE
            probability
        }
        pry_predict_matrix <- function(compiled, X, center, engine, chunk_size = 512L, max_prediction_mb = 64) {
            pry_predict_assert(isTRUE(compiled$validation_state$passed), "NATIVE_AND_TWO_CAUSE_VALIDATION_REQUIRED_BEFORE_PRODUCTION")
            pry_predict_core(compiled, X, center, engine, chunk_size, max_prediction_mb)
        }
        pry_scalar_R_probability <- function(hazard) {
            pry_predict_assert(is.matrix(hazard) && ncol(hazard) == 2L && all(is.finite(hazard) & hazard >= 0), 
                "R_REFERENCE_HAZARD_INVALID")
            survival <- 1
            cif <- compensation <- c(0, 0)
            for (i in seq_len(nrow(hazard))) {
                a <- hazard[i, 1L]
                b <- hazard[i, 2L]
                if (a == 0 && b == 0) next
                total <- a + b
                if (is.finite(total)) {
                  share <- c(a, b)/total
                  mass <- -expm1(-total)
                  next_survival <- survival * exp(-total)
                } else {
                  scale <- max(a, b)
                  share <- (c(a, b)/scale)/sum(c(a, b)/scale)
                  mass <- 1
                  next_survival <- 0
                }
                increment <- survival * mass * share
                for (j in 1:2) {
                  y <- increment[j] - compensation[j]
                  t <- cif[j] + y
                  compensation[j] <- (t - cif[j]) - y
                  cif[j] <- t
                }
                survival <- next_survival
            }
            c(disease = cif[1L], death = cif[2L], survival = survival)
        }
        pry_validation_probe <- function(X, center, ethnicity, extra_counterfactual = list()) {
            pry_predict_assert(is.matrix(X) && "B" %in% colnames(X) && nrow(X) == length(center) && length(center) == 
                length(ethnicity) && !anyNA(center) && !anyNA(ethnicity), "VALIDATION_SCOPE_METADATA_REQUIRED")
            cells <- paste(center, ethnicity, sep = "\r")
            rows <- sort(unique(c(match(unique(cells), cells), as.integer(round(seq(1L, nrow(X), length.out = min(16L, 
                nrow(X))))))))
            original <- X[rows, , drop = FALSE]
            zero <- one <- original
            zero[, "B"] <- 0
            one[, "B"] <- 1
            make <- function(xx, label, kind) list(X = xx, center = as.character(center[rows]), ethnicity = as.character(ethnicity[rows]), 
                label = label, kind = kind)
            c(list(observed = make(original, "FIXED_OBSERVED_ROWS", "observed"), B0 = make(zero, "FIXED_B0_COUNTERFACTUAL_ROWS", 
                "counterfactual"), B1 = make(one, "FIXED_B1_COUNTERFACTUAL_ROWS", "counterfactual")), extra_counterfactual)
        }
        pry_validate_prediction <- function(compiled, cases, engine, expected_ethnicity, native_rows_per_call = 8L, 
            require_counterfactual = TRUE) {
            compiled$validation_state$passed <- FALSE
            pry_predict_assert(is.list(cases) && length(cases) > 0L && is.numeric(native_rows_per_call), "EXPLICIT_BOUNDED_VALIDATION_CASES_REQUIRED")
            pry_predict_assert(length(native_rows_per_call) == 1L && native_rows_per_call >= 1L && native_rows_per_call == 
                as.integer(native_rows_per_call) && native_rows_per_call <= 16L, "BOUNDED_NATIVE_CURVE_BATCH_MAX16")
            tolerance <- compiled$tolerances
            errors <- c(native_lp = compiled$reference_lp_max_abs, native_cumhaz = 0, compiled_vs_native_CIF = 0, 
                C_vs_scalar_R = 0, probability_mass = 0, time_monotonicity = 0)
            coverage <- list()
            rows_checked <- 0L
            kinds <- character()
            for (case in cases) {
                X <- case$X
                center <- as.character(case$center)
                ethnicity <- as.character(case$ethnicity)
                pry_predict_assert(is.matrix(X) && "B" %in% colnames(X) && nrow(X) == length(center) && length(ethnicity) == 
                  nrow(X) && !anyNA(ethnicity) && case$kind %in% c("observed", "counterfactual"), "FIXED_VALIDATION_ROW_METADATA_REQUIRED")
                kinds <- c(kinds, case$kind)
                coverage[[length(coverage) + 1L]] <- data.frame(center = center, ethnicity = ethnicity, B = X[, 
                  "B"], stringsAsFactors = FALSE)
                fast <- pry_predict_core(compiled, X, center, engine, chunk_size = native_rows_per_call)
                errors["probability_mass"] <- max(errors["probability_mass"], abs(rowSums(fast) - 1))
                for (g in unique(center)) {
                  group_rows <- which(center == g)
                  for (start in seq.int(1L, length(group_rows), by = native_rows_per_call)) {
                    ii <- group_rows[seq.int(start, min(length(group_rows), start + native_rows_per_call - 1L))]
                    native <- vector("list", 2L)
                    names(native) <- names(compiled$pieces)
                    for (cause in names(compiled$pieces)) {
                      piece <- compiled$pieces[[cause]]
                      prepared <- pry_predict_design(piece, X[ii, , drop = FALSE], center[ii], compiled$transform_X)
                      native_lp <- pry_predict_native_lp(piece, prepared$newdata)
                      errors["native_lp"] <- max(errors["native_lp"], pry_predict_close(prepared$lp, native_lp, 
                        tolerance$native_lp_abs, 0, "NATIVE_LP_POSITIONAL_EQUIVALENCE"))
                      native[[cause]] <- pry_predict_native_curves(piece, prepared$newdata)
                      for (k in seq_along(ii)) {
                        curve <- native[[cause]][[k]]
                        ref <- compiled$baselines[[cause]][[g]]
                        keep <- curve$time <= compiled$tau
                        expected <- c(0, ref$cumulative)[findInterval(curve$time[keep], ref$time) + 1L] * exp(prepared$lp[k] - 
                          compiled$anchors[[cause]][[g]])
                        errors["native_cumhaz"] <- max(errors["native_cumhaz"], pry_predict_close(curve$cumulative[keep], 
                          expected, tolerance$cumhaz_abs, tolerance$cumhaz_rel, "NATIVE_ANCHORED_CUMHAZ_EQUIVALENCE"))
                      }
                    }
                    for (k in seq_along(ii)) {
                      merged <- pry_predict_merge_curves(lapply(native, `[[`, k), compiled$tau)
                      native_probability <- pry_predict_kernel(engine, merged$hazard, matrix(c(1, 1), 1L), tolerance)
                      scalar <- matrix(pry_scalar_R_probability(merged$hazard), nrow = 1L)
                      actual <- matrix(as.numeric(fast[ii[k], ]), nrow = 1L)
                      errors["compiled_vs_native_CIF"] <- max(errors["compiled_vs_native_CIF"], pry_predict_close(actual, 
                        native_probability, tolerance$probability_abs, tolerance$probability_rel, "NATIVE_TWO_CAUSE_CIF_EQUIVALENCE"))
                      errors["C_vs_scalar_R"] <- max(errors["C_vs_scalar_R"], pry_predict_close(native_probability, 
                        scalar, tolerance$scalar_R_probability_abs, 0, "FROZEN_C_KERNEL_VS_SCALAR_R_REFERENCE"))
                    }
                  }
                }
                before <- pry_predict_core(compiled, X, center, engine, chunk_size = native_rows_per_call, tau = 0)
                for (time in c(0.25, 0.5, 0.75, 1) * compiled$tau) {
                  after <- pry_predict_core(compiled, X, center, engine, chunk_size = native_rows_per_call, tau = time)
                  delta <- after - before
                  errors["time_monotonicity"] <- max(errors["time_monotonicity"], -delta[, "disease"], -delta[, 
                    "death"], delta[, "survival"])
                  before <- after
                }
                rows_checked <- rows_checked + nrow(X)
            }
            scope <- do.call(rbind, coverage)
            pry_predict_assert(setequal(unique(scope$center), compiled$center_levels) && length(expected_ethnicity) > 
                0L && setequal(unique(scope$ethnicity), as.character(expected_ethnicity)), "EVERY_PRESENT_CENTER_AND_ETHNICITY_MUST_BE_VALIDATED")
            pry_predict_assert("observed" %in% kinds && (!require_counterfactual || "counterfactual" %in% kinds), 
                "OBSERVED_AND_COUNTERFACTUAL_VALIDATION_REQUIRED")
            pairs <- unique(scope[, c("center", "ethnicity"), drop = FALSE])
            for (k in seq_len(nrow(pairs))) pry_predict_assert(all(c(0, 1) %in% scope$B[scope$center == pairs$center[k] & 
                scope$ethnicity == pairs$ethnicity[k]]), "B0_AND_B1_REQUIRED_FOR_EVERY_VALIDATED_CENTER_ETHNICITY_CELL")
            pry_predict_assert(errors["time_monotonicity"] <= tolerance$time_monotonicity_abs && errors["probability_mass"] <= 
                tolerance$probability_mass_abs, "TIME_MONOTONICITY_OR_MASS_FAILURE")
            pry_predict_assert(identical(lapply(compiled$pieces, function(z) class(z$fit)), compiled$original_fit_classes), 
                "NATIVE_COXPH_PENAL_CLASSES_MUST_REMAIN_UNCHANGED")
            receipt <- list(schema = "SELECTED44_PROTEIN_RIDGE_NATIVE_PREDICTION_VALIDATION_V1", passed = TRUE, 
                rows_checked = rows_checked, fixed_cases = lapply(cases, function(z) list(label = z$label, kind = z$kind, 
                  rows = nrow(z$X))), centers = compiled$center_levels, ethnicities = as.character(expected_ethnicity), 
                B_interventions = c(0, 1), max_abs = as.list(errors), tolerances = tolerance, native_rows_per_call = native_rows_per_call, 
                native_interface = "survfit(coxph.penal,newdata=explicit,se.fit=FALSE,stype=2,ctype=2)", baseline_reference = "ACTUAL_FINITE_ROW_PER_CENTER;exp(LP_new-LP_anchor);same_saved_transform", 
                compiled_kernel = compiled$kernel_spec$method, original_classes = compiled$original_fit_classes, 
                zero_event_supported_centers = lapply(compiled$pieces, function(p) as.character(p$center_support$center[p$center_support$N > 
                  0 & p$center_support[[paste0(p$cause, "_events")]] == 0])), no_class_deletion = TRUE, no_clipping = TRUE, 
                no_NA_fill = TRUE, fits_performed = 0L, uses_RNG = FALSE)
            compiled$validation_state$receipt <- receipt
            compiled$validation_state$passed <- TRUE
            receipt
        }
    }), envir = context)
}

