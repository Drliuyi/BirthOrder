# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

pry_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

pry_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

pry_models <- c("S44_PROT_REF11", "S44_PROT_PLUS15")

pry_matrix <- function(X) {
    pry_assert(is.matrix(X) && is.numeric(X) && nrow(X) > 0L && ncol(X) > 0L && all(is.finite(X)) && !is.null(colnames(X)) && 
        !anyNA(colnames(X)) && !anyDuplicated(colnames(X)), "PRY_FINITE_ORIGINAL_LABELLED_MATRIX_REQUIRED")
    if ("(Intercept)" %in% colnames(X)) {
        pry_assert(all(X[, "(Intercept)"] == 1), "PRY_ORIGINAL_INTERCEPT_CHANGED")
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    X
}

pry_pool_X <- function(X, encoding = NULL) {
    X <- pry_matrix(X)
    asian <- "ethnicityAsian or Asian British"
    chinese <- "ethnicityChinese"
    pooled <- "ethnicityAsian_or_Chinese_pooled"
    raw <- all(c(asian, chinese) %in% colnames(X))
    mode <- if (raw) 
        "ORIGINAL_ASIAN_CHINESE_DUMMIES"
    else "ALREADY_APPROVED_ASIAN_CHINESE_POOL"
    if (!is.null(encoding)) 
        pry_assert(identical(mode, encoding), "PRY_FROZEN_INPUT_ETHNICITY_ENCODING_DRIFT")
    if (raw) {
        pry_assert(!pooled %in% colnames(X) && all(X[, asian] %in% 0:1) && all(X[, chinese] %in% 0:1) && all(X[, 
            asian] + X[, chinese] <= 1), "PRY_ORIGINAL_ASIAN_CHINESE_DISJOINT_DUMMIES_REQUIRED")
        Z <- X[, colnames(X) != chinese, drop = FALSE]
        Z[, asian] <- X[, asian] + X[, chinese]
        colnames(Z)[colnames(Z) == asian] <- pooled
    }
    else {
        pry_assert(pooled %in% colnames(X) && !any(c(asian, chinese) %in% colnames(X)), "PRY_EXACT_APPROVED_ETHNICITY_POOL_REQUIRED")
        Z <- X
    }
    list(X = Z, encoding = mode)
}

pry_weight_scale <- function(X, weights) {
    mu <- colSums(X * (weights/sum(weights)))
    sd <- sqrt(colSums(sweep(X, 2L, mu, `-`)^2 * (weights/sum(weights))))
    pry_assert(all(is.finite(mu)) && all(is.finite(sd) & sd > 0), "PRY_B0_CONSTANT_OR_UNSCALABLE_BLOCK_COLUMN")
    list(center = mu, scale = sd, rule = "B0_SELECTION_WEIGHTED_POPULATION_SD;FROZEN_FOR_ALL_B0_B1_100;RIDGE_SCALE_FALSE")
}

pry_check_centers <- function(X, center, columns) {
    pry_assert(length(center) == nrow(X) && !anyNA(center), "PRY_ORIGINAL_ACTUAL_CENTER_REQUIRED")
    center <- as.character(center)
    for (column in columns) {
        label <- sub("^center", "", column)
        pry_assert(all(X[, column] %in% 0:1) && identical(as.numeric(X[, column]), as.numeric(center == label)), 
            paste("PRY_ORIGINAL_CENTER_DUMMY_ACTUAL_BINDING_CHANGED", column))
    }
    center
}

pry_freeze_reference <- function(X_B0, center_B0, weights_B0, model_id, ethnicity_columns = NULL, body_columns = NULL, 
    df_accuracy = 0.1) {
    pry_assert(model_id %in% pry_models, "PRY_EXACT_TWO_PROTEIN_MODELS_ONLY")
    original <- pry_matrix(X_B0)
    pooled <- pry_pool_X(original)
    X <- pooled$X
    w <- as.numeric(weights_B0)
    pry_assert(length(w) == nrow(X) && all(is.finite(w) & w > 0), "PRY_EXACT_POSITIVE_B0_SELECTION_WEIGHTS")
    center_columns <- grep("^centerc[0-9]+$", colnames(X), value = TRUE)
    pry_assert(length(center_columns) > 0L && !length(setdiff(grep("^center", colnames(X), value = TRUE), center_columns)), 
        "PRY_FROZEN_CENTER_DUMMY_MAPPING_REQUIRED")
    center <- pry_check_centers(X, center_B0, center_columns)
    levels <- if (is.factor(center_B0)) 
        levels(center_B0)
    else sort(unique(center))
    pry_assert(all(center %in% levels) && length(levels) >= 2L, "PRY_FROZEN_CENTER_LEVELS_REQUIRED")
    if (is.null(ethnicity_columns)) 
        ethnicity_columns <- grep("^ethnicity", colnames(X), value = TRUE)
    if (is.null(body_columns)) 
        body_columns <- grep("^v5__(weight|fat_mass|fat_free_mass)__ns[0-9]+$", colnames(X), value = TRUE)
    pry_assert(length(ethnicity_columns) > 0L && length(body_columns) > 0L && !anyDuplicated(c(ethnicity_columns, 
        body_columns, center_columns)) && all(c(ethnicity_columns, body_columns) %in% colnames(X)), "PRY_EXACT_DISJOINT_ETHNIC_BODY_CENTER_BLOCKS")
    pry_assert(all(X[, ethnicity_columns, drop = FALSE] %in% 0:1) && all(rowSums(X[, ethnicity_columns, drop = FALSE]) <= 
        1L), "PRY_FROZEN_TREATMENT_CONTRASTS_NOT_INTEGER_ETHNICITY")
    pry_assert(all(colSums(X[, ethnicity_columns, drop = FALSE]) > 0) && any(rowSums(X[, ethnicity_columns, drop = FALSE]) == 
        0), "PRY_B0_ETHNIC_LEVEL_OR_REFERENCE_ABSENT")
    for (id in c("weight", "fat_mass", "fat_free_mass")) pry_assert(any(startsWith(body_columns, paste0("v5__", 
        id, "__ns"))), "PRY_WHOLE_THREE_VARIABLE_BODY_BLOCK_REQUIRED")
    other <- setdiff(colnames(X), c(center_columns, ethnicity_columns, body_columns))
    pry_assert(all(c("B", "age") %in% other), "PRY_B_AND_LINEAR_AGE_REMAIN_UNPENALIZED")
    ref <- list(schema = "SELECTED44_PROTEIN_RIDGE_B0_REFERENCE_V1", model_id = model_id, original_columns = colnames(original), 
        pooled_columns = colnames(X), ethnicity_encoding = pooled$encoding, center_columns = center_columns, center_levels = levels, 
        ethnicity_columns = ethnicity_columns, body_columns = body_columns, other_columns = other, ethnic = pry_weight_scale(X[, 
            ethnicity_columns, drop = FALSE], w), body = pry_weight_scale(X[, body_columns, drop = FALSE], w), 
        ethnicity_k = length(ethnicity_columns), body_k = length(body_columns), df_accuracy = df_accuracy, df_rule = "EACH_MODEL_ENDPOINT_CAUSE_DRAW_RECALIBRATES_DF_K_OVER_2", 
        control = list(iter.max = 100L, outer.max = 40L, eps = 1e-09), B0_design_hash = pry_hash(original), B0_weights_hash = pry_hash(w), 
        whole_function_space_preserved = TRUE, no_column_truncation = TRUE)
    pry_assert(is.numeric(df_accuracy) && length(df_accuracy) == 1L && df_accuracy > 0 && df_accuracy <= 0.1, "PRY_DF_ACCURACY_MUST_BE_PREFROZEN_AT_MOST_OFFICIAL_DEFAULT")
    ref$reference_sha256 <- pry_hash(ref)
    ref
}

pry_check_reference <- function(reference) {
    r <- reference
    wanted <- r$reference_sha256
    r$reference_sha256 <- NULL
    pry_assert(identical(reference$schema, "SELECTED44_PROTEIN_RIDGE_B0_REFERENCE_V1") && identical(pry_hash(r), 
        wanted), "PRY_FROZEN_B0_REFERENCE_CHANGED")
}

pry_screen_design <- function(X, center, weights, reference) {
    pry_check_reference(reference)
    original <- pry_matrix(X)
    pry_assert(identical(colnames(original), reference$original_columns), "PRY_ORIGINAL_X_COLUMN_ORDER_DRIFT")
    z <- pry_pool_X(original, reference$ethnicity_encoding)$X
    pry_assert(identical(colnames(z), reference$pooled_columns), "PRY_POOLED_X_COLUMN_ORDER_DRIFT")
    center <- pry_check_centers(z, center, reference$center_columns)
    pry_assert(all(center %in% reference$center_levels) && length(weights) == nrow(z) && all(is.finite(weights) & 
        weights > 0), "PRY_CENTER_WEIGHT_BINDING_INVALID")
    E <- z[, reference$ethnicity_columns, drop = FALSE]
    pry_assert(all(E %in% 0:1) && all(rowSums(E) <= 1), "PRY_ETHNIC_CONTRAST_ENCODING_CHANGED")
    span <- vapply(seq_len(ncol(z)), function(j) diff(range(z[, j])), 0)
    names(span) <- colnames(z)
    constant <- setdiff(names(span)[span == 0], reference$center_columns)
    absent <- reference$ethnicity_columns[colSums(E) == 0]
    if (!any(rowSums(E) == 0)) 
        absent <- c(absent, "FROZEN_REFERENCE_LEVEL")
    list(passed = !length(absent) && !length(constant), n = nrow(z), absent_ethnic_levels = absent, constant_noncenter_columns = constant, 
        observed_centers = reference$center_levels[reference$center_levels %in% center], ethnicity_k_unchanged = ncol(E), 
        target_df_unchanged = reference$ethnicity_k/2, X = z, original = original, center = center)
}

pry_preflight_fixed100 <- function(reference, design_provider) {
    pry_assert(is.function(design_provider), "PRY_FIXED_DRAW_PROVIDER_REQUIRED_NO_NEW_RESAMPLING")
    rows <- lapply(0:100, function(b) {
        z <- design_provider(b)
        s <- pry_screen_design(z$X, z$center, z$weights, reference)
        list(b = b, passed = s$passed, n = s$n, absent_ethnic_levels = s$absent_ethnic_levels, constant_noncenter_columns = s$constant_noncenter_columns, 
            ethnicity_k = s$ethnicity_k_unchanged, target_df = s$target_df_unchanged)
    })
    list(schema = "SELECTED44_PROTEIN_RIDGE_FIXED100_DESIGN_PREFLIGHT_V1", reference_sha256 = reference$reference_sha256, 
        passed = all(vapply(rows, `[[`, TRUE, "passed")), draws = rows, fits = 0L, new_draws = 0L, effective_df_attainability = "NATIVE_PER_FIT_OUTER_CONTROL_MUST_PASS;NOT_CLAIMED_BY_ZERO_FIT_DESIGN_SCREEN")
}

pry_build_transform <- function(screen, weights, reference, body_penalty) {
    B <- sweep(sweep(screen$X[, reference$body_columns, drop = FALSE], 2L, reference$body$center, `-`), 2L, reference$body$scale, 
        `/`)
    k <- ncol(B)
    if (isTRUE(body_penalty)) {
        R <- NULL
        permutation <- seq_len(k)
        forward <- inverse <- diag(k)
        witness <- 0
    }
    else {
        pry_assert(nrow(B) >= k, "PRY_BODY_QR_TOO_FEW_ROWS_NO_TRUNCATION")
        q <- qr(B * sqrt(weights/sum(weights)), LAPACK = TRUE)
        R <- qr.R(q)[seq_len(k), seq_len(k), drop = FALSE]
        permutation <- q$pivot
        pry_assert(all(is.finite(R)) && all(diag(R) != 0), "PRY_BODY_QR_NOT_INVERTIBLE_NO_TRUNCATION")
        P <- diag(k)[, permutation, drop = FALSE]
        forward <- P %*% backsolve(R, diag(k))
        inverse <- R %*% t(P)
        witness <- max(abs(forward %*% inverse - diag(k)))
        pry_assert(all(is.finite(forward)) && all(is.finite(inverse)) && witness <= 1e-08, "PRY_BODY_QR_INVERSE_WITNESS_FAILED_NO_TRUNCATION")
    }
    list(schema = "SELECTED44_PROTEIN_RIDGE_DESIGN_TRANSFORM_V1", reference = reference, body_penalty = isTRUE(body_penalty), 
        body_forward = forward, body_inverse = inverse, body_QR_R = R, body_column_permutation = permutation, body_inverse_max_abs_error = witness, 
        qr_definition = if (body_penalty) "NO_QR;PENALIZE_B0_STANDARDIZED_ORIGINAL_BODY_BASIS" else "sqrt(w/sum(w))*B0_standardized_original_BODY;LAPACK_pivoted_QR;ALL_COLUMNS", 
        other_safe = setNames(sprintf(".pry_x%03d", seq_along(reference$other_columns)), reference$other_columns), 
        body_safe = sprintf(".pry_q%03d", seq_len(k)), fallback_penalty_basis = "B0_STANDARDIZED_ORIGINAL_NS_BASIS_NOT_QR_COEFFICIENTS", 
        no_truncation = TRUE, no_penalty_mixing = TRUE)
}

pry_check_qr_equivalence <- function(X, time, status, center, weights, reference, endpoint, cause = c("disease", 
    "death"), prediction_X = NULL, prediction_center = NULL, tolerance = 1e-10) {
    pry_runtime_guard()
    cause <- match.arg(cause)
    s <- pry_screen_design(X, center, weights, reference)
    pry_assert(s$passed, "PRY_EQUIVALENCE_FROZEN_DESIGN_INVALID")
    tr <- pry_build_transform(s, weights, reference, FALSE)
    event <- as.integer(status == if (cause == "disease") 1L else 2L)
    raw <- as.data.frame(s$X[, setdiff(reference$pooled_columns, reference$center_columns), drop = FALSE])
    raw_names <- sprintf(".raw%03d", seq_len(ncol(raw)))
    names(raw) <- raw_names
    raw$.pry_center <- factor(s$center, levels = reference$center_levels)
    transformed <- pry_transform_X(X, tr, center)
    make_fit <- function(work, terms) {
        work$.pry_time <- as.numeric(time)
        work$.pry_event <- event
        work$.pry_weight <- as.numeric(weights)
        e <- new.env(parent = asNamespace("survival"))
        e$strata <- survival::strata
        f <- stats::reformulate(c(terms, "strata(.pry_center)"), response = "survival::Surv(.pry_time,.pry_event)", 
            env = e)
        warnings <- character()
        fit <- withCallingHandlers(survival::coxph(f, data = work, weights = .pry_weight, ties = "efron", robust = FALSE, 
            singular.ok = FALSE, na.action = stats::na.fail, x = TRUE, y = TRUE, model = TRUE, control = do.call(survival::coxph.control, 
                reference$control)), warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
        list(fit = fit, warnings = warnings)
    }
    a <- make_fit(raw, raw_names)
    z <- make_fit(transformed, c(unname(tr$other_safe), tr$body_safe, ".pry_E"))
    if (is.null(prediction_X)) {
        take <- unique(unlist(lapply(s$observed_centers, function(g) {
            k <- which(s$center == g)
            k[unique(round(seq(1L, length(k), length.out = min(16L, length(k)))))]
        })))
        prediction_X <- s$original[take, , drop = FALSE]
        prediction_center <- s$center[take]
    }
    pp <- pry_pool_X(prediction_X, reference$ethnicity_encoding)$X
    nd0 <- as.data.frame(pp[, setdiff(reference$pooled_columns, reference$center_columns), drop = FALSE])
    names(nd0) <- raw_names
    nd0$.pry_center <- factor(prediction_center, levels = reference$center_levels)
    nd1 <- pry_transform_X(prediction_X, tr, prediction_center)
    lp0 <- as.numeric(stats::predict(a$fit, newdata = nd0, type = "lp", reference = "zero"))
    lp1 <- as.numeric(stats::predict(z$fit, newdata = nd1, type = "lp", reference = "zero"))
    delta <- lp1 - lp0
    centered_delta <- delta
    for (g in unique(as.character(prediction_center))) {
        k <- as.character(prediction_center) == g
        centered_delta[k] <- delta[k] - delta[which(k)[1L]]
    }
    hazards <- vapply(seq_len(nrow(prediction_X)), function(j) {
        h0 <- survival::survfit(a$fit, newdata = nd0[j, , drop = FALSE], se.fit = FALSE, stype = 2L, ctype = 2L)
        h1 <- survival::survfit(z$fit, newdata = nd1[j, , drop = FALSE], se.fit = FALSE, stype = 2L, ctype = 2L)
        pry_assert(identical(as.numeric(h0$time), as.numeric(h1$time)), "PRY_QR_EQUIVALENCE_NATIVE_EVENT_GRID_CHANGED")
        max(abs(as.numeric(h0$cumhaz) - as.numeric(h1$cumhaz)), 0)
    }, 0)
    values <- c(lp_after_center_anchor = max(abs(centered_delta)), cumulative_hazard = max(hazards), final_loglik = abs(tail(a$fit$loglik, 
        1L) - tail(z$fit$loglik, 1L)))
    list(schema = "SELECTED44_PROTEIN_UNPENALIZED_QR_EQUIVALENCE_V1", endpoint = endpoint, cause = cause, passed = !length(c(a$warnings, 
        z$warnings)) && all(is.finite(values)) && all(values <= tolerance), differences = values, tolerance = tolerance, 
        tolerance_origin = "FROZEN_PROJECT_NATIVE_LP_BASELINE_1E_MINUS10_NO_POSTHOC_RELAXATION", original_warnings = a$warnings, 
        transformed_warnings = z$warnings, reference_sha256 = reference$reference_sha256, transform = tr, old_fixed_center_to_strata_equivalence_claimed = FALSE, 
        prediction_rows = nrow(prediction_X), fits = 2L, integrations = 0L, duplicated_rows_unchanged = TRUE)
}

pry_transform_X <- function(X, transform, center) {
    ref <- transform$reference
    pry_check_reference(ref)
    original <- pry_matrix(X)
    pry_assert(identical(colnames(original), ref$original_columns), "PRY_PREDICTION_ORIGINAL_COLUMNS_CHANGED")
    z <- pry_pool_X(original, ref$ethnicity_encoding)$X
    center <- pry_check_centers(z, center, ref$center_columns)
    pry_assert(all(center %in% ref$center_levels), "PRY_PREDICTION_UNKNOWN_FROZEN_CENTER")
    out <- as.data.frame(z[, ref$other_columns, drop = FALSE])
    names(out) <- unname(transform$other_safe)
    E <- sweep(sweep(z[, ref$ethnicity_columns, drop = FALSE], 2L, ref$ethnic$center, `-`), 2L, ref$ethnic$scale, 
        `/`)
    B <- sweep(sweep(z[, ref$body_columns, drop = FALSE], 2L, ref$body$center, `-`), 2L, ref$body$scale, `/`)
    if (isTRUE(transform$body_penalty)) 
        out$.pry_B <- I(B)
    else {
        Q <- B %*% transform$body_forward
        for (j in seq_len(ncol(Q))) out[[transform$body_safe[j]]] <- Q[, j]
    }
    out$.pry_E <- I(E)
    out$.pry_center <- factor(center, levels = ref$center_levels)
    pry_assert(all(is.finite(as.matrix(out[, setdiff(names(out), ".pry_center"), drop = FALSE]))), "PRY_TRANSFORM_NONFINITE")
    out
}

pry_runtime_guard <- function() {
    pry_assert(identical(Sys.getenv("PRY_EXECUTION_HOST"), "macmini") && Sys.info()[["sysname"]] == "Darwin" && 
        as.character(getRversion()) == "4.5.3" && as.character(utils::packageVersion("survival")) == "3.8.12", 
        "PRY_FITS_ONLY_ON_REVIEWED_MINI_R453_SURVIVAL3812")
}

pry_fit_piece <- function(X, time, status, center, weights, reference, model_id = reference$model_id, b, endpoint, 
    cause = c("disease", "death"), body_penalty = FALSE, ethnic_df_fraction = 1/2, body_df_fraction = 1/2) {
    pry_runtime_guard()
    cause <- match.arg(cause)
    pry_assert(is.logical(body_penalty) && length(body_penalty) == 1L && !is.na(body_penalty), "PRY_EXPLICIT_PRIMARY_OR_SOLE_AUTHORIZED_BODY_FALLBACK_REQUIRED")
    pry_assert(model_id %in% pry_models && identical(model_id, reference$model_id) && b %in% 0:100 && length(b) == 
        1L && b == as.integer(b) && endpoint %in% c("IHD", "T2D"), "PRY_EXACT_TWO_PROTEIN_FIXED_SCOPE")
    allowed <- c(1/3, 1/2, 2/3)
    valid_fraction <- function(value) is.numeric(value) && length(value) == 1L && is.finite(value) && value %in% 
        allowed
    pry_assert(valid_fraction(ethnic_df_fraction) && valid_fraction(body_df_fraction), "PRY_DF_FRACTIONS_ONLY_ONE_THIRD_ONE_HALF_TWO_THIRDS")
    changed <- c(ethnicity = ethnic_df_fraction != 1/2, body = body_df_fraction != 1/2)
    pry_assert(sum(changed) <= 1L && (b == 0L || !any(changed)), "PRY_B0_SINGLE_BLOCK_SENSITIVITY_ONLY_NUMBERED_DRAWS_FIXED_HALF")
    pry_assert(body_penalty || !changed[["body"]], "PRY_BODY_DF_SENSITIVITY_REQUIRES_BODY_PENALTY")
    primary_setting <- !any(changed)
    sensitivity <- list(schema = "SELECTED44_PROTEIN_RIDGE_B0_PENALTY_SENSITIVITY_V1", enabled = !primary_setting, 
        changed_block = if (primary_setting) "NONE" else names(changed)[changed], ethnic_df_fraction = ethnic_df_fraction, 
        body_df_fraction = body_df_fraction, ethnicity_k = reference$ethnicity_k, body_k = reference$body_k, body_penalty = body_penalty, 
        primary_gate_eligible = primary_setting, rule = "PRIMARY_ALL_202_HALF;B0_ONLY_ONE_ACTIVE_BLOCK_ONE_THIRD_OR_TWO_THIRDS;NO_K_CHANGE")
    s <- pry_screen_design(X, center, weights, reference)
    pry_assert(s$passed, paste("PRY_FIXED_ENCODING_ABSENT_OR_CONSTANT_REQUIRES_REVIEW_NOT_DELETE", paste(c(s$absent_ethnic_levels, 
        s$constant_noncenter_columns), collapse = ",")))
    pry_assert(length(time) == s$n && length(status) == s$n && all(is.finite(time) & time > 0) && !anyNA(status) && 
        all(status %in% 0:2), "PRY_ORIGINAL_FULLFOLLOWUP_RISK_REQUIRED")
    event <- as.integer(status == if (cause == "disease") 1L else 2L)
    pry_assert(sum(event) > 0L, "PRY_GLOBAL_CAUSE_HAS_NO_EVENTS")
    support <- do.call(rbind, lapply(reference$center_levels, function(g) {
        k <- s$center == g
        data.frame(center = g, N = sum(k), max_followup = if (any(k)) 
            max(time[k])
        else NA_real_, at_risk_at_tau = sum(time[k] >= 11), disease_events = sum(status[k] == 1L), death_events = sum(status[k] == 
            2L), supported_at_tau = any(time[k] >= 11), stringsAsFactors = FALSE)
    }))
    pry_assert(all(support$supported_at_tau[support$N > 0L]), "PRY_PRESENT_CENTER_TAU11_UNSUPPORTED_NO_EXTEND")
    transform <- pry_build_transform(s, weights, reference, body_penalty)
    work <- pry_transform_X(s$original, transform, s$center)
    work$.pry_time <- as.numeric(time)
    work$.pry_event <- event
    work$.pry_weight <- as.numeric(weights)
    environment <- new.env(parent = asNamespace("survival"))
    environment$ridge <- survival::ridge
    environment$strata <- survival::strata
    ordinary <- unname(transform$other_safe)
    if (!body_penalty) 
        ordinary <- c(ordinary, transform$body_safe)
    terms <- c(ordinary, sprintf("ridge(.pry_E, df = %.17g, eps = %.17g, scale = FALSE)", reference$ethnicity_k * 
        ethnic_df_fraction, reference$df_accuracy))
    if (body_penalty) 
        terms <- c(terms, sprintf("ridge(.pry_B, df = %.17g, eps = %.17g, scale = FALSE)", reference$body_k * body_df_fraction, 
            reference$df_accuracy))
    formula <- stats::reformulate(c(terms, "strata(.pry_center)"), response = "survival::Surv(.pry_time, .pry_event)", 
        env = environment)
    control <- do.call(survival::coxph.control, reference$control)
    warnings <- character()
    fit <- withCallingHandlers(survival::coxph(formula, data = work, weights = .pry_weight, ties = "efron", robust = FALSE, 
        singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE, control = control), 
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    pry_assert(inherits(fit, "coxph.penal") && inherits(fit, "coxph") && is.null(fit$naive.var) && !length(fit$call$cluster) && 
        identical(fit$method, "efron"), "PRY_NATIVE_PENALIZED_CLASS_AND_NO_CLUSTER_REQUIRED")
    pry_assert(fit$n == s$n && identical(as.numeric(fit$y[, "time"]), as.numeric(time)) && identical(as.numeric(fit$y[, 
        "status"]), as.numeric(event)) && identical(as.numeric(stats::model.weights(fit$model)), as.numeric(weights)), 
        "PRY_DUPLICATED_ROWS_TIME_EVENT_WEIGHT_DRIFT")
    pen <- which(fit$pterms > 0)
    targets <- c(ethnicity = reference$ethnicity_k * ethnic_df_fraction, if (body_penalty) c(body = reference$body_k * 
        body_df_fraction))
    block_k <- c(ethnicity = reference$ethnicity_k, if (body_penalty) c(body = reference$body_k))
    pry_assert(length(pen) == length(targets) && length(fit$history) == length(targets), "PRY_EXACT_SEPARATE_PENALTY_BLOCKS")
    penalty <- lapply(seq_along(pen), function(j) {
        h <- fit$history[[j]]
        history <- as.matrix(h$history)
        theta <- if (nrow(history)) 
            as.numeric(tail(history[, 1L], 1L))
        else NA_real_
        achieved <- as.numeric(fit$df[pen[j]])
        list(block = names(targets)[j], k = unname(block_k[j]), target_df = unname(targets[j]), actual_df = achieved, 
            theta = theta, next_theta = h$theta, outer_done = isTRUE(h$done), history = h, passed = isTRUE(h$done) && 
                is.finite(theta) && theta >= 0 && is.finite(achieved) && abs(achieved - targets[j]) < reference$df_accuracy)
    })
    names(penalty) <- names(targets)
    finite <- all(is.finite(stats::coef(fit))) && all(is.finite(stats::vcov(fit))) && all(is.finite(fit$linear.predictors)) && 
        all(is.finite(fit$loglik))
    passed <- finite && !length(warnings) && all(vapply(penalty, `[[`, TRUE, "passed")) && fit$iter[1L] < control$outer.max
    anchors <- match(s$observed_centers, s$center)
    mapping <- list(original_columns = reference$original_columns, pooled_columns = reference$pooled_columns, center_columns_removed_for_strata = reference$center_columns, 
        other_safe = transform$other_safe, ethnicity_columns = reference$ethnicity_columns, body_columns = reference$body_columns, 
        native_design_columns = colnames(fit$x), native_coefficient_names = names(stats::coef(fit)), coefficient_alignment = "POSITIONAL_FIT_X_MODEL_MATRIX_ORDER;RIDGE_MATRIX_COEF_NAMES_MAY_REPEAT", 
        original_body_beta = if (body_penalty) "beta_std / B0_sd" else "body_forward %*% beta_QR / B0_sd", lp_origin_shift = "Subtract frozen block means on BOTH fitting/prediction; use native fitted-coordinate baseline")
    list(schema = "SELECTED44_PROTEIN_RIDGE_FIT_PIECE_V1", fit = fit, transform = transform, mapping = mapping, 
        reference = reference, center_levels = reference$center_levels, observed_center_levels = s$observed_centers, 
        center_support = support, reference_X = s$original[anchors, , drop = FALSE], reference_center = s$center[anchors], 
        model_id = model_id, b = as.integer(b), endpoint = endpoint, cause = cause, method_id = if (body_penalty) "PROTEIN_ETHNIC_BODY_RIDGE_CENTER_STRATA_V1" else "PROTEIN_ETHNIC_RIDGE_CENTER_STRATA_QR_V1", 
        sensitivity_setting = sensitivity, primary_gate_eligible = primary_setting, diagnostics = list(passed = passed, 
            finite = finite, warnings = warnings, penalties = penalty, targets = as.list(targets), sensitivity_setting = sensitivity, 
            primary_gate_eligible = primary_setting, iterations = fit$iter, control = control, n = s$n, events = sum(event), 
            weight_sum = sum(weights), original_X_hash = pry_hash(s$original), weights_hash = pry_hash(as.numeric(weights)), 
            time_status_hash = pry_hash(list(time = time, status = status)), robust = FALSE, cluster_argument_used = FALSE, 
            duplicates_collapsed = FALSE, inference = "EXTERNAL_FROZEN_FAMILY_BOOTSTRAP_NOT_MODEL_COEFFICIENT_SE", 
            R = R.version.string, survival = as.character(utils::packageVersion("survival"))))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        pry_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        pry_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        pry_models <- c("S44_PROT_REF11", "S44_PROT_PLUS15")
        pry_matrix <- function(X) {
            pry_assert(is.matrix(X) && is.numeric(X) && nrow(X) > 0L && ncol(X) > 0L && all(is.finite(X)) && !is.null(colnames(X)) && 
                !anyNA(colnames(X)) && !anyDuplicated(colnames(X)), "PRY_FINITE_ORIGINAL_LABELLED_MATRIX_REQUIRED")
            if ("(Intercept)" %in% colnames(X)) {
                pry_assert(all(X[, "(Intercept)"] == 1), "PRY_ORIGINAL_INTERCEPT_CHANGED")
                X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            }
            X
        }
        pry_pool_X <- function(X, encoding = NULL) {
            X <- pry_matrix(X)
            asian <- "ethnicityAsian or Asian British"
            chinese <- "ethnicityChinese"
            pooled <- "ethnicityAsian_or_Chinese_pooled"
            raw <- all(c(asian, chinese) %in% colnames(X))
            mode <- if (raw) "ORIGINAL_ASIAN_CHINESE_DUMMIES" else "ALREADY_APPROVED_ASIAN_CHINESE_POOL"
            if (!is.null(encoding)) pry_assert(identical(mode, encoding), "PRY_FROZEN_INPUT_ETHNICITY_ENCODING_DRIFT")
            if (raw) {
                pry_assert(!pooled %in% colnames(X) && all(X[, asian] %in% 0:1) && all(X[, chinese] %in% 0:1) && 
                  all(X[, asian] + X[, chinese] <= 1), "PRY_ORIGINAL_ASIAN_CHINESE_DISJOINT_DUMMIES_REQUIRED")
                Z <- X[, colnames(X) != chinese, drop = FALSE]
                Z[, asian] <- X[, asian] + X[, chinese]
                colnames(Z)[colnames(Z) == asian] <- pooled
            } else {
                pry_assert(pooled %in% colnames(X) && !any(c(asian, chinese) %in% colnames(X)), "PRY_EXACT_APPROVED_ETHNICITY_POOL_REQUIRED")
                Z <- X
            }
            list(X = Z, encoding = mode)
        }
        pry_weight_scale <- function(X, weights) {
            mu <- colSums(X * (weights/sum(weights)))
            sd <- sqrt(colSums(sweep(X, 2L, mu, `-`)^2 * (weights/sum(weights))))
            pry_assert(all(is.finite(mu)) && all(is.finite(sd) & sd > 0), "PRY_B0_CONSTANT_OR_UNSCALABLE_BLOCK_COLUMN")
            list(center = mu, scale = sd, rule = "B0_SELECTION_WEIGHTED_POPULATION_SD;FROZEN_FOR_ALL_B0_B1_100;RIDGE_SCALE_FALSE")
        }
        pry_check_centers <- function(X, center, columns) {
            pry_assert(length(center) == nrow(X) && !anyNA(center), "PRY_ORIGINAL_ACTUAL_CENTER_REQUIRED")
            center <- as.character(center)
            for (column in columns) {
                label <- sub("^center", "", column)
                pry_assert(all(X[, column] %in% 0:1) && identical(as.numeric(X[, column]), as.numeric(center == 
                  label)), paste("PRY_ORIGINAL_CENTER_DUMMY_ACTUAL_BINDING_CHANGED", column))
            }
            center
        }
        pry_freeze_reference <- function(X_B0, center_B0, weights_B0, model_id, ethnicity_columns = NULL, body_columns = NULL, 
            df_accuracy = 0.1) {
            pry_assert(model_id %in% pry_models, "PRY_EXACT_TWO_PROTEIN_MODELS_ONLY")
            original <- pry_matrix(X_B0)
            pooled <- pry_pool_X(original)
            X <- pooled$X
            w <- as.numeric(weights_B0)
            pry_assert(length(w) == nrow(X) && all(is.finite(w) & w > 0), "PRY_EXACT_POSITIVE_B0_SELECTION_WEIGHTS")
            center_columns <- grep("^centerc[0-9]+$", colnames(X), value = TRUE)
            pry_assert(length(center_columns) > 0L && !length(setdiff(grep("^center", colnames(X), value = TRUE), 
                center_columns)), "PRY_FROZEN_CENTER_DUMMY_MAPPING_REQUIRED")
            center <- pry_check_centers(X, center_B0, center_columns)
            levels <- if (is.factor(center_B0)) levels(center_B0) else sort(unique(center))
            pry_assert(all(center %in% levels) && length(levels) >= 2L, "PRY_FROZEN_CENTER_LEVELS_REQUIRED")
            if (is.null(ethnicity_columns)) ethnicity_columns <- grep("^ethnicity", colnames(X), value = TRUE)
            if (is.null(body_columns)) body_columns <- grep("^v5__(weight|fat_mass|fat_free_mass)__ns[0-9]+$", 
                colnames(X), value = TRUE)
            pry_assert(length(ethnicity_columns) > 0L && length(body_columns) > 0L && !anyDuplicated(c(ethnicity_columns, 
                body_columns, center_columns)) && all(c(ethnicity_columns, body_columns) %in% colnames(X)), "PRY_EXACT_DISJOINT_ETHNIC_BODY_CENTER_BLOCKS")
            pry_assert(all(X[, ethnicity_columns, drop = FALSE] %in% 0:1) && all(rowSums(X[, ethnicity_columns, 
                drop = FALSE]) <= 1L), "PRY_FROZEN_TREATMENT_CONTRASTS_NOT_INTEGER_ETHNICITY")
            pry_assert(all(colSums(X[, ethnicity_columns, drop = FALSE]) > 0) && any(rowSums(X[, ethnicity_columns, 
                drop = FALSE]) == 0), "PRY_B0_ETHNIC_LEVEL_OR_REFERENCE_ABSENT")
            for (id in c("weight", "fat_mass", "fat_free_mass")) pry_assert(any(startsWith(body_columns, paste0("v5__", 
                id, "__ns"))), "PRY_WHOLE_THREE_VARIABLE_BODY_BLOCK_REQUIRED")
            other <- setdiff(colnames(X), c(center_columns, ethnicity_columns, body_columns))
            pry_assert(all(c("B", "age") %in% other), "PRY_B_AND_LINEAR_AGE_REMAIN_UNPENALIZED")
            ref <- list(schema = "SELECTED44_PROTEIN_RIDGE_B0_REFERENCE_V1", model_id = model_id, original_columns = colnames(original), 
                pooled_columns = colnames(X), ethnicity_encoding = pooled$encoding, center_columns = center_columns, 
                center_levels = levels, ethnicity_columns = ethnicity_columns, body_columns = body_columns, other_columns = other, 
                ethnic = pry_weight_scale(X[, ethnicity_columns, drop = FALSE], w), body = pry_weight_scale(X[, 
                  body_columns, drop = FALSE], w), ethnicity_k = length(ethnicity_columns), body_k = length(body_columns), 
                df_accuracy = df_accuracy, df_rule = "EACH_MODEL_ENDPOINT_CAUSE_DRAW_RECALIBRATES_DF_K_OVER_2", 
                control = list(iter.max = 100L, outer.max = 40L, eps = 1e-09), B0_design_hash = pry_hash(original), 
                B0_weights_hash = pry_hash(w), whole_function_space_preserved = TRUE, no_column_truncation = TRUE)
            pry_assert(is.numeric(df_accuracy) && length(df_accuracy) == 1L && df_accuracy > 0 && df_accuracy <= 
                0.1, "PRY_DF_ACCURACY_MUST_BE_PREFROZEN_AT_MOST_OFFICIAL_DEFAULT")
            ref$reference_sha256 <- pry_hash(ref)
            ref
        }
        pry_check_reference <- function(reference) {
            r <- reference
            wanted <- r$reference_sha256
            r$reference_sha256 <- NULL
            pry_assert(identical(reference$schema, "SELECTED44_PROTEIN_RIDGE_B0_REFERENCE_V1") && identical(pry_hash(r), 
                wanted), "PRY_FROZEN_B0_REFERENCE_CHANGED")
        }
        pry_screen_design <- function(X, center, weights, reference) {
            pry_check_reference(reference)
            original <- pry_matrix(X)
            pry_assert(identical(colnames(original), reference$original_columns), "PRY_ORIGINAL_X_COLUMN_ORDER_DRIFT")
            z <- pry_pool_X(original, reference$ethnicity_encoding)$X
            pry_assert(identical(colnames(z), reference$pooled_columns), "PRY_POOLED_X_COLUMN_ORDER_DRIFT")
            center <- pry_check_centers(z, center, reference$center_columns)
            pry_assert(all(center %in% reference$center_levels) && length(weights) == nrow(z) && all(is.finite(weights) & 
                weights > 0), "PRY_CENTER_WEIGHT_BINDING_INVALID")
            E <- z[, reference$ethnicity_columns, drop = FALSE]
            pry_assert(all(E %in% 0:1) && all(rowSums(E) <= 1), "PRY_ETHNIC_CONTRAST_ENCODING_CHANGED")
            span <- vapply(seq_len(ncol(z)), function(j) diff(range(z[, j])), 0)
            names(span) <- colnames(z)
            constant <- setdiff(names(span)[span == 0], reference$center_columns)
            absent <- reference$ethnicity_columns[colSums(E) == 0]
            if (!any(rowSums(E) == 0)) absent <- c(absent, "FROZEN_REFERENCE_LEVEL")
            list(passed = !length(absent) && !length(constant), n = nrow(z), absent_ethnic_levels = absent, constant_noncenter_columns = constant, 
                observed_centers = reference$center_levels[reference$center_levels %in% center], ethnicity_k_unchanged = ncol(E), 
                target_df_unchanged = reference$ethnicity_k/2, X = z, original = original, center = center)
        }
        pry_preflight_fixed100 <- function(reference, design_provider) {
            pry_assert(is.function(design_provider), "PRY_FIXED_DRAW_PROVIDER_REQUIRED_NO_NEW_RESAMPLING")
            rows <- lapply(0:100, function(b) {
                z <- design_provider(b)
                s <- pry_screen_design(z$X, z$center, z$weights, reference)
                list(b = b, passed = s$passed, n = s$n, absent_ethnic_levels = s$absent_ethnic_levels, constant_noncenter_columns = s$constant_noncenter_columns, 
                  ethnicity_k = s$ethnicity_k_unchanged, target_df = s$target_df_unchanged)
            })
            list(schema = "SELECTED44_PROTEIN_RIDGE_FIXED100_DESIGN_PREFLIGHT_V1", reference_sha256 = reference$reference_sha256, 
                passed = all(vapply(rows, `[[`, TRUE, "passed")), draws = rows, fits = 0L, new_draws = 0L, effective_df_attainability = "NATIVE_PER_FIT_OUTER_CONTROL_MUST_PASS;NOT_CLAIMED_BY_ZERO_FIT_DESIGN_SCREEN")
        }
        pry_build_transform <- function(screen, weights, reference, body_penalty) {
            B <- sweep(sweep(screen$X[, reference$body_columns, drop = FALSE], 2L, reference$body$center, `-`), 
                2L, reference$body$scale, `/`)
            k <- ncol(B)
            if (isTRUE(body_penalty)) {
                R <- NULL
                permutation <- seq_len(k)
                forward <- inverse <- diag(k)
                witness <- 0
            } else {
                pry_assert(nrow(B) >= k, "PRY_BODY_QR_TOO_FEW_ROWS_NO_TRUNCATION")
                q <- qr(B * sqrt(weights/sum(weights)), LAPACK = TRUE)
                R <- qr.R(q)[seq_len(k), seq_len(k), drop = FALSE]
                permutation <- q$pivot
                pry_assert(all(is.finite(R)) && all(diag(R) != 0), "PRY_BODY_QR_NOT_INVERTIBLE_NO_TRUNCATION")
                P <- diag(k)[, permutation, drop = FALSE]
                forward <- P %*% backsolve(R, diag(k))
                inverse <- R %*% t(P)
                witness <- max(abs(forward %*% inverse - diag(k)))
                pry_assert(all(is.finite(forward)) && all(is.finite(inverse)) && witness <= 1e-08, "PRY_BODY_QR_INVERSE_WITNESS_FAILED_NO_TRUNCATION")
            }
            list(schema = "SELECTED44_PROTEIN_RIDGE_DESIGN_TRANSFORM_V1", reference = reference, body_penalty = isTRUE(body_penalty), 
                body_forward = forward, body_inverse = inverse, body_QR_R = R, body_column_permutation = permutation, 
                body_inverse_max_abs_error = witness, qr_definition = if (body_penalty) "NO_QR;PENALIZE_B0_STANDARDIZED_ORIGINAL_BODY_BASIS" else "sqrt(w/sum(w))*B0_standardized_original_BODY;LAPACK_pivoted_QR;ALL_COLUMNS", 
                other_safe = setNames(sprintf(".pry_x%03d", seq_along(reference$other_columns)), reference$other_columns), 
                body_safe = sprintf(".pry_q%03d", seq_len(k)), fallback_penalty_basis = "B0_STANDARDIZED_ORIGINAL_NS_BASIS_NOT_QR_COEFFICIENTS", 
                no_truncation = TRUE, no_penalty_mixing = TRUE)
        }
        pry_check_qr_equivalence <- function(X, time, status, center, weights, reference, endpoint, cause = c("disease", 
            "death"), prediction_X = NULL, prediction_center = NULL, tolerance = 1e-10) {
            pry_runtime_guard()
            cause <- match.arg(cause)
            s <- pry_screen_design(X, center, weights, reference)
            pry_assert(s$passed, "PRY_EQUIVALENCE_FROZEN_DESIGN_INVALID")
            tr <- pry_build_transform(s, weights, reference, FALSE)
            event <- as.integer(status == if (cause == "disease") 1L else 2L)
            raw <- as.data.frame(s$X[, setdiff(reference$pooled_columns, reference$center_columns), drop = FALSE])
            raw_names <- sprintf(".raw%03d", seq_len(ncol(raw)))
            names(raw) <- raw_names
            raw$.pry_center <- factor(s$center, levels = reference$center_levels)
            transformed <- pry_transform_X(X, tr, center)
            make_fit <- function(work, terms) {
                work$.pry_time <- as.numeric(time)
                work$.pry_event <- event
                work$.pry_weight <- as.numeric(weights)
                e <- new.env(parent = asNamespace("survival"))
                e$strata <- survival::strata
                f <- stats::reformulate(c(terms, "strata(.pry_center)"), response = "survival::Surv(.pry_time,.pry_event)", 
                  env = e)
                warnings <- character()
                fit <- withCallingHandlers(survival::coxph(f, data = work, weights = .pry_weight, ties = "efron", 
                  robust = FALSE, singular.ok = FALSE, na.action = stats::na.fail, x = TRUE, y = TRUE, model = TRUE, 
                  control = do.call(survival::coxph.control, reference$control)), warning = function(w) {
                  warnings <<- c(warnings, conditionMessage(w))
                  invokeRestart("muffleWarning")
                })
                list(fit = fit, warnings = warnings)
            }
            a <- make_fit(raw, raw_names)
            z <- make_fit(transformed, c(unname(tr$other_safe), tr$body_safe, ".pry_E"))
            if (is.null(prediction_X)) {
                take <- unique(unlist(lapply(s$observed_centers, function(g) {
                  k <- which(s$center == g)
                  k[unique(round(seq(1L, length(k), length.out = min(16L, length(k)))))]
                })))
                prediction_X <- s$original[take, , drop = FALSE]
                prediction_center <- s$center[take]
            }
            pp <- pry_pool_X(prediction_X, reference$ethnicity_encoding)$X
            nd0 <- as.data.frame(pp[, setdiff(reference$pooled_columns, reference$center_columns), drop = FALSE])
            names(nd0) <- raw_names
            nd0$.pry_center <- factor(prediction_center, levels = reference$center_levels)
            nd1 <- pry_transform_X(prediction_X, tr, prediction_center)
            lp0 <- as.numeric(stats::predict(a$fit, newdata = nd0, type = "lp", reference = "zero"))
            lp1 <- as.numeric(stats::predict(z$fit, newdata = nd1, type = "lp", reference = "zero"))
            delta <- lp1 - lp0
            centered_delta <- delta
            for (g in unique(as.character(prediction_center))) {
                k <- as.character(prediction_center) == g
                centered_delta[k] <- delta[k] - delta[which(k)[1L]]
            }
            hazards <- vapply(seq_len(nrow(prediction_X)), function(j) {
                h0 <- survival::survfit(a$fit, newdata = nd0[j, , drop = FALSE], se.fit = FALSE, stype = 2L, ctype = 2L)
                h1 <- survival::survfit(z$fit, newdata = nd1[j, , drop = FALSE], se.fit = FALSE, stype = 2L, ctype = 2L)
                pry_assert(identical(as.numeric(h0$time), as.numeric(h1$time)), "PRY_QR_EQUIVALENCE_NATIVE_EVENT_GRID_CHANGED")
                max(abs(as.numeric(h0$cumhaz) - as.numeric(h1$cumhaz)), 0)
            }, 0)
            values <- c(lp_after_center_anchor = max(abs(centered_delta)), cumulative_hazard = max(hazards), final_loglik = abs(tail(a$fit$loglik, 
                1L) - tail(z$fit$loglik, 1L)))
            list(schema = "SELECTED44_PROTEIN_UNPENALIZED_QR_EQUIVALENCE_V1", endpoint = endpoint, cause = cause, 
                passed = !length(c(a$warnings, z$warnings)) && all(is.finite(values)) && all(values <= tolerance), 
                differences = values, tolerance = tolerance, tolerance_origin = "FROZEN_PROJECT_NATIVE_LP_BASELINE_1E_MINUS10_NO_POSTHOC_RELAXATION", 
                original_warnings = a$warnings, transformed_warnings = z$warnings, reference_sha256 = reference$reference_sha256, 
                transform = tr, old_fixed_center_to_strata_equivalence_claimed = FALSE, prediction_rows = nrow(prediction_X), 
                fits = 2L, integrations = 0L, duplicated_rows_unchanged = TRUE)
        }
        pry_transform_X <- function(X, transform, center) {
            ref <- transform$reference
            pry_check_reference(ref)
            original <- pry_matrix(X)
            pry_assert(identical(colnames(original), ref$original_columns), "PRY_PREDICTION_ORIGINAL_COLUMNS_CHANGED")
            z <- pry_pool_X(original, ref$ethnicity_encoding)$X
            center <- pry_check_centers(z, center, ref$center_columns)
            pry_assert(all(center %in% ref$center_levels), "PRY_PREDICTION_UNKNOWN_FROZEN_CENTER")
            out <- as.data.frame(z[, ref$other_columns, drop = FALSE])
            names(out) <- unname(transform$other_safe)
            E <- sweep(sweep(z[, ref$ethnicity_columns, drop = FALSE], 2L, ref$ethnic$center, `-`), 2L, ref$ethnic$scale, 
                `/`)
            B <- sweep(sweep(z[, ref$body_columns, drop = FALSE], 2L, ref$body$center, `-`), 2L, ref$body$scale, 
                `/`)
            if (isTRUE(transform$body_penalty)) out$.pry_B <- I(B) else {
                Q <- B %*% transform$body_forward
                for (j in seq_len(ncol(Q))) out[[transform$body_safe[j]]] <- Q[, j]
            }
            out$.pry_E <- I(E)
            out$.pry_center <- factor(center, levels = ref$center_levels)
            pry_assert(all(is.finite(as.matrix(out[, setdiff(names(out), ".pry_center"), drop = FALSE]))), "PRY_TRANSFORM_NONFINITE")
            out
        }
        pry_runtime_guard <- function() {
            pry_assert(identical(Sys.getenv("PRY_EXECUTION_HOST"), "macmini") && Sys.info()[["sysname"]] == "Darwin" && 
                as.character(getRversion()) == "4.5.3" && as.character(utils::packageVersion("survival")) == "3.8.12", 
                "PRY_FITS_ONLY_ON_REVIEWED_MINI_R453_SURVIVAL3812")
        }
        pry_fit_piece <- function(X, time, status, center, weights, reference, model_id = reference$model_id, b, 
            endpoint, cause = c("disease", "death"), body_penalty = FALSE, ethnic_df_fraction = 1/2, body_df_fraction = 1/2) {
            pry_runtime_guard()
            cause <- match.arg(cause)
            pry_assert(is.logical(body_penalty) && length(body_penalty) == 1L && !is.na(body_penalty), "PRY_EXPLICIT_PRIMARY_OR_SOLE_AUTHORIZED_BODY_FALLBACK_REQUIRED")
            pry_assert(model_id %in% pry_models && identical(model_id, reference$model_id) && b %in% 0:100 && length(b) == 
                1L && b == as.integer(b) && endpoint %in% c("IHD", "T2D"), "PRY_EXACT_TWO_PROTEIN_FIXED_SCOPE")
            allowed <- c(1/3, 1/2, 2/3)
            valid_fraction <- function(value) is.numeric(value) && length(value) == 1L && is.finite(value) && value %in% 
                allowed
            pry_assert(valid_fraction(ethnic_df_fraction) && valid_fraction(body_df_fraction), "PRY_DF_FRACTIONS_ONLY_ONE_THIRD_ONE_HALF_TWO_THIRDS")
            changed <- c(ethnicity = ethnic_df_fraction != 1/2, body = body_df_fraction != 1/2)
            pry_assert(sum(changed) <= 1L && (b == 0L || !any(changed)), "PRY_B0_SINGLE_BLOCK_SENSITIVITY_ONLY_NUMBERED_DRAWS_FIXED_HALF")
            pry_assert(body_penalty || !changed[["body"]], "PRY_BODY_DF_SENSITIVITY_REQUIRES_BODY_PENALTY")
            primary_setting <- !any(changed)
            sensitivity <- list(schema = "SELECTED44_PROTEIN_RIDGE_B0_PENALTY_SENSITIVITY_V1", enabled = !primary_setting, 
                changed_block = if (primary_setting) "NONE" else names(changed)[changed], ethnic_df_fraction = ethnic_df_fraction, 
                body_df_fraction = body_df_fraction, ethnicity_k = reference$ethnicity_k, body_k = reference$body_k, 
                body_penalty = body_penalty, primary_gate_eligible = primary_setting, rule = "PRIMARY_ALL_202_HALF;B0_ONLY_ONE_ACTIVE_BLOCK_ONE_THIRD_OR_TWO_THIRDS;NO_K_CHANGE")
            s <- pry_screen_design(X, center, weights, reference)
            pry_assert(s$passed, paste("PRY_FIXED_ENCODING_ABSENT_OR_CONSTANT_REQUIRES_REVIEW_NOT_DELETE", paste(c(s$absent_ethnic_levels, 
                s$constant_noncenter_columns), collapse = ",")))
            pry_assert(length(time) == s$n && length(status) == s$n && all(is.finite(time) & time > 0) && !anyNA(status) && 
                all(status %in% 0:2), "PRY_ORIGINAL_FULLFOLLOWUP_RISK_REQUIRED")
            event <- as.integer(status == if (cause == "disease") 1L else 2L)
            pry_assert(sum(event) > 0L, "PRY_GLOBAL_CAUSE_HAS_NO_EVENTS")
            support <- do.call(rbind, lapply(reference$center_levels, function(g) {
                k <- s$center == g
                data.frame(center = g, N = sum(k), max_followup = if (any(k)) max(time[k]) else NA_real_, at_risk_at_tau = sum(time[k] >= 
                  11), disease_events = sum(status[k] == 1L), death_events = sum(status[k] == 2L), supported_at_tau = any(time[k] >= 
                  11), stringsAsFactors = FALSE)
            }))
            pry_assert(all(support$supported_at_tau[support$N > 0L]), "PRY_PRESENT_CENTER_TAU11_UNSUPPORTED_NO_EXTEND")
            transform <- pry_build_transform(s, weights, reference, body_penalty)
            work <- pry_transform_X(s$original, transform, s$center)
            work$.pry_time <- as.numeric(time)
            work$.pry_event <- event
            work$.pry_weight <- as.numeric(weights)
            environment <- new.env(parent = asNamespace("survival"))
            environment$ridge <- survival::ridge
            environment$strata <- survival::strata
            ordinary <- unname(transform$other_safe)
            if (!body_penalty) ordinary <- c(ordinary, transform$body_safe)
            terms <- c(ordinary, sprintf("ridge(.pry_E, df = %.17g, eps = %.17g, scale = FALSE)", reference$ethnicity_k * 
                ethnic_df_fraction, reference$df_accuracy))
            if (body_penalty) terms <- c(terms, sprintf("ridge(.pry_B, df = %.17g, eps = %.17g, scale = FALSE)", 
                reference$body_k * body_df_fraction, reference$df_accuracy))
            formula <- stats::reformulate(c(terms, "strata(.pry_center)"), response = "survival::Surv(.pry_time, .pry_event)", 
                env = environment)
            control <- do.call(survival::coxph.control, reference$control)
            warnings <- character()
            fit <- withCallingHandlers(survival::coxph(formula, data = work, weights = .pry_weight, ties = "efron", 
                robust = FALSE, singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE, 
                control = control), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
            pry_assert(inherits(fit, "coxph.penal") && inherits(fit, "coxph") && is.null(fit$naive.var) && !length(fit$call$cluster) && 
                identical(fit$method, "efron"), "PRY_NATIVE_PENALIZED_CLASS_AND_NO_CLUSTER_REQUIRED")
            pry_assert(fit$n == s$n && identical(as.numeric(fit$y[, "time"]), as.numeric(time)) && identical(as.numeric(fit$y[, 
                "status"]), as.numeric(event)) && identical(as.numeric(stats::model.weights(fit$model)), as.numeric(weights)), 
                "PRY_DUPLICATED_ROWS_TIME_EVENT_WEIGHT_DRIFT")
            pen <- which(fit$pterms > 0)
            targets <- c(ethnicity = reference$ethnicity_k * ethnic_df_fraction, if (body_penalty) c(body = reference$body_k * 
                body_df_fraction))
            block_k <- c(ethnicity = reference$ethnicity_k, if (body_penalty) c(body = reference$body_k))
            pry_assert(length(pen) == length(targets) && length(fit$history) == length(targets), "PRY_EXACT_SEPARATE_PENALTY_BLOCKS")
            penalty <- lapply(seq_along(pen), function(j) {
                h <- fit$history[[j]]
                history <- as.matrix(h$history)
                theta <- if (nrow(history)) as.numeric(tail(history[, 1L], 1L)) else NA_real_
                achieved <- as.numeric(fit$df[pen[j]])
                list(block = names(targets)[j], k = unname(block_k[j]), target_df = unname(targets[j]), actual_df = achieved, 
                  theta = theta, next_theta = h$theta, outer_done = isTRUE(h$done), history = h, passed = isTRUE(h$done) && 
                    is.finite(theta) && theta >= 0 && is.finite(achieved) && abs(achieved - targets[j]) < reference$df_accuracy)
            })
            names(penalty) <- names(targets)
            finite <- all(is.finite(stats::coef(fit))) && all(is.finite(stats::vcov(fit))) && all(is.finite(fit$linear.predictors)) && 
                all(is.finite(fit$loglik))
            passed <- finite && !length(warnings) && all(vapply(penalty, `[[`, TRUE, "passed")) && fit$iter[1L] < 
                control$outer.max
            anchors <- match(s$observed_centers, s$center)
            mapping <- list(original_columns = reference$original_columns, pooled_columns = reference$pooled_columns, 
                center_columns_removed_for_strata = reference$center_columns, other_safe = transform$other_safe, 
                ethnicity_columns = reference$ethnicity_columns, body_columns = reference$body_columns, native_design_columns = colnames(fit$x), 
                native_coefficient_names = names(stats::coef(fit)), coefficient_alignment = "POSITIONAL_FIT_X_MODEL_MATRIX_ORDER;RIDGE_MATRIX_COEF_NAMES_MAY_REPEAT", 
                original_body_beta = if (body_penalty) "beta_std / B0_sd" else "body_forward %*% beta_QR / B0_sd", 
                lp_origin_shift = "Subtract frozen block means on BOTH fitting/prediction; use native fitted-coordinate baseline")
            list(schema = "SELECTED44_PROTEIN_RIDGE_FIT_PIECE_V1", fit = fit, transform = transform, mapping = mapping, 
                reference = reference, center_levels = reference$center_levels, observed_center_levels = s$observed_centers, 
                center_support = support, reference_X = s$original[anchors, , drop = FALSE], reference_center = s$center[anchors], 
                model_id = model_id, b = as.integer(b), endpoint = endpoint, cause = cause, method_id = if (body_penalty) "PROTEIN_ETHNIC_BODY_RIDGE_CENTER_STRATA_V1" else "PROTEIN_ETHNIC_RIDGE_CENTER_STRATA_QR_V1", 
                sensitivity_setting = sensitivity, primary_gate_eligible = primary_setting, diagnostics = list(passed = passed, 
                  finite = finite, warnings = warnings, penalties = penalty, targets = as.list(targets), sensitivity_setting = sensitivity, 
                  primary_gate_eligible = primary_setting, iterations = fit$iter, control = control, n = s$n, events = sum(event), 
                  weight_sum = sum(weights), original_X_hash = pry_hash(s$original), weights_hash = pry_hash(as.numeric(weights)), 
                  time_status_hash = pry_hash(list(time = time, status = status)), robust = FALSE, cluster_argument_used = FALSE, 
                  duplicates_collapsed = FALSE, inference = "EXTERNAL_FROZEN_FAMILY_BOOTSTRAP_NOT_MODEL_COEFFICIENT_SE", 
                  R = R.version.string, survival = as.character(utils::packageVersion("survival"))))
        }
    }), envir = context)
}

