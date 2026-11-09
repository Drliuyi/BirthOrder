# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

nmr23435_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

nmr23435_id <- function() ".nmrA_UNTIED_23435_linear"

nmr23435_revision <- function() "NMR23435_LOG_RAW_RATIO_A_V3"

nmr23435_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

nmr23435_check_positive <- function(y) {
    nmr23435_assert(is.numeric(y) && length(y) > 0L && all(is.finite(y)), "NMR23435_FINITE_OBSERVED_RAW_REQUIRED")
    nmr23435_assert(all(y > 0), paste("NMR23435_OBSERVED_NONPOSITIVE_REQUIRES_ENDPOINT_STATE_NO_EPSILON", "zero_n=", 
        sum(y == 0), "negative_n=", sum(y < 0), "min=", format(min(y), digits = 17), "max=", format(max(y), digits = 17)))
    invisible(list(N = length(y), min = min(y), max = max(y), zero_n = 0L, negative_n = 0L, strict_positive = TRUE))
}

nmr23435_conditional_fit <- function(data, y, base, previous, weights, engine, id, lower = -Inf, upper = Inf, rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED") {
    if (!identical(id, nmr23435_id())) 
        return(.nmr23435_original_conditional_fit(data, y, base, previous, weights, engine, id, lower, upper, rule))
    support <- nmr23435_check_positive(y)
    nmr23435_assert(length(weights) == length(y) && all(is.finite(weights) & weights > 0) && is.finite(sum(weights)), 
        "NMR23435_EXACT_POSITIVE_WEIGHTS_REQUIRED")
    nmr23435_assert(length(lower) == 1L && length(upper) == 1L && !is.na(lower) && !is.na(upper) && upper > 0 && 
        all(y >= lower & y <= upper), "NMR23435_ORIGINAL_PHYSICAL_SUPPORT_CHANGED")
    latent <- log(y)
    f <- .nmr23435_original_conditional_fit(data, latent, base, previous, weights, engine, id, -Inf, Inf, "EMPIRICAL_RESIDUAL_UNTRUNCATED")
    center <- sum(weights * y)/sum(weights)
    scale <- sqrt(sum(weights * (y - center)^2)/sum(weights))
    if (scale == 0) 
        scale <- 1
    f$predictor_center <- center
    f$predictor_scale <- scale
    f$lower <- max(0, lower)
    f$upper <- upper
    f$nmr23435_transform <- list(revision = nmr23435_revision(), response = "log(original_raw_ratio)", inverse = "exp(latent)", 
        raw_support = support, original_raw_hash = nmr23435_hash(y), exact_weight_hash = nmr23435_hash(as.numeric(weights)), 
        original_raw_predictor_center = center, original_raw_predictor_scale = scale, internal_key = id, original_raw_lower = lower, 
        original_raw_upper = upper, same_B_C_H_and_preceding_NMR_conditioning = TRUE, subsequent_predictor_is_original_raw_ratio = TRUE, 
        observed_design_hash = nmr23435_hash(engine$c39_a_design(data, base, previous)$X), no_epsilon_no_clipping_no_retry_no_reseed = TRUE)
    f
}

nmr23435_conditional_draw <- function(fit, data, u, engine, audit = NULL, lower = fit$lower, upper = fit$upper) {
    if (!identical(fit$id, nmr23435_id())) 
        return(.nmr23435_original_conditional_draw(fit, data, u, engine, audit, lower, upper))
    tag <- fit$nmr23435_transform
    nmr23435_assert(identical(tag$revision, nmr23435_revision()) && identical(tag$internal_key, fit$id) && identical(fit$predictor_center, 
        tag$original_raw_predictor_center) && identical(fit$predictor_scale, tag$original_raw_predictor_scale), 
        "NMR23435_LOG_MARKER_OR_ORIGINAL_RAW_PREDICTOR_DRIFT")
    X <- engine$c39_a_design(data, NULL, NULL, fit$design_reference)$X
    n <- nrow(X)
    lo <- pmax(0, rep_len(lower, n))
    hi <- rep_len(upper, n)
    nmr23435_assert(all(!is.na(lo) & !is.na(hi) & hi > 0 & hi >= lo), "NMR23435_INVALID_RAW_BOUND_NO_REPAIR")
    loglo <- rep(-Inf, n)
    positive <- lo > 0
    loglo[positive] <- log(lo[positive])
    loghi <- log(hi)
    f <- fit
    f$rule <- if (any(is.finite(loglo)) || any(is.finite(loghi))) 
        "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
    else "EMPIRICAL_RESIDUAL_UNTRUNCATED"
    z <- engine$c39_residual_draw(f, X, u, loglo, loghi, FALSE, fit$id, audit)
    value <- exp(z)
    nmr23435_assert(all(is.finite(z) & is.finite(value) & value > 0 & value >= lo & value <= hi), "NMR23435_EXP_RAW_SUPPORT_OR_NUMERICAL_FAILURE_NO_CLIPPING_NO_RETRY")
    value
}

nmr23435_install <- function(methods) {
    nmr23435_assert(is.environment(methods) && !identical(methods, .GlobalEnv) && exists(".nmr23591_installed", 
        methods, inherits = FALSE) && is.function(methods$nmrA_conditional_fit) && is.function(methods$nmrA_conditional_draw), 
        "NMR23435_INSTALL_AFTER_EXISTING_23591_MODULE_IN_ISOLATED_METHODS")
    if (exists(".nmr23435_installed", methods, inherits = FALSE)) {
        nmr23435_assert(identical(methods$.nmr23435_installed, nmr23435_revision()) && identical(methods$nmrA_conditional_fit, 
            methods$nmr23435_conditional_fit) && identical(methods$nmrA_conditional_draw, methods$nmr23435_conditional_draw), 
            "NMR23435_INSTALL_REVISION_DRIFT")
        return(invisible(methods))
    }
    methods$.nmr23435_original_conditional_fit <- methods$nmrA_conditional_fit
    methods$.nmr23435_original_conditional_draw <- methods$nmrA_conditional_draw
    src <- environment(nmr23435_install)
    for (name in c("nmr23435_assert", "nmr23435_id", "nmr23435_revision", "nmr23435_hash", "nmr23435_check_positive", 
        "nmr23435_conditional_fit", "nmr23435_conditional_draw", "nmr23435_observed_prefix", "nmr23435_revise_saved_A", 
        "nmr23435_distribution_receipt")) {
        f <- get(name, src, inherits = FALSE)
        environment(f) <- methods
        methods[[name]] <- f
    }
    methods$nmrA_conditional_fit <- methods$nmr23435_conditional_fit
    methods$nmrA_conditional_draw <- methods$nmr23435_conditional_draw
    methods$.nmr23435_installed <- nmr23435_revision()
    invisible(methods)
}

nmr23435_observed_prefix <- function(A, d, methods, target_step = nmr23435_id()) {
    original <- d
    raw <- methods$nmrA_raw_values(d, A$members[A$nmr_ids])
    if (is.null(target_step)) 
        stop_before <- length(A$extension) + 1L
    else if (is.numeric(target_step)) {
        nmr23435_assert(length(target_step) == 1L && is.finite(target_step) && target_step == floor(target_step) && 
            target_step >= 1L && target_step <= length(A$extension), "NMR23435_ACTUAL_TARGET_STEP_REQUIRED")
        stop_before <- as.integer(target_step)
    }
    else {
        ii <- which(vapply(A$extension, function(s) identical(s$type, "scalar") && identical(s$fit$id, target_step), 
            logical(1)))
        nmr23435_assert(length(ii) == 1L, "NMR23435_UNIQUE_ACTUAL_TARGET_STEP_REQUIRED")
        stop_before <- ii
    }
    for (i in seq_along(A$extension)) {
        if (i == stop_before) 
            break
        step <- A$extension[[i]]
        if (identical(step$type, "scalar")) {
            d[[step$key]] <- methods$nmrA_predictor(step$fit, raw[[paste0("F", step$block$field)]])
        }
        else {
            enc <- methods$nmrA_block_encode(step$block, raw, nrow(d))
            d[[step$mask_key]] <- factor(enc$mask, levels = step$status$levels)
            d[[step$log_key]] <- 0
            for (key in step$ilr_keys) d[[key]] <- 0
            for (pattern in names(step$patterns)) {
                pf <- step$patterns[[pattern]]
                ii <- enc$mask == pattern
                if (!any(ii) || !length(pf$active)) 
                  next
                z <- methods$nmrA_ilr(enc$parts[ii, pf$active, drop = FALSE])
                for (j in seq_along(pf$ilr)) {
                  key <- names(pf$ilr)[j]
                  d[[key]][ii] <- methods$nmrA_predictor(pf$ilr[[key]], z[, j])
                }
                if (!is.null(pf$total)) 
                  d[[step$log_key]][ii] <- methods$nmrA_predictor(pf$total, log(enc$total[ii]))
            }
            for (key in names(step$auxiliary)) {
                af <- step$auxiliary[[key]]
                d[[af$column]] <- methods$nmrA_predictor(af$fit, enc$aux[[key]])
            }
        }
    }
    nmr23435_assert(identical(original, d[names(original)]), "NMR23435_OBSERVED_RAW_COLUMNS_CHANGED")
    d
}

nmr23435_revise_saved_A <- function(A, d, w, engine, methods, provenance = list()) {
    nmr23435_assert(identical(methods$.nmr23435_installed, nmr23435_revision()) && is.data.frame(d) && nrow(d) == 
        65936L && length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)) && A$uniform_dim == 49L && 
        length(A$nmr_ids) == 29L && !is.null(A$base11$body_pmm), "NMR23435_EXACT_SAVED_PMM_NMR_B0_SAMPLE_WEIGHT_WIDTH_REQUIRED")
    target <- which(vapply(A$extension, function(s) identical(s$type, "scalar") && identical(s$fit$id, nmr23435_id()), 
        logical(1)))
    nmr23435_assert(length(target) == 1L, "NMR23435_EXACTLY_ONE_TARGET_SCALAR_REQUIRED")
    step <- A$extension[[target]]
    old <- step$fit
    y <- d[[A$members[["nmr__nmr_23435"]]$column]]
    support <- nmr23435_check_positive(y)
    binding <- nmr23435_hash(list(eid = if ("eid" %in% names(d)) as.character(d$eid), raw = d[vapply(A$members, 
        `[[`, "", "column")], base = d[A$base], weights = as.numeric(w)))
    if (!is.null(old$nmr23435_transform)) {
        nmr23435_assert(identical(A$nmr23435_revision$observed_binding_hash, binding), "NMR23435_REVISED_CACHE_SAMPLE_OR_WEIGHTS_CHANGED")
        return(A)
    }
    before <- methods$nmrA_rng()
    on.exit(nmr23435_assert(identical(before, methods$nmrA_rng()), "NMR23435_ONLY_TARGET_REVISION_CHANGED_RNG"), 
        add = TRUE)
    observed <- nmr23435_observed_prefix(A, d, methods, nmr23435_id())
    X <- engine$c39_a_design(observed, NULL, NULL, old$design_reference)$X
    terms <- attr(stats::terms(old$design_reference$formula), "term.labels")
    nmr23435_assert(all(A$base %in% terms), "NMR23435_ORIGINAL_B_C_CONDITIONING_MISSING")
    previous <- setdiff(terms, A$base)
    new <- methods$nmrA_conditional_fit(observed, y, A$base, previous, w, engine, old$id, old$lower, old$upper, 
        old$rule)
    nmr23435_assert(identical(colnames(X), new$design_reference$columns) && identical(as.numeric(X), as.numeric(engine$c39_a_design(observed, 
        A$base, previous)$X)), "NMR23435_LOG_RESPONSE_FIT_CHANGED_ACTUAL_PRECEDING_DESIGN")
    new$design_reference <- old$design_reference
    new$predictor_center <- old$predictor_center
    new$predictor_scale <- old$predictor_scale
    new$nmr23435_transform$original_raw_predictor_center <- old$predictor_center
    new$nmr23435_transform$original_raw_predictor_scale <- old$predictor_scale
    nmr23435_assert(identical(methods$nmrA_predictor(old, y), methods$nmrA_predictor(new, y)), "NMR23435_SUBSEQUENT_OBSERVED_RAW_RATIO_REPRESENTATION_CHANGED")
    revised <- A
    revised$extension[[target]]$fit <- new
    nmr23435_assert(identical(A$extension[-target], revised$extension[-target]) && identical(A$base11, revised$base11) && 
        identical(A$members, revised$members) && identical(A$uniform_columns, revised$uniform_columns) && identical(step$key, 
        revised$extension[[target]]$key), "NMR23435_OTHER_A_OBJECTS_BODY_PREFIX_MEMBERS_KEY_OR_UNIFORMS_CHANGED")
    rest_before <- A
    rest_after <- revised
    rest_before$extension[[target]]$fit <- NULL
    rest_after$extension[[target]]$fit <- NULL
    nmr23435_assert(identical(rest_before, rest_after), "NMR23435_NON_TARGET_A_OBJECT_CHANGED")
    later <- if (target < length(A$extension)) 
        seq.int(target + 1L, length(A$extension))
    else integer()
    revised$nmr23435_revision <- list(revision = nmr23435_revision(), actual_target_step = target, target_key = step$key, 
        original_target_fit_hash = nmr23435_hash(old), revised_target_fit_hash = nmr23435_hash(new), observed_binding_hash = binding, 
        observed_support = support, other_extension_hash = nmr23435_hash(A$extension[-target]), other_conditionals_not_refitted = TRUE, 
        new_fitted_conditionals = 1L, BODY_prefix_identical = TRUE, original_raw_predictor_center_scale_and_key_identical = TRUE, 
        later_actual_steps = later, later_observed_representation = "EXACT_SAME_STANDARDIZED_RAW_RATIO_NOT_LOG", 
        generated_affected_chain_must_rerun = TRUE, original_raw_Y_columns_unchanged = TRUE, no_change_to_members_N_weights_K49_stream_or_K50_estimand = TRUE, 
        provenance = provenance)
    revised
}

nmr23435_distribution_receipt <- function(diagnostic, A, d, w, methods, tail_review = list()) {
    observed <- methods$nmr_observed_summary(A, d, w)
    extract <- function(s) s$raw[s$raw$name %in% c("F23435", "F23591"), , drop = FALSE]
    summaries <- c(list(observed = extract(observed)), lapply(diagnostic$summaries, extract))
    failed <- !isTRUE(diagnostic$diagnostic_complete) || !isTRUE(diagnostic$physical_pass) || identical(diagnostic$engineering_pass, 
        FALSE) || isTRUE(tail_review$severe_distortion)
    list(schema = "NMR23435_LOG_A_REAL_TAIL_COMPOSITION_RECEIPT_V3", observed = TRUE, synthetic = FALSE, revision = nmr23435_revision(), 
        status = if (failed) "FAILED_REAL_DIAGNOSTIC_OR_SEVERE_TAIL_REVIEW" else "REQUIRES_EXPLICIT_TAIL_REVIEW_NOT_AUTOMATIC_FULL_B0_PASS", 
        weighted_raw_23435_23591_moments_q99_q999 = summaries, physical_and_original_composition_pass = isTRUE(diagnostic$physical_pass), 
        complete_original_NMR_diagnostic = isTRUE(diagnostic$diagnostic_complete), original_existing_engineering_screen = diagnostic$engineering_pass, 
        composition_relations = lapply(diagnostic$summaries, function(s) s[c("parts", "physical_proportions", "reported_closure", 
            "reported_discrepancy", "physical_failures")]), tail_review = tail_review, no_new_percentage_difference_gate = TRUE, 
        inverse_exp_support_not_substitute_for_real_tail_or_composition_checks = TRUE, complete_generated_affected_chain_required = TRUE, 
        B1_to_B100_unchanged_and_not_started_by_this_module = TRUE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        nmr23435_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        nmr23435_id <- function() ".nmrA_UNTIED_23435_linear"
        nmr23435_revision <- function() "NMR23435_LOG_RAW_RATIO_A_V3"
        nmr23435_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        nmr23435_check_positive <- function(y) {
            nmr23435_assert(is.numeric(y) && length(y) > 0L && all(is.finite(y)), "NMR23435_FINITE_OBSERVED_RAW_REQUIRED")
            nmr23435_assert(all(y > 0), paste("NMR23435_OBSERVED_NONPOSITIVE_REQUIRES_ENDPOINT_STATE_NO_EPSILON", 
                "zero_n=", sum(y == 0), "negative_n=", sum(y < 0), "min=", format(min(y), digits = 17), "max=", 
                format(max(y), digits = 17)))
            invisible(list(N = length(y), min = min(y), max = max(y), zero_n = 0L, negative_n = 0L, strict_positive = TRUE))
        }
        nmr23435_conditional_fit <- function(data, y, base, previous, weights, engine, id, lower = -Inf, upper = Inf, 
            rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED") {
            if (!identical(id, nmr23435_id())) return(.nmr23435_original_conditional_fit(data, y, base, previous, 
                weights, engine, id, lower, upper, rule))
            support <- nmr23435_check_positive(y)
            nmr23435_assert(length(weights) == length(y) && all(is.finite(weights) & weights > 0) && is.finite(sum(weights)), 
                "NMR23435_EXACT_POSITIVE_WEIGHTS_REQUIRED")
            nmr23435_assert(length(lower) == 1L && length(upper) == 1L && !is.na(lower) && !is.na(upper) && upper > 
                0 && all(y >= lower & y <= upper), "NMR23435_ORIGINAL_PHYSICAL_SUPPORT_CHANGED")
            latent <- log(y)
            f <- .nmr23435_original_conditional_fit(data, latent, base, previous, weights, engine, id, -Inf, Inf, 
                "EMPIRICAL_RESIDUAL_UNTRUNCATED")
            center <- sum(weights * y)/sum(weights)
            scale <- sqrt(sum(weights * (y - center)^2)/sum(weights))
            if (scale == 0) scale <- 1
            f$predictor_center <- center
            f$predictor_scale <- scale
            f$lower <- max(0, lower)
            f$upper <- upper
            f$nmr23435_transform <- list(revision = nmr23435_revision(), response = "log(original_raw_ratio)", 
                inverse = "exp(latent)", raw_support = support, original_raw_hash = nmr23435_hash(y), exact_weight_hash = nmr23435_hash(as.numeric(weights)), 
                original_raw_predictor_center = center, original_raw_predictor_scale = scale, internal_key = id, 
                original_raw_lower = lower, original_raw_upper = upper, same_B_C_H_and_preceding_NMR_conditioning = TRUE, 
                subsequent_predictor_is_original_raw_ratio = TRUE, observed_design_hash = nmr23435_hash(engine$c39_a_design(data, 
                  base, previous)$X), no_epsilon_no_clipping_no_retry_no_reseed = TRUE)
            f
        }
        nmr23435_conditional_draw <- function(fit, data, u, engine, audit = NULL, lower = fit$lower, upper = fit$upper) {
            if (!identical(fit$id, nmr23435_id())) return(.nmr23435_original_conditional_draw(fit, data, u, engine, 
                audit, lower, upper))
            tag <- fit$nmr23435_transform
            nmr23435_assert(identical(tag$revision, nmr23435_revision()) && identical(tag$internal_key, fit$id) && 
                identical(fit$predictor_center, tag$original_raw_predictor_center) && identical(fit$predictor_scale, 
                tag$original_raw_predictor_scale), "NMR23435_LOG_MARKER_OR_ORIGINAL_RAW_PREDICTOR_DRIFT")
            X <- engine$c39_a_design(data, NULL, NULL, fit$design_reference)$X
            n <- nrow(X)
            lo <- pmax(0, rep_len(lower, n))
            hi <- rep_len(upper, n)
            nmr23435_assert(all(!is.na(lo) & !is.na(hi) & hi > 0 & hi >= lo), "NMR23435_INVALID_RAW_BOUND_NO_REPAIR")
            loglo <- rep(-Inf, n)
            positive <- lo > 0
            loglo[positive] <- log(lo[positive])
            loghi <- log(hi)
            f <- fit
            f$rule <- if (any(is.finite(loglo)) || any(is.finite(loghi))) "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL" else "EMPIRICAL_RESIDUAL_UNTRUNCATED"
            z <- engine$c39_residual_draw(f, X, u, loglo, loghi, FALSE, fit$id, audit)
            value <- exp(z)
            nmr23435_assert(all(is.finite(z) & is.finite(value) & value > 0 & value >= lo & value <= hi), "NMR23435_EXP_RAW_SUPPORT_OR_NUMERICAL_FAILURE_NO_CLIPPING_NO_RETRY")
            value
        }
        nmr23435_install <- function(methods) {
            nmr23435_assert(is.environment(methods) && !identical(methods, .GlobalEnv) && exists(".nmr23591_installed", 
                methods, inherits = FALSE) && is.function(methods$nmrA_conditional_fit) && is.function(methods$nmrA_conditional_draw), 
                "NMR23435_INSTALL_AFTER_EXISTING_23591_MODULE_IN_ISOLATED_METHODS")
            if (exists(".nmr23435_installed", methods, inherits = FALSE)) {
                nmr23435_assert(identical(methods$.nmr23435_installed, nmr23435_revision()) && identical(methods$nmrA_conditional_fit, 
                  methods$nmr23435_conditional_fit) && identical(methods$nmrA_conditional_draw, methods$nmr23435_conditional_draw), 
                  "NMR23435_INSTALL_REVISION_DRIFT")
                return(invisible(methods))
            }
            methods$.nmr23435_original_conditional_fit <- methods$nmrA_conditional_fit
            methods$.nmr23435_original_conditional_draw <- methods$nmrA_conditional_draw
            src <- environment(nmr23435_install)
            for (name in c("nmr23435_assert", "nmr23435_id", "nmr23435_revision", "nmr23435_hash", "nmr23435_check_positive", 
                "nmr23435_conditional_fit", "nmr23435_conditional_draw", "nmr23435_observed_prefix", "nmr23435_revise_saved_A", 
                "nmr23435_distribution_receipt")) {
                f <- get(name, src, inherits = FALSE)
                environment(f) <- methods
                methods[[name]] <- f
            }
            methods$nmrA_conditional_fit <- methods$nmr23435_conditional_fit
            methods$nmrA_conditional_draw <- methods$nmr23435_conditional_draw
            methods$.nmr23435_installed <- nmr23435_revision()
            invisible(methods)
        }
        nmr23435_observed_prefix <- function(A, d, methods, target_step = nmr23435_id()) {
            original <- d
            raw <- methods$nmrA_raw_values(d, A$members[A$nmr_ids])
            if (is.null(target_step)) stop_before <- length(A$extension) + 1L else if (is.numeric(target_step)) {
                nmr23435_assert(length(target_step) == 1L && is.finite(target_step) && target_step == floor(target_step) && 
                  target_step >= 1L && target_step <= length(A$extension), "NMR23435_ACTUAL_TARGET_STEP_REQUIRED")
                stop_before <- as.integer(target_step)
            } else {
                ii <- which(vapply(A$extension, function(s) identical(s$type, "scalar") && identical(s$fit$id, 
                  target_step), logical(1)))
                nmr23435_assert(length(ii) == 1L, "NMR23435_UNIQUE_ACTUAL_TARGET_STEP_REQUIRED")
                stop_before <- ii
            }
            for (i in seq_along(A$extension)) {
                if (i == stop_before) break
                step <- A$extension[[i]]
                if (identical(step$type, "scalar")) {
                  d[[step$key]] <- methods$nmrA_predictor(step$fit, raw[[paste0("F", step$block$field)]])
                } else {
                  enc <- methods$nmrA_block_encode(step$block, raw, nrow(d))
                  d[[step$mask_key]] <- factor(enc$mask, levels = step$status$levels)
                  d[[step$log_key]] <- 0
                  for (key in step$ilr_keys) d[[key]] <- 0
                  for (pattern in names(step$patterns)) {
                    pf <- step$patterns[[pattern]]
                    ii <- enc$mask == pattern
                    if (!any(ii) || !length(pf$active)) next
                    z <- methods$nmrA_ilr(enc$parts[ii, pf$active, drop = FALSE])
                    for (j in seq_along(pf$ilr)) {
                      key <- names(pf$ilr)[j]
                      d[[key]][ii] <- methods$nmrA_predictor(pf$ilr[[key]], z[, j])
                    }
                    if (!is.null(pf$total)) d[[step$log_key]][ii] <- methods$nmrA_predictor(pf$total, log(enc$total[ii]))
                  }
                  for (key in names(step$auxiliary)) {
                    af <- step$auxiliary[[key]]
                    d[[af$column]] <- methods$nmrA_predictor(af$fit, enc$aux[[key]])
                  }
                }
            }
            nmr23435_assert(identical(original, d[names(original)]), "NMR23435_OBSERVED_RAW_COLUMNS_CHANGED")
            d
        }
        nmr23435_revise_saved_A <- function(A, d, w, engine, methods, provenance = list()) {
            nmr23435_assert(identical(methods$.nmr23435_installed, nmr23435_revision()) && is.data.frame(d) && 
                nrow(d) == 65936L && length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)) && 
                A$uniform_dim == 49L && length(A$nmr_ids) == 29L && !is.null(A$base11$body_pmm), "NMR23435_EXACT_SAVED_PMM_NMR_B0_SAMPLE_WEIGHT_WIDTH_REQUIRED")
            target <- which(vapply(A$extension, function(s) identical(s$type, "scalar") && identical(s$fit$id, 
                nmr23435_id()), logical(1)))
            nmr23435_assert(length(target) == 1L, "NMR23435_EXACTLY_ONE_TARGET_SCALAR_REQUIRED")
            step <- A$extension[[target]]
            old <- step$fit
            y <- d[[A$members[["nmr__nmr_23435"]]$column]]
            support <- nmr23435_check_positive(y)
            binding <- nmr23435_hash(list(eid = if ("eid" %in% names(d)) as.character(d$eid), raw = d[vapply(A$members, 
                `[[`, "", "column")], base = d[A$base], weights = as.numeric(w)))
            if (!is.null(old$nmr23435_transform)) {
                nmr23435_assert(identical(A$nmr23435_revision$observed_binding_hash, binding), "NMR23435_REVISED_CACHE_SAMPLE_OR_WEIGHTS_CHANGED")
                return(A)
            }
            before <- methods$nmrA_rng()
            on.exit(nmr23435_assert(identical(before, methods$nmrA_rng()), "NMR23435_ONLY_TARGET_REVISION_CHANGED_RNG"), 
                add = TRUE)
            observed <- nmr23435_observed_prefix(A, d, methods, nmr23435_id())
            X <- engine$c39_a_design(observed, NULL, NULL, old$design_reference)$X
            terms <- attr(stats::terms(old$design_reference$formula), "term.labels")
            nmr23435_assert(all(A$base %in% terms), "NMR23435_ORIGINAL_B_C_CONDITIONING_MISSING")
            previous <- setdiff(terms, A$base)
            new <- methods$nmrA_conditional_fit(observed, y, A$base, previous, w, engine, old$id, old$lower, old$upper, 
                old$rule)
            nmr23435_assert(identical(colnames(X), new$design_reference$columns) && identical(as.numeric(X), as.numeric(engine$c39_a_design(observed, 
                A$base, previous)$X)), "NMR23435_LOG_RESPONSE_FIT_CHANGED_ACTUAL_PRECEDING_DESIGN")
            new$design_reference <- old$design_reference
            new$predictor_center <- old$predictor_center
            new$predictor_scale <- old$predictor_scale
            new$nmr23435_transform$original_raw_predictor_center <- old$predictor_center
            new$nmr23435_transform$original_raw_predictor_scale <- old$predictor_scale
            nmr23435_assert(identical(methods$nmrA_predictor(old, y), methods$nmrA_predictor(new, y)), "NMR23435_SUBSEQUENT_OBSERVED_RAW_RATIO_REPRESENTATION_CHANGED")
            revised <- A
            revised$extension[[target]]$fit <- new
            nmr23435_assert(identical(A$extension[-target], revised$extension[-target]) && identical(A$base11, 
                revised$base11) && identical(A$members, revised$members) && identical(A$uniform_columns, revised$uniform_columns) && 
                identical(step$key, revised$extension[[target]]$key), "NMR23435_OTHER_A_OBJECTS_BODY_PREFIX_MEMBERS_KEY_OR_UNIFORMS_CHANGED")
            rest_before <- A
            rest_after <- revised
            rest_before$extension[[target]]$fit <- NULL
            rest_after$extension[[target]]$fit <- NULL
            nmr23435_assert(identical(rest_before, rest_after), "NMR23435_NON_TARGET_A_OBJECT_CHANGED")
            later <- if (target < length(A$extension)) seq.int(target + 1L, length(A$extension)) else integer()
            revised$nmr23435_revision <- list(revision = nmr23435_revision(), actual_target_step = target, target_key = step$key, 
                original_target_fit_hash = nmr23435_hash(old), revised_target_fit_hash = nmr23435_hash(new), observed_binding_hash = binding, 
                observed_support = support, other_extension_hash = nmr23435_hash(A$extension[-target]), other_conditionals_not_refitted = TRUE, 
                new_fitted_conditionals = 1L, BODY_prefix_identical = TRUE, original_raw_predictor_center_scale_and_key_identical = TRUE, 
                later_actual_steps = later, later_observed_representation = "EXACT_SAME_STANDARDIZED_RAW_RATIO_NOT_LOG", 
                generated_affected_chain_must_rerun = TRUE, original_raw_Y_columns_unchanged = TRUE, no_change_to_members_N_weights_K49_stream_or_K50_estimand = TRUE, 
                provenance = provenance)
            revised
        }
        nmr23435_distribution_receipt <- function(diagnostic, A, d, w, methods, tail_review = list()) {
            observed <- methods$nmr_observed_summary(A, d, w)
            extract <- function(s) s$raw[s$raw$name %in% c("F23435", "F23591"), , drop = FALSE]
            summaries <- c(list(observed = extract(observed)), lapply(diagnostic$summaries, extract))
            failed <- !isTRUE(diagnostic$diagnostic_complete) || !isTRUE(diagnostic$physical_pass) || identical(diagnostic$engineering_pass, 
                FALSE) || isTRUE(tail_review$severe_distortion)
            list(schema = "NMR23435_LOG_A_REAL_TAIL_COMPOSITION_RECEIPT_V3", observed = TRUE, synthetic = FALSE, 
                revision = nmr23435_revision(), status = if (failed) "FAILED_REAL_DIAGNOSTIC_OR_SEVERE_TAIL_REVIEW" else "REQUIRES_EXPLICIT_TAIL_REVIEW_NOT_AUTOMATIC_FULL_B0_PASS", 
                weighted_raw_23435_23591_moments_q99_q999 = summaries, physical_and_original_composition_pass = isTRUE(diagnostic$physical_pass), 
                complete_original_NMR_diagnostic = isTRUE(diagnostic$diagnostic_complete), original_existing_engineering_screen = diagnostic$engineering_pass, 
                composition_relations = lapply(diagnostic$summaries, function(s) s[c("parts", "physical_proportions", 
                  "reported_closure", "reported_discrepancy", "physical_failures")]), tail_review = tail_review, 
                no_new_percentage_difference_gate = TRUE, inverse_exp_support_not_substitute_for_real_tail_or_composition_checks = TRUE, 
                complete_generated_affected_chain_required = TRUE, B1_to_B100_unchanged_and_not_started_by_this_module = TRUE)
        }
    }), envir = context)
}

