# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

nmr23591_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

nmr23591_id <- function() ".nmrA_UNTIED_23591_linear"

nmr23591_revision <- function() "NMR23591_LOGIT_PERCENT_A_V2"

nmr23591_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

nmr23591_check_percent <- function(y) {
    nmr23591_assert(is.numeric(y) && length(y) > 0L && all(is.finite(y)), "NMR23591_FINITE_RAW_PERCENT_REQUIRED")
    nmr23591_assert(all(y > 0 & y < 100), paste("NMR23591_ENDPOINT_STATE_MODEL_REQUIRED_NO_EPSILON_NO_REPAIR", 
        "zero_n=", sum(y == 0), "hundred_n=", sum(y == 100), "outside_n=", sum(y < 0 | y > 100), "min=", format(min(y), 
            digits = 17), "max=", format(max(y), digits = 17)))
    invisible(list(N = length(y), min = min(y), max = max(y), zero_n = 0L, hundred_n = 0L, strictly_inside_0_100 = TRUE))
}

nmr23591_conditional_fit <- function(data, y, base, previous, weights, engine, id, lower = -Inf, upper = Inf, rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED") {
    if (!identical(id, nmr23591_id())) 
        return(.nmr23591_original_conditional_fit(data, y, base, previous, weights, engine, id, lower, upper, rule))
    support <- nmr23591_check_percent(y)
    latent <- stats::qlogis(y/100)
    nmr23591_assert(all(is.finite(latent)), "NMR23591_NONFINITE_EXACT_LOGIT_NO_EPSILON")
    f <- .nmr23591_original_conditional_fit(data, latent, base, previous, weights, engine, id, -Inf, Inf, "EMPIRICAL_RESIDUAL_UNTRUNCATED")
    center <- sum(weights * y)/sum(weights)
    scale <- sqrt(sum(weights * (y - center)^2)/sum(weights))
    if (scale == 0) 
        scale <- 1
    f$predictor_center <- center
    f$predictor_scale <- scale
    f$lower <- max(0, lower)
    f$upper <- min(100, upper)
    f$nmr23591_transform <- list(revision = nmr23591_revision(), response = "logit(raw_percent/100)", inverse = "100*plogis(latent)", 
        raw_support = support, original_raw_hash = nmr23591_hash(y), original_raw_predictor_center = center, original_raw_predictor_scale = scale, 
        internal_key = id, unchanged_B_C_H_and_preceding_NMR_predictors = TRUE, subsequent_predictor_is_original_percent = TRUE, 
        observed_design_hash = nmr23591_hash(engine$c39_a_design(data, base, previous)$X), exact_weights_hash = nmr23591_hash(as.numeric(weights)), 
        no_epsilon_no_clipping_no_retry_no_reseed = TRUE)
    f
}

nmr23591_conditional_draw <- function(fit, data, u, engine, audit = NULL, lower = fit$lower, upper = fit$upper) {
    if (!identical(fit$id, nmr23591_id())) 
        return(.nmr23591_original_conditional_draw(fit, data, u, engine, audit, lower, upper))
    tag <- fit$nmr23591_transform
    nmr23591_assert(identical(tag$revision, nmr23591_revision()) && identical(tag$internal_key, fit$id) && identical(fit$predictor_center, 
        tag$original_raw_predictor_center) && identical(fit$predictor_scale, tag$original_raw_predictor_scale), 
        "NMR23591_REQUIRED_LOGIT_TRANSFORM_OR_PERCENT_PREDICTOR_DRIFT")
    X <- engine$c39_a_design(data, NULL, NULL, fit$design_reference)$X
    n <- nrow(X)
    lo <- pmax(0, rep_len(lower, n))
    hi <- pmin(100, rep_len(upper, n))
    nmr23591_assert(all(!is.na(lo) & !is.na(hi) & lo < hi), "NMR23591_INVALID_RAW_SUPPORT_NO_REPAIR")
    latent_lo <- ifelse(lo <= 0, -Inf, stats::qlogis(lo/100))
    latent_hi <- ifelse(hi >= 100, Inf, stats::qlogis(hi/100))
    f <- fit
    f$rule <- if (any(is.finite(latent_lo)) || any(is.finite(latent_hi))) 
        "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
    else "EMPIRICAL_RESIDUAL_UNTRUNCATED"
    z <- engine$c39_residual_draw(f, X, u, latent_lo, latent_hi, FALSE, fit$id, audit)
    value <- 100 * stats::plogis(z)
    nmr23591_assert(all(is.finite(z) & is.finite(value) & value > 0 & value < 100 & value >= lo & value <= hi), 
        "NMR23591_INVERSE_LOGIT_RAW_SUPPORT_FAILURE_NO_CLIPPING_NO_RETRY")
    value
}

nmr23591_install <- function(methods) {
    nmr23591_assert(is.environment(methods) && !identical(methods, .GlobalEnv) && is.function(methods$nmrA_conditional_fit) && 
        is.function(methods$nmrA_conditional_draw), "NMR23591_ISOLATED_ORIGINAL_NMR_METHODS_REQUIRED")
    if (exists(".nmr23591_installed", methods, inherits = FALSE)) {
        nmr23591_assert(identical(methods$.nmr23591_installed, nmr23591_revision()), "NMR23591_INSTALL_REVISION_DRIFT")
        return(invisible(methods))
    }
    methods$.nmr23591_original_conditional_fit <- methods$nmrA_conditional_fit
    methods$.nmr23591_original_conditional_draw <- methods$nmrA_conditional_draw
    src <- environment(nmr23591_install)
    for (name in c("nmr23591_assert", "nmr23591_id", "nmr23591_revision", "nmr23591_hash", "nmr23591_check_percent", 
        "nmr23591_conditional_fit", "nmr23591_conditional_draw", "nmr23591_observed_prefix", "nmr23591_revise_saved_A", 
        "nmr23591_distribution_receipt")) {
        f <- get(name, src, inherits = FALSE)
        environment(f) <- methods
        methods[[name]] <- f
    }
    methods$nmrA_conditional_fit <- methods$nmr23591_conditional_fit
    methods$nmrA_conditional_draw <- methods$nmr23591_conditional_draw
    methods$.nmr23591_installed <- nmr23591_revision()
    invisible(methods)
}

nmr23591_observed_prefix <- function(A, d, methods) {
    original <- d
    raw <- methods$nmrA_raw_values(d, A$members[A$nmr_ids])
    found <- FALSE
    for (step in A$extension) {
        if (identical(step$type, "scalar")) {
            if (identical(step$fit$id, nmr23591_id())) {
                found <- TRUE
                break
            }
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
    nmr23591_assert(found && identical(original, d[names(original)]), "NMR23591_OBSERVED_PREFIX_OR_RAW_DATA_DRIFT")
    d
}

nmr23591_revise_saved_A <- function(A, d, w, engine, methods, provenance = list()) {
    nmr23591_assert(identical(methods$.nmr23591_installed, nmr23591_revision()) && is.data.frame(d) && length(w) == 
        nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)) && A$uniform_dim == 49L && length(A$nmr_ids) == 
        29L, "NMR23591_SAVED_SELECTED29_WEIGHTED_A_REQUIRED")
    target <- which(vapply(A$extension, function(s) identical(s$type, "scalar") && identical(s$fit$id, nmr23591_id()), 
        logical(1)))
    nmr23591_assert(length(target) == 1L, "NMR23591_EXACTLY_ONE_SAVED_SCALAR_REQUIRED")
    step <- A$extension[[target]]
    old <- step$fit
    member <- A$members[["nmr__nmr_23591"]]
    y <- d[[member$column]]
    observed_support <- nmr23591_check_percent(y)
    binding <- nmr23591_hash(list(eid = if ("eid" %in% names(d)) as.character(d$eid), raw = d[vapply(A$members, 
        `[[`, "", "column")], base = d[A$base], weights = as.numeric(w)))
    if (!is.null(old$nmr23591_transform)) {
        nmr23591_assert(identical(A$nmr23591_revision$observed_binding_hash, binding), "NMR23591_REVISED_CACHE_OBSERVED_SAMPLE_WEIGHT_CHANGED")
        return(A)
    }
    before <- methods$nmrA_rng()
    on.exit(nmr23591_assert(identical(before, methods$nmrA_rng()), "NMR23591_REVISION_CHANGED_RNG"), add = TRUE)
    observed <- nmr23591_observed_prefix(A, d, methods)
    X <- engine$c39_a_design(observed, NULL, NULL, old$design_reference)$X
    terms <- attr(stats::terms(old$design_reference$formula), "term.labels")
    nmr23591_assert(all(A$base %in% terms), "NMR23591_ORIGINAL_BASE_PREDICTORS_MISSING")
    previous <- setdiff(terms, A$base)
    new <- methods$nmrA_conditional_fit(observed, y, A$base, previous, w, engine, old$id, old$lower, old$upper, 
        old$rule)
    nmr23591_assert(identical(colnames(X), new$design_reference$columns) && identical(as.numeric(X), as.numeric(engine$c39_a_design(observed, 
        A$base, previous)$X)), "NMR23591_NEW_LATENT_FIT_CONDITIONING_DESIGN_CHANGED")
    new$design_reference <- old$design_reference
    new$predictor_center <- old$predictor_center
    new$predictor_scale <- old$predictor_scale
    new$nmr23591_transform$original_raw_predictor_center <- old$predictor_center
    new$nmr23591_transform$original_raw_predictor_scale <- old$predictor_scale
    revised <- A
    revised$extension[[target]]$fit <- new
    nmr23591_assert(identical(A$extension[-target], revised$extension[-target]) && identical(A$base11, revised$base11) && 
        identical(A$members, revised$members) && identical(A$uniform_columns, revised$uniform_columns) && identical(step$key, 
        revised$extension[[target]]$key), "NMR23591_UNRELATED_FITS_MEMBERS_PREFIX_OR_UNIFORMS_CHANGED")
    revised$nmr23591_revision <- list(revision = nmr23591_revision(), target_id = old$id, target_key = step$key, 
        original_target_fit_hash = nmr23591_hash(old), revised_target_fit_hash = nmr23591_hash(new), observed_binding_hash = binding, 
        observed_support = observed_support, unchanged_other_extension_hash = nmr23591_hash(A$extension[-target]), 
        other_28_raw_fits_not_refitted = TRUE, exact_raw_predictor_center_scale_retained = TRUE, new_fitted_conditionals = 1L, 
        generated_affected_chain_must_rerun = TRUE, observed_raw_Y_columns_unchanged = TRUE, provenance = provenance)
    revised
}

nmr23591_distribution_receipt <- function(diagnostic, A, d, w, methods) {
    observed <- methods$nmr_observed_summary(A, d, w)
    extract <- function(summary) summary$raw[summary$raw$name == "F23591", , drop = FALSE]
    tables <- c(list(observed = extract(observed)), lapply(diagnostic$summaries, extract))
    list(schema = "NMR23591_REAL_DISTRIBUTION_AND_COMPOSITION_REVIEW_V2", observed = TRUE, synthetic = FALSE, revision = nmr23591_revision(), 
        diagnostic_complete = isTRUE(diagnostic$diagnostic_complete), target_percent_weighted_moments_q99_q999 = tables, 
        strict_percent_support = lapply(tables, function(t) if (nrow(t)) all(t$min > 0 & t$max < 100) else FALSE), 
        full_original_NMR_physical_and_composition_pass = isTRUE(diagnostic$physical_pass), composition_relations = lapply(diagnostic$summaries, 
            function(s) s[c("parts", "physical_proportions", "reported_closure", "reported_discrepancy", "physical_failures")]), 
        inverse_transform_does_not_replace_real_distribution_checks = TRUE, no_new_percent_difference_gate = TRUE, 
        tail_changes_are_reported_not_hidden = TRUE, original_full_diagnostic_required = TRUE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        nmr23591_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        nmr23591_id <- function() ".nmrA_UNTIED_23591_linear"
        nmr23591_revision <- function() "NMR23591_LOGIT_PERCENT_A_V2"
        nmr23591_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        nmr23591_check_percent <- function(y) {
            nmr23591_assert(is.numeric(y) && length(y) > 0L && all(is.finite(y)), "NMR23591_FINITE_RAW_PERCENT_REQUIRED")
            nmr23591_assert(all(y > 0 & y < 100), paste("NMR23591_ENDPOINT_STATE_MODEL_REQUIRED_NO_EPSILON_NO_REPAIR", 
                "zero_n=", sum(y == 0), "hundred_n=", sum(y == 100), "outside_n=", sum(y < 0 | y > 100), "min=", 
                format(min(y), digits = 17), "max=", format(max(y), digits = 17)))
            invisible(list(N = length(y), min = min(y), max = max(y), zero_n = 0L, hundred_n = 0L, strictly_inside_0_100 = TRUE))
        }
        nmr23591_conditional_fit <- function(data, y, base, previous, weights, engine, id, lower = -Inf, upper = Inf, 
            rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED") {
            if (!identical(id, nmr23591_id())) return(.nmr23591_original_conditional_fit(data, y, base, previous, 
                weights, engine, id, lower, upper, rule))
            support <- nmr23591_check_percent(y)
            latent <- stats::qlogis(y/100)
            nmr23591_assert(all(is.finite(latent)), "NMR23591_NONFINITE_EXACT_LOGIT_NO_EPSILON")
            f <- .nmr23591_original_conditional_fit(data, latent, base, previous, weights, engine, id, -Inf, Inf, 
                "EMPIRICAL_RESIDUAL_UNTRUNCATED")
            center <- sum(weights * y)/sum(weights)
            scale <- sqrt(sum(weights * (y - center)^2)/sum(weights))
            if (scale == 0) scale <- 1
            f$predictor_center <- center
            f$predictor_scale <- scale
            f$lower <- max(0, lower)
            f$upper <- min(100, upper)
            f$nmr23591_transform <- list(revision = nmr23591_revision(), response = "logit(raw_percent/100)", inverse = "100*plogis(latent)", 
                raw_support = support, original_raw_hash = nmr23591_hash(y), original_raw_predictor_center = center, 
                original_raw_predictor_scale = scale, internal_key = id, unchanged_B_C_H_and_preceding_NMR_predictors = TRUE, 
                subsequent_predictor_is_original_percent = TRUE, observed_design_hash = nmr23591_hash(engine$c39_a_design(data, 
                  base, previous)$X), exact_weights_hash = nmr23591_hash(as.numeric(weights)), no_epsilon_no_clipping_no_retry_no_reseed = TRUE)
            f
        }
        nmr23591_conditional_draw <- function(fit, data, u, engine, audit = NULL, lower = fit$lower, upper = fit$upper) {
            if (!identical(fit$id, nmr23591_id())) return(.nmr23591_original_conditional_draw(fit, data, u, engine, 
                audit, lower, upper))
            tag <- fit$nmr23591_transform
            nmr23591_assert(identical(tag$revision, nmr23591_revision()) && identical(tag$internal_key, fit$id) && 
                identical(fit$predictor_center, tag$original_raw_predictor_center) && identical(fit$predictor_scale, 
                tag$original_raw_predictor_scale), "NMR23591_REQUIRED_LOGIT_TRANSFORM_OR_PERCENT_PREDICTOR_DRIFT")
            X <- engine$c39_a_design(data, NULL, NULL, fit$design_reference)$X
            n <- nrow(X)
            lo <- pmax(0, rep_len(lower, n))
            hi <- pmin(100, rep_len(upper, n))
            nmr23591_assert(all(!is.na(lo) & !is.na(hi) & lo < hi), "NMR23591_INVALID_RAW_SUPPORT_NO_REPAIR")
            latent_lo <- ifelse(lo <= 0, -Inf, stats::qlogis(lo/100))
            latent_hi <- ifelse(hi >= 100, Inf, stats::qlogis(hi/100))
            f <- fit
            f$rule <- if (any(is.finite(latent_lo)) || any(is.finite(latent_hi))) "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL" else "EMPIRICAL_RESIDUAL_UNTRUNCATED"
            z <- engine$c39_residual_draw(f, X, u, latent_lo, latent_hi, FALSE, fit$id, audit)
            value <- 100 * stats::plogis(z)
            nmr23591_assert(all(is.finite(z) & is.finite(value) & value > 0 & value < 100 & value >= lo & value <= 
                hi), "NMR23591_INVERSE_LOGIT_RAW_SUPPORT_FAILURE_NO_CLIPPING_NO_RETRY")
            value
        }
        nmr23591_install <- function(methods) {
            nmr23591_assert(is.environment(methods) && !identical(methods, .GlobalEnv) && is.function(methods$nmrA_conditional_fit) && 
                is.function(methods$nmrA_conditional_draw), "NMR23591_ISOLATED_ORIGINAL_NMR_METHODS_REQUIRED")
            if (exists(".nmr23591_installed", methods, inherits = FALSE)) {
                nmr23591_assert(identical(methods$.nmr23591_installed, nmr23591_revision()), "NMR23591_INSTALL_REVISION_DRIFT")
                return(invisible(methods))
            }
            methods$.nmr23591_original_conditional_fit <- methods$nmrA_conditional_fit
            methods$.nmr23591_original_conditional_draw <- methods$nmrA_conditional_draw
            src <- environment(nmr23591_install)
            for (name in c("nmr23591_assert", "nmr23591_id", "nmr23591_revision", "nmr23591_hash", "nmr23591_check_percent", 
                "nmr23591_conditional_fit", "nmr23591_conditional_draw", "nmr23591_observed_prefix", "nmr23591_revise_saved_A", 
                "nmr23591_distribution_receipt")) {
                f <- get(name, src, inherits = FALSE)
                environment(f) <- methods
                methods[[name]] <- f
            }
            methods$nmrA_conditional_fit <- methods$nmr23591_conditional_fit
            methods$nmrA_conditional_draw <- methods$nmr23591_conditional_draw
            methods$.nmr23591_installed <- nmr23591_revision()
            invisible(methods)
        }
        nmr23591_observed_prefix <- function(A, d, methods) {
            original <- d
            raw <- methods$nmrA_raw_values(d, A$members[A$nmr_ids])
            found <- FALSE
            for (step in A$extension) {
                if (identical(step$type, "scalar")) {
                  if (identical(step$fit$id, nmr23591_id())) {
                    found <- TRUE
                    break
                  }
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
            nmr23591_assert(found && identical(original, d[names(original)]), "NMR23591_OBSERVED_PREFIX_OR_RAW_DATA_DRIFT")
            d
        }
        nmr23591_revise_saved_A <- function(A, d, w, engine, methods, provenance = list()) {
            nmr23591_assert(identical(methods$.nmr23591_installed, nmr23591_revision()) && is.data.frame(d) && 
                length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)) && A$uniform_dim == 49L && 
                length(A$nmr_ids) == 29L, "NMR23591_SAVED_SELECTED29_WEIGHTED_A_REQUIRED")
            target <- which(vapply(A$extension, function(s) identical(s$type, "scalar") && identical(s$fit$id, 
                nmr23591_id()), logical(1)))
            nmr23591_assert(length(target) == 1L, "NMR23591_EXACTLY_ONE_SAVED_SCALAR_REQUIRED")
            step <- A$extension[[target]]
            old <- step$fit
            member <- A$members[["nmr__nmr_23591"]]
            y <- d[[member$column]]
            observed_support <- nmr23591_check_percent(y)
            binding <- nmr23591_hash(list(eid = if ("eid" %in% names(d)) as.character(d$eid), raw = d[vapply(A$members, 
                `[[`, "", "column")], base = d[A$base], weights = as.numeric(w)))
            if (!is.null(old$nmr23591_transform)) {
                nmr23591_assert(identical(A$nmr23591_revision$observed_binding_hash, binding), "NMR23591_REVISED_CACHE_OBSERVED_SAMPLE_WEIGHT_CHANGED")
                return(A)
            }
            before <- methods$nmrA_rng()
            on.exit(nmr23591_assert(identical(before, methods$nmrA_rng()), "NMR23591_REVISION_CHANGED_RNG"), add = TRUE)
            observed <- nmr23591_observed_prefix(A, d, methods)
            X <- engine$c39_a_design(observed, NULL, NULL, old$design_reference)$X
            terms <- attr(stats::terms(old$design_reference$formula), "term.labels")
            nmr23591_assert(all(A$base %in% terms), "NMR23591_ORIGINAL_BASE_PREDICTORS_MISSING")
            previous <- setdiff(terms, A$base)
            new <- methods$nmrA_conditional_fit(observed, y, A$base, previous, w, engine, old$id, old$lower, old$upper, 
                old$rule)
            nmr23591_assert(identical(colnames(X), new$design_reference$columns) && identical(as.numeric(X), as.numeric(engine$c39_a_design(observed, 
                A$base, previous)$X)), "NMR23591_NEW_LATENT_FIT_CONDITIONING_DESIGN_CHANGED")
            new$design_reference <- old$design_reference
            new$predictor_center <- old$predictor_center
            new$predictor_scale <- old$predictor_scale
            new$nmr23591_transform$original_raw_predictor_center <- old$predictor_center
            new$nmr23591_transform$original_raw_predictor_scale <- old$predictor_scale
            revised <- A
            revised$extension[[target]]$fit <- new
            nmr23591_assert(identical(A$extension[-target], revised$extension[-target]) && identical(A$base11, 
                revised$base11) && identical(A$members, revised$members) && identical(A$uniform_columns, revised$uniform_columns) && 
                identical(step$key, revised$extension[[target]]$key), "NMR23591_UNRELATED_FITS_MEMBERS_PREFIX_OR_UNIFORMS_CHANGED")
            revised$nmr23591_revision <- list(revision = nmr23591_revision(), target_id = old$id, target_key = step$key, 
                original_target_fit_hash = nmr23591_hash(old), revised_target_fit_hash = nmr23591_hash(new), observed_binding_hash = binding, 
                observed_support = observed_support, unchanged_other_extension_hash = nmr23591_hash(A$extension[-target]), 
                other_28_raw_fits_not_refitted = TRUE, exact_raw_predictor_center_scale_retained = TRUE, new_fitted_conditionals = 1L, 
                generated_affected_chain_must_rerun = TRUE, observed_raw_Y_columns_unchanged = TRUE, provenance = provenance)
            revised
        }
        nmr23591_distribution_receipt <- function(diagnostic, A, d, w, methods) {
            observed <- methods$nmr_observed_summary(A, d, w)
            extract <- function(summary) summary$raw[summary$raw$name == "F23591", , drop = FALSE]
            tables <- c(list(observed = extract(observed)), lapply(diagnostic$summaries, extract))
            list(schema = "NMR23591_REAL_DISTRIBUTION_AND_COMPOSITION_REVIEW_V2", observed = TRUE, synthetic = FALSE, 
                revision = nmr23591_revision(), diagnostic_complete = isTRUE(diagnostic$diagnostic_complete), target_percent_weighted_moments_q99_q999 = tables, 
                strict_percent_support = lapply(tables, function(t) if (nrow(t)) all(t$min > 0 & t$max < 100) else FALSE), 
                full_original_NMR_physical_and_composition_pass = isTRUE(diagnostic$physical_pass), composition_relations = lapply(diagnostic$summaries, 
                  function(s) s[c("parts", "physical_proportions", "reported_closure", "reported_discrepancy", 
                    "physical_failures")]), inverse_transform_does_not_replace_real_distribution_checks = TRUE, 
                no_new_percent_difference_gate = TRUE, tail_changes_are_reported_not_hidden = TRUE, original_full_diagnostic_required = TRUE)
        }
    }), envir = context)
}

