# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

p2_policy_spec <- function() paste(c("FOUR_PROTEIN_SELECTIVE_WARNING_RECOVERY_PHASE2_20261031_V1", "scope=fcrl6,gzma,gzmh,klrd1;frozen_actual_draw_and_retained_A_only;no_A_fit;no_new_draw", 
    "source=adopted_current_V4_exact_endpoint_cause_formula_Efron_weights_basis_tau11", "clean_cause=original_native;warned_cause=original_plus_one_native_tighter_warm_retry_only", 
    "allowed_warning=fully_index_mapped_infinite_coefficient_warning_only_for_centerc_numeric_or_ethnicity_nuisance_terms", 
    "warning_support=binary_dummy_observed_in_both_groups_and_zero_cause_events_in_at_least_one_group;otherwise_explicit_unresolved_support_review", 
    "forbidden=core_B_or_protein_warning;unknown_warning;singular;changed_formula;changed_weights;fallback_after_failure", 
    "estimate=original_native_original_family;warm_family=diagnostic_only_no_IE_sign_P_selection", "required=finite_coef_vcov;positive_core_information_and_covariance;max_core_beta_movement_le_0.01_original_SE", 
    "required=all_original_native_time_event_design_weight_center_LP_baseline_CIF_support_guards", "required=legal_observed_and_R00_R10_R11_three_state_probabilities_original_and_warm;frozen_closure_guards", 
    "individual_max_movement_1e-9=DIAGNOSTIC_ONLY_KEEP_FAILED_FLAG_NOT_ACCEPTANCE_NOT_POSTHOC_NEW_CUTOFF", "known_125_singular=explicit_ledger_skip_no_blind_refit;new_native_failure=retain_failure_no_fallback", 
    "old_VALID=never_refit_or_overwrite;all_new_warned=WARNING_QUALIFIED_CONDITIONAL_ESTIMATE", "cache=private_G_prediction_only_compact_original_and_warm_no_full_subject_model_matrices", 
    "full_failure_models=at_most_first_new_failure_per_worker_and_explicit_byte_cap", "restart=skip_only_hash_verified_current_run_sealed_outputs;partial_artifacts_never_overwritten", 
    "claims=no_new_P_q_or_FDR;no_automatic_formal_or_finite_joint_MLE_certificate"), collapse = "\n")

p2_method_review_sha256 <- function() birth_private_setting("PRIVATE_HASH_1")

p2_rebase <- function(path) {
    rw_assert(is.character(path) && length(path) == 1L && nzchar(path), "P2_EXACT_PATH_REQUIRED")
    prefix <- birth_private_setting("PRIVATE_INPUT")
    input <- Sys.getenv("FF_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH"))
    if (startsWith(path, prefix)) 
        file.path(input, substring(path, nchar(prefix) + 1L))
    else path
}

p2_key <- function(id, bid, endpoint) paste(id, as.integer(bid), endpoint, sep = "|")

p2_warning_map <- function(warnings, reference, column) {
    z <- rw_warning_map(warnings, reference, column)
    for (i in seq_along(z$warnings)) {
        q <- z$warnings[[i]]
        terms <- unlist(q$actual_terms, use.names = FALSE)
        q$recognized_center_ethnicity_nuisance_only <- isTRUE(q$eligible_non_target_infinite_warning) && length(terms) > 
            0L && all(grepl("^centerc[0-9]+$|^ethnicity", terms))
        z$warnings[[i]] <- q
    }
    z$qualified <- !length(z$warnings) || all(vapply(z$warnings, function(q) isTRUE(q$recognized_center_ethnicity_nuisance_only), 
        logical(1L)))
    z
}

p2_clone_fit <- function(state, column, private_dir) {
    fn <- rw_clone_fit(state, column, private_dir)
    env <- environment(fn)
    n <- 0L
    visit <- function(z) {
        if (is.call(z)) {
            if (identical(z[[1L]], as.name("if")) && identical(z[[2L]], TRUE)) {
                z[[2L]] <- quote(ff_ce_bad_warning(got$warnings))
                n <<- n + 1L
            }
            for (i in seq_along(z)) if (i != 1L) 
                z[[i]] <- visit(z[[i]])
        }
        z
    }
    body(fn) <- visit(body(fn))
    rw_assert(n == 1L, "P2_SELECTIVE_RETRY_AST_COUNT_CHANGED")
    env$rw_original_eligible <- function(warnings, finite, reference) {
        finite && isTRUE(p2_warning_map(warnings, reference, column)$qualified)
    }
    env$rw_retry_gate <- function(warnings, rec, reference) {
        map <- p2_warning_map(warnings, reference, column)
        rw_assert(isTRUE(map$qualified), "P2_UNRECOGNIZED_OR_CORE_WARM_WARNING_UNRESOLVED")
        state$warning_gates[[length(state$warning_gates) + 1L]] <- list(cause = rec$cause, native_gate = "CE_UNRESOLVED_NATIVE_WARNING_AFTER_SINGLE_RETRY", 
            native_gate_passed = FALSE, warning_mapping = map, disposition = "SEPARATELY_AUTHORIZED_CONDITIONAL_WARNING_POLICY_NOT_NATIVE_GATE_PASS")
        invisible(TRUE)
    }
    env$rw_native <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL) {
        key <- paste(if (is.null(init)) 
            "original"
        else "warm", c("disease", "death")[cause], sep = "_")
        rw_assert(is.null(state$models[[key]]), "P2_ADDITIONAL_NATIVE_RETRY_FORBIDDEN")
        state$fits_attempted <- state$fits_attempted + 1L
        result <- tryCatch(ff_ce_native(work, safe_names, cause, stratified, cluster_bound, control, init), error = function(e) {
            state$native_errors[[key]] <- list(message = conditionMessage(e), class = class(e), singular = grepl("singular", 
                conditionMessage(e), ignore.case = TRUE))
            stop(e)
        })
        state$models[[key]] <- result$value
        state$warnings[[key]] <- result$warnings
        state$controls[[key]] <- control
        result
    }
    state$clone_body_sha256 <- digest::digest(paste(deparse(body(fn), width.cutoff = 500L), collapse = "\n"), algo = "sha256", 
        serialize = FALSE)
    fn
}

p2_diagnostics <- function(state, column, endpoint, w) {
    lapply(stats::setNames(c("disease", "death"), c("disease", "death")), function(cause) {
        dd <- state$designs[[cause]]
        core <- rw_core(dd$reference, column)
        out <- list(cause = cause, core_terms = as.list(dd$reference$columns[core]), core_internal_terms = as.list(dd$reference$safe_names[core]))
        for (family in c("original", "warm")) {
            key <- paste(family, cause, sep = "_")
            m <- state$models[[key]]
            alias <- FALSE
            if (family == "warm" && is.null(m) && !is.null(state$models[[paste0("original_", cause)]]) && !ff_ce_bad_warning(state$warnings[[paste0("original_", 
                cause)]])) {
                m <- state$models[[paste0("original_", cause)]]
                key <- paste0("original_", cause)
                alias <- TRUE
            }
            if (is.null(m)) 
                next
            V <- stats::vcov(m)
            co <- stats::coef(m)
            map <- p2_warning_map(state$warnings[[key]], dd$reference, column)
            wix <- unique(unlist(lapply(map$warnings, function(q) unlist(q$parsed_indices))))
            warning_support <- lapply(wix, function(k) {
                x <- dd$X[, k]
                nz <- x != 0
                event <- endpoint$status == match(cause, c("disease", "death"))
                list(term = dd$reference$columns[k], binary = all(x %in% 0:1), N_nonzero = sum(nz), N_zero = sum(!nz), 
                  cause_events_nonzero = sum(event & nz), cause_events_zero = sum(event & !nz), qualified = all(x %in% 
                    0:1) && any(nz) && any(!nz) && (sum(event & nz) == 0L || sum(event & !nz) == 0L))
            })
            info <- tryCatch(rw_information(m, core), error = function(e) list(positive = FALSE, error = conditionMessage(e)))
            out[[family]] <- list(core_coefficients = as.list(stats::setNames(unname(co[core]), dd$reference$columns[core])), 
                core_SE = as.list(stats::setNames(sqrt(diag(V))[core], dd$reference$columns[core])), core_covariance = unname(V[core, 
                  core, drop = FALSE]), information = info, finite_coefficients_and_covariance = all(is.finite(co)) && 
                  all(is.finite(V)), warning_map = map, warnings = as.list(state$warnings[[key]]), controls = state$controls[[key]], 
                warning_support = warning_support, warning_support_qualified = !length(warning_support) || all(vapply(warning_support, 
                  function(z) isTRUE(z$qualified), logical(1L))), clean_original_reused_without_extra_warm_fit = alias, 
                iterations = m$iter, final_loglik = tail(m$loglik, 1L), core_and_warning_support = rw_support(dd, 
                  endpoint, w, unique(c(core, wix))))
        }
        if (!is.null(out$original) && !is.null(out$warm)) {
            co0 <- unlist(out$original$core_coefficients)
            co1 <- unlist(out$warm$core_coefficients)
            se0 <- unlist(out$original$core_SE)
            delta <- abs(co1 - co0)/se0
            out$core_beta_max_movement_over_original_SE <- max(delta)
            out$core_stability_passed <- all(is.finite(delta)) && all(se0 > 0) && max(delta) <= 0.01
            out$core_positive_information_both <- isTRUE(out$original$information$positive) && isTRUE(out$warm$information$positive)
        }
        out
    })
}

p2_compact_family <- function(object) {
    out <- object
    out$fit <- lapply(object$fit, function(m) {
        keep <- c("coefficients", "var", "naive.var", "means", "method", "n", "nevent", "terms", "xlevels", "contrasts", 
            "assign", "loglik", "iter")
        z <- m[intersect(keep, names(m))]
        z$original_native_fit_means <- z$means
        z$means <- stats::setNames(rep(0, length(z$means)), names(z$means))
        z$prediction_cache_zero_means_encoding_for_reference_zero_only <- TRUE
        e <- new.env(parent = baseenv())
        e$strata <- survival::strata
        if (!is.null(z$terms)) 
            environment(z$terms) <- e
        class(z) <- class(m)
        z
    })
    out$endpoint <- object$endpoint[setdiff(names(object$endpoint), c("time", "status"))]
    out$binding <- list(weight_hash = ff_a_hash(format(object$binding$weight, digits = 17, scientific = TRUE)), 
        cluster = NULL)
    out$numerical_audit <- NULL
    out$numerical_retry <- NULL
    out$PH <- NULL
    out$engine$baseline_public_validation <- lapply(out$engine$baseline_public_validation, function(cause) lapply(cause, 
        function(b) b[setdiff(names(b), c("time", "cumulative"))]))
    out$prediction_only_compact_cache <- TRUE
    out$native_baselines_may_not_be_recomputed_from_stripped_models <- TRUE
    out$compact_cache_supports_only_native_LP_reference_zero <- TRUE
    out$coefficients <- do.call(rbind, lapply(seq_along(out$fit), function(j) data.frame(cause = names(out$fit)[j], 
        term = out$cause_references[[j]]$columns, internal_term = out$cause_references[[j]]$safe_names, coefficient = unname(stats::coef(out$fit[[j]])), 
        SE = sqrt(diag(stats::vcov(out$fit[[j]]))), row.names = NULL)))
    out
}

p2_import_phase1 <- function(item, q, group, policy_sha) {
    fields <- c("phase1_result", "phase1_complete", "phase1_adjudication")
    for (field in fields) rw_check(p2_rebase(item[[paste0(field, "_path")]]), item[[paste0(field, "_sha256")]])
    result <- jsonlite::fromJSON(p2_rebase(item$phase1_result_path), simplifyVector = FALSE)
    seal <- jsonlite::fromJSON(p2_rebase(item$phase1_complete_path), simplifyVector = FALSE)
    adjudication <- jsonlite::fromJSON(p2_rebase(item$phase1_adjudication_path), simplifyVector = FALSE)
    rw_assert(isTRUE(seal$complete) && identical(seal$result_sha256, item$phase1_result_sha256) && identical(seal$case_id, 
        result$case_id), "P2_PHASE1_RESULT_SEAL_CHANGED")
    rw_assert(identical(adjudication$schema, "FOUR_PROTEIN_PHASE1_POSTPROCESS_ADJUDICATION_V1"), "P2_PHASE1_ADJUDICATION_SCHEMA_REQUIRED")
    entry <- Filter(function(z) identical(z$case_id, result$case_id), adjudication$entries)
    rw_assert(length(entry) == 1L && isTRUE(entry[[1L]]$point_recovery_accepted) && isTRUE(entry[[1L]]$postrun_task_A_draw_model_source_hashes_match) && 
        identical(entry[[1L]]$status, "WARNING_QUALIFIED_POINT_CONDITIONAL_RECOVERY"), "P2_PHASE1_ADJUDICATION_NOT_ACCEPTED")
    rw_assert(identical(entry[[1L]]$result_sha256, item$phase1_result_sha256) && identical(entry[[1L]]$complete_sha256, 
        item$phase1_complete_sha256), "P2_PHASE1_ADJUDICATION_FILE_BINDING_CHANGED")
    rw_assert(identical(result$indicator_id, q$indicator_id) && as.integer(result$bootstrap_id) == as.integer(q$bootstrap_id) && 
        identical(result$endpoint_id, q$endpoint_id) && identical(result$original_status, q$original_status) && 
        identical(result$original_reason, q$original_reason), "P2_PHASE1_CASE_OR_ORIGINAL_STATE_CHANGED")
    for (k in c("method_hash", "source_hash", "input_hash", "actual_draw_hash", "target_hash", "sample_hash", "selection_weight_hash")) rw_assert(identical(result[[k]], 
        q[[k]]), paste("P2_PHASE1_FROZEN_BINDING_CHANGED", k))
    rw_assert(identical(result$draw_private_sha256, group$draw_rds_sha256) && identical(result$A_descriptor_sha256, 
        group$a_descriptors_sha256), "P2_PHASE1_A_OR_DRAW_SOURCE_CHANGED")
    rw_assert(length(seal$private_models) == 4L, "P2_PHASE1_FOUR_RETAINED_NATIVE_MODELS_REQUIRED")
    for (model in seal$private_models) rw_check(p2_rebase(model$path), model$sha256)
    rw_assert(isTRUE(result$native_validation$original$passed) && isTRUE(result$native_validation$warm$passed), 
        "P2_PHASE1_NATIVE_VALIDATION_FAILED")
    diagnostics <- result$coefficient_diagnostics
    rw_assert(length(diagnostics) == 2L && all(vapply(diagnostics, function(d) isTRUE(d$core_stability_descriptive_rule_passed) && 
        isTRUE(d$core_positive_information_both), logical(1L))), "P2_PHASE1_CORE_RULE_NOT_MET")
    for (d in diagnostics) for (family in c("original", "warm")) {
        map <- d[[family]]$warning_map
        for (warning in map$warnings) {
            terms <- unlist(warning$actual_terms, use.names = FALSE)
            rw_assert(isTRUE(warning$eligible_non_target_infinite_warning) && length(terms) > 0L && all(grepl("^centerc[0-9]+$|^ethnicity", 
                terms)), "P2_PHASE1_WARNING_CLASS_UNRECOGNIZED")
            for (term in terms) {
                support <- Filter(function(s) identical(s$term, term), d[[family]]$core_and_warning_event_weight_support)
                rw_assert(length(support) == 1L, "P2_PHASE1_WARNING_SUPPORT_NOT_RETAINED")
                s <- support[[1L]]
                event_prefix <- if (identical(d$cause, "disease")) 
                  "disease"
                else "death"
                rw_assert(s$minimum %in% 0:1 && s$maximum %in% 0:1 && s$N_nonzero > 0 && s$N_zero > 0 && (s[[paste0(event_prefix, 
                  "_events_nonzero")]] == 0 || s[[paste0(event_prefix, "_events_zero")]] == 0), "P2_PHASE1_WARNING_SUPPORT_REQUIRES_SEPARATE_REVIEW")
            }
        }
    }
    rw_assert(length(result$predictions) == 4L && all(vapply(result$predictions, function(p) isTRUE(p$original$legal) && 
        isTRUE(p$warm$legal), logical(1L))), "P2_PHASE1_PROBABILITY_LEGALITY_FAILED")
    q$original_effects <- result$original_effects
    q$warm_effects <- result$warm_effects
    for (k in c("R00", "R10", "R11", "TE", "DE", "IE", "PM")) q[[k]] <- q$original_effects[[k]]
    q$coefficient_diagnostics <- diagnostics
    q$predictions <- result$predictions
    q$native_validation <- result$native_validation
    q$event_weight_support <- result$event_weight_support
    q$individual_probability_1e9_diagnostic_passed <- all(vapply(result$predictions, function(p) isTRUE(p$individual_probability_stability_passed), 
        logical(1L)))
    q$individual_probability_1e9_is_acceptance_gate <- FALSE
    q$individual_probability_1e9_failed_flag_was_not_renamed_or_relaxed <- TRUE
    q$status <- "WARNING_QUALIFIED_CONDITIONAL_ESTIMATE"
    q$reason <- "AUTHENTICATED_ADJUDICATED_PHASE1_REUSE_NO_Y_REFIT_OR_REPREDICTION"
    q$estimate_family <- "NATIVE_ORIGINAL_FAMILY_NO_IE_SIGN_P_SELECTION"
    q$phase1_reused_without_Y_refit <- TRUE
    q$phase1_source <- item[grep("^phase1_", names(item))]
    q$platform_fit_binding_difference_not_claimed_bitwise <- TRUE
    q$phase1_native_warning_gate_exceptions <- result$native_warning_gate_exceptions
    cache <- list(schema = "P2_PRIVATE_PHASE1_REFERENCE_POSTPROCESS_CACHE_V1", recovery_policy_sha256 = policy_sha, 
        phase1_source = q$phase1_source, private_models = seal$private_models, original_effects = q$original_effects, 
        warm_effects = q$warm_effects, predictions = q$predictions, diagnostics = q$coefficient_diagnostics, Y_refits_or_repredictions = 0L, 
        task_seal_sha256 = group$task_seal_sha256, draw_rds_sha256 = group$draw_rds_sha256, A_descriptors_sha256 = group$a_descriptors_sha256)
    list(result = q, cache = cache)
}

p2_verified <- function(directory, identity) {
    sp <- file.path(directory, "COMPLETE.json")
    rw_assert(file.exists(sp), paste("P2_PARTIAL_ARTIFACT_PRESERVED_NO_AUTOMATIC_OVERWRITE", directory))
    seal <- jsonlite::fromJSON(sp, simplifyVector = FALSE)
    rw_assert(isTRUE(seal$complete) && all(vapply(names(identity), function(k) identical(seal[[k]], identity[[k]]), 
        logical(1L))), "P2_EXISTING_SEAL_IDENTITY_CHANGED")
    for (name in names(seal$artifacts_sha256)) rw_check(file.path(directory, name), seal$artifacts_sha256[[name]])
    jsonlite::fromJSON(file.path(directory, "RESULT.json"), simplifyVector = FALSE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        p2_policy_spec <- function() paste(c("FOUR_PROTEIN_SELECTIVE_WARNING_RECOVERY_PHASE2_20261031_V1", "scope=fcrl6,gzma,gzmh,klrd1;frozen_actual_draw_and_retained_A_only;no_A_fit;no_new_draw", 
            "source=adopted_current_V4_exact_endpoint_cause_formula_Efron_weights_basis_tau11", "clean_cause=original_native;warned_cause=original_plus_one_native_tighter_warm_retry_only", 
            "allowed_warning=fully_index_mapped_infinite_coefficient_warning_only_for_centerc_numeric_or_ethnicity_nuisance_terms", 
            "warning_support=binary_dummy_observed_in_both_groups_and_zero_cause_events_in_at_least_one_group;otherwise_explicit_unresolved_support_review", 
            "forbidden=core_B_or_protein_warning;unknown_warning;singular;changed_formula;changed_weights;fallback_after_failure", 
            "estimate=original_native_original_family;warm_family=diagnostic_only_no_IE_sign_P_selection", "required=finite_coef_vcov;positive_core_information_and_covariance;max_core_beta_movement_le_0.01_original_SE", 
            "required=all_original_native_time_event_design_weight_center_LP_baseline_CIF_support_guards", "required=legal_observed_and_R00_R10_R11_three_state_probabilities_original_and_warm;frozen_closure_guards", 
            "individual_max_movement_1e-9=DIAGNOSTIC_ONLY_KEEP_FAILED_FLAG_NOT_ACCEPTANCE_NOT_POSTHOC_NEW_CUTOFF", 
            "known_125_singular=explicit_ledger_skip_no_blind_refit;new_native_failure=retain_failure_no_fallback", 
            "old_VALID=never_refit_or_overwrite;all_new_warned=WARNING_QUALIFIED_CONDITIONAL_ESTIMATE", "cache=private_G_prediction_only_compact_original_and_warm_no_full_subject_model_matrices", 
            "full_failure_models=at_most_first_new_failure_per_worker_and_explicit_byte_cap", "restart=skip_only_hash_verified_current_run_sealed_outputs;partial_artifacts_never_overwritten", 
            "claims=no_new_P_q_or_FDR;no_automatic_formal_or_finite_joint_MLE_certificate"), collapse = "\n")
        p2_method_review_sha256 <- function() birth_private_setting("PRIVATE_HASH_4")
        p2_rebase <- function(path) {
            rw_assert(is.character(path) && length(path) == 1L && nzchar(path), "P2_EXACT_PATH_REQUIRED")
            prefix <- birth_private_setting("PRIVATE_INPUT")
            input <- Sys.getenv("FF_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH"))
            if (startsWith(path, prefix)) file.path(input, substring(path, nchar(prefix) + 1L)) else path
        }
        p2_key <- function(id, bid, endpoint) paste(id, as.integer(bid), endpoint, sep = "|")
        p2_warning_map <- function(warnings, reference, column) {
            z <- rw_warning_map(warnings, reference, column)
            for (i in seq_along(z$warnings)) {
                q <- z$warnings[[i]]
                terms <- unlist(q$actual_terms, use.names = FALSE)
                q$recognized_center_ethnicity_nuisance_only <- isTRUE(q$eligible_non_target_infinite_warning) && 
                  length(terms) > 0L && all(grepl("^centerc[0-9]+$|^ethnicity", terms))
                z$warnings[[i]] <- q
            }
            z$qualified <- !length(z$warnings) || all(vapply(z$warnings, function(q) isTRUE(q$recognized_center_ethnicity_nuisance_only), 
                logical(1L)))
            z
        }
        p2_clone_fit <- function(state, column, private_dir) {
            fn <- rw_clone_fit(state, column, private_dir)
            env <- environment(fn)
            n <- 0L
            visit <- function(z) {
                if (is.call(z)) {
                  if (identical(z[[1L]], as.name("if")) && identical(z[[2L]], TRUE)) {
                    z[[2L]] <- quote(ff_ce_bad_warning(got$warnings))
                    n <<- n + 1L
                  }
                  for (i in seq_along(z)) if (i != 1L) z[[i]] <- visit(z[[i]])
                }
                z
            }
            body(fn) <- visit(body(fn))
            rw_assert(n == 1L, "P2_SELECTIVE_RETRY_AST_COUNT_CHANGED")
            env$rw_original_eligible <- function(warnings, finite, reference) {
                finite && isTRUE(p2_warning_map(warnings, reference, column)$qualified)
            }
            env$rw_retry_gate <- function(warnings, rec, reference) {
                map <- p2_warning_map(warnings, reference, column)
                rw_assert(isTRUE(map$qualified), "P2_UNRECOGNIZED_OR_CORE_WARM_WARNING_UNRESOLVED")
                state$warning_gates[[length(state$warning_gates) + 1L]] <- list(cause = rec$cause, native_gate = "CE_UNRESOLVED_NATIVE_WARNING_AFTER_SINGLE_RETRY", 
                  native_gate_passed = FALSE, warning_mapping = map, disposition = "SEPARATELY_AUTHORIZED_CONDITIONAL_WARNING_POLICY_NOT_NATIVE_GATE_PASS")
                invisible(TRUE)
            }
            env$rw_native <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL) {
                key <- paste(if (is.null(init)) "original" else "warm", c("disease", "death")[cause], sep = "_")
                rw_assert(is.null(state$models[[key]]), "P2_ADDITIONAL_NATIVE_RETRY_FORBIDDEN")
                state$fits_attempted <- state$fits_attempted + 1L
                result <- tryCatch(ff_ce_native(work, safe_names, cause, stratified, cluster_bound, control, init), 
                  error = function(e) {
                    state$native_errors[[key]] <- list(message = conditionMessage(e), class = class(e), singular = grepl("singular", 
                      conditionMessage(e), ignore.case = TRUE))
                    stop(e)
                  })
                state$models[[key]] <- result$value
                state$warnings[[key]] <- result$warnings
                state$controls[[key]] <- control
                result
            }
            state$clone_body_sha256 <- digest::digest(paste(deparse(body(fn), width.cutoff = 500L), collapse = "\n"), 
                algo = "sha256", serialize = FALSE)
            fn
        }
        p2_diagnostics <- function(state, column, endpoint, w) {
            lapply(stats::setNames(c("disease", "death"), c("disease", "death")), function(cause) {
                dd <- state$designs[[cause]]
                core <- rw_core(dd$reference, column)
                out <- list(cause = cause, core_terms = as.list(dd$reference$columns[core]), core_internal_terms = as.list(dd$reference$safe_names[core]))
                for (family in c("original", "warm")) {
                  key <- paste(family, cause, sep = "_")
                  m <- state$models[[key]]
                  alias <- FALSE
                  if (family == "warm" && is.null(m) && !is.null(state$models[[paste0("original_", cause)]]) && 
                    !ff_ce_bad_warning(state$warnings[[paste0("original_", cause)]])) {
                    m <- state$models[[paste0("original_", cause)]]
                    key <- paste0("original_", cause)
                    alias <- TRUE
                  }
                  if (is.null(m)) next
                  V <- stats::vcov(m)
                  co <- stats::coef(m)
                  map <- p2_warning_map(state$warnings[[key]], dd$reference, column)
                  wix <- unique(unlist(lapply(map$warnings, function(q) unlist(q$parsed_indices))))
                  warning_support <- lapply(wix, function(k) {
                    x <- dd$X[, k]
                    nz <- x != 0
                    event <- endpoint$status == match(cause, c("disease", "death"))
                    list(term = dd$reference$columns[k], binary = all(x %in% 0:1), N_nonzero = sum(nz), N_zero = sum(!nz), 
                      cause_events_nonzero = sum(event & nz), cause_events_zero = sum(event & !nz), qualified = all(x %in% 
                        0:1) && any(nz) && any(!nz) && (sum(event & nz) == 0L || sum(event & !nz) == 0L))
                  })
                  info <- tryCatch(rw_information(m, core), error = function(e) list(positive = FALSE, error = conditionMessage(e)))
                  out[[family]] <- list(core_coefficients = as.list(stats::setNames(unname(co[core]), dd$reference$columns[core])), 
                    core_SE = as.list(stats::setNames(sqrt(diag(V))[core], dd$reference$columns[core])), core_covariance = unname(V[core, 
                      core, drop = FALSE]), information = info, finite_coefficients_and_covariance = all(is.finite(co)) && 
                      all(is.finite(V)), warning_map = map, warnings = as.list(state$warnings[[key]]), controls = state$controls[[key]], 
                    warning_support = warning_support, warning_support_qualified = !length(warning_support) || 
                      all(vapply(warning_support, function(z) isTRUE(z$qualified), logical(1L))), clean_original_reused_without_extra_warm_fit = alias, 
                    iterations = m$iter, final_loglik = tail(m$loglik, 1L), core_and_warning_support = rw_support(dd, 
                      endpoint, w, unique(c(core, wix))))
                }
                if (!is.null(out$original) && !is.null(out$warm)) {
                  co0 <- unlist(out$original$core_coefficients)
                  co1 <- unlist(out$warm$core_coefficients)
                  se0 <- unlist(out$original$core_SE)
                  delta <- abs(co1 - co0)/se0
                  out$core_beta_max_movement_over_original_SE <- max(delta)
                  out$core_stability_passed <- all(is.finite(delta)) && all(se0 > 0) && max(delta) <= 0.01
                  out$core_positive_information_both <- isTRUE(out$original$information$positive) && isTRUE(out$warm$information$positive)
                }
                out
            })
        }
        p2_compact_family <- function(object) {
            out <- object
            out$fit <- lapply(object$fit, function(m) {
                keep <- c("coefficients", "var", "naive.var", "means", "method", "n", "nevent", "terms", "xlevels", 
                  "contrasts", "assign", "loglik", "iter")
                z <- m[intersect(keep, names(m))]
                z$original_native_fit_means <- z$means
                z$means <- stats::setNames(rep(0, length(z$means)), names(z$means))
                z$prediction_cache_zero_means_encoding_for_reference_zero_only <- TRUE
                e <- new.env(parent = baseenv())
                e$strata <- survival::strata
                if (!is.null(z$terms)) environment(z$terms) <- e
                class(z) <- class(m)
                z
            })
            out$endpoint <- object$endpoint[setdiff(names(object$endpoint), c("time", "status"))]
            out$binding <- list(weight_hash = ff_a_hash(format(object$binding$weight, digits = 17, scientific = TRUE)), 
                cluster = NULL)
            out$numerical_audit <- NULL
            out$numerical_retry <- NULL
            out$PH <- NULL
            out$engine$baseline_public_validation <- lapply(out$engine$baseline_public_validation, function(cause) lapply(cause, 
                function(b) b[setdiff(names(b), c("time", "cumulative"))]))
            out$prediction_only_compact_cache <- TRUE
            out$native_baselines_may_not_be_recomputed_from_stripped_models <- TRUE
            out$compact_cache_supports_only_native_LP_reference_zero <- TRUE
            out$coefficients <- do.call(rbind, lapply(seq_along(out$fit), function(j) data.frame(cause = names(out$fit)[j], 
                term = out$cause_references[[j]]$columns, internal_term = out$cause_references[[j]]$safe_names, 
                coefficient = unname(stats::coef(out$fit[[j]])), SE = sqrt(diag(stats::vcov(out$fit[[j]]))), row.names = NULL)))
            out
        }
        p2_import_phase1 <- function(item, q, group, policy_sha) {
            fields <- c("phase1_result", "phase1_complete", "phase1_adjudication")
            for (field in fields) rw_check(p2_rebase(item[[paste0(field, "_path")]]), item[[paste0(field, "_sha256")]])
            result <- jsonlite::fromJSON(p2_rebase(item$phase1_result_path), simplifyVector = FALSE)
            seal <- jsonlite::fromJSON(p2_rebase(item$phase1_complete_path), simplifyVector = FALSE)
            adjudication <- jsonlite::fromJSON(p2_rebase(item$phase1_adjudication_path), simplifyVector = FALSE)
            rw_assert(isTRUE(seal$complete) && identical(seal$result_sha256, item$phase1_result_sha256) && identical(seal$case_id, 
                result$case_id), "P2_PHASE1_RESULT_SEAL_CHANGED")
            rw_assert(identical(adjudication$schema, "FOUR_PROTEIN_PHASE1_POSTPROCESS_ADJUDICATION_V1"), "P2_PHASE1_ADJUDICATION_SCHEMA_REQUIRED")
            entry <- Filter(function(z) identical(z$case_id, result$case_id), adjudication$entries)
            rw_assert(length(entry) == 1L && isTRUE(entry[[1L]]$point_recovery_accepted) && isTRUE(entry[[1L]]$postrun_task_A_draw_model_source_hashes_match) && 
                identical(entry[[1L]]$status, "WARNING_QUALIFIED_POINT_CONDITIONAL_RECOVERY"), "P2_PHASE1_ADJUDICATION_NOT_ACCEPTED")
            rw_assert(identical(entry[[1L]]$result_sha256, item$phase1_result_sha256) && identical(entry[[1L]]$complete_sha256, 
                item$phase1_complete_sha256), "P2_PHASE1_ADJUDICATION_FILE_BINDING_CHANGED")
            rw_assert(identical(result$indicator_id, q$indicator_id) && as.integer(result$bootstrap_id) == as.integer(q$bootstrap_id) && 
                identical(result$endpoint_id, q$endpoint_id) && identical(result$original_status, q$original_status) && 
                identical(result$original_reason, q$original_reason), "P2_PHASE1_CASE_OR_ORIGINAL_STATE_CHANGED")
            for (k in c("method_hash", "source_hash", "input_hash", "actual_draw_hash", "target_hash", "sample_hash", 
                "selection_weight_hash")) rw_assert(identical(result[[k]], q[[k]]), paste("P2_PHASE1_FROZEN_BINDING_CHANGED", 
                k))
            rw_assert(identical(result$draw_private_sha256, group$draw_rds_sha256) && identical(result$A_descriptor_sha256, 
                group$a_descriptors_sha256), "P2_PHASE1_A_OR_DRAW_SOURCE_CHANGED")
            rw_assert(length(seal$private_models) == 4L, "P2_PHASE1_FOUR_RETAINED_NATIVE_MODELS_REQUIRED")
            for (model in seal$private_models) rw_check(p2_rebase(model$path), model$sha256)
            rw_assert(isTRUE(result$native_validation$original$passed) && isTRUE(result$native_validation$warm$passed), 
                "P2_PHASE1_NATIVE_VALIDATION_FAILED")
            diagnostics <- result$coefficient_diagnostics
            rw_assert(length(diagnostics) == 2L && all(vapply(diagnostics, function(d) isTRUE(d$core_stability_descriptive_rule_passed) && 
                isTRUE(d$core_positive_information_both), logical(1L))), "P2_PHASE1_CORE_RULE_NOT_MET")
            for (d in diagnostics) for (family in c("original", "warm")) {
                map <- d[[family]]$warning_map
                for (warning in map$warnings) {
                  terms <- unlist(warning$actual_terms, use.names = FALSE)
                  rw_assert(isTRUE(warning$eligible_non_target_infinite_warning) && length(terms) > 0L && all(grepl("^centerc[0-9]+$|^ethnicity", 
                    terms)), "P2_PHASE1_WARNING_CLASS_UNRECOGNIZED")
                  for (term in terms) {
                    support <- Filter(function(s) identical(s$term, term), d[[family]]$core_and_warning_event_weight_support)
                    rw_assert(length(support) == 1L, "P2_PHASE1_WARNING_SUPPORT_NOT_RETAINED")
                    s <- support[[1L]]
                    event_prefix <- if (identical(d$cause, "disease")) "disease" else "death"
                    rw_assert(s$minimum %in% 0:1 && s$maximum %in% 0:1 && s$N_nonzero > 0 && s$N_zero > 0 && (s[[paste0(event_prefix, 
                      "_events_nonzero")]] == 0 || s[[paste0(event_prefix, "_events_zero")]] == 0), "P2_PHASE1_WARNING_SUPPORT_REQUIRES_SEPARATE_REVIEW")
                  }
                }
            }
            rw_assert(length(result$predictions) == 4L && all(vapply(result$predictions, function(p) isTRUE(p$original$legal) && 
                isTRUE(p$warm$legal), logical(1L))), "P2_PHASE1_PROBABILITY_LEGALITY_FAILED")
            q$original_effects <- result$original_effects
            q$warm_effects <- result$warm_effects
            for (k in c("R00", "R10", "R11", "TE", "DE", "IE", "PM")) q[[k]] <- q$original_effects[[k]]
            q$coefficient_diagnostics <- diagnostics
            q$predictions <- result$predictions
            q$native_validation <- result$native_validation
            q$event_weight_support <- result$event_weight_support
            q$individual_probability_1e9_diagnostic_passed <- all(vapply(result$predictions, function(p) isTRUE(p$individual_probability_stability_passed), 
                logical(1L)))
            q$individual_probability_1e9_is_acceptance_gate <- FALSE
            q$individual_probability_1e9_failed_flag_was_not_renamed_or_relaxed <- TRUE
            q$status <- "WARNING_QUALIFIED_CONDITIONAL_ESTIMATE"
            q$reason <- "AUTHENTICATED_ADJUDICATED_PHASE1_REUSE_NO_Y_REFIT_OR_REPREDICTION"
            q$estimate_family <- "NATIVE_ORIGINAL_FAMILY_NO_IE_SIGN_P_SELECTION"
            q$phase1_reused_without_Y_refit <- TRUE
            q$phase1_source <- item[grep("^phase1_", names(item))]
            q$platform_fit_binding_difference_not_claimed_bitwise <- TRUE
            q$phase1_native_warning_gate_exceptions <- result$native_warning_gate_exceptions
            cache <- list(schema = "P2_PRIVATE_PHASE1_REFERENCE_POSTPROCESS_CACHE_V1", recovery_policy_sha256 = policy_sha, 
                phase1_source = q$phase1_source, private_models = seal$private_models, original_effects = q$original_effects, 
                warm_effects = q$warm_effects, predictions = q$predictions, diagnostics = q$coefficient_diagnostics, 
                Y_refits_or_repredictions = 0L, task_seal_sha256 = group$task_seal_sha256, draw_rds_sha256 = group$draw_rds_sha256, 
                A_descriptors_sha256 = group$a_descriptors_sha256)
            list(result = q, cache = cache)
        }
        p2_verified <- function(directory, identity) {
            sp <- file.path(directory, "COMPLETE.json")
            rw_assert(file.exists(sp), paste("P2_PARTIAL_ARTIFACT_PRESERVED_NO_AUTOMATIC_OVERWRITE", directory))
            seal <- jsonlite::fromJSON(sp, simplifyVector = FALSE)
            rw_assert(isTRUE(seal$complete) && all(vapply(names(identity), function(k) identical(seal[[k]], identity[[k]]), 
                logical(1L))), "P2_EXISTING_SEAL_IDENTITY_CHANGED")
            for (name in names(seal$artifacts_sha256)) rw_check(file.path(directory, name), seal$artifacts_sha256[[name]])
            jsonlite::fromJSON(file.path(directory, "RESULT.json"), simplifyVector = FALSE)
        }
    }), envir = context)
}

