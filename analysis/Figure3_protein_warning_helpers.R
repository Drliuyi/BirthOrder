# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

rw_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

rw_json <- function(value, path) {
    rw_assert(!file.exists(path), paste("PRESERVE_EXISTING_RECOVERY_ARTIFACT", path))
    jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
}

rw_sha <- function(path) {
    rw_assert(file.exists(path) && !dir.exists(path), paste("RECOVERY_FILE_MISSING", path))
    digest::digest(file = path, algo = "sha256")
}

rw_check <- function(path, sha) {
    rw_assert(is.character(sha) && length(sha) == 1L && grepl("^[a-f0-9]{64}$", sha), "RECOVERY_EXACT_SHA256_REQUIRED")
    rw_assert(identical(rw_sha(path), sha), paste("RECOVERY_SOURCE_SHA256_CHANGED", path))
}

rw_inside <- function(path, parent) {
    path <- normalizePath(path, mustWork = TRUE)
    parent <- normalizePath(parent, mustWork = TRUE)
    startsWith(paste0(path, "/"), paste0(parent, "/"))
}

rw_newdir <- function(path) {
    rw_assert(!file.exists(path) && dir.create(path, recursive = TRUE, showWarnings = FALSE), paste("RECOVERY_NEW_ISOLATED_DIRECTORY_REQUIRED", 
        path))
    normalizePath(path, mustWork = TRUE)
}

rw_core <- function(reference, column) {
    ix <- which(reference$columns == "B" | grepl(column, reference$columns, fixed = TRUE))
    rw_assert(sum(reference$columns[ix] == "B") == 1L && length(ix) >= 2L, "RECOVERY_B_AND_PROTEIN_SPLINE_CORE_MAPPING_REQUIRED")
    ix
}

rw_warning_map <- function(warnings, reference, column) {
    core <- rw_core(reference, column)
    rows <- lapply(warnings, function(w) {
        pattern <- "(?i)variable[s]?[[:space:]]+([0-9]+(?:[[:space:],]+[0-9]+)*)"
        hit <- regmatches(w, regexec(pattern, w, perl = TRUE))[[1L]]
        indices <- if (length(hit) >= 2L) 
            as.integer(regmatches(hit[2L], gregexpr("[0-9]+", hit[2L], perl = TRUE))[[1L]])
        else integer()
        eligible_text <- grepl("^Loglik converged before variable.*coefficient may be infinite[.]?[[:space:]]*$", 
            trimws(w), ignore.case = TRUE)
        mapped <- length(indices) > 0L && !anyNA(indices) && all(indices %in% seq_along(reference$columns))
        list(message = w, parsed_indices = as.list(indices), completely_mapped = mapped, actual_terms = if (mapped) as.list(reference$columns[indices]) else list(), 
            internal_terms = if (mapped) as.list(reference$safe_names[indices]) else list(), core_warning = if (mapped) any(indices %in% 
                core) else NA, eligible_non_target_infinite_warning = eligible_text && mapped && !any(indices %in% 
                core))
    })
    list(warnings = rows, all_non_target_infinite = !length(rows) || all(vapply(rows, function(z) isTRUE(z$eligible_non_target_infinite_warning), 
        logical(1L))), any_core_warning = any(vapply(rows, function(z) isTRUE(z$core_warning), logical(1L))), any_unmapped_warning = any(vapply(rows, 
        function(z) !isTRUE(z$completely_mapped), logical(1L))))
}

rw_information <- function(fit, core) {
    detail <- survival::coxph.detail(fit, riskmat = FALSE)
    p <- length(stats::coef(fit))
    information <- if (p == 1L) 
        matrix(sum(detail$imat), 1L)
    else apply(detail$imat, c(1L, 2L), sum)
    core_information <- information[core, core, drop = FALSE]
    V <- stats::vcov(fit)[core, core, drop = FALSE]
    eig_I <- eigen((core_information + t(core_information))/2, symmetric = TRUE, only.values = TRUE)$values
    eig_V <- eigen((V + t(V))/2, symmetric = TRUE, only.values = TRUE)$values
    profile_information <- solve(V)
    eig_profile <- eigen((profile_information + t(profile_information))/2, symmetric = TRUE, only.values = TRUE)$values
    list(core_information_min_eigenvalue = min(eig_I), core_covariance_min_eigenvalue = min(eig_V), inverse_core_covariance_min_eigenvalue = min(eig_profile), 
        positive = all(is.finite(c(eig_I, eig_V, eig_profile))) && min(eig_I) > 0 && min(eig_V) > 0 && min(eig_profile) > 
            0, finite_joint_MLE_certified = FALSE, interpretation = "CORE_BLOCK_AND_PROFILE_CURVATURE_NOT_FINITE_JOINT_MLE_CERTIFICATE")
}

rw_support <- function(dd, endpoint, w, indices) {
    lapply(indices, function(k) {
        x <- dd$X[, k]
        nz <- x != 0
        event <- endpoint$status
        list(term = dd$reference$columns[k], minimum = min(x), maximum = max(x), N_nonzero = sum(nz), N_zero = sum(!nz), 
            disease_events_nonzero = sum(nz & event == 1L), disease_events_zero = sum(!nz & event == 1L), death_events_nonzero = sum(nz & 
                event == 2L), death_events_zero = sum(!nz & event == 2L), weight_sum_nonzero = sum(w[nz]), disease_event_weight_nonzero = sum(w[nz & 
                event == 1L]), death_event_weight_nonzero = sum(w[nz & event == 2L]))
    })
}

rw_replace <- function(tree, from, to, expected = 1L) {
    n <- 0L
    visit <- function(z) {
        if (identical(z, from)) {
            n <<- n + 1L
            return(to)
        }
        if (is.call(z)) 
            for (i in seq_along(z)) z[[i]] <- visit(z[[i]])
        z
    }
    answer <- visit(tree)
    rw_assert(n == expected, paste("RECOVERY_NATIVE_AST_PATCH_COUNT_CHANGED", n, expected))
    answer
}

rw_clone_fit <- function(state, column, private_dir) {
    native <- ff_ce_fit
    env <- new.env(parent = environment(native))
    env$rw_original_eligible <- function(warnings, finite, reference) {
        map <- rw_warning_map(warnings, reference, column)
        rw_assert(!map$any_core_warning, "RECOVERY_CORE_INFINITE_WARNING_UNRESOLVED")
        finite && map$all_non_target_infinite
    }
    env$rw_retry_gate <- function(warnings, rec, reference) {
        map <- rw_warning_map(warnings, reference, column)
        rw_assert(map$all_non_target_infinite && !map$any_core_warning, "RECOVERY_SINGLE_WARM_NATIVE_WARNING_UNRESOLVED")
        state$warning_gates[[length(state$warning_gates) + 1L]] <- list(cause = rec$cause, native_gate = "CE_UNRESOLVED_NATIVE_WARNING_AFTER_SINGLE_RETRY", 
            native_gate_passed = FALSE, warning_mapping = map, disposition = "NON_TARGET_WARNING_DIAGNOSTIC_CONTINUATION_NOT_VALIDITY")
        invisible(TRUE)
    }
    env$rw_native <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL) {
        family <- if (is.null(init)) 
            "original"
        else "warm"
        cause_name <- c("disease", "death")[cause]
        key <- paste(family, cause_name, sep = "_")
        rw_assert(is.null(state$models[[key]]), "RECOVERY_MORE_THAN_ONE_ORIGINAL_OR_WARM_FIT_FORBIDDEN")
        state$fits_attempted <- state$fits_attempted + 1L
        result <- tryCatch(ff_ce_native(work, safe_names, cause, stratified, cluster_bound, control, init), error = function(e) {
            state$native_errors[[key]] <- list(message = conditionMessage(e), condition_class = class(e), native_singular = grepl("singular", 
                conditionMessage(e), ignore.case = TRUE))
            stop(e)
        })
        state$models[[key]] <- result$value
        state$warnings[[key]] <- result$warnings
        state$controls[[key]] <- control
        path <- file.path(private_dir, paste0(key, ".private.rds"))
        rw_assert(!file.exists(path), "RECOVERY_MODEL_ALREADY_EXISTS")
        saveRDS(result$value, path, compress = FALSE, version = 2)
        state$model_files[[key]] <- list(path = path, sha256 = rw_sha(path))
        result
    }
    tree <- body(native)
    tree <- rw_replace(tree, quote(if (ff_ce_bad_warning(got$warnings)) {
        eligible <- length(got$warnings) > 0L && finite && all(grepl("Loglik converged before variable.*coefficient may be infinite", 
            got$warnings))
        if (!eligible) ff_ce_fail("CE_UNRESOLVED_ORIGINAL_NATIVE_WARNING_NO_ELIGIBLE_RETRY", rec)
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
        if (ff_ce_bad_warning(got$warnings)) ff_ce_fail("CE_UNRESOLVED_NATIVE_WARNING_AFTER_SINGLE_RETRY", rec)
    }), quote(if (TRUE) {
        eligible <- rw_original_eligible(got$warnings, finite, dd$reference)
        if (!eligible) ff_ce_fail("CE_UNRESOLVED_ORIGINAL_NATIVE_WARNING_NO_ELIGIBLE_RETRY", rec)
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
        if (ff_ce_bad_warning(got$warnings)) rw_retry_gate(got$warnings, rec, dd$reference)
    }))
    rename <- function(z) {
        if (is.call(z)) {
            if (identical(z[[1L]], as.name("ff_ce_native"))) 
                z[[1L]] <- as.name("rw_native")
            for (i in seq_along(z)) if (i != 1L) 
                z[[i]] <- rename(z[[i]])
        }
        z
    }
    tree <- rename(tree)
    body(native) <- tree
    environment(native) <- env
    state$clone_body_sha256 <- digest::digest(paste(deparse(tree, width.cutoff = 500L), collapse = "\n"), algo = "sha256", 
        serialize = FALSE)
    native
}

rw_rebuild_family <- function(template, models, X, risk, w, centres, tau) {
    out <- template
    out$fit <- models
    endpoint <- ff_y_risk(risk, tau)
    original <- ff_y_matrix(X, template$matrix_reference)
    changed <- ff_y_matrix(ff_ce_columns(original$X, template$center_reference$columns))
    designs <- list(disease = if (template$stratified[1L]) changed else original, death = changed)
    out$engine$baseline_public_validation <- vector("list", 2L)
    names(out$engine$baseline_public_validation) <- names(models)
    lp_max <- numeric(2L)
    probe <- unique(as.integer(round(seq(1, endpoint$n, length.out = min(16L, endpoint$n)))))
    for (j in 1:2) {
        m <- models[[j]]
        dd <- designs[[j]]
        mf <- stats::model.frame(m)
        ff_y_assert(inherits(m, "coxph") && !inherits(m, "coxphms") && identical(m$method, "efron") && all(is.finite(stats::coef(m))) && 
            all(is.finite(stats::vcov(m))), "RECOVERY_NATIVE_NONFINITE_OR_WRONG_ESTIMATOR")
        ff_y_assert(m$n == endpoint$n && identical(as.numeric(m$y[, "time"]), as.numeric(endpoint$time)) && identical(as.numeric(m$y[, 
            "status"]), as.numeric(endpoint$status == j)) && identical(as.numeric(stats::model.weights(mf)), as.numeric(w)) && 
            isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), "RECOVERY_ORIGINAL_FROZEN_TIME_EVENT_WEIGHT_OR_DESIGN_CHANGED")
        if (template$stratified[j]) 
            ff_y_assert(identical(as.character(mf[["strata(.ce_center)"]]), centres$group), "RECOVERY_ORIGINAL_NATIVE_STRATUM_BINDING_CHANGED")
        out$engine$baseline_public_validation[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, template$stratified[j], 
            template$center_reference$observed_levels, template$center_support, j)
        nd <- dd$work[probe, , drop = FALSE]
        nd$.ce_center <- factor(centres$group[probe], levels = centres$reference$observed_levels)
        lp <- ff_y_capture(as.numeric(predict(m, newdata = nd, type = "lp", reference = "zero")))
        ff_y_assert(!ff_ce_bad_warning(lp$warnings), "RECOVERY_ORIGINAL_NATIVE_LP_WARNING")
        lp_max[j] <- max(abs(lp$value - as.numeric(dd$X[probe, , drop = FALSE] %*% stats::coef(m))))
        ff_y_assert(lp_max[j] <= 1e-10, "RECOVERY_ORIGINAL_NATIVE_LP_MATRIX_BINDING_FAILED")
    }
    baselines <- out$engine$baseline_public_validation
    out$engine$by_center <- stats::setNames(lapply(centres$reference$observed_levels, function(g) ff_ce_merge_baselines(list(if (out$stratified[1L]) baselines$disease[[g]] else baselines$disease$GLOBAL, 
        baselines$death[[g]]))), centres$reference$observed_levels)
    all_baselines <- unlist(baselines, recursive = FALSE)
    out$native_validation <- list(passed = FALSE, lp_max_abs = max(lp_max), basehaz_survfit_max_abs = max(vapply(all_baselines, 
        `[[`, numeric(1L), "max_abs_basehaz_vs_survfit")))
    curve <- lapply(c(0, 0.25, 0.5, 0.75, 1) * tau, function(t) ff_ce_predict_at(out, original$X[probe, , drop = FALSE], 
        t))
    violation <- 0
    for (k in 2:length(curve)) {
        delta <- curve[[k]] - curve[[k - 1L]]
        violation <- max(violation, -delta[, "disease"], -delta[, "death"], delta[, "survival"])
    }
    ff_y_assert(violation <= 1e-12, "RECOVERY_ORIGINAL_CIF_MONOTONICITY_FAILED")
    out$native_validation$monotonicity_violation <- violation
    out$native_validation$three_state_closure <- max(vapply(curve, function(z) max(abs(rowSums(z) - 1)), numeric(1L)))
    out$native_validation$passed <- TRUE
    out
}

rw_probability_pair <- function(original, warm, X, w, tau) {
    p <- ff_ce_predict_y_matrix(original, X, tau)
    q <- ff_ce_predict_y_matrix(warm, X, tau)
    summarize <- function(z) list(rows = nrow(z), state_min = as.list(apply(z, 2L, min)), state_max = as.list(apply(z, 
        2L, max)), max_abs_closure = max(abs(rowSums(z) - 1)), nonfinite = sum(!is.finite(z)), outside_unit_interval = sum(z < 
        0 | z > 1), closure_legality_tolerance = 1e-08, legal = all(is.finite(z)) && all(z >= 0 & z <= 1) && max(abs(rowSums(z) - 
        1)) < 1e-08)
    maxdiff <- max(abs(p - q))
    list(original = summarize(p), warm = summarize(q), max_abs_individual_state_movement = maxdiff, state_max_abs_individual_movement = as.list(apply(abs(p - 
        q), 2L, max)), cif_tolerance_probability_units = 1e-09, stability_rule_is_descriptive_not_original_V4_acceptance_policy = TRUE, 
        individual_probability_stability_passed = is.finite(maxdiff) && maxdiff <= 1e-09, original_weighted_disease_risk = weighted.mean(p[, 
            "disease"], w), warm_weighted_disease_risk = weighted.mean(q[, "disease"], w))
}

rw_effects <- function(risk) {
    TE <- unname(risk[["R11"]] - risk[["R00"]])
    DE <- unname(risk[["R10"]] - risk[["R00"]])
    IE <- unname(risk[["R11"]] - risk[["R10"]])
    list(R00 = unname(risk[["R00"]]), R10 = unname(risk[["R10"]]), R11 = unname(risk[["R11"]]), TE = TE, DE = DE, 
        IE = IE, PM = if (TE != 0) IE/TE else NA_real_, closure = TE - DE - IE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        rw_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        rw_json <- function(value, path) {
            rw_assert(!file.exists(path), paste("PRESERVE_EXISTING_RECOVERY_ARTIFACT", path))
            jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
        }
        rw_sha <- function(path) {
            rw_assert(file.exists(path) && !dir.exists(path), paste("RECOVERY_FILE_MISSING", path))
            digest::digest(file = path, algo = "sha256")
        }
        rw_check <- function(path, sha) {
            rw_assert(is.character(sha) && length(sha) == 1L && grepl("^[a-f0-9]{64}$", sha), "RECOVERY_EXACT_SHA256_REQUIRED")
            rw_assert(identical(rw_sha(path), sha), paste("RECOVERY_SOURCE_SHA256_CHANGED", path))
        }
        rw_inside <- function(path, parent) {
            path <- normalizePath(path, mustWork = TRUE)
            parent <- normalizePath(parent, mustWork = TRUE)
            startsWith(paste0(path, "/"), paste0(parent, "/"))
        }
        rw_newdir <- function(path) {
            rw_assert(!file.exists(path) && dir.create(path, recursive = TRUE, showWarnings = FALSE), paste("RECOVERY_NEW_ISOLATED_DIRECTORY_REQUIRED", 
                path))
            normalizePath(path, mustWork = TRUE)
        }
        rw_core <- function(reference, column) {
            ix <- which(reference$columns == "B" | grepl(column, reference$columns, fixed = TRUE))
            rw_assert(sum(reference$columns[ix] == "B") == 1L && length(ix) >= 2L, "RECOVERY_B_AND_PROTEIN_SPLINE_CORE_MAPPING_REQUIRED")
            ix
        }
        rw_warning_map <- function(warnings, reference, column) {
            core <- rw_core(reference, column)
            rows <- lapply(warnings, function(w) {
                pattern <- "(?i)variable[s]?[[:space:]]+([0-9]+(?:[[:space:],]+[0-9]+)*)"
                hit <- regmatches(w, regexec(pattern, w, perl = TRUE))[[1L]]
                indices <- if (length(hit) >= 2L) as.integer(regmatches(hit[2L], gregexpr("[0-9]+", hit[2L], perl = TRUE))[[1L]]) else integer()
                eligible_text <- grepl("^Loglik converged before variable.*coefficient may be infinite[.]?[[:space:]]*$", 
                  trimws(w), ignore.case = TRUE)
                mapped <- length(indices) > 0L && !anyNA(indices) && all(indices %in% seq_along(reference$columns))
                list(message = w, parsed_indices = as.list(indices), completely_mapped = mapped, actual_terms = if (mapped) as.list(reference$columns[indices]) else list(), 
                  internal_terms = if (mapped) as.list(reference$safe_names[indices]) else list(), core_warning = if (mapped) any(indices %in% 
                    core) else NA, eligible_non_target_infinite_warning = eligible_text && mapped && !any(indices %in% 
                    core))
            })
            list(warnings = rows, all_non_target_infinite = !length(rows) || all(vapply(rows, function(z) isTRUE(z$eligible_non_target_infinite_warning), 
                logical(1L))), any_core_warning = any(vapply(rows, function(z) isTRUE(z$core_warning), logical(1L))), 
                any_unmapped_warning = any(vapply(rows, function(z) !isTRUE(z$completely_mapped), logical(1L))))
        }
        rw_information <- function(fit, core) {
            detail <- survival::coxph.detail(fit, riskmat = FALSE)
            p <- length(stats::coef(fit))
            information <- if (p == 1L) matrix(sum(detail$imat), 1L) else apply(detail$imat, c(1L, 2L), sum)
            core_information <- information[core, core, drop = FALSE]
            V <- stats::vcov(fit)[core, core, drop = FALSE]
            eig_I <- eigen((core_information + t(core_information))/2, symmetric = TRUE, only.values = TRUE)$values
            eig_V <- eigen((V + t(V))/2, symmetric = TRUE, only.values = TRUE)$values
            profile_information <- solve(V)
            eig_profile <- eigen((profile_information + t(profile_information))/2, symmetric = TRUE, only.values = TRUE)$values
            list(core_information_min_eigenvalue = min(eig_I), core_covariance_min_eigenvalue = min(eig_V), inverse_core_covariance_min_eigenvalue = min(eig_profile), 
                positive = all(is.finite(c(eig_I, eig_V, eig_profile))) && min(eig_I) > 0 && min(eig_V) > 0 && 
                  min(eig_profile) > 0, finite_joint_MLE_certified = FALSE, interpretation = "CORE_BLOCK_AND_PROFILE_CURVATURE_NOT_FINITE_JOINT_MLE_CERTIFICATE")
        }
        rw_support <- function(dd, endpoint, w, indices) {
            lapply(indices, function(k) {
                x <- dd$X[, k]
                nz <- x != 0
                event <- endpoint$status
                list(term = dd$reference$columns[k], minimum = min(x), maximum = max(x), N_nonzero = sum(nz), N_zero = sum(!nz), 
                  disease_events_nonzero = sum(nz & event == 1L), disease_events_zero = sum(!nz & event == 1L), 
                  death_events_nonzero = sum(nz & event == 2L), death_events_zero = sum(!nz & event == 2L), weight_sum_nonzero = sum(w[nz]), 
                  disease_event_weight_nonzero = sum(w[nz & event == 1L]), death_event_weight_nonzero = sum(w[nz & 
                    event == 2L]))
            })
        }
        rw_replace <- function(tree, from, to, expected = 1L) {
            n <- 0L
            visit <- function(z) {
                if (identical(z, from)) {
                  n <<- n + 1L
                  return(to)
                }
                if (is.call(z)) for (i in seq_along(z)) z[[i]] <- visit(z[[i]])
                z
            }
            answer <- visit(tree)
            rw_assert(n == expected, paste("RECOVERY_NATIVE_AST_PATCH_COUNT_CHANGED", n, expected))
            answer
        }
        rw_clone_fit <- function(state, column, private_dir) {
            native <- ff_ce_fit
            env <- new.env(parent = environment(native))
            env$rw_original_eligible <- function(warnings, finite, reference) {
                map <- rw_warning_map(warnings, reference, column)
                rw_assert(!map$any_core_warning, "RECOVERY_CORE_INFINITE_WARNING_UNRESOLVED")
                finite && map$all_non_target_infinite
            }
            env$rw_retry_gate <- function(warnings, rec, reference) {
                map <- rw_warning_map(warnings, reference, column)
                rw_assert(map$all_non_target_infinite && !map$any_core_warning, "RECOVERY_SINGLE_WARM_NATIVE_WARNING_UNRESOLVED")
                state$warning_gates[[length(state$warning_gates) + 1L]] <- list(cause = rec$cause, native_gate = "CE_UNRESOLVED_NATIVE_WARNING_AFTER_SINGLE_RETRY", 
                  native_gate_passed = FALSE, warning_mapping = map, disposition = "NON_TARGET_WARNING_DIAGNOSTIC_CONTINUATION_NOT_VALIDITY")
                invisible(TRUE)
            }
            env$rw_native <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL) {
                family <- if (is.null(init)) "original" else "warm"
                cause_name <- c("disease", "death")[cause]
                key <- paste(family, cause_name, sep = "_")
                rw_assert(is.null(state$models[[key]]), "RECOVERY_MORE_THAN_ONE_ORIGINAL_OR_WARM_FIT_FORBIDDEN")
                state$fits_attempted <- state$fits_attempted + 1L
                result <- tryCatch(ff_ce_native(work, safe_names, cause, stratified, cluster_bound, control, init), 
                  error = function(e) {
                    state$native_errors[[key]] <- list(message = conditionMessage(e), condition_class = class(e), 
                      native_singular = grepl("singular", conditionMessage(e), ignore.case = TRUE))
                    stop(e)
                  })
                state$models[[key]] <- result$value
                state$warnings[[key]] <- result$warnings
                state$controls[[key]] <- control
                path <- file.path(private_dir, paste0(key, ".private.rds"))
                rw_assert(!file.exists(path), "RECOVERY_MODEL_ALREADY_EXISTS")
                saveRDS(result$value, path, compress = FALSE, version = 2)
                state$model_files[[key]] <- list(path = path, sha256 = rw_sha(path))
                result
            }
            tree <- body(native)
            tree <- rw_replace(tree, quote(if (ff_ce_bad_warning(got$warnings)) {
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
            }), quote(if (TRUE) {
                eligible <- rw_original_eligible(got$warnings, finite, dd$reference)
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
                if (ff_ce_bad_warning(got$warnings)) rw_retry_gate(got$warnings, rec, dd$reference)
            }))
            rename <- function(z) {
                if (is.call(z)) {
                  if (identical(z[[1L]], as.name("ff_ce_native"))) z[[1L]] <- as.name("rw_native")
                  for (i in seq_along(z)) if (i != 1L) z[[i]] <- rename(z[[i]])
                }
                z
            }
            tree <- rename(tree)
            body(native) <- tree
            environment(native) <- env
            state$clone_body_sha256 <- digest::digest(paste(deparse(tree, width.cutoff = 500L), collapse = "\n"), 
                algo = "sha256", serialize = FALSE)
            native
        }
        rw_rebuild_family <- function(template, models, X, risk, w, centres, tau) {
            out <- template
            out$fit <- models
            endpoint <- ff_y_risk(risk, tau)
            original <- ff_y_matrix(X, template$matrix_reference)
            changed <- ff_y_matrix(ff_ce_columns(original$X, template$center_reference$columns))
            designs <- list(disease = if (template$stratified[1L]) changed else original, death = changed)
            out$engine$baseline_public_validation <- vector("list", 2L)
            names(out$engine$baseline_public_validation) <- names(models)
            lp_max <- numeric(2L)
            probe <- unique(as.integer(round(seq(1, endpoint$n, length.out = min(16L, endpoint$n)))))
            for (j in 1:2) {
                m <- models[[j]]
                dd <- designs[[j]]
                mf <- stats::model.frame(m)
                ff_y_assert(inherits(m, "coxph") && !inherits(m, "coxphms") && identical(m$method, "efron") && 
                  all(is.finite(stats::coef(m))) && all(is.finite(stats::vcov(m))), "RECOVERY_NATIVE_NONFINITE_OR_WRONG_ESTIMATOR")
                ff_y_assert(m$n == endpoint$n && identical(as.numeric(m$y[, "time"]), as.numeric(endpoint$time)) && 
                  identical(as.numeric(m$y[, "status"]), as.numeric(endpoint$status == j)) && identical(as.numeric(stats::model.weights(mf)), 
                  as.numeric(w)) && isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), 
                  "RECOVERY_ORIGINAL_FROZEN_TIME_EVENT_WEIGHT_OR_DESIGN_CHANGED")
                if (template$stratified[j]) ff_y_assert(identical(as.character(mf[["strata(.ce_center)"]]), centres$group), 
                  "RECOVERY_ORIGINAL_NATIVE_STRATUM_BINDING_CHANGED")
                out$engine$baseline_public_validation[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, template$stratified[j], 
                  template$center_reference$observed_levels, template$center_support, j)
                nd <- dd$work[probe, , drop = FALSE]
                nd$.ce_center <- factor(centres$group[probe], levels = centres$reference$observed_levels)
                lp <- ff_y_capture(as.numeric(predict(m, newdata = nd, type = "lp", reference = "zero")))
                ff_y_assert(!ff_ce_bad_warning(lp$warnings), "RECOVERY_ORIGINAL_NATIVE_LP_WARNING")
                lp_max[j] <- max(abs(lp$value - as.numeric(dd$X[probe, , drop = FALSE] %*% stats::coef(m))))
                ff_y_assert(lp_max[j] <= 1e-10, "RECOVERY_ORIGINAL_NATIVE_LP_MATRIX_BINDING_FAILED")
            }
            baselines <- out$engine$baseline_public_validation
            out$engine$by_center <- stats::setNames(lapply(centres$reference$observed_levels, function(g) ff_ce_merge_baselines(list(if (out$stratified[1L]) baselines$disease[[g]] else baselines$disease$GLOBAL, 
                baselines$death[[g]]))), centres$reference$observed_levels)
            all_baselines <- unlist(baselines, recursive = FALSE)
            out$native_validation <- list(passed = FALSE, lp_max_abs = max(lp_max), basehaz_survfit_max_abs = max(vapply(all_baselines, 
                `[[`, numeric(1L), "max_abs_basehaz_vs_survfit")))
            curve <- lapply(c(0, 0.25, 0.5, 0.75, 1) * tau, function(t) ff_ce_predict_at(out, original$X[probe, 
                , drop = FALSE], t))
            violation <- 0
            for (k in 2:length(curve)) {
                delta <- curve[[k]] - curve[[k - 1L]]
                violation <- max(violation, -delta[, "disease"], -delta[, "death"], delta[, "survival"])
            }
            ff_y_assert(violation <= 1e-12, "RECOVERY_ORIGINAL_CIF_MONOTONICITY_FAILED")
            out$native_validation$monotonicity_violation <- violation
            out$native_validation$three_state_closure <- max(vapply(curve, function(z) max(abs(rowSums(z) - 1)), 
                numeric(1L)))
            out$native_validation$passed <- TRUE
            out
        }
        rw_probability_pair <- function(original, warm, X, w, tau) {
            p <- ff_ce_predict_y_matrix(original, X, tau)
            q <- ff_ce_predict_y_matrix(warm, X, tau)
            summarize <- function(z) list(rows = nrow(z), state_min = as.list(apply(z, 2L, min)), state_max = as.list(apply(z, 
                2L, max)), max_abs_closure = max(abs(rowSums(z) - 1)), nonfinite = sum(!is.finite(z)), outside_unit_interval = sum(z < 
                0 | z > 1), closure_legality_tolerance = 1e-08, legal = all(is.finite(z)) && all(z >= 0 & z <= 
                1) && max(abs(rowSums(z) - 1)) < 1e-08)
            maxdiff <- max(abs(p - q))
            list(original = summarize(p), warm = summarize(q), max_abs_individual_state_movement = maxdiff, state_max_abs_individual_movement = as.list(apply(abs(p - 
                q), 2L, max)), cif_tolerance_probability_units = 1e-09, stability_rule_is_descriptive_not_original_V4_acceptance_policy = TRUE, 
                individual_probability_stability_passed = is.finite(maxdiff) && maxdiff <= 1e-09, original_weighted_disease_risk = weighted.mean(p[, 
                  "disease"], w), warm_weighted_disease_risk = weighted.mean(q[, "disease"], w))
        }
        rw_effects <- function(risk) {
            TE <- unname(risk[["R11"]] - risk[["R00"]])
            DE <- unname(risk[["R10"]] - risk[["R00"]])
            IE <- unname(risk[["R11"]] - risk[["R10"]])
            list(R00 = unname(risk[["R00"]]), R10 = unname(risk[["R10"]]), R11 = unname(risk[["R11"]]), TE = TE, 
                DE = DE, IE = IE, PM = if (TE != 0) IE/TE else NA_real_, closure = TE - DE - IE)
        }
    }), envir = context)
}

