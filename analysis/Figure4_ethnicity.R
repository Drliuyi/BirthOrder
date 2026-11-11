# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

c39_uy_assert <- function(ok, msg) if (!isTRUE(ok)) stop(paste0("CORE39_UY_", msg), call. = FALSE)

c39_uy_pool_columns <- function(X) {
    a <- "ethnicityAsian or Asian British"
    c <- "ethnicityChinese"
    c39_uy_assert(is.matrix(X) && all(c(a, c) %in% colnames(X)), "ETHNICITY_COLUMNS_REQUIRED")
    c39_uy_assert(all(X[, a] %in% 0:1) && all(X[, c] %in% 0:1) && all(X[, a] + X[, c] <= 1), "ETHNICITY_DISJOINT_REQUIRED")
    z <- X[, colnames(X) != c, drop = FALSE]
    z[, a] <- X[, a] + X[, c]
    colnames(z)[colnames(z) == a] <- "ethnicityAsian_or_Chinese_pooled"
    z
}

c39_uy_replace <- function(fn, old, new, expected = 1L) {
    count <- 0L
    walk <- function(x) {
        if (identical(x, old)) {
            count <<- count + 1L
            return(new)
        }
        if (is.call(x)) 
            return(as.call(lapply(as.list(x), walk)))
        if (is.pairlist(x)) 
            return(as.pairlist(lapply(x, walk)))
        x
    }
    body(fn) <- walk(body(fn))
    c39_uy_assert(count == expected, "PINNED_FUNCTION_EXPRESSION_DRIFT")
    fn
}

c39_uy_death_cache <- function(old_y, provenance) {
    c39_uy_assert(inherits(old_y, "ff_ce_matrix_fit") && isTRUE(old_y$native_validation$passed) && identical(as.logical(old_y$stratified), 
        c(FALSE, TRUE)) && old_y$endpoint_id %in% c("250.2", "411"), "ACCEPTED_OLD_Y_REQUIRED")
    c39_uy_assert(is.list(provenance) && isTRUE(provenance$sha256_verified) && is.character(provenance$sha256) && 
        grepl("^[a-f0-9]{64}$", provenance$sha256), "CALLER_VERIFIED_CACHE_PROVENANCE_REQUIRED")
    audit <- old_y$numerical_audit$death
    final_warnings <- if (isTRUE(audit$retry_attempted)) 
        audit$retry_warnings
    else audit$original_warnings
    c39_uy_assert(is.list(audit) && !any(grepl("infinite|converg|singular|iteration", unlist(final_warnings), ignore.case = TRUE)), 
        "UNRESOLVED_DEATH_CACHE_WARNING")
    list(fit = old_y$fit$death, reference = old_y$cause_references$death, center_reference = old_y$center_reference, 
        endpoint_id = old_y$endpoint_id, baseline = old_y$engine$baseline_public_validation$death, source_numerical_audit = audit, 
        provenance = provenance, source_method_hash = old_y$method_spec_hash)
}

c39_uy_native_or_reuse <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL, cache = NULL) {
    if (is.null(cache)) 
        return(ff_ce_native(work, safe_names, cause, stratified, cluster_bound, control, init))
    c39_uy_assert(identical(cause, 2L) && isTRUE(stratified) && isTRUE(cluster_bound) && is.null(init), "DEATH_ONLY_CACHE_REUSE")
    m <- cache$fit
    mf <- stats::model.frame(m)
    c39_uy_assert(inherits(m, "coxph") && !inherits(m, "coxph.penal") && identical(m$method, "efron") && all(is.finite(stats::coef(m))) && 
        all(is.finite(stats::vcov(m))) && !is.null(m$naive.var), "NATIVE_CACHE_ESTIMATOR_INVALID")
    c39_uy_assert(identical(names(stats::coef(m)), safe_names) && identical(cache$reference$safe_names, safe_names) && 
        identical(dim(m$x), dim(as.matrix(work[, safe_names, drop = FALSE]))) && identical(colnames(m$x), safe_names) && 
        identical(as.numeric(m$x), as.numeric(as.matrix(work[, safe_names, drop = FALSE]))) && identical(as.numeric(m$y[, 
        "time"]), as.numeric(work$.ff_time)) && identical(as.numeric(m$y[, "status"]), as.numeric(work$.ff_death)) && 
        identical(as.numeric(stats::model.weights(mf)), as.numeric(work$.ff_weight)) && identical(as.character(stats::model.extract(mf, 
        "cluster")), as.character(work$.ff_cluster)) && identical(as.character(mf[["strata(.ce_center)"]]), as.character(work$.ce_center)) && 
        identical(levels(work$.ce_center), cache$center_reference$observed_levels), "DEATH_CACHE_EXACT_BINDING_MISMATCH")
    c39_uy_assert(identical(attr(m$terms, "term.labels"), c(safe_names, "strata(.ce_center)")) && !length(m$offset) && 
        !length(m$na.action) && is.null(m$frail), "DEATH_CACHE_TERM_DRIFT")
    list(value = m, warnings = character(), cache_reuse = list(reused = TRUE, optimizer_called = FALSE, exact_X_time_status_weight_cluster_center_binding = TRUE, 
        provenance = cache$provenance, source_method_hash = cache$source_method_hash, source_numerical_audit = cache$source_numerical_audit))
}

c39_uy_verify_death_baseline <- function(current, cache) {
    if (is.null(cache)) 
        return(invisible(TRUE))
    c39_uy_assert(identical(names(current), names(cache$baseline)) && all(vapply(names(current), function(g) identical(current[[g]]$time, 
        cache$baseline[[g]]$time) && identical(current[[g]]$cumulative, cache$baseline[[g]]$cumulative), FALSE)), 
        "DEATH_CACHE_BASELINE_CHANGED")
    invisible(TRUE)
}

c39_uy_install <- function(engine) {
    c39_uy_assert(is.environment(engine) && is.function(engine$ff_ce_fit) && is.function(engine$fd_fp_prepare), 
        "FROZEN_ENGINE_REQUIRED")
    old_spec <- engine$ff_ce_method_spec()
    c39_uy_assert(identical(old_spec$id, "CENTER_ETHNICITY_V4_20261028"), "INSTALL_ONCE_BEFORE_OTHER_METHOD_OVERLAYS")
    helpers <- c("c39_uy_assert", "c39_uy_pool_columns", "c39_uy_native_or_reuse", "c39_uy_verify_death_baseline")
    for (n in helpers) {
        f <- get(n, environment())
        environment(f) <- engine
        engine[[n]] <- f
    }
    spec_text <- paste(c("CORE39_Y_ASIAN_CHINESE_POOL_V1_20261101", "scope=250.2,411;all_B0_to_B100;no_failure_conditional_pooling", 
        "disease=Asian_Chinese_pool+unchanged_center_fixed_covariates", "death=unchanged_Asian_Chinese_pool+unchanged_center_strata", 
        "A+input_ethnicity+sample+draw+selection_weights+Efron+robust_cluster+tau11=UNCHANGED", "death_cache_reuse=accepted_native_only_exact_X_time_status_weights_cluster_center_baseline", 
        "disease_cache_reuse=FORBIDDEN;old_whole_Y_cache=INCOMPATIBLE", paste0("parent_method_sha256=", old_spec$sha256), 
        "strict_warning_gate_and_one_native_numerical_retry=UNCHANGED"), collapse = "\n")
    method <- list(id = "CORE39_Y_ASIAN_CHINESE_POOL_V1_20261101", canonical = spec_text, sha256 = digest::digest(spec_text, 
        algo = "sha256", serialize = FALSE))
    engine$ff_ce_method_spec <- local({
        value <- method
        function() value
    })
    f <- engine$ff_ce_fit
    f <- c39_uy_replace(f, quote(c("250", "250.2", "411", "411.3", "411.4")), quote(c("250.2", "411")))
    f <- c39_uy_replace(f, quote(if (stratified[1L]) changed else original), quote(ff_y_matrix(c39_uy_pool_columns(original$X))))
    f <- c39_uy_replace(f, quote(got <- ff_ce_native(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), 
        ctrl0)), quote(got <- c39_uy_native_or_reuse(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), 
        ctrl0, cache = if (j == 2L) death_cache else NULL)))
    f <- c39_uy_replace(f, quote(emit(c(rec, list(stage = "ORIGINAL_ATTEMPT_PRESERVED")))), quote({
        if (!is.null(got$cache_reuse)) {
            rec$cache_reuse <- got$cache_reuse
            rec$original_control <- NULL
        }
        emit(c(rec, list(stage = if (is.null(got$cache_reuse)) "ORIGINAL_ATTEMPT_PRESERVED" else "REUSED_ACCEPTED_DEATH_MODEL_PRESERVED")))
    }))
    f <- c39_uy_replace(f, quote(baselines[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, stratified[j], centers$reference$observed_levels, 
        support, j)), quote({
        baselines[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, stratified[j], centers$reference$observed_levels, 
            support, j)
        if (j == 2L) c39_uy_verify_death_baseline(baselines[[j]], death_cache)
    }))
    formals(f) <- c(formals(f), alist(death_cache = NULL))
    body(f) <- substitute({
        if (!is.null(death_cache)) c39_uy_assert(identical(death_cache$endpoint_id, endpoint_id), "DEATH_CACHE_ENDPOINT_MISMATCH")
        BODY
    }, list(BODY = body(f)))
    environment(f) <- engine
    engine$ff_ce_fit <- f
    p <- engine$ff_ce_predict_at
    p <- c39_uy_replace(p, quote(if (object$stratified[1L]) changed else original$X), quote(c39_uy_pool_columns(original$X)))
    body(p) <- substitute({
        c39_uy_assert(identical(object$method_spec_hash, ff_ce_method_spec()$sha256) && identical(as.logical(object$stratified), 
            c(FALSE, TRUE)), "INCOMPATIBLE_Y_CACHE")
        BODY
    }, list(BODY = body(p)))
    environment(p) <- engine
    engine$ff_ce_predict_at <- p
    p <- engine$ff_ce_predict_y_matrix
    environment(p) <- engine
    engine$ff_ce_predict_y_matrix <- p
    engine$fd_fp_prepare <- eval(quote(function(X, compiled, v4) {
        o <- compiled$object
        c39_uy_assert(identical(o$method_spec_hash, v4$ff_ce_method_spec()$sha256) && identical(as.logical(o$stratified), 
            c(FALSE, TRUE)), "INCOMPATIBLE_FAST_Y_CACHE")
        X <- fd_fp_matrix(X, o$matrix_reference)
        centers <- v4$ff_ce_centers(X, o$center_reference)
        changed <- v4$ff_ce_columns(X, o$center_reference$columns)
        list(designs = list(fd_fp_matrix(c39_uy_pool_columns(X), o$cause_references$disease), fd_fp_matrix(changed, 
            o$cause_references$death)), group = centers$group)
    }), envir = engine)
    invisible(method)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        c39_uy_assert <- function(ok, msg) if (!isTRUE(ok)) stop(paste0("CORE39_UY_", msg), call. = FALSE)
        c39_uy_pool_columns <- function(X) {
            a <- "ethnicityAsian or Asian British"
            c <- "ethnicityChinese"
            c39_uy_assert(is.matrix(X) && all(c(a, c) %in% colnames(X)), "ETHNICITY_COLUMNS_REQUIRED")
            c39_uy_assert(all(X[, a] %in% 0:1) && all(X[, c] %in% 0:1) && all(X[, a] + X[, c] <= 1), "ETHNICITY_DISJOINT_REQUIRED")
            z <- X[, colnames(X) != c, drop = FALSE]
            z[, a] <- X[, a] + X[, c]
            colnames(z)[colnames(z) == a] <- "ethnicityAsian_or_Chinese_pooled"
            z
        }
        c39_uy_replace <- function(fn, old, new, expected = 1L) {
            count <- 0L
            walk <- function(x) {
                if (identical(x, old)) {
                  count <<- count + 1L
                  return(new)
                }
                if (is.call(x)) return(as.call(lapply(as.list(x), walk)))
                if (is.pairlist(x)) return(as.pairlist(lapply(x, walk)))
                x
            }
            body(fn) <- walk(body(fn))
            c39_uy_assert(count == expected, "PINNED_FUNCTION_EXPRESSION_DRIFT")
            fn
        }
        c39_uy_death_cache <- function(old_y, provenance) {
            c39_uy_assert(inherits(old_y, "ff_ce_matrix_fit") && isTRUE(old_y$native_validation$passed) && identical(as.logical(old_y$stratified), 
                c(FALSE, TRUE)) && old_y$endpoint_id %in% c("250.2", "411"), "ACCEPTED_OLD_Y_REQUIRED")
            c39_uy_assert(is.list(provenance) && isTRUE(provenance$sha256_verified) && is.character(provenance$sha256) && 
                grepl("^[a-f0-9]{64}$", provenance$sha256), "CALLER_VERIFIED_CACHE_PROVENANCE_REQUIRED")
            audit <- old_y$numerical_audit$death
            final_warnings <- if (isTRUE(audit$retry_attempted)) audit$retry_warnings else audit$original_warnings
            c39_uy_assert(is.list(audit) && !any(grepl("infinite|converg|singular|iteration", unlist(final_warnings), 
                ignore.case = TRUE)), "UNRESOLVED_DEATH_CACHE_WARNING")
            list(fit = old_y$fit$death, reference = old_y$cause_references$death, center_reference = old_y$center_reference, 
                endpoint_id = old_y$endpoint_id, baseline = old_y$engine$baseline_public_validation$death, source_numerical_audit = audit, 
                provenance = provenance, source_method_hash = old_y$method_spec_hash)
        }
        c39_uy_native_or_reuse <- function(work, safe_names, cause, stratified, cluster_bound, control, init = NULL, 
            cache = NULL) {
            if (is.null(cache)) return(ff_ce_native(work, safe_names, cause, stratified, cluster_bound, control, 
                init))
            c39_uy_assert(identical(cause, 2L) && isTRUE(stratified) && isTRUE(cluster_bound) && is.null(init), 
                "DEATH_ONLY_CACHE_REUSE")
            m <- cache$fit
            mf <- stats::model.frame(m)
            c39_uy_assert(inherits(m, "coxph") && !inherits(m, "coxph.penal") && identical(m$method, "efron") && 
                all(is.finite(stats::coef(m))) && all(is.finite(stats::vcov(m))) && !is.null(m$naive.var), "NATIVE_CACHE_ESTIMATOR_INVALID")
            c39_uy_assert(identical(names(stats::coef(m)), safe_names) && identical(cache$reference$safe_names, 
                safe_names) && identical(dim(m$x), dim(as.matrix(work[, safe_names, drop = FALSE]))) && identical(colnames(m$x), 
                safe_names) && identical(as.numeric(m$x), as.numeric(as.matrix(work[, safe_names, drop = FALSE]))) && 
                identical(as.numeric(m$y[, "time"]), as.numeric(work$.ff_time)) && identical(as.numeric(m$y[, "status"]), 
                as.numeric(work$.ff_death)) && identical(as.numeric(stats::model.weights(mf)), as.numeric(work$.ff_weight)) && 
                identical(as.character(stats::model.extract(mf, "cluster")), as.character(work$.ff_cluster)) && 
                identical(as.character(mf[["strata(.ce_center)"]]), as.character(work$.ce_center)) && identical(levels(work$.ce_center), 
                cache$center_reference$observed_levels), "DEATH_CACHE_EXACT_BINDING_MISMATCH")
            c39_uy_assert(identical(attr(m$terms, "term.labels"), c(safe_names, "strata(.ce_center)")) && !length(m$offset) && 
                !length(m$na.action) && is.null(m$frail), "DEATH_CACHE_TERM_DRIFT")
            list(value = m, warnings = character(), cache_reuse = list(reused = TRUE, optimizer_called = FALSE, 
                exact_X_time_status_weight_cluster_center_binding = TRUE, provenance = cache$provenance, source_method_hash = cache$source_method_hash, 
                source_numerical_audit = cache$source_numerical_audit))
        }
        c39_uy_verify_death_baseline <- function(current, cache) {
            if (is.null(cache)) return(invisible(TRUE))
            c39_uy_assert(identical(names(current), names(cache$baseline)) && all(vapply(names(current), function(g) identical(current[[g]]$time, 
                cache$baseline[[g]]$time) && identical(current[[g]]$cumulative, cache$baseline[[g]]$cumulative), 
                FALSE)), "DEATH_CACHE_BASELINE_CHANGED")
            invisible(TRUE)
        }
        c39_uy_install <- function(engine) {
            c39_uy_assert(is.environment(engine) && is.function(engine$ff_ce_fit) && is.function(engine$fd_fp_prepare), 
                "FROZEN_ENGINE_REQUIRED")
            old_spec <- engine$ff_ce_method_spec()
            c39_uy_assert(identical(old_spec$id, "CENTER_ETHNICITY_V4_20261028"), "INSTALL_ONCE_BEFORE_OTHER_METHOD_OVERLAYS")
            helpers <- c("c39_uy_assert", "c39_uy_pool_columns", "c39_uy_native_or_reuse", "c39_uy_verify_death_baseline")
            for (n in helpers) {
                f <- get(n, environment())
                environment(f) <- engine
                engine[[n]] <- f
            }
            spec_text <- paste(c("CORE39_Y_ASIAN_CHINESE_POOL_V1_20261101", "scope=250.2,411;all_B0_to_B100;no_failure_conditional_pooling", 
                "disease=Asian_Chinese_pool+unchanged_center_fixed_covariates", "death=unchanged_Asian_Chinese_pool+unchanged_center_strata", 
                "A+input_ethnicity+sample+draw+selection_weights+Efron+robust_cluster+tau11=UNCHANGED", "death_cache_reuse=accepted_native_only_exact_X_time_status_weights_cluster_center_baseline", 
                "disease_cache_reuse=FORBIDDEN;old_whole_Y_cache=INCOMPATIBLE", paste0("parent_method_sha256=", 
                  old_spec$sha256), "strict_warning_gate_and_one_native_numerical_retry=UNCHANGED"), collapse = "\n")
            method <- list(id = "CORE39_Y_ASIAN_CHINESE_POOL_V1_20261101", canonical = spec_text, sha256 = digest::digest(spec_text, 
                algo = "sha256", serialize = FALSE))
            engine$ff_ce_method_spec <- local({
                value <- method
                function() value
            })
            f <- engine$ff_ce_fit
            f <- c39_uy_replace(f, quote(c("250", "250.2", "411", "411.3", "411.4")), quote(c("250.2", "411")))
            f <- c39_uy_replace(f, quote(if (stratified[1L]) changed else original), quote(ff_y_matrix(c39_uy_pool_columns(original$X))))
            f <- c39_uy_replace(f, quote(got <- ff_ce_native(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), 
                ctrl0)), quote(got <- c39_uy_native_or_reuse(work, dd$reference$safe_names, j, stratified[j], !is.null(cluster), 
                ctrl0, cache = if (j == 2L) death_cache else NULL)))
            f <- c39_uy_replace(f, quote(emit(c(rec, list(stage = "ORIGINAL_ATTEMPT_PRESERVED")))), quote({
                if (!is.null(got$cache_reuse)) {
                  rec$cache_reuse <- got$cache_reuse
                  rec$original_control <- NULL
                }
                emit(c(rec, list(stage = if (is.null(got$cache_reuse)) "ORIGINAL_ATTEMPT_PRESERVED" else "REUSED_ACCEPTED_DEATH_MODEL_PRESERVED")))
            }))
            f <- c39_uy_replace(f, quote(baselines[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, stratified[j], 
                centers$reference$observed_levels, support, j)), quote({
                baselines[[j]] <- ff_ce_baselines(m, dd$reference$safe_names, stratified[j], centers$reference$observed_levels, 
                  support, j)
                if (j == 2L) c39_uy_verify_death_baseline(baselines[[j]], death_cache)
            }))
            formals(f) <- c(formals(f), alist(death_cache = NULL))
            body(f) <- substitute({
                if (!is.null(death_cache)) c39_uy_assert(identical(death_cache$endpoint_id, endpoint_id), "DEATH_CACHE_ENDPOINT_MISMATCH")
                BODY
            }, list(BODY = body(f)))
            environment(f) <- engine
            engine$ff_ce_fit <- f
            p <- engine$ff_ce_predict_at
            p <- c39_uy_replace(p, quote(if (object$stratified[1L]) changed else original$X), quote(c39_uy_pool_columns(original$X)))
            body(p) <- substitute({
                c39_uy_assert(identical(object$method_spec_hash, ff_ce_method_spec()$sha256) && identical(as.logical(object$stratified), 
                  c(FALSE, TRUE)), "INCOMPATIBLE_Y_CACHE")
                BODY
            }, list(BODY = body(p)))
            environment(p) <- engine
            engine$ff_ce_predict_at <- p
            p <- engine$ff_ce_predict_y_matrix
            environment(p) <- engine
            engine$ff_ce_predict_y_matrix <- p
            engine$fd_fp_prepare <- eval(quote(function(X, compiled, v4) {
                o <- compiled$object
                c39_uy_assert(identical(o$method_spec_hash, v4$ff_ce_method_spec()$sha256) && identical(as.logical(o$stratified), 
                  c(FALSE, TRUE)), "INCOMPATIBLE_FAST_Y_CACHE")
                X <- fd_fp_matrix(X, o$matrix_reference)
                centers <- v4$ff_ce_centers(X, o$center_reference)
                changed <- v4$ff_ce_columns(X, o$center_reference$columns)
                list(designs = list(fd_fp_matrix(c39_uy_pool_columns(X), o$cause_references$disease), fd_fp_matrix(changed, 
                  o$cause_references$death)), group = centers$group)
            }), envir = engine)
            invisible(method)
        }
    }), envir = context)
}

