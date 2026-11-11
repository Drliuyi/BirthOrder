# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

fd_fp_assert <- function(ok, code) {
    if (!isTRUE(ok)) 
        stop(paste0("FD_FAST_", code), call. = FALSE)
    invisible(TRUE)
}

fd_fp_matrix <- function(X, reference) {
    fd_fp_assert(is.matrix(X) && is.numeric(X) && nrow(X) >= 1L && all(is.finite(X)) && !is.null(colnames(X)) && 
        !anyNA(colnames(X)) && !anyDuplicated(colnames(X)), "INVALID_DESIGN")
    if ("(Intercept)" %in% colnames(X)) {
        fd_fp_assert(all(X[, "(Intercept)"] == 1), "INVALID_INTERCEPT")
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    fd_fp_assert(identical(colnames(X), reference$columns) && all(c("B", "age") %in% colnames(X)) && all(X[, "B"] >= 
        0 & X[, "B"] <= 1), "EXACT_DESIGN_BINDING_FAILED")
    X
}

fd_fp_compile <- function(object, v4) {
    fd_fp_assert(inherits(object, "ff_ce_matrix_fit") && identical(as.numeric(object$tau_report), 11) && identical(object$ties, 
        "efron") && object$endpoint_id %in% c("250.2", "411") && identical(object$method_spec_hash, v4$ff_ce_method_spec()$sha256) && 
        isTRUE(object$native_validation$passed), "UNVERIFIED_FROZEN_V4_FIT")
    fd_fp_assert(identical(names(object$fit), c("disease", "death")) && identical(as.logical(object$stratified), 
        c(FALSE, TRUE)), "CAUSE_STRATA_DRIFT")
    beta <- lapply(seq_along(object$fit), function(j) {
        f <- object$fit[[j]]
        r <- object$cause_references[[j]]
        expected <- c(r$safe_names, if (object$stratified[j]) "strata(.ce_center)")
        fd_fp_assert(inherits(f, "coxph") && !inherits(f, "coxphms") && !inherits(f, "coxph.penal") && identical(f$method, 
            "efron") && identical(names(stats::coef(f)), r$safe_names) && length(r$columns) == length(r$safe_names) && 
            !anyDuplicated(r$columns) && !anyDuplicated(r$safe_names) && all(is.finite(stats::coef(f))) && identical(attr(f$terms, 
            "term.labels"), expected) && all(attr(f$terms, "order") == 1L) && !length(attr(f$terms, "offset")) && 
            !length(f$offset) && !length(f$na.action) && is.null(f$frail), "COEFFICIENT_TERM_OR_OFFSET_MAPPING_UNSUPPORTED")
        stats::coef(f)
    })
    levels <- object$center_reference$observed_levels
    fd_fp_assert(identical(names(object$engine$by_center), levels), "BASELINE_CENTER_BINDING")
    hazards <- lapply(levels, function(g) {
        support <- object$center_support[object$center_support$center_design_label == g, , drop = FALSE]
        fd_fp_assert(nrow(support) == 1L && support$max_followup >= 11, "CENTER_TAU11_UNSUPPORTED")
        eng <- object$engine$by_center[[g]]
        fd_fp_assert(is.numeric(eng$time) && all(is.finite(eng$time)) && !is.unsorted(eng$time, strictly = TRUE) && 
            is.matrix(eng$hazard) && ncol(eng$hazard) == 2L && nrow(eng$hazard) == length(eng$time) && all(is.finite(eng$hazard) & 
            eng$hazard >= 0), "INVALID_NATIVE_BASELINE_ENGINE")
        eng$hazard[eng$time <= 11, , drop = FALSE]
    })
    names(hazards) <- levels
    list(object = object, beta = beta, hazards = hazards)
}

fd_fp_prepare <- function(X, compiled, v4) {
    o <- compiled$object
    X <- fd_fp_matrix(X, o$matrix_reference)
    centers <- v4$ff_ce_centers(X, o$center_reference)
    changed <- v4$ff_ce_columns(X, o$center_reference$columns)
    list(designs = list(fd_fp_matrix(X, o$cause_references$disease), fd_fp_matrix(changed, o$cause_references$death)), 
        group = centers$group)
}

fd_fp_predict_matrix <- function(compiled, X, v4, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
    64)) {
    fd_fp_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1L && length(max_prediction_mb) == 
        1L && is.finite(max_prediction_mb) && max_prediction_mb > 0, "INVALID_MEMORY_BUDGET")
    prepared <- fd_fp_prepare(X, compiled, v4)
    step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
    p <- matrix(NA_real_, nrow(X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
    for (g in unique(prepared$group)) {
        rows <- which(prepared$group == g)
        h <- compiled$hazards[[g]]
        for (start in seq.int(1L, length(rows), by = step)) {
            ii <- rows[seq.int(start, min(length(rows), start + step - 1L))]
            lp <- vapply(1:2, function(j) as.numeric(prepared$designs[[j]][ii, , drop = FALSE] %*% compiled$beta[[j]]), 
                numeric(length(ii)))
            rr <- exp(matrix(lp, nrow = length(ii), ncol = 2L))
            fd_fp_assert(all(is.finite(rr)) && all(rr > 0), "INVALID_RR_NO_CLIPPING")
            if (is.function(v4$ff_aj_status) && isTRUE(v4$ff_aj_status()$ready)) {
                p[ii, ] <- v4$ff_aj_backend(h, rr)
            }
            else {
                fd_fp_assert(!isTRUE(getOption("ff.y.require_compiled_aj", FALSE)), "COMPILED_AJ_REQUIRED")
                p[ii, ] <- v4$ff_y_aj(h, cbind(rr, rr), rep(1L, nrow(h)))
            }
            v4$ff_y_probability_check(p[ii, , drop = FALSE])
        }
    }
    v4$ff_y_probability_check(p)
    attr(p, "warnings") <- character()
    attr(p, "chunk_size") <- step
    attr(p, "product.limit") <- TRUE
    attr(p, "risk_range") <- apply(p, 2L, range)
    p
}

fd_attach_fast_prediction <- function(z, ctx, validation_data = NULL, cache_entries = 3L, lp_tolerance = 1e-10, 
    cif_tolerance = 1e-12) {
    fd_fp_assert(is.list(z) && length(z$fits) > 0L && all(names(z$fits) %in% c("250.2", "411")), "MISSING_FROZEN_FITS")
    fd_fp_assert(length(cache_entries) == 1L && is.finite(cache_entries) && cache_entries >= 2L && cache_entries == 
        as.integer(cache_entries), "INVALID_CACHE_BOUND")
    fd_fp_assert(length(lp_tolerance) == 1L && is.finite(lp_tolerance) && lp_tolerance >= 0 && lp_tolerance <= 
        1e-10 && length(cif_tolerance) == 1L && is.finite(cif_tolerance) && cif_tolerance >= 0 && cif_tolerance <= 
        1e-12, "TOLERANCE_MAY_ONLY_BE_TIGHTENED")
    v4 <- ctx$v4
    members <- z$pop$members
    base <- z$pop$base
    reference <- z$reference
    fits <- z$fits
    compiled <- lapply(fits, fd_fp_compile, v4 = v4)
    design_fun <- fd_design
    native_fun <- v4$ff_ce_predict_y_matrix
    original <- function(d, endpoint) {
        if (is.null(fits[[endpoint]])) 
            return(rep(0, nrow(d)))
        X <- design_fun(d, members, base, reference)$X
        native_fun(fits[[endpoint]], X, 11)[, "disease"]
    }
    discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
    design_terms <- stats::terms(stats::reformulate(unique(c(base, discrete))))
    design_columns <- unique(c(all.vars(design_terms), vapply(members, `[[`, "", "column")))
    labels <- attr(design_terms, "term.labels")
    b_reuse <- "B" %in% labels && all(vapply(setdiff(labels, "B"), function(term) !"B" %in% all.vars(stats::reformulate(term)), 
        logical(1L))) && !"B" %in% vapply(members, `[[`, "", "column")
    state <- new.env(parent = emptyenv())
    state$cache <- list()
    state$ready <- FALSE
    state$stats <- list(design_builds = 0L, design_cache_hits = 0L, b_column_reuses = 0L, predictions = 0L, validations = 0L)
    state$validation <- list()
    get_design <- function(d) {
        fd_fp_assert(is.data.frame(d) && nrow(d) >= 1L && all(design_columns %in% names(d)), "INPUT_COLUMNS")
        reuse <- b_reuse && is.numeric(d$B) && all(is.finite(d$B) & d$B >= 0 & d$B <= 1)
        keys <- if (reuse) 
            setdiff(design_columns, "B")
        else design_columns
        key <- d[, keys, drop = FALSE]
        hit <- which(vapply(state$cache, function(entry) identical(entry$reuse, reuse) && identical(entry$key, 
            key), logical(1L)))
        if (length(hit)) {
            entry <- state$cache[[hit[1L]]]
            state$stats$design_cache_hits <- state$stats$design_cache_hits + 1L
            X <- entry$X
            if (reuse && !identical(as.numeric(X[, "B"]), as.numeric(d$B))) {
                X[, "B"] <- d$B
                state$stats$b_column_reuses <- state$stats$b_column_reuses + 1L
            }
            return(X)
        }
        X <- design_fun(d, members, base, reference)$X
        if (reuse) 
            fd_fp_assert("B" %in% colnames(X) && identical(as.numeric(X[, "B"]), as.numeric(d$B)), "B_COLUMN_NOT_IDENTITY")
        state$stats$design_builds <- state$stats$design_builds + 1L
        state$cache <- c(list(list(key = key, reuse = reuse, X = X)), state$cache)
        if (length(state$cache) > cache_entries) 
            state$cache <- state$cache[seq_len(cache_entries)]
        X
    }
    validate <- function(data, label = "EXPLICIT_SAME_INPUT_VALIDATION") {
        state$ready <- FALSE
        samples <- if (is.data.frame(data)) 
            list(data)
        else data
        fd_fp_assert(is.list(samples) && length(samples) >= 1L, "VALIDATION_INPUT_REQUIRED")
        lp_error <- cif_error <- design_error <- 0
        checked_rows <- 0L
        for (d in samples) {
            X <- get_design(d)
            expected_X <- design_fun(d, members, base, reference)$X
            fd_fp_assert(identical(dimnames(X), dimnames(expected_X)) && identical(dim(X), dim(expected_X)), "CACHED_DESIGN_DIMENSION_DRIFT")
            design_error <- max(design_error, abs(X - expected_X))
            fd_fp_assert(design_error == 0, "CACHED_DESIGN_NOT_EXACT")
            for (endpoint in names(compiled)) {
                cfit <- compiled[[endpoint]]
                prepared <- fd_fp_prepare(X, cfit, v4)
                for (g in unique(prepared$group)) for (j in 1:2) {
                  ii <- which(prepared$group == g)
                  nd <- as.data.frame(prepared$designs[[j]][ii, , drop = FALSE])
                  names(nd) <- cfit$object$cause_references[[j]]$safe_names
                  nd$.ce_center <- factor(g, levels = cfit$object$center_reference$observed_levels)
                  native <- v4$ff_y_capture(as.numeric(stats::predict(cfit$object$fit[[j]], newdata = nd, type = "lp", 
                    reference = "zero")))
                  fd_fp_assert(!length(native$warnings), "NATIVE_LP_VALIDATION_WARNING")
                  fast <- as.numeric(prepared$designs[[j]][ii, , drop = FALSE] %*% cfit$beta[[j]])
                  lp_error <- max(lp_error, abs(native$value - fast))
                }
                fast <- fd_fp_predict_matrix(cfit, X, v4)
                native <- native_fun(fits[[endpoint]], expected_X, 11)
                fd_fp_assert(!length(attr(native, "warnings")), "NATIVE_CIF_VALIDATION_WARNING")
                cif_error <- max(cif_error, abs(native - fast))
                fd_fp_assert(is.finite(lp_error) && lp_error <= lp_tolerance && is.finite(cif_error) && cif_error <= 
                  cif_tolerance, "NATIVE_EQUIVALENCE_FAILED")
            }
            checked_rows <- checked_rows + nrow(d)
        }
        receipt <- list(label = label, rows = checked_rows, endpoints = names(compiled), design_max_abs = design_error, 
            lp_max_abs = lp_error, cif_all_states_max_abs = cif_error, lp_tolerance = lp_tolerance, cif_tolerance = cif_tolerance, 
            passed = TRUE, original_fits_baselines_weights_unchanged = TRUE, consumes_RNG = FALSE)
        state$validation[[length(state$validation) + 1L]] <- receipt
        state$stats$validations <- state$stats$validations + 1L
        state$ready <- TRUE
        receipt
    }
    if (is.null(validation_data)) {
        d <- z$pop$data
        fd_fp_assert(is.data.frame(d) && nrow(d) > 0L && "center" %in% names(d), "DEFAULT_VALIDATION_DATA_REQUIRED")
        rows <- sort(unique(c(match(unique(d$center), d$center), as.integer(round(seq(1L, nrow(d), length.out = min(32L, 
            nrow(d))))))))
        probe <- d[rows, , drop = FALSE]
        zero <- one <- probe
        zero$B <- 0
        one$B <- 1
        validation_data <- list(observed = probe, R00 = zero, R10 = one)
    }
    validate(validation_data, "ATTACHMENT_SAME_INPUT_VALIDATION")
    state$cache <- list()
    state$stats$design_builds <- state$stats$design_cache_hits <- state$stats$b_column_reuses <- 0L
    z$predict_fun_original <- original
    z$predict_fun <- function(d, endpoint) {
        fd_fp_assert(isTRUE(state$ready), "PREDICTOR_NOT_VALIDATED")
        fd_fp_assert(length(endpoint) == 1L && endpoint %in% c("250.2", "411"), "ENDPOINT_SCOPE")
        if (is.null(compiled[[endpoint]])) 
            return(rep(0, nrow(d)))
        state$stats$predictions <- state$stats$predictions + 1L
        fd_fp_predict_matrix(compiled[[endpoint]], get_design(d), v4)[, "disease"]
    }
    z$validate_prediction <- validate
    z$prediction_validation <- function() state$validation
    z$prediction_cache_stats <- function() c(state$stats, list(cache_entries = length(state$cache), cache_limit = cache_entries, 
        numeric_additive_B_reuse = b_reuse, cache_bytes = as.numeric(object.size(state$cache))))
    z
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        fd_fp_assert <- function(ok, code) {
            if (!isTRUE(ok)) stop(paste0("FD_FAST_", code), call. = FALSE)
            invisible(TRUE)
        }
        fd_fp_matrix <- function(X, reference) {
            fd_fp_assert(is.matrix(X) && is.numeric(X) && nrow(X) >= 1L && all(is.finite(X)) && !is.null(colnames(X)) && 
                !anyNA(colnames(X)) && !anyDuplicated(colnames(X)), "INVALID_DESIGN")
            if ("(Intercept)" %in% colnames(X)) {
                fd_fp_assert(all(X[, "(Intercept)"] == 1), "INVALID_INTERCEPT")
                X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            }
            fd_fp_assert(identical(colnames(X), reference$columns) && all(c("B", "age") %in% colnames(X)) && all(X[, 
                "B"] >= 0 & X[, "B"] <= 1), "EXACT_DESIGN_BINDING_FAILED")
            X
        }
        fd_fp_compile <- function(object, v4) {
            fd_fp_assert(inherits(object, "ff_ce_matrix_fit") && identical(as.numeric(object$tau_report), 11) && 
                identical(object$ties, "efron") && object$endpoint_id %in% c("250.2", "411") && identical(object$method_spec_hash, 
                v4$ff_ce_method_spec()$sha256) && isTRUE(object$native_validation$passed), "UNVERIFIED_FROZEN_V4_FIT")
            fd_fp_assert(identical(names(object$fit), c("disease", "death")) && identical(as.logical(object$stratified), 
                c(FALSE, TRUE)), "CAUSE_STRATA_DRIFT")
            beta <- lapply(seq_along(object$fit), function(j) {
                f <- object$fit[[j]]
                r <- object$cause_references[[j]]
                expected <- c(r$safe_names, if (object$stratified[j]) "strata(.ce_center)")
                fd_fp_assert(inherits(f, "coxph") && !inherits(f, "coxphms") && !inherits(f, "coxph.penal") && 
                  identical(f$method, "efron") && identical(names(stats::coef(f)), r$safe_names) && length(r$columns) == 
                  length(r$safe_names) && !anyDuplicated(r$columns) && !anyDuplicated(r$safe_names) && all(is.finite(stats::coef(f))) && 
                  identical(attr(f$terms, "term.labels"), expected) && all(attr(f$terms, "order") == 1L) && !length(attr(f$terms, 
                  "offset")) && !length(f$offset) && !length(f$na.action) && is.null(f$frail), "COEFFICIENT_TERM_OR_OFFSET_MAPPING_UNSUPPORTED")
                stats::coef(f)
            })
            levels <- object$center_reference$observed_levels
            fd_fp_assert(identical(names(object$engine$by_center), levels), "BASELINE_CENTER_BINDING")
            hazards <- lapply(levels, function(g) {
                support <- object$center_support[object$center_support$center_design_label == g, , drop = FALSE]
                fd_fp_assert(nrow(support) == 1L && support$max_followup >= 11, "CENTER_TAU11_UNSUPPORTED")
                eng <- object$engine$by_center[[g]]
                fd_fp_assert(is.numeric(eng$time) && all(is.finite(eng$time)) && !is.unsorted(eng$time, strictly = TRUE) && 
                  is.matrix(eng$hazard) && ncol(eng$hazard) == 2L && nrow(eng$hazard) == length(eng$time) && all(is.finite(eng$hazard) & 
                  eng$hazard >= 0), "INVALID_NATIVE_BASELINE_ENGINE")
                eng$hazard[eng$time <= 11, , drop = FALSE]
            })
            names(hazards) <- levels
            list(object = object, beta = beta, hazards = hazards)
        }
        fd_fp_prepare <- function(X, compiled, v4) {
            o <- compiled$object
            X <- fd_fp_matrix(X, o$matrix_reference)
            centers <- v4$ff_ce_centers(X, o$center_reference)
            changed <- v4$ff_ce_columns(X, o$center_reference$columns)
            list(designs = list(fd_fp_matrix(X, o$cause_references$disease), fd_fp_matrix(changed, o$cause_references$death)), 
                group = centers$group)
        }
        fd_fp_predict_matrix <- function(compiled, X, v4, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
            64)) {
            fd_fp_assert(length(chunk_size) == 1L && is.finite(chunk_size) && chunk_size >= 1L && length(max_prediction_mb) == 
                1L && is.finite(max_prediction_mb) && max_prediction_mb > 0, "INVALID_MEMORY_BUDGET")
            prepared <- fd_fp_prepare(X, compiled, v4)
            step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
            p <- matrix(NA_real_, nrow(X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
            for (g in unique(prepared$group)) {
                rows <- which(prepared$group == g)
                h <- compiled$hazards[[g]]
                for (start in seq.int(1L, length(rows), by = step)) {
                  ii <- rows[seq.int(start, min(length(rows), start + step - 1L))]
                  lp <- vapply(1:2, function(j) as.numeric(prepared$designs[[j]][ii, , drop = FALSE] %*% compiled$beta[[j]]), 
                    numeric(length(ii)))
                  rr <- exp(matrix(lp, nrow = length(ii), ncol = 2L))
                  fd_fp_assert(all(is.finite(rr)) && all(rr > 0), "INVALID_RR_NO_CLIPPING")
                  if (is.function(v4$ff_aj_status) && isTRUE(v4$ff_aj_status()$ready)) {
                    p[ii, ] <- v4$ff_aj_backend(h, rr)
                  } else {
                    fd_fp_assert(!isTRUE(getOption("ff.y.require_compiled_aj", FALSE)), "COMPILED_AJ_REQUIRED")
                    p[ii, ] <- v4$ff_y_aj(h, cbind(rr, rr), rep(1L, nrow(h)))
                  }
                  v4$ff_y_probability_check(p[ii, , drop = FALSE])
                }
            }
            v4$ff_y_probability_check(p)
            attr(p, "warnings") <- character()
            attr(p, "chunk_size") <- step
            attr(p, "product.limit") <- TRUE
            attr(p, "risk_range") <- apply(p, 2L, range)
            p
        }
        fd_attach_fast_prediction <- function(z, ctx, validation_data = NULL, cache_entries = 3L, lp_tolerance = 1e-10, 
            cif_tolerance = 1e-12) {
            fd_fp_assert(is.list(z) && length(z$fits) > 0L && all(names(z$fits) %in% c("250.2", "411")), "MISSING_FROZEN_FITS")
            fd_fp_assert(length(cache_entries) == 1L && is.finite(cache_entries) && cache_entries >= 2L && cache_entries == 
                as.integer(cache_entries), "INVALID_CACHE_BOUND")
            fd_fp_assert(length(lp_tolerance) == 1L && is.finite(lp_tolerance) && lp_tolerance >= 0 && lp_tolerance <= 
                1e-10 && length(cif_tolerance) == 1L && is.finite(cif_tolerance) && cif_tolerance >= 0 && cif_tolerance <= 
                1e-12, "TOLERANCE_MAY_ONLY_BE_TIGHTENED")
            v4 <- ctx$v4
            members <- z$pop$members
            base <- z$pop$base
            reference <- z$reference
            fits <- z$fits
            compiled <- lapply(fits, fd_fp_compile, v4 = v4)
            design_fun <- fd_design
            native_fun <- v4$ff_ce_predict_y_matrix
            original <- function(d, endpoint) {
                if (is.null(fits[[endpoint]])) return(rep(0, nrow(d)))
                X <- design_fun(d, members, base, reference)$X
                native_fun(fits[[endpoint]], X, 11)[, "disease"]
            }
            discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
            design_terms <- stats::terms(stats::reformulate(unique(c(base, discrete))))
            design_columns <- unique(c(all.vars(design_terms), vapply(members, `[[`, "", "column")))
            labels <- attr(design_terms, "term.labels")
            b_reuse <- "B" %in% labels && all(vapply(setdiff(labels, "B"), function(term) !"B" %in% all.vars(stats::reformulate(term)), 
                logical(1L))) && !"B" %in% vapply(members, `[[`, "", "column")
            state <- new.env(parent = emptyenv())
            state$cache <- list()
            state$ready <- FALSE
            state$stats <- list(design_builds = 0L, design_cache_hits = 0L, b_column_reuses = 0L, predictions = 0L, 
                validations = 0L)
            state$validation <- list()
            get_design <- function(d) {
                fd_fp_assert(is.data.frame(d) && nrow(d) >= 1L && all(design_columns %in% names(d)), "INPUT_COLUMNS")
                reuse <- b_reuse && is.numeric(d$B) && all(is.finite(d$B) & d$B >= 0 & d$B <= 1)
                keys <- if (reuse) setdiff(design_columns, "B") else design_columns
                key <- d[, keys, drop = FALSE]
                hit <- which(vapply(state$cache, function(entry) identical(entry$reuse, reuse) && identical(entry$key, 
                  key), logical(1L)))
                if (length(hit)) {
                  entry <- state$cache[[hit[1L]]]
                  state$stats$design_cache_hits <- state$stats$design_cache_hits + 1L
                  X <- entry$X
                  if (reuse && !identical(as.numeric(X[, "B"]), as.numeric(d$B))) {
                    X[, "B"] <- d$B
                    state$stats$b_column_reuses <- state$stats$b_column_reuses + 1L
                  }
                  return(X)
                }
                X <- design_fun(d, members, base, reference)$X
                if (reuse) fd_fp_assert("B" %in% colnames(X) && identical(as.numeric(X[, "B"]), as.numeric(d$B)), 
                  "B_COLUMN_NOT_IDENTITY")
                state$stats$design_builds <- state$stats$design_builds + 1L
                state$cache <- c(list(list(key = key, reuse = reuse, X = X)), state$cache)
                if (length(state$cache) > cache_entries) state$cache <- state$cache[seq_len(cache_entries)]
                X
            }
            validate <- function(data, label = "EXPLICIT_SAME_INPUT_VALIDATION") {
                state$ready <- FALSE
                samples <- if (is.data.frame(data)) list(data) else data
                fd_fp_assert(is.list(samples) && length(samples) >= 1L, "VALIDATION_INPUT_REQUIRED")
                lp_error <- cif_error <- design_error <- 0
                checked_rows <- 0L
                for (d in samples) {
                  X <- get_design(d)
                  expected_X <- design_fun(d, members, base, reference)$X
                  fd_fp_assert(identical(dimnames(X), dimnames(expected_X)) && identical(dim(X), dim(expected_X)), 
                    "CACHED_DESIGN_DIMENSION_DRIFT")
                  design_error <- max(design_error, abs(X - expected_X))
                  fd_fp_assert(design_error == 0, "CACHED_DESIGN_NOT_EXACT")
                  for (endpoint in names(compiled)) {
                    cfit <- compiled[[endpoint]]
                    prepared <- fd_fp_prepare(X, cfit, v4)
                    for (g in unique(prepared$group)) for (j in 1:2) {
                      ii <- which(prepared$group == g)
                      nd <- as.data.frame(prepared$designs[[j]][ii, , drop = FALSE])
                      names(nd) <- cfit$object$cause_references[[j]]$safe_names
                      nd$.ce_center <- factor(g, levels = cfit$object$center_reference$observed_levels)
                      native <- v4$ff_y_capture(as.numeric(stats::predict(cfit$object$fit[[j]], newdata = nd, type = "lp", 
                        reference = "zero")))
                      fd_fp_assert(!length(native$warnings), "NATIVE_LP_VALIDATION_WARNING")
                      fast <- as.numeric(prepared$designs[[j]][ii, , drop = FALSE] %*% cfit$beta[[j]])
                      lp_error <- max(lp_error, abs(native$value - fast))
                    }
                    fast <- fd_fp_predict_matrix(cfit, X, v4)
                    native <- native_fun(fits[[endpoint]], expected_X, 11)
                    fd_fp_assert(!length(attr(native, "warnings")), "NATIVE_CIF_VALIDATION_WARNING")
                    cif_error <- max(cif_error, abs(native - fast))
                    fd_fp_assert(is.finite(lp_error) && lp_error <= lp_tolerance && is.finite(cif_error) && cif_error <= 
                      cif_tolerance, "NATIVE_EQUIVALENCE_FAILED")
                  }
                  checked_rows <- checked_rows + nrow(d)
                }
                receipt <- list(label = label, rows = checked_rows, endpoints = names(compiled), design_max_abs = design_error, 
                  lp_max_abs = lp_error, cif_all_states_max_abs = cif_error, lp_tolerance = lp_tolerance, cif_tolerance = cif_tolerance, 
                  passed = TRUE, original_fits_baselines_weights_unchanged = TRUE, consumes_RNG = FALSE)
                state$validation[[length(state$validation) + 1L]] <- receipt
                state$stats$validations <- state$stats$validations + 1L
                state$ready <- TRUE
                receipt
            }
            if (is.null(validation_data)) {
                d <- z$pop$data
                fd_fp_assert(is.data.frame(d) && nrow(d) > 0L && "center" %in% names(d), "DEFAULT_VALIDATION_DATA_REQUIRED")
                rows <- sort(unique(c(match(unique(d$center), d$center), as.integer(round(seq(1L, nrow(d), length.out = min(32L, 
                  nrow(d))))))))
                probe <- d[rows, , drop = FALSE]
                zero <- one <- probe
                zero$B <- 0
                one$B <- 1
                validation_data <- list(observed = probe, R00 = zero, R10 = one)
            }
            validate(validation_data, "ATTACHMENT_SAME_INPUT_VALIDATION")
            state$cache <- list()
            state$stats$design_builds <- state$stats$design_cache_hits <- state$stats$b_column_reuses <- 0L
            z$predict_fun_original <- original
            z$predict_fun <- function(d, endpoint) {
                fd_fp_assert(isTRUE(state$ready), "PREDICTOR_NOT_VALIDATED")
                fd_fp_assert(length(endpoint) == 1L && endpoint %in% c("250.2", "411"), "ENDPOINT_SCOPE")
                if (is.null(compiled[[endpoint]])) return(rep(0, nrow(d)))
                state$stats$predictions <- state$stats$predictions + 1L
                fd_fp_predict_matrix(compiled[[endpoint]], get_design(d), v4)[, "disease"]
            }
            z$validate_prediction <- validate
            z$prediction_validation <- function() state$validation
            z$prediction_cache_stats <- function() c(state$stats, list(cache_entries = length(state$cache), cache_limit = cache_entries, 
                numeric_additive_B_reuse = b_reuse, cache_bytes = as.numeric(object.size(state$cache))))
            z
        }
    }), envir = context)
}

