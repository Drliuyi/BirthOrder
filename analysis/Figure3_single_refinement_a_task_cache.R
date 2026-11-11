# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

b16_a_cache_assert <- function(ok, message) {
    if (!isTRUE(ok)) 
        stop(message, call. = FALSE)
}

b16_a_cache_digest <- function(x) {
    digest::digest(serialize(x, NULL, version = 2), algo = "sha256", serialize = FALSE)
}

b16_a_cache_is_hash <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[a-f0-9]{64}$", x)
}

b16_a_parameter_hash <- function(params) {
    if (is.data.frame(params)) 
        return(b16_a_cache_digest(params))
    fields <- c("alpha", "mediator", "mediators", "selection_coefficients", "selection_formula", "base_terms", 
        "factor_levels", "columns", "spline_reference")
    b16_a_cache_digest(stats::setNames(lapply(fields, function(k) params[[k]]), fields))
}

b16_a_cache_function_hash <- function(ctx) {
    native <- attr(ff_a_scenarios, "b16_native_a")
    if (is.null(native)) 
        native <- ff_a_scenarios
    selection <- attr(ff_a_selection, "b16_native_selection")
    if (is.null(selection)) 
        selection <- ff_a_selection
    else {
        b16_a_cache_assert(is.function(selection) && identical(digest::digest(paste(c(deparse(formals(selection), 
            width.cutoff = 500L), deparse(body(selection), width.cutoff = 500L)), collapse = "\n"), algo = "sha256", 
            serialize = FALSE), birth_private_setting("PRIVATE_HASH_1")), "B16_A_CACHE_FROZEN_SELECTION_FUNCTION_CHANGED")
    }
    fns <- list(native, selection, ff_a_category_probability, ctx$engine$b3_missing, ctx$engine$fp_a_path, ctx$engine$b3_design, 
        ctx$discrete$c40_discrete_fit)
    b16_a_cache_digest(lapply(fns, function(f) {
        if (!is.function(f)) 
            return(NULL)
        list(formals = formals(f), body = deparse(body(f), width.cutoff = 500L))
    }))
}

b16_a_cache_binding <- function(ctx, id, bid, pr, d, base, params) {
    ca <- ctx$features[[id]]
    column <- ca$column
    b16_a_cache_assert(is.data.frame(d) && nrow(d) == length(pr) && length(pr) > 0L && identical(as.character(d$eid), 
        as.character(ctx$d_full$eid[pr])), "B16_A_CACHE_ACTUAL_PARENT_ROWS_CHANGED")
    complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
    selection_formula <- if (is.null(params$selection_formula)) 
        ff_a_selection_formula()
    else stats::as.formula(paste(params$selection_formula, collapse = " "))
    vars <- unique(c("eid", base, column, setdiff(all.vars(selection_formula), "complete")))
    b16_a_cache_assert(all(vars %in% names(d)), "B16_A_CACHE_REQUIRED_INPUT_MISSING")
    b16_a_cache_assert(identical(d[, vars, drop = FALSE], ctx$d_full[pr, vars, drop = FALSE]), "B16_A_CACHE_FROZEN_INPUT_VALUES_CHANGED")
    feature_fields <- c("column", "candidate_id", "kind", "prepared_transform", "prepared_unit", "inverse_transform", 
        "physical_lower_bound", "physical_upper_bound", "technical_covariates")
    ans <- list(indicator_id = as.character(id), bootstrap_id = as.character(bid), actual_draw_hash = ctx$actual_draw_hash, 
        target_hash = ff_a_hash(pr), sample_hash = ff_a_hash(pr[complete]), input_lock_sha256 = ctx$input_lock_sha256, 
        source_input_lock_sha256 = ctx$old_input_lock_sha256, source_prepared_sha256 = ctx$old_prepared_inputs_sha256, 
        base_hash = b16_a_cache_digest(as.character(base)), feature_contract_hash = b16_a_cache_digest(stats::setNames(lapply(feature_fields, 
            function(k) ca[[k]]), feature_fields)), factor_levels_hash = b16_a_cache_digest(lapply(d[, vars, drop = FALSE], 
            function(x) if (is.factor(x)) list(levels = levels(x), ordered = is.ordered(x)) else NULL)), a_input_values_hash = b16_a_cache_digest(d[, 
            vars, drop = FALSE]), a_function_hash = b16_a_cache_function_hash(ctx), A_parameter_hash = b16_a_parameter_hash(params))
    b16_a_cache_assert(all(vapply(ans[-c(1L, 2L)], b16_a_cache_is_hash, logical(1))), "B16_A_CACHE_BINDING_HASH_MISSING")
    ans
}

b16_a_cache_proof <- function(desc) {
    p <- desc$a_group_proof_path
    if (is.null(p) && is.null(desc$proof_receipt_sha256)) 
        return(NULL)
    b16_a_cache_assert(is.character(p) && length(p) == 1L && b16_a_cache_is_hash(desc$proof_receipt_sha256) && 
        b16_a_cache_is_hash(desc$a_group_id), "B16_A_CACHE_PROOF_DESCRIPTOR_INCOMPLETE")
    if (exists("ff_p_rebase", mode = "function")) 
        p <- ff_p_rebase(p)
    b16_a_cache_assert(file.exists(p) && file.info(p)$size <= 1024^2, "B16_A_CACHE_PROOF_FILE_MISSING_OR_TOO_LARGE")
    b16_a_cache_assert(identical(ff_a_sha(p), desc$proof_receipt_sha256), "B16_A_CACHE_PROOF_SHA_CHANGED")
    receipt <- jsonlite::fromJSON(p, simplifyVector = FALSE)
    b16_a_cache_assert(identical(ff_a_sha(p), desc$proof_receipt_sha256), "B16_A_CACHE_PROOF_SHA_CHANGED")
    proof <- receipt$task_cache
    if (is.null(proof) || !isTRUE(proof$certified_for_task_cache)) 
        return(NULL)
    b16_a_cache_assert(identical(proof$schema, "B16_A_TASK_CACHE_BINDING_V1") && identical(proof$a_group_id, desc$a_group_id) && 
        is.list(proof$binding), "B16_A_CACHE_PROOF_SCHEMA_OR_GROUP_CHANGED")
    proof
}

b16_a_cache_rng <- function() {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    else NULL
}

b16_a_weight_output_check <- function(ctx, answer, params, proof, expected, actual, d, column) {
    if (is.null(expected) || identical(expected, actual)) 
        return(list(mode = "BITWISE_EQUAL_OUTPUT", max_relative_difference = 0, expected_hash = expected, actual_hash = actual))
    b16_a_cache_assert(identical(proof$weight_proof_kind, "RECONSTRUCTED_SEALED_SELECTION_COEFFICIENTS") && identical(answer$ledger$selection_origin, 
        "RECONSTRUCTED_SEALED_FROZEN_SELECTION_COEFFICIENTS"), "B16_RETAINED_WEIGHT_VECTOR_CHANGED_NOT_ROUNDOFF")
    lock <- jsonlite::fromJSON(file.path(ctx$root, "config/METHOD_LOCK.json"), simplifyVector = FALSE)
    tol <- if (is.null(lock$tolerance$weights)) 
        lock$tolerance
    else lock$tolerance$weights
    b16_a_cache_assert(identical(tol$version, "NUMERIC_OUTPUT_RELATIVE_V1") && identical(as.numeric(tol$relative), 
        1e-12) && identical(as.numeric(tol$absolute), 0), "B16_NUMERIC_TOLERANCE_NOT_LOCKED")
    co <- unlist(params$selection_coefficients, use.names = TRUE)
    b16_a_cache_assert(length(co) > 0L && all(is.finite(co)), "B16_NUMERIC_REFERENCE_COEFFICIENTS_MISSING")
    dd <- d
    dd$complete <- answer$mi$complete
    f <- if (is.null(params$selection_formula)) 
        ff_a_selection_formula()
    else stats::as.formula(paste(params$selection_formula, collapse = " "))
    X <- stats::model.matrix(f, stats::model.frame(f, dd, na.action = stats::na.fail, drop.unused.levels = FALSE))
    if (is.null(names(co)) && length(co) == ncol(X)) 
        names(co) <- colnames(X)
    b16_a_cache_assert(identical(colnames(X), names(co)) && all(is.finite(X)), "B16_NUMERIC_REFERENCE_INPUT_OR_COEFFICIENT_ORDER_CHANGED")
    eta <- rep(0, nrow(X))
    for (j in seq_len(ncol(X))) eta <- eta + X[, j] * co[j]
    reference <- 1/stats::binomial()$linkinv(eta)[answer$mi$complete]
    actual_values <- answer$mi$w
    b16_a_cache_assert(length(reference) == length(actual_values) && all(is.finite(reference) & reference > 0), 
        "B16_NUMERIC_REFERENCE_WEIGHT_INVALID")
    relative <- abs(actual_values - reference)/pmax(abs(actual_values), abs(reference))
    b16_a_cache_assert(all(is.finite(relative)) && all(relative <= tol$relative), "B16_A_WEIGHT_NUMERIC_RELATIVE_TOLERANCE_EXCEEDED")
    list(mode = "SEALED_COEFFICIENT_INDEPENDENT_NUMERIC_RECONSTRUCTION", max_relative_difference = max(relative), 
        relative_tolerance = tol$relative, expected_hash = expected, actual_hash = actual, unretained_historical_vector_claimed_equal = FALSE)
}

b16_a_cache_call <- function(native, ctx, id, bid, ep, pr, d, base) {
    state <- ctx$b16_a_task_cache
    if (!is.environment(state) || !isTRUE(state$active)) 
        return(native(ctx, id, bid, ep, pr, d, base))
    b16_a_cache_assert(identical(state$task, c(as.character(id), as.character(bid))), "B16_A_CACHE_TASK_CHANGED")
    desc <- ff_a_descriptor(ctx, id, as.integer(bid), ep)
    proof <- b16_a_cache_proof(desc)
    if (is.null(proof)) 
        return(native(ctx, id, bid, ep, pr, d, base))
    b16_a_cache_assert(ep %in% as.character(unlist(proof$endpoint_ids)) && identical(as.character(desc$indicator_id), 
        as.character(id)) && identical(as.character(desc$bootstrap_id), as.character(bid)) && identical(as.character(desc$endpoint_id), 
        as.character(ep)), "B16_A_CACHE_DESCRIPTOR_ENDPOINT_CHANGED")
    b16_a_cache_assert(identical(desc$draw_hash, ctx$actual_draw_hash) && identical(desc$target_hash, ff_a_hash(pr)) && 
        identical(desc$sample_hash, ff_a_hash(pr[stats::complete.cases(d[, unique(c(base, ctx$features[[id]]$column)), 
            drop = FALSE])])), "B16_A_CACHE_DRAW_TARGET_OR_SAMPLE_CHANGED")
    params <- ff_a_read(desc, ep, id, bid)
    binding <- b16_a_cache_binding(ctx, id, bid, pr, d, base, params)
    b16_a_cache_assert(identical(desc$source_input_lock_sha256, binding$source_input_lock_sha256) && identical(desc$source_prepared_sha256, 
        binding$source_prepared_sha256), "B16_A_CACHE_DESCRIPTOR_INPUT_LOCK_CHANGED")
    b16_a_cache_assert(identical(names(binding), names(proof$binding)) && all(vapply(names(binding), function(k) identical(binding[[k]], 
        proof$binding[[k]]), logical(1))), "B16_A_CACHE_CERTIFIED_BINDING_CHANGED")
    weight_kind <- proof$weight_proof_kind
    allowed <- c("RETAINED_INDIVIDUAL_VECTOR", "RECONSTRUCTED_SEALED_SELECTION_COEFFICIENTS", "AUDITED_MISSING_FROZEN_REDERIVATION")
    b16_a_cache_assert(weight_kind %in% allowed && b16_a_cache_is_hash(proof$weight_proof_sha256), "B16_A_CACHE_WEIGHT_PROOF_MISSING")
    expected <- proof$reconstructed_weight_hash
    if (identical(weight_kind, "AUDITED_MISSING_FROZEN_REDERIVATION")) {
        b16_a_cache_assert(isTRUE(desc$missing_object) || isTRUE(desc$selection_object_missing) || isTRUE(proof$A_nonreusable_under_unretained_weights), 
            "B16_A_CACHE_MISSING_OBJECT_NOT_AUDITED")
    }
    else b16_a_cache_assert(b16_a_cache_is_hash(expected), "B16_A_CACHE_EXPECTED_WEIGHT_HASH_MISSING")
    history <- desc$selection_weight_hash
    if (is.null(history) || (length(history) == 1L && (is.na(history) || !nzchar(history)))) 
        history <- NULL
    b16_a_cache_assert(identical(history, proof$historical_selection_weight_hash), "B16_A_CACHE_HISTORICAL_WEIGHT_PROOF_CHANGED")
    policy <- list(binding = binding, weight_kind = weight_kind, expected = expected, weight_proof_sha256 = proof$weight_proof_sha256, 
        history = history, allow_recompute = isTRUE(ctx$allow_missing_a_recompute), missing_object = isTRUE(desc$missing_object), 
        selection_object_missing = isTRUE(desc$selection_object_missing), mediator_column = desc$mediator_column)
    key <- desc$a_group_id
    signature <- b16_a_cache_digest(policy)
    hit <- exists(key, state$entries, inherits = FALSE)
    if (hit) {
        item <- get(key, state$entries, inherits = FALSE)
        b16_a_cache_assert(identical(item$signature, signature) && identical(item$proof_sha256, desc$proof_receipt_sha256), 
            "B16_A_CACHE_GROUP_OR_WEIGHT_CHANGED")
        if (!is.null(item$error)) 
            stop(item$error)
        answer <- item$value
    }
    else {
        before <- b16_a_cache_rng()
        result <- tryCatch(list(value = native(ctx, id, bid, ep, pr, d, base)), error = function(e) list(error = e))
        if (!identical(before, b16_a_cache_rng())) {
            error <- simpleError("B16_A_CACHE_NATIVE_RNG_CHANGED_SHARING_REFUSED")
            assign(key, list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, error = error), state$entries)
            stop(error)
        }
        if (!is.null(result$error)) {
            assign(key, list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, error = result$error), 
                state$entries)
            stop(result$error)
        }
        answer <- result$value
        validation <- tryCatch({
            actual_weight <- ff_a_hash(format(answer$mi$w, digits = 17, scientific = TRUE))
            b16_a_cache_assert(is.numeric(answer$mi$w) && all(is.finite(answer$mi$w) & answer$mi$w > 0) && identical(actual_weight, 
                answer$ledger$selection_weight_hash), "B16_A_CACHE_ACTUAL_WEIGHT_CHANGED")
            comparison <- b16_a_weight_output_check(ctx, answer, params, proof, expected, actual_weight, d, ctx$features[[id]]$column)
            answer$ledger$A_weight_output_comparison <- comparison$mode
            answer$ledger$A_weight_output_max_relative_difference <- comparison$max_relative_difference
            answer$ledger$A_weight_certifier_output_hash <- if (is.null(expected)) 
                NA_character_
            else expected
            b16_a_cache_assert(identical(answer$mi$complete, stats::complete.cases(d[, unique(c(base, ctx$features[[id]]$column)), 
                drop = FALSE])) && identical(answer$parent_rows, pr[answer$mi$complete]), "B16_A_CACHE_NATIVE_COMPLETE_SAMPLE_CHANGED")
            NULL
        }, error = function(e) e)
        if (!is.null(validation)) {
            assign(key, list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, error = validation), 
                state$entries)
            stop(validation)
        }
        item <- list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, value = answer, source_endpoint = ep, 
            parameter_endpoint = if (is.null(desc$parameter_endpoint)) ep else desc$parameter_endpoint)
        assign(key, item, state$entries)
    }
    answer$ledger$endpoint_id <- as.character(ep)
    answer$ledger$A_source_endpoint <- as.character(item$parameter_endpoint)
    answer$ledger$A_task_cache <- if (hit) 
        "HIT_CERTIFIED_IDENTICAL_A_GROUP"
    else "MISS_FROZEN_A_CALLED_ONCE"
    answer$ledger$A_task_cache_group <- key
    answer$ledger$A_task_cache_source_endpoint <- as.character(item$source_endpoint)
    answer$ledger$A_task_cache_proof_sha256 <- desc$proof_receipt_sha256
    answer$ledger$A_task_cache_validated_source_path <- if (is.null(desc$path)) 
        NA_character_
    else desc$path
    answer$ledger$A_task_cache_validated_source_sha256 <- if (is.null(desc$sha256)) 
        NA_character_
    else desc$sha256
    answer
}

b16_install_a_task_cache <- function(envir = parent.frame()) {
    native_a <- get("ff_a_scenarios", envir = envir)
    native_compute <- get("ff_compute", envir = envir)
    b16_a_cache_assert(is.null(attr(native_a, "b16_native_a")), "B16_A_CACHE_ALREADY_INSTALLED")
    wrapped_a <- function(ctx, id, bid, ep, pr, d, base) b16_a_cache_call(native_a, ctx, id, bid, ep, pr, d, base)
    attr(wrapped_a, "b16_native_a") <- native_a
    wrapped_compute <- function(ctx, command) {
        state <- new.env(parent = emptyenv())
        state$entries <- new.env(parent = emptyenv())
        state$active <- TRUE
        state$task <- c(as.character(command$indicator_id), as.character(command$bootstrap_id))
        ctx$b16_a_task_cache <- state
        on.exit({
            state$active <- FALSE
            rm(list = ls(state$entries, all.names = TRUE), envir = state$entries)
        }, add = TRUE)
        native_compute(ctx, command)
    }
    assign("ff_a_scenarios", wrapped_a, envir = envir)
    assign("ff_compute", wrapped_compute, envir = envir)
    invisible(TRUE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        b16_a_cache_assert <- function(ok, message) {
            if (!isTRUE(ok)) stop(message, call. = FALSE)
        }
        b16_a_cache_digest <- function(x) {
            digest::digest(serialize(x, NULL, version = 2), algo = "sha256", serialize = FALSE)
        }
        b16_a_cache_is_hash <- function(x) {
            is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[a-f0-9]{64}$", x)
        }
        b16_a_parameter_hash <- function(params) {
            if (is.data.frame(params)) return(b16_a_cache_digest(params))
            fields <- c("alpha", "mediator", "mediators", "selection_coefficients", "selection_formula", "base_terms", 
                "factor_levels", "columns", "spline_reference")
            b16_a_cache_digest(stats::setNames(lapply(fields, function(k) params[[k]]), fields))
        }
        b16_a_cache_function_hash <- function(ctx) {
            native <- attr(ff_a_scenarios, "b16_native_a")
            if (is.null(native)) native <- ff_a_scenarios
            selection <- attr(ff_a_selection, "b16_native_selection")
            if (is.null(selection)) selection <- ff_a_selection else {
                b16_a_cache_assert(is.function(selection) && identical(digest::digest(paste(c(deparse(formals(selection), 
                  width.cutoff = 500L), deparse(body(selection), width.cutoff = 500L)), collapse = "\n"), algo = "sha256", 
                  serialize = FALSE), birth_private_setting("PRIVATE_HASH_2")), "B16_A_CACHE_FROZEN_SELECTION_FUNCTION_CHANGED")
            }
            fns <- list(native, selection, ff_a_category_probability, ctx$engine$b3_missing, ctx$engine$fp_a_path, 
                ctx$engine$b3_design, ctx$discrete$c40_discrete_fit)
            b16_a_cache_digest(lapply(fns, function(f) {
                if (!is.function(f)) return(NULL)
                list(formals = formals(f), body = deparse(body(f), width.cutoff = 500L))
            }))
        }
        b16_a_cache_binding <- function(ctx, id, bid, pr, d, base, params) {
            ca <- ctx$features[[id]]
            column <- ca$column
            b16_a_cache_assert(is.data.frame(d) && nrow(d) == length(pr) && length(pr) > 0L && identical(as.character(d$eid), 
                as.character(ctx$d_full$eid[pr])), "B16_A_CACHE_ACTUAL_PARENT_ROWS_CHANGED")
            complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
            selection_formula <- if (is.null(params$selection_formula)) ff_a_selection_formula() else stats::as.formula(paste(params$selection_formula, 
                collapse = " "))
            vars <- unique(c("eid", base, column, setdiff(all.vars(selection_formula), "complete")))
            b16_a_cache_assert(all(vars %in% names(d)), "B16_A_CACHE_REQUIRED_INPUT_MISSING")
            b16_a_cache_assert(identical(d[, vars, drop = FALSE], ctx$d_full[pr, vars, drop = FALSE]), "B16_A_CACHE_FROZEN_INPUT_VALUES_CHANGED")
            feature_fields <- c("column", "candidate_id", "kind", "prepared_transform", "prepared_unit", "inverse_transform", 
                "physical_lower_bound", "physical_upper_bound", "technical_covariates")
            ans <- list(indicator_id = as.character(id), bootstrap_id = as.character(bid), actual_draw_hash = ctx$actual_draw_hash, 
                target_hash = ff_a_hash(pr), sample_hash = ff_a_hash(pr[complete]), input_lock_sha256 = ctx$input_lock_sha256, 
                source_input_lock_sha256 = ctx$old_input_lock_sha256, source_prepared_sha256 = ctx$old_prepared_inputs_sha256, 
                base_hash = b16_a_cache_digest(as.character(base)), feature_contract_hash = b16_a_cache_digest(stats::setNames(lapply(feature_fields, 
                  function(k) ca[[k]]), feature_fields)), factor_levels_hash = b16_a_cache_digest(lapply(d[, vars, 
                  drop = FALSE], function(x) if (is.factor(x)) list(levels = levels(x), ordered = is.ordered(x)) else NULL)), 
                a_input_values_hash = b16_a_cache_digest(d[, vars, drop = FALSE]), a_function_hash = b16_a_cache_function_hash(ctx), 
                A_parameter_hash = b16_a_parameter_hash(params))
            b16_a_cache_assert(all(vapply(ans[-c(1L, 2L)], b16_a_cache_is_hash, logical(1))), "B16_A_CACHE_BINDING_HASH_MISSING")
            ans
        }
        b16_a_cache_proof <- function(desc) {
            p <- desc$a_group_proof_path
            if (is.null(p) && is.null(desc$proof_receipt_sha256)) return(NULL)
            b16_a_cache_assert(is.character(p) && length(p) == 1L && b16_a_cache_is_hash(desc$proof_receipt_sha256) && 
                b16_a_cache_is_hash(desc$a_group_id), "B16_A_CACHE_PROOF_DESCRIPTOR_INCOMPLETE")
            if (exists("ff_p_rebase", mode = "function")) p <- ff_p_rebase(p)
            b16_a_cache_assert(file.exists(p) && file.info(p)$size <= 1024^2, "B16_A_CACHE_PROOF_FILE_MISSING_OR_TOO_LARGE")
            b16_a_cache_assert(identical(ff_a_sha(p), desc$proof_receipt_sha256), "B16_A_CACHE_PROOF_SHA_CHANGED")
            receipt <- jsonlite::fromJSON(p, simplifyVector = FALSE)
            b16_a_cache_assert(identical(ff_a_sha(p), desc$proof_receipt_sha256), "B16_A_CACHE_PROOF_SHA_CHANGED")
            proof <- receipt$task_cache
            if (is.null(proof) || !isTRUE(proof$certified_for_task_cache)) return(NULL)
            b16_a_cache_assert(identical(proof$schema, "B16_A_TASK_CACHE_BINDING_V1") && identical(proof$a_group_id, 
                desc$a_group_id) && is.list(proof$binding), "B16_A_CACHE_PROOF_SCHEMA_OR_GROUP_CHANGED")
            proof
        }
        b16_a_cache_rng <- function() {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, 
                inherits = FALSE) else NULL
        }
        b16_a_weight_output_check <- function(ctx, answer, params, proof, expected, actual, d, column) {
            if (is.null(expected) || identical(expected, actual)) return(list(mode = "BITWISE_EQUAL_OUTPUT", max_relative_difference = 0, 
                expected_hash = expected, actual_hash = actual))
            b16_a_cache_assert(identical(proof$weight_proof_kind, "RECONSTRUCTED_SEALED_SELECTION_COEFFICIENTS") && 
                identical(answer$ledger$selection_origin, "RECONSTRUCTED_SEALED_FROZEN_SELECTION_COEFFICIENTS"), 
                "B16_RETAINED_WEIGHT_VECTOR_CHANGED_NOT_ROUNDOFF")
            lock <- jsonlite::fromJSON(file.path(ctx$root, "config/METHOD_LOCK.json"), simplifyVector = FALSE)
            tol <- if (is.null(lock$tolerance$weights)) lock$tolerance else lock$tolerance$weights
            b16_a_cache_assert(identical(tol$version, "NUMERIC_OUTPUT_RELATIVE_V1") && identical(as.numeric(tol$relative), 
                1e-12) && identical(as.numeric(tol$absolute), 0), "B16_NUMERIC_TOLERANCE_NOT_LOCKED")
            co <- unlist(params$selection_coefficients, use.names = TRUE)
            b16_a_cache_assert(length(co) > 0L && all(is.finite(co)), "B16_NUMERIC_REFERENCE_COEFFICIENTS_MISSING")
            dd <- d
            dd$complete <- answer$mi$complete
            f <- if (is.null(params$selection_formula)) ff_a_selection_formula() else stats::as.formula(paste(params$selection_formula, 
                collapse = " "))
            X <- stats::model.matrix(f, stats::model.frame(f, dd, na.action = stats::na.fail, drop.unused.levels = FALSE))
            if (is.null(names(co)) && length(co) == ncol(X)) names(co) <- colnames(X)
            b16_a_cache_assert(identical(colnames(X), names(co)) && all(is.finite(X)), "B16_NUMERIC_REFERENCE_INPUT_OR_COEFFICIENT_ORDER_CHANGED")
            eta <- rep(0, nrow(X))
            for (j in seq_len(ncol(X))) eta <- eta + X[, j] * co[j]
            reference <- 1/stats::binomial()$linkinv(eta)[answer$mi$complete]
            actual_values <- answer$mi$w
            b16_a_cache_assert(length(reference) == length(actual_values) && all(is.finite(reference) & reference > 
                0), "B16_NUMERIC_REFERENCE_WEIGHT_INVALID")
            relative <- abs(actual_values - reference)/pmax(abs(actual_values), abs(reference))
            b16_a_cache_assert(all(is.finite(relative)) && all(relative <= tol$relative), "B16_A_WEIGHT_NUMERIC_RELATIVE_TOLERANCE_EXCEEDED")
            list(mode = "SEALED_COEFFICIENT_INDEPENDENT_NUMERIC_RECONSTRUCTION", max_relative_difference = max(relative), 
                relative_tolerance = tol$relative, expected_hash = expected, actual_hash = actual, unretained_historical_vector_claimed_equal = FALSE)
        }
        b16_a_cache_call <- function(native, ctx, id, bid, ep, pr, d, base) {
            state <- ctx$b16_a_task_cache
            if (!is.environment(state) || !isTRUE(state$active)) return(native(ctx, id, bid, ep, pr, d, base))
            b16_a_cache_assert(identical(state$task, c(as.character(id), as.character(bid))), "B16_A_CACHE_TASK_CHANGED")
            desc <- ff_a_descriptor(ctx, id, as.integer(bid), ep)
            proof <- b16_a_cache_proof(desc)
            if (is.null(proof)) return(native(ctx, id, bid, ep, pr, d, base))
            b16_a_cache_assert(ep %in% as.character(unlist(proof$endpoint_ids)) && identical(as.character(desc$indicator_id), 
                as.character(id)) && identical(as.character(desc$bootstrap_id), as.character(bid)) && identical(as.character(desc$endpoint_id), 
                as.character(ep)), "B16_A_CACHE_DESCRIPTOR_ENDPOINT_CHANGED")
            b16_a_cache_assert(identical(desc$draw_hash, ctx$actual_draw_hash) && identical(desc$target_hash, ff_a_hash(pr)) && 
                identical(desc$sample_hash, ff_a_hash(pr[stats::complete.cases(d[, unique(c(base, ctx$features[[id]]$column)), 
                  drop = FALSE])])), "B16_A_CACHE_DRAW_TARGET_OR_SAMPLE_CHANGED")
            params <- ff_a_read(desc, ep, id, bid)
            binding <- b16_a_cache_binding(ctx, id, bid, pr, d, base, params)
            b16_a_cache_assert(identical(desc$source_input_lock_sha256, binding$source_input_lock_sha256) && identical(desc$source_prepared_sha256, 
                binding$source_prepared_sha256), "B16_A_CACHE_DESCRIPTOR_INPUT_LOCK_CHANGED")
            b16_a_cache_assert(identical(names(binding), names(proof$binding)) && all(vapply(names(binding), function(k) identical(binding[[k]], 
                proof$binding[[k]]), logical(1))), "B16_A_CACHE_CERTIFIED_BINDING_CHANGED")
            weight_kind <- proof$weight_proof_kind
            allowed <- c("RETAINED_INDIVIDUAL_VECTOR", "RECONSTRUCTED_SEALED_SELECTION_COEFFICIENTS", "AUDITED_MISSING_FROZEN_REDERIVATION")
            b16_a_cache_assert(weight_kind %in% allowed && b16_a_cache_is_hash(proof$weight_proof_sha256), "B16_A_CACHE_WEIGHT_PROOF_MISSING")
            expected <- proof$reconstructed_weight_hash
            if (identical(weight_kind, "AUDITED_MISSING_FROZEN_REDERIVATION")) {
                b16_a_cache_assert(isTRUE(desc$missing_object) || isTRUE(desc$selection_object_missing) || isTRUE(proof$A_nonreusable_under_unretained_weights), 
                  "B16_A_CACHE_MISSING_OBJECT_NOT_AUDITED")
            } else b16_a_cache_assert(b16_a_cache_is_hash(expected), "B16_A_CACHE_EXPECTED_WEIGHT_HASH_MISSING")
            history <- desc$selection_weight_hash
            if (is.null(history) || (length(history) == 1L && (is.na(history) || !nzchar(history)))) history <- NULL
            b16_a_cache_assert(identical(history, proof$historical_selection_weight_hash), "B16_A_CACHE_HISTORICAL_WEIGHT_PROOF_CHANGED")
            policy <- list(binding = binding, weight_kind = weight_kind, expected = expected, weight_proof_sha256 = proof$weight_proof_sha256, 
                history = history, allow_recompute = isTRUE(ctx$allow_missing_a_recompute), missing_object = isTRUE(desc$missing_object), 
                selection_object_missing = isTRUE(desc$selection_object_missing), mediator_column = desc$mediator_column)
            key <- desc$a_group_id
            signature <- b16_a_cache_digest(policy)
            hit <- exists(key, state$entries, inherits = FALSE)
            if (hit) {
                item <- get(key, state$entries, inherits = FALSE)
                b16_a_cache_assert(identical(item$signature, signature) && identical(item$proof_sha256, desc$proof_receipt_sha256), 
                  "B16_A_CACHE_GROUP_OR_WEIGHT_CHANGED")
                if (!is.null(item$error)) stop(item$error)
                answer <- item$value
            } else {
                before <- b16_a_cache_rng()
                result <- tryCatch(list(value = native(ctx, id, bid, ep, pr, d, base)), error = function(e) list(error = e))
                if (!identical(before, b16_a_cache_rng())) {
                  error <- simpleError("B16_A_CACHE_NATIVE_RNG_CHANGED_SHARING_REFUSED")
                  assign(key, list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, error = error), 
                    state$entries)
                  stop(error)
                }
                if (!is.null(result$error)) {
                  assign(key, list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, error = result$error), 
                    state$entries)
                  stop(result$error)
                }
                answer <- result$value
                validation <- tryCatch({
                  actual_weight <- ff_a_hash(format(answer$mi$w, digits = 17, scientific = TRUE))
                  b16_a_cache_assert(is.numeric(answer$mi$w) && all(is.finite(answer$mi$w) & answer$mi$w > 0) && 
                    identical(actual_weight, answer$ledger$selection_weight_hash), "B16_A_CACHE_ACTUAL_WEIGHT_CHANGED")
                  comparison <- b16_a_weight_output_check(ctx, answer, params, proof, expected, actual_weight, 
                    d, ctx$features[[id]]$column)
                  answer$ledger$A_weight_output_comparison <- comparison$mode
                  answer$ledger$A_weight_output_max_relative_difference <- comparison$max_relative_difference
                  answer$ledger$A_weight_certifier_output_hash <- if (is.null(expected)) NA_character_ else expected
                  b16_a_cache_assert(identical(answer$mi$complete, stats::complete.cases(d[, unique(c(base, ctx$features[[id]]$column)), 
                    drop = FALSE])) && identical(answer$parent_rows, pr[answer$mi$complete]), "B16_A_CACHE_NATIVE_COMPLETE_SAMPLE_CHANGED")
                  NULL
                }, error = function(e) e)
                if (!is.null(validation)) {
                  assign(key, list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, error = validation), 
                    state$entries)
                  stop(validation)
                }
                item <- list(signature = signature, proof_sha256 = desc$proof_receipt_sha256, value = answer, source_endpoint = ep, 
                  parameter_endpoint = if (is.null(desc$parameter_endpoint)) ep else desc$parameter_endpoint)
                assign(key, item, state$entries)
            }
            answer$ledger$endpoint_id <- as.character(ep)
            answer$ledger$A_source_endpoint <- as.character(item$parameter_endpoint)
            answer$ledger$A_task_cache <- if (hit) "HIT_CERTIFIED_IDENTICAL_A_GROUP" else "MISS_FROZEN_A_CALLED_ONCE"
            answer$ledger$A_task_cache_group <- key
            answer$ledger$A_task_cache_source_endpoint <- as.character(item$source_endpoint)
            answer$ledger$A_task_cache_proof_sha256 <- desc$proof_receipt_sha256
            answer$ledger$A_task_cache_validated_source_path <- if (is.null(desc$path)) NA_character_ else desc$path
            answer$ledger$A_task_cache_validated_source_sha256 <- if (is.null(desc$sha256)) NA_character_ else desc$sha256
            answer
        }
        b16_install_a_task_cache <- function(envir = parent.frame()) {
            native_a <- get("ff_a_scenarios", envir = envir)
            native_compute <- get("ff_compute", envir = envir)
            b16_a_cache_assert(is.null(attr(native_a, "b16_native_a")), "B16_A_CACHE_ALREADY_INSTALLED")
            wrapped_a <- function(ctx, id, bid, ep, pr, d, base) b16_a_cache_call(native_a, ctx, id, bid, ep, pr, 
                d, base)
            attr(wrapped_a, "b16_native_a") <- native_a
            wrapped_compute <- function(ctx, command) {
                state <- new.env(parent = emptyenv())
                state$entries <- new.env(parent = emptyenv())
                state$active <- TRUE
                state$task <- c(as.character(command$indicator_id), as.character(command$bootstrap_id))
                ctx$b16_a_task_cache <- state
                on.exit({
                  state$active <- FALSE
                  rm(list = ls(state$entries, all.names = TRUE), envir = state$entries)
                }, add = TRUE)
                native_compute(ctx, command)
            }
            assign("ff_a_scenarios", wrapped_a, envir = envir)
            assign("ff_compute", wrapped_compute, envir = envir)
            invisible(TRUE)
        }
        if (exists("ff_a_scenarios", mode = "function") && exists("ff_compute", mode = "function")) b16_install_a_task_cache(environment())
    }), envir = context)
}

