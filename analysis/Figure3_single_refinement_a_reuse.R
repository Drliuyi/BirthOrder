# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ff_a_assert <- function(ok, code) if (!isTRUE(ok)) stop(code, call. = FALSE)

ff_a_hash <- function(x) digest::digest(paste0(nchar(as.character(x)), ":", as.character(x), "\n", collapse = ""), 
    algo = "sha256", serialize = FALSE)

ff_a_sha <- function(p) {
    ff_a_assert(file.exists(p) && !dir.exists(p), paste("FF_A_SOURCE_FILE_MISSING", p))
    digest::digest(file = p, algo = "sha256")
}

ff_a_check <- function(p, h) {
    ff_a_assert(is.character(h) && length(h) == 1L && grepl("^[a-f0-9]{64}$", h), "FF_A_EXACT_SOURCE_SHA256_REQUIRED")
    ff_a_assert(identical(ff_a_sha(p), h), paste("FF_A_SOURCE_SHA256_CHANGED", p))
}

ff_a_descriptor <- function(ctx, id, bid, endpoint) {
    if (is.function(ctx$a_source_for)) 
        out <- ctx$a_source_for(id, bid, endpoint)
    else {
        key <- paste(id, bid, endpoint, sep = "|")
        index <- ctx$a_reuse_index
        if (is.environment(index)) 
            out <- get0(key, index, inherits = FALSE)
        else out <- index[[key]]
    }
    ff_a_assert(is.list(out), paste("FF_A_SOURCE_NOT_AUDITED", id, bid, endpoint))
    ff_a_assert(!is.null(out$format) || (isTRUE(out$missing_object) && nzchar(out$reason)), "FF_A_SOURCE_DESCRIPTOR_INCOMPLETE")
    out
}

ff_a_read <- function(desc, endpoint, id = NULL, bid = NULL) {
    if (isTRUE(desc$missing_object)) 
        return(NULL)
    parameter_endpoint <- if (is.null(desc$parameter_endpoint)) 
        endpoint
    else desc$parameter_endpoint
    if (identical(desc$format, "c40_shard_gzip")) {
        ff_a_check(desc$path, desc$sha256)
        ff_a_assert(file.info(desc$path)$size <= 16 * 1024^2, "FF_A_C40_SHARD_COMPRESSED_LIMIT")
        con <- gzfile(desc$path, "rb")
        limit <- 64L * 1024L * 1024L
        bytes <- tryCatch(readBin(con, "raw", n = limit + 1L), finally = close(con))
        ff_a_assert(length(bytes) <= limit, "FF_A_C40_SHARD_DECOMPRESSED_LIMIT")
        shard <- jsonlite::fromJSON(rawToChar(bytes), simplifyVector = FALSE)
        ff_a_assert(identical(shard$schema, "C40_FIVE_ENDPOINT_AGGREGATE_SHARD_V2") && is.list(shard$bundles) && 
            length(shard$bundles) >= 1L && length(shard$bundles) <= 25L, "FF_A_C40_SHARD_SCHEMA_CHANGED")
        units <- Filter(function(u) identical(u$indicator_id, id) && identical(as.integer(u$bootstrap_id), as.integer(bid)), 
            shard$bundles)
        ff_a_assert(length(units) == 1L, "FF_A_C40_NUMBERED_UNIT_MISSING_OR_DUPLICATED")
        unit <- units[[1L]]
        ff_a_assert(identical(unit$draw$draw_hash, desc$draw_hash), "FF_A_C40_ACTUAL_DRAW_CHANGED")
        rows <- Filter(function(r) identical(as.character(r$endpoint_id), as.character(endpoint)), unit$rows)
        same_or_audited_missing <- function(value, expected) {
            if (is.character(value) && length(value) == 1L && !is.na(value) && grepl("^[a-f0-9]{64}$", value)) 
                identical(value, expected)
            else isTRUE(desc$binding_reconstruction_required)
        }
        ff_a_assert(length(rows) == 1L && same_or_audited_missing(rows[[1L]]$target_hash, desc$target_hash) && 
            same_or_audited_missing(rows[[1L]]$sample_hash, desc$sample_hash), "FF_A_C40_TARGET_SAMPLE_PROOF_CHANGED")
        ff_a_assert(is.list(unit$files) && "PARAMETERS.json" %in% names(unit$files), "FF_A_C40_SEALED_PARAMETER_ARTIFACT_ABSENT")
        z <- unit$files[["PARAMETERS.json"]][[as.character(parameter_endpoint)]]
        ff_a_check(desc$path, desc$sha256)
        return(z)
    }
    else if (!is.null(desc$archive)) {
        ff_a_check(desc$archive, desc$archive_sha256)
        member <- desc$member
        ff_a_assert(is.character(member) && length(member) == 1L && nzchar(member) && !grepl("(^/|(^|/)[.][.](/|$))", 
            member), "FF_A_UNSAFE_ARCHIVE_MEMBER")
        ff_a_assert(is.character(desc$member_sha256) && grepl("^[a-f0-9]{64}$", desc$member_sha256), "FF_A_MEMBER_SHA256_REQUIRED")
        limit <- 2L * 1024L * 1024L
        con <- pipe(paste("tar -xOf", shQuote(desc$archive), shQuote(member)), "rb")
        bytes <- tryCatch(readBin(con, "raw", n = limit + 1L), finally = close(con))
        ff_a_assert(length(bytes) <= limit, "FF_A_AGGREGATE_CHECKPOINT_EXCEEDS_2_MIB")
        ff_a_assert(identical(digest::digest(bytes, algo = "sha256", serialize = FALSE), desc$member_sha256), "FF_A_ARCHIVE_MEMBER_CHANGED")
        ff_a_assert(identical(desc$format, "c40_json"), "FF_A_ONLY_JSON_ARCHIVE_CHECKPOINT_SUPPORTED")
        z <- jsonlite::fromJSON(rawToChar(bytes), simplifyVector = TRUE)
        ff_a_check(desc$archive, desc$archive_sha256)
    }
    else {
        ff_a_check(desc$path, desc$sha256)
        z <- switch(desc$format, v5_rds = readRDS(desc$path), b18_rds = readRDS(desc$path), c40_json = jsonlite::fromJSON(desc$path, 
            simplifyVector = TRUE), legacy_coefficients_csv = read.csv(desc$path, stringsAsFactors = FALSE, check.names = FALSE), 
            stop("FF_A_UNKNOWN_CHECKPOINT_FORMAT", call. = FALSE))
        ff_a_check(desc$path, desc$sha256)
    }
    if (desc$format %in% c("v5_rds", "c40_json")) 
        z <- z[[as.character(parameter_endpoint)]]
    z
}

ff_a_selection_formula <- function() stats::reformulate(c("B", "splines::ns(age, df=3)", "sex", "center", "assessment_year", 
    "S", "splines::ns(aux_t_250_2, knots=c(1,3), Boundary.knots=c(0,5))", "aux_event_250_2", "aux_death_250_2", 
    "splines::ns(aux_t_411, knots=c(1,3), Boundary.knots=c(0,5))", "aux_event_411", "aux_death_411"), response = "complete")

ff_a_selection <- function(ctx, d, column, base, params, desc) {
    complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
    ff_a_assert(any(complete), "NO_COMPLETE_RECORDS")
    coefficients <- params$selection_coefficients
    A_refit_required <- FALSE
    coefficient_max_abs_difference <- NA_real_
    if (all(complete)) {
        ff_a_assert(is.null(coefficients) || !length(coefficients), "FF_A_ALL_COMPLETE_CONFLICTS_WITH_SAVED_SELECTION_MODEL")
        mi <- list(complete = complete, w = rep(1, nrow(d)), fit = NULL, warnings = character(), diagnostics = data.frame())
        origin <- "FROZEN_ALL_COMPLETE_UNIT_SELECTION_WEIGHTS"
    }
    else if (!is.null(coefficients) && length(coefficients)) {
        coefficients <- unlist(coefficients, use.names = TRUE)
        dd <- d
        dd$complete <- complete
        f <- if (is.null(params$selection_formula)) 
            ff_a_selection_formula()
        else stats::as.formula(paste(params$selection_formula, collapse = " "))
        X <- stats::model.matrix(f, stats::model.frame(f, dd, na.action = stats::na.fail, drop.unused.levels = FALSE))
        if (is.null(names(coefficients)) && length(coefficients) == ncol(X)) 
            names(coefficients) <- colnames(X)
        ff_a_assert(identical(colnames(X), names(coefficients)) && all(is.finite(coefficients)), "FF_A_SAVED_SELECTION_DESIGN_OR_COEFFICIENT_DRIFT")
        pi <- stats::binomial()$linkinv(as.numeric(X %*% coefficients))
        ff_a_assert(all(is.finite(pi) & pi > 0 & pi <= 1), "FF_A_SELECTION_POSITIVITY_FAILURE")
        mi <- list(complete = complete, w = 1/pi[complete], fit = NULL, warnings = character(), diagnostics = data.frame(), 
            saved_coefficients = coefficients)
        origin <- "RECONSTRUCTED_SEALED_FROZEN_SELECTION_COEFFICIENTS"
    }
    else {
        missing_confirmed <- isTRUE(desc$selection_object_missing) || isTRUE(desc$missing_object) || identical(desc$format, 
            "legacy_coefficients_csv") || is.null(params)
        ff_a_assert(missing_confirmed, "FF_A_SELECTION_CHECKPOINT_MISSING_NOT_AUDITED")
        ff_a_assert(isTRUE(ctx$allow_missing_a_recompute), "FF_A_MISSING_SELECTION_RECOMPUTE_NOT_AUTHORIZED")
        ff_a_assert(is.function(ctx$engine$b3_missing), "FF_A_FROZEN_SELECTION_FUNCTION_MISSING")
        mi <- ctx$engine$b3_missing(d, column, base)
        ff_a_assert(identical(mi$complete, complete), "FF_A_RECOMPUTED_SELECTION_SAMPLE_CHANGED")
        origin <- "RECOMPUTED_FROZEN_SELECTION_ONLY_PERSISTED_OBJECT_MISSING"
        A_refit_required <- TRUE
        origin <- "A_NONREUSABLE_SELECTION_FIT_AND_WEIGHT_VECTOR_NOT_RETAINED"
    }
    ff_a_assert(length(mi$w) == sum(complete) && all(is.finite(mi$w) & mi$w > 0), "FF_A_INVALID_SELECTION_WEIGHT_VECTOR")
    hash <- ff_a_hash(format(mi$w, digits = 17, scientific = TRUE))
    if (!is.null(desc$selection_weight_hash) && !is.na(desc$selection_weight_hash) && nzchar(desc$selection_weight_hash) && 
        !identical(hash, desc$selection_weight_hash)) {
        reason <- "A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED"
        ff_a_assert(!isTRUE(desc$individual_weight_vector_retained), "FF_A_ACTUALLY_RETAINED_SELECTION_WEIGHT_VECTOR_MISMATCH_INPUT_DRIFT")
        ff_a_assert(isTRUE(ctx$allow_missing_a_recompute) && is.function(ctx$engine$b3_missing), paste("FF_A_SELECTION_WEIGHT_HASH_MISMATCH", 
            reason))
        rebuilt <- ctx$engine$b3_missing(d, column, base)
        ff_a_assert(identical(rebuilt$complete, complete), paste("FF_A_REBUILT_SELECTION_SAMPLE_CHANGED", reason))
        ff_a_assert(!is.null(rebuilt$fit) && identical(attr(stats::terms(rebuilt$fit), "term.labels"), attr(stats::terms(ff_a_selection_formula()), 
            "term.labels")), "FF_A_REBUILT_SELECTION_FROZEN_FORMULA_CHANGED")
        if (!is.null(coefficients) && length(coefficients)) {
            newcoef <- stats::coef(rebuilt$fit)
            ff_a_assert(identical(names(newcoef), names(coefficients)) && all(is.finite(newcoef)), "FF_A_REBUILT_SELECTION_COEFFICIENT_DESIGN_CHANGED")
            coefficient_max_abs_difference <- max(abs(as.numeric(newcoef) - as.numeric(coefficients)))
        }
        rebuilt_hash <- ff_a_hash(format(rebuilt$w, digits = 17, scientific = TRUE))
        mi <- rebuilt
        hash <- rebuilt_hash
        A_refit_required <- TRUE
        origin <- reason
    }
    list(mi = mi, origin = origin, weight_hash = hash, old_weight_hash = if (is.null(desc$selection_weight_hash)) NA_character_ else desc$selection_weight_hash, 
        A_refit_required = A_refit_required, saved_selection_coefficient_max_abs_difference = coefficient_max_abs_difference)
}

ff_a_numeric <- function(x, expected_names = NULL) {
    ans <- unlist(x, use.names = TRUE)
    if (!is.null(expected_names)) {
        ff_a_assert(length(ans) == length(expected_names), "FF_A_PARAMETER_LENGTH_CHANGED")
        if (is.null(names(ans))) 
            names(ans) <- expected_names
        ff_a_assert(identical(names(ans), expected_names), "FF_A_PARAMETER_COLUMNS_CHANGED")
    }
    ff_a_assert(is.numeric(ans) && all(is.finite(ans)), "FF_A_NONFINITE_SAVED_PARAMETER")
    ans
}

ff_a_category_probability <- function(X, saved, levels) {
    mu <- ff_a_numeric(saved$conditioning_mean, colnames(X))
    scale <- ff_a_numeric(saved$conditioning_scale, colnames(X))
    ff_a_assert(all(scale > 0), "FF_A_INVALID_SAVED_CONDITIONING_SCALE")
    Z <- sweep(sweep(X, 2L, mu, "-"), 2L, scale, "/")
    coef <- saved$coefficients
    if (identical(saved$estimator, "WEIGHTED_BINOMIAL_ML_GLM_FIT")) {
        ff_a_assert(length(levels) == 2L, "FF_A_BINARY_LEVEL_DRIFT")
        beta <- ff_a_numeric(coef, colnames(X))
        p <- stats::plogis(as.numeric(Z %*% beta))
        ans <- cbind(1 - p, p)
    }
    else {
        ff_a_assert(identical(saved$estimator, "WEIGHTED_NOMINAL_MULTINOMIAL_ML_NNET_DECAY_ZERO"), "FF_A_UNKNOWN_SAVED_MEDIATOR_ESTIMATOR")
        beta <- if (is.matrix(coef)) 
            coef
        else if (is.data.frame(coef)) 
            as.matrix(coef)
        else if (is.list(coef) && length(coef) == length(levels) - 1L) 
            do.call(rbind, lapply(coef, unlist))
        else as.numeric(unlist(coef))
        if (length(levels) == 2L && is.null(dim(beta))) 
            beta <- matrix(beta, nrow = 1L)
        ff_a_assert(is.matrix(beta) && identical(dim(beta), c(length(levels) - 1L, ncol(X))) && all(is.finite(beta)), 
            "FF_A_MULTINOMIAL_COEFFICIENT_DIMENSION_DRIFT")
        eta <- cbind(0, Z %*% t(beta))
        eta <- eta - apply(eta, 1L, max)
        p <- exp(eta)
        ans <- p/rowSums(p)
    }
    colnames(ans) <- levels
    ff_a_assert(all(is.finite(ans) & ans >= 0 & ans <= 1) && max(abs(rowSums(ans) - 1)) < 1e-10, "FF_A_INVALID_CATEGORY_PROBABILITIES_NO_CLIPPING")
    ans
}

ff_a_scenarios <- function(ctx, id, bid, ep, pr, d, base) {
    ff_a_assert(is.numeric(bid) && length(bid) == 1L && bid == as.integer(bid) && bid %in% 0:1000, "FF_A_FIXED_ID_ZERO_TO_1000_REQUIRED")
    ff_a_assert(is.data.frame(d) && length(pr) == nrow(d) && length(pr) > 0L && identical(as.character(d$eid), 
        as.character(ctx$d_full$eid[pr])), "FF_A_ACTUAL_PARENT_DRAW_ROW_BINDING_MISMATCH")
    ca <- ctx$features[[id]]
    ff_a_assert(is.list(ca) && ca$column %in% names(d), "FF_A_FROZEN_FEATURE_MISSING")
    column <- ca$column
    desc <- ff_a_descriptor(ctx, id, as.integer(bid), ep)
    target_hash <- ff_a_hash(pr)
    binding_reconstructed <- FALSE
    if (isTRUE(desc$binding_reconstruction_required)) {
        ff_a_assert(identical(desc$source_input_lock_sha256, ctx$old_input_lock_sha256) && identical(desc$source_prepared_sha256, 
            ctx$old_prepared_inputs_sha256), "FF_A_MISSING_BINDING_NOT_BOUND_TO_FROZEN_OLD_INPUTS")
        missing_hash <- function(h) is.null(h) || length(h) != 1L || is.na(h) || !is.character(h) || !grepl("^[a-f0-9]{64}$", 
            h)
        if (missing_hash(desc$target_hash)) 
            desc$target_hash <- target_hash
        if (missing_hash(desc$sample_hash)) 
            desc$sample_hash <- ff_a_hash(pr[stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])])
        binding_reconstructed <- TRUE
    }
    ff_a_assert(is.character(desc$target_hash) && identical(desc$target_hash, target_hash), "FF_A_ACTUAL_TARGET_HASH_CHANGED")
    complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
    sample_hash <- ff_a_hash(pr[complete])
    ff_a_assert(is.character(desc$sample_hash) && identical(desc$sample_hash, sample_hash), "FF_A_ACTUAL_COMPLETE_SAMPLE_HASH_CHANGED")
    ff_a_assert(is.character(desc$draw_hash) && length(desc$draw_hash) == 1L && grepl("^[a-f0-9]{64}$", desc$draw_hash), 
        "FF_A_ACTUAL_FULL_PARENT_DRAW_PROOF_REQUIRED")
    ff_a_assert(is.character(ctx$actual_draw_hash) && length(ctx$actual_draw_hash) == 1L && identical(desc$draw_hash, 
        ctx$actual_draw_hash), "FF_A_ACTUAL_DRAW_HASH_CHANGED_OR_NOT_BOUND")
    params <- ff_a_read(desc, ep, id, bid)
    if (!is.null(params$base_terms)) 
        ff_a_assert(identical(as.character(unlist(params$base_terms)), as.character(base)), "FF_A_SAVED_BASE_FORMULA_CHANGED")
    if (!is.null(params$factor_levels)) 
        for (name in names(params$factor_levels)) ff_a_assert(name %in% names(d) && is.factor(d[[name]]) && identical(levels(d[[name]]), 
            as.character(unlist(params$factor_levels[[name]]))), "FF_A_SAVED_FACTOR_LEVEL_DICTIONARY_CHANGED")
    selection <- ff_a_selection(ctx, d, column, base, params, desc)
    mi <- selection$mi
    dd <- d[mi$complete, , drop = FALSE]
    factor_mediator <- is.factor(dd[[column]])
    ledger <- data.frame(indicator_id = id, endpoint_id = ep, bootstrap_id = as.integer(bid), source_path = if (!is.null(desc$archive)) 
        desc$archive
    else if (is.null(desc$path)) 
        NA_character_
    else desc$path, source_member = if (is.null(desc$member)) 
        NA_character_
    else desc$member, source_sha256 = if (!is.null(desc$archive)) 
        desc$archive_sha256
    else if (is.null(desc$sha256)) 
        NA_character_
    else desc$sha256, source_format = if (is.null(desc$format)) 
        "AUDITED_MISSING"
    else desc$format, A_source_endpoint = if (is.null(desc$parameter_endpoint)) 
        ep
    else desc$parameter_endpoint, source_missing_reason = if (is.null(desc$reason)) 
        ""
    else desc$reason, draw_hash = desc$draw_hash, target_hash = target_hash, sample_hash = sample_hash, draw_proof_kind = if (is.null(desc$draw_proof_kind)) 
        "AUDITED_ACTUAL_DRAW_BINDING"
    else desc$draw_proof_kind, draw_artifact = if (is.null(desc$draw_artifact)) 
        NA_character_
    else desc$draw_artifact, old_failed_cell_missing_sample_binding_reconstructed = binding_reconstructed, selection_weight_hash = selection$weight_hash, 
        selection_origin = selection$origin, old_selection_weight_hash = selection$old_weight_hash, selection_and_A_rebuilt_not_claimed_bitwise_old_weights = selection$A_refit_required, 
        saved_selection_coefficient_max_abs_difference = selection$saved_selection_coefficient_max_abs_difference, 
        N_target = nrow(d), N_complete = nrow(dd), A_origin = NA_character_, simulation_rule = NA_character_, reason = "", 
        stringsAsFactors = FALSE)
    out <- list(mi = mi, data = dd, column = column, parent_rows = pr[mi$complete], ledger = ledger, reference_basis = NULL, 
        alpha = NULL, scenarios = NULL, g0 = NULL, g1 = NULL, levels = NULL)
    if (!factor_mediator) {
        alpha <- if (isTRUE(selection$A_refit_required)) 
            NULL
        else if (is.data.frame(params) && identical(desc$format, "legacy_coefficients_csv")) {
            term <- paste0(if (is.null(desc$mediator_column)) 
                column
            else desc$mediator_column, "__a_path")
            q <- params[params$term == term, , drop = FALSE]
            if (nrow(q) == 1L) 
                stats::setNames(q$coefficient, column)
            else NULL
        }
        else params$alpha
        if (!is.null(alpha) && length(alpha)) {
            alpha <- unlist(alpha, use.names = TRUE)
            old_column <- if (is.null(desc$mediator_column)) 
                column
            else desc$mediator_column
            ff_a_assert(old_column %in% names(alpha), "FF_A_SAVED_ALPHA_COLUMN_MISSING")
            alpha <- stats::setNames(as.numeric(alpha[[old_column]]), column)
            ff_a_assert(is.finite(alpha), "FF_A_SAVED_ALPHA_NONFINITE")
            out$ledger$A_origin <- "REUSED_NUMBER_SPECIFIC_SAVED_ALPHA"
        }
        else {
            ff_a_assert(isTRUE(ctx$allow_missing_a_recompute), "FF_A_MISSING_ALPHA_RECOMPUTE_NOT_AUTHORIZED")
            ff_a_assert(is.function(ctx$engine$fp_a_path), "FF_A_FROZEN_A_FUNCTION_MISSING")
            aa <- ctx$engine$fp_a_path(d, column, base, mi, ctx, paste0("single_", ca$candidate_id), as.integer(bid))
            alpha <- aa$alpha
            out$ledger$A_origin <- "RECOMPUTED_NUMBER_SPECIFIC_A_OBJECT_MISSING"
            out$ledger$reason <- if (isTRUE(selection$A_refit_required)) 
                selection$origin
            else "OLD_CHECKPOINT_HAS_NO_ALPHA_FOR_THIS_NUMBER_AND_SAMPLE"
        }
        m0 <- dd[[column]] - dd$B * alpha[[column]]
        m1 <- m0 + alpha[[column]]
        lo <- if (is.null(ca$physical_lower_bound)) 
            -Inf
        else ca$physical_lower_bound
        hi <- if (is.null(ca$physical_upper_bound)) 
            Inf
        else ca$physical_upper_bound
        ff_a_assert(all(is.finite(m0)) && all(is.finite(m1)) && all(m0 >= lo & m0 <= hi) && all(m1 >= lo & m1 <= 
            hi), "SCENARIO_OUTSIDE_PHYSICAL_DOMAIN_NO_CLIPPING")
        reference <- params$spline_reference
        if (!is.null(reference) && !column %in% names(reference) && !is.null(desc$mediator_column) && desc$mediator_column %in% 
            names(reference)) 
            names(reference)[names(reference) == desc$mediator_column] <- column
        ff_a_assert(is.function(ctx$engine$b3_design), "FF_A_FROZEN_BASIS_FUNCTION_MISSING")
        basis <- ctx$engine$b3_design(dd, column, base, reference)
        out$alpha <- alpha
        out$scenarios <- list(M0 = m0, M1 = m1)
        out$M0 <- m0
        out$M1 <- m1
        out$reference_basis <- basis$reference
        out$ledger$simulation_rule <- "FROZEN_EMPIRICAL_RESIDUAL_PRESERVING_ADDITIVE_SHIFT_NO_MC_DRAW"
    }
    else {
        levels <- base::levels(dd[[column]])
        ff_a_assert(length(levels) >= 2L, "FF_A_DEFINED_CATEGORY_LEVELS_MISSING")
        X <- stats::model.matrix(stats::reformulate(base), dd)
        if (!is.null(params$columns)) 
            ff_a_assert(identical(head(as.character(unlist(params$columns)), ncol(X)), colnames(X)), "FF_A_SAVED_CATEGORY_BASE_DESIGN_CHANGED")
        X0 <- X
        X0[, "B"] <- 0
        X1 <- X
        X1[, "B"] <- 1
        saved <- if (isTRUE(selection$A_refit_required)) 
            NULL
        else params$mediator
        if (!is.null(saved) && length(saved)) {
            ff_a_assert(identical(as.character(unlist(saved$levels)), levels) && identical(as.character(saved$reference), 
                levels[1L]), "FF_A_SAVED_CATEGORY_LEVEL_REFERENCE_DRIFT")
            out$g0 <- ff_a_category_probability(X0, saved, levels)
            out$g1 <- ff_a_category_probability(X1, saved, levels)
            out$ledger$A_origin <- "REUSED_NUMBER_SPECIFIC_CONDITIONED_ML_COEFFICIENTS"
        }
        else {
            ff_a_assert(isTRUE(ctx$allow_missing_a_recompute), "FF_A_MISSING_CATEGORY_A_RECOMPUTE_NOT_AUTHORIZED")
            ff_a_assert(is.function(ctx$discrete$c40_discrete_fit), "FF_A_FROZEN_DISCRETE_FITTER_MISSING")
            kind <- if (is.null(ca$kind)) 
                if (length(levels) == 2L) 
                  "binary"
                else "categorical"
            else ca$kind
            fit <- ctx$discrete$c40_discrete_fit(X, dd[[column]], mi$w, kind, ctx$engine)
            out$g0 <- fit$predict(X0)
            out$g1 <- fit$predict(X1)
            out$ledger$A_origin <- "RECOMPUTED_NUMBER_SPECIFIC_CATEGORY_A_OBJECT_MISSING"
            out$ledger$reason <- if (isTRUE(selection$A_refit_required)) 
                selection$origin
            else "OLD_CHECKPOINT_HAS_NO_CATEGORY_A_FOR_THIS_NUMBER_AND_SAMPLE"
        }
        out$levels <- levels
        out$reference_basis <- list(levels = levels, reference = levels[1L], contrast = "FROZEN_FACTOR_TREATMENT_DUMMIES_NO_ORDINAL_CONVERSION")
        out$ledger$simulation_rule <- "FROZEN_EXACT_ALL_CATEGORY_PROBABILITY_ENUMERATION_NO_MC_DRAW"
    }
    out
}

ff_a_bundle <- ff_a_scenarios

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ff_a_assert <- function(ok, code) if (!isTRUE(ok)) stop(code, call. = FALSE)
        ff_a_hash <- function(x) digest::digest(paste0(nchar(as.character(x)), ":", as.character(x), "\n", collapse = ""), 
            algo = "sha256", serialize = FALSE)
        ff_a_sha <- function(p) {
            ff_a_assert(file.exists(p) && !dir.exists(p), paste("FF_A_SOURCE_FILE_MISSING", p))
            digest::digest(file = p, algo = "sha256")
        }
        ff_a_check <- function(p, h) {
            ff_a_assert(is.character(h) && length(h) == 1L && grepl("^[a-f0-9]{64}$", h), "FF_A_EXACT_SOURCE_SHA256_REQUIRED")
            ff_a_assert(identical(ff_a_sha(p), h), paste("FF_A_SOURCE_SHA256_CHANGED", p))
        }
        ff_a_descriptor <- function(ctx, id, bid, endpoint) {
            if (is.function(ctx$a_source_for)) out <- ctx$a_source_for(id, bid, endpoint) else {
                key <- paste(id, bid, endpoint, sep = "|")
                index <- ctx$a_reuse_index
                if (is.environment(index)) out <- get0(key, index, inherits = FALSE) else out <- index[[key]]
            }
            ff_a_assert(is.list(out), paste("FF_A_SOURCE_NOT_AUDITED", id, bid, endpoint))
            ff_a_assert(!is.null(out$format) || (isTRUE(out$missing_object) && nzchar(out$reason)), "FF_A_SOURCE_DESCRIPTOR_INCOMPLETE")
            out
        }
        ff_a_read <- function(desc, endpoint, id = NULL, bid = NULL) {
            if (isTRUE(desc$missing_object)) return(NULL)
            parameter_endpoint <- if (is.null(desc$parameter_endpoint)) endpoint else desc$parameter_endpoint
            if (identical(desc$format, "c40_shard_gzip")) {
                ff_a_check(desc$path, desc$sha256)
                ff_a_assert(file.info(desc$path)$size <= 16 * 1024^2, "FF_A_C40_SHARD_COMPRESSED_LIMIT")
                con <- gzfile(desc$path, "rb")
                limit <- 64L * 1024L * 1024L
                bytes <- tryCatch(readBin(con, "raw", n = limit + 1L), finally = close(con))
                ff_a_assert(length(bytes) <= limit, "FF_A_C40_SHARD_DECOMPRESSED_LIMIT")
                shard <- jsonlite::fromJSON(rawToChar(bytes), simplifyVector = FALSE)
                ff_a_assert(identical(shard$schema, "C40_FIVE_ENDPOINT_AGGREGATE_SHARD_V2") && is.list(shard$bundles) && 
                  length(shard$bundles) >= 1L && length(shard$bundles) <= 25L, "FF_A_C40_SHARD_SCHEMA_CHANGED")
                units <- Filter(function(u) identical(u$indicator_id, id) && identical(as.integer(u$bootstrap_id), 
                  as.integer(bid)), shard$bundles)
                ff_a_assert(length(units) == 1L, "FF_A_C40_NUMBERED_UNIT_MISSING_OR_DUPLICATED")
                unit <- units[[1L]]
                ff_a_assert(identical(unit$draw$draw_hash, desc$draw_hash), "FF_A_C40_ACTUAL_DRAW_CHANGED")
                rows <- Filter(function(r) identical(as.character(r$endpoint_id), as.character(endpoint)), unit$rows)
                same_or_audited_missing <- function(value, expected) {
                  if (is.character(value) && length(value) == 1L && !is.na(value) && grepl("^[a-f0-9]{64}$", value)) identical(value, 
                    expected) else isTRUE(desc$binding_reconstruction_required)
                }
                ff_a_assert(length(rows) == 1L && same_or_audited_missing(rows[[1L]]$target_hash, desc$target_hash) && 
                  same_or_audited_missing(rows[[1L]]$sample_hash, desc$sample_hash), "FF_A_C40_TARGET_SAMPLE_PROOF_CHANGED")
                ff_a_assert(is.list(unit$files) && "PARAMETERS.json" %in% names(unit$files), "FF_A_C40_SEALED_PARAMETER_ARTIFACT_ABSENT")
                z <- unit$files[["PARAMETERS.json"]][[as.character(parameter_endpoint)]]
                ff_a_check(desc$path, desc$sha256)
                return(z)
            } else if (!is.null(desc$archive)) {
                ff_a_check(desc$archive, desc$archive_sha256)
                member <- desc$member
                ff_a_assert(is.character(member) && length(member) == 1L && nzchar(member) && !grepl("(^/|(^|/)[.][.](/|$))", 
                  member), "FF_A_UNSAFE_ARCHIVE_MEMBER")
                ff_a_assert(is.character(desc$member_sha256) && grepl("^[a-f0-9]{64}$", desc$member_sha256), "FF_A_MEMBER_SHA256_REQUIRED")
                limit <- 2L * 1024L * 1024L
                con <- pipe(paste("tar -xOf", shQuote(desc$archive), shQuote(member)), "rb")
                bytes <- tryCatch(readBin(con, "raw", n = limit + 1L), finally = close(con))
                ff_a_assert(length(bytes) <= limit, "FF_A_AGGREGATE_CHECKPOINT_EXCEEDS_2_MIB")
                ff_a_assert(identical(digest::digest(bytes, algo = "sha256", serialize = FALSE), desc$member_sha256), 
                  "FF_A_ARCHIVE_MEMBER_CHANGED")
                ff_a_assert(identical(desc$format, "c40_json"), "FF_A_ONLY_JSON_ARCHIVE_CHECKPOINT_SUPPORTED")
                z <- jsonlite::fromJSON(rawToChar(bytes), simplifyVector = TRUE)
                ff_a_check(desc$archive, desc$archive_sha256)
            } else {
                ff_a_check(desc$path, desc$sha256)
                z <- switch(desc$format, v5_rds = readRDS(desc$path), b18_rds = readRDS(desc$path), c40_json = jsonlite::fromJSON(desc$path, 
                  simplifyVector = TRUE), legacy_coefficients_csv = read.csv(desc$path, stringsAsFactors = FALSE, 
                  check.names = FALSE), stop("FF_A_UNKNOWN_CHECKPOINT_FORMAT", call. = FALSE))
                ff_a_check(desc$path, desc$sha256)
            }
            if (desc$format %in% c("v5_rds", "c40_json")) z <- z[[as.character(parameter_endpoint)]]
            z
        }
        ff_a_selection_formula <- function() stats::reformulate(c("B", "splines::ns(age, df=3)", "sex", "center", 
            "assessment_year", "S", "splines::ns(aux_t_250_2, knots=c(1,3), Boundary.knots=c(0,5))", "aux_event_250_2", 
            "aux_death_250_2", "splines::ns(aux_t_411, knots=c(1,3), Boundary.knots=c(0,5))", "aux_event_411", 
            "aux_death_411"), response = "complete")
        ff_a_selection <- function(ctx, d, column, base, params, desc) {
            complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
            ff_a_assert(any(complete), "NO_COMPLETE_RECORDS")
            coefficients <- params$selection_coefficients
            A_refit_required <- FALSE
            coefficient_max_abs_difference <- NA_real_
            if (all(complete)) {
                ff_a_assert(is.null(coefficients) || !length(coefficients), "FF_A_ALL_COMPLETE_CONFLICTS_WITH_SAVED_SELECTION_MODEL")
                mi <- list(complete = complete, w = rep(1, nrow(d)), fit = NULL, warnings = character(), diagnostics = data.frame())
                origin <- "FROZEN_ALL_COMPLETE_UNIT_SELECTION_WEIGHTS"
            } else if (!is.null(coefficients) && length(coefficients)) {
                coefficients <- unlist(coefficients, use.names = TRUE)
                dd <- d
                dd$complete <- complete
                f <- if (is.null(params$selection_formula)) ff_a_selection_formula() else stats::as.formula(paste(params$selection_formula, 
                  collapse = " "))
                X <- stats::model.matrix(f, stats::model.frame(f, dd, na.action = stats::na.fail, drop.unused.levels = FALSE))
                if (is.null(names(coefficients)) && length(coefficients) == ncol(X)) names(coefficients) <- colnames(X)
                ff_a_assert(identical(colnames(X), names(coefficients)) && all(is.finite(coefficients)), "FF_A_SAVED_SELECTION_DESIGN_OR_COEFFICIENT_DRIFT")
                pi <- stats::binomial()$linkinv(as.numeric(X %*% coefficients))
                ff_a_assert(all(is.finite(pi) & pi > 0 & pi <= 1), "FF_A_SELECTION_POSITIVITY_FAILURE")
                mi <- list(complete = complete, w = 1/pi[complete], fit = NULL, warnings = character(), diagnostics = data.frame(), 
                  saved_coefficients = coefficients)
                origin <- "RECONSTRUCTED_SEALED_FROZEN_SELECTION_COEFFICIENTS"
            } else {
                missing_confirmed <- isTRUE(desc$selection_object_missing) || isTRUE(desc$missing_object) || identical(desc$format, 
                  "legacy_coefficients_csv") || is.null(params)
                ff_a_assert(missing_confirmed, "FF_A_SELECTION_CHECKPOINT_MISSING_NOT_AUDITED")
                ff_a_assert(isTRUE(ctx$allow_missing_a_recompute), "FF_A_MISSING_SELECTION_RECOMPUTE_NOT_AUTHORIZED")
                ff_a_assert(is.function(ctx$engine$b3_missing), "FF_A_FROZEN_SELECTION_FUNCTION_MISSING")
                mi <- ctx$engine$b3_missing(d, column, base)
                ff_a_assert(identical(mi$complete, complete), "FF_A_RECOMPUTED_SELECTION_SAMPLE_CHANGED")
                origin <- "RECOMPUTED_FROZEN_SELECTION_ONLY_PERSISTED_OBJECT_MISSING"
                A_refit_required <- TRUE
                origin <- "A_NONREUSABLE_SELECTION_FIT_AND_WEIGHT_VECTOR_NOT_RETAINED"
            }
            ff_a_assert(length(mi$w) == sum(complete) && all(is.finite(mi$w) & mi$w > 0), "FF_A_INVALID_SELECTION_WEIGHT_VECTOR")
            hash <- ff_a_hash(format(mi$w, digits = 17, scientific = TRUE))
            if (!is.null(desc$selection_weight_hash) && !is.na(desc$selection_weight_hash) && nzchar(desc$selection_weight_hash) && 
                !identical(hash, desc$selection_weight_hash)) {
                reason <- "A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED"
                ff_a_assert(!isTRUE(desc$individual_weight_vector_retained), "FF_A_ACTUALLY_RETAINED_SELECTION_WEIGHT_VECTOR_MISMATCH_INPUT_DRIFT")
                ff_a_assert(isTRUE(ctx$allow_missing_a_recompute) && is.function(ctx$engine$b3_missing), paste("FF_A_SELECTION_WEIGHT_HASH_MISMATCH", 
                  reason))
                rebuilt <- ctx$engine$b3_missing(d, column, base)
                ff_a_assert(identical(rebuilt$complete, complete), paste("FF_A_REBUILT_SELECTION_SAMPLE_CHANGED", 
                  reason))
                ff_a_assert(!is.null(rebuilt$fit) && identical(attr(stats::terms(rebuilt$fit), "term.labels"), 
                  attr(stats::terms(ff_a_selection_formula()), "term.labels")), "FF_A_REBUILT_SELECTION_FROZEN_FORMULA_CHANGED")
                if (!is.null(coefficients) && length(coefficients)) {
                  newcoef <- stats::coef(rebuilt$fit)
                  ff_a_assert(identical(names(newcoef), names(coefficients)) && all(is.finite(newcoef)), "FF_A_REBUILT_SELECTION_COEFFICIENT_DESIGN_CHANGED")
                  coefficient_max_abs_difference <- max(abs(as.numeric(newcoef) - as.numeric(coefficients)))
                }
                rebuilt_hash <- ff_a_hash(format(rebuilt$w, digits = 17, scientific = TRUE))
                mi <- rebuilt
                hash <- rebuilt_hash
                A_refit_required <- TRUE
                origin <- reason
            }
            list(mi = mi, origin = origin, weight_hash = hash, old_weight_hash = if (is.null(desc$selection_weight_hash)) NA_character_ else desc$selection_weight_hash, 
                A_refit_required = A_refit_required, saved_selection_coefficient_max_abs_difference = coefficient_max_abs_difference)
        }
        ff_a_numeric <- function(x, expected_names = NULL) {
            ans <- unlist(x, use.names = TRUE)
            if (!is.null(expected_names)) {
                ff_a_assert(length(ans) == length(expected_names), "FF_A_PARAMETER_LENGTH_CHANGED")
                if (is.null(names(ans))) names(ans) <- expected_names
                ff_a_assert(identical(names(ans), expected_names), "FF_A_PARAMETER_COLUMNS_CHANGED")
            }
            ff_a_assert(is.numeric(ans) && all(is.finite(ans)), "FF_A_NONFINITE_SAVED_PARAMETER")
            ans
        }
        ff_a_category_probability <- function(X, saved, levels) {
            mu <- ff_a_numeric(saved$conditioning_mean, colnames(X))
            scale <- ff_a_numeric(saved$conditioning_scale, colnames(X))
            ff_a_assert(all(scale > 0), "FF_A_INVALID_SAVED_CONDITIONING_SCALE")
            Z <- sweep(sweep(X, 2L, mu, "-"), 2L, scale, "/")
            coef <- saved$coefficients
            if (identical(saved$estimator, "WEIGHTED_BINOMIAL_ML_GLM_FIT")) {
                ff_a_assert(length(levels) == 2L, "FF_A_BINARY_LEVEL_DRIFT")
                beta <- ff_a_numeric(coef, colnames(X))
                p <- stats::plogis(as.numeric(Z %*% beta))
                ans <- cbind(1 - p, p)
            } else {
                ff_a_assert(identical(saved$estimator, "WEIGHTED_NOMINAL_MULTINOMIAL_ML_NNET_DECAY_ZERO"), "FF_A_UNKNOWN_SAVED_MEDIATOR_ESTIMATOR")
                beta <- if (is.matrix(coef)) coef else if (is.data.frame(coef)) as.matrix(coef) else if (is.list(coef) && 
                  length(coef) == length(levels) - 1L) do.call(rbind, lapply(coef, unlist)) else as.numeric(unlist(coef))
                if (length(levels) == 2L && is.null(dim(beta))) beta <- matrix(beta, nrow = 1L)
                ff_a_assert(is.matrix(beta) && identical(dim(beta), c(length(levels) - 1L, ncol(X))) && all(is.finite(beta)), 
                  "FF_A_MULTINOMIAL_COEFFICIENT_DIMENSION_DRIFT")
                eta <- cbind(0, Z %*% t(beta))
                eta <- eta - apply(eta, 1L, max)
                p <- exp(eta)
                ans <- p/rowSums(p)
            }
            colnames(ans) <- levels
            ff_a_assert(all(is.finite(ans) & ans >= 0 & ans <= 1) && max(abs(rowSums(ans) - 1)) < 1e-10, "FF_A_INVALID_CATEGORY_PROBABILITIES_NO_CLIPPING")
            ans
        }
        ff_a_scenarios <- function(ctx, id, bid, ep, pr, d, base) {
            ff_a_assert(is.numeric(bid) && length(bid) == 1L && bid == as.integer(bid) && bid %in% 0:1000, "FF_A_FIXED_ID_ZERO_TO_1000_REQUIRED")
            ff_a_assert(is.data.frame(d) && length(pr) == nrow(d) && length(pr) > 0L && identical(as.character(d$eid), 
                as.character(ctx$d_full$eid[pr])), "FF_A_ACTUAL_PARENT_DRAW_ROW_BINDING_MISMATCH")
            ca <- ctx$features[[id]]
            ff_a_assert(is.list(ca) && ca$column %in% names(d), "FF_A_FROZEN_FEATURE_MISSING")
            column <- ca$column
            desc <- ff_a_descriptor(ctx, id, as.integer(bid), ep)
            target_hash <- ff_a_hash(pr)
            binding_reconstructed <- FALSE
            if (isTRUE(desc$binding_reconstruction_required)) {
                ff_a_assert(identical(desc$source_input_lock_sha256, ctx$old_input_lock_sha256) && identical(desc$source_prepared_sha256, 
                  ctx$old_prepared_inputs_sha256), "FF_A_MISSING_BINDING_NOT_BOUND_TO_FROZEN_OLD_INPUTS")
                missing_hash <- function(h) is.null(h) || length(h) != 1L || is.na(h) || !is.character(h) || !grepl("^[a-f0-9]{64}$", 
                  h)
                if (missing_hash(desc$target_hash)) desc$target_hash <- target_hash
                if (missing_hash(desc$sample_hash)) desc$sample_hash <- ff_a_hash(pr[stats::complete.cases(d[, 
                  unique(c(base, column)), drop = FALSE])])
                binding_reconstructed <- TRUE
            }
            ff_a_assert(is.character(desc$target_hash) && identical(desc$target_hash, target_hash), "FF_A_ACTUAL_TARGET_HASH_CHANGED")
            complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
            sample_hash <- ff_a_hash(pr[complete])
            ff_a_assert(is.character(desc$sample_hash) && identical(desc$sample_hash, sample_hash), "FF_A_ACTUAL_COMPLETE_SAMPLE_HASH_CHANGED")
            ff_a_assert(is.character(desc$draw_hash) && length(desc$draw_hash) == 1L && grepl("^[a-f0-9]{64}$", 
                desc$draw_hash), "FF_A_ACTUAL_FULL_PARENT_DRAW_PROOF_REQUIRED")
            ff_a_assert(is.character(ctx$actual_draw_hash) && length(ctx$actual_draw_hash) == 1L && identical(desc$draw_hash, 
                ctx$actual_draw_hash), "FF_A_ACTUAL_DRAW_HASH_CHANGED_OR_NOT_BOUND")
            params <- ff_a_read(desc, ep, id, bid)
            if (!is.null(params$base_terms)) ff_a_assert(identical(as.character(unlist(params$base_terms)), as.character(base)), 
                "FF_A_SAVED_BASE_FORMULA_CHANGED")
            if (!is.null(params$factor_levels)) for (name in names(params$factor_levels)) ff_a_assert(name %in% 
                names(d) && is.factor(d[[name]]) && identical(levels(d[[name]]), as.character(unlist(params$factor_levels[[name]]))), 
                "FF_A_SAVED_FACTOR_LEVEL_DICTIONARY_CHANGED")
            selection <- ff_a_selection(ctx, d, column, base, params, desc)
            mi <- selection$mi
            dd <- d[mi$complete, , drop = FALSE]
            factor_mediator <- is.factor(dd[[column]])
            ledger <- data.frame(indicator_id = id, endpoint_id = ep, bootstrap_id = as.integer(bid), source_path = if (!is.null(desc$archive)) desc$archive else if (is.null(desc$path)) NA_character_ else desc$path, 
                source_member = if (is.null(desc$member)) NA_character_ else desc$member, source_sha256 = if (!is.null(desc$archive)) desc$archive_sha256 else if (is.null(desc$sha256)) NA_character_ else desc$sha256, 
                source_format = if (is.null(desc$format)) "AUDITED_MISSING" else desc$format, A_source_endpoint = if (is.null(desc$parameter_endpoint)) ep else desc$parameter_endpoint, 
                source_missing_reason = if (is.null(desc$reason)) "" else desc$reason, draw_hash = desc$draw_hash, 
                target_hash = target_hash, sample_hash = sample_hash, draw_proof_kind = if (is.null(desc$draw_proof_kind)) "AUDITED_ACTUAL_DRAW_BINDING" else desc$draw_proof_kind, 
                draw_artifact = if (is.null(desc$draw_artifact)) NA_character_ else desc$draw_artifact, old_failed_cell_missing_sample_binding_reconstructed = binding_reconstructed, 
                selection_weight_hash = selection$weight_hash, selection_origin = selection$origin, old_selection_weight_hash = selection$old_weight_hash, 
                selection_and_A_rebuilt_not_claimed_bitwise_old_weights = selection$A_refit_required, saved_selection_coefficient_max_abs_difference = selection$saved_selection_coefficient_max_abs_difference, 
                N_target = nrow(d), N_complete = nrow(dd), A_origin = NA_character_, simulation_rule = NA_character_, 
                reason = "", stringsAsFactors = FALSE)
            out <- list(mi = mi, data = dd, column = column, parent_rows = pr[mi$complete], ledger = ledger, reference_basis = NULL, 
                alpha = NULL, scenarios = NULL, g0 = NULL, g1 = NULL, levels = NULL)
            if (!factor_mediator) {
                alpha <- if (isTRUE(selection$A_refit_required)) NULL else if (is.data.frame(params) && identical(desc$format, 
                  "legacy_coefficients_csv")) {
                  term <- paste0(if (is.null(desc$mediator_column)) column else desc$mediator_column, "__a_path")
                  q <- params[params$term == term, , drop = FALSE]
                  if (nrow(q) == 1L) stats::setNames(q$coefficient, column) else NULL
                } else params$alpha
                if (!is.null(alpha) && length(alpha)) {
                  alpha <- unlist(alpha, use.names = TRUE)
                  old_column <- if (is.null(desc$mediator_column)) column else desc$mediator_column
                  ff_a_assert(old_column %in% names(alpha), "FF_A_SAVED_ALPHA_COLUMN_MISSING")
                  alpha <- stats::setNames(as.numeric(alpha[[old_column]]), column)
                  ff_a_assert(is.finite(alpha), "FF_A_SAVED_ALPHA_NONFINITE")
                  out$ledger$A_origin <- "REUSED_NUMBER_SPECIFIC_SAVED_ALPHA"
                } else {
                  ff_a_assert(isTRUE(ctx$allow_missing_a_recompute), "FF_A_MISSING_ALPHA_RECOMPUTE_NOT_AUTHORIZED")
                  ff_a_assert(is.function(ctx$engine$fp_a_path), "FF_A_FROZEN_A_FUNCTION_MISSING")
                  aa <- ctx$engine$fp_a_path(d, column, base, mi, ctx, paste0("single_", ca$candidate_id), as.integer(bid))
                  alpha <- aa$alpha
                  out$ledger$A_origin <- "RECOMPUTED_NUMBER_SPECIFIC_A_OBJECT_MISSING"
                  out$ledger$reason <- if (isTRUE(selection$A_refit_required)) selection$origin else "OLD_CHECKPOINT_HAS_NO_ALPHA_FOR_THIS_NUMBER_AND_SAMPLE"
                }
                m0 <- dd[[column]] - dd$B * alpha[[column]]
                m1 <- m0 + alpha[[column]]
                lo <- if (is.null(ca$physical_lower_bound)) -Inf else ca$physical_lower_bound
                hi <- if (is.null(ca$physical_upper_bound)) Inf else ca$physical_upper_bound
                ff_a_assert(all(is.finite(m0)) && all(is.finite(m1)) && all(m0 >= lo & m0 <= hi) && all(m1 >= lo & 
                  m1 <= hi), "SCENARIO_OUTSIDE_PHYSICAL_DOMAIN_NO_CLIPPING")
                reference <- params$spline_reference
                if (!is.null(reference) && !column %in% names(reference) && !is.null(desc$mediator_column) && desc$mediator_column %in% 
                  names(reference)) names(reference)[names(reference) == desc$mediator_column] <- column
                ff_a_assert(is.function(ctx$engine$b3_design), "FF_A_FROZEN_BASIS_FUNCTION_MISSING")
                basis <- ctx$engine$b3_design(dd, column, base, reference)
                out$alpha <- alpha
                out$scenarios <- list(M0 = m0, M1 = m1)
                out$M0 <- m0
                out$M1 <- m1
                out$reference_basis <- basis$reference
                out$ledger$simulation_rule <- "FROZEN_EMPIRICAL_RESIDUAL_PRESERVING_ADDITIVE_SHIFT_NO_MC_DRAW"
            } else {
                levels <- base::levels(dd[[column]])
                ff_a_assert(length(levels) >= 2L, "FF_A_DEFINED_CATEGORY_LEVELS_MISSING")
                X <- stats::model.matrix(stats::reformulate(base), dd)
                if (!is.null(params$columns)) ff_a_assert(identical(head(as.character(unlist(params$columns)), 
                  ncol(X)), colnames(X)), "FF_A_SAVED_CATEGORY_BASE_DESIGN_CHANGED")
                X0 <- X
                X0[, "B"] <- 0
                X1 <- X
                X1[, "B"] <- 1
                saved <- if (isTRUE(selection$A_refit_required)) NULL else params$mediator
                if (!is.null(saved) && length(saved)) {
                  ff_a_assert(identical(as.character(unlist(saved$levels)), levels) && identical(as.character(saved$reference), 
                    levels[1L]), "FF_A_SAVED_CATEGORY_LEVEL_REFERENCE_DRIFT")
                  out$g0 <- ff_a_category_probability(X0, saved, levels)
                  out$g1 <- ff_a_category_probability(X1, saved, levels)
                  out$ledger$A_origin <- "REUSED_NUMBER_SPECIFIC_CONDITIONED_ML_COEFFICIENTS"
                } else {
                  ff_a_assert(isTRUE(ctx$allow_missing_a_recompute), "FF_A_MISSING_CATEGORY_A_RECOMPUTE_NOT_AUTHORIZED")
                  ff_a_assert(is.function(ctx$discrete$c40_discrete_fit), "FF_A_FROZEN_DISCRETE_FITTER_MISSING")
                  kind <- if (is.null(ca$kind)) if (length(levels) == 2L) "binary" else "categorical" else ca$kind
                  fit <- ctx$discrete$c40_discrete_fit(X, dd[[column]], mi$w, kind, ctx$engine)
                  out$g0 <- fit$predict(X0)
                  out$g1 <- fit$predict(X1)
                  out$ledger$A_origin <- "RECOMPUTED_NUMBER_SPECIFIC_CATEGORY_A_OBJECT_MISSING"
                  out$ledger$reason <- if (isTRUE(selection$A_refit_required)) selection$origin else "OLD_CHECKPOINT_HAS_NO_CATEGORY_A_FOR_THIS_NUMBER_AND_SAMPLE"
                }
                out$levels <- levels
                out$reference_basis <- list(levels = levels, reference = levels[1L], contrast = "FROZEN_FACTOR_TREATMENT_DUMMIES_NO_ORDINAL_CONVERSION")
                out$ledger$simulation_rule <- "FROZEN_EXACT_ALL_CATEGORY_PROBABILITY_ENUMERATION_NO_MC_DRAW"
            }
            out
        }
        ff_a_bundle <- ff_a_scenarios
    }), envir = context)
}

