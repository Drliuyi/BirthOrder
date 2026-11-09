# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

new121_y_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

new121_y_engine <- function(old_root, source_parent = dirname(old_root), numerical_backend = NULL, build_dir = NULL) {
    for (pkg in c("survival", "digest", "jsonlite", "splines")) new121_y_assert(requireNamespace(pkg, quietly = TRUE), 
        paste("INSTALLED_PACKAGE_REQUIRED", pkg))
    old_root <- normalizePath(old_root, mustWork = TRUE)
    source_parent <- normalizePath(source_parent, mustWork = TRUE)
    message("Y source root: ", old_root, "; shared source parent: ", source_parent)
    lock_path <- file.path(old_root, "config", "SOURCE_LOCK.json")
    new121_y_assert(file.exists(lock_path), paste("MISSING_SOURCE_LOCK", lock_path))
    lock <- jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
    wanted <- c("f/30.outcome.R", "f/61.center_strata.R", "R/fast_predict.R")
    engine <- new.env(parent = globalenv())
    source_audit <- list()
    for (suffix in wanted) {
        entries <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)
        new121_y_assert(length(entries) == 1L, paste("PINNED_Y_SOURCE_NOT_UNIQUE", suffix))
        entry <- entries[[1L]]
        path <- file.path(source_parent, entry$path)
        new121_y_assert(file.exists(path), paste("MISSING_Y_SOURCE", path))
        sha <- digest::digest(file = path, algo = "sha256")
        new121_y_assert(identical(sha, entry$sha256), paste("PINNED_Y_SOURCE_CHANGED", path))
        sys.source(path, envir = engine)
        source_audit[[suffix]] <- list(path = path, sha256 = sha)
    }
    overlay <- file.path(old_root, "repairs", "exponential_fast_b100_v1")
    helpers <- c(file.path(old_root, "R", "core39_models.R"), file.path(overlay, "y_model", "ethnicity_pool.R"), 
        file.path(overlay, "numerics", "bounded_aj.R"))
    for (path in helpers) {
        new121_y_assert(file.exists(path), paste("MISSING_APPROVED_HELPER", path))
        sys.source(path, envir = engine)
        source_audit[[basename(path)]] <- list(path = path, sha256 = digest::digest(file = path, algo = "sha256"))
    }
    engine$c39_uy_install(engine)
    if (is.null(numerical_backend)) {
        new121_y_assert(is.character(build_dir) && length(build_dir) == 1L && nzchar(build_dir), "EXPLICIT_NEW_NUMERICAL_BUILD_DIR_OR_ALREADY_LOADED_BACKEND_REQUIRED")
        message("Y worker numerical build root: ", build_dir)
        numerical_backend <- engine$c39_num_initialize(file.path(overlay, "numerics"), build_dir)
    }
    new121_y_assert(is.list(numerical_backend) && !is.null(numerical_backend$symbol) && identical(numerical_backend$manifest$spec$computation_precision, 
        "DOUBLE_ONLY_BINARY64"), "APPROVED_ALREADY_LOADED_DOUBLE_EXPONENTIAL_BACKEND_REQUIRED")
    engine$c39_num_install(engine, numerical_backend)
    engine$new121_source_audit <- source_audit
    engine$new121_numeric_backend <- numerical_backend
    engine
}

new121_y_exact_alias_plan <- function(X, mediator_design_columns) {
    new121_y_assert(is.matrix(X) && is.numeric(X) && all(is.finite(X)) && !is.null(colnames(X)) && !anyDuplicated(colnames(X)), 
        "FINITE_NAMED_FULL_X_REQUIRED")
    if ("(Intercept)" %in% colnames(X)) {
        new121_y_assert(all(X[, "(Intercept)"] == 1), "INVALID_INTERCEPT")
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    new121_y_assert(all(mediator_design_columns %in% colnames(X)) && !anyDuplicated(mediator_design_columns) && 
        !any(c("B", "age") %in% mediator_design_columns), "EXPLICIT_MEDIATION_ONLY_COLUMN_OWNERSHIP_REQUIRED")
    base_columns <- setdiff(colnames(X), mediator_design_columns)
    augmented <- cbind(`(Intercept)` = rep(1, nrow(X)), X)
    base <- c("(Intercept)", base_columns)
    new121_y_assert(qr(augmented[, base, drop = FALSE], tol = 1e-11)$rank == length(base), "BASE_COVARIATES_UNIDENTIFIABLE_NO_SILENT_BASE_CATEGORY_REMOVAL")
    ordered <- c(base, mediator_design_columns)
    q <- qr(augmented[, ordered, drop = FALSE], tol = 1e-11)
    keep <- ordered[q$pivot[seq_len(q$rank)]]
    new121_y_assert(all(base %in% keep), "QR_WOULD_REMOVE_PROTECTED_COVARIATE")
    dropped <- setdiff(ordered, keep)
    witness <- list()
    if (length(dropped)) {
        kept_X <- augmented[, keep, drop = FALSE]
        coefficients <- qr.coef(qr(kept_X, tol = 1e-13), augmented[, dropped, drop = FALSE])
        coefficients <- matrix(coefficients, nrow = length(keep), dimnames = list(keep, dropped))
        for (column in dropped) {
            coefficient <- coefficients[, column]
            coefficient[setdiff(base, "(Intercept)")] <- 0
            fitted <- as.numeric(kept_X %*% coefficient)
            bound <- 128 * .Machine$double.eps * pmax(1, abs(augmented[, column]) + as.numeric(abs(kept_X) %*% 
                abs(coefficient)))
            residual <- abs(augmented[, column] - fitted)
            new121_y_assert(all(is.finite(coefficient)) && all(residual <= bound), paste("NOT_AN_EXACT_MEDIATION_ALGEBRAIC_ALIAS_NO_HIGH_CORRELATION_DROP", 
                column))
            witness[[column]] <- list(coefficients = coefficient, fit_max_abs_residual = max(residual), roundoff_bound_max = max(bound), 
                rule = "MEDIATION_PLUS_INTERCEPT_IDENTITY_AT_ROUNDOFF_CHECKED_EVERY_PREDICTION")
        }
    }
    list(input_columns = colnames(X), keep_columns = setdiff(keep, "(Intercept)"), dropped_columns = dropped, aliases = witness, 
        base_columns = base_columns, approximate_dimension_reduction = FALSE)
}

new121_y_apply_alias_plan <- function(X, plan) {
    if ("(Intercept)" %in% colnames(X)) {
        new121_y_assert(all(X[, "(Intercept)"] == 1), "INVALID_PREDICTION_INTERCEPT")
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    new121_y_assert(is.matrix(X) && is.numeric(X) && all(is.finite(X)) && identical(colnames(X), plan$input_columns), 
        "FULL_PREDICTION_X_COLUMN_DRIFT")
    augmented <- cbind(`(Intercept)` = rep(1, nrow(X)), X)
    for (column in names(plan$aliases)) {
        coefficient <- plan$aliases[[column]]$coefficients
        kept_X <- augmented[, names(coefficient), drop = FALSE]
        residual <- abs(augmented[, column] - as.numeric(kept_X %*% coefficient))
        bound <- 128 * .Machine$double.eps * pmax(1, abs(augmented[, column]) + as.numeric(abs(kept_X) %*% abs(coefficient)))
        new121_y_assert(all(residual <= bound), paste("ALGEBRAIC_REPRESENTATION_IDENTITY_FAILED_ON_GENERATED_DATA", 
            column))
    }
    X[, plan$keep_columns, drop = FALSE]
}

new121_fit_y <- function(engine, fullX, risks, weights, cluster, mediator_design_columns, centre_actual = NULL, 
    audit_callback = NULL, diagnostics = FALSE) {
    new121_y_assert(is.environment(engine) && is.function(engine$ff_ce_fit) && isTRUE(engine$ff_aj_status()$ready) && 
        identical(engine$ff_aj_status()$computation_precision, "DOUBLE_ONLY_BINARY64"), "INSTALLED_DOUBLE_EXPONENTIAL_Y_ENGINE_REQUIRED")
    new121_y_assert(is.list(risks) && setequal(names(risks), c("IHD", "T2D")), "PAIRED_IHD_T2D_RISK_FRAMES_REQUIRED")
    new121_y_assert(length(cluster) == nrow(fullX) && !anyNA(cluster), "EXPLICIT_ROBUST_FAMILY_CLUSTER_REQUIRED")
    plan <- new121_y_exact_alias_plan(fullX, mediator_design_columns)
    X <- new121_y_apply_alias_plan(fullX, plan)
    fits <- list()
    for (endpoint in c("IHD", "T2D")) {
        if (is.function(audit_callback)) 
            audit_callback(list(stage = "Y_FITTING", endpoint = endpoint))
        fits[[endpoint]] <- engine$ff_ce_fit(X = X, risk = risks[[endpoint]], w = weights, cluster = cluster, tau_report = 11, 
            ties = "efron", diagnostics = diagnostics, endpoint_id = if (endpoint == "IHD") 
                "411"
            else "250.2", centre_actual = centre_actual, audit_callback = if (is.function(audit_callback)) 
                local({
                  current_endpoint <- endpoint
                  function(record) audit_callback(c(list(endpoint = current_endpoint), record))
                })
            else NULL)
    }
    compiled <- lapply(fits, engine$fd_fp_compile, v4 = engine)
    list(fits = fits, compiled = compiled, alias_plan = plan, method = list(Y = engine$ff_ce_method_spec(), numeric = engine$c39_num_spec(), 
        sources = engine$new121_source_audit, fit_scope = "FRESH_CURRENT_PANEL_BOOTSTRAP_NO_OLD_COEFFICIENT_REUSE", 
        tau = 11, fit_followup_cap = Inf, ties = "efron", robust_family_cluster = TRUE, warnings_policy = "INHERITED_CLASSIFIED_WARNINGS_ONE_ELIGIBLE_NATIVE_NUMERICAL_RETRY_NO_PENALTY", 
        dropped_exact_aliases = plan$dropped_columns))
}

new121_y_design <- function(data, members, base, reference = NULL, engine) {
    engine$c39_design(data, members, base, reference)
}

new121_y_design_factory <- function(data, members, base, domains, engine, cache_entries = 4L, reference = NULL) {
    new121_y_assert(cache_entries >= 1L && cache_entries == as.integer(cache_entries), "POSITIVE_CACHE_BOUND_REQUIRED")
    columns <- vapply(members, `[[`, "", "column")
    new121_y_assert(setequal(unlist(domains, use.names = FALSE), columns) && !anyDuplicated(unlist(domains, use.names = FALSE)) && 
        !length(intersect(base, columns)), "EXACT_EXPANDED_COORDINATE_DOMAIN_PARTITION_REQUIRED")
    initial <- engine$c39_design(data, members, base, reference)
    reference <- initial$reference
    continuous <- Filter(function(m) m$type == "continuous", members)
    discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
    form <- stats::reformulate(unique(c(base, discrete)))
    mm <- stats::model.matrix(form, data)
    labels <- attr(stats::terms(form), "term.labels")
    assignment <- attr(mm, "assign")
    owner <- setNames(rep("BASE", ncol(initial$X)), colnames(initial$X))
    for (g in names(domains)) for (column in domains[[g]]) {
        member <- members[[match(column, columns)]]
        owned <- if (member$type == "continuous") 
            reference[[member$indicator_id]]$column_names
        else colnames(mm)[assignment == match(column, labels)]
        new121_y_assert(all(owned %in% names(owner)), "UNMAPPED_EXPANDED_MEDIATOR_BASIS")
        owner[owned] <- g
    }
    state <- new.env(parent = emptyenv())
    state$cache <- list()
    state$builds <- 0L
    state$hits <- 0L
    basis <- function(d, key, vars, build) {
        values <- d[, vars, drop = FALSE]
        entries <- state$cache[[key]]
        hit <- which(vapply(entries, function(e) identical(e$values, values), logical(1)))
        if (length(hit)) {
            state$hits <- state$hits + 1L
            return(entries[[hit[1L]]]$X)
        }
        X <- build()
        state$cache[[key]] <- head(c(list(list(values = values, X = X)), entries), cache_entries)
        state$builds <- state$builds + 1L
        X
    }
    design <- function(b, M, C) {
        new121_y_assert(length(b) == 1L && b %in% 0:1 && nrow(M) == nrow(C), "SCALAR_EXPOSURE_AND_MATCHED_ROWS_REQUIRED")
        d <- C
        d$B <- b
        for (column in columns) d[[column]] <- M[[column]]
        for (column in intersect(names(d), names(data))) if (is.factor(data[[column]])) {
            ref <- data[[column]]
            new121_y_assert(all(as.character(d[[column]]) %in% levels(ref)), paste("UNKNOWN_FROZEN_FACTOR_LEVEL", 
                column))
            d[[column]] <- factor(as.character(d[[column]]), levels = levels(ref), ordered = is.ordered(ref))
            if (!is.null(attr(ref, "contrasts"))) 
                attr(d[[column]], "contrasts") <- attr(ref, "contrasts")
        }
        X <- basis(d, "BASE", base, function() stats::model.matrix(stats::reformulate(base), d))
        for (column in discrete) {
            part <- basis(d, column, column, function() {
                q <- stats::model.matrix(stats::reformulate(column), d)
                q[, colnames(q) != "(Intercept)", drop = FALSE]
            })
            X <- cbind(X, part)
        }
        for (member in continuous) {
            column <- member$column
            ref <- reference[[member$indicator_id]]
            X <- cbind(X, basis(d, column, column, function() engine$c39_basis(d[[column]], ref, column)))
        }
        new121_y_assert(setequal(colnames(X), initial$fullcolnames) && all(is.finite(X)), "CACHED_FULL_X_DRIFT")
        X[, initial$fullcolnames, drop = FALSE]
    }
    attr(design, "column_owner") <- owner
    attr(design, "reference") <- reference
    attr(design, "fullX") <- initial$X
    attr(design, "stats") <- function() list(basis_builds = state$builds, basis_cache_hits = state$hits)
    design
}

new121_y_prediction <- function(fit_bundle, design_fn, engine, domains = NULL, cache_entries = 4L) {
    new121_y_assert(is.function(design_fn) && is.list(fit_bundle$compiled) && setequal(names(fit_bundle$compiled), 
        c("IHD", "T2D")), "FRESH_PAIRED_FITS_AND_DESIGN_CALLBACK_REQUIRED")
    new121_y_assert(cache_entries >= 1L && cache_entries == as.integer(cache_entries), "POSITIVE_CACHE_BOUND_REQUIRED")
    compiled <- fit_bundle$compiled
    plan <- fit_bundle$alias_plan
    owner <- attr(design_fn, "column_owner")
    if (is.null(owner)) 
        owner <- setNames(rep("ALL", length(plan$keep_columns)), plan$keep_columns)
    owner <- c(owner, ethnicityAsian_or_Chinese_pooled = "BASE")
    state <- new.env(parent = emptyenv())
    state$cache <- list()
    state$predictions <- state$builds <- state$hits <- 0L
    lp <- function(X, beta, key) {
        owners <- owner[colnames(X)]
        owners[is.na(owners)] <- "BASE"
        total <- numeric(nrow(X))
        for (g in unique(unname(owners))) {
            ii <- which(owners == g)
            part <- X[, ii, drop = FALSE]
            name <- paste(key, g, sep = "|")
            entries <- state$cache[[name]]
            hit <- which(vapply(entries, function(e) identical(e$X, part), logical(1)))
            if (length(hit)) {
                z <- entries[[hit[1L]]]$lp
                state$hits <- state$hits + 1L
            }
            else {
                z <- as.numeric(part %*% beta[ii])
                state$cache[[name]] <- head(c(list(list(X = part, lp = z)), entries), cache_entries)
                state$builds <- state$builds + 1L
            }
            total <- total + z
        }
        total
    }
    predict_all <- function(b, M, C, validate = FALSE, weights = NULL) {
        X <- new121_y_apply_alias_plan(design_fn(b = b, M = M, C = C), plan)
        probabilities <- list()
        diffs <- numeric()
        for (endpoint in c("IHD", "T2D")) {
            cfit <- compiled[[endpoint]]
            prepared <- engine$fd_fp_prepare(X, cfit, engine)
            lps <- vapply(1:2, function(j) lp(prepared$designs[[j]], cfit$beta[[j]], paste(endpoint, j, sep = "|")), 
                numeric(nrow(X)))
            rr <- exp(matrix(lps, nrow = nrow(X), ncol = 2L))
            new121_y_assert(all(is.finite(rr) & rr > 0), "NONFINITE_RISK_MULTIPLIER_NO_CLIPPING")
            p <- matrix(NA_real_, nrow(X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
            for (g in unique(prepared$group)) {
                ii <- which(prepared$group == g)
                p[ii, ] <- engine$ff_aj_backend(cfit$hazards[[g]], rr[ii, , drop = FALSE])
            }
            engine$ff_y_probability_check(p)
            probabilities[[endpoint]] <- p
            if (validate) {
                native <- engine$ff_ce_predict_y_matrix(cfit$object, X, 11)
                new121_y_assert(!engine$ff_ce_bad_warning(attr(native, "warnings")), "GENUINE_NATIVE_PREDICTION_WARNING")
                diffs <- c(diffs, abs(p - native))
            }
        }
        state$predictions <- state$predictions + 1L
        risk <- cbind(IHD = probabilities$IHD[, "disease"], T2D = probabilities$T2D[, "disease"])
        validation <- list(performed = validate, engineering_probability_reference = 1e-06, maximum_absolute_probability_difference = if (length(diffs)) max(diffs) else NULL, 
            mean_absolute_probability_difference = if (length(diffs)) mean(diffs) else NULL, bitwise_agreement_required = FALSE)
        if (validate && !is.null(weights)) {
            new121_y_assert(is.numeric(weights) && length(weights) == nrow(X) && all(is.finite(weights) & weights >= 
                0) && sum(weights) > 0, "INVALID_BLOCK_COMPARISON_WEIGHTS")
            validation$weighted_mean_absolute_probability_difference <- sum(diffs * rep(weights, 6L))/(sum(weights) * 
                6)
        }
        if (validate) 
            validation$within_engineering_reference <- max(diffs) <= 1e-06
        list(risk = risk, all_states = probabilities, validation = validation)
    }
    list(predict = function(b, M, C) predict_all(b, M, C)$risk, validate = function(b, M, C, weights = NULL) predict_all(b, 
        M, C, TRUE, weights), stats = function() c(list(predictions = state$predictions, lp_builds = state$builds, 
        lp_cache_hits = state$hits, bounded_entries_per_piece = cache_entries), if (is.function(attr(design_fn, 
        "stats"))) attr(design_fn, "stats")() else list()), design = design_fn, alias_plan = plan)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        new121_y_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        new121_y_engine <- function(old_root, source_parent = dirname(old_root), numerical_backend = NULL, build_dir = NULL) {
            for (pkg in c("survival", "digest", "jsonlite", "splines")) new121_y_assert(requireNamespace(pkg, quietly = TRUE), 
                paste("INSTALLED_PACKAGE_REQUIRED", pkg))
            old_root <- normalizePath(old_root, mustWork = TRUE)
            source_parent <- normalizePath(source_parent, mustWork = TRUE)
            message("Y source root: ", old_root, "; shared source parent: ", source_parent)
            lock_path <- file.path(old_root, "config", "SOURCE_LOCK.json")
            new121_y_assert(file.exists(lock_path), paste("MISSING_SOURCE_LOCK", lock_path))
            lock <- jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
            wanted <- c("f/30.outcome.R", "f/61.center_strata.R", "R/fast_predict.R")
            engine <- new.env(parent = globalenv())
            source_audit <- list()
            for (suffix in wanted) {
                entries <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)
                new121_y_assert(length(entries) == 1L, paste("PINNED_Y_SOURCE_NOT_UNIQUE", suffix))
                entry <- entries[[1L]]
                path <- file.path(source_parent, entry$path)
                new121_y_assert(file.exists(path), paste("MISSING_Y_SOURCE", path))
                sha <- digest::digest(file = path, algo = "sha256")
                new121_y_assert(identical(sha, entry$sha256), paste("PINNED_Y_SOURCE_CHANGED", path))
                sys.source(path, envir = engine)
                source_audit[[suffix]] <- list(path = path, sha256 = sha)
            }
            overlay <- file.path(old_root, "repairs", "exponential_fast_b100_v1")
            helpers <- c(file.path(old_root, "R", "core39_models.R"), file.path(overlay, "y_model", "ethnicity_pool.R"), 
                file.path(overlay, "numerics", "bounded_aj.R"))
            for (path in helpers) {
                new121_y_assert(file.exists(path), paste("MISSING_APPROVED_HELPER", path))
                sys.source(path, envir = engine)
                source_audit[[basename(path)]] <- list(path = path, sha256 = digest::digest(file = path, algo = "sha256"))
            }
            engine$c39_uy_install(engine)
            if (is.null(numerical_backend)) {
                new121_y_assert(is.character(build_dir) && length(build_dir) == 1L && nzchar(build_dir), "EXPLICIT_NEW_NUMERICAL_BUILD_DIR_OR_ALREADY_LOADED_BACKEND_REQUIRED")
                message("Y worker numerical build root: ", build_dir)
                numerical_backend <- engine$c39_num_initialize(file.path(overlay, "numerics"), build_dir)
            }
            new121_y_assert(is.list(numerical_backend) && !is.null(numerical_backend$symbol) && identical(numerical_backend$manifest$spec$computation_precision, 
                "DOUBLE_ONLY_BINARY64"), "APPROVED_ALREADY_LOADED_DOUBLE_EXPONENTIAL_BACKEND_REQUIRED")
            engine$c39_num_install(engine, numerical_backend)
            engine$new121_source_audit <- source_audit
            engine$new121_numeric_backend <- numerical_backend
            engine
        }
        new121_y_exact_alias_plan <- function(X, mediator_design_columns) {
            new121_y_assert(is.matrix(X) && is.numeric(X) && all(is.finite(X)) && !is.null(colnames(X)) && !anyDuplicated(colnames(X)), 
                "FINITE_NAMED_FULL_X_REQUIRED")
            if ("(Intercept)" %in% colnames(X)) {
                new121_y_assert(all(X[, "(Intercept)"] == 1), "INVALID_INTERCEPT")
                X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            }
            new121_y_assert(all(mediator_design_columns %in% colnames(X)) && !anyDuplicated(mediator_design_columns) && 
                !any(c("B", "age") %in% mediator_design_columns), "EXPLICIT_MEDIATION_ONLY_COLUMN_OWNERSHIP_REQUIRED")
            base_columns <- setdiff(colnames(X), mediator_design_columns)
            augmented <- cbind(`(Intercept)` = rep(1, nrow(X)), X)
            base <- c("(Intercept)", base_columns)
            new121_y_assert(qr(augmented[, base, drop = FALSE], tol = 1e-11)$rank == length(base), "BASE_COVARIATES_UNIDENTIFIABLE_NO_SILENT_BASE_CATEGORY_REMOVAL")
            ordered <- c(base, mediator_design_columns)
            q <- qr(augmented[, ordered, drop = FALSE], tol = 1e-11)
            keep <- ordered[q$pivot[seq_len(q$rank)]]
            new121_y_assert(all(base %in% keep), "QR_WOULD_REMOVE_PROTECTED_COVARIATE")
            dropped <- setdiff(ordered, keep)
            witness <- list()
            if (length(dropped)) {
                kept_X <- augmented[, keep, drop = FALSE]
                coefficients <- qr.coef(qr(kept_X, tol = 1e-13), augmented[, dropped, drop = FALSE])
                coefficients <- matrix(coefficients, nrow = length(keep), dimnames = list(keep, dropped))
                for (column in dropped) {
                  coefficient <- coefficients[, column]
                  coefficient[setdiff(base, "(Intercept)")] <- 0
                  fitted <- as.numeric(kept_X %*% coefficient)
                  bound <- 128 * .Machine$double.eps * pmax(1, abs(augmented[, column]) + as.numeric(abs(kept_X) %*% 
                    abs(coefficient)))
                  residual <- abs(augmented[, column] - fitted)
                  new121_y_assert(all(is.finite(coefficient)) && all(residual <= bound), paste("NOT_AN_EXACT_MEDIATION_ALGEBRAIC_ALIAS_NO_HIGH_CORRELATION_DROP", 
                    column))
                  witness[[column]] <- list(coefficients = coefficient, fit_max_abs_residual = max(residual), roundoff_bound_max = max(bound), 
                    rule = "MEDIATION_PLUS_INTERCEPT_IDENTITY_AT_ROUNDOFF_CHECKED_EVERY_PREDICTION")
                }
            }
            list(input_columns = colnames(X), keep_columns = setdiff(keep, "(Intercept)"), dropped_columns = dropped, 
                aliases = witness, base_columns = base_columns, approximate_dimension_reduction = FALSE)
        }
        new121_y_apply_alias_plan <- function(X, plan) {
            if ("(Intercept)" %in% colnames(X)) {
                new121_y_assert(all(X[, "(Intercept)"] == 1), "INVALID_PREDICTION_INTERCEPT")
                X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            }
            new121_y_assert(is.matrix(X) && is.numeric(X) && all(is.finite(X)) && identical(colnames(X), plan$input_columns), 
                "FULL_PREDICTION_X_COLUMN_DRIFT")
            augmented <- cbind(`(Intercept)` = rep(1, nrow(X)), X)
            for (column in names(plan$aliases)) {
                coefficient <- plan$aliases[[column]]$coefficients
                kept_X <- augmented[, names(coefficient), drop = FALSE]
                residual <- abs(augmented[, column] - as.numeric(kept_X %*% coefficient))
                bound <- 128 * .Machine$double.eps * pmax(1, abs(augmented[, column]) + as.numeric(abs(kept_X) %*% 
                  abs(coefficient)))
                new121_y_assert(all(residual <= bound), paste("ALGEBRAIC_REPRESENTATION_IDENTITY_FAILED_ON_GENERATED_DATA", 
                  column))
            }
            X[, plan$keep_columns, drop = FALSE]
        }
        new121_fit_y <- function(engine, fullX, risks, weights, cluster, mediator_design_columns, centre_actual = NULL, 
            audit_callback = NULL, diagnostics = FALSE) {
            new121_y_assert(is.environment(engine) && is.function(engine$ff_ce_fit) && isTRUE(engine$ff_aj_status()$ready) && 
                identical(engine$ff_aj_status()$computation_precision, "DOUBLE_ONLY_BINARY64"), "INSTALLED_DOUBLE_EXPONENTIAL_Y_ENGINE_REQUIRED")
            new121_y_assert(is.list(risks) && setequal(names(risks), c("IHD", "T2D")), "PAIRED_IHD_T2D_RISK_FRAMES_REQUIRED")
            new121_y_assert(length(cluster) == nrow(fullX) && !anyNA(cluster), "EXPLICIT_ROBUST_FAMILY_CLUSTER_REQUIRED")
            plan <- new121_y_exact_alias_plan(fullX, mediator_design_columns)
            X <- new121_y_apply_alias_plan(fullX, plan)
            fits <- list()
            for (endpoint in c("IHD", "T2D")) {
                if (is.function(audit_callback)) audit_callback(list(stage = "Y_FITTING", endpoint = endpoint))
                fits[[endpoint]] <- engine$ff_ce_fit(X = X, risk = risks[[endpoint]], w = weights, cluster = cluster, 
                  tau_report = 11, ties = "efron", diagnostics = diagnostics, endpoint_id = if (endpoint == "IHD") "411" else "250.2", 
                  centre_actual = centre_actual, audit_callback = if (is.function(audit_callback)) local({
                    current_endpoint <- endpoint
                    function(record) audit_callback(c(list(endpoint = current_endpoint), record))
                  }) else NULL)
            }
            compiled <- lapply(fits, engine$fd_fp_compile, v4 = engine)
            list(fits = fits, compiled = compiled, alias_plan = plan, method = list(Y = engine$ff_ce_method_spec(), 
                numeric = engine$c39_num_spec(), sources = engine$new121_source_audit, fit_scope = "FRESH_CURRENT_PANEL_BOOTSTRAP_NO_OLD_COEFFICIENT_REUSE", 
                tau = 11, fit_followup_cap = Inf, ties = "efron", robust_family_cluster = TRUE, warnings_policy = "INHERITED_CLASSIFIED_WARNINGS_ONE_ELIGIBLE_NATIVE_NUMERICAL_RETRY_NO_PENALTY", 
                dropped_exact_aliases = plan$dropped_columns))
        }
        new121_y_design <- function(data, members, base, reference = NULL, engine) {
            engine$c39_design(data, members, base, reference)
        }
        new121_y_design_factory <- function(data, members, base, domains, engine, cache_entries = 4L, reference = NULL) {
            new121_y_assert(cache_entries >= 1L && cache_entries == as.integer(cache_entries), "POSITIVE_CACHE_BOUND_REQUIRED")
            columns <- vapply(members, `[[`, "", "column")
            new121_y_assert(setequal(unlist(domains, use.names = FALSE), columns) && !anyDuplicated(unlist(domains, 
                use.names = FALSE)) && !length(intersect(base, columns)), "EXACT_EXPANDED_COORDINATE_DOMAIN_PARTITION_REQUIRED")
            initial <- engine$c39_design(data, members, base, reference)
            reference <- initial$reference
            continuous <- Filter(function(m) m$type == "continuous", members)
            discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
            form <- stats::reformulate(unique(c(base, discrete)))
            mm <- stats::model.matrix(form, data)
            labels <- attr(stats::terms(form), "term.labels")
            assignment <- attr(mm, "assign")
            owner <- setNames(rep("BASE", ncol(initial$X)), colnames(initial$X))
            for (g in names(domains)) for (column in domains[[g]]) {
                member <- members[[match(column, columns)]]
                owned <- if (member$type == "continuous") reference[[member$indicator_id]]$column_names else colnames(mm)[assignment == 
                  match(column, labels)]
                new121_y_assert(all(owned %in% names(owner)), "UNMAPPED_EXPANDED_MEDIATOR_BASIS")
                owner[owned] <- g
            }
            state <- new.env(parent = emptyenv())
            state$cache <- list()
            state$builds <- 0L
            state$hits <- 0L
            basis <- function(d, key, vars, build) {
                values <- d[, vars, drop = FALSE]
                entries <- state$cache[[key]]
                hit <- which(vapply(entries, function(e) identical(e$values, values), logical(1)))
                if (length(hit)) {
                  state$hits <- state$hits + 1L
                  return(entries[[hit[1L]]]$X)
                }
                X <- build()
                state$cache[[key]] <- head(c(list(list(values = values, X = X)), entries), cache_entries)
                state$builds <- state$builds + 1L
                X
            }
            design <- function(b, M, C) {
                new121_y_assert(length(b) == 1L && b %in% 0:1 && nrow(M) == nrow(C), "SCALAR_EXPOSURE_AND_MATCHED_ROWS_REQUIRED")
                d <- C
                d$B <- b
                for (column in columns) d[[column]] <- M[[column]]
                for (column in intersect(names(d), names(data))) if (is.factor(data[[column]])) {
                  ref <- data[[column]]
                  new121_y_assert(all(as.character(d[[column]]) %in% levels(ref)), paste("UNKNOWN_FROZEN_FACTOR_LEVEL", 
                    column))
                  d[[column]] <- factor(as.character(d[[column]]), levels = levels(ref), ordered = is.ordered(ref))
                  if (!is.null(attr(ref, "contrasts"))) attr(d[[column]], "contrasts") <- attr(ref, "contrasts")
                }
                X <- basis(d, "BASE", base, function() stats::model.matrix(stats::reformulate(base), d))
                for (column in discrete) {
                  part <- basis(d, column, column, function() {
                    q <- stats::model.matrix(stats::reformulate(column), d)
                    q[, colnames(q) != "(Intercept)", drop = FALSE]
                  })
                  X <- cbind(X, part)
                }
                for (member in continuous) {
                  column <- member$column
                  ref <- reference[[member$indicator_id]]
                  X <- cbind(X, basis(d, column, column, function() engine$c39_basis(d[[column]], ref, column)))
                }
                new121_y_assert(setequal(colnames(X), initial$fullcolnames) && all(is.finite(X)), "CACHED_FULL_X_DRIFT")
                X[, initial$fullcolnames, drop = FALSE]
            }
            attr(design, "column_owner") <- owner
            attr(design, "reference") <- reference
            attr(design, "fullX") <- initial$X
            attr(design, "stats") <- function() list(basis_builds = state$builds, basis_cache_hits = state$hits)
            design
        }
        new121_y_prediction <- function(fit_bundle, design_fn, engine, domains = NULL, cache_entries = 4L) {
            new121_y_assert(is.function(design_fn) && is.list(fit_bundle$compiled) && setequal(names(fit_bundle$compiled), 
                c("IHD", "T2D")), "FRESH_PAIRED_FITS_AND_DESIGN_CALLBACK_REQUIRED")
            new121_y_assert(cache_entries >= 1L && cache_entries == as.integer(cache_entries), "POSITIVE_CACHE_BOUND_REQUIRED")
            compiled <- fit_bundle$compiled
            plan <- fit_bundle$alias_plan
            owner <- attr(design_fn, "column_owner")
            if (is.null(owner)) owner <- setNames(rep("ALL", length(plan$keep_columns)), plan$keep_columns)
            owner <- c(owner, ethnicityAsian_or_Chinese_pooled = "BASE")
            state <- new.env(parent = emptyenv())
            state$cache <- list()
            state$predictions <- state$builds <- state$hits <- 0L
            lp <- function(X, beta, key) {
                owners <- owner[colnames(X)]
                owners[is.na(owners)] <- "BASE"
                total <- numeric(nrow(X))
                for (g in unique(unname(owners))) {
                  ii <- which(owners == g)
                  part <- X[, ii, drop = FALSE]
                  name <- paste(key, g, sep = "|")
                  entries <- state$cache[[name]]
                  hit <- which(vapply(entries, function(e) identical(e$X, part), logical(1)))
                  if (length(hit)) {
                    z <- entries[[hit[1L]]]$lp
                    state$hits <- state$hits + 1L
                  } else {
                    z <- as.numeric(part %*% beta[ii])
                    state$cache[[name]] <- head(c(list(list(X = part, lp = z)), entries), cache_entries)
                    state$builds <- state$builds + 1L
                  }
                  total <- total + z
                }
                total
            }
            predict_all <- function(b, M, C, validate = FALSE, weights = NULL) {
                X <- new121_y_apply_alias_plan(design_fn(b = b, M = M, C = C), plan)
                probabilities <- list()
                diffs <- numeric()
                for (endpoint in c("IHD", "T2D")) {
                  cfit <- compiled[[endpoint]]
                  prepared <- engine$fd_fp_prepare(X, cfit, engine)
                  lps <- vapply(1:2, function(j) lp(prepared$designs[[j]], cfit$beta[[j]], paste(endpoint, j, sep = "|")), 
                    numeric(nrow(X)))
                  rr <- exp(matrix(lps, nrow = nrow(X), ncol = 2L))
                  new121_y_assert(all(is.finite(rr) & rr > 0), "NONFINITE_RISK_MULTIPLIER_NO_CLIPPING")
                  p <- matrix(NA_real_, nrow(X), 3L, dimnames = list(NULL, c("disease", "death", "survival")))
                  for (g in unique(prepared$group)) {
                    ii <- which(prepared$group == g)
                    p[ii, ] <- engine$ff_aj_backend(cfit$hazards[[g]], rr[ii, , drop = FALSE])
                  }
                  engine$ff_y_probability_check(p)
                  probabilities[[endpoint]] <- p
                  if (validate) {
                    native <- engine$ff_ce_predict_y_matrix(cfit$object, X, 11)
                    new121_y_assert(!engine$ff_ce_bad_warning(attr(native, "warnings")), "GENUINE_NATIVE_PREDICTION_WARNING")
                    diffs <- c(diffs, abs(p - native))
                  }
                }
                state$predictions <- state$predictions + 1L
                risk <- cbind(IHD = probabilities$IHD[, "disease"], T2D = probabilities$T2D[, "disease"])
                validation <- list(performed = validate, engineering_probability_reference = 1e-06, maximum_absolute_probability_difference = if (length(diffs)) max(diffs) else NULL, 
                  mean_absolute_probability_difference = if (length(diffs)) mean(diffs) else NULL, bitwise_agreement_required = FALSE)
                if (validate && !is.null(weights)) {
                  new121_y_assert(is.numeric(weights) && length(weights) == nrow(X) && all(is.finite(weights) & 
                    weights >= 0) && sum(weights) > 0, "INVALID_BLOCK_COMPARISON_WEIGHTS")
                  validation$weighted_mean_absolute_probability_difference <- sum(diffs * rep(weights, 6L))/(sum(weights) * 
                    6)
                }
                if (validate) validation$within_engineering_reference <- max(diffs) <= 1e-06
                list(risk = risk, all_states = probabilities, validation = validation)
            }
            list(predict = function(b, M, C) predict_all(b, M, C)$risk, validate = function(b, M, C, weights = NULL) predict_all(b, 
                M, C, TRUE, weights), stats = function() c(list(predictions = state$predictions, lp_builds = state$builds, 
                lp_cache_hits = state$hits, bounded_entries_per_piece = cache_entries), if (is.function(attr(design_fn, 
                "stats"))) attr(design_fn, "stats")() else list()), design = design_fn, alias_plan = plan)
        }
    }), envir = context)
}

