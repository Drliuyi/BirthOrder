# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

s44yc_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

s44yc_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

s44yc_design <- function(design, members, weights, transform = NULL) {
    X <- attr(design, "fullX")
    reference <- attr(design, "reference")
    s44yc_assert(is.function(design) && is.matrix(X) && all(is.finite(X)) && length(weights) == nrow(X) && all(is.finite(weights) & 
        weights > 0), "Y_CENTER_EXACT_OBSERVED_DESIGN_AND_WEIGHTS_REQUIRED")
    continuous <- Filter(function(m) identical(m$type, "continuous"), members)
    columns <- unname(unlist(lapply(continuous, function(m) reference[[m$indicator_id]]$column_names), use.names = FALSE))
    owner <- attr(design, "column_owner")
    s44yc_assert(length(columns) > 0L && !anyDuplicated(columns) && all(columns %in% colnames(X)) && all(owner[columns] != 
        "BASE") && !any(grepl("^((Intercept)|B$|age$|sex|center|ethnicity|assessment_year|S$)", columns)), "Y_CENTER_ONLY_CONTINUOUS_MEDIATOR_BASIS_ALLOWED")
    offset <- setNames(colSums(X[, columns, drop = FALSE] * (weights/sum(weights))), columns)
    implementation <- s44yc_hash(list(formals(s44yc_design), body(s44yc_design)))
    if (is.null(transform)) 
        transform <- list(schema = "SELECTED44_Y_CONTINUOUS_NS_ORIGIN_SHIFT_V1", columns = columns, offset = offset, 
            original_fullX_hash = s44yc_hash(X), original_reference_hash = s44yc_hash(reference), implementation_hash = implementation, 
            formula = "X_centered[,columns]=X_original[,columns]-frozen_observed_weighted_mean", unchanged_columns = setdiff(colnames(X), 
                columns), coefficient_mapping = "beta_original=beta_centered; H0_original=H0_centered*exp(-sum(offset*beta)); evaluate in centered coordinates", 
            same_ns_function_space = TRUE, no_column_removed = TRUE, no_factor_center_ethnicity_B_C_changed = TRUE, 
            no_penalty_no_clipping_no_outcome_selection = TRUE)
    s44yc_assert(identical(transform$schema, "SELECTED44_Y_CONTINUOUS_NS_ORIGIN_SHIFT_V1") && identical(transform$columns, 
        columns) && identical(transform$original_fullX_hash, s44yc_hash(X)) && identical(transform$original_reference_hash, 
        s44yc_hash(reference)) && identical(transform$implementation_hash, implementation) && identical(as.numeric(transform$offset), 
        as.numeric(offset)), "Y_CENTER_FROZEN_TRANSFORM_DESIGN_WEIGHT_DRIFT")
    apply <- function(xx) {
        s44yc_assert(identical(colnames(xx), colnames(X)) && all(is.finite(xx)), "Y_CENTER_PREDICTION_COLUMNS_OR_FINITE_DRIFT")
        xx[, columns] <- sweep(xx[, columns, drop = FALSE], 2L, transform$offset, `-`)
        s44yc_assert(all(is.finite(xx)), "Y_CENTER_NONFINITE_TRANSFORM_NO_CLIPPING")
        xx
    }
    result <- function(b, M, C) apply(design(b, M, C))
    attributes(result) <- attributes(design)
    attr(result, "fullX") <- apply(X)
    attr(result, "Y_equivalent_transform") <- transform
    result
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        s44yc_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        s44yc_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        s44yc_design <- function(design, members, weights, transform = NULL) {
            X <- attr(design, "fullX")
            reference <- attr(design, "reference")
            s44yc_assert(is.function(design) && is.matrix(X) && all(is.finite(X)) && length(weights) == nrow(X) && 
                all(is.finite(weights) & weights > 0), "Y_CENTER_EXACT_OBSERVED_DESIGN_AND_WEIGHTS_REQUIRED")
            continuous <- Filter(function(m) identical(m$type, "continuous"), members)
            columns <- unname(unlist(lapply(continuous, function(m) reference[[m$indicator_id]]$column_names), 
                use.names = FALSE))
            owner <- attr(design, "column_owner")
            s44yc_assert(length(columns) > 0L && !anyDuplicated(columns) && all(columns %in% colnames(X)) && all(owner[columns] != 
                "BASE") && !any(grepl("^((Intercept)|B$|age$|sex|center|ethnicity|assessment_year|S$)", columns)), 
                "Y_CENTER_ONLY_CONTINUOUS_MEDIATOR_BASIS_ALLOWED")
            offset <- setNames(colSums(X[, columns, drop = FALSE] * (weights/sum(weights))), columns)
            implementation <- s44yc_hash(list(formals(s44yc_design), body(s44yc_design)))
            if (is.null(transform)) transform <- list(schema = "SELECTED44_Y_CONTINUOUS_NS_ORIGIN_SHIFT_V1", columns = columns, 
                offset = offset, original_fullX_hash = s44yc_hash(X), original_reference_hash = s44yc_hash(reference), 
                implementation_hash = implementation, formula = "X_centered[,columns]=X_original[,columns]-frozen_observed_weighted_mean", 
                unchanged_columns = setdiff(colnames(X), columns), coefficient_mapping = "beta_original=beta_centered; H0_original=H0_centered*exp(-sum(offset*beta)); evaluate in centered coordinates", 
                same_ns_function_space = TRUE, no_column_removed = TRUE, no_factor_center_ethnicity_B_C_changed = TRUE, 
                no_penalty_no_clipping_no_outcome_selection = TRUE)
            s44yc_assert(identical(transform$schema, "SELECTED44_Y_CONTINUOUS_NS_ORIGIN_SHIFT_V1") && identical(transform$columns, 
                columns) && identical(transform$original_fullX_hash, s44yc_hash(X)) && identical(transform$original_reference_hash, 
                s44yc_hash(reference)) && identical(transform$implementation_hash, implementation) && identical(as.numeric(transform$offset), 
                as.numeric(offset)), "Y_CENTER_FROZEN_TRANSFORM_DESIGN_WEIGHT_DRIFT")
            apply <- function(xx) {
                s44yc_assert(identical(colnames(xx), colnames(X)) && all(is.finite(xx)), "Y_CENTER_PREDICTION_COLUMNS_OR_FINITE_DRIFT")
                xx[, columns] <- sweep(xx[, columns, drop = FALSE], 2L, transform$offset, `-`)
                s44yc_assert(all(is.finite(xx)), "Y_CENTER_NONFINITE_TRANSFORM_NO_CLIPPING")
                xx
            }
            result <- function(b, M, C) apply(design(b, M, C))
            attributes(result) <- attributes(design)
            attr(result, "fullX") <- apply(X)
            attr(result, "Y_equivalent_transform") <- transform
            result
        }
    }), envir = context)
}

