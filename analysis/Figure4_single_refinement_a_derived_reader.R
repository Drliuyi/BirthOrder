# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

b16_a_derived_function_sha <- function(fn) digest::digest(paste(c(deparse(formals(fn), width.cutoff = 500L), deparse(body(fn), 
    width.cutoff = 500L)), collapse = "\n"), algo = "sha256", serialize = FALSE)

b16_a_derived_install <- function(envir = parent.frame(), transport_path) {
    stopifnot(file.exists(transport_path))
    transport_sha <- ff_a_sha(transport_path)
    original_read <- get("ff_a_read", envir = envir)
    original_selection <- get("ff_a_selection", envir = envir)
    stopifnot(is.function(original_read), is.function(original_selection), is.null(attr(original_selection, "b16_native_selection")), 
        identical(b16_a_derived_function_sha(original_selection), b16_a_derived_native_sha))
    reader <- local({
        native <- original_read
        pin <- transport_sha
        source_path <- transport_path
        function(desc, endpoint, id = NULL, bid = NULL) {
            if (!identical(desc$format, "derived_group_rds")) 
                return(native(desc, endpoint, id, bid))
            ff_a_assert(!isTRUE(desc$missing_object) && isTRUE(desc$individual_weight_vector_retained), "DERIVED_A_COMPLETE_RETAINED_OBJECT_REQUIRED")
            ff_a_assert(identical(ff_a_sha(source_path), pin) && identical(desc$derivation_transport_sha256, pin), 
                "DERIVED_A_TRANSPORT_SOURCE_CHANGED")
            p <- ff_p_rebase(desc$path)
            pp <- ff_p_rebase(desc$a_group_proof_path)
            ff_a_check(p, desc$sha256)
            ff_a_check(pp, desc$proof_receipt_sha256)
            ff_a_assert(file.info(p)$size <= 32 * 1024^2 && file.info(pp)$size <= 1024^2, "DERIVED_A_OBJECT_SIZE_LIMIT")
            z <- readRDS(p)
            proof <- jsonlite::fromJSON(pp, simplifyVector = FALSE)
            ff_a_check(p, desc$sha256)
            ff_a_check(pp, desc$proof_receipt_sha256)
            ff_a_assert(identical(z$schema, "DERIVED_A_GROUP_RDS_V1") && identical(proof$schema, "DERIVED_A_GROUP_RECEIPT_V1") && 
                isTRUE(proof$complete) && identical(z$a_group_id, desc$a_group_id) && identical(proof$a_group_id, 
                desc$a_group_id) && identical(z$indicator_id, id) && identical(as.integer(z$bootstrap_id), as.integer(bid)) && 
                endpoint %in% as.character(unlist(proof$endpoint_ids)), "DERIVED_A_ID_OR_SCHEMA_CHANGED")
            for (k in c("actual_draw_hash", "target_hash", "sample_hash", "source_input_lock_sha256", "source_prepared_sha256")) {
                dk <- if (k == "actual_draw_hash") 
                  "draw_hash"
                else k
                ff_a_assert(identical(z$binding[[k]], desc[[dk]]) && identical(z$binding[[k]], proof$binding[[k]]), 
                  "DERIVED_A_BINDING_CHANGED")
            }
            ff_a_assert(identical(proof$derivation_transport_sha256, pin) && identical(proof$parameters_sha256, 
                desc$sha256) && identical(b16_a_parameter_hash(z$parameters), z$binding$A_parameter_hash) && identical(z$binding, 
                proof$binding), "DERIVED_A_PARAMETERS_OR_PROOF_CHANGED")
            w <- z$selection_weights
            ff_a_assert(is.numeric(w) && length(w) > 0L && length(w) == proof$N_complete && proof$N_complete <= 
                proof$N_target && all(is.finite(w) & w > 0) && identical(ff_a_hash(format(w, digits = 17, scientific = TRUE)), 
                proof$stored_selection_weight_hash) && identical(proof$stored_selection_weight_hash, desc$derived_selection_weight_hash), 
                "DERIVED_A_RETAINED_WEIGHT_CHANGED")
            params <- z$parameters
            attr(params, "b16_derived_envelope") <- list(binding = z$binding, weights = w, proof = proof, desc = desc)
            params
        }
    })
    selection <- local({
        native <- original_selection
        function(ctx, d, column, base, params, desc) {
            if (!identical(desc$format, "derived_group_rds")) 
                return(native(ctx, d, column, base, params, desc))
            z <- attr(params, "b16_derived_envelope")
            ff_a_assert(is.list(z), "DERIVED_A_READER_MUST_VERIFY_FIRST")
            ff_a_assert(!anyNA(ctx$d_full$eid) && !anyDuplicated(as.character(ctx$d_full$eid)) && nrow(d) == z$proof$N_target, 
                "DERIVED_A_PARENT_IDS_OR_TARGET_LENGTH_CHANGED")
            pr <- match(as.character(d$eid), as.character(ctx$d_full$eid))
            ff_a_assert(!anyNA(pr) && identical(column, ctx$features[[desc$indicator_id]]$column), "DERIVED_A_FROZEN_PARENT_OR_COLUMN_CHANGED")
            binding <- b16_a_cache_binding(ctx, desc$indicator_id, as.integer(desc$bootstrap_id), pr, d, base, 
                params)
            ff_a_assert(identical(binding, z$binding), "DERIVED_A_ACTUAL_DRAW_SAMPLE_INPUT_OR_METHOD_CHANGED")
            complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
            ff_a_assert(length(z$weights) == sum(complete), "DERIVED_A_WEIGHT_SAMPLE_LENGTH_CHANGED")
            list(mi = list(complete = complete, w = z$weights, fit = NULL, warnings = character(), diagnostics = data.frame(), 
                saved_coefficients = params$selection_coefficients), origin = "REUSED_NEW_DERIVED_RETAINED_SELECTION_VECTOR", 
                weight_hash = z$proof$stored_selection_weight_hash, old_weight_hash = NA_character_, A_refit_required = FALSE, 
                saved_selection_coefficient_max_abs_difference = NA_real_)
        }
    })
    attr(selection, "b16_native_selection") <- original_selection
    attr(selection, "b16_transport_sha256") <- transport_sha
    assign("ff_a_read", reader, envir = envir)
    assign("ff_a_selection", selection, envir = envir)
    invisible(list(native_selection_sha256 = b16_a_derived_function_sha(original_selection), transport_sha256 = transport_sha))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        b16_a_derived_native_sha <- birth_private_setting("PRIVATE_HASH_1")
        b16_a_derived_function_sha <- function(fn) digest::digest(paste(c(deparse(formals(fn), width.cutoff = 500L), 
            deparse(body(fn), width.cutoff = 500L)), collapse = "\n"), algo = "sha256", serialize = FALSE)
        b16_a_derived_install <- function(envir = parent.frame(), transport_path) {
            stopifnot(file.exists(transport_path))
            transport_sha <- ff_a_sha(transport_path)
            original_read <- get("ff_a_read", envir = envir)
            original_selection <- get("ff_a_selection", envir = envir)
            stopifnot(is.function(original_read), is.function(original_selection), is.null(attr(original_selection, 
                "b16_native_selection")), identical(b16_a_derived_function_sha(original_selection), b16_a_derived_native_sha))
            reader <- local({
                native <- original_read
                pin <- transport_sha
                source_path <- transport_path
                function(desc, endpoint, id = NULL, bid = NULL) {
                  if (!identical(desc$format, "derived_group_rds")) return(native(desc, endpoint, id, bid))
                  ff_a_assert(!isTRUE(desc$missing_object) && isTRUE(desc$individual_weight_vector_retained), "DERIVED_A_COMPLETE_RETAINED_OBJECT_REQUIRED")
                  ff_a_assert(identical(ff_a_sha(source_path), pin) && identical(desc$derivation_transport_sha256, 
                    pin), "DERIVED_A_TRANSPORT_SOURCE_CHANGED")
                  p <- ff_p_rebase(desc$path)
                  pp <- ff_p_rebase(desc$a_group_proof_path)
                  ff_a_check(p, desc$sha256)
                  ff_a_check(pp, desc$proof_receipt_sha256)
                  ff_a_assert(file.info(p)$size <= 32 * 1024^2 && file.info(pp)$size <= 1024^2, "DERIVED_A_OBJECT_SIZE_LIMIT")
                  z <- readRDS(p)
                  proof <- jsonlite::fromJSON(pp, simplifyVector = FALSE)
                  ff_a_check(p, desc$sha256)
                  ff_a_check(pp, desc$proof_receipt_sha256)
                  ff_a_assert(identical(z$schema, "DERIVED_A_GROUP_RDS_V1") && identical(proof$schema, "DERIVED_A_GROUP_RECEIPT_V1") && 
                    isTRUE(proof$complete) && identical(z$a_group_id, desc$a_group_id) && identical(proof$a_group_id, 
                    desc$a_group_id) && identical(z$indicator_id, id) && identical(as.integer(z$bootstrap_id), 
                    as.integer(bid)) && endpoint %in% as.character(unlist(proof$endpoint_ids)), "DERIVED_A_ID_OR_SCHEMA_CHANGED")
                  for (k in c("actual_draw_hash", "target_hash", "sample_hash", "source_input_lock_sha256", "source_prepared_sha256")) {
                    dk <- if (k == "actual_draw_hash") "draw_hash" else k
                    ff_a_assert(identical(z$binding[[k]], desc[[dk]]) && identical(z$binding[[k]], proof$binding[[k]]), 
                      "DERIVED_A_BINDING_CHANGED")
                  }
                  ff_a_assert(identical(proof$derivation_transport_sha256, pin) && identical(proof$parameters_sha256, 
                    desc$sha256) && identical(b16_a_parameter_hash(z$parameters), z$binding$A_parameter_hash) && 
                    identical(z$binding, proof$binding), "DERIVED_A_PARAMETERS_OR_PROOF_CHANGED")
                  w <- z$selection_weights
                  ff_a_assert(is.numeric(w) && length(w) > 0L && length(w) == proof$N_complete && proof$N_complete <= 
                    proof$N_target && all(is.finite(w) & w > 0) && identical(ff_a_hash(format(w, digits = 17, scientific = TRUE)), 
                    proof$stored_selection_weight_hash) && identical(proof$stored_selection_weight_hash, desc$derived_selection_weight_hash), 
                    "DERIVED_A_RETAINED_WEIGHT_CHANGED")
                  params <- z$parameters
                  attr(params, "b16_derived_envelope") <- list(binding = z$binding, weights = w, proof = proof, 
                    desc = desc)
                  params
                }
            })
            selection <- local({
                native <- original_selection
                function(ctx, d, column, base, params, desc) {
                  if (!identical(desc$format, "derived_group_rds")) return(native(ctx, d, column, base, params, 
                    desc))
                  z <- attr(params, "b16_derived_envelope")
                  ff_a_assert(is.list(z), "DERIVED_A_READER_MUST_VERIFY_FIRST")
                  ff_a_assert(!anyNA(ctx$d_full$eid) && !anyDuplicated(as.character(ctx$d_full$eid)) && nrow(d) == 
                    z$proof$N_target, "DERIVED_A_PARENT_IDS_OR_TARGET_LENGTH_CHANGED")
                  pr <- match(as.character(d$eid), as.character(ctx$d_full$eid))
                  ff_a_assert(!anyNA(pr) && identical(column, ctx$features[[desc$indicator_id]]$column), "DERIVED_A_FROZEN_PARENT_OR_COLUMN_CHANGED")
                  binding <- b16_a_cache_binding(ctx, desc$indicator_id, as.integer(desc$bootstrap_id), pr, d, 
                    base, params)
                  ff_a_assert(identical(binding, z$binding), "DERIVED_A_ACTUAL_DRAW_SAMPLE_INPUT_OR_METHOD_CHANGED")
                  complete <- stats::complete.cases(d[, unique(c(base, column)), drop = FALSE])
                  ff_a_assert(length(z$weights) == sum(complete), "DERIVED_A_WEIGHT_SAMPLE_LENGTH_CHANGED")
                  list(mi = list(complete = complete, w = z$weights, fit = NULL, warnings = character(), diagnostics = data.frame(), 
                    saved_coefficients = params$selection_coefficients), origin = "REUSED_NEW_DERIVED_RETAINED_SELECTION_VECTOR", 
                    weight_hash = z$proof$stored_selection_weight_hash, old_weight_hash = NA_character_, A_refit_required = FALSE, 
                    saved_selection_coefficient_max_abs_difference = NA_real_)
                }
            })
            attr(selection, "b16_native_selection") <- original_selection
            attr(selection, "b16_transport_sha256") <- transport_sha
            assign("ff_a_read", reader, envir = envir)
            assign("ff_a_selection", selection, envir = envir)
            invisible(list(native_selection_sha256 = b16_a_derived_function_sha(original_selection), transport_sha256 = transport_sha))
        }
    }), envir = context)
}

