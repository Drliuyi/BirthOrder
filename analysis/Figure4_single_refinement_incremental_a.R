# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

b9_a_put <- function(x, p) {
    ff_a_assert(!file.exists(p), paste("ADD9_PRESERVE_EXISTING_ARTIFACT", p))
    jsonlite::write_json(x, p, auto_unbox = TRUE, digits = 17, null = "null", na = "null")
}

b9_a_canonical <- function(p) {
    prefix <- paste0(ff_p_env("FF_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH")), "/")
    if (startsWith(p, prefix)) 
        paste0(birth_private_setting("PRIVATE_INPUT"), substring(p, nchar(prefix) + 1L))
    else p
}

b9_a_private_path <- function(ctx, relative) {
    ff_a_assert(is.character(relative) && length(relative) == 1L && nzchar(relative) && !grepl("(^/|(^|/)[.][.](/|$))", 
        relative), "ADD9_PRIVATE_RELATIVE_PATH_REQUIRED")
    path <- file.path(ctx$run, relative)
    ff_a_assert(startsWith(normalizePath(path, mustWork = TRUE), paste0(normalizePath(file.path(ctx$run, "private"), 
        mustWork = TRUE), "/")), "ADD9_POLICY_OUTSIDE_NEW_PRIVATE")
    path
}

b9_a_initialize <- function(ctx, root) {
    r <- jsonlite::fromJSON(file.path(root, "config/A_EXTENSION_READINESS.json"), simplifyVector = FALSE)
    ff_a_assert(!isTRUE(r$ready) && identical(r$mode, "INCREMENTAL_TASK_BOUND_REUSE_OR_CERTIFIED_DERIVE") && identical(r$status, 
        "PREPARED") && identical(r$policy, "NO_DEFAULT_MISSING_OBJECT_NO_BLIND_RECOMPUTATION"), "ADD9_INCREMENTAL_POLICY_REQUIRED_NOT_OLD16_CERTIFICATE")
    p <- b9_a_private_path(ctx, r$source_policy_path)
    cp <- b9_a_private_path(ctx, r$source_policy_complete_path)
    cert <- jsonlite::fromJSON(cp, simplifyVector = FALSE)
    policy <- jsonlite::fromJSON(p, simplifyVector = FALSE)
    ff_a_assert(identical(cert$schema, "ADD9_INCREMENTAL_A_POLICY_COMPLETE_V1") && isTRUE(cert$complete) && cert$models_run == 
        0L && !isTRUE(cert$global_A_completion_claimed) && identical(cert$run_id, ctx$run_id) && identical(cert$runtime_sha256, 
        ff_p_sha(file.path(root, "config/runtime.json"))) && identical(cert$policy_sha256, ff_a_sha(p)), "ADD9_SOURCE_POLICY_SEAL_CHANGED")
    for (rel in names(cert$producer_source_sha256)) ff_a_check(file.path(root, rel), cert$producer_source_sha256[[rel]])
    ff_a_assert(identical(policy$schema, "ADD9_INCREMENTAL_A_SOURCE_POLICY_V1") && policy$models_run == 0L && setequal(unlist(policy$indicator_ids), 
        ctx$scope$indicator_id) && identical(as.integer(unlist(policy$bootstrap_ids)), 0:1000) && setequal(unlist(policy$endpoint_ids), 
        unlist(ctx$ff_config$endpoints)) && identical(policy$source_input_lock_sha256, ctx$old_input_lock_sha256) && 
        identical(policy$source_prepared_sha256, ctx$old_prepared_inputs_sha256) && identical(policy$source_fullfollowup_input_lock_sha256, 
        ctx$input_lock_sha256), "ADD9_SOURCE_POLICY_SCOPE_OR_INPUT_CHANGED")
    sourcecode <- ctx$ff_config$source_fullfollowup$code_paths[[ff_p_host()]]
    oldcfg <- ff_p_cfg(sourcecode)
    source_path <- function(path) {
        aliases <- list(`@SOURCE_FULLFOLLOWUP_CODE_ROOT@` = sourcecode, `@V5_CODE_ROOT@` = oldcfg$previous_code_paths[[ff_p_host()]])
        for (alias in names(aliases)) if (startsWith(path, paste0(alias, "/"))) 
            return(file.path(aliases[[alias]], substring(path, nchar(alias) + 2L)))
        ff_p_rebase(path)
    }
    for (path in names(policy$source_evidence_sha256)) ff_a_check(source_path(path), policy$source_evidence_sha256[[path]])
    state <- new.env(parent = emptyenv())
    state$policy <- policy
    state$policy_sha256 <- cert$policy_sha256
    state$policy_complete_sha256 <- ff_a_sha(cp)
    state$rows <- new.env(parent = emptyenv())
    state$draws <- new.env(parent = emptyenv())
    state$draw_rng <- new.env(parent = emptyenv())
    state$source_equivalence <- new.env(parent = emptyenv())
    state$transport_sha256 <- ff_a_sha(file.path(root, "f/53.a_derived_reader.R"))
    ctx$a_legacy_index_count_before_release <- length(ctx$a_reuse_index)
    keep <- vapply(ctx$a_reuse_index, function(d) d$indicator_id %in% ctx$scope$indicator_id, FALSE)
    state$legacy <- ctx$a_reuse_index[keep]
    ctx$a_reuse_index <- NULL
    ff_a_assert(length(state$legacy) == 9L * 101L * 5L, "ADD9_ORIGINAL_BASE_INDEX_EXACT_COVERAGE_REQUIRED")
    ctx$add9_a <- state
    ctx$a_extension_index <- list()
    ctx$a_index_sha256 <- cert$policy_sha256
    ctx$a_base_index_sha256 <- policy$original_A_source_index_sha256
    ctx$a_source_for <- local({
        s <- state
        function(id, bid, endpoint) {
            z <- get0(paste(id, bid, endpoint, sep = "|"), s$rows, inherits = FALSE)
            ff_a_assert(is.list(z), "ADD9_A_GROUP_NOT_YET_AUDITED_NO_BLIND_REFIT")
            z
        }
    })
    ctx
}

b9_a_source_parent <- function(ctx, id, desc, base) {
    path <- if (is.null(desc$archive)) 
        desc$path
    else desc$archive
    if (startsWith(path, paste0(ctx$old_run, "/"))) 
        return(list(compatible = TRUE, kind = "SAME_FROZEN_V5_CONTEXT", source_input_lock_sha256 = ctx$old_input_lock_sha256))
    ff_a_assert(identical(desc$format, "b18_rds"), "ADD9_UNREVIEWED_DONOR_SOURCE_NAMESPACE")
    key <- paste(desc$format, id, sep = "|")
    state <- ctx$add9_a
    if (exists(key, state$source_equivalence, inherits = FALSE)) 
        return(get(key, state$source_equivalence))
    ref <- state$policy$parent_sources[[desc$format]]
    ff_a_assert(is.list(ref), "ADD9_SOURCE_PARENT_METADATA_REQUIRED")
    pp <- ff_p_rebase(ref$path)
    lp <- ff_p_rebase(ref$input_lock_path)
    ff_a_check(lp, ref$input_lock_sha256)
    ff_a_check(pp, ref$sha256)
    src <- readRDS(pp)
    ff_a_check(pp, ref$sha256)
    ca <- ctx$features[[id]]
    column <- ca$column
    vars <- unique(c("eid", base, column, setdiff(all.vars(ff_a_selection_formula()), "complete")))
    sca <- Filter(function(x) identical(x$feature_id, id), src$cfg$features)
    same <- length(sca) == 1L && all(vars %in% names(src$d_full)) && identical(ctx$d_full[, vars, drop = FALSE], 
        src$d_full[, vars, drop = FALSE])
    if (same) 
        for (k in c("column", "prepared_transform", "prepared_unit", "kind")) if (!is.null(sca[[1L]][[k]])) 
            same <- same && identical(sca[[1L]][[k]], ca[[k]])
    if (same) {
        target <- src$eligibility[[sca[[1L]]$target_key]]
        same <- all(vapply(as.character(unlist(ctx$ff_config$endpoints)), function(e) identical(target & src$endpoints[[e]]$eligibility, 
            ctx$eligibility[[id]] & ctx$endpoints[[e]]$eligibility), FALSE))
    }
    answer <- list(compatible = same, kind = if (same) "EXACT_SOURCE_PARENT_VALUES_FACTORS_FEATURE_AND_FIVE_TARGETS" else "NOT_CERTIFIED_FOR_TASK_SOURCE_PARENT_VALUES_FEATURE_OR_TARGETS_DIFFER", 
        source = ref)
    assign(key, answer, state$source_equivalence)
    rm(src)
    invisible(gc(FALSE))
    answer
}

b9_a_candidates <- function(ctx, id, bid) {
    s <- ctx$add9_a
    out <- list()
    checks <- list()
    if (bid <= 100L) 
        for (e in as.character(unlist(ctx$ff_config$endpoints))) {
            d <- s$legacy[[paste(id, bid, e, sep = "|")]]
            ff_a_assert(is.list(d), "ADD9_BASE_CELL_NOT_IN_FROZEN_INDEX")
            out[[length(out) + 1L]] <- d
        }
    tasks <- s$policy$donor_tasks[[paste(id, bid, sep = "|")]]
    for (t in tasks) {
        dirs <- vapply(t$output_dirs, ff_p_rebase, "")
        present <- dirs[dir.exists(dirs)]
        ff_a_assert(length(present) <= 1L, "ADD9_AMBIGUOUS_DONOR_DIRECTORIES")
        if (!length(present)) {
            checks[[length(checks) + 1L]] <- list(task = t, status = "EXACT_CATALOG_PATH_NOT_PRESENT_NOT_GLOBAL_ABSENCE")
            next
        }
        directory <- present[[1L]]
        v5 <- identical(t$kind, "V5_UNIT")
        if (identical(t$kind, "COMPLETION_UNIT")) {
            sp <- file.path(directory, "COMPLETE.json")
            if (!file.exists(sp)) {
                checks[[length(checks) + 1L]] <- list(task = t, status = "COMPLETION_UNIT_UNSEALED_NOT_CERTIFIED_FOR_TASK")
                next
            }
            sh <- ff_a_sha(sp)
            seal <- jsonlite::fromJSON(sp, simplifyVector = FALSE)
            ff_a_assert(isTRUE(seal$complete) && identical(seal$indicator_id, id) && seal$bootstrap_id == bid, 
                "ADD9_COMPLETION_UNIT_SEAL_CHANGED")
            rp <- file.path(directory, "RESULTS.json")
            ff_a_check(rp, seal$sha256)
            result <- jsonlite::fromJSON(rp, simplifyVector = FALSE)
            ff_a_check(sp, sh)
            for (e in as.character(unlist(ctx$ff_config$endpoints))) {
                rr <- result$rows[[e]]
                if (is.null(rr$new_A_checkpoint_sha256)) 
                  next
                ap <- file.path(ff_p_rebase(t$completion_run), "private/new_A", sprintf("%s_B%04d_%s.rds", id, 
                  bid, e))
                ff_a_check(ap, rr$new_A_checkpoint_sha256)
                out[[length(out) + 1L]] <- list(indicator_id = id, bootstrap_id = bid, endpoint_id = e, format = "add9_completion_retained_rds", 
                  path = ap, sha256 = rr$new_A_checkpoint_sha256, sealed_result = rr, unit_seal_path = b9_a_canonical(sp), 
                  unit_seal_sha256 = sh, unit_results_sha256 = seal$sha256)
            }
            checks[[length(checks) + 1L]] <- list(task = t, status = "COMPLETION_RETAINED_A_CANDIDATES_INSPECTED", 
                unit_seal_sha256 = sh, results_sha256 = seal$sha256)
            next
        }
        sealpath <- file.path(directory, if (v5) 
            "V5_COMPLETE.json"
        else "COMPLETE.json")
        if (!file.exists(sealpath)) {
            checks[[length(checks) + 1L]] <- list(task = t, status = "INCOMPLETE_NO_SEAL_NOT_CERTIFIED_FOR_TASK")
            next
        }
        sealsha <- ff_a_sha(sealpath)
        seal <- jsonlite::fromJSON(sealpath, simplifyVector = FALSE)
        ff_a_assert(isTRUE(seal$complete), "ADD9_DONOR_INCOMPLETE_SEAL")
        arts <- seal$artifacts_sha256
        dp <- file.path(directory, "DRAW.json")
        ff_a_assert(is.list(arts) && !is.null(arts[["DRAW.json"]]), "ADD9_DONOR_DRAW_SEAL_REQUIRED")
        ff_a_check(dp, arts[["DRAW.json"]])
        dr <- jsonlite::fromJSON(dp, simplifyVector = FALSE)
        ff_a_check(sealpath, sealsha)
        check <- list(task = t, seal_path = b9_a_canonical(sealpath), seal_sha256 = sealsha, draw = dr)
        if (is.null(arts[["PARAMETERS.rds"]])) {
            check$status <- "SEALED_PARAMETER_NOT_RETAINED_NOT_CERTIFIED_FOR_TASK"
            checks[[length(checks) + 1L]] <- check
            next
        }
        pp <- file.path(directory, "PARAMETERS.rds")
        ff_a_check(pp, arts[["PARAMETERS.rds"]])
        eps <- if (v5) 
            as.character(unlist(ctx$ff_config$endpoints))
        else t$endpoint_id
        for (e in eps) out[[length(out) + 1L]] <- list(indicator_id = id, bootstrap_id = bid, endpoint_id = e, 
            format = if (v5) "v5_rds" else "b18_rds", path = pp, sha256 = arts[["PARAMETERS.rds"]], draw_hash = dr$draw_hash, 
            target_hash = dr$target_hash, sample_hash = dr$sample_hash, binding_reconstruction_required = TRUE, 
            source_input_lock_sha256 = ctx$old_input_lock_sha256, source_prepared_sha256 = ctx$old_prepared_inputs_sha256, 
            individual_weight_vector_retained = FALSE, selection_weight_hash = NULL, catalog_seal_path = b9_a_canonical(sealpath), 
            catalog_seal_sha256 = sealsha)
        check$status <- "SEALED_CANDIDATE_REQUIRES_ACTUAL_TASK_BINDING"
        checks[[length(checks) + 1L]] <- check
    }
    list(descriptors = out, catalog_checks = checks)
}

b9_a_completion_donor <- function(ctx, id, bid, ds, pr, d, base) {
    ff_a_check(ds$path, ds$sha256)
    aa <- readRDS(ds$path)
    ff_a_check(ds$path, ds$sha256)
    rr <- ds$sealed_result
    ca <- ctx$features[[id]]
    ff_a_assert(identical(rr$indicator_id, id) && rr$bootstrap_id == bid && identical(rr$method_hash, ctx$frozen_A_method_hash) && 
        identical(rr$source_hash, ctx$frozen_A_source_hash) && identical(rr$input_hash, ctx$input_lock_sha256) && 
        identical(aa$input_hash, ctx$input_lock_sha256) && identical(aa$draw$draw_hash, ctx$actual_draw_hash) && 
        identical(rr$draw$draw_hash, ctx$actual_draw_hash), "ADD9_COMPLETION_DONOR_METHOD_INPUT_OR_ACTUAL_DRAW_CHANGED")
    cc <- complete.cases(d[, unique(c(base, ca$column)), drop = FALSE])
    dd <- d[cc, , drop = FALSE]
    wh <- ff_a_hash(format(aa$selection_weight, digits = 17, scientific = TRUE))
    ff_a_assert(identical(aa$target_hash, ff_a_hash(pr)) && identical(aa$sample_hash, ff_a_hash(pr[cc])) && identical(rr$target_hash, 
        aa$target_hash) && identical(rr$sample_hash, aa$sample_hash) && identical(rr$selection_weight_hash, wh) && 
        identical(as.character(aa$ledger$selection_weight_hash), wh) && is.numeric(aa$selection_weight) && length(aa$selection_weight) == 
        sum(cc) && all(is.finite(aa$selection_weight) & aa$selection_weight > 0), "ADD9_COMPLETION_RETAINED_WEIGHT_OR_SAMPLE_CHANGED")
    if (!all(cc)) {
        co <- aa$selection_coefficients
        ff_a_assert(is.numeric(co) && length(co) > 0L && all(is.finite(co)), "ADD9_COMPLETION_SELECTION_COEFFICIENTS_REQUIRED")
        sd <- d
        sd$complete <- cc
        X <- model.matrix(ff_a_selection_formula(), model.frame(ff_a_selection_formula(), sd, na.action = na.fail, 
            drop.unused.levels = FALSE))
        ff_a_assert(identical(names(co), colnames(X)), "ADD9_COMPLETION_SELECTION_DESIGN_CHANGED")
        eta <- rep(0, nrow(X))
        for (j in seq_len(ncol(X))) eta <- eta + X[, j] * co[j]
        reference <- 1/binomial()$linkinv(eta)[cc]
        ff_a_assert(all(is.finite(reference) & reference > 0) && all(abs(reference - aa$selection_weight)/pmax(abs(reference), 
            abs(aa$selection_weight)) <= 1e-12), "ADD9_COMPLETION_SELECTION_NUMERIC_CHECK_FAILED")
    }
    else ff_a_assert(identical(aa$selection_weight, rep(1, sum(cc))), "ADD9_COMPLETION_ALL_COMPLETE_WEIGHTS_CHANGED")
    alpha <- unlist(aa$alpha, use.names = TRUE)
    ff_a_assert(ca$column %in% names(alpha) && length(alpha) == 1L && all(is.finite(alpha)), "ADD9_COMPLETION_ALPHA_MISSING_OR_CHANGED")
    params <- list(alpha = alpha, mediator = NULL, mediators = ca$column, selection_coefficients = aa$selection_coefficients, 
        selection_formula = deparse(ff_a_selection_formula()), base_terms = base, factor_levels = lapply(d[, vapply(d, 
            is.factor, FALSE), drop = FALSE], levels), columns = colnames(model.matrix(reformulate(base), dd)), 
        spline_reference = aa$reference_basis)
    m0 <- dd[[ca$column]] - dd$B * alpha[[ca$column]]
    m1 <- m0 + alpha[[ca$column]]
    lo <- if (is.null(ca$physical_lower_bound)) 
        -Inf
    else ca$physical_lower_bound
    hi <- if (is.null(ca$physical_upper_bound)) 
        Inf
    else ca$physical_upper_bound
    ff_a_assert(all(is.finite(m0) & is.finite(m1) & m0 >= lo & m0 <= hi & m1 >= lo & m1 <= hi), "SCENARIO_OUTSIDE_PHYSICAL_DOMAIN_NO_CLIPPING")
    list(params = params, answer = list(mi = list(complete = cc, w = aa$selection_weight), parent_rows = pr[cc], 
        M0 = m0, M1 = m1), descriptor = ds)
}

b9_a_rederivation_reason <- function(message) {
    if (message %in% c("FF_A_MISSING_SELECTION_RECOMPUTE_NOT_AUTHORIZED", "FF_A_MISSING_ALPHA_RECOMPUTE_NOT_AUTHORIZED", 
        "FF_A_MISSING_CATEGORY_A_RECOMPUTE_NOT_AUTHORIZED")) 
        return("AUDITED_MISSING_FROZEN_SELECTION_OR_A_OBJECT")
    if (identical(message, "FF_A_SELECTION_WEIGHT_HASH_MISMATCH A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED")) 
        return("A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED")
    NULL
}

b9_a_compare_donors <- function(reference, answer) {
    ff_a_assert(identical(reference$parent_rows, answer$parent_rows) && identical(reference$mi$complete, answer$mi$complete) && 
        length(reference$mi$w) == length(answer$mi$w) && all(abs(reference$mi$w - answer$mi$w)/pmax(abs(reference$mi$w), 
        abs(answer$mi$w)) <= 1e-12) && max(abs(reference$M0 - answer$M0), abs(reference$M1 - answer$M1)) < 1e-12, 
        "ADD9_CERTIFIED_DONOR_WEIGHT_OR_SCENARIO_DIVERGENCE")
    TRUE
}

b9_a_materialize <- function(ctx, id, bid, ep, pr, d, base) {
    s <- ctx$add9_a
    ff_a_assert(is.environment(s), "ADD9_INCREMENTAL_STATE_REQUIRED")
    key <- paste(id, bid, ep, sep = "|")
    if (exists(key, s$rows, inherits = FALSE)) 
        return(invisible(TRUE))
    endpoints <- as.character(unlist(ctx$ff_config$endpoints))
    ca <- ctx$features[[id]]
    ff_a_assert(id %in% ctx$scope$indicator_id && bid %in% 0:1000 && is.numeric(d[[ca$column]]), "ADD9_EXACT_CONTINUOUS_TASK_REQUIRED")
    masks <- lapply(endpoints, function(e) ctx$eligibility[[id]] & ctx$endpoints[[e]]$eligibility)
    ff_a_assert(all(vapply(masks, function(x) identical(x, masks[[1L]]), FALSE)), "ADD9_FIVE_ENDPOINT_A_TARGETS_DIFFER_NO_SHARING")
    avars <- unique(c("eid", base, ca$column, setdiff(all.vars(ff_a_selection_formula()), "complete")))
    d <- d[, avars, drop = FALSE]
    pre <- b16_a_cache_binding(ctx, id, bid, pr, d, base, NULL)
    gid <- ff_a_hash(c("A_GROUP_V1", pre$source_input_lock_sha256, id, as.character(bid), pre$actual_draw_hash, 
        pre$target_hash, pre$sample_hash, pre$base_hash, pre$feature_contract_hash))
    parent <- file.path(ctx$run, "private/a_incremental")
    dir.create(parent, showWarnings = FALSE, recursive = TRUE)
    group <- file.path(parent, gid)
    completepath <- file.path(group, "COMPLETE.json")
    if (file.exists(completepath)) {
        seal <- jsonlite::fromJSON(completepath, simplifyVector = FALSE)
        ff_a_assert(identical(seal$a_group_id, gid) && identical(seal$source_policy_sha256, s$policy_sha256), "ADD9_EXISTING_GROUP_POLICY_OR_BINDING_CHANGED")
        ff_a_check(file.path(group, "DESCRIPTORS.json"), seal$descriptors_sha256)
        rows <- jsonlite::fromJSON(file.path(group, "DESCRIPTORS.json"), simplifyVector = FALSE)$rows
        ff_a_assert(length(rows) == 5L && setequal(vapply(rows, function(z) z$endpoint_id, ""), endpoints), "ADD9_EXISTING_GROUP_FIVE_DESCRIPTORS_REQUIRED")
        for (z in rows) {
            ff_a_check(ff_p_rebase(z$path), z$sha256)
            ff_a_check(ff_p_rebase(z$a_group_proof_path), z$proof_receipt_sha256)
            ff_a_read(z, z$endpoint_id, id, bid)
            assign(paste(id, bid, z$endpoint_id, sep = "|"), z, s$rows)
        }
        return(invisible(TRUE))
    }
    ff_a_assert(!dir.exists(group) && dir.create(group), "ADD9_PARTIAL_OR_CLAIMED_A_GROUP_PRESERVED_REVIEW_REQUIRED")
    b9_a_put(list(schema = "ADD9_A_GROUP_OWNER_V1", a_group_id = gid, indicator_id = id, bootstrap_id = bid, host = ff_p_host(), 
        pid = Sys.getpid(), source_policy_sha256 = s$policy_sha256), file.path(group, "OWNER.json"))
    tryCatch({
        ff_w_stage(ctx, "A_DONOR_AUDIT", ep)
        catalog <- b9_a_candidates(ctx, id, bid)
        audits <- list()
        donors <- list()
        native <- attr(ff_a_scenarios, "b16_native_a")
        ff_a_assert(is.function(native), "ADD9_FROZEN_A_FUNCTION_REQUIRED")
        for (i in seq_along(catalog$descriptors)) {
            ds <- catalog$descriptors[[i]]
            e <- as.character(ds$endpoint_id)
            check <- list(descriptor = ds)
            if (identical(ds$format, "add9_completion_retained_rds")) {
                donor <- b9_a_completion_donor(ctx, id, bid, ds, pr, d, base)
                check$status <- "CERTIFIED_COMPLETION_ACTUAL_RETAINED_VECTOR_DONOR"
                check$A_parameter_hash <- b16_a_parameter_hash(donor$params)
                if (length(donors)) {
                  check$certified_donor_numeric_outputs_match <- b9_a_compare_donors(donors[[1L]]$answer, donor$answer)
                  donor$answer <- NULL
                }
                donors[[length(donors) + 1L]] <- donor
                audits[[i]] <- check
                next
            }
            if (isTRUE(ds$missing_object)) {
                ff_a_assert(is.character(ds$reason) && nzchar(ds$reason), "ADD9_MISSING_SOURCE_REASON_REQUIRED")
                check$status <- "AUDITED_MISSING_FROZEN_OBJECT"
                audits[[i]] <- check
                next
            }
            params <- ff_a_read(ds, e, id, bid)
            source <- b9_a_source_parent(ctx, id, ds, base)
            check$source_parent <- source
            if (!isTRUE(source$compatible)) {
                check$status <- "NOT_CERTIFIED_FOR_TASK"
                check$reason <- source$kind
                audits[[i]] <- check
                next
            }
            ff_a_assert(identical(ds$source_input_lock_sha256, ctx$old_input_lock_sha256) && identical(ds$source_prepared_sha256, 
                ctx$old_prepared_inputs_sha256), "ADD9_DONOR_FROZEN_INPUT_DRIFT")
            if (!identical(ds$draw_hash, ctx$actual_draw_hash)) {
                check$status <- "NOT_CERTIFIED_FOR_TASK"
                check$reason <- "SEALED_DONOR_ACTUAL_DRAW_DIFFERS"
                audits[[i]] <- check
                next
            }
            for (k in c("target_hash", "sample_hash")) {
                if (b16_a_cache_is_hash(ds[[k]])) 
                  ff_a_assert(identical(ds[[k]], pre[[k]]), "ADD9_DONOR_TARGET_OR_SAMPLE_DRIFT")
                else {
                  ff_a_assert(isTRUE(ds$binding_reconstruction_required), "ADD9_DONOR_MISSING_BINDING_NOT_AUDITED")
                  ds[[k]] <- pre[[k]]
                }
            }
            probe <- ctx
            probe$allow_missing_a_recompute <- FALSE
            probe$a_source_for <- local({
                z <- ds
                function(id, bid, endpoint) z
            })
            answer <- tryCatch(native(probe, id, bid, e, pr, d, base), error = function(err) err)
            if (inherits(answer, "error") && identical(conditionMessage(answer), "FF_A_SELECTION_WEIGHT_HASH_MISMATCH A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED")) {
                ff_a_assert(!isTRUE(ds$individual_weight_vector_retained) && is.list(params) && length(params$selection_coefficients) > 
                  0L, "ADD9_ROUNDOFF_REQUIRES_SEALED_SELECTION_COEFFICIENTS")
                numeric_ds <- ds
                numeric_ds$selection_weight_hash <- NULL
                probe$a_source_for <- local({
                  z <- numeric_ds
                  function(id, bid, endpoint) z
                })
                answer <- native(probe, id, bid, e, pr, d, base)
                numeric_check <- b16_a_weight_output_check(ctx, answer, params, list(weight_proof_kind = "RECONSTRUCTED_SEALED_SELECTION_COEFFICIENTS"), 
                  ds$selection_weight_hash, ff_a_hash(format(answer$mi$w, digits = 17, scientific = TRUE)), d, 
                  ca$column)
                check$numeric_reconstruction_check <- numeric_check
                check$historical_weight_identity_claimed <- FALSE
            }
            if (inherits(answer, "error")) {
                reason <- b9_a_rederivation_reason(conditionMessage(answer))
                ff_a_assert(!is.null(reason), conditionMessage(answer))
                check$status <- "NOT_CERTIFIED_FOR_TASK"
                check$reason <- reason
            }
            else {
                ff_a_assert(!is.null(params) && is.list(params), "ADD9_VALID_A_PARAMETER_OBJECT_REQUIRED")
                check$status <- "CERTIFIED_CURRENT_TASK_DONOR"
                check$A_parameter_hash <- b16_a_parameter_hash(params)
                if (length(donors)) 
                  check$certified_donor_numeric_outputs_match <- b9_a_compare_donors(donors[[1L]]$answer, answer)
                donors[[length(donors) + 1L]] <- list(params = params, answer = if (!length(donors)) answer else NULL, 
                  descriptor = ds)
            }
            audits[[i]] <- check
        }
        reason <- if (length(donors)) 
            "REUSED_CERTIFIED_CURRENT_TASK_DONOR"
        else "AUDITED_MISSING_OR_NOT_CERTIFIED_FOR_CURRENT_TASK_UNCHANGED_FROZEN_DERIVATION"
        dr <- get0(as.character(bid), s$draws, inherits = FALSE)
        ff_a_assert(is.list(dr) && identical(dr$draw_hash, ctx$actual_draw_hash) && b16_a_cache_is_hash(dr$actual_draw_private_sha256), 
            "ADD9_ACTUAL_RETAINED_DRAW_REQUIRED")
        sourceproof <- list(schema = "ADD9_A_TASK_SOURCE_AUDIT_V1", models_run = 0L, a_group_id = gid, indicator_id = id, 
            bootstrap_id = bid, endpoint_ids = as.list(endpoints), binding = pre, source_policy_sha256 = s$policy_sha256, 
            source_policy_complete_sha256 = s$policy_complete_sha256, source_disposition = reason, existing_donor_count = length(donors), 
            source_checks = audits, catalog_checks = catalog$catalog_checks, upstream_policy = s$policy$upstream_compatibility[[id]], 
            actual_draw_receipt = dr, global_storage_absence_claim = FALSE)
        sp <- file.path(group, "SOURCE_PROOF.json")
        b9_a_put(sourceproof, sp)
        selection_fits <- A_fits <- 0L
        if (length(donors)) {
            params <- donors[[1L]]$params
            reference <- donors[[1L]]$answer
            mi <- reference$mi
        }
        else {
            ff_a_assert(identical(Sys.getenv("FF_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "ADD9_A_DERIVATION_REQUIRES_EXPLICIT_MANUAL_START")
            rng <- get0(as.character(bid), s$draw_rng, inherits = FALSE)
            ff_a_assert(is.integer(rng) && length(rng) > 0L, "ADD9_FROZEN_POST_DRAW_RNG_STATE_REQUIRED")
            assign(".Random.seed", rng, envir = globalenv())
            ff_w_stage(ctx, "A_FROZEN_DERIVATION", ep)
            mi <- ctx$engine$b3_missing(d, ca$column, base)
            selection_fits <- as.integer(!all(mi$complete))
            ff_a_assert(identical(mi$complete, complete.cases(d[, unique(c(base, ca$column)), drop = FALSE])), 
                "ADD9_FROZEN_SELECTION_SAMPLE_CHANGED")
            fit <- ctx$engine$fp_a_path(d, ca$column, base, mi, ctx, paste0("single_", ca$candidate_id), as.integer(bid))
            A_fits <- 1L
            dd <- d[mi$complete, , drop = FALSE]
            X <- model.matrix(reformulate(base), dd)
            fd <- d[, unique(c(base, ca$column)), drop = FALSE]
            params <- list(alpha = fit$alpha, mediator = NULL, mediators = ca$column, selection_coefficients = if (is.null(mi$fit)) NULL else coef(mi$fit), 
                selection_formula = deparse(ff_a_selection_formula()), base_terms = base, factor_levels = lapply(fd[, 
                  vapply(fd, is.factor, FALSE), drop = FALSE], levels), columns = colnames(X), spline_reference = ctx$engine$b3_design(dd, 
                  ca$column, base)$reference)
            physical <- ctx$engine$fp_physical(d, ca$column, mi, fit$alpha, ctx)
            ff_a_assert(all(physical$status == "VALID"), "ADD9_FROZEN_A_PHYSICAL_DIAGNOSTIC_FAILED")
            reference <- list(mi = mi, parent_rows = pr[mi$complete], M0 = dd[[ca$column]] - dd$B * params$alpha[[ca$column]])
            reference$M1 <- reference$M0 + params$alpha[[ca$column]]
        }
        binding <- b16_a_cache_binding(ctx, id, bid, pr, d, base, params)
        weight_hash <- ff_a_hash(format(mi$w, digits = 17, scientific = TRUE))
        rp <- file.path(group, "PARAMETERS.private.rds")
        envelope <- list(schema = "DERIVED_A_GROUP_RDS_V1", a_group_id = gid, indicator_id = id, bootstrap_id = as.integer(bid), 
            binding = binding, parameters = params, selection_weights = mi$w)
        ff_a_assert(!file.exists(rp), "ADD9_PRESERVE_EXISTING_A_VECTOR")
        saveRDS(envelope, rp, compress = FALSE, version = 2)
        proof <- list(schema = "DERIVED_A_GROUP_RECEIPT_V1", complete = TRUE, ready = FALSE, a_group_id = gid, 
            indicator_id = id, bootstrap_id = as.integer(bid), endpoint_ids = as.list(endpoints), binding = binding, 
            source_proof_path = b9_a_canonical(sp), source_proof_sha256 = ff_a_sha(sp), source_disposition = reason, 
            derivation_transport_sha256 = s$transport_sha256, parameters_sha256 = ff_a_sha(rp), stored_selection_weight_hash = weight_hash, 
            historical_selection_weight_hash = NULL, historical_weight_identity_claimed = FALSE, weight_origin = reason, 
            weight_normalization = "NONE", weight_trimming = "NONE", weight_stabilization = "NONE", A_fits = A_fits, 
            selection_fits = selection_fits, Y_fits = 0L, N_target = nrow(d), N_complete = length(mi$w), host = ff_p_host(), 
            task_cache = list(schema = "B16_A_TASK_CACHE_BINDING_V1", certified_for_task_cache = TRUE, a_group_id = gid, 
                endpoint_ids = as.list(endpoints), binding = binding, weight_proof_kind = "RETAINED_INDIVIDUAL_VECTOR", 
                weight_proof_sha256 = b16_a_cache_digest(list(weight_hash = weight_hash, parameters_sha256 = ff_a_sha(rp))), 
                reconstructed_weight_hash = weight_hash, historical_selection_weight_hash = NULL))
        pp <- file.path(group, "PROOF.json")
        b9_a_put(proof, pp)
        ff_w_stage(ctx, "A_RETAINED_CACHE_VALIDATION", ep)
        desc <- list(format = "derived_group_rds", path = b9_a_canonical(rp), sha256 = ff_a_sha(rp), a_group_id = gid, 
            a_group_proof_path = b9_a_canonical(pp), proof_receipt_sha256 = ff_a_sha(pp), derivation_transport_sha256 = s$transport_sha256, 
            derived_selection_weight_hash = weight_hash, indicator_id = id, bootstrap_id = as.integer(bid), draw_hash = binding$actual_draw_hash, 
            target_hash = binding$target_hash, sample_hash = binding$sample_hash, source_input_lock_sha256 = binding$source_input_lock_sha256, 
            source_prepared_sha256 = binding$source_prepared_sha256, missing_object = FALSE, selection_object_missing = FALSE, 
            individual_weight_vector_retained = TRUE, selection_weight_hash = NULL, reason = reason)
        probe <- ctx
        probe$a_source_for <- local({
            z <- desc
            function(id, bid, endpoint) {
                z$endpoint_id <- endpoint
                z
            }
        })
        replay <- native(probe, id, bid, ep, pr, d, base)
        ff_a_assert(identical(replay$mi$w, mi$w) && identical(replay$parent_rows, reference$parent_rows), "ADD9_RETAINED_A_WEIGHT_OR_SAMPLE_ROUNDTRIP_CHANGED")
        difference <- max(abs(replay$M0 - reference$M0), abs(replay$M1 - reference$M1))
        ff_a_assert(is.finite(difference) && difference < 1e-12, "ADD9_RETAINED_A_SCENARIO_ROUNDTRIP_CHANGED")
        rows <- lapply(endpoints, function(e) {
            z <- desc
            z$endpoint_id <- e
            z
        })
        dp <- file.path(group, "DESCRIPTORS.json")
        b9_a_put(list(rows = rows), dp)
        b9_a_put(list(schema = "ADD9_A_GROUP_COMPLETE_V1", complete = TRUE, a_group_id = gid, source_policy_sha256 = s$policy_sha256, 
            descriptors_sha256 = ff_a_sha(dp), A_fits = A_fits, selection_fits = selection_fits, Y_fits = 0L, source_disposition = reason, 
            roundtrip_max_abs_difference = difference, actual_weight_bitwise_roundtrip = TRUE), completepath)
        for (z in rows) assign(paste(id, bid, z$endpoint_id, sep = "|"), z, s$rows)
        invisible(TRUE)
    }, error = function(error) {
        fp <- file.path(group, "FAILED.json")
        if (!file.exists(fp)) 
            b9_a_put(list(a_group_id = gid, indicator_id = id, bootstrap_id = bid, reason = conditionMessage(error), 
                outputs_retained = TRUE, automatic_recompute_authorized = FALSE), fp)
        stop(error)
    })
}

b9_a_install <- function(envir = parent.frame()) {
    cached <- get("ff_a_scenarios", envir = envir)
    native <- attr(cached, "b16_native_a")
    ff_a_assert(is.function(native), "ADD9_INSTALL_REQUIRES_FROZEN_F52_CACHE")
    wrapper <- function(ctx, id, bid, ep, pr, d, base) {
        b9_a_materialize(ctx, id, as.integer(bid), ep, pr, d, base)
        ff_w_stage(ctx, "A_RETAINED_CACHE_VALIDATION", ep)
        cached(ctx, id, bid, ep, pr, d, base)
    }
    attr(wrapper, "b16_native_a") <- native
    assign("ff_a_scenarios", wrapper, envir = envir)
    invisible(TRUE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        b9_a_put <- function(x, p) {
            ff_a_assert(!file.exists(p), paste("ADD9_PRESERVE_EXISTING_ARTIFACT", p))
            jsonlite::write_json(x, p, auto_unbox = TRUE, digits = 17, null = "null", na = "null")
        }
        b9_a_canonical <- function(p) {
            prefix <- paste0(ff_p_env("FF_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH")), "/")
            if (startsWith(p, prefix)) paste0(birth_private_setting("PRIVATE_INPUT"), substring(p, nchar(prefix) + 
                1L)) else p
        }
        b9_a_private_path <- function(ctx, relative) {
            ff_a_assert(is.character(relative) && length(relative) == 1L && nzchar(relative) && !grepl("(^/|(^|/)[.][.](/|$))", 
                relative), "ADD9_PRIVATE_RELATIVE_PATH_REQUIRED")
            path <- file.path(ctx$run, relative)
            ff_a_assert(startsWith(normalizePath(path, mustWork = TRUE), paste0(normalizePath(file.path(ctx$run, 
                "private"), mustWork = TRUE), "/")), "ADD9_POLICY_OUTSIDE_NEW_PRIVATE")
            path
        }
        b9_a_initialize <- function(ctx, root) {
            r <- jsonlite::fromJSON(file.path(root, "config/A_EXTENSION_READINESS.json"), simplifyVector = FALSE)
            ff_a_assert(!isTRUE(r$ready) && identical(r$mode, "INCREMENTAL_TASK_BOUND_REUSE_OR_CERTIFIED_DERIVE") && 
                identical(r$status, "PREPARED") && identical(r$policy, "NO_DEFAULT_MISSING_OBJECT_NO_BLIND_RECOMPUTATION"), 
                "ADD9_INCREMENTAL_POLICY_REQUIRED_NOT_OLD16_CERTIFICATE")
            p <- b9_a_private_path(ctx, r$source_policy_path)
            cp <- b9_a_private_path(ctx, r$source_policy_complete_path)
            cert <- jsonlite::fromJSON(cp, simplifyVector = FALSE)
            policy <- jsonlite::fromJSON(p, simplifyVector = FALSE)
            ff_a_assert(identical(cert$schema, "ADD9_INCREMENTAL_A_POLICY_COMPLETE_V1") && isTRUE(cert$complete) && 
                cert$models_run == 0L && !isTRUE(cert$global_A_completion_claimed) && identical(cert$run_id, ctx$run_id) && 
                identical(cert$runtime_sha256, ff_p_sha(file.path(root, "config/runtime.json"))) && identical(cert$policy_sha256, 
                ff_a_sha(p)), "ADD9_SOURCE_POLICY_SEAL_CHANGED")
            for (rel in names(cert$producer_source_sha256)) ff_a_check(file.path(root, rel), cert$producer_source_sha256[[rel]])
            ff_a_assert(identical(policy$schema, "ADD9_INCREMENTAL_A_SOURCE_POLICY_V1") && policy$models_run == 
                0L && setequal(unlist(policy$indicator_ids), ctx$scope$indicator_id) && identical(as.integer(unlist(policy$bootstrap_ids)), 
                0:1000) && setequal(unlist(policy$endpoint_ids), unlist(ctx$ff_config$endpoints)) && identical(policy$source_input_lock_sha256, 
                ctx$old_input_lock_sha256) && identical(policy$source_prepared_sha256, ctx$old_prepared_inputs_sha256) && 
                identical(policy$source_fullfollowup_input_lock_sha256, ctx$input_lock_sha256), "ADD9_SOURCE_POLICY_SCOPE_OR_INPUT_CHANGED")
            sourcecode <- ctx$ff_config$source_fullfollowup$code_paths[[ff_p_host()]]
            oldcfg <- ff_p_cfg(sourcecode)
            source_path <- function(path) {
                aliases <- list(`@SOURCE_FULLFOLLOWUP_CODE_ROOT@` = sourcecode, `@V5_CODE_ROOT@` = oldcfg$previous_code_paths[[ff_p_host()]])
                for (alias in names(aliases)) if (startsWith(path, paste0(alias, "/"))) return(file.path(aliases[[alias]], 
                  substring(path, nchar(alias) + 2L)))
                ff_p_rebase(path)
            }
            for (path in names(policy$source_evidence_sha256)) ff_a_check(source_path(path), policy$source_evidence_sha256[[path]])
            state <- new.env(parent = emptyenv())
            state$policy <- policy
            state$policy_sha256 <- cert$policy_sha256
            state$policy_complete_sha256 <- ff_a_sha(cp)
            state$rows <- new.env(parent = emptyenv())
            state$draws <- new.env(parent = emptyenv())
            state$draw_rng <- new.env(parent = emptyenv())
            state$source_equivalence <- new.env(parent = emptyenv())
            state$transport_sha256 <- ff_a_sha(file.path(root, "f/53.a_derived_reader.R"))
            ctx$a_legacy_index_count_before_release <- length(ctx$a_reuse_index)
            keep <- vapply(ctx$a_reuse_index, function(d) d$indicator_id %in% ctx$scope$indicator_id, FALSE)
            state$legacy <- ctx$a_reuse_index[keep]
            ctx$a_reuse_index <- NULL
            ff_a_assert(length(state$legacy) == 9L * 101L * 5L, "ADD9_ORIGINAL_BASE_INDEX_EXACT_COVERAGE_REQUIRED")
            ctx$add9_a <- state
            ctx$a_extension_index <- list()
            ctx$a_index_sha256 <- cert$policy_sha256
            ctx$a_base_index_sha256 <- policy$original_A_source_index_sha256
            ctx$a_source_for <- local({
                s <- state
                function(id, bid, endpoint) {
                  z <- get0(paste(id, bid, endpoint, sep = "|"), s$rows, inherits = FALSE)
                  ff_a_assert(is.list(z), "ADD9_A_GROUP_NOT_YET_AUDITED_NO_BLIND_REFIT")
                  z
                }
            })
            ctx
        }
        b9_a_source_parent <- function(ctx, id, desc, base) {
            path <- if (is.null(desc$archive)) desc$path else desc$archive
            if (startsWith(path, paste0(ctx$old_run, "/"))) return(list(compatible = TRUE, kind = "SAME_FROZEN_V5_CONTEXT", 
                source_input_lock_sha256 = ctx$old_input_lock_sha256))
            ff_a_assert(identical(desc$format, "b18_rds"), "ADD9_UNREVIEWED_DONOR_SOURCE_NAMESPACE")
            key <- paste(desc$format, id, sep = "|")
            state <- ctx$add9_a
            if (exists(key, state$source_equivalence, inherits = FALSE)) return(get(key, state$source_equivalence))
            ref <- state$policy$parent_sources[[desc$format]]
            ff_a_assert(is.list(ref), "ADD9_SOURCE_PARENT_METADATA_REQUIRED")
            pp <- ff_p_rebase(ref$path)
            lp <- ff_p_rebase(ref$input_lock_path)
            ff_a_check(lp, ref$input_lock_sha256)
            ff_a_check(pp, ref$sha256)
            src <- readRDS(pp)
            ff_a_check(pp, ref$sha256)
            ca <- ctx$features[[id]]
            column <- ca$column
            vars <- unique(c("eid", base, column, setdiff(all.vars(ff_a_selection_formula()), "complete")))
            sca <- Filter(function(x) identical(x$feature_id, id), src$cfg$features)
            same <- length(sca) == 1L && all(vars %in% names(src$d_full)) && identical(ctx$d_full[, vars, drop = FALSE], 
                src$d_full[, vars, drop = FALSE])
            if (same) for (k in c("column", "prepared_transform", "prepared_unit", "kind")) if (!is.null(sca[[1L]][[k]])) same <- same && 
                identical(sca[[1L]][[k]], ca[[k]])
            if (same) {
                target <- src$eligibility[[sca[[1L]]$target_key]]
                same <- all(vapply(as.character(unlist(ctx$ff_config$endpoints)), function(e) identical(target & 
                  src$endpoints[[e]]$eligibility, ctx$eligibility[[id]] & ctx$endpoints[[e]]$eligibility), FALSE))
            }
            answer <- list(compatible = same, kind = if (same) "EXACT_SOURCE_PARENT_VALUES_FACTORS_FEATURE_AND_FIVE_TARGETS" else "NOT_CERTIFIED_FOR_TASK_SOURCE_PARENT_VALUES_FEATURE_OR_TARGETS_DIFFER", 
                source = ref)
            assign(key, answer, state$source_equivalence)
            rm(src)
            invisible(gc(FALSE))
            answer
        }
        b9_a_candidates <- function(ctx, id, bid) {
            s <- ctx$add9_a
            out <- list()
            checks <- list()
            if (bid <= 100L) for (e in as.character(unlist(ctx$ff_config$endpoints))) {
                d <- s$legacy[[paste(id, bid, e, sep = "|")]]
                ff_a_assert(is.list(d), "ADD9_BASE_CELL_NOT_IN_FROZEN_INDEX")
                out[[length(out) + 1L]] <- d
            }
            tasks <- s$policy$donor_tasks[[paste(id, bid, sep = "|")]]
            for (t in tasks) {
                dirs <- vapply(t$output_dirs, ff_p_rebase, "")
                present <- dirs[dir.exists(dirs)]
                ff_a_assert(length(present) <= 1L, "ADD9_AMBIGUOUS_DONOR_DIRECTORIES")
                if (!length(present)) {
                  checks[[length(checks) + 1L]] <- list(task = t, status = "EXACT_CATALOG_PATH_NOT_PRESENT_NOT_GLOBAL_ABSENCE")
                  next
                }
                directory <- present[[1L]]
                v5 <- identical(t$kind, "V5_UNIT")
                if (identical(t$kind, "COMPLETION_UNIT")) {
                  sp <- file.path(directory, "COMPLETE.json")
                  if (!file.exists(sp)) {
                    checks[[length(checks) + 1L]] <- list(task = t, status = "COMPLETION_UNIT_UNSEALED_NOT_CERTIFIED_FOR_TASK")
                    next
                  }
                  sh <- ff_a_sha(sp)
                  seal <- jsonlite::fromJSON(sp, simplifyVector = FALSE)
                  ff_a_assert(isTRUE(seal$complete) && identical(seal$indicator_id, id) && seal$bootstrap_id == 
                    bid, "ADD9_COMPLETION_UNIT_SEAL_CHANGED")
                  rp <- file.path(directory, "RESULTS.json")
                  ff_a_check(rp, seal$sha256)
                  result <- jsonlite::fromJSON(rp, simplifyVector = FALSE)
                  ff_a_check(sp, sh)
                  for (e in as.character(unlist(ctx$ff_config$endpoints))) {
                    rr <- result$rows[[e]]
                    if (is.null(rr$new_A_checkpoint_sha256)) next
                    ap <- file.path(ff_p_rebase(t$completion_run), "private/new_A", sprintf("%s_B%04d_%s.rds", 
                      id, bid, e))
                    ff_a_check(ap, rr$new_A_checkpoint_sha256)
                    out[[length(out) + 1L]] <- list(indicator_id = id, bootstrap_id = bid, endpoint_id = e, format = "add9_completion_retained_rds", 
                      path = ap, sha256 = rr$new_A_checkpoint_sha256, sealed_result = rr, unit_seal_path = b9_a_canonical(sp), 
                      unit_seal_sha256 = sh, unit_results_sha256 = seal$sha256)
                  }
                  checks[[length(checks) + 1L]] <- list(task = t, status = "COMPLETION_RETAINED_A_CANDIDATES_INSPECTED", 
                    unit_seal_sha256 = sh, results_sha256 = seal$sha256)
                  next
                }
                sealpath <- file.path(directory, if (v5) "V5_COMPLETE.json" else "COMPLETE.json")
                if (!file.exists(sealpath)) {
                  checks[[length(checks) + 1L]] <- list(task = t, status = "INCOMPLETE_NO_SEAL_NOT_CERTIFIED_FOR_TASK")
                  next
                }
                sealsha <- ff_a_sha(sealpath)
                seal <- jsonlite::fromJSON(sealpath, simplifyVector = FALSE)
                ff_a_assert(isTRUE(seal$complete), "ADD9_DONOR_INCOMPLETE_SEAL")
                arts <- seal$artifacts_sha256
                dp <- file.path(directory, "DRAW.json")
                ff_a_assert(is.list(arts) && !is.null(arts[["DRAW.json"]]), "ADD9_DONOR_DRAW_SEAL_REQUIRED")
                ff_a_check(dp, arts[["DRAW.json"]])
                dr <- jsonlite::fromJSON(dp, simplifyVector = FALSE)
                ff_a_check(sealpath, sealsha)
                check <- list(task = t, seal_path = b9_a_canonical(sealpath), seal_sha256 = sealsha, draw = dr)
                if (is.null(arts[["PARAMETERS.rds"]])) {
                  check$status <- "SEALED_PARAMETER_NOT_RETAINED_NOT_CERTIFIED_FOR_TASK"
                  checks[[length(checks) + 1L]] <- check
                  next
                }
                pp <- file.path(directory, "PARAMETERS.rds")
                ff_a_check(pp, arts[["PARAMETERS.rds"]])
                eps <- if (v5) as.character(unlist(ctx$ff_config$endpoints)) else t$endpoint_id
                for (e in eps) out[[length(out) + 1L]] <- list(indicator_id = id, bootstrap_id = bid, endpoint_id = e, 
                  format = if (v5) "v5_rds" else "b18_rds", path = pp, sha256 = arts[["PARAMETERS.rds"]], draw_hash = dr$draw_hash, 
                  target_hash = dr$target_hash, sample_hash = dr$sample_hash, binding_reconstruction_required = TRUE, 
                  source_input_lock_sha256 = ctx$old_input_lock_sha256, source_prepared_sha256 = ctx$old_prepared_inputs_sha256, 
                  individual_weight_vector_retained = FALSE, selection_weight_hash = NULL, catalog_seal_path = b9_a_canonical(sealpath), 
                  catalog_seal_sha256 = sealsha)
                check$status <- "SEALED_CANDIDATE_REQUIRES_ACTUAL_TASK_BINDING"
                checks[[length(checks) + 1L]] <- check
            }
            list(descriptors = out, catalog_checks = checks)
        }
        b9_a_completion_donor <- function(ctx, id, bid, ds, pr, d, base) {
            ff_a_check(ds$path, ds$sha256)
            aa <- readRDS(ds$path)
            ff_a_check(ds$path, ds$sha256)
            rr <- ds$sealed_result
            ca <- ctx$features[[id]]
            ff_a_assert(identical(rr$indicator_id, id) && rr$bootstrap_id == bid && identical(rr$method_hash, ctx$frozen_A_method_hash) && 
                identical(rr$source_hash, ctx$frozen_A_source_hash) && identical(rr$input_hash, ctx$input_lock_sha256) && 
                identical(aa$input_hash, ctx$input_lock_sha256) && identical(aa$draw$draw_hash, ctx$actual_draw_hash) && 
                identical(rr$draw$draw_hash, ctx$actual_draw_hash), "ADD9_COMPLETION_DONOR_METHOD_INPUT_OR_ACTUAL_DRAW_CHANGED")
            cc <- complete.cases(d[, unique(c(base, ca$column)), drop = FALSE])
            dd <- d[cc, , drop = FALSE]
            wh <- ff_a_hash(format(aa$selection_weight, digits = 17, scientific = TRUE))
            ff_a_assert(identical(aa$target_hash, ff_a_hash(pr)) && identical(aa$sample_hash, ff_a_hash(pr[cc])) && 
                identical(rr$target_hash, aa$target_hash) && identical(rr$sample_hash, aa$sample_hash) && identical(rr$selection_weight_hash, 
                wh) && identical(as.character(aa$ledger$selection_weight_hash), wh) && is.numeric(aa$selection_weight) && 
                length(aa$selection_weight) == sum(cc) && all(is.finite(aa$selection_weight) & aa$selection_weight > 
                0), "ADD9_COMPLETION_RETAINED_WEIGHT_OR_SAMPLE_CHANGED")
            if (!all(cc)) {
                co <- aa$selection_coefficients
                ff_a_assert(is.numeric(co) && length(co) > 0L && all(is.finite(co)), "ADD9_COMPLETION_SELECTION_COEFFICIENTS_REQUIRED")
                sd <- d
                sd$complete <- cc
                X <- model.matrix(ff_a_selection_formula(), model.frame(ff_a_selection_formula(), sd, na.action = na.fail, 
                  drop.unused.levels = FALSE))
                ff_a_assert(identical(names(co), colnames(X)), "ADD9_COMPLETION_SELECTION_DESIGN_CHANGED")
                eta <- rep(0, nrow(X))
                for (j in seq_len(ncol(X))) eta <- eta + X[, j] * co[j]
                reference <- 1/binomial()$linkinv(eta)[cc]
                ff_a_assert(all(is.finite(reference) & reference > 0) && all(abs(reference - aa$selection_weight)/pmax(abs(reference), 
                  abs(aa$selection_weight)) <= 1e-12), "ADD9_COMPLETION_SELECTION_NUMERIC_CHECK_FAILED")
            } else ff_a_assert(identical(aa$selection_weight, rep(1, sum(cc))), "ADD9_COMPLETION_ALL_COMPLETE_WEIGHTS_CHANGED")
            alpha <- unlist(aa$alpha, use.names = TRUE)
            ff_a_assert(ca$column %in% names(alpha) && length(alpha) == 1L && all(is.finite(alpha)), "ADD9_COMPLETION_ALPHA_MISSING_OR_CHANGED")
            params <- list(alpha = alpha, mediator = NULL, mediators = ca$column, selection_coefficients = aa$selection_coefficients, 
                selection_formula = deparse(ff_a_selection_formula()), base_terms = base, factor_levels = lapply(d[, 
                  vapply(d, is.factor, FALSE), drop = FALSE], levels), columns = colnames(model.matrix(reformulate(base), 
                  dd)), spline_reference = aa$reference_basis)
            m0 <- dd[[ca$column]] - dd$B * alpha[[ca$column]]
            m1 <- m0 + alpha[[ca$column]]
            lo <- if (is.null(ca$physical_lower_bound)) -Inf else ca$physical_lower_bound
            hi <- if (is.null(ca$physical_upper_bound)) Inf else ca$physical_upper_bound
            ff_a_assert(all(is.finite(m0) & is.finite(m1) & m0 >= lo & m0 <= hi & m1 >= lo & m1 <= hi), "SCENARIO_OUTSIDE_PHYSICAL_DOMAIN_NO_CLIPPING")
            list(params = params, answer = list(mi = list(complete = cc, w = aa$selection_weight), parent_rows = pr[cc], 
                M0 = m0, M1 = m1), descriptor = ds)
        }
        b9_a_rederivation_reason <- function(message) {
            if (message %in% c("FF_A_MISSING_SELECTION_RECOMPUTE_NOT_AUTHORIZED", "FF_A_MISSING_ALPHA_RECOMPUTE_NOT_AUTHORIZED", 
                "FF_A_MISSING_CATEGORY_A_RECOMPUTE_NOT_AUTHORIZED")) return("AUDITED_MISSING_FROZEN_SELECTION_OR_A_OBJECT")
            if (identical(message, "FF_A_SELECTION_WEIGHT_HASH_MISMATCH A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED")) return("A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED")
            NULL
        }
        b9_a_compare_donors <- function(reference, answer) {
            ff_a_assert(identical(reference$parent_rows, answer$parent_rows) && identical(reference$mi$complete, 
                answer$mi$complete) && length(reference$mi$w) == length(answer$mi$w) && all(abs(reference$mi$w - 
                answer$mi$w)/pmax(abs(reference$mi$w), abs(answer$mi$w)) <= 1e-12) && max(abs(reference$M0 - answer$M0), 
                abs(reference$M1 - answer$M1)) < 1e-12, "ADD9_CERTIFIED_DONOR_WEIGHT_OR_SCENARIO_DIVERGENCE")
            TRUE
        }
        b9_a_materialize <- function(ctx, id, bid, ep, pr, d, base) {
            s <- ctx$add9_a
            ff_a_assert(is.environment(s), "ADD9_INCREMENTAL_STATE_REQUIRED")
            key <- paste(id, bid, ep, sep = "|")
            if (exists(key, s$rows, inherits = FALSE)) return(invisible(TRUE))
            endpoints <- as.character(unlist(ctx$ff_config$endpoints))
            ca <- ctx$features[[id]]
            ff_a_assert(id %in% ctx$scope$indicator_id && bid %in% 0:1000 && is.numeric(d[[ca$column]]), "ADD9_EXACT_CONTINUOUS_TASK_REQUIRED")
            masks <- lapply(endpoints, function(e) ctx$eligibility[[id]] & ctx$endpoints[[e]]$eligibility)
            ff_a_assert(all(vapply(masks, function(x) identical(x, masks[[1L]]), FALSE)), "ADD9_FIVE_ENDPOINT_A_TARGETS_DIFFER_NO_SHARING")
            avars <- unique(c("eid", base, ca$column, setdiff(all.vars(ff_a_selection_formula()), "complete")))
            d <- d[, avars, drop = FALSE]
            pre <- b16_a_cache_binding(ctx, id, bid, pr, d, base, NULL)
            gid <- ff_a_hash(c("A_GROUP_V1", pre$source_input_lock_sha256, id, as.character(bid), pre$actual_draw_hash, 
                pre$target_hash, pre$sample_hash, pre$base_hash, pre$feature_contract_hash))
            parent <- file.path(ctx$run, "private/a_incremental")
            dir.create(parent, showWarnings = FALSE, recursive = TRUE)
            group <- file.path(parent, gid)
            completepath <- file.path(group, "COMPLETE.json")
            if (file.exists(completepath)) {
                seal <- jsonlite::fromJSON(completepath, simplifyVector = FALSE)
                ff_a_assert(identical(seal$a_group_id, gid) && identical(seal$source_policy_sha256, s$policy_sha256), 
                  "ADD9_EXISTING_GROUP_POLICY_OR_BINDING_CHANGED")
                ff_a_check(file.path(group, "DESCRIPTORS.json"), seal$descriptors_sha256)
                rows <- jsonlite::fromJSON(file.path(group, "DESCRIPTORS.json"), simplifyVector = FALSE)$rows
                ff_a_assert(length(rows) == 5L && setequal(vapply(rows, function(z) z$endpoint_id, ""), endpoints), 
                  "ADD9_EXISTING_GROUP_FIVE_DESCRIPTORS_REQUIRED")
                for (z in rows) {
                  ff_a_check(ff_p_rebase(z$path), z$sha256)
                  ff_a_check(ff_p_rebase(z$a_group_proof_path), z$proof_receipt_sha256)
                  ff_a_read(z, z$endpoint_id, id, bid)
                  assign(paste(id, bid, z$endpoint_id, sep = "|"), z, s$rows)
                }
                return(invisible(TRUE))
            }
            ff_a_assert(!dir.exists(group) && dir.create(group), "ADD9_PARTIAL_OR_CLAIMED_A_GROUP_PRESERVED_REVIEW_REQUIRED")
            b9_a_put(list(schema = "ADD9_A_GROUP_OWNER_V1", a_group_id = gid, indicator_id = id, bootstrap_id = bid, 
                host = ff_p_host(), pid = Sys.getpid(), source_policy_sha256 = s$policy_sha256), file.path(group, 
                "OWNER.json"))
            tryCatch({
                ff_w_stage(ctx, "A_DONOR_AUDIT", ep)
                catalog <- b9_a_candidates(ctx, id, bid)
                audits <- list()
                donors <- list()
                native <- attr(ff_a_scenarios, "b16_native_a")
                ff_a_assert(is.function(native), "ADD9_FROZEN_A_FUNCTION_REQUIRED")
                for (i in seq_along(catalog$descriptors)) {
                  ds <- catalog$descriptors[[i]]
                  e <- as.character(ds$endpoint_id)
                  check <- list(descriptor = ds)
                  if (identical(ds$format, "add9_completion_retained_rds")) {
                    donor <- b9_a_completion_donor(ctx, id, bid, ds, pr, d, base)
                    check$status <- "CERTIFIED_COMPLETION_ACTUAL_RETAINED_VECTOR_DONOR"
                    check$A_parameter_hash <- b16_a_parameter_hash(donor$params)
                    if (length(donors)) {
                      check$certified_donor_numeric_outputs_match <- b9_a_compare_donors(donors[[1L]]$answer, donor$answer)
                      donor$answer <- NULL
                    }
                    donors[[length(donors) + 1L]] <- donor
                    audits[[i]] <- check
                    next
                  }
                  if (isTRUE(ds$missing_object)) {
                    ff_a_assert(is.character(ds$reason) && nzchar(ds$reason), "ADD9_MISSING_SOURCE_REASON_REQUIRED")
                    check$status <- "AUDITED_MISSING_FROZEN_OBJECT"
                    audits[[i]] <- check
                    next
                  }
                  params <- ff_a_read(ds, e, id, bid)
                  source <- b9_a_source_parent(ctx, id, ds, base)
                  check$source_parent <- source
                  if (!isTRUE(source$compatible)) {
                    check$status <- "NOT_CERTIFIED_FOR_TASK"
                    check$reason <- source$kind
                    audits[[i]] <- check
                    next
                  }
                  ff_a_assert(identical(ds$source_input_lock_sha256, ctx$old_input_lock_sha256) && identical(ds$source_prepared_sha256, 
                    ctx$old_prepared_inputs_sha256), "ADD9_DONOR_FROZEN_INPUT_DRIFT")
                  if (!identical(ds$draw_hash, ctx$actual_draw_hash)) {
                    check$status <- "NOT_CERTIFIED_FOR_TASK"
                    check$reason <- "SEALED_DONOR_ACTUAL_DRAW_DIFFERS"
                    audits[[i]] <- check
                    next
                  }
                  for (k in c("target_hash", "sample_hash")) {
                    if (b16_a_cache_is_hash(ds[[k]])) ff_a_assert(identical(ds[[k]], pre[[k]]), "ADD9_DONOR_TARGET_OR_SAMPLE_DRIFT") else {
                      ff_a_assert(isTRUE(ds$binding_reconstruction_required), "ADD9_DONOR_MISSING_BINDING_NOT_AUDITED")
                      ds[[k]] <- pre[[k]]
                    }
                  }
                  probe <- ctx
                  probe$allow_missing_a_recompute <- FALSE
                  probe$a_source_for <- local({
                    z <- ds
                    function(id, bid, endpoint) z
                  })
                  answer <- tryCatch(native(probe, id, bid, e, pr, d, base), error = function(err) err)
                  if (inherits(answer, "error") && identical(conditionMessage(answer), "FF_A_SELECTION_WEIGHT_HASH_MISMATCH A_NONREUSABLE_COEFFICIENT_ONLY_WEIGHT_VECTOR_NOT_RETAINED")) {
                    ff_a_assert(!isTRUE(ds$individual_weight_vector_retained) && is.list(params) && length(params$selection_coefficients) > 
                      0L, "ADD9_ROUNDOFF_REQUIRES_SEALED_SELECTION_COEFFICIENTS")
                    numeric_ds <- ds
                    numeric_ds$selection_weight_hash <- NULL
                    probe$a_source_for <- local({
                      z <- numeric_ds
                      function(id, bid, endpoint) z
                    })
                    answer <- native(probe, id, bid, e, pr, d, base)
                    numeric_check <- b16_a_weight_output_check(ctx, answer, params, list(weight_proof_kind = "RECONSTRUCTED_SEALED_SELECTION_COEFFICIENTS"), 
                      ds$selection_weight_hash, ff_a_hash(format(answer$mi$w, digits = 17, scientific = TRUE)), 
                      d, ca$column)
                    check$numeric_reconstruction_check <- numeric_check
                    check$historical_weight_identity_claimed <- FALSE
                  }
                  if (inherits(answer, "error")) {
                    reason <- b9_a_rederivation_reason(conditionMessage(answer))
                    ff_a_assert(!is.null(reason), conditionMessage(answer))
                    check$status <- "NOT_CERTIFIED_FOR_TASK"
                    check$reason <- reason
                  } else {
                    ff_a_assert(!is.null(params) && is.list(params), "ADD9_VALID_A_PARAMETER_OBJECT_REQUIRED")
                    check$status <- "CERTIFIED_CURRENT_TASK_DONOR"
                    check$A_parameter_hash <- b16_a_parameter_hash(params)
                    if (length(donors)) check$certified_donor_numeric_outputs_match <- b9_a_compare_donors(donors[[1L]]$answer, 
                      answer)
                    donors[[length(donors) + 1L]] <- list(params = params, answer = if (!length(donors)) answer else NULL, 
                      descriptor = ds)
                  }
                  audits[[i]] <- check
                }
                reason <- if (length(donors)) "REUSED_CERTIFIED_CURRENT_TASK_DONOR" else "AUDITED_MISSING_OR_NOT_CERTIFIED_FOR_CURRENT_TASK_UNCHANGED_FROZEN_DERIVATION"
                dr <- get0(as.character(bid), s$draws, inherits = FALSE)
                ff_a_assert(is.list(dr) && identical(dr$draw_hash, ctx$actual_draw_hash) && b16_a_cache_is_hash(dr$actual_draw_private_sha256), 
                  "ADD9_ACTUAL_RETAINED_DRAW_REQUIRED")
                sourceproof <- list(schema = "ADD9_A_TASK_SOURCE_AUDIT_V1", models_run = 0L, a_group_id = gid, 
                  indicator_id = id, bootstrap_id = bid, endpoint_ids = as.list(endpoints), binding = pre, source_policy_sha256 = s$policy_sha256, 
                  source_policy_complete_sha256 = s$policy_complete_sha256, source_disposition = reason, existing_donor_count = length(donors), 
                  source_checks = audits, catalog_checks = catalog$catalog_checks, upstream_policy = s$policy$upstream_compatibility[[id]], 
                  actual_draw_receipt = dr, global_storage_absence_claim = FALSE)
                sp <- file.path(group, "SOURCE_PROOF.json")
                b9_a_put(sourceproof, sp)
                selection_fits <- A_fits <- 0L
                if (length(donors)) {
                  params <- donors[[1L]]$params
                  reference <- donors[[1L]]$answer
                  mi <- reference$mi
                } else {
                  ff_a_assert(identical(Sys.getenv("FF_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "ADD9_A_DERIVATION_REQUIRES_EXPLICIT_MANUAL_START")
                  rng <- get0(as.character(bid), s$draw_rng, inherits = FALSE)
                  ff_a_assert(is.integer(rng) && length(rng) > 0L, "ADD9_FROZEN_POST_DRAW_RNG_STATE_REQUIRED")
                  assign(".Random.seed", rng, envir = globalenv())
                  ff_w_stage(ctx, "A_FROZEN_DERIVATION", ep)
                  mi <- ctx$engine$b3_missing(d, ca$column, base)
                  selection_fits <- as.integer(!all(mi$complete))
                  ff_a_assert(identical(mi$complete, complete.cases(d[, unique(c(base, ca$column)), drop = FALSE])), 
                    "ADD9_FROZEN_SELECTION_SAMPLE_CHANGED")
                  fit <- ctx$engine$fp_a_path(d, ca$column, base, mi, ctx, paste0("single_", ca$candidate_id), 
                    as.integer(bid))
                  A_fits <- 1L
                  dd <- d[mi$complete, , drop = FALSE]
                  X <- model.matrix(reformulate(base), dd)
                  fd <- d[, unique(c(base, ca$column)), drop = FALSE]
                  params <- list(alpha = fit$alpha, mediator = NULL, mediators = ca$column, selection_coefficients = if (is.null(mi$fit)) NULL else coef(mi$fit), 
                    selection_formula = deparse(ff_a_selection_formula()), base_terms = base, factor_levels = lapply(fd[, 
                      vapply(fd, is.factor, FALSE), drop = FALSE], levels), columns = colnames(X), spline_reference = ctx$engine$b3_design(dd, 
                      ca$column, base)$reference)
                  physical <- ctx$engine$fp_physical(d, ca$column, mi, fit$alpha, ctx)
                  ff_a_assert(all(physical$status == "VALID"), "ADD9_FROZEN_A_PHYSICAL_DIAGNOSTIC_FAILED")
                  reference <- list(mi = mi, parent_rows = pr[mi$complete], M0 = dd[[ca$column]] - dd$B * params$alpha[[ca$column]])
                  reference$M1 <- reference$M0 + params$alpha[[ca$column]]
                }
                binding <- b16_a_cache_binding(ctx, id, bid, pr, d, base, params)
                weight_hash <- ff_a_hash(format(mi$w, digits = 17, scientific = TRUE))
                rp <- file.path(group, "PARAMETERS.private.rds")
                envelope <- list(schema = "DERIVED_A_GROUP_RDS_V1", a_group_id = gid, indicator_id = id, bootstrap_id = as.integer(bid), 
                  binding = binding, parameters = params, selection_weights = mi$w)
                ff_a_assert(!file.exists(rp), "ADD9_PRESERVE_EXISTING_A_VECTOR")
                saveRDS(envelope, rp, compress = FALSE, version = 2)
                proof <- list(schema = "DERIVED_A_GROUP_RECEIPT_V1", complete = TRUE, ready = FALSE, a_group_id = gid, 
                  indicator_id = id, bootstrap_id = as.integer(bid), endpoint_ids = as.list(endpoints), binding = binding, 
                  source_proof_path = b9_a_canonical(sp), source_proof_sha256 = ff_a_sha(sp), source_disposition = reason, 
                  derivation_transport_sha256 = s$transport_sha256, parameters_sha256 = ff_a_sha(rp), stored_selection_weight_hash = weight_hash, 
                  historical_selection_weight_hash = NULL, historical_weight_identity_claimed = FALSE, weight_origin = reason, 
                  weight_normalization = "NONE", weight_trimming = "NONE", weight_stabilization = "NONE", A_fits = A_fits, 
                  selection_fits = selection_fits, Y_fits = 0L, N_target = nrow(d), N_complete = length(mi$w), 
                  host = ff_p_host(), task_cache = list(schema = "B16_A_TASK_CACHE_BINDING_V1", certified_for_task_cache = TRUE, 
                    a_group_id = gid, endpoint_ids = as.list(endpoints), binding = binding, weight_proof_kind = "RETAINED_INDIVIDUAL_VECTOR", 
                    weight_proof_sha256 = b16_a_cache_digest(list(weight_hash = weight_hash, parameters_sha256 = ff_a_sha(rp))), 
                    reconstructed_weight_hash = weight_hash, historical_selection_weight_hash = NULL))
                pp <- file.path(group, "PROOF.json")
                b9_a_put(proof, pp)
                ff_w_stage(ctx, "A_RETAINED_CACHE_VALIDATION", ep)
                desc <- list(format = "derived_group_rds", path = b9_a_canonical(rp), sha256 = ff_a_sha(rp), a_group_id = gid, 
                  a_group_proof_path = b9_a_canonical(pp), proof_receipt_sha256 = ff_a_sha(pp), derivation_transport_sha256 = s$transport_sha256, 
                  derived_selection_weight_hash = weight_hash, indicator_id = id, bootstrap_id = as.integer(bid), 
                  draw_hash = binding$actual_draw_hash, target_hash = binding$target_hash, sample_hash = binding$sample_hash, 
                  source_input_lock_sha256 = binding$source_input_lock_sha256, source_prepared_sha256 = binding$source_prepared_sha256, 
                  missing_object = FALSE, selection_object_missing = FALSE, individual_weight_vector_retained = TRUE, 
                  selection_weight_hash = NULL, reason = reason)
                probe <- ctx
                probe$a_source_for <- local({
                  z <- desc
                  function(id, bid, endpoint) {
                    z$endpoint_id <- endpoint
                    z
                  }
                })
                replay <- native(probe, id, bid, ep, pr, d, base)
                ff_a_assert(identical(replay$mi$w, mi$w) && identical(replay$parent_rows, reference$parent_rows), 
                  "ADD9_RETAINED_A_WEIGHT_OR_SAMPLE_ROUNDTRIP_CHANGED")
                difference <- max(abs(replay$M0 - reference$M0), abs(replay$M1 - reference$M1))
                ff_a_assert(is.finite(difference) && difference < 1e-12, "ADD9_RETAINED_A_SCENARIO_ROUNDTRIP_CHANGED")
                rows <- lapply(endpoints, function(e) {
                  z <- desc
                  z$endpoint_id <- e
                  z
                })
                dp <- file.path(group, "DESCRIPTORS.json")
                b9_a_put(list(rows = rows), dp)
                b9_a_put(list(schema = "ADD9_A_GROUP_COMPLETE_V1", complete = TRUE, a_group_id = gid, source_policy_sha256 = s$policy_sha256, 
                  descriptors_sha256 = ff_a_sha(dp), A_fits = A_fits, selection_fits = selection_fits, Y_fits = 0L, 
                  source_disposition = reason, roundtrip_max_abs_difference = difference, actual_weight_bitwise_roundtrip = TRUE), 
                  completepath)
                for (z in rows) assign(paste(id, bid, z$endpoint_id, sep = "|"), z, s$rows)
                invisible(TRUE)
            }, error = function(error) {
                fp <- file.path(group, "FAILED.json")
                if (!file.exists(fp)) b9_a_put(list(a_group_id = gid, indicator_id = id, bootstrap_id = bid, reason = conditionMessage(error), 
                  outputs_retained = TRUE, automatic_recompute_authorized = FALSE), fp)
                stop(error)
            })
        }
        b9_a_install <- function(envir = parent.frame()) {
            cached <- get("ff_a_scenarios", envir = envir)
            native <- attr(cached, "b16_native_a")
            ff_a_assert(is.function(native), "ADD9_INSTALL_REQUIRES_FROZEN_F52_CACHE")
            wrapper <- function(ctx, id, bid, ep, pr, d, base) {
                b9_a_materialize(ctx, id, as.integer(bid), ep, pr, d, base)
                ff_w_stage(ctx, "A_RETAINED_CACHE_VALIDATION", ep)
                cached(ctx, id, bid, ep, pr, d, base)
            }
            attr(wrapper, "b16_native_a") <- native
            assign("ff_a_scenarios", wrapper, envir = envir)
            invisible(TRUE)
        }
    }), envir = context)
}

