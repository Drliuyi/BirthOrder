# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

c39_assert <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)

c39_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

c39_sha <- function(path) {
    c39_assert(file.exists(path) && !dir.exists(path), paste("CORE39_REQUIRED_FILE", path))
    digest::digest(file = path, algo = "sha256")
}

c39_read <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)

c39_json <- function(x, path) {
    c39_assert(!file.exists(path), paste("CORE39_PRESERVE_EXISTING", path))
    tmp <- paste0(path, ".pending_", Sys.getpid())
    c39_assert(!file.exists(tmp), "CORE39_PENDING_EXISTS")
    jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
    c39_assert(!file.exists(path) && file.rename(tmp, path), "CORE39_ATOMIC_JSON_FAILED")
}

c39_rds <- function(x, path) {
    c39_assert(!file.exists(path), paste("CORE39_PRESERVE_EXISTING", path))
    tmp <- paste0(path, ".pending_", Sys.getpid())
    c39_assert(!file.exists(tmp), "CORE39_PENDING_EXISTS")
    saveRDS(x, tmp, compress = FALSE, version = 2L)
    c39_assert(!file.exists(path) && file.rename(tmp, path), "CORE39_ATOMIC_RDS_FAILED")
}

c39_run <- function(run, create = FALSE) {
    c39_assert(!grepl("(^|/)[.][.](/|$)", run) && grepl("/UKB_results/birth/border_jointdomain_decomposition_20261101_v2/[^/]*CORE39[^/]*$", 
        run) && (startsWith(run, birth_private_setting("PRIVATE_INPUT")) || startsWith(run, birth_private_setting("PRIVATE_INPUT"))), 
        "CORE39_ISOLATED_SHARED_RUN_REQUIRED")
    if (create && !dir.exists(run)) 
        c39_assert(dir.create(run, recursive = TRUE), "CORE39_RUN_CREATE_FAILED")
    normalizePath(run, mustWork = TRUE)
}

c39_config <- function(root) {
    p <- file.path(root, "config/CORE39_DESIGN.json")
    cfg <- c39_read(p)
    cfg$domains <- lapply(cfg$domains, function(x) as.character(unlist(x, use.names = FALSE)))
    cfg$generation_order <- unlist(cfg$domains, use.names = FALSE)
    c39_assert(length(cfg$domains) == 15L && length(cfg$generation_order) == 39L && !anyDuplicated(cfg$generation_order) && 
        isTRUE(cfg$user_confirmation_received), "CORE39_EXACT_15_39_REQUIRED")
    c39_assert(identical(as.integer(unlist(cfg$bootstrap$ids)), c(0L, 100L)), "CORE39_B0_PLUS_FIXED_B100_ONLY")
    cfg$config_hash <- c39_sha(p)
    cfg
}

c39_rng <- function(code) {
    kind <- RNGkind()
    had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (had) 
        old <- get(".Random.seed", .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(kind))
        if (had) assign(".Random.seed", old, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", 
            envir = .GlobalEnv)
    })
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    force(code)
}

c39_selection_engine <- function(text) {
    e <- new.env(parent = globalenv())
    for (n in names(text)) e[[n]] <- eval(parse(text = text[[n]]), envir = e)
    e
}

c39_unwrap_monitor <- function(f) {
    for (i in 1:4) {
        body_text <- paste(deparse(body(f)), collapse = " ")
        if (!grepl("v5_frame_endpoint", body_text, fixed = TRUE)) 
            break
        e <- environment(f)
        c39_assert(exists("fn", e, inherits = FALSE) && is.function(e$fn), "CORE39_UNKNOWN_SELECTION_WRAPPER")
        f <- e$fn
    }
    c39_assert(!grepl("v5_frame_endpoint", paste(deparse(body(f)), collapse = " "), fixed = TRUE), "CORE39_WRAPPER_NOT_REMOVED")
    f
}

c39_members <- function(ctx, cfg) {
    ans <- lapply(cfg$generation_order, function(id) {
        ca <- ctx$features[[id]]
        s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
        c39_assert(is.list(ca) && nrow(s) == 1L && identical(ctx$variables[[id]]$status, "AVAILABLE_FROZEN"), paste("CORE39_MEMBER_NOT_FROZEN", 
            id))
        type <- as.character(s$type)
        x <- ctx$d_full[[ca$column]]
        c39_assert(type %in% c("continuous", "binary", "categorical", "ordinal", "nominal") && !is.null(x), paste("CORE39_MEMBER_TYPE_OR_COLUMN", 
            id))
        if (type == "continuous") {
            allowed <- identical(ca$prepared_transform, "identity") || (identical(id, "clinical__tg") && identical(ca$prepared_transform, 
                "log") && identical(ca$prepared_unit, "log(mmol/L)") && identical(ca$inverse_transform, "exp") && 
                is.null(ca$physical_lower_bound) && is.null(ca$physical_upper_bound))
            c39_assert(allowed, paste("CORE39_UNREVIEWED_TRANSFORM", id))
        }
        list(indicator_id = id, column = ca$column, type = type, levels = if (type != "continuous") levels(x) else NULL, 
            lower = if (is.null(ca$physical_lower_bound)) -Inf else ca$physical_lower_bound, upper = if (is.null(ca$physical_upper_bound)) Inf else ca$physical_upper_bound, 
            prepared_transform = ca$prepared_transform, prepared_unit = ca$prepared_unit, raw_unit = ca$raw_unit, 
            inverse_transform = ca$inverse_transform, technical_covariates = ca$technical_covariates)
    })
    names(ans) <- cfg$generation_order
    c39_assert(!anyDuplicated(vapply(ans, `[[`, "", "column")), "CORE39_DUPLICATE_COORDINATE")
    ans
}

c39_selection <- function(engine, d, columns, base) {
    z <- engine$b3_missing(d, columns, base)
    c39_assert(identical(z$complete, complete.cases(d[, unique(c(base, columns)), drop = FALSE])), "CORE39_SELECTION_MASK_DRIFT")
    c39_assert(length(z$w) == sum(z$complete) && all(is.finite(z$w) & z$w > 0), "CORE39_SELECTION_WEIGHT_INVALID")
    z$design_weights <- rep(1, length(z$w))
    z$selection_weights <- as.numeric(z$w)
    z$selection_coefficients <- if (is.null(z$fit)) 
        NULL
    else stats::coef(z$fit)
    z$selection_formula <- if (is.null(z$fit)) 
        NULL
    else paste(deparse(stats::formula(z$fit)), collapse = " ")
    z$fit <- NULL
    z
}

c39_counts <- function(d, risks, complete) {
    list(N_target = nrow(d), N_joint_fit = sum(complete), N_standardization = sum(complete), N_joint39_missing = sum(!complete), 
        joint39_missing_fraction = mean(!complete), N_families = length(unique(d$cluster_id)), endpoints = lapply(risks, 
            function(r) list(target_events = sum(r$event), target_competing_deaths = sum(r$death_event), joint_fit_events = sum(r$event[complete]), 
                joint_fit_competing_deaths = sum(r$death_event[complete]), followup_max = max(r$exit), followup_not_capped_at_11 = TRUE)))
}

c39_member_audit <- function(d, members, complete) {
    lapply(members, function(m) {
        x <- d[[m$column]]
        y <- x[complete]
        list(indicator_id = m$indicator_id, column = m$column, type = m$type, lower = m$lower, upper = m$upper, 
            prepared_transform = m$prepared_transform, missing_target = sum(is.na(x)), observed_target = sum(!is.na(x)), 
            zero_on_prepared_scale_target = sum(x == 0, na.rm = TRUE), true_raw_zero_target = if (identical(m$prepared_transform, 
                "identity")) sum(x == 0, na.rm = TRUE) else NULL, half_value_target = sum(x == 0.5, na.rm = TRUE), 
            negative_target = if (is.numeric(x)) sum(x < 0, na.rm = TRUE) else NULL, complete_support = if (is.factor(y)) as.list(table(y)) else as.list(stats::quantile(y, 
                c(0, 0.01, 0.5, 0.99, 1))), observed_range_is_not_physical_bound = TRUE)
    })
}

c39_attach_backend <- function(ctx, engine) {
    for (n in ls(ctx$v4, all.names = TRUE)) if (startsWith(n, "ff_") && is.function(ctx$v4[[n]])) 
        engine[[n]] <- ctx$v4[[n]]
    invisible(engine)
}

c39_backend_load <- function(root, run, engine, expected) {
    v4root <- Sys.getenv("FD_V4_ROOT", file.path(dirname(root), "border_fullfollowup_b1000_add9_20261030_v1"))
    decision <- c39_read(file.path(v4root, "config/V4_METHOD_DECISION.json"))
    for (n in names(decision$scientific_files_sha256)) c39_assert(identical(c39_sha(file.path(v4root, n)), decision$scientific_files_sha256[[n]]), 
        paste("CORE39_V4_SOURCE_CHANGED", n))
    c39_assert(identical(decision$method_spec_sha256, expected), "CORE39_METHOD_CHANGED")
    env <- new.env(parent = globalenv())
    for (n in c("10.prepare.R", "20.a_reuse.R", "30.outcome.R", "31.aj_backend.R", "61.center_strata.R")) sys.source(file.path(v4root, 
        "f", n), envir = env)
    c39_assert(identical(env$ff_ce_method_spec()$sha256, expected), "CORE39_V4_METHOD_CODE_CHANGED")
    cfg <- c39_read(file.path(v4root, "config/runtime.json"))
    host <- Sys.getenv("FF_HOST", "winpc")
    work <- Sys.getenv("FD_V4_WORK", cfg$work_paths[[host]])
    runtime <- list(R = R.version$version.string, platform = R.version$platform, arch = R.version$arch, pointer_size = .Machine$sizeof.pointer, 
        extension = .Platform$dynlib.ext, isolated_build_contract = "R_CMD_SHLIB_V2_EXPLICIT_RUNTIME_COMPILER_PATH", 
        compiler_paths = Sys.which(c("cc", "clang", "gcc", "arm64-apple-darwin20.0.0-clang")))
    build <- file.path(work, "software", paste0("aj_", substr(c39_sha(file.path(v4root, "src/aj.c")), 1, 16), "_", 
        substr(digest::digest(runtime, algo = "sha256"), 1, 12)))
    c39_assert(file.exists(file.path(build, "BUILD_MANIFEST.rds")) && file.exists(file.path(build, paste0("ff_aj", 
        .Platform$dynlib.ext))), "CORE39_EXISTING_AJ_BINARY_REQUIRED")
    backend <- env$ff_aj_initialize(v4root, work, require_compiled = TRUE)
    z <- list(v4 = env, v4_decision = decision, backend = backend)
    c39_attach_backend(z, engine)
    z
}

c39_prepare <- function(root, run, engine) {
    c39_assert(Sys.info()[["sysname"]] != "Darwin" && identical(Sys.getenv("FF_HOST"), "winpc"), "CORE39_PREPARATION_WINPC_ONLY")
    root <- normalizePath(root, mustWork = TRUE)
    run <- c39_run(run, TRUE)
    cfg <- c39_config(root)
    private <- file.path(run, "private")
    if (!dir.exists(private)) 
        dir.create(private)
    lockpath <- file.path(private, "INPUT_MANIFEST.json")
    c39_assert(!file.exists(lockpath), "CORE39_INPUT_ALREADY_PREPARED")
    old <- file.path(dirname(root), "border_fine_domain_behavior_joint_20261031_v1")
    loaderpath <- file.path(old, "f/10.inputs.R")
    loaderhash <- c39_read(file.path(root, "config/PROPOSED_DESIGN.json"))$source_files$old_loader$sha256
    c39_assert(identical(c39_sha(loaderpath), loaderhash), "CORE39_ORIGINAL_LOADER_CHANGED")
    loader <- new.env(parent = engine)
    sys.source(loaderpath, envir = loader)
    ctx <- loader$fd_load(old, list(), run)
    members <- c39_members(ctx, cfg)
    base <- unique(c(as.character(unlist(ctx$cfg$base_covariates)), unlist(lapply(members, `[[`, "technical_covariates"))))
    c39_assert(identical(base, as.character(unlist(cfg$base))), "CORE39_ORIGINAL_BASE_CHANGED")
    ep <- ctx$endpoints[c("250.2", "411")]
    c39_assert(all(vapply(ep, function(x) identical(x$status, "AVAILABLE_FROZEN"), FALSE)), "CORE39_ENDPOINT_UNAVAILABLE")
    shared <- ep[[1]]$eligibility & ep[[2]]$eligibility
    c39_assert(is.logical(shared) && !anyNA(shared) && length(shared) == nrow(ctx$d_full), "CORE39_INVALID_TARGET_MASK")
    for (e in names(ep)) c39_assert(is.logical(ep[[e]]$incident_eligibility) && length(ep[[e]]$incident_eligibility) == 
        length(shared) && !anyNA(ep[[e]]$incident_eligibility) && all(ep[[e]]$incident_eligibility[shared]), paste("CORE39_ORIGINAL_INCIDENT_RISKSET_NOT_COVERING_TARGET", 
        e))
    for (id in cfg$generation_order) c39_assert(all(ctx$eligibility[[id]][shared]), paste("CORE39_MEMBER_APPLICABILITY_TARGET_CONFLICT", 
        id))
    rows <- which(shared)
    c39_assert(length(rows) > 0L, "CORE39_EMPTY_TARGET")
    columns <- vapply(members, `[[`, "", "column")
    keep <- unique(c("eid", "cluster_id", base, columns, grep("^aux_", names(ctx$d_full), value = TRUE)))
    d <- ctx$d_full[rows, keep, drop = FALSE]
    rownames(d) <- NULL
    c39_assert(all(as.numeric(as.character(d$S)) %in% 2:8) && all(is.finite(d$B) & d$B >= 0 & d$B <= 1) && !anyNA(d$cluster_id) && 
        !anyDuplicated(d$eid), "CORE39_ORIGINAL_ELIGIBILITY_INVALID")
    risks <- lapply(ep, function(x) x$risk[rows, , drop = FALSE])
    for (r in risks) c39_assert(identical(as.character(r$eid), as.character(d$eid)) && all(is.finite(r$exit) & 
        r$exit > 0 & r$entry == 0 & r$event %in% 0:1 & r$death_event %in% 0:1 & r$event + r$death_event <= 1), 
        "CORE39_RISK_ALIGNMENT_INVALID")
    funcs <- c("b3_missing", "b3_formula", "b3_capture", "b3_assert", "b3_boot_index")
    texts <- setNames(lapply(funcs, function(n) paste(deparse(c39_unwrap_monitor(ctx$engine[[n]]), width.cutoff = 500L), 
        collapse = "\n")), funcs)
    selection_engine <- c39_selection_engine(texts)
    selection <- c39_selection(selection_engine, d, columns, base)
    counts <- c39_counts(d, risks, selection$complete)
    counts$N_missing_base_covariates <- sum(!complete.cases(d[, base, drop = FALSE]))
    engine$fd_support(d[selection$complete, , drop = FALSE], members)
    domains <- lapply(cfg$domains, function(ids) vapply(members[ids], `[[`, "", "column"))
    model <- list(members = members, base = base, continuous_rule = cfg$continuous_rule, count_rule = cfg$count_rule)
    object <- list(schema = "CORE39_FROZEN_ACTUAL_INPUTS_V21", data = d, risks = risks, members = members, model = model, 
        domains = domains, parent_rows = rows, parent_frame_hash = ctx$parent_frame_hash, source_input_hash = ctx$input_hash, 
        source_hash = ctx$source_hash, method_hash = ctx$v4_decision$method_spec_sha256, config_hash = cfg$config_hash, 
        target_hash = c39_hash(as.character(d$eid)), cluster_membership_hash = c39_hash(paste(d$eid, d$cluster_id, 
            sep = ":")), selection_functions = texts, selection_source_hash = c39_hash(texts), B0_selection = selection, 
        counts = counts)
    inputpath <- file.path(private, "CORE39_INPUTS.private.rds")
    c39_rds(object, inputpath)
    drawsdir <- file.path(private, "draws")
    c39_assert(!dir.exists(drawsdir) && dir.create(drawsdir), "CORE39_DRAW_DIRECTORY_EXISTS")
    clusters <- unique(as.character(d$cluster_id))
    ci <- match(as.character(d$cluster_id), clusters)
    sizes <- tabulate(ci, length(clusters))
    draws <- c39_rng(lapply(0:100, function(b) {
        ix <- selection_engine$b3_boot_index(d, as.integer(b), as.integer(cfg$bootstrap$seed))
        multiplicity <- tabulate(ci[ix], length(clusters))/sizes
        c39_assert(all(is.finite(multiplicity) & multiplicity >= 0 & multiplicity == floor(multiplicity)) && sum(multiplicity) == 
            length(clusters), "CORE39_WHOLE_CLUSTER_DRAW_INVALID")
        p <- file.path(drawsdir, sprintf("bootstrap_%04d.private.rds", b))
        draw <- list(bootstrap_id = as.integer(b), indices = as.integer(ix), cluster_multiplicity = as.integer(multiplicity), 
            target_hash = object$target_hash, cluster_membership_hash = object$cluster_membership_hash, seed = as.integer(cfg$bootstrap$seed) + 
                1009L * as.integer(b), RNG_kind = RNGkind())
        c39_rds(draw, p)
        if (b%%100L == 0L) {
            message("CORE39_ACTUAL_DRAW_FROZEN ", b, "/100")
            flush.console()
        }
        list(b = as.integer(b), seed = draw$seed, sampling_hash = c39_hash(draw$indices), index_path = p, index_sha256 = c39_sha(p), 
            index_relative_path = file.path("private/draws", basename(p)), N_resampled = length(ix), multiplicity_hash = c39_hash(draw$cluster_multiplicity), 
            target_hash = object$target_hash)
    }))
    drawmanifest <- file.path(private, "DRAW_MANIFEST.json")
    c39_json(list(schema = "CORE39_ACTUAL_CLUSTER_DRAWS_V21", unit = "Frozen cluster_id family", mechanism = cfg$bootstrap$mechanism, 
        input_sha256 = c39_sha(inputpath), target_hash = object$target_hash, cluster_membership_hash = object$cluster_membership_hash, 
        RNG_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"), fixed_attempts = 100L, draws = draws), drawmanifest)
    manifest <- list(schema = object$schema, input_path = inputpath, input_relative_path = "private/CORE39_INPUTS.private.rds", 
        input_hash = c39_sha(inputpath), files = list(list(path = inputpath, sha256 = c39_sha(inputpath))), config_hash = cfg$config_hash, 
        source_input_hash = ctx$input_hash, source_hash = ctx$source_hash, method_hash = object$method_hash, selection_source_hash = object$selection_source_hash, 
        target_hash = object$target_hash, cluster_membership_hash = object$cluster_membership_hash, draw_manifest_path = drawmanifest, 
        draw_manifest_hash = c39_sha(drawmanifest), design_weight_rule = cfg$design_weights, selection_rule = cfg$selection, 
        standardization = cfg$standardization, frozen_scale_notes = cfg$frozen_scale_notes, missingness_assumptions = cfg$missingness_assumptions, 
        counts = counts, selection = selection[c("diagnostics", "distribution", "warnings", "selection_coefficients", 
            "selection_formula")], member_audit = c39_member_audit(d, members, selection$complete), frozen_feature_metadata = ctx$features[c(cfg$generation_order, 
            "whr", "blood_pressure__MAP", "weight")], frozen_variable_metadata = lapply(ctx$variables[cfg$generation_order], 
            function(x) x[intersect(names(x), c("indicator_id", "status", "column", "field", "instance", "encoding", 
                "eligible", "levels", "reference", "ordered", "source_definition", "target_rule", "invalid_missing_codes", 
                "measurement_instance", "source_label_retained"))]), endpoint_censoring = "Original frozen per-endpoint risk, no cross-disease censoring added", 
        A_Y_models_fitted = 0L)
    c39_json(manifest, lockpath)
    c39_attach_backend(ctx, engine)
    list(input_manifest = lockpath, draw_manifest = drawmanifest, input_hash = manifest$input_hash, counts = counts, 
        config_hash = cfg$config_hash, source_hash = ctx$source_hash, method_hash = object$method_hash)
}

c39_load <- function(root, run, engine) {
    c39_assert(identical(Sys.getenv("FD_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "CORE39_REAL_AUTHORITY_REQUIRED")
    if (Sys.info()[["sysname"]] == "Darwin") 
        c39_assert(identical(Sys.getenv("FD_ALLOW_REVIEWED_MAC_COMPUTE"), "YES") && identical(Sys.getenv("FF_HOST"), 
            "macmini"), "CORE39_NO_PARTICIPANT_DATA_ON_CONTROL_MAC")
    run <- c39_run(run)
    cfg <- c39_config(root)
    manifest <- c39_read(file.path(run, "private/INPUT_MANIFEST.json"))
    c39_assert(identical(manifest$config_hash, cfg$config_hash), "CORE39_CONFIG_CHANGED")
    p <- file.path(run, manifest$input_relative_path)
    c39_assert(identical(c39_sha(p), manifest$input_hash), "CORE39_FROZEN_INPUT_CHANGED")
    ctx <- readRDS(p)
    c39_assert(identical(ctx$config_hash, cfg$config_hash) && identical(c39_hash(ctx$selection_functions), manifest$selection_source_hash), 
        "CORE39_INPUT_BINDING_CHANGED")
    dp <- file.path(run, "private/DRAW_MANIFEST.json")
    c39_assert(identical(c39_sha(dp), manifest$draw_manifest_hash), "CORE39_DRAW_MANIFEST_CHANGED")
    dm <- c39_read(dp)
    c39_assert(identical(dm$input_sha256, manifest$input_hash) && length(dm$draws) == 101L && dm$fixed_attempts == 
        100L && identical(vapply(dm$draws, function(x) as.integer(x$b), 0L), 0:100), "CORE39_DRAW_INPUT_MISMATCH")
    ctx$selection_engine <- c39_selection_engine(ctx$selection_functions)
    ctx$run <- run
    ctx$cfg <- cfg
    ctx$input_hash <- manifest$input_hash
    ctx$manifest <- manifest
    ctx$draw_manifest <- dm
    ctx$engine <- engine
    backend <- c39_backend_load(root, run, engine, ctx$method_hash)
    ctx$v4 <- backend$v4
    ctx$backend <- backend$backend
    ctx
}

c39_task <- function(ctx, b) {
    c39_assert(length(b) == 1L && is.finite(b) && b == as.integer(b) && b %in% 0:100, "CORE39_FIXED_ID_REQUIRED")
    rec <- ctx$draw_manifest$draws[[as.integer(b) + 1L]]
    c39_assert(rec$b == b, "CORE39_DRAW_ID_MISMATCH")
    p <- file.path(ctx$run, rec$index_relative_path)
    c39_assert(identical(c39_sha(p), rec$index_sha256), "CORE39_ACTUAL_INDEX_CHANGED")
    actual <- readRDS(p)
    ix <- actual$indices
    c39_assert(actual$bootstrap_id == b && identical(actual$target_hash, ctx$target_hash) && identical(actual$cluster_membership_hash, 
        ctx$cluster_membership_hash) && identical(c39_hash(ix), rec$sampling_hash) && identical(c39_hash(actual$cluster_multiplicity), 
        rec$multiplicity_hash) && all(is.finite(ix) & ix >= 1L & ix <= nrow(ctx$data)), "CORE39_ACTUAL_DRAW_BINDING_INVALID")
    d <- ctx$data[ix, , drop = FALSE]
    cols <- vapply(ctx$members, `[[`, "", "column")
    mi <- if (b == 0L) 
        ctx$B0_selection
    else c39_selection(ctx$selection_engine, d, cols, ctx$model$base)
    c39_assert(length(mi$complete) == nrow(d), "CORE39_DRAW_SELECTION_LENGTH_MISMATCH")
    complete <- mi$complete
    fit <- d[complete, , drop = FALSE]
    risks <- lapply(ctx$risks, function(r) r[ix[complete], , drop = FALSE])
    w <- as.numeric(mi$design_weights * mi$selection_weights)
    ctx$engine$fd_support(fit, ctx$members)
    support <- local({
        members <- ctx$members
        eng <- ctx$engine
        function(m, C = NULL, b = NULL) {
            eng$fd_support(m, members)
            dbp <- members[["blood_pressure__DBP"]]$column
            sbp <- members[["blood_pressure__SBP"]]$column
            c39_assert(all(m[[sbp]] > m[[dbp]]), "CORE39_SBP_NOT_GREATER_THAN_DBP")
            TRUE
        }
    })
    fit <- fit[, unique(c(ctx$model$base, cols)), drop = FALSE]
    prepared <- list(fit_data = fit, risks = risks, fit_weights = w, target = fit, standardization_weights = w, 
        design_weights = mi$design_weights, selection_weights = mi$selection_weights, cluster = as.character(d$cluster_id[complete]), 
        joint_target_identity_checked = TRUE, selection_weight_process_checked = TRUE, validate_joint_support = support, 
        actual_draw_proof = c(rec, list(actual_index_sha256 = c39_sha(p), index_path_host = p, new_common_target = TRUE, 
            old_indices_reused = FALSE)), selection = mi, counts = c39_counts(d, lapply(ctx$risks, function(r) r[ix, 
            , drop = FALSE]), complete))
    binding <- list(input_hash = ctx$input_hash, config_hash = ctx$config_hash, method_hash = ctx$method_hash, 
        bootstrap_id = as.integer(b), sampling_hash = rec$sampling_hash, actual_index_sha256 = rec$index_sha256, 
        target_hash = c39_hash(as.character(d$eid)), sample_hash = c39_hash(as.character(d$eid[complete])), design_weight_hash = c39_hash(mi$design_weights), 
        selection_weight_hash = c39_hash(mi$selection_weights), weight_hash = c39_hash(w), generation_order = cols, 
        domain_partition_sha256 = c39_hash(ctx$domains), allocation_order = names(ctx$domains), mc_draws = 50L, 
        mc_seed = as.integer(20261101L + 1009L * b + 100003L))
    list(prepared = prepared, model = ctx$model, domains = ctx$domains, binding = binding, members = ctx$members, 
        base = ctx$model$base, counts = prepared$counts, actual_draw_proof = prepared$actual_draw_proof)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        c39_assert <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
        c39_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        c39_sha <- function(path) {
            c39_assert(file.exists(path) && !dir.exists(path), paste("CORE39_REQUIRED_FILE", path))
            digest::digest(file = path, algo = "sha256")
        }
        c39_read <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
        c39_json <- function(x, path) {
            c39_assert(!file.exists(path), paste("CORE39_PRESERVE_EXISTING", path))
            tmp <- paste0(path, ".pending_", Sys.getpid())
            c39_assert(!file.exists(tmp), "CORE39_PENDING_EXISTS")
            jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
            c39_assert(!file.exists(path) && file.rename(tmp, path), "CORE39_ATOMIC_JSON_FAILED")
        }
        c39_rds <- function(x, path) {
            c39_assert(!file.exists(path), paste("CORE39_PRESERVE_EXISTING", path))
            tmp <- paste0(path, ".pending_", Sys.getpid())
            c39_assert(!file.exists(tmp), "CORE39_PENDING_EXISTS")
            saveRDS(x, tmp, compress = FALSE, version = 2L)
            c39_assert(!file.exists(path) && file.rename(tmp, path), "CORE39_ATOMIC_RDS_FAILED")
        }
        c39_run <- function(run, create = FALSE) {
            c39_assert(!grepl("(^|/)[.][.](/|$)", run) && grepl("/UKB_results/birth/border_jointdomain_decomposition_20261101_v2/[^/]*CORE39[^/]*$", 
                run) && (startsWith(run, birth_private_setting("PRIVATE_INPUT")) || startsWith(run, birth_private_setting("PRIVATE_INPUT"))), 
                "CORE39_ISOLATED_SHARED_RUN_REQUIRED")
            if (create && !dir.exists(run)) c39_assert(dir.create(run, recursive = TRUE), "CORE39_RUN_CREATE_FAILED")
            normalizePath(run, mustWork = TRUE)
        }
        c39_config <- function(root) {
            p <- file.path(root, "config/CORE39_DESIGN.json")
            cfg <- c39_read(p)
            cfg$domains <- lapply(cfg$domains, function(x) as.character(unlist(x, use.names = FALSE)))
            cfg$generation_order <- unlist(cfg$domains, use.names = FALSE)
            c39_assert(length(cfg$domains) == 15L && length(cfg$generation_order) == 39L && !anyDuplicated(cfg$generation_order) && 
                isTRUE(cfg$user_confirmation_received), "CORE39_EXACT_15_39_REQUIRED")
            c39_assert(identical(as.integer(unlist(cfg$bootstrap$ids)), c(0L, 100L)), "CORE39_B0_PLUS_FIXED_B100_ONLY")
            cfg$config_hash <- c39_sha(p)
            cfg
        }
        c39_rng <- function(code) {
            kind <- RNGkind()
            had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
            if (had) old <- get(".Random.seed", .GlobalEnv)
            on.exit({
                do.call(RNGkind, as.list(kind))
                if (had) assign(".Random.seed", old, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", 
                  envir = .GlobalEnv)
            })
            RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
            force(code)
        }
        c39_selection_engine <- function(text) {
            e <- new.env(parent = globalenv())
            for (n in names(text)) e[[n]] <- eval(parse(text = text[[n]]), envir = e)
            e
        }
        c39_unwrap_monitor <- function(f) {
            for (i in 1:4) {
                body_text <- paste(deparse(body(f)), collapse = " ")
                if (!grepl("v5_frame_endpoint", body_text, fixed = TRUE)) break
                e <- environment(f)
                c39_assert(exists("fn", e, inherits = FALSE) && is.function(e$fn), "CORE39_UNKNOWN_SELECTION_WRAPPER")
                f <- e$fn
            }
            c39_assert(!grepl("v5_frame_endpoint", paste(deparse(body(f)), collapse = " "), fixed = TRUE), "CORE39_WRAPPER_NOT_REMOVED")
            f
        }
        c39_members <- function(ctx, cfg) {
            ans <- lapply(cfg$generation_order, function(id) {
                ca <- ctx$features[[id]]
                s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
                c39_assert(is.list(ca) && nrow(s) == 1L && identical(ctx$variables[[id]]$status, "AVAILABLE_FROZEN"), 
                  paste("CORE39_MEMBER_NOT_FROZEN", id))
                type <- as.character(s$type)
                x <- ctx$d_full[[ca$column]]
                c39_assert(type %in% c("continuous", "binary", "categorical", "ordinal", "nominal") && !is.null(x), 
                  paste("CORE39_MEMBER_TYPE_OR_COLUMN", id))
                if (type == "continuous") {
                  allowed <- identical(ca$prepared_transform, "identity") || (identical(id, "clinical__tg") && 
                    identical(ca$prepared_transform, "log") && identical(ca$prepared_unit, "log(mmol/L)") && identical(ca$inverse_transform, 
                    "exp") && is.null(ca$physical_lower_bound) && is.null(ca$physical_upper_bound))
                  c39_assert(allowed, paste("CORE39_UNREVIEWED_TRANSFORM", id))
                }
                list(indicator_id = id, column = ca$column, type = type, levels = if (type != "continuous") levels(x) else NULL, 
                  lower = if (is.null(ca$physical_lower_bound)) -Inf else ca$physical_lower_bound, upper = if (is.null(ca$physical_upper_bound)) Inf else ca$physical_upper_bound, 
                  prepared_transform = ca$prepared_transform, prepared_unit = ca$prepared_unit, raw_unit = ca$raw_unit, 
                  inverse_transform = ca$inverse_transform, technical_covariates = ca$technical_covariates)
            })
            names(ans) <- cfg$generation_order
            c39_assert(!anyDuplicated(vapply(ans, `[[`, "", "column")), "CORE39_DUPLICATE_COORDINATE")
            ans
        }
        c39_selection <- function(engine, d, columns, base) {
            z <- engine$b3_missing(d, columns, base)
            c39_assert(identical(z$complete, complete.cases(d[, unique(c(base, columns)), drop = FALSE])), "CORE39_SELECTION_MASK_DRIFT")
            c39_assert(length(z$w) == sum(z$complete) && all(is.finite(z$w) & z$w > 0), "CORE39_SELECTION_WEIGHT_INVALID")
            z$design_weights <- rep(1, length(z$w))
            z$selection_weights <- as.numeric(z$w)
            z$selection_coefficients <- if (is.null(z$fit)) NULL else stats::coef(z$fit)
            z$selection_formula <- if (is.null(z$fit)) NULL else paste(deparse(stats::formula(z$fit)), collapse = " ")
            z$fit <- NULL
            z
        }
        c39_counts <- function(d, risks, complete) {
            list(N_target = nrow(d), N_joint_fit = sum(complete), N_standardization = sum(complete), N_joint39_missing = sum(!complete), 
                joint39_missing_fraction = mean(!complete), N_families = length(unique(d$cluster_id)), endpoints = lapply(risks, 
                  function(r) list(target_events = sum(r$event), target_competing_deaths = sum(r$death_event), 
                    joint_fit_events = sum(r$event[complete]), joint_fit_competing_deaths = sum(r$death_event[complete]), 
                    followup_max = max(r$exit), followup_not_capped_at_11 = TRUE)))
        }
        c39_member_audit <- function(d, members, complete) {
            lapply(members, function(m) {
                x <- d[[m$column]]
                y <- x[complete]
                list(indicator_id = m$indicator_id, column = m$column, type = m$type, lower = m$lower, upper = m$upper, 
                  prepared_transform = m$prepared_transform, missing_target = sum(is.na(x)), observed_target = sum(!is.na(x)), 
                  zero_on_prepared_scale_target = sum(x == 0, na.rm = TRUE), true_raw_zero_target = if (identical(m$prepared_transform, 
                    "identity")) sum(x == 0, na.rm = TRUE) else NULL, half_value_target = sum(x == 0.5, na.rm = TRUE), 
                  negative_target = if (is.numeric(x)) sum(x < 0, na.rm = TRUE) else NULL, complete_support = if (is.factor(y)) as.list(table(y)) else as.list(stats::quantile(y, 
                    c(0, 0.01, 0.5, 0.99, 1))), observed_range_is_not_physical_bound = TRUE)
            })
        }
        c39_attach_backend <- function(ctx, engine) {
            for (n in ls(ctx$v4, all.names = TRUE)) if (startsWith(n, "ff_") && is.function(ctx$v4[[n]])) engine[[n]] <- ctx$v4[[n]]
            invisible(engine)
        }
        c39_backend_load <- function(root, run, engine, expected) {
            v4root <- Sys.getenv("FD_V4_ROOT", file.path(dirname(root), "border_fullfollowup_b1000_add9_20261030_v1"))
            decision <- c39_read(file.path(v4root, "config/V4_METHOD_DECISION.json"))
            for (n in names(decision$scientific_files_sha256)) c39_assert(identical(c39_sha(file.path(v4root, n)), 
                decision$scientific_files_sha256[[n]]), paste("CORE39_V4_SOURCE_CHANGED", n))
            c39_assert(identical(decision$method_spec_sha256, expected), "CORE39_METHOD_CHANGED")
            env <- new.env(parent = globalenv())
            for (n in c("10.prepare.R", "20.a_reuse.R", "30.outcome.R", "31.aj_backend.R", "61.center_strata.R")) sys.source(file.path(v4root, 
                "f", n), envir = env)
            c39_assert(identical(env$ff_ce_method_spec()$sha256, expected), "CORE39_V4_METHOD_CODE_CHANGED")
            cfg <- c39_read(file.path(v4root, "config/runtime.json"))
            host <- Sys.getenv("FF_HOST", "winpc")
            work <- Sys.getenv("FD_V4_WORK", cfg$work_paths[[host]])
            runtime <- list(R = R.version$version.string, platform = R.version$platform, arch = R.version$arch, 
                pointer_size = .Machine$sizeof.pointer, extension = .Platform$dynlib.ext, isolated_build_contract = "R_CMD_SHLIB_V2_EXPLICIT_RUNTIME_COMPILER_PATH", 
                compiler_paths = Sys.which(c("cc", "clang", "gcc", "arm64-apple-darwin20.0.0-clang")))
            build <- file.path(work, "software", paste0("aj_", substr(c39_sha(file.path(v4root, "src/aj.c")), 1, 
                16), "_", substr(digest::digest(runtime, algo = "sha256"), 1, 12)))
            c39_assert(file.exists(file.path(build, "BUILD_MANIFEST.rds")) && file.exists(file.path(build, paste0("ff_aj", 
                .Platform$dynlib.ext))), "CORE39_EXISTING_AJ_BINARY_REQUIRED")
            backend <- env$ff_aj_initialize(v4root, work, require_compiled = TRUE)
            z <- list(v4 = env, v4_decision = decision, backend = backend)
            c39_attach_backend(z, engine)
            z
        }
        c39_prepare <- function(root, run, engine) {
            c39_assert(Sys.info()[["sysname"]] != "Darwin" && identical(Sys.getenv("FF_HOST"), "winpc"), "CORE39_PREPARATION_WINPC_ONLY")
            root <- normalizePath(root, mustWork = TRUE)
            run <- c39_run(run, TRUE)
            cfg <- c39_config(root)
            private <- file.path(run, "private")
            if (!dir.exists(private)) dir.create(private)
            lockpath <- file.path(private, "INPUT_MANIFEST.json")
            c39_assert(!file.exists(lockpath), "CORE39_INPUT_ALREADY_PREPARED")
            old <- file.path(dirname(root), "border_fine_domain_behavior_joint_20261031_v1")
            loaderpath <- file.path(old, "f/10.inputs.R")
            loaderhash <- c39_read(file.path(root, "config/PROPOSED_DESIGN.json"))$source_files$old_loader$sha256
            c39_assert(identical(c39_sha(loaderpath), loaderhash), "CORE39_ORIGINAL_LOADER_CHANGED")
            loader <- new.env(parent = engine)
            sys.source(loaderpath, envir = loader)
            ctx <- loader$fd_load(old, list(), run)
            members <- c39_members(ctx, cfg)
            base <- unique(c(as.character(unlist(ctx$cfg$base_covariates)), unlist(lapply(members, `[[`, "technical_covariates"))))
            c39_assert(identical(base, as.character(unlist(cfg$base))), "CORE39_ORIGINAL_BASE_CHANGED")
            ep <- ctx$endpoints[c("250.2", "411")]
            c39_assert(all(vapply(ep, function(x) identical(x$status, "AVAILABLE_FROZEN"), FALSE)), "CORE39_ENDPOINT_UNAVAILABLE")
            shared <- ep[[1]]$eligibility & ep[[2]]$eligibility
            c39_assert(is.logical(shared) && !anyNA(shared) && length(shared) == nrow(ctx$d_full), "CORE39_INVALID_TARGET_MASK")
            for (e in names(ep)) c39_assert(is.logical(ep[[e]]$incident_eligibility) && length(ep[[e]]$incident_eligibility) == 
                length(shared) && !anyNA(ep[[e]]$incident_eligibility) && all(ep[[e]]$incident_eligibility[shared]), 
                paste("CORE39_ORIGINAL_INCIDENT_RISKSET_NOT_COVERING_TARGET", e))
            for (id in cfg$generation_order) c39_assert(all(ctx$eligibility[[id]][shared]), paste("CORE39_MEMBER_APPLICABILITY_TARGET_CONFLICT", 
                id))
            rows <- which(shared)
            c39_assert(length(rows) > 0L, "CORE39_EMPTY_TARGET")
            columns <- vapply(members, `[[`, "", "column")
            keep <- unique(c("eid", "cluster_id", base, columns, grep("^aux_", names(ctx$d_full), value = TRUE)))
            d <- ctx$d_full[rows, keep, drop = FALSE]
            rownames(d) <- NULL
            c39_assert(all(as.numeric(as.character(d$S)) %in% 2:8) && all(is.finite(d$B) & d$B >= 0 & d$B <= 1) && 
                !anyNA(d$cluster_id) && !anyDuplicated(d$eid), "CORE39_ORIGINAL_ELIGIBILITY_INVALID")
            risks <- lapply(ep, function(x) x$risk[rows, , drop = FALSE])
            for (r in risks) c39_assert(identical(as.character(r$eid), as.character(d$eid)) && all(is.finite(r$exit) & 
                r$exit > 0 & r$entry == 0 & r$event %in% 0:1 & r$death_event %in% 0:1 & r$event + r$death_event <= 
                1), "CORE39_RISK_ALIGNMENT_INVALID")
            funcs <- c("b3_missing", "b3_formula", "b3_capture", "b3_assert", "b3_boot_index")
            texts <- setNames(lapply(funcs, function(n) paste(deparse(c39_unwrap_monitor(ctx$engine[[n]]), width.cutoff = 500L), 
                collapse = "\n")), funcs)
            selection_engine <- c39_selection_engine(texts)
            selection <- c39_selection(selection_engine, d, columns, base)
            counts <- c39_counts(d, risks, selection$complete)
            counts$N_missing_base_covariates <- sum(!complete.cases(d[, base, drop = FALSE]))
            engine$fd_support(d[selection$complete, , drop = FALSE], members)
            domains <- lapply(cfg$domains, function(ids) vapply(members[ids], `[[`, "", "column"))
            model <- list(members = members, base = base, continuous_rule = cfg$continuous_rule, count_rule = cfg$count_rule)
            object <- list(schema = "CORE39_FROZEN_ACTUAL_INPUTS_V21", data = d, risks = risks, members = members, 
                model = model, domains = domains, parent_rows = rows, parent_frame_hash = ctx$parent_frame_hash, 
                source_input_hash = ctx$input_hash, source_hash = ctx$source_hash, method_hash = ctx$v4_decision$method_spec_sha256, 
                config_hash = cfg$config_hash, target_hash = c39_hash(as.character(d$eid)), cluster_membership_hash = c39_hash(paste(d$eid, 
                  d$cluster_id, sep = ":")), selection_functions = texts, selection_source_hash = c39_hash(texts), 
                B0_selection = selection, counts = counts)
            inputpath <- file.path(private, "CORE39_INPUTS.private.rds")
            c39_rds(object, inputpath)
            drawsdir <- file.path(private, "draws")
            c39_assert(!dir.exists(drawsdir) && dir.create(drawsdir), "CORE39_DRAW_DIRECTORY_EXISTS")
            clusters <- unique(as.character(d$cluster_id))
            ci <- match(as.character(d$cluster_id), clusters)
            sizes <- tabulate(ci, length(clusters))
            draws <- c39_rng(lapply(0:100, function(b) {
                ix <- selection_engine$b3_boot_index(d, as.integer(b), as.integer(cfg$bootstrap$seed))
                multiplicity <- tabulate(ci[ix], length(clusters))/sizes
                c39_assert(all(is.finite(multiplicity) & multiplicity >= 0 & multiplicity == floor(multiplicity)) && 
                  sum(multiplicity) == length(clusters), "CORE39_WHOLE_CLUSTER_DRAW_INVALID")
                p <- file.path(drawsdir, sprintf("bootstrap_%04d.private.rds", b))
                draw <- list(bootstrap_id = as.integer(b), indices = as.integer(ix), cluster_multiplicity = as.integer(multiplicity), 
                  target_hash = object$target_hash, cluster_membership_hash = object$cluster_membership_hash, seed = as.integer(cfg$bootstrap$seed) + 
                    1009L * as.integer(b), RNG_kind = RNGkind())
                c39_rds(draw, p)
                if (b%%100L == 0L) {
                  message("CORE39_ACTUAL_DRAW_FROZEN ", b, "/100")
                  flush.console()
                }
                list(b = as.integer(b), seed = draw$seed, sampling_hash = c39_hash(draw$indices), index_path = p, 
                  index_sha256 = c39_sha(p), index_relative_path = file.path("private/draws", basename(p)), N_resampled = length(ix), 
                  multiplicity_hash = c39_hash(draw$cluster_multiplicity), target_hash = object$target_hash)
            }))
            drawmanifest <- file.path(private, "DRAW_MANIFEST.json")
            c39_json(list(schema = "CORE39_ACTUAL_CLUSTER_DRAWS_V21", unit = "Frozen cluster_id family", mechanism = cfg$bootstrap$mechanism, 
                input_sha256 = c39_sha(inputpath), target_hash = object$target_hash, cluster_membership_hash = object$cluster_membership_hash, 
                RNG_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"), fixed_attempts = 100L, draws = draws), 
                drawmanifest)
            manifest <- list(schema = object$schema, input_path = inputpath, input_relative_path = "private/CORE39_INPUTS.private.rds", 
                input_hash = c39_sha(inputpath), files = list(list(path = inputpath, sha256 = c39_sha(inputpath))), 
                config_hash = cfg$config_hash, source_input_hash = ctx$input_hash, source_hash = ctx$source_hash, 
                method_hash = object$method_hash, selection_source_hash = object$selection_source_hash, target_hash = object$target_hash, 
                cluster_membership_hash = object$cluster_membership_hash, draw_manifest_path = drawmanifest, draw_manifest_hash = c39_sha(drawmanifest), 
                design_weight_rule = cfg$design_weights, selection_rule = cfg$selection, standardization = cfg$standardization, 
                frozen_scale_notes = cfg$frozen_scale_notes, missingness_assumptions = cfg$missingness_assumptions, 
                counts = counts, selection = selection[c("diagnostics", "distribution", "warnings", "selection_coefficients", 
                  "selection_formula")], member_audit = c39_member_audit(d, members, selection$complete), frozen_feature_metadata = ctx$features[c(cfg$generation_order, 
                  "whr", "blood_pressure__MAP", "weight")], frozen_variable_metadata = lapply(ctx$variables[cfg$generation_order], 
                  function(x) x[intersect(names(x), c("indicator_id", "status", "column", "field", "instance", 
                    "encoding", "eligible", "levels", "reference", "ordered", "source_definition", "target_rule", 
                    "invalid_missing_codes", "measurement_instance", "source_label_retained"))]), endpoint_censoring = "Original frozen per-endpoint risk, no cross-disease censoring added", 
                A_Y_models_fitted = 0L)
            c39_json(manifest, lockpath)
            c39_attach_backend(ctx, engine)
            list(input_manifest = lockpath, draw_manifest = drawmanifest, input_hash = manifest$input_hash, counts = counts, 
                config_hash = cfg$config_hash, source_hash = ctx$source_hash, method_hash = object$method_hash)
        }
        c39_load <- function(root, run, engine) {
            c39_assert(identical(Sys.getenv("FD_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "CORE39_REAL_AUTHORITY_REQUIRED")
            if (Sys.info()[["sysname"]] == "Darwin") c39_assert(identical(Sys.getenv("FD_ALLOW_REVIEWED_MAC_COMPUTE"), 
                "YES") && identical(Sys.getenv("FF_HOST"), "macmini"), "CORE39_NO_PARTICIPANT_DATA_ON_CONTROL_MAC")
            run <- c39_run(run)
            cfg <- c39_config(root)
            manifest <- c39_read(file.path(run, "private/INPUT_MANIFEST.json"))
            c39_assert(identical(manifest$config_hash, cfg$config_hash), "CORE39_CONFIG_CHANGED")
            p <- file.path(run, manifest$input_relative_path)
            c39_assert(identical(c39_sha(p), manifest$input_hash), "CORE39_FROZEN_INPUT_CHANGED")
            ctx <- readRDS(p)
            c39_assert(identical(ctx$config_hash, cfg$config_hash) && identical(c39_hash(ctx$selection_functions), 
                manifest$selection_source_hash), "CORE39_INPUT_BINDING_CHANGED")
            dp <- file.path(run, "private/DRAW_MANIFEST.json")
            c39_assert(identical(c39_sha(dp), manifest$draw_manifest_hash), "CORE39_DRAW_MANIFEST_CHANGED")
            dm <- c39_read(dp)
            c39_assert(identical(dm$input_sha256, manifest$input_hash) && length(dm$draws) == 101L && dm$fixed_attempts == 
                100L && identical(vapply(dm$draws, function(x) as.integer(x$b), 0L), 0:100), "CORE39_DRAW_INPUT_MISMATCH")
            ctx$selection_engine <- c39_selection_engine(ctx$selection_functions)
            ctx$run <- run
            ctx$cfg <- cfg
            ctx$input_hash <- manifest$input_hash
            ctx$manifest <- manifest
            ctx$draw_manifest <- dm
            ctx$engine <- engine
            backend <- c39_backend_load(root, run, engine, ctx$method_hash)
            ctx$v4 <- backend$v4
            ctx$backend <- backend$backend
            ctx
        }
        c39_task <- function(ctx, b) {
            c39_assert(length(b) == 1L && is.finite(b) && b == as.integer(b) && b %in% 0:100, "CORE39_FIXED_ID_REQUIRED")
            rec <- ctx$draw_manifest$draws[[as.integer(b) + 1L]]
            c39_assert(rec$b == b, "CORE39_DRAW_ID_MISMATCH")
            p <- file.path(ctx$run, rec$index_relative_path)
            c39_assert(identical(c39_sha(p), rec$index_sha256), "CORE39_ACTUAL_INDEX_CHANGED")
            actual <- readRDS(p)
            ix <- actual$indices
            c39_assert(actual$bootstrap_id == b && identical(actual$target_hash, ctx$target_hash) && identical(actual$cluster_membership_hash, 
                ctx$cluster_membership_hash) && identical(c39_hash(ix), rec$sampling_hash) && identical(c39_hash(actual$cluster_multiplicity), 
                rec$multiplicity_hash) && all(is.finite(ix) & ix >= 1L & ix <= nrow(ctx$data)), "CORE39_ACTUAL_DRAW_BINDING_INVALID")
            d <- ctx$data[ix, , drop = FALSE]
            cols <- vapply(ctx$members, `[[`, "", "column")
            mi <- if (b == 0L) ctx$B0_selection else c39_selection(ctx$selection_engine, d, cols, ctx$model$base)
            c39_assert(length(mi$complete) == nrow(d), "CORE39_DRAW_SELECTION_LENGTH_MISMATCH")
            complete <- mi$complete
            fit <- d[complete, , drop = FALSE]
            risks <- lapply(ctx$risks, function(r) r[ix[complete], , drop = FALSE])
            w <- as.numeric(mi$design_weights * mi$selection_weights)
            ctx$engine$fd_support(fit, ctx$members)
            support <- local({
                members <- ctx$members
                eng <- ctx$engine
                function(m, C = NULL, b = NULL) {
                  eng$fd_support(m, members)
                  dbp <- members[["blood_pressure__DBP"]]$column
                  sbp <- members[["blood_pressure__SBP"]]$column
                  c39_assert(all(m[[sbp]] > m[[dbp]]), "CORE39_SBP_NOT_GREATER_THAN_DBP")
                  TRUE
                }
            })
            fit <- fit[, unique(c(ctx$model$base, cols)), drop = FALSE]
            prepared <- list(fit_data = fit, risks = risks, fit_weights = w, target = fit, standardization_weights = w, 
                design_weights = mi$design_weights, selection_weights = mi$selection_weights, cluster = as.character(d$cluster_id[complete]), 
                joint_target_identity_checked = TRUE, selection_weight_process_checked = TRUE, validate_joint_support = support, 
                actual_draw_proof = c(rec, list(actual_index_sha256 = c39_sha(p), index_path_host = p, new_common_target = TRUE, 
                  old_indices_reused = FALSE)), selection = mi, counts = c39_counts(d, lapply(ctx$risks, function(r) r[ix, 
                  , drop = FALSE]), complete))
            binding <- list(input_hash = ctx$input_hash, config_hash = ctx$config_hash, method_hash = ctx$method_hash, 
                bootstrap_id = as.integer(b), sampling_hash = rec$sampling_hash, actual_index_sha256 = rec$index_sha256, 
                target_hash = c39_hash(as.character(d$eid)), sample_hash = c39_hash(as.character(d$eid[complete])), 
                design_weight_hash = c39_hash(mi$design_weights), selection_weight_hash = c39_hash(mi$selection_weights), 
                weight_hash = c39_hash(w), generation_order = cols, domain_partition_sha256 = c39_hash(ctx$domains), 
                allocation_order = names(ctx$domains), mc_draws = 50L, mc_seed = as.integer(20261101L + 1009L * 
                  b + 100003L))
            list(prepared = prepared, model = ctx$model, domains = ctx$domains, binding = binding, members = ctx$members, 
                base = ctx$model$base, counts = prepared$counts, actual_draw_proof = prepared$actual_draw_proof)
        }
    }), envir = context)
}

