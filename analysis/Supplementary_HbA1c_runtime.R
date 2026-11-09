# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

gs_rt_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

gs_rt_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

gs_rt_sha <- function(path) digest::digest(file = path, algo = "sha256")

gs_rt_norm <- function(path, must = TRUE) normalizePath(path, winslash = "/", mustWork = must)

gs_rt_inside <- function(path, root) {
    p <- tolower(gs_rt_norm(path, FALSE))
    r <- tolower(gs_rt_norm(root))
    identical(p, r) || startsWith(paste0(p, "/"), paste0(r, "/"))
}

gs_rt_task_name <- function() "border_glycaemia_sensitivity_BASE11_B0_20261106_v1"

gs_rt_task_boundary <- function(path) {
    root <- gs_rt_norm(path)
    while (!identical(basename(root), gs_rt_task_name()) && !identical(dirname(root), root)) root <- dirname(root)
    gs_rt_assert(identical(basename(root), gs_rt_task_name()), "EXACT_NAMED_TASK_ANCESTOR_REQUIRED")
    root
}

gs_rt_pins <- function() {
    common <- c(`vendor/core39/config/SOURCE_LOCK.json` = birth_private_setting("PRIVATE_HASH_1"), `vendor/core39/R/core39_models.R` = birth_private_setting("PRIVATE_HASH_2"), 
        `vendor/core39/repairs/exponential_fast_b100_v1/y_model/ethnicity_pool.R` = birth_private_setting("PRIVATE_HASH_3"), 
        `vendor/core39/repairs/exponential_fast_b100_v1/numerics/bounded_aj.R` = birth_private_setting("PRIVATE_HASH_4"), 
        `vendor/core39/repairs/exponential_fast_b100_v1/numerics/piecewise_exponential.c` = birth_private_setting("PRIVATE_HASH_5"))
    current <- c(`R/worker.R` = birth_private_setting("PRIVATE_HASH_6"), `vendor/prior121/R/product_domain_kernel.R` = birth_private_setting("PRIVATE_HASH_7"), 
        `ops/precompute_B0_Y.R` = birth_private_setting("PRIVATE_HASH_8"), `R/y_center_continuous.R` = birth_private_setting("PRIVATE_HASH_9"), 
        `R/product_reverse_isolation.R` = birth_private_setting("PRIVATE_HASH_10"), `R/product_kernel_runtime.R` = birth_private_setting("PRIVATE_HASH_11"), 
        `R/selected_nmr_composition_A.R` = birth_private_setting("PRIVATE_HASH_12"), `R/selected_nmr_composition_diagnostics.R` = birth_private_setting("PRIVATE_HASH_13"), 
        `R/selected_nmr_blockspec.R` = birth_private_setting("PRIVATE_HASH_14"), `R/body_tqe_diagnostics.R` = birth_private_setting("PRIVATE_HASH_15"), 
        `vendor/core39/repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R` = birth_private_setting("PRIVATE_HASH_16"), 
        `R/body_tqe_A.R` = birth_private_setting("PRIVATE_HASH_17"), `R/body_group_pmm_A_v2.R` = birth_private_setting("PRIVATE_HASH_18"), 
        `R/nmr23591_logit_A_v2.R` = birth_private_setting("PRIVATE_HASH_19"))
    y <- c(`vendor/prior121/R/y_prediction_adapter.R` = birth_private_setting("PRIVATE_HASH_20"))
    list(current = c(common, current), original_Y = c(common, y))
}

gs_rt_source_closure <- function(source_root, y_root) {
    pins <- gs_rt_pins()
    roots <- list(current = gs_rt_norm(source_root), original_Y = gs_rt_norm(y_root))
    evidence <- list()
    order <- 0L
    add <- function(scope, relative, expected, kind) {
        gs_rt_assert(is.character(relative) && length(relative) == 1L && !grepl("(^|/)\\.\\.(/|$)|^[A-Za-z]:|^/", 
            relative), "SOURCE_RELATIVE_PATH_REQUIRED")
        path <- gs_rt_norm(file.path(roots[[scope]], relative))
        gs_rt_assert(gs_rt_inside(path, roots[[scope]]) && file.exists(path), "PINNED_SOURCE_BOUNDARY_OR_FILE_MISSING")
        actual <- gs_rt_sha(path)
        gs_rt_assert(identical(actual, expected), paste("PINNED_SOURCE_CHANGED", scope, relative))
        order <<- order + 1L
        evidence[[order]] <<- list(scope = scope, relative_path = relative, path = path, sha256 = actual, kind = kind)
    }
    for (scope in names(pins)) {
        for (relative in names(pins[[scope]])) add(scope, relative, pins[[scope]][[relative]], "STATIC_BYTE_PIN")
        lock <- jsonlite::fromJSON(file.path(roots[[scope]], "vendor/core39/config/SOURCE_LOCK.json"), simplifyVector = FALSE)
        wanted <- if (identical(scope, "current")) 
            c("f/00.core.R", "f/20.joint.R")
        else c("f/30.outcome.R", "f/61.center_strata.R", "R/fast_predict.R")
        for (suffix in wanted) {
            entries <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)
            gs_rt_assert(length(entries) == 1L, paste("UNIQUE_PINNED_SOURCE_REQUIRED", scope, suffix))
            add(scope, file.path("vendor", entries[[1L]]$path), entries[[1L]]$sha256, "BYTE_PINNED_SOURCE_LOCK_DEPENDENCY")
        }
    }
    list(schema = "GS_BASE11_B0_RUNTIME_SOURCE_CLOSURE_V1", roots = roots, files = evidence, source_only_reuse = TRUE, 
        old_fitted_A_Y_or_cache_reads = 0L, closure_hash = gs_rt_hash(lapply(evidence, function(z) z[c("scope", 
            "relative_path", "sha256", "kind")])))
}

gs_rt_no_old_A_options <- function() {
    forbidden <- c("body_tqe.upstream_A", "body_tqe.prefix_A", "body_pmm.old_A11", "body_pmm.prefix_A")
    gs_rt_assert(all(vapply(forbidden, function(n) is.null(getOption(n)), logical(1))), "FRESH_BASE11_REQUIRES_NO_OLD_OR_PREFIX_A_OPTIONS")
    invisible(TRUE)
}

gs_rt_new_path <- function(runtime, out, name) {
    allowed <- if (!is.null(runtime$output_root)) 
        runtime$output_root
    else runtime$task_root
    gs_rt_assert(dir.exists(out) && gs_rt_inside(out, allowed) && !gs_rt_inside(out, runtime$source_root) && !gs_rt_inside(out, 
        runtime$y_root), "NEW_TASK_PRIVATE_OUTPUT_BOUNDARY_REQUIRED")
    gs_rt_assert(identical(basename(name), name) && nzchar(name), "ONE_OUTPUT_BASENAME_REQUIRED")
    path <- file.path(out, name)
    gs_rt_assert(!file.exists(path), paste("PRESERVE_EXISTING_RUNTIME_OUTPUT", name))
    path
}

gs_rt_rds <- function(runtime, x, out, name) {
    path <- gs_rt_new_path(runtime, out, name)
    saveRDS(x, path, version = 2L)
    list(path = basename(path), sha256 = gs_rt_sha(path), bytes = unname(file.info(path)$size))
}

gs_rt_json <- function(runtime, x, out, name) {
    path <- gs_rt_new_path(runtime, out, name)
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null", na = "null", matrix = "rowmajor")
    list(path = basename(path), sha256 = gs_rt_sha(path), bytes = unname(file.info(path)$size))
}

gs_load_runtime <- function(source_root, y_root, state_root, output_root = Sys.getenv("GS44_OUTPUT_ROOT", unset = "")) {
    for (p in c("digest", "jsonlite", "survival", "splines", "nnet")) gs_rt_assert(requireNamespace(p, quietly = TRUE), 
        paste("INSTALLED_PACKAGE_REQUIRED_NO_INSTALL", p))
    gs_rt_assert(identical(unname(Sys.info()[["sysname"]]), "Linux") && identical(as.character(getRversion()), 
        "4.5.3"), "AUTHORIZED_WINPC_WSL_LINUX_R453_RUNTIME_ONLY")
    source_root <- gs_rt_norm(source_root)
    y_root <- gs_rt_norm(y_root)
    state_root <- gs_rt_norm(state_root)
    gs_rt_assert(dir.exists(state_root), "NEW_TASK_STATE_ROOT_REQUIRED")
    task_root <- gs_rt_task_boundary(state_root)
    gs_rt_assert(!gs_rt_inside(state_root, source_root) && !gs_rt_inside(state_root, y_root), "ISOLATED_NAMED_NEW_TASK_STATE_REQUIRED")
    if (!nzchar(output_root)) 
        output_root <- task_root
    output_root <- gs_rt_norm(output_root)
    gs_rt_task_boundary(output_root)
    gs_rt_assert(dir.exists(output_root) && !gs_rt_inside(output_root, source_root) && !gs_rt_inside(output_root, 
        y_root), "EXPLICIT_NEW_TASK_OUTPUT_ROOT_REQUIRED")
    gs_rt_no_old_A_options()
    evidence <- gs_rt_source_closure(source_root, y_root)
    methods <- new.env(parent = globalenv())
    sourced <- list()
    src <- function(root, relative, envir) {
        path <- gs_rt_norm(file.path(root, relative))
        sys.source(path, envir = envir)
        sourced[[length(sourced) + 1L]] <<- list(path = path, sha256 = gs_rt_sha(path))
    }
    src(source_root, "vendor/prior121/R/product_domain_kernel.R", methods)
    src(y_root, "vendor/prior121/R/y_prediction_adapter.R", methods)
    src(source_root, "ops/precompute_B0_Y.R", methods)
    src(source_root, "R/y_center_continuous.R", methods)
    for (name in c("product_reverse_isolation.R", "product_kernel_runtime.R", "selected_nmr_composition_A.R", "selected_nmr_composition_diagnostics.R", 
        "selected_nmr_blockspec.R", "body_tqe_diagnostics.R")) src(source_root, file.path("R", name), methods)
    methods$nmr_install_product_runtime(methods)
    build_parent <- file.path(state_root, "native_double")
    if (!dir.exists(build_parent)) 
        gs_rt_assert(dir.create(build_parent), "PRIVATE_NATIVE_PARENT_CREATE_FAILED")
    build <- file.path(build_parent, paste0("fresh_B0_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%dT%H%M%OS6", 
        tz = "UTC")))
    gs_rt_assert(gs_rt_inside(build, task_root) && !file.exists(build), "NEW_ISOLATED_NATIVE_BUILD_REQUIRED")
    engine <- methods$new121_y_engine(file.path(y_root, "vendor/core39"), source_parent = file.path(y_root, "vendor"), 
        build_dir = build)
    lock <- jsonlite::fromJSON(file.path(source_root, "vendor/core39/config/SOURCE_LOCK.json"), simplifyVector = FALSE)
    for (suffix in c("f/00.core.R", "f/20.joint.R")) {
        entry <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)[[1L]]
        src(source_root, file.path("vendor", entry$path), engine)
    }
    engine$c39_install_models(engine)
    src(source_root, "vendor/core39/repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R", engine)
    engine$c39_mixeda_install(engine)
    src(source_root, "R/body_tqe_A.R", engine)
    engine$body_tqe_install(engine)
    src(source_root, "R/body_group_pmm_A_v2.R", engine)
    engine$body_pmm_install(engine)
    src(source_root, "R/nmr23591_logit_A_v2.R", methods)
    methods$nmr23591_install(methods)
    gs_rt_assert(identical(engine$body_pmm_revision(), "SELECTED44_GROUPED_BODY_PREDICTIVE_MEAN_CELLS_V2") && identical(engine$ff_ce_method_spec()$id, 
        "CORE39_Y_ASIAN_CHINESE_POOL_V1_20261101"), "ORIGINAL_BASE_A_Y_METHOD_ID_CHANGED")
    gs_rt_assert(identical(engine$new121_numeric_backend$manifest$source_sha256, gs_rt_pins()$original_Y[["vendor/core39/repairs/exponential_fast_b100_v1/numerics/piecewise_exponential.c"]]), 
        "FRESH_NATIVE_SOURCE_CHANGED")
    evidence$explicit_source_order <- sourced
    evidence$engine_internal_Y_source_order <- engine$new121_source_audit
    evidence$native <- list(build_dir = build, manifest = engine$new121_numeric_backend$manifest, build_manifest_sha256 = gs_rt_sha(file.path(build, 
        "BUILD_MANIFEST.rds")), fresh_host_local_build = TRUE)
    evidence$software <- list(R = R.version.string, platform = R.version$platform, packages = setNames(lapply(c("survival", 
        "splines", "digest", "jsonlite", "nnet"), function(p) as.character(utils::packageVersion(p))), c("survival", 
        "splines", "digest", "jsonlite", "nnet")))
    message("B0 runtime source: ", source_root, "; original Y: ", y_root, "; private state: ", state_root, "; private output: ", 
        output_root)
    list(engine = engine, methods = methods, source_pin_evidence = evidence, source_root = source_root, y_root = y_root, 
        state_root = state_root, task_root = task_root, output_root = output_root)
}

gs_rt_validate_inputs <- function(runtime, d, members, w) {
    gs_rt_no_old_A_options()
    gs_rt_assert(is.data.frame(d) && nrow(d) > 1L && length(members) == 11L && identical(names(members), runtime$engine$body_tqe_ids()) && 
        all(c("B", "center", "cluster_id") %in% names(d)) && length(w) == nrow(d) && all(is.finite(w) & w > 0) && 
        is.finite(sum(w)) && !anyNA(d$cluster_id), "EXACT_BASE11_CURRENT_SUBSET_POSITIVE_WEIGHTS_FAMILY_CLUSTER_REQUIRED")
    runtime$engine$c39_support(d, members)
    invisible(TRUE)
}

gs_validate_generated_A <- function(runtime, A, d, w, seed, out) {
    gs_rt_validate_inputs(runtime, d, A$members, w)
    engine <- runtime$engine
    methods <- runtime$methods
    gs_rt_assert(identical(A$body_pmm$revision, engine$body_pmm_revision()) && is.null(A$body_pmm$reused_old_provenance) && 
        length(A$body_pmm$reused_old_ids) == 0L && identical(names(A$models), names(A$members)), "FRESH_ALL_ELEVEN_A_REQUIRED_NO_OLD_CLINICAL_OR_UPSTREAM_REUSE")
    body <- engine$body_pmm_diagnostics(A, d, w, engine, methods, seed, width = 11L, K = 3L)
    before <- methods$nmrA_rng()
    on.exit({
        do.call(RNGkind, as.list(before$kind))
        if (before$present) assign(".Random.seed", before$seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, 
            inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
    chunks <- vector("list", 3L)
    cols <- vapply(Filter(function(m) identical(m$type, "continuous"), A$members), `[[`, "", "column")
    for (k in seq_len(3L)) {
        U <- matrix(stats::runif(nrow(d) * 11L), nrow(d), byrow = TRUE)
        draw <- engine$c39_draw_a(A, d, d$B, U)
        engine$c39_support(draw, A$members)
        chunks[[k]] <- draw[, cols, drop = FALSE]
    }
    make <- function(raw) list(raw = raw, selected_spearman = data.frame(), physical_pass = TRUE)
    comparison <- methods$nmr_compare_summaries(make(methods$nmrD_table(as.list(d[cols]), w, "SELECTED_ORIGINAL_CONTINUOUS")), 
        make(methods$nmrD_table(as.list(do.call(rbind, chunks)), rep(w/3L, 3L), "SELECTED_ORIGINAL_CONTINUOUS")))
    clinical <- list(schema = "SELECTED_ACTUAL_OBSERVED_B_ENGINEERING_V1", actual_rows = nrow(d), Ksmall = 3L, 
        observed = TRUE, synthetic = FALSE, comparison = comparison, passed = comparison$engineering_pass, official_MC_integrations = 0L, 
        fitted_models = 0L, physical_support_checked_every_generated_full_world = TRUE)
    record <- list(schema = "GS_BASE11_B0_FRESH_A_DISTRIBUTION_VALIDATION_V1", passed = isTRUE(body$passed) && 
        isTRUE(clinical$passed), actual_rows = nrow(d), weight_hash = gs_rt_hash(as.numeric(w)), observed = TRUE, 
        synthetic = FALSE, body = body, clinical = clinical, no_refit_or_law_changes = TRUE, no_old_A_reuse = TRUE)
    artifact <- gs_rt_json(runtime, record, out, "ACTUAL_A_DISTRIBUTION_VALIDATION.json")
    gs_rt_assert(record$passed, "ACTUAL_BASE11_A_DISTRIBUTION_OR_PHYSICAL_SUPPORT_FAILURE")
    c(record, list(artifacts = list(artifact)))
}

gs_rt_probe_rows <- function(d) {
    gs_rt_assert(is.data.frame(d) && nrow(d) > 0L && "center" %in% names(d) && !anyNA(d$center), "ACTUAL_CENTER_ROWS_REQUIRED")
    unique(c(match(unique(as.character(d$center)), as.character(d$center)), as.integer(round(seq(1, nrow(d), length.out = min(16L, 
        nrow(d)))))))
}

gs_rt_validate_y <- function(runtime, y, design_fn, predictor, d, members, base, domains, w) {
    engine <- runtime$engine
    methods <- runtime$methods
    take <- gs_rt_probe_rows(d)
    C <- d[take, , drop = FALSE]
    cols <- vapply(members, `[[`, "", "column")
    M <- C[, cols, drop = FALSE]
    records <- list()
    max_mass <- 0
    for (b in 0:1) {
        probe <- C
        probe$B <- b
        X <- design_fn(b, M, C)
        original <- engine$c39_design(probe, members, base, attr(design_fn, "reference"))$X
        gs_rt_assert(identical(dimnames(X), dimnames(original)) && identical(dim(X), dim(original)) && max(abs(X - 
            original)) == 0, "CACHED_ORIGINAL_Y_DESIGN_NOT_EXACT")
        X <- methods$new121_y_apply_alias_plan(X, y$alias_plan)
        check <- predictor$validate(b, M, C, weights = w[take])
        gs_rt_assert(isTRUE(check$validation$within_engineering_reference), "ORIGINAL_ADAPTER_NATIVE_CIF_ENGINEERING_REFERENCE_EXCEEDED")
        for (endpoint in c("IHD", "T2D")) {
            cf <- y$compiled[[endpoint]]
            prepared <- engine$fd_fp_prepare(X, cf, engine)
            lp_error <- 0
            for (g in unique(prepared$group)) for (j in 1:2) {
                ii <- which(prepared$group == g)
                nd <- as.data.frame(prepared$designs[[j]][ii, , drop = FALSE])
                names(nd) <- cf$object$cause_references[[j]]$safe_names
                nd$.ce_center <- factor(g, levels = cf$object$center_reference$observed_levels)
                got <- engine$ff_y_capture(as.numeric(stats::predict(cf$object$fit[[j]], newdata = nd, type = "lp", 
                  reference = "zero")))
                gs_rt_assert(!length(got$warnings), "NATIVE_LP_VALIDATION_WARNING")
                fast_lp <- as.numeric(prepared$designs[[j]][ii, , drop = FALSE] %*% cf$beta[[j]])
                lp_error <- max(lp_error, abs(got$value - fast_lp))
            }
            fast <- engine$fd_fp_predict_matrix(cf, X, engine)
            native <- engine$ff_ce_predict_y_matrix(cf$object, X, 11)
            gs_rt_assert(!length(attr(native, "warnings")), "NATIVE_CIF_VALIDATION_WARNING")
            cif_error <- max(abs(fast - native))
            mass <- max(abs(rowSums(check$all_states[[endpoint]]) - 1))
            max_mass <- max(max_mass, mass)
            gs_rt_assert(is.finite(lp_error) && lp_error <= 1e-10 && is.finite(cif_error) && cif_error <= 1e-12 && 
                is.finite(mass) && mass <= 1e-06, "ORIGINAL_FAST_NATIVE_LP_CIF_OR_MASS_VALIDATION_FAILED")
            rr <- matrix(exp(vapply(1:2, function(j) as.numeric(prepared$designs[[j]][1L, , drop = FALSE] %*% cf$beta[[j]]), 
                numeric(1))), 1L, 2L)
            h <- cf$hazards[[prepared$group[1L]]]
            slow <- engine$c39_num_quadrature(h, rr)
            slow_error <- max(abs(slow - fast[1L, , drop = FALSE]))
            gs_rt_assert(is.finite(slow_error) && slow_error <= 1e-10, "ORIGINAL_SLOW_QUADRATURE_REFERENCE_FAILED")
            records[[length(records) + 1L]] <- list(exposure = b, endpoint = endpoint, actual_rows = nrow(C), design_max_abs = 0, 
                lp_max_abs = lp_error, cif_all_states_max_abs = cif_error, lp_tolerance = 1e-10, cif_tolerance = 1e-12, 
                slow_entire_hazard_intervals = nrow(h), slow_actual_rows = 1L, slow_reference = "ORIGINAL_c39_num_quadrature_COMPLETE_ACTUAL_HAZARD_GRID", 
                slow_max_abs = slow_error, slow_tolerance = 1e-10, adapter = check$validation, max_mass_residual = mass, 
                passed = TRUE)
        }
    }
    list(schema = "GS_BASE11_B0_ACTUAL_Y_FAST_NATIVE_SLOW_VALIDATION_V1", passed = TRUE, observed = TRUE, synthetic = FALSE, 
        actual_rows = length(take), actual_centers = length(unique(as.character(C$center))), exposures = 0:1, endpoints = c("IHD", 
            "T2D"), maximum_mass_residual = max_mass, genuine_probability_check = TRUE, four_slow_rowpaths = TRUE, 
        complete_baseline_grid = TRUE, no_new_fits = TRUE, consumes_RNG = FALSE, checks = records)
}

gs_fit_and_validate_y <- function(runtime, d, members, base, domains, risks, w, out, progress) {
    gs_rt_validate_inputs(runtime, d, members, w)
    gs_rt_assert(is.function(progress) && is.list(risks) && setequal(names(risks), c("IHD", "T2D")) && all(vapply(risks, 
        nrow, integer(1)) == nrow(d)), "EXACT_CURRENT_PAIRED_RISK_FRAMES_AND_PROGRESS_REQUIRED")
    engine <- runtime$engine
    methods <- runtime$methods
    progress("FRESH_ORIGINAL_BASE11_Y_DESIGN")
    design_fn <- methods$new121_y_design_factory(d, members, base, domains, engine, cache_entries = 2L)
    reference <- gs_rt_rds(runtime, attr(design_fn, "reference"), out, "Y_DESIGN_REFERENCE.private.rds")
    owner <- attr(design_fn, "column_owner")
    mediator_cols <- names(owner)[owner != "BASE"]
    audits <- list()
    progress("FIT_FOUR_FRESH_FAMILY_ROBUST_CAUSE_SPECIFIC_Y_BRANCHES")
    y <- methods$new121_fit_y(engine, attr(design_fn, "fullX"), risks, w, d$cluster_id, mediator_cols, centre_actual = d$center, 
        diagnostics = FALSE, audit_callback = function(record) {
            audits[[length(audits) + 1L]] <<- gs_rt_rds(runtime, record, out, sprintf("Y_NATIVE_AUDIT_%04d.private.rds", 
                length(audits) + 1L))
            progress("FRESH_Y_NATIVE_AUDIT", record[intersect(c("stage", "endpoint", "cause", "warnings", "warning", 
                "status", "reason"), names(record))])
        })
    gs_rt_assert(identical(y$method$tau, 11) && identical(y$method$fit_followup_cap, Inf) && identical(y$method$ties, 
        "efron") && isTRUE(y$method$robust_family_cluster), "ORIGINAL_FULLFOLLOWUP_TAU11_EFRON_FAMILY_METHOD_CHANGED")
    for (ep in c("IHD", "T2D")) {
        fit <- y$fits[[ep]]
        expected_risk <- engine$ff_y_risk(risks[[ep]], 11)
        gs_rt_assert(identical(names(fit$fit), c("disease", "death")) && identical(as.logical(fit$stratified), 
            c(FALSE, TRUE)) && identical(as.numeric(fit$endpoint$time), as.numeric(expected_risk$time)) && identical(as.numeric(fit$endpoint$status), 
            as.numeric(expected_risk$status)) && all(vapply(fit$fit, function(f) !is.null(f$naive.var) && f$n == 
            nrow(d) && identical(f$method, "efron"), logical(1))), "FRESH_Y_BRANCH_CLUSTER_FULLFOLLOWUP_OR_CAUSE_STRATA_CHANGED")
    }
    yr <- gs_rt_rds(runtime, y, out, "FRESH_Y.private.rds")
    manifest <- list(schema = "GS_BASE11_B0_FRESH_Y_MANIFEST_V1", complete = TRUE, fresh_four_branch_fit = TRUE, 
        old_A_Y_cache_reads = 0L, file = yr, reference = reference, private_audits = audits, source_closure_hash = runtime$source_pin_evidence$closure_hash, 
        model = "S44_BASE11", bootstrap_id = 0L, rows = nrow(d), weights_hash = gs_rt_hash(as.numeric(w)), cluster_hash = gs_rt_hash(d$cluster_id), 
        risk_hash = gs_rt_hash(risks), members_hash = gs_rt_hash(members), base = base, domains = domains, method = y$method, 
        fullX_hash = gs_rt_hash(attr(design_fn, "fullX")))
    mr <- gs_rt_json(runtime, manifest, out, "FRESH_Y_MANIFEST.json")
    predictor <- methods$new121_y_prediction(y, design_fn, engine, domains, cache_entries = 2L)
    progress("VALIDATE_ACTUAL_Y_FAST_NATIVE_AND_SLOW_ORIGINAL_REFERENCES")
    validation <- gs_rt_validate_y(runtime, y, design_fn, predictor, d, members, base, domains, w)
    vr <- gs_rt_json(runtime, validation, out, "FIRST_REAL_RISK_BLOCK.json")
    progress("FRESH_ORIGINAL_SCIENCE_Y_READY", list(actual_validation_rows = validation$actual_rows, four_fresh_branches = TRUE))
    list(y = y, design_fn = design_fn, predictor = predictor, predict_risk = predictor$predict, validation = validation, 
        artifacts = list(fit = yr, reference = reference, manifest = mr, validation = vr, audits = audits))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        gs_rt_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        gs_rt_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        gs_rt_sha <- function(path) digest::digest(file = path, algo = "sha256")
        gs_rt_norm <- function(path, must = TRUE) normalizePath(path, winslash = "/", mustWork = must)
        gs_rt_inside <- function(path, root) {
            p <- tolower(gs_rt_norm(path, FALSE))
            r <- tolower(gs_rt_norm(root))
            identical(p, r) || startsWith(paste0(p, "/"), paste0(r, "/"))
        }
        gs_rt_task_name <- function() "border_glycaemia_sensitivity_BASE11_B0_20261106_v1"
        gs_rt_task_boundary <- function(path) {
            root <- gs_rt_norm(path)
            while (!identical(basename(root), gs_rt_task_name()) && !identical(dirname(root), root)) root <- dirname(root)
            gs_rt_assert(identical(basename(root), gs_rt_task_name()), "EXACT_NAMED_TASK_ANCESTOR_REQUIRED")
            root
        }
        gs_rt_pins <- function() {
            common <- c(`vendor/core39/config/SOURCE_LOCK.json` = birth_private_setting("PRIVATE_HASH_21"), `vendor/core39/R/core39_models.R` = birth_private_setting("PRIVATE_HASH_22"), 
                `vendor/core39/repairs/exponential_fast_b100_v1/y_model/ethnicity_pool.R` = birth_private_setting("PRIVATE_HASH_23"), 
                `vendor/core39/repairs/exponential_fast_b100_v1/numerics/bounded_aj.R` = birth_private_setting("PRIVATE_HASH_24"), 
                `vendor/core39/repairs/exponential_fast_b100_v1/numerics/piecewise_exponential.c` = birth_private_setting("PRIVATE_HASH_25"))
            current <- c(`R/worker.R` = birth_private_setting("PRIVATE_HASH_26"), `vendor/prior121/R/product_domain_kernel.R` = birth_private_setting("PRIVATE_HASH_27"), 
                `ops/precompute_B0_Y.R` = birth_private_setting("PRIVATE_HASH_28"), `R/y_center_continuous.R` = birth_private_setting("PRIVATE_HASH_29"), 
                `R/product_reverse_isolation.R` = birth_private_setting("PRIVATE_HASH_30"), `R/product_kernel_runtime.R` = birth_private_setting("PRIVATE_HASH_31"), 
                `R/selected_nmr_composition_A.R` = birth_private_setting("PRIVATE_HASH_32"), `R/selected_nmr_composition_diagnostics.R` = birth_private_setting("PRIVATE_HASH_33"), 
                `R/selected_nmr_blockspec.R` = birth_private_setting("PRIVATE_HASH_34"), `R/body_tqe_diagnostics.R` = birth_private_setting("PRIVATE_HASH_35"), 
                `vendor/core39/repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R` = birth_private_setting("PRIVATE_HASH_36"), 
                `R/body_tqe_A.R` = birth_private_setting("PRIVATE_HASH_37"), `R/body_group_pmm_A_v2.R` = birth_private_setting("PRIVATE_HASH_38"), 
                `R/nmr23591_logit_A_v2.R` = birth_private_setting("PRIVATE_HASH_39"))
            y <- c(`vendor/prior121/R/y_prediction_adapter.R` = birth_private_setting("PRIVATE_HASH_40"))
            list(current = c(common, current), original_Y = c(common, y))
        }
        gs_rt_source_closure <- function(source_root, y_root) {
            pins <- gs_rt_pins()
            roots <- list(current = gs_rt_norm(source_root), original_Y = gs_rt_norm(y_root))
            evidence <- list()
            order <- 0L
            add <- function(scope, relative, expected, kind) {
                gs_rt_assert(is.character(relative) && length(relative) == 1L && !grepl("(^|/)\\.\\.(/|$)|^[A-Za-z]:|^/", 
                  relative), "SOURCE_RELATIVE_PATH_REQUIRED")
                path <- gs_rt_norm(file.path(roots[[scope]], relative))
                gs_rt_assert(gs_rt_inside(path, roots[[scope]]) && file.exists(path), "PINNED_SOURCE_BOUNDARY_OR_FILE_MISSING")
                actual <- gs_rt_sha(path)
                gs_rt_assert(identical(actual, expected), paste("PINNED_SOURCE_CHANGED", scope, relative))
                order <<- order + 1L
                evidence[[order]] <<- list(scope = scope, relative_path = relative, path = path, sha256 = actual, 
                  kind = kind)
            }
            for (scope in names(pins)) {
                for (relative in names(pins[[scope]])) add(scope, relative, pins[[scope]][[relative]], "STATIC_BYTE_PIN")
                lock <- jsonlite::fromJSON(file.path(roots[[scope]], "vendor/core39/config/SOURCE_LOCK.json"), 
                  simplifyVector = FALSE)
                wanted <- if (identical(scope, "current")) c("f/00.core.R", "f/20.joint.R") else c("f/30.outcome.R", 
                  "f/61.center_strata.R", "R/fast_predict.R")
                for (suffix in wanted) {
                  entries <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)
                  gs_rt_assert(length(entries) == 1L, paste("UNIQUE_PINNED_SOURCE_REQUIRED", scope, suffix))
                  add(scope, file.path("vendor", entries[[1L]]$path), entries[[1L]]$sha256, "BYTE_PINNED_SOURCE_LOCK_DEPENDENCY")
                }
            }
            list(schema = "GS_BASE11_B0_RUNTIME_SOURCE_CLOSURE_V1", roots = roots, files = evidence, source_only_reuse = TRUE, 
                old_fitted_A_Y_or_cache_reads = 0L, closure_hash = gs_rt_hash(lapply(evidence, function(z) z[c("scope", 
                  "relative_path", "sha256", "kind")])))
        }
        gs_rt_no_old_A_options <- function() {
            forbidden <- c("body_tqe.upstream_A", "body_tqe.prefix_A", "body_pmm.old_A11", "body_pmm.prefix_A")
            gs_rt_assert(all(vapply(forbidden, function(n) is.null(getOption(n)), logical(1))), "FRESH_BASE11_REQUIRES_NO_OLD_OR_PREFIX_A_OPTIONS")
            invisible(TRUE)
        }
        gs_rt_new_path <- function(runtime, out, name) {
            allowed <- if (!is.null(runtime$output_root)) runtime$output_root else runtime$task_root
            gs_rt_assert(dir.exists(out) && gs_rt_inside(out, allowed) && !gs_rt_inside(out, runtime$source_root) && 
                !gs_rt_inside(out, runtime$y_root), "NEW_TASK_PRIVATE_OUTPUT_BOUNDARY_REQUIRED")
            gs_rt_assert(identical(basename(name), name) && nzchar(name), "ONE_OUTPUT_BASENAME_REQUIRED")
            path <- file.path(out, name)
            gs_rt_assert(!file.exists(path), paste("PRESERVE_EXISTING_RUNTIME_OUTPUT", name))
            path
        }
        gs_rt_rds <- function(runtime, x, out, name) {
            path <- gs_rt_new_path(runtime, out, name)
            saveRDS(x, path, version = 2L)
            list(path = basename(path), sha256 = gs_rt_sha(path), bytes = unname(file.info(path)$size))
        }
        gs_rt_json <- function(runtime, x, out, name) {
            path <- gs_rt_new_path(runtime, out, name)
            jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null", na = "null", 
                matrix = "rowmajor")
            list(path = basename(path), sha256 = gs_rt_sha(path), bytes = unname(file.info(path)$size))
        }
        gs_load_runtime <- function(source_root, y_root, state_root, output_root = Sys.getenv("GS44_OUTPUT_ROOT", 
            unset = "")) {
            for (p in c("digest", "jsonlite", "survival", "splines", "nnet")) gs_rt_assert(requireNamespace(p, 
                quietly = TRUE), paste("INSTALLED_PACKAGE_REQUIRED_NO_INSTALL", p))
            gs_rt_assert(identical(unname(Sys.info()[["sysname"]]), "Linux") && identical(as.character(getRversion()), 
                "4.5.3"), "AUTHORIZED_WINPC_WSL_LINUX_R453_RUNTIME_ONLY")
            source_root <- gs_rt_norm(source_root)
            y_root <- gs_rt_norm(y_root)
            state_root <- gs_rt_norm(state_root)
            gs_rt_assert(dir.exists(state_root), "NEW_TASK_STATE_ROOT_REQUIRED")
            task_root <- gs_rt_task_boundary(state_root)
            gs_rt_assert(!gs_rt_inside(state_root, source_root) && !gs_rt_inside(state_root, y_root), "ISOLATED_NAMED_NEW_TASK_STATE_REQUIRED")
            if (!nzchar(output_root)) output_root <- task_root
            output_root <- gs_rt_norm(output_root)
            gs_rt_task_boundary(output_root)
            gs_rt_assert(dir.exists(output_root) && !gs_rt_inside(output_root, source_root) && !gs_rt_inside(output_root, 
                y_root), "EXPLICIT_NEW_TASK_OUTPUT_ROOT_REQUIRED")
            gs_rt_no_old_A_options()
            evidence <- gs_rt_source_closure(source_root, y_root)
            methods <- new.env(parent = globalenv())
            sourced <- list()
            src <- function(root, relative, envir) {
                path <- gs_rt_norm(file.path(root, relative))
                sys.source(path, envir = envir)
                sourced[[length(sourced) + 1L]] <<- list(path = path, sha256 = gs_rt_sha(path))
            }
            src(source_root, "vendor/prior121/R/product_domain_kernel.R", methods)
            src(y_root, "vendor/prior121/R/y_prediction_adapter.R", methods)
            src(source_root, "ops/precompute_B0_Y.R", methods)
            src(source_root, "R/y_center_continuous.R", methods)
            for (name in c("product_reverse_isolation.R", "product_kernel_runtime.R", "selected_nmr_composition_A.R", 
                "selected_nmr_composition_diagnostics.R", "selected_nmr_blockspec.R", "body_tqe_diagnostics.R")) src(source_root, 
                file.path("R", name), methods)
            methods$nmr_install_product_runtime(methods)
            build_parent <- file.path(state_root, "native_double")
            if (!dir.exists(build_parent)) gs_rt_assert(dir.create(build_parent), "PRIVATE_NATIVE_PARENT_CREATE_FAILED")
            build <- file.path(build_parent, paste0("fresh_B0_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%dT%H%M%OS6", 
                tz = "UTC")))
            gs_rt_assert(gs_rt_inside(build, task_root) && !file.exists(build), "NEW_ISOLATED_NATIVE_BUILD_REQUIRED")
            engine <- methods$new121_y_engine(file.path(y_root, "vendor/core39"), source_parent = file.path(y_root, 
                "vendor"), build_dir = build)
            lock <- jsonlite::fromJSON(file.path(source_root, "vendor/core39/config/SOURCE_LOCK.json"), simplifyVector = FALSE)
            for (suffix in c("f/00.core.R", "f/20.joint.R")) {
                entry <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)[[1L]]
                src(source_root, file.path("vendor", entry$path), engine)
            }
            engine$c39_install_models(engine)
            src(source_root, "vendor/core39/repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R", engine)
            engine$c39_mixeda_install(engine)
            src(source_root, "R/body_tqe_A.R", engine)
            engine$body_tqe_install(engine)
            src(source_root, "R/body_group_pmm_A_v2.R", engine)
            engine$body_pmm_install(engine)
            src(source_root, "R/nmr23591_logit_A_v2.R", methods)
            methods$nmr23591_install(methods)
            gs_rt_assert(identical(engine$body_pmm_revision(), "SELECTED44_GROUPED_BODY_PREDICTIVE_MEAN_CELLS_V2") && 
                identical(engine$ff_ce_method_spec()$id, "CORE39_Y_ASIAN_CHINESE_POOL_V1_20261101"), "ORIGINAL_BASE_A_Y_METHOD_ID_CHANGED")
            gs_rt_assert(identical(engine$new121_numeric_backend$manifest$source_sha256, gs_rt_pins()$original_Y[["vendor/core39/repairs/exponential_fast_b100_v1/numerics/piecewise_exponential.c"]]), 
                "FRESH_NATIVE_SOURCE_CHANGED")
            evidence$explicit_source_order <- sourced
            evidence$engine_internal_Y_source_order <- engine$new121_source_audit
            evidence$native <- list(build_dir = build, manifest = engine$new121_numeric_backend$manifest, build_manifest_sha256 = gs_rt_sha(file.path(build, 
                "BUILD_MANIFEST.rds")), fresh_host_local_build = TRUE)
            evidence$software <- list(R = R.version.string, platform = R.version$platform, packages = setNames(lapply(c("survival", 
                "splines", "digest", "jsonlite", "nnet"), function(p) as.character(utils::packageVersion(p))), 
                c("survival", "splines", "digest", "jsonlite", "nnet")))
            message("B0 runtime source: ", source_root, "; original Y: ", y_root, "; private state: ", state_root, 
                "; private output: ", output_root)
            list(engine = engine, methods = methods, source_pin_evidence = evidence, source_root = source_root, 
                y_root = y_root, state_root = state_root, task_root = task_root, output_root = output_root)
        }
        gs_rt_validate_inputs <- function(runtime, d, members, w) {
            gs_rt_no_old_A_options()
            gs_rt_assert(is.data.frame(d) && nrow(d) > 1L && length(members) == 11L && identical(names(members), 
                runtime$engine$body_tqe_ids()) && all(c("B", "center", "cluster_id") %in% names(d)) && length(w) == 
                nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)) && !anyNA(d$cluster_id), "EXACT_BASE11_CURRENT_SUBSET_POSITIVE_WEIGHTS_FAMILY_CLUSTER_REQUIRED")
            runtime$engine$c39_support(d, members)
            invisible(TRUE)
        }
        gs_validate_generated_A <- function(runtime, A, d, w, seed, out) {
            gs_rt_validate_inputs(runtime, d, A$members, w)
            engine <- runtime$engine
            methods <- runtime$methods
            gs_rt_assert(identical(A$body_pmm$revision, engine$body_pmm_revision()) && is.null(A$body_pmm$reused_old_provenance) && 
                length(A$body_pmm$reused_old_ids) == 0L && identical(names(A$models), names(A$members)), "FRESH_ALL_ELEVEN_A_REQUIRED_NO_OLD_CLINICAL_OR_UPSTREAM_REUSE")
            body <- engine$body_pmm_diagnostics(A, d, w, engine, methods, seed, width = 11L, K = 3L)
            before <- methods$nmrA_rng()
            on.exit({
                do.call(RNGkind, as.list(before$kind))
                if (before$present) assign(".Random.seed", before$seed, .GlobalEnv) else if (exists(".Random.seed", 
                  .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            }, add = TRUE)
            set.seed(as.integer(seed))
            chunks <- vector("list", 3L)
            cols <- vapply(Filter(function(m) identical(m$type, "continuous"), A$members), `[[`, "", "column")
            for (k in seq_len(3L)) {
                U <- matrix(stats::runif(nrow(d) * 11L), nrow(d), byrow = TRUE)
                draw <- engine$c39_draw_a(A, d, d$B, U)
                engine$c39_support(draw, A$members)
                chunks[[k]] <- draw[, cols, drop = FALSE]
            }
            make <- function(raw) list(raw = raw, selected_spearman = data.frame(), physical_pass = TRUE)
            comparison <- methods$nmr_compare_summaries(make(methods$nmrD_table(as.list(d[cols]), w, "SELECTED_ORIGINAL_CONTINUOUS")), 
                make(methods$nmrD_table(as.list(do.call(rbind, chunks)), rep(w/3L, 3L), "SELECTED_ORIGINAL_CONTINUOUS")))
            clinical <- list(schema = "SELECTED_ACTUAL_OBSERVED_B_ENGINEERING_V1", actual_rows = nrow(d), Ksmall = 3L, 
                observed = TRUE, synthetic = FALSE, comparison = comparison, passed = comparison$engineering_pass, 
                official_MC_integrations = 0L, fitted_models = 0L, physical_support_checked_every_generated_full_world = TRUE)
            record <- list(schema = "GS_BASE11_B0_FRESH_A_DISTRIBUTION_VALIDATION_V1", passed = isTRUE(body$passed) && 
                isTRUE(clinical$passed), actual_rows = nrow(d), weight_hash = gs_rt_hash(as.numeric(w)), observed = TRUE, 
                synthetic = FALSE, body = body, clinical = clinical, no_refit_or_law_changes = TRUE, no_old_A_reuse = TRUE)
            artifact <- gs_rt_json(runtime, record, out, "ACTUAL_A_DISTRIBUTION_VALIDATION.json")
            gs_rt_assert(record$passed, "ACTUAL_BASE11_A_DISTRIBUTION_OR_PHYSICAL_SUPPORT_FAILURE")
            c(record, list(artifacts = list(artifact)))
        }
        gs_rt_probe_rows <- function(d) {
            gs_rt_assert(is.data.frame(d) && nrow(d) > 0L && "center" %in% names(d) && !anyNA(d$center), "ACTUAL_CENTER_ROWS_REQUIRED")
            unique(c(match(unique(as.character(d$center)), as.character(d$center)), as.integer(round(seq(1, nrow(d), 
                length.out = min(16L, nrow(d)))))))
        }
        gs_rt_validate_y <- function(runtime, y, design_fn, predictor, d, members, base, domains, w) {
            engine <- runtime$engine
            methods <- runtime$methods
            take <- gs_rt_probe_rows(d)
            C <- d[take, , drop = FALSE]
            cols <- vapply(members, `[[`, "", "column")
            M <- C[, cols, drop = FALSE]
            records <- list()
            max_mass <- 0
            for (b in 0:1) {
                probe <- C
                probe$B <- b
                X <- design_fn(b, M, C)
                original <- engine$c39_design(probe, members, base, attr(design_fn, "reference"))$X
                gs_rt_assert(identical(dimnames(X), dimnames(original)) && identical(dim(X), dim(original)) && 
                  max(abs(X - original)) == 0, "CACHED_ORIGINAL_Y_DESIGN_NOT_EXACT")
                X <- methods$new121_y_apply_alias_plan(X, y$alias_plan)
                check <- predictor$validate(b, M, C, weights = w[take])
                gs_rt_assert(isTRUE(check$validation$within_engineering_reference), "ORIGINAL_ADAPTER_NATIVE_CIF_ENGINEERING_REFERENCE_EXCEEDED")
                for (endpoint in c("IHD", "T2D")) {
                  cf <- y$compiled[[endpoint]]
                  prepared <- engine$fd_fp_prepare(X, cf, engine)
                  lp_error <- 0
                  for (g in unique(prepared$group)) for (j in 1:2) {
                    ii <- which(prepared$group == g)
                    nd <- as.data.frame(prepared$designs[[j]][ii, , drop = FALSE])
                    names(nd) <- cf$object$cause_references[[j]]$safe_names
                    nd$.ce_center <- factor(g, levels = cf$object$center_reference$observed_levels)
                    got <- engine$ff_y_capture(as.numeric(stats::predict(cf$object$fit[[j]], newdata = nd, type = "lp", 
                      reference = "zero")))
                    gs_rt_assert(!length(got$warnings), "NATIVE_LP_VALIDATION_WARNING")
                    fast_lp <- as.numeric(prepared$designs[[j]][ii, , drop = FALSE] %*% cf$beta[[j]])
                    lp_error <- max(lp_error, abs(got$value - fast_lp))
                  }
                  fast <- engine$fd_fp_predict_matrix(cf, X, engine)
                  native <- engine$ff_ce_predict_y_matrix(cf$object, X, 11)
                  gs_rt_assert(!length(attr(native, "warnings")), "NATIVE_CIF_VALIDATION_WARNING")
                  cif_error <- max(abs(fast - native))
                  mass <- max(abs(rowSums(check$all_states[[endpoint]]) - 1))
                  max_mass <- max(max_mass, mass)
                  gs_rt_assert(is.finite(lp_error) && lp_error <= 1e-10 && is.finite(cif_error) && cif_error <= 
                    1e-12 && is.finite(mass) && mass <= 1e-06, "ORIGINAL_FAST_NATIVE_LP_CIF_OR_MASS_VALIDATION_FAILED")
                  rr <- matrix(exp(vapply(1:2, function(j) as.numeric(prepared$designs[[j]][1L, , drop = FALSE] %*% 
                    cf$beta[[j]]), numeric(1))), 1L, 2L)
                  h <- cf$hazards[[prepared$group[1L]]]
                  slow <- engine$c39_num_quadrature(h, rr)
                  slow_error <- max(abs(slow - fast[1L, , drop = FALSE]))
                  gs_rt_assert(is.finite(slow_error) && slow_error <= 1e-10, "ORIGINAL_SLOW_QUADRATURE_REFERENCE_FAILED")
                  records[[length(records) + 1L]] <- list(exposure = b, endpoint = endpoint, actual_rows = nrow(C), 
                    design_max_abs = 0, lp_max_abs = lp_error, cif_all_states_max_abs = cif_error, lp_tolerance = 1e-10, 
                    cif_tolerance = 1e-12, slow_entire_hazard_intervals = nrow(h), slow_actual_rows = 1L, slow_reference = "ORIGINAL_c39_num_quadrature_COMPLETE_ACTUAL_HAZARD_GRID", 
                    slow_max_abs = slow_error, slow_tolerance = 1e-10, adapter = check$validation, max_mass_residual = mass, 
                    passed = TRUE)
                }
            }
            list(schema = "GS_BASE11_B0_ACTUAL_Y_FAST_NATIVE_SLOW_VALIDATION_V1", passed = TRUE, observed = TRUE, 
                synthetic = FALSE, actual_rows = length(take), actual_centers = length(unique(as.character(C$center))), 
                exposures = 0:1, endpoints = c("IHD", "T2D"), maximum_mass_residual = max_mass, genuine_probability_check = TRUE, 
                four_slow_rowpaths = TRUE, complete_baseline_grid = TRUE, no_new_fits = TRUE, consumes_RNG = FALSE, 
                checks = records)
        }
        gs_fit_and_validate_y <- function(runtime, d, members, base, domains, risks, w, out, progress) {
            gs_rt_validate_inputs(runtime, d, members, w)
            gs_rt_assert(is.function(progress) && is.list(risks) && setequal(names(risks), c("IHD", "T2D")) && 
                all(vapply(risks, nrow, integer(1)) == nrow(d)), "EXACT_CURRENT_PAIRED_RISK_FRAMES_AND_PROGRESS_REQUIRED")
            engine <- runtime$engine
            methods <- runtime$methods
            progress("FRESH_ORIGINAL_BASE11_Y_DESIGN")
            design_fn <- methods$new121_y_design_factory(d, members, base, domains, engine, cache_entries = 2L)
            reference <- gs_rt_rds(runtime, attr(design_fn, "reference"), out, "Y_DESIGN_REFERENCE.private.rds")
            owner <- attr(design_fn, "column_owner")
            mediator_cols <- names(owner)[owner != "BASE"]
            audits <- list()
            progress("FIT_FOUR_FRESH_FAMILY_ROBUST_CAUSE_SPECIFIC_Y_BRANCHES")
            y <- methods$new121_fit_y(engine, attr(design_fn, "fullX"), risks, w, d$cluster_id, mediator_cols, 
                centre_actual = d$center, diagnostics = FALSE, audit_callback = function(record) {
                  audits[[length(audits) + 1L]] <<- gs_rt_rds(runtime, record, out, sprintf("Y_NATIVE_AUDIT_%04d.private.rds", 
                    length(audits) + 1L))
                  progress("FRESH_Y_NATIVE_AUDIT", record[intersect(c("stage", "endpoint", "cause", "warnings", 
                    "warning", "status", "reason"), names(record))])
                })
            gs_rt_assert(identical(y$method$tau, 11) && identical(y$method$fit_followup_cap, Inf) && identical(y$method$ties, 
                "efron") && isTRUE(y$method$robust_family_cluster), "ORIGINAL_FULLFOLLOWUP_TAU11_EFRON_FAMILY_METHOD_CHANGED")
            for (ep in c("IHD", "T2D")) {
                fit <- y$fits[[ep]]
                expected_risk <- engine$ff_y_risk(risks[[ep]], 11)
                gs_rt_assert(identical(names(fit$fit), c("disease", "death")) && identical(as.logical(fit$stratified), 
                  c(FALSE, TRUE)) && identical(as.numeric(fit$endpoint$time), as.numeric(expected_risk$time)) && 
                  identical(as.numeric(fit$endpoint$status), as.numeric(expected_risk$status)) && all(vapply(fit$fit, 
                  function(f) !is.null(f$naive.var) && f$n == nrow(d) && identical(f$method, "efron"), logical(1))), 
                  "FRESH_Y_BRANCH_CLUSTER_FULLFOLLOWUP_OR_CAUSE_STRATA_CHANGED")
            }
            yr <- gs_rt_rds(runtime, y, out, "FRESH_Y.private.rds")
            manifest <- list(schema = "GS_BASE11_B0_FRESH_Y_MANIFEST_V1", complete = TRUE, fresh_four_branch_fit = TRUE, 
                old_A_Y_cache_reads = 0L, file = yr, reference = reference, private_audits = audits, source_closure_hash = runtime$source_pin_evidence$closure_hash, 
                model = "S44_BASE11", bootstrap_id = 0L, rows = nrow(d), weights_hash = gs_rt_hash(as.numeric(w)), 
                cluster_hash = gs_rt_hash(d$cluster_id), risk_hash = gs_rt_hash(risks), members_hash = gs_rt_hash(members), 
                base = base, domains = domains, method = y$method, fullX_hash = gs_rt_hash(attr(design_fn, "fullX")))
            mr <- gs_rt_json(runtime, manifest, out, "FRESH_Y_MANIFEST.json")
            predictor <- methods$new121_y_prediction(y, design_fn, engine, domains, cache_entries = 2L)
            progress("VALIDATE_ACTUAL_Y_FAST_NATIVE_AND_SLOW_ORIGINAL_REFERENCES")
            validation <- gs_rt_validate_y(runtime, y, design_fn, predictor, d, members, base, domains, w)
            vr <- gs_rt_json(runtime, validation, out, "FIRST_REAL_RISK_BLOCK.json")
            progress("FRESH_ORIGINAL_SCIENCE_Y_READY", list(actual_validation_rows = validation$actual_rows, four_fresh_branches = TRUE))
            list(y = y, design_fn = design_fn, predictor = predictor, predict_risk = predictor$predict, validation = validation, 
                artifacts = list(fit = yr, reference = reference, manifest = mr, validation = vr, audits = audits))
        }
    }), envir = context)
}

