# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

gs100_worker_task_name <- function() "border_glycaemia_sensitivity_BASE11_B100_20261106_v1"

gs100_worker_runtime_sha <- function() birth_private_setting("PRIVATE_HASH_1")

gs100_worker_mkdir <- function(path, root, exclusive = FALSE) {
    root <- gs100_norm(root)
    gs100_assert(gs100_inside(gs100_norm(path, FALSE), root) && !identical(path, root), "WORKER_DIRECTORY_OUTSIDE_EXACT_ATTEMPT_BOUNDARY")
    if (exclusive) 
        gs100_assert(!file.exists(path), "PRESERVE_EXISTING_WORKER_DIRECTORY")
    if (!dir.exists(path)) 
        gs100_assert(dir.create(path, recursive = TRUE, showWarnings = FALSE, mode = "0700"), "WORKER_DIRECTORY_CREATION_FAILED")
    actual <- gs100_norm(path)
    gs100_assert(gs100_inside(actual, root), "WORKER_DIRECTORY_SYMLINK_ESCAPE_REJECTED")
    actual
}

gs100_worker_file <- function(record, root) {
    gs100_assert(is.list(record) && is.character(record$path) && length(record$path) == 1L && is.character(record$sha256) && 
        length(record$sha256) == 1L, "OWN_CHECKPOINT_FILE_RECORD_REQUIRED")
    path <- gs100_norm(record$path)
    gs100_assert(gs100_inside(path, gs100_norm(root)) && !dir.exists(path) && identical(gs100_sha(path), record$sha256), 
        "OWN_CHECKPOINT_BOUNDARY_OR_BYTES_CHANGED")
    path
}

gs100_worker_rng <- function() list(kind = RNGkind(), present = exists(".Random.seed", .GlobalEnv, inherits = FALSE), 
    seed = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL)

gs100_worker_restore_rng <- function(rng) {
    gs100_assert(is.list(rng) && length(rng$kind) == 3L && is.logical(rng$present) && length(rng$present) == 1L && 
        !is.na(rng$present), "OWN_A_POSTFIT_RNG_RECORD_REQUIRED")
    do.call(RNGkind, as.list(rng$kind))
    if (rng$present) {
        gs100_assert(is.integer(rng$seed) && !anyNA(rng$seed), "OWN_A_POSTFIT_RNG_STATE_INVALID")
        assign(".Random.seed", rng$seed, .GlobalEnv)
    }
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) 
        rm(".Random.seed", envir = .GlobalEnv)
    invisible(NULL)
}

gs100_worker_failure_class <- function(error) {
    scientific <- inherits(error, "ff_ce_error") || grepl("singular|infinite|non.?converg|rank_deficient|absent|unestimable|cox|(^|[^[:alnum:]_])(ff_)?ce_|logistic_failed|continuous_fit_failed|support_empty|draw_illegal|catastrophic_distortion|negative_biochemical_part|physical_failure|probability_mass_failure|nonfinite_risk_multiplier", 
        conditionMessage(error), ignore.case = TRUE)
    if (scientific) 
        "STATISTICAL_FAILURE"
    else "TECHNICAL_FAILURE"
}

gs100_worker_a_check <- function(runtime, A, d, members, base, w) {
    engine <- runtime$engine
    gs100_assert(identical(A$members, members) && identical(names(A$models), names(members)) && identical(A$body_pmm$revision, 
        engine$body_pmm_revision()) && length(A$body_pmm$reused_old_ids) == 0L && is.null(A$body_pmm$reused_old_provenance) && 
        identical(A$body_pmm$prefix_binding, engine$body_tqe_binding(d, members, base, w)), "ALL_ELEVEN_A_MUST_BE_FRESH_AND_BOUND_TO_THIS_SAMPLE_WEIGHTS")
    invisible(TRUE)
}

gs100_worker_y_check <- function(runtime, y, d, risks) {
    gs100_assert(identical(y$method$tau, 11) && identical(y$method$fit_followup_cap, Inf) && identical(y$method$ties, 
        "efron") && isTRUE(y$method$robust_family_cluster), "ORIGINAL_FULLFOLLOWUP_TAU11_EFRON_FAMILY_METHOD_CHANGED")
    for (endpoint in c("IHD", "T2D")) {
        fit <- y$fits[[endpoint]]
        expected <- runtime$engine$ff_y_risk(risks[[endpoint]], 11)
        gs100_assert(identical(names(fit$fit), c("disease", "death")) && identical(as.logical(fit$stratified), 
            c(FALSE, TRUE)) && identical(as.numeric(fit$endpoint$time), as.numeric(expected$time)) && identical(as.numeric(fit$endpoint$status), 
            as.numeric(expected$status)) && all(vapply(fit$fit, function(f) !is.null(f$naive.var) && f$n == nrow(d) && 
            identical(f$method, "efron"), logical(1))), "FRESH_Y_BRANCH_CLUSTER_FULLFOLLOWUP_OR_CAUSE_STRATA_CHANGED")
    }
    invisible(TRUE)
}

gs100_worker_completed_phase <- function(receipts, phase, identity_hash, private) {
    paths <- sort(list.files(receipts, pattern = paste0("^FRESH_", phase, "_[0-9]{4}[.]json$"), full.names = TRUE))
    if (!length(paths)) 
        return(NULL)
    records <- lapply(paths, function(path) {
        gs100_assert(gs100_inside(gs100_norm(path), gs100_norm(receipts)), "OWN_PHASE_MANIFEST_ESCAPE_REJECTED")
        rec <- jsonlite::fromJSON(path, simplifyVector = FALSE)
        gs100_assert(isTRUE(rec$complete) && identical(rec$phase, phase) && identical(rec$identity_hash, identity_hash), 
            "OWN_PHASE_CONFIG_SAMPLE_WEIGHT_SOURCE_BINDING_CHANGED")
        for (file in rec$files) gs100_worker_file(file, private)
        rec
    })
    records[[length(records)]]
}

gs100_worker_kernel_audit <- function(audit) {
    if (!is.null(audit$reverse_failure)) 
        audit$reverse_failure <- list(allocation_order = "reverse", failed_stage = "REVERSE_CUMULATIVE_REPLACEMENT", 
            reason = "REVERSE_STAGE_FAILURE_REVIEW_PRIVATE_DETAIL", disposition = "ALL_REVERSE_ESTIMATES_DROPPED_FORWARD_CONTINUES")
    audit
}

gs100_worker_validate_integral <- function(integral, allocation_order, b, seed) {
    audit <- integral$audit
    orders <- gs100_vec(audit$returned_orders)
    gs100_assert(is.list(audit) && audit$n_mc == 50L && audit$uniform_dim == 11L && audit$block_size == 512L && 
        audit$seed == seed && audit$binding$bootstrap_id == b && identical(orders, if (isTRUE(audit$reverse_failure_recorded)) 
        "forward"
    else c("forward", "reverse")) && identical(names(integral$orders), orders) && identical(gs100_vec(integral$orders$forward), 
        allocation_order) && isTRUE(audit$paired_uniforms_at_exposure_0_1) && isTRUE(audit$independent_domain_streams) && 
        isTRUE(audit$shared_product_endpoints_both_orders), "COMPLETE_PAIRED_WIDTH11_K50_FORWARD_REQUIRED")
    if ("reverse" %in% orders) 
        gs100_assert(identical(gs100_vec(integral$orders$reverse), rev(allocation_order)), "REVERSE_ALLOCATION_ORDER_CHANGED")
    gs100_assert(is.data.frame(integral$components) && is.data.frame(integral$paired) && is.data.frame(integral$scenarios) && 
        all(is.finite(integral$components$effect)) && all(is.finite(integral$paired$Delta)) && all(is.finite(integral$scenarios$risk) & 
        integral$scenarios$risk >= -1e-12 & integral$scenarios$risk <= 1 + 1e-12), "ACTUAL_INTEGRATION_COMPONENTS_OR_RISKS_INVALID")
    for (order in orders) for (endpoint in c("IHD", "T2D")) {
        x <- integral$components[integral$components$allocation_order == order & integral$components$endpoint == 
            endpoint, , drop = FALSE]
        pick <- function(kind) x$effect[x$component_type == kind]
        gs100_assert(nrow(x) == 8L && length(pick("TE")) == 1L && length(pick("DE")) == 1L && length(pick("IE_joint")) == 
            1L && length(pick("DEP")) == 1L && length(pick("IIE")) == 4L && setequal(x$domain_id[x$component_type == 
            "IIE"], allocation_order) && abs(pick("TE") - pick("DE") - pick("IE_joint")) <= 1e-12 && abs(sum(pick("IIE")) + 
            pick("DEP") - pick("IE_joint")) <= 1e-12, "ACTUAL_SAME_VECTOR_DECOMPOSITION_CLOSURE_FAILED")
    }
    list(forward = TRUE, reverse = "reverse" %in% orders)
}

gs100_worker_main <- function(args = commandArgs(TRUE)) {
    if (!(length(args) %in% 2:3) || (length(args) == 3L && !identical(args[[3L]], "--resume"))) 
        stop("USAGE: worker_B100.R /absolute/config.json b [--resume]", call. = FALSE)
    for (p in c("jsonlite", "digest", "survival", "nnet")) if (!requireNamespace(p, quietly = TRUE)) 
        stop(paste("MISSING_REQUIRED_PACKAGE", p), call. = FALSE)
    config_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
    cfg <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
    adapter <- Sys.getenv("GS100_ADAPTER_ROOT", unset = if (is.null(cfg$adapter_root)) 
        ""
    else cfg$adapter_root)
    if (!nzchar(adapter) || !file.exists(file.path(adapter, "R", "prepare_B100.R"))) 
        stop("MISSING_B100_PREPARATION_MODULE", call. = FALSE)
    sys.source(file.path(adapter, "R", "prepare_B100.R"), envir = environment(gs100_worker_main))
    cfg <- gs100_resolve_config(cfg)
    gs100_assert(identical(gs100_norm(adapter), cfg$adapter_root), "ADAPTER_RESOLUTION_MISMATCH")
    b <- suppressWarnings(as.integer(args[[2L]]))
    gs100_assert(length(b) == 1L && !is.na(b) && b >= 1L && b <= 100L && identical(as.character(b), args[[2L]]), 
        "ONLY_EXACT_BOOTSTRAP_IDS_1_THROUGH_100_ALLOWED")
    resume <- length(args) == 3L
    gs100_assert(identical(unname(Sys.info()[["sysname"]]), "Linux") && identical(unname(Sys.info()[["nodename"]]), 
        "DESKTOP-5G16RQU") && identical(as.character(getRversion()), "4.5.3"), "REVIEWED_WINPC_WSL_R453_ONLY")
    gs100_assert(dir.exists(cfg$output_root) && dir.exists(cfg$state_root), "PREPARED_NEW_OUTPUT_AND_STATE_ROOTS_REQUIRED")
    runtime_path <- file.path(cfg$adapter_root, "vendor", "B0_model_runtime.R")
    gs100_assert(file.exists(runtime_path) && identical(gs100_sha(runtime_path), gs100_worker_runtime_sha()), "EXACT_SUCCESSFUL_B0_SOURCE_RUNTIME_BYTES_REQUIRED")
    worker_path <- file.path(cfg$adapter_root, "R", "worker_B100.R")
    prepare_path <- file.path(cfg$adapter_root, "R", "prepare_B100.R")
    seed <- as.integer(20261103L + 1009L * b)
    task_parent <- gs100_worker_mkdir(file.path(cfg$output_root, "tasks"), cfg$output_root)
    task <- gs100_worker_mkdir(file.path(task_parent, sprintf("B%04d", b)), task_parent)
    attempt <- file.path(task, "attempt_0001")
    state_parent <- gs100_worker_mkdir(file.path(cfg$state_root, "workers"), cfg$state_root)
    state_b <- gs100_worker_mkdir(file.path(state_parent, sprintf("B%04d", b)), state_parent)
    state <- file.path(state_b, "attempt_0001")
    if (!resume) {
        attempt <- gs100_worker_mkdir(attempt, task, exclusive = TRUE)
        state <- gs100_worker_mkdir(state, state_b, exclusive = TRUE)
        for (name in c("private", "artifacts", "receipts", "progress", "invocations")) gs100_worker_mkdir(file.path(attempt, 
            name), attempt, exclusive = TRUE)
    }
    else {
        gs100_assert(dir.exists(attempt) && dir.exists(state) && all(dir.exists(file.path(attempt, c("private", 
            "artifacts", "receipts", "progress", "invocations")))), "RESUME_REQUIRES_OWN_EXISTING_COMPLETE_ATTEMPT_LAYOUT")
        attempt <- gs100_norm(attempt)
        state <- gs100_norm(state)
        gs100_assert(gs100_inside(attempt, task) && gs100_inside(state, state_b), "OWN_ATTEMPT_SYMLINK_ESCAPE_REJECTED")
    }
    gs100_assert(!file.exists(file.path(attempt, "RESULT.json")), "TERMINAL_SUCCESS_OR_STATISTICAL_FAILURE_MUST_NOT_BE_RETRIED")
    claim <- list(schema = "GS100_EXCLUSIVE_ATTEMPT_CLAIM_V1", task = gs100_worker_task_name(), b = b, seed = seed, 
        config_sha256 = gs100_sha(config_path), resolved_config_hash = gs100_hash(cfg), code = list(worker = gs100_sha(worker_path), 
            prepare = gs100_sha(prepare_path), runtime = gs100_sha(runtime_path)), attempt_root = attempt, state_root = state)
    claim_path <- file.path(attempt, "receipts", "ATTEMPT_CLAIM.json")
    if (resume) {
        gs100_assert(file.exists(claim_path) && identical(gs100_hash(jsonlite::fromJSON(claim_path, simplifyVector = FALSE)), 
            gs100_hash(jsonlite::fromJSON(jsonlite::toJSON(claim, auto_unbox = TRUE, null = "null"), simplifyVector = FALSE))), 
            "OWN_ATTEMPT_CONFIG_CODE_PATH_OR_SEED_CHANGED")
    }
    else gs100_new_json(claim, claim_path)
    existing <- list.files(file.path(attempt, "invocations"), pattern = "^launch_[0-9]{4}$")
    launch_number <- if (length(existing)) 
        max(as.integer(sub("launch_", "", existing, fixed = TRUE))) + 1L
    else 1L
    gs100_assert(launch_number <= 9999L && (resume || launch_number == 1L), "OWN_INVOCATION_NUMBER_INVALID")
    launch_id <- sprintf("launch_%04d", launch_number)
    invocation <- gs100_worker_mkdir(file.path(attempt, "invocations", launch_id), attempt, exclusive = TRUE)
    private <- gs100_worker_mkdir(file.path(attempt, "private", launch_id), attempt, exclusive = TRUE)
    artifacts <- gs100_worker_mkdir(file.path(attempt, "artifacts", launch_id), attempt, exclusive = TRUE)
    receipts <- file.path(attempt, "receipts")
    stage <- "CLAIMED"
    event <- 0L
    last_integrating <- -Inf
    warnings <- list()
    progress <- function(new_stage, details = list()) {
        gs100_assert(is.character(new_stage) && length(new_stage) == 1L && grepl("^[A-Za-z][A-Za-z0-9_]*$", new_stage), 
            "SAFE_AGGREGATE_STAGE_NAME_REQUIRED")
        stage <<- new_stage
        if (identical(new_stage, "integrating") && proc.time()[["elapsed"]] - last_integrating < 20) 
            return(invisible(NULL))
        if (identical(new_stage, "integrating")) 
            last_integrating <<- proc.time()[["elapsed"]]
        allow <- c("N_target", "N_fit", "conditional_index", "total", "mc", "n_mc", "rows_complete", "n_target", 
            "endpoint", "cause", "status", "b", "fresh_A", "fresh_Y", "passed", "reverse_failure_recorded")
        details <- details[intersect(names(details), allow)]
        event <<- event + 1L
        record <- list(schema = "GLYCAEMIA_BASE11_B100_PROGRESS_V1", utc = gs100_utc(), pid = Sys.getpid(), b = b, 
            launch = launch_id, stage = new_stage, details = details)
        gs100_new_json(record, file.path(attempt, "progress", sprintf("%s_%05d.json", launch_id, event)))
        cat(jsonlite::toJSON(record, auto_unbox = TRUE), "\n")
        flush.console()
        invisible(NULL)
    }
    started <- proc.time()[["elapsed"]]
    old_options <- options(body_tqe.upstream_A = NULL, body_tqe.prefix_A = NULL, body_pmm.old_A11 = NULL, body_pmm.prefix_A = NULL)
    on.exit(options(old_options), add = TRUE)
    result <- tryCatch(withCallingHandlers({
        gs100_new_json(list(utc = gs100_utc(), pid = Sys.getpid(), b = b, seed = seed, resume = resume, claim = claim, 
            roots = cfg[c("input_path", "adapter_root", "source_root", "y_root", "state_root", "output_root")], 
            old_fitted_A_Y_or_prediction_cache_reads = 0L, automatic_retry = FALSE, participant_rows_in_reports = 0L), 
            file.path(invocation, "INVOCATION.json"))
        progress("LOADING_OWN_FRESH_DRAW_AND_WEIGHTS")
        prepared <- gs100_load_draw(cfg, b)
        z <- prepared$z
        d <- prepared$d
        w <- prepared$w
        risks <- prepared$risks
        gs100_assert(prepared$manifest$b == b && prepared$manifest$seed == seed, "EXACT_CURRENT_BOOTSTRAP_NUMBER_AND_ORIGINAL_SEED_REQUIRED")
        rt <- new.env(parent = environment())
        sys.source(runtime_path, envir = rt)
        rt$gs_rt_task_name <- gs100_worker_task_name
        progress("LOADING_PINNED_SCIENTIFIC_RUNTIME")
        runtime <- rt$gs_load_runtime(cfg$source_root, cfg$y_root, state, attempt)
        gs100_new_json(runtime$source_pin_evidence, file.path(invocation, "SCIENTIFIC_SOURCE_PINS.json"))
        engine <- runtime$engine
        methods <- runtime$methods
        members <- z$members
        base <- gs100_vec(z$base)
        model_spec <- z$model_specifications$S44_BASE11
        allocation_order <- gs100_vec(model_spec$allocation_order)
        domains <- lapply(z$domains[allocation_order], gs100_vec)
        model_cols <- vapply(members, `[[`, "", "column")
        gs100_assert(length(members) == 11L && length(domains) == 4L && identical(names(members), gs100_vec(model_spec$members)) && 
            identical(allocation_order, c("SMOKING", "BODY_COMPOSITION", "CLINICAL_LIPIDS", "GLYCAEMIA")) && setequal(unlist(domains, 
            use.names = FALSE), model_cols) && !anyDuplicated(unlist(domains, use.names = FALSE)) && !isTRUE(model_spec$total_only), 
            "EXACT_BASE11_FOUR_DOMAIN_SCIENTIFIC_CONTRACT_CHANGED")
        rt$gs_rt_validate_inputs(runtime, d, members, w)
        rng_spec <- list(kind = "L'Ecuyer-CMRG", normal = "Inversion", sample = "Rejection", mapping = "ROW_MAJOR_WIDTH11_COMMON_JOINT", 
            independent_domain_streams = TRUE, seed = seed)
        binding <- list(panel_id = "S44_BASE11_HBA1C_LT48_SENSITIVITY_B100", bootstrap_id = b, sample_hash = prepared$manifest$fit_sample_hash, 
            population_hash = z$target_hash, model_spec_hash = gs100_hash(model_spec), input_hash = prepared$files$input$sha256, 
            rng_spec_hash = gs100_hash(rng_spec))
        identity <- list(schema = "GLYCAEMIA_BASE11_B100_OWN_ATTEMPT_IDENTITY_V1", claim = claim, binding = binding, 
            rng_spec = rng_spec, draw_files = prepared$files, draw_manifest_sha256 = gs100_sha(file.path(gs100_draw_folder(cfg, 
                b), "DRAW_MANIFEST.json")), sampling_hash = gs100_hash(as.integer(prepared$index$indices)), sample_hash = gs100_hash(as.character(d$eid)), 
            data_hash = gs100_hash(d), risk_hash = gs100_hash(risks), weight_hash = gs100_hash(as.numeric(w)), 
            pair_binding_hash = gs100_hash(prepared$selection$pair_binding), cluster_hash = gs100_hash(d$cluster_id), 
            selection_source_hash = z$selection_source_hash, members_hash = gs100_hash(members), base = base, domains = domains, 
            source_closure_hash = runtime$source_pin_evidence$closure_hash, source_pins_hash = gs100_hash(rt$gs_rt_pins()), 
            R = R.version.string, platform = R.version$platform, packages = setNames(as.list(vapply(c("jsonlite", 
                "digest", "survival", "splines", "nnet"), function(p) as.character(utils::packageVersion(p)), "")), 
                c("jsonlite", "digest", "survival", "splines", "nnet")))
        identity_hash <- gs100_hash(identity)
        identity_path <- file.path(receipts, "MODEL_BINDING.private.rds")
        if (resume && file.exists(identity_path)) {
            previous <- readRDS(identity_path)
            gs100_assert(identical(previous$identity_hash, identity_hash) && identical(previous$identity, identity), 
                "OWN_RESUME_SOURCE_CONFIG_SAMPLE_RISK_WEIGHT_BINDING_CHANGED")
        }
        else {
            gs100_assert(!file.exists(identity_path), "PRESERVE_EXISTING_MODEL_BINDING")
            gs100_new_rds(list(identity = identity, identity_hash = identity_hash), identity_path)
        }
        gs100_new_json(list(binding = binding, identity_hash = identity_hash, rng_spec = rng_spec, model_spec = model_spec, 
            N_target = prepared$manifest$N_target, N_fit = nrow(d), new_draw = prepared$files, no_old_A_Y_or_prediction_reuse = TRUE), 
            file.path(invocation, "MODEL_BINDING.json"))
        RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
        set.seed(seed)
        A_record <- if (resume) 
            gs100_worker_completed_phase(receipts, "A", identity_hash, file.path(attempt, "private"))
        else NULL
        reused_A <- !is.null(A_record)
        if (is.null(A_record)) {
            gs100_assert(is.function(engine$body_pmm_fit_a) && identical(engine$c39_fit_a, engine$body_pmm_fit_a) && 
                all(c("prefix_fit", "old_fit") %in% names(formals(engine$body_pmm_fit_a))), "EXACT_PMM_FRESH_ALL_A_API_REQUIRED")
            cb_options <- options(c39.a_checkpoint_callback = function(id, fit, index) progress("FRESH_A_CONDITIONAL", 
                list(conditional_index = index, total = 11L)))
            progress("FITTING_ALL_ELEVEN_FRESH_A")
            A <- tryCatch(engine$body_pmm_fit_a(d, members, base, w, prefix_fit = NULL, old_fit = NULL), finally = options(cb_options))
            gs100_worker_a_check(runtime, A, d, members, base, w)
            a_file <- gs100_new_rds(A, file.path(private, "FRESH_A.private.rds"))
            rng_file <- gs100_new_rds(gs100_worker_rng(), file.path(private, "FRESH_A_POSTFIT_RNG.private.rds"))
            A_record <- list(phase = "A", complete = TRUE, identity_hash = identity_hash, fresh_all_eleven = TRUE, 
                old_or_prefix_reuse = FALSE, files = list(A = a_file, rng = rng_file))
            gs100_new_json(A_record, file.path(receipts, sprintf("FRESH_A_%04d.json", launch_number)))
        }
        else {
            progress("RESTORING_OWN_EXACT_SAME_ATTEMPT_A")
            A <- readRDS(gs100_worker_file(A_record$files$A, file.path(attempt, "private")))
            gs100_worker_restore_rng(readRDS(gs100_worker_file(A_record$files$rng, file.path(attempt, "private"))))
        }
        if (reused_A) 
            gs100_worker_a_check(runtime, A, d, members, base, w)
        progress("BUILDING_ORIGINAL_BASE11_Y_DESIGN")
        design_fn <- methods$new121_y_design_factory(d, members, base, domains, engine, cache_entries = 2L)
        fullX_hash <- gs100_hash(attr(design_fn, "fullX"))
        Y_record <- if (resume) 
            gs100_worker_completed_phase(receipts, "Y", identity_hash, file.path(attempt, "private"))
        else NULL
        reused_Y <- !is.null(Y_record)
        if (is.null(Y_record)) {
            reference <- gs100_new_rds(attr(design_fn, "reference"), file.path(private, "Y_DESIGN_REFERENCE.private.rds"))
            owner <- attr(design_fn, "column_owner")
            mediator_cols <- names(owner)[owner != "BASE"]
            audits <- list()
            progress("FITTING_FOUR_FRESH_FAMILY_ROBUST_CAUSE_SPECIFIC_Y")
            y <- methods$new121_fit_y(engine, attr(design_fn, "fullX"), risks, w, d$cluster_id, mediator_cols, 
                centre_actual = d$center, diagnostics = FALSE, audit_callback = function(record) {
                  audits[[length(audits) + 1L]] <<- gs100_new_rds(record, file.path(private, sprintf("Y_NATIVE_AUDIT_%04d.private.rds", 
                    length(audits) + 1L)))
                  progress("FRESH_Y_NATIVE_AUDIT", record[intersect(c("endpoint", "cause", "status"), names(record))])
                })
            gs100_worker_y_check(runtime, y, d, risks)
            y_file <- gs100_new_rds(y, file.path(private, "FRESH_Y.private.rds"))
            Y_record <- list(phase = "Y", complete = TRUE, identity_hash = identity_hash, A_sha256 = A_record$files$A$sha256, 
                fullX_hash = fullX_hash, fresh_four_branches = TRUE, full_followup = TRUE, risk_horizon = 11L, 
                files = c(list(Y = y_file, reference = reference), setNames(audits, sprintf("audit_%04d", seq_along(audits)))))
            gs100_new_json(Y_record, file.path(receipts, sprintf("FRESH_Y_%04d.json", launch_number)))
        }
        else {
            gs100_assert(identical(Y_record$A_sha256, A_record$files$A$sha256) && identical(Y_record$fullX_hash, 
                fullX_hash), "OWN_Y_A_OR_ACTUAL_DESIGN_BINDING_CHANGED")
            progress("RESTORING_OWN_EXACT_SAME_ATTEMPT_Y")
            y <- readRDS(gs100_worker_file(Y_record$files$Y, file.path(attempt, "private")))
            gs100_worker_y_check(runtime, y, d, risks)
        }
        predictor <- methods$new121_y_prediction(y, design_fn, engine, domains, cache_entries = 2L)
        validation <- NULL
        predict_risk <- function(b, M, C) {
            if (is.null(validation)) {
                take <- seq_len(min(16L, nrow(C)))
                check <- predictor$validate(b, M[take, , drop = FALSE], C[take, , drop = FALSE])
                mass <- vapply(check$all_states, function(p) max(abs(rowSums(p) - 1)), 0)
                gs100_assert(all(is.finite(mass) & mass <= 1e-06), "REAL_FIRST_BLOCK_PROBABILITY_MASS_FAILURE")
                validation <<- c(check$validation, list(max_mass_residual = max(mass), genuine_probability_check = TRUE, 
                  actual_rows = length(take), endpoints = c("IHD", "T2D"), passed = TRUE, source = "ORIGINAL_B_GT0_FIRST_REAL_16_ROW_GUARD", 
                  no_standalone_numeric_probe = TRUE))
                gs100_new_json(validation, file.path(artifacts, "FIRST_REAL_RISK_BLOCK.json"))
            }
            predictor$predict(b, M, C)
        }
        progress("ORIGINAL_FIRST_REAL_RISK_BLOCK_GUARD")
        take <- seq_len(min(16L, nrow(d)))
        invisible(predict_risk(0, d[take, model_cols, drop = FALSE], d[take, , drop = FALSE]))
        gs100_new_json(list(passed = TRUE, b = b, identity_hash = identity_hash, fresh_all_A = TRUE, fresh_four_Y = TRUE, 
            validation = validation, K50_started = FALSE, no_B0_K3_distribution_diagnostic = TRUE, no_slow_quadrature_probe = TRUE), 
            file.path(invocation, "PREINTEGRATION_GATE_PASSED.json"))
        draw_joint <- function(b, C, uniforms, requested_domain = NULL) engine$c39_draw_a(A, C, b, uniforms[, seq_along(members), 
            drop = FALSE])[, model_cols, drop = FALSE]
        support <- function(M, C, context) {
            requested <- Filter(function(m) m$column %in% names(M), members)
            gs100_assert(length(requested) > 0L, "EMPTY_REQUESTED_MEDIATOR_SUPPORT")
            engine$c39_support(M, requested)
            TRUE
        }
        checkpoint <- NULL
        cp_paths <- list.files(receipts, pattern = "^checkpoint_MC[0-9]{2}_launch_[0-9]{4}[.]json$", full.names = TRUE)
        if (resume && length(cp_paths)) {
            cp_records <- lapply(cp_paths, function(path) {
                gs100_assert(gs100_inside(gs100_norm(path), gs100_norm(receipts)), "OWN_MC_MANIFEST_ESCAPE_REJECTED")
                rec <- jsonlite::fromJSON(path, simplifyVector = FALSE)
                gs100_assert(isTRUE(rec$complete) && identical(rec$identity_hash, identity_hash) && identical(rec$A_sha256, 
                  A_record$files$A$sha256) && identical(rec$Y_sha256, Y_record$files$Y$sha256) && is.numeric(rec$completed_mc) && 
                  length(rec$completed_mc) == 1L && rec$completed_mc >= 1L && rec$completed_mc <= 50L && rec$completed_mc == 
                  floor(rec$completed_mc), "OWN_MC_CONFIG_A_Y_SAMPLE_WEIGHT_SOURCE_BINDING_CHANGED")
                gs100_worker_file(rec$file, file.path(attempt, "private"))
                rec
            })
            rec <- cp_records[[which.max(vapply(cp_records, function(x) as.integer(x$completed_mc), 1L))]]
            checkpoint <- readRDS(gs100_worker_file(rec$file, file.path(attempt, "private")))
            gs100_assert(identical(as.integer(checkpoint$completed_mc), as.integer(rec$completed_mc)), "OWN_MC_MANIFEST_PAYLOAD_PASS_ACCOUNTING_CHANGED")
            progress("RESUMING_OWN_COMMITTED_MC_PASS", list(mc = rec$completed_mc, n_mc = 50L))
        }
        save_checkpoint <- function(state) {
            k <- as.integer(state$completed_mc)
            gs100_assert(length(k) == 1L && !is.na(k) && k >= 1L && k <= 50L, "MC_CHECKPOINT_NUMBER_INVALID")
            file <- gs100_new_rds(state, file.path(private, sprintf("MC_%02d.private.rds", k)))
            gs100_new_json(list(complete = TRUE, b = b, completed_mc = k, identity_hash = identity_hash, A_sha256 = A_record$files$A$sha256, 
                Y_sha256 = Y_record$files$Y$sha256, file = file), file.path(receipts, sprintf("checkpoint_MC%02d_%s.json", 
                k, launch_id)))
        }
        progress("INTEGRATING_K50_FORWARD_REVERSE_DEP")
        integral <- methods$nmr_product_stream(d, w, domains, draw_joint, predict_risk, support, uniform_dim = 11L, 
            n_mc = 50L, block_size = 512L, seed = seed, allocation_order = allocation_order, include_reverse = TRUE, 
            endpoints = c("IHD", "T2D"), binding = binding, checkpoint = checkpoint, save_checkpoint = save_checkpoint, 
            progress = function(v) progress(v$stage, v[setdiff(names(v), "stage")]))
        order_success <- gs100_worker_validate_integral(integral, allocation_order, b, seed)
        for (name in c("components", "paired", "scenarios")) {
            path <- file.path(artifacts, paste0(toupper(name), ".csv"))
            gs100_assert(!file.exists(path), "PRESERVE_EXISTING_AGGREGATE_RESULT_TABLE")
            utils::write.csv(integral[[name]], path, row.names = FALSE)
        }
        integral_file <- gs100_new_rds(integral, file.path(private, "INTEGRATION.private.rds"))
        progress("REPLICATE_COMPLETE", list(fresh_A = TRUE, fresh_Y = TRUE, passed = TRUE))
        list(schema = "GLYCAEMIA_SENSITIVITY_BASE11_B100_REPLICATE_V1", status = "SUCCESS", utc = gs100_utc(), 
            b = b, model_id = "S44_BASE11", seed = seed, binding = binding, identity_hash = identity_hash, N = prepared$manifest[c("N_target", 
                "N_frozen_target", "N_fit", "N_standardization", "endpoints", "IPW")], components = integral$components, 
            paired = integral$paired, scenarios = integral$scenarios, stage = "REPLICATE_COMPLETE", reason = "COMPLETE_FORWARD_K50", 
            safe_code = "SUCCESS", order_success = order_success, audit = list(kernel = gs100_worker_kernel_audit(integral$audit), 
                fresh_all_A = TRUE, fresh_four_Y = TRUE, original_seed = seed, rng_spec = rng_spec, first_real_risk_block = validation, 
                no_B0_K3_distribution_diagnostic = TRUE, no_slow_quadrature_probe = TRUE, old_A_Y_or_prediction_cache_reuse = FALSE, 
                own_same_attempt_A_resume = reused_A, own_same_attempt_Y_resume = reused_Y, native_prediction = predictor$stats(), 
                selection_refits = 1L, bootstrap_sampling_performed = TRUE, full_followup_preserved = TRUE, report_horizon_years = 11L), 
            artifacts = list(A = A_record$files$A, Y = Y_record$files$Y, integration = integral_file), participant_rows_in_reports = 0L)
    }, warning = function(warning) {
        warnings[[length(warnings) + 1L]] <<- list(stage = stage, message = conditionMessage(warning))
        invokeRestart("muffleWarning")
    }), error = function(error) {
        gs100_new_rds(list(stage = stage, class = class(error), message = conditionMessage(error)), file.path(private, 
            "FAILURE_DETAIL.private.rds"))
        status <- gs100_worker_failure_class(error)
        list(schema = "GLYCAEMIA_SENSITIVITY_BASE11_B100_REPLICATE_V1", status = status, utc = gs100_utc(), b = b, 
            model_id = "S44_BASE11", seed = seed, stage = stage, reason = if (identical(status, "STATISTICAL_FAILURE")) "BOOTSTRAP_STATISTICAL_STAGE_FAILURE" else "WORKER_TECHNICAL_STAGE_FAILURE", 
            safe_code = paste0(status, "_REVIEW_PRIVATE_DETAIL"), automatic_retry = FALSE, resume_possible = identical(status, 
                "TECHNICAL_FAILURE"), participant_rows_in_reports = 0L)
    })
    gs100_new_rds(warnings, file.path(private, "WARNINGS.private.rds"))
    result$warning_count <- length(warnings)
    result$elapsed_seconds <- proc.time()[["elapsed"]] - started
    result$launch <- launch_id
    result$resume_requested <- resume
    gs100_new_json(result, file.path(invocation, "RESULT.json"))
    if (result$status %in% c("SUCCESS", "STATISTICAL_FAILURE")) 
        gs100_new_json(result, file.path(attempt, "RESULT.json"))
    if (identical(result$status, "SUCCESS")) 
        0L
    else if (identical(result$status, "STATISTICAL_FAILURE")) 
        2L
    else 1L
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        gs100_worker_task_name <- function() "border_glycaemia_sensitivity_BASE11_B100_20261106_v1"
        gs100_worker_runtime_sha <- function() birth_private_setting("PRIVATE_HASH_2")
        gs100_worker_mkdir <- function(path, root, exclusive = FALSE) {
            root <- gs100_norm(root)
            gs100_assert(gs100_inside(gs100_norm(path, FALSE), root) && !identical(path, root), "WORKER_DIRECTORY_OUTSIDE_EXACT_ATTEMPT_BOUNDARY")
            if (exclusive) gs100_assert(!file.exists(path), "PRESERVE_EXISTING_WORKER_DIRECTORY")
            if (!dir.exists(path)) gs100_assert(dir.create(path, recursive = TRUE, showWarnings = FALSE, mode = "0700"), 
                "WORKER_DIRECTORY_CREATION_FAILED")
            actual <- gs100_norm(path)
            gs100_assert(gs100_inside(actual, root), "WORKER_DIRECTORY_SYMLINK_ESCAPE_REJECTED")
            actual
        }
        gs100_worker_file <- function(record, root) {
            gs100_assert(is.list(record) && is.character(record$path) && length(record$path) == 1L && is.character(record$sha256) && 
                length(record$sha256) == 1L, "OWN_CHECKPOINT_FILE_RECORD_REQUIRED")
            path <- gs100_norm(record$path)
            gs100_assert(gs100_inside(path, gs100_norm(root)) && !dir.exists(path) && identical(gs100_sha(path), 
                record$sha256), "OWN_CHECKPOINT_BOUNDARY_OR_BYTES_CHANGED")
            path
        }
        gs100_worker_rng <- function() list(kind = RNGkind(), present = exists(".Random.seed", .GlobalEnv, inherits = FALSE), 
            seed = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL)
        gs100_worker_restore_rng <- function(rng) {
            gs100_assert(is.list(rng) && length(rng$kind) == 3L && is.logical(rng$present) && length(rng$present) == 
                1L && !is.na(rng$present), "OWN_A_POSTFIT_RNG_RECORD_REQUIRED")
            do.call(RNGkind, as.list(rng$kind))
            if (rng$present) {
                gs100_assert(is.integer(rng$seed) && !anyNA(rng$seed), "OWN_A_POSTFIT_RNG_STATE_INVALID")
                assign(".Random.seed", rng$seed, .GlobalEnv)
            } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            invisible(NULL)
        }
        gs100_worker_failure_class <- function(error) {
            scientific <- inherits(error, "ff_ce_error") || grepl("singular|infinite|non.?converg|rank_deficient|absent|unestimable|cox|(^|[^[:alnum:]_])(ff_)?ce_|logistic_failed|continuous_fit_failed|support_empty|draw_illegal|catastrophic_distortion|negative_biochemical_part|physical_failure|probability_mass_failure|nonfinite_risk_multiplier", 
                conditionMessage(error), ignore.case = TRUE)
            if (scientific) "STATISTICAL_FAILURE" else "TECHNICAL_FAILURE"
        }
        gs100_worker_a_check <- function(runtime, A, d, members, base, w) {
            engine <- runtime$engine
            gs100_assert(identical(A$members, members) && identical(names(A$models), names(members)) && identical(A$body_pmm$revision, 
                engine$body_pmm_revision()) && length(A$body_pmm$reused_old_ids) == 0L && is.null(A$body_pmm$reused_old_provenance) && 
                identical(A$body_pmm$prefix_binding, engine$body_tqe_binding(d, members, base, w)), "ALL_ELEVEN_A_MUST_BE_FRESH_AND_BOUND_TO_THIS_SAMPLE_WEIGHTS")
            invisible(TRUE)
        }
        gs100_worker_y_check <- function(runtime, y, d, risks) {
            gs100_assert(identical(y$method$tau, 11) && identical(y$method$fit_followup_cap, Inf) && identical(y$method$ties, 
                "efron") && isTRUE(y$method$robust_family_cluster), "ORIGINAL_FULLFOLLOWUP_TAU11_EFRON_FAMILY_METHOD_CHANGED")
            for (endpoint in c("IHD", "T2D")) {
                fit <- y$fits[[endpoint]]
                expected <- runtime$engine$ff_y_risk(risks[[endpoint]], 11)
                gs100_assert(identical(names(fit$fit), c("disease", "death")) && identical(as.logical(fit$stratified), 
                  c(FALSE, TRUE)) && identical(as.numeric(fit$endpoint$time), as.numeric(expected$time)) && identical(as.numeric(fit$endpoint$status), 
                  as.numeric(expected$status)) && all(vapply(fit$fit, function(f) !is.null(f$naive.var) && f$n == 
                  nrow(d) && identical(f$method, "efron"), logical(1))), "FRESH_Y_BRANCH_CLUSTER_FULLFOLLOWUP_OR_CAUSE_STRATA_CHANGED")
            }
            invisible(TRUE)
        }
        gs100_worker_completed_phase <- function(receipts, phase, identity_hash, private) {
            paths <- sort(list.files(receipts, pattern = paste0("^FRESH_", phase, "_[0-9]{4}[.]json$"), full.names = TRUE))
            if (!length(paths)) return(NULL)
            records <- lapply(paths, function(path) {
                gs100_assert(gs100_inside(gs100_norm(path), gs100_norm(receipts)), "OWN_PHASE_MANIFEST_ESCAPE_REJECTED")
                rec <- jsonlite::fromJSON(path, simplifyVector = FALSE)
                gs100_assert(isTRUE(rec$complete) && identical(rec$phase, phase) && identical(rec$identity_hash, 
                  identity_hash), "OWN_PHASE_CONFIG_SAMPLE_WEIGHT_SOURCE_BINDING_CHANGED")
                for (file in rec$files) gs100_worker_file(file, private)
                rec
            })
            records[[length(records)]]
        }
        gs100_worker_kernel_audit <- function(audit) {
            if (!is.null(audit$reverse_failure)) audit$reverse_failure <- list(allocation_order = "reverse", failed_stage = "REVERSE_CUMULATIVE_REPLACEMENT", 
                reason = "REVERSE_STAGE_FAILURE_REVIEW_PRIVATE_DETAIL", disposition = "ALL_REVERSE_ESTIMATES_DROPPED_FORWARD_CONTINUES")
            audit
        }
        gs100_worker_validate_integral <- function(integral, allocation_order, b, seed) {
            audit <- integral$audit
            orders <- gs100_vec(audit$returned_orders)
            gs100_assert(is.list(audit) && audit$n_mc == 50L && audit$uniform_dim == 11L && audit$block_size == 
                512L && audit$seed == seed && audit$binding$bootstrap_id == b && identical(orders, if (isTRUE(audit$reverse_failure_recorded)) "forward" else c("forward", 
                "reverse")) && identical(names(integral$orders), orders) && identical(gs100_vec(integral$orders$forward), 
                allocation_order) && isTRUE(audit$paired_uniforms_at_exposure_0_1) && isTRUE(audit$independent_domain_streams) && 
                isTRUE(audit$shared_product_endpoints_both_orders), "COMPLETE_PAIRED_WIDTH11_K50_FORWARD_REQUIRED")
            if ("reverse" %in% orders) gs100_assert(identical(gs100_vec(integral$orders$reverse), rev(allocation_order)), 
                "REVERSE_ALLOCATION_ORDER_CHANGED")
            gs100_assert(is.data.frame(integral$components) && is.data.frame(integral$paired) && is.data.frame(integral$scenarios) && 
                all(is.finite(integral$components$effect)) && all(is.finite(integral$paired$Delta)) && all(is.finite(integral$scenarios$risk) & 
                integral$scenarios$risk >= -1e-12 & integral$scenarios$risk <= 1 + 1e-12), "ACTUAL_INTEGRATION_COMPONENTS_OR_RISKS_INVALID")
            for (order in orders) for (endpoint in c("IHD", "T2D")) {
                x <- integral$components[integral$components$allocation_order == order & integral$components$endpoint == 
                  endpoint, , drop = FALSE]
                pick <- function(kind) x$effect[x$component_type == kind]
                gs100_assert(nrow(x) == 8L && length(pick("TE")) == 1L && length(pick("DE")) == 1L && length(pick("IE_joint")) == 
                  1L && length(pick("DEP")) == 1L && length(pick("IIE")) == 4L && setequal(x$domain_id[x$component_type == 
                  "IIE"], allocation_order) && abs(pick("TE") - pick("DE") - pick("IE_joint")) <= 1e-12 && abs(sum(pick("IIE")) + 
                  pick("DEP") - pick("IE_joint")) <= 1e-12, "ACTUAL_SAME_VECTOR_DECOMPOSITION_CLOSURE_FAILED")
            }
            list(forward = TRUE, reverse = "reverse" %in% orders)
        }
        gs100_worker_main <- function(args = commandArgs(TRUE)) {
            if (!(length(args) %in% 2:3) || (length(args) == 3L && !identical(args[[3L]], "--resume"))) stop("USAGE: worker_B100.R /absolute/config.json b [--resume]", 
                call. = FALSE)
            for (p in c("jsonlite", "digest", "survival", "nnet")) if (!requireNamespace(p, quietly = TRUE)) stop(paste("MISSING_REQUIRED_PACKAGE", 
                p), call. = FALSE)
            config_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
            cfg <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
            adapter <- Sys.getenv("GS100_ADAPTER_ROOT", unset = if (is.null(cfg$adapter_root)) "" else cfg$adapter_root)
            if (!nzchar(adapter) || !file.exists(file.path(adapter, "R", "prepare_B100.R"))) stop("MISSING_B100_PREPARATION_MODULE", 
                call. = FALSE)
            sys.source(file.path(adapter, "R", "prepare_B100.R"), envir = environment(gs100_worker_main))
            cfg <- gs100_resolve_config(cfg)
            gs100_assert(identical(gs100_norm(adapter), cfg$adapter_root), "ADAPTER_RESOLUTION_MISMATCH")
            b <- suppressWarnings(as.integer(args[[2L]]))
            gs100_assert(length(b) == 1L && !is.na(b) && b >= 1L && b <= 100L && identical(as.character(b), args[[2L]]), 
                "ONLY_EXACT_BOOTSTRAP_IDS_1_THROUGH_100_ALLOWED")
            resume <- length(args) == 3L
            gs100_assert(identical(unname(Sys.info()[["sysname"]]), "Linux") && identical(unname(Sys.info()[["nodename"]]), 
                "DESKTOP-5G16RQU") && identical(as.character(getRversion()), "4.5.3"), "REVIEWED_WINPC_WSL_R453_ONLY")
            gs100_assert(dir.exists(cfg$output_root) && dir.exists(cfg$state_root), "PREPARED_NEW_OUTPUT_AND_STATE_ROOTS_REQUIRED")
            runtime_path <- file.path(cfg$adapter_root, "vendor", "B0_model_runtime.R")
            gs100_assert(file.exists(runtime_path) && identical(gs100_sha(runtime_path), gs100_worker_runtime_sha()), 
                "EXACT_SUCCESSFUL_B0_SOURCE_RUNTIME_BYTES_REQUIRED")
            worker_path <- file.path(cfg$adapter_root, "R", "worker_B100.R")
            prepare_path <- file.path(cfg$adapter_root, "R", "prepare_B100.R")
            seed <- as.integer(20261103L + 1009L * b)
            task_parent <- gs100_worker_mkdir(file.path(cfg$output_root, "tasks"), cfg$output_root)
            task <- gs100_worker_mkdir(file.path(task_parent, sprintf("B%04d", b)), task_parent)
            attempt <- file.path(task, "attempt_0001")
            state_parent <- gs100_worker_mkdir(file.path(cfg$state_root, "workers"), cfg$state_root)
            state_b <- gs100_worker_mkdir(file.path(state_parent, sprintf("B%04d", b)), state_parent)
            state <- file.path(state_b, "attempt_0001")
            if (!resume) {
                attempt <- gs100_worker_mkdir(attempt, task, exclusive = TRUE)
                state <- gs100_worker_mkdir(state, state_b, exclusive = TRUE)
                for (name in c("private", "artifacts", "receipts", "progress", "invocations")) gs100_worker_mkdir(file.path(attempt, 
                  name), attempt, exclusive = TRUE)
            } else {
                gs100_assert(dir.exists(attempt) && dir.exists(state) && all(dir.exists(file.path(attempt, c("private", 
                  "artifacts", "receipts", "progress", "invocations")))), "RESUME_REQUIRES_OWN_EXISTING_COMPLETE_ATTEMPT_LAYOUT")
                attempt <- gs100_norm(attempt)
                state <- gs100_norm(state)
                gs100_assert(gs100_inside(attempt, task) && gs100_inside(state, state_b), "OWN_ATTEMPT_SYMLINK_ESCAPE_REJECTED")
            }
            gs100_assert(!file.exists(file.path(attempt, "RESULT.json")), "TERMINAL_SUCCESS_OR_STATISTICAL_FAILURE_MUST_NOT_BE_RETRIED")
            claim <- list(schema = "GS100_EXCLUSIVE_ATTEMPT_CLAIM_V1", task = gs100_worker_task_name(), b = b, 
                seed = seed, config_sha256 = gs100_sha(config_path), resolved_config_hash = gs100_hash(cfg), code = list(worker = gs100_sha(worker_path), 
                  prepare = gs100_sha(prepare_path), runtime = gs100_sha(runtime_path)), attempt_root = attempt, 
                state_root = state)
            claim_path <- file.path(attempt, "receipts", "ATTEMPT_CLAIM.json")
            if (resume) {
                gs100_assert(file.exists(claim_path) && identical(gs100_hash(jsonlite::fromJSON(claim_path, simplifyVector = FALSE)), 
                  gs100_hash(jsonlite::fromJSON(jsonlite::toJSON(claim, auto_unbox = TRUE, null = "null"), simplifyVector = FALSE))), 
                  "OWN_ATTEMPT_CONFIG_CODE_PATH_OR_SEED_CHANGED")
            } else gs100_new_json(claim, claim_path)
            existing <- list.files(file.path(attempt, "invocations"), pattern = "^launch_[0-9]{4}$")
            launch_number <- if (length(existing)) max(as.integer(sub("launch_", "", existing, fixed = TRUE))) + 
                1L else 1L
            gs100_assert(launch_number <= 9999L && (resume || launch_number == 1L), "OWN_INVOCATION_NUMBER_INVALID")
            launch_id <- sprintf("launch_%04d", launch_number)
            invocation <- gs100_worker_mkdir(file.path(attempt, "invocations", launch_id), attempt, exclusive = TRUE)
            private <- gs100_worker_mkdir(file.path(attempt, "private", launch_id), attempt, exclusive = TRUE)
            artifacts <- gs100_worker_mkdir(file.path(attempt, "artifacts", launch_id), attempt, exclusive = TRUE)
            receipts <- file.path(attempt, "receipts")
            stage <- "CLAIMED"
            event <- 0L
            last_integrating <- -Inf
            warnings <- list()
            progress <- function(new_stage, details = list()) {
                gs100_assert(is.character(new_stage) && length(new_stage) == 1L && grepl("^[A-Za-z][A-Za-z0-9_]*$", 
                  new_stage), "SAFE_AGGREGATE_STAGE_NAME_REQUIRED")
                stage <<- new_stage
                if (identical(new_stage, "integrating") && proc.time()[["elapsed"]] - last_integrating < 20) return(invisible(NULL))
                if (identical(new_stage, "integrating")) last_integrating <<- proc.time()[["elapsed"]]
                allow <- c("N_target", "N_fit", "conditional_index", "total", "mc", "n_mc", "rows_complete", "n_target", 
                  "endpoint", "cause", "status", "b", "fresh_A", "fresh_Y", "passed", "reverse_failure_recorded")
                details <- details[intersect(names(details), allow)]
                event <<- event + 1L
                record <- list(schema = "GLYCAEMIA_BASE11_B100_PROGRESS_V1", utc = gs100_utc(), pid = Sys.getpid(), 
                  b = b, launch = launch_id, stage = new_stage, details = details)
                gs100_new_json(record, file.path(attempt, "progress", sprintf("%s_%05d.json", launch_id, event)))
                cat(jsonlite::toJSON(record, auto_unbox = TRUE), "\n")
                flush.console()
                invisible(NULL)
            }
            started <- proc.time()[["elapsed"]]
            old_options <- options(body_tqe.upstream_A = NULL, body_tqe.prefix_A = NULL, body_pmm.old_A11 = NULL, 
                body_pmm.prefix_A = NULL)
            on.exit(options(old_options), add = TRUE)
            result <- tryCatch(withCallingHandlers({
                gs100_new_json(list(utc = gs100_utc(), pid = Sys.getpid(), b = b, seed = seed, resume = resume, 
                  claim = claim, roots = cfg[c("input_path", "adapter_root", "source_root", "y_root", "state_root", 
                    "output_root")], old_fitted_A_Y_or_prediction_cache_reads = 0L, automatic_retry = FALSE, participant_rows_in_reports = 0L), 
                  file.path(invocation, "INVOCATION.json"))
                progress("LOADING_OWN_FRESH_DRAW_AND_WEIGHTS")
                prepared <- gs100_load_draw(cfg, b)
                z <- prepared$z
                d <- prepared$d
                w <- prepared$w
                risks <- prepared$risks
                gs100_assert(prepared$manifest$b == b && prepared$manifest$seed == seed, "EXACT_CURRENT_BOOTSTRAP_NUMBER_AND_ORIGINAL_SEED_REQUIRED")
                rt <- new.env(parent = environment())
                sys.source(runtime_path, envir = rt)
                rt$gs_rt_task_name <- gs100_worker_task_name
                progress("LOADING_PINNED_SCIENTIFIC_RUNTIME")
                runtime <- rt$gs_load_runtime(cfg$source_root, cfg$y_root, state, attempt)
                gs100_new_json(runtime$source_pin_evidence, file.path(invocation, "SCIENTIFIC_SOURCE_PINS.json"))
                engine <- runtime$engine
                methods <- runtime$methods
                members <- z$members
                base <- gs100_vec(z$base)
                model_spec <- z$model_specifications$S44_BASE11
                allocation_order <- gs100_vec(model_spec$allocation_order)
                domains <- lapply(z$domains[allocation_order], gs100_vec)
                model_cols <- vapply(members, `[[`, "", "column")
                gs100_assert(length(members) == 11L && length(domains) == 4L && identical(names(members), gs100_vec(model_spec$members)) && 
                  identical(allocation_order, c("SMOKING", "BODY_COMPOSITION", "CLINICAL_LIPIDS", "GLYCAEMIA")) && 
                  setequal(unlist(domains, use.names = FALSE), model_cols) && !anyDuplicated(unlist(domains, use.names = FALSE)) && 
                  !isTRUE(model_spec$total_only), "EXACT_BASE11_FOUR_DOMAIN_SCIENTIFIC_CONTRACT_CHANGED")
                rt$gs_rt_validate_inputs(runtime, d, members, w)
                rng_spec <- list(kind = "L'Ecuyer-CMRG", normal = "Inversion", sample = "Rejection", mapping = "ROW_MAJOR_WIDTH11_COMMON_JOINT", 
                  independent_domain_streams = TRUE, seed = seed)
                binding <- list(panel_id = "S44_BASE11_HBA1C_LT48_SENSITIVITY_B100", bootstrap_id = b, sample_hash = prepared$manifest$fit_sample_hash, 
                  population_hash = z$target_hash, model_spec_hash = gs100_hash(model_spec), input_hash = prepared$files$input$sha256, 
                  rng_spec_hash = gs100_hash(rng_spec))
                identity <- list(schema = "GLYCAEMIA_BASE11_B100_OWN_ATTEMPT_IDENTITY_V1", claim = claim, binding = binding, 
                  rng_spec = rng_spec, draw_files = prepared$files, draw_manifest_sha256 = gs100_sha(file.path(gs100_draw_folder(cfg, 
                    b), "DRAW_MANIFEST.json")), sampling_hash = gs100_hash(as.integer(prepared$index$indices)), 
                  sample_hash = gs100_hash(as.character(d$eid)), data_hash = gs100_hash(d), risk_hash = gs100_hash(risks), 
                  weight_hash = gs100_hash(as.numeric(w)), pair_binding_hash = gs100_hash(prepared$selection$pair_binding), 
                  cluster_hash = gs100_hash(d$cluster_id), selection_source_hash = z$selection_source_hash, members_hash = gs100_hash(members), 
                  base = base, domains = domains, source_closure_hash = runtime$source_pin_evidence$closure_hash, 
                  source_pins_hash = gs100_hash(rt$gs_rt_pins()), R = R.version.string, platform = R.version$platform, 
                  packages = setNames(as.list(vapply(c("jsonlite", "digest", "survival", "splines", "nnet"), function(p) as.character(utils::packageVersion(p)), 
                    "")), c("jsonlite", "digest", "survival", "splines", "nnet")))
                identity_hash <- gs100_hash(identity)
                identity_path <- file.path(receipts, "MODEL_BINDING.private.rds")
                if (resume && file.exists(identity_path)) {
                  previous <- readRDS(identity_path)
                  gs100_assert(identical(previous$identity_hash, identity_hash) && identical(previous$identity, 
                    identity), "OWN_RESUME_SOURCE_CONFIG_SAMPLE_RISK_WEIGHT_BINDING_CHANGED")
                } else {
                  gs100_assert(!file.exists(identity_path), "PRESERVE_EXISTING_MODEL_BINDING")
                  gs100_new_rds(list(identity = identity, identity_hash = identity_hash), identity_path)
                }
                gs100_new_json(list(binding = binding, identity_hash = identity_hash, rng_spec = rng_spec, model_spec = model_spec, 
                  N_target = prepared$manifest$N_target, N_fit = nrow(d), new_draw = prepared$files, no_old_A_Y_or_prediction_reuse = TRUE), 
                  file.path(invocation, "MODEL_BINDING.json"))
                RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
                set.seed(seed)
                A_record <- if (resume) gs100_worker_completed_phase(receipts, "A", identity_hash, file.path(attempt, 
                  "private")) else NULL
                reused_A <- !is.null(A_record)
                if (is.null(A_record)) {
                  gs100_assert(is.function(engine$body_pmm_fit_a) && identical(engine$c39_fit_a, engine$body_pmm_fit_a) && 
                    all(c("prefix_fit", "old_fit") %in% names(formals(engine$body_pmm_fit_a))), "EXACT_PMM_FRESH_ALL_A_API_REQUIRED")
                  cb_options <- options(c39.a_checkpoint_callback = function(id, fit, index) progress("FRESH_A_CONDITIONAL", 
                    list(conditional_index = index, total = 11L)))
                  progress("FITTING_ALL_ELEVEN_FRESH_A")
                  A <- tryCatch(engine$body_pmm_fit_a(d, members, base, w, prefix_fit = NULL, old_fit = NULL), 
                    finally = options(cb_options))
                  gs100_worker_a_check(runtime, A, d, members, base, w)
                  a_file <- gs100_new_rds(A, file.path(private, "FRESH_A.private.rds"))
                  rng_file <- gs100_new_rds(gs100_worker_rng(), file.path(private, "FRESH_A_POSTFIT_RNG.private.rds"))
                  A_record <- list(phase = "A", complete = TRUE, identity_hash = identity_hash, fresh_all_eleven = TRUE, 
                    old_or_prefix_reuse = FALSE, files = list(A = a_file, rng = rng_file))
                  gs100_new_json(A_record, file.path(receipts, sprintf("FRESH_A_%04d.json", launch_number)))
                } else {
                  progress("RESTORING_OWN_EXACT_SAME_ATTEMPT_A")
                  A <- readRDS(gs100_worker_file(A_record$files$A, file.path(attempt, "private")))
                  gs100_worker_restore_rng(readRDS(gs100_worker_file(A_record$files$rng, file.path(attempt, "private"))))
                }
                if (reused_A) gs100_worker_a_check(runtime, A, d, members, base, w)
                progress("BUILDING_ORIGINAL_BASE11_Y_DESIGN")
                design_fn <- methods$new121_y_design_factory(d, members, base, domains, engine, cache_entries = 2L)
                fullX_hash <- gs100_hash(attr(design_fn, "fullX"))
                Y_record <- if (resume) gs100_worker_completed_phase(receipts, "Y", identity_hash, file.path(attempt, 
                  "private")) else NULL
                reused_Y <- !is.null(Y_record)
                if (is.null(Y_record)) {
                  reference <- gs100_new_rds(attr(design_fn, "reference"), file.path(private, "Y_DESIGN_REFERENCE.private.rds"))
                  owner <- attr(design_fn, "column_owner")
                  mediator_cols <- names(owner)[owner != "BASE"]
                  audits <- list()
                  progress("FITTING_FOUR_FRESH_FAMILY_ROBUST_CAUSE_SPECIFIC_Y")
                  y <- methods$new121_fit_y(engine, attr(design_fn, "fullX"), risks, w, d$cluster_id, mediator_cols, 
                    centre_actual = d$center, diagnostics = FALSE, audit_callback = function(record) {
                      audits[[length(audits) + 1L]] <<- gs100_new_rds(record, file.path(private, sprintf("Y_NATIVE_AUDIT_%04d.private.rds", 
                        length(audits) + 1L)))
                      progress("FRESH_Y_NATIVE_AUDIT", record[intersect(c("endpoint", "cause", "status"), names(record))])
                    })
                  gs100_worker_y_check(runtime, y, d, risks)
                  y_file <- gs100_new_rds(y, file.path(private, "FRESH_Y.private.rds"))
                  Y_record <- list(phase = "Y", complete = TRUE, identity_hash = identity_hash, A_sha256 = A_record$files$A$sha256, 
                    fullX_hash = fullX_hash, fresh_four_branches = TRUE, full_followup = TRUE, risk_horizon = 11L, 
                    files = c(list(Y = y_file, reference = reference), setNames(audits, sprintf("audit_%04d", seq_along(audits)))))
                  gs100_new_json(Y_record, file.path(receipts, sprintf("FRESH_Y_%04d.json", launch_number)))
                } else {
                  gs100_assert(identical(Y_record$A_sha256, A_record$files$A$sha256) && identical(Y_record$fullX_hash, 
                    fullX_hash), "OWN_Y_A_OR_ACTUAL_DESIGN_BINDING_CHANGED")
                  progress("RESTORING_OWN_EXACT_SAME_ATTEMPT_Y")
                  y <- readRDS(gs100_worker_file(Y_record$files$Y, file.path(attempt, "private")))
                  gs100_worker_y_check(runtime, y, d, risks)
                }
                predictor <- methods$new121_y_prediction(y, design_fn, engine, domains, cache_entries = 2L)
                validation <- NULL
                predict_risk <- function(b, M, C) {
                  if (is.null(validation)) {
                    take <- seq_len(min(16L, nrow(C)))
                    check <- predictor$validate(b, M[take, , drop = FALSE], C[take, , drop = FALSE])
                    mass <- vapply(check$all_states, function(p) max(abs(rowSums(p) - 1)), 0)
                    gs100_assert(all(is.finite(mass) & mass <= 1e-06), "REAL_FIRST_BLOCK_PROBABILITY_MASS_FAILURE")
                    validation <<- c(check$validation, list(max_mass_residual = max(mass), genuine_probability_check = TRUE, 
                      actual_rows = length(take), endpoints = c("IHD", "T2D"), passed = TRUE, source = "ORIGINAL_B_GT0_FIRST_REAL_16_ROW_GUARD", 
                      no_standalone_numeric_probe = TRUE))
                    gs100_new_json(validation, file.path(artifacts, "FIRST_REAL_RISK_BLOCK.json"))
                  }
                  predictor$predict(b, M, C)
                }
                progress("ORIGINAL_FIRST_REAL_RISK_BLOCK_GUARD")
                take <- seq_len(min(16L, nrow(d)))
                invisible(predict_risk(0, d[take, model_cols, drop = FALSE], d[take, , drop = FALSE]))
                gs100_new_json(list(passed = TRUE, b = b, identity_hash = identity_hash, fresh_all_A = TRUE, fresh_four_Y = TRUE, 
                  validation = validation, K50_started = FALSE, no_B0_K3_distribution_diagnostic = TRUE, no_slow_quadrature_probe = TRUE), 
                  file.path(invocation, "PREINTEGRATION_GATE_PASSED.json"))
                draw_joint <- function(b, C, uniforms, requested_domain = NULL) engine$c39_draw_a(A, C, b, uniforms[, 
                  seq_along(members), drop = FALSE])[, model_cols, drop = FALSE]
                support <- function(M, C, context) {
                  requested <- Filter(function(m) m$column %in% names(M), members)
                  gs100_assert(length(requested) > 0L, "EMPTY_REQUESTED_MEDIATOR_SUPPORT")
                  engine$c39_support(M, requested)
                  TRUE
                }
                checkpoint <- NULL
                cp_paths <- list.files(receipts, pattern = "^checkpoint_MC[0-9]{2}_launch_[0-9]{4}[.]json$", full.names = TRUE)
                if (resume && length(cp_paths)) {
                  cp_records <- lapply(cp_paths, function(path) {
                    gs100_assert(gs100_inside(gs100_norm(path), gs100_norm(receipts)), "OWN_MC_MANIFEST_ESCAPE_REJECTED")
                    rec <- jsonlite::fromJSON(path, simplifyVector = FALSE)
                    gs100_assert(isTRUE(rec$complete) && identical(rec$identity_hash, identity_hash) && identical(rec$A_sha256, 
                      A_record$files$A$sha256) && identical(rec$Y_sha256, Y_record$files$Y$sha256) && is.numeric(rec$completed_mc) && 
                      length(rec$completed_mc) == 1L && rec$completed_mc >= 1L && rec$completed_mc <= 50L && rec$completed_mc == 
                      floor(rec$completed_mc), "OWN_MC_CONFIG_A_Y_SAMPLE_WEIGHT_SOURCE_BINDING_CHANGED")
                    gs100_worker_file(rec$file, file.path(attempt, "private"))
                    rec
                  })
                  rec <- cp_records[[which.max(vapply(cp_records, function(x) as.integer(x$completed_mc), 1L))]]
                  checkpoint <- readRDS(gs100_worker_file(rec$file, file.path(attempt, "private")))
                  gs100_assert(identical(as.integer(checkpoint$completed_mc), as.integer(rec$completed_mc)), "OWN_MC_MANIFEST_PAYLOAD_PASS_ACCOUNTING_CHANGED")
                  progress("RESUMING_OWN_COMMITTED_MC_PASS", list(mc = rec$completed_mc, n_mc = 50L))
                }
                save_checkpoint <- function(state) {
                  k <- as.integer(state$completed_mc)
                  gs100_assert(length(k) == 1L && !is.na(k) && k >= 1L && k <= 50L, "MC_CHECKPOINT_NUMBER_INVALID")
                  file <- gs100_new_rds(state, file.path(private, sprintf("MC_%02d.private.rds", k)))
                  gs100_new_json(list(complete = TRUE, b = b, completed_mc = k, identity_hash = identity_hash, 
                    A_sha256 = A_record$files$A$sha256, Y_sha256 = Y_record$files$Y$sha256, file = file), file.path(receipts, 
                    sprintf("checkpoint_MC%02d_%s.json", k, launch_id)))
                }
                progress("INTEGRATING_K50_FORWARD_REVERSE_DEP")
                integral <- methods$nmr_product_stream(d, w, domains, draw_joint, predict_risk, support, uniform_dim = 11L, 
                  n_mc = 50L, block_size = 512L, seed = seed, allocation_order = allocation_order, include_reverse = TRUE, 
                  endpoints = c("IHD", "T2D"), binding = binding, checkpoint = checkpoint, save_checkpoint = save_checkpoint, 
                  progress = function(v) progress(v$stage, v[setdiff(names(v), "stage")]))
                order_success <- gs100_worker_validate_integral(integral, allocation_order, b, seed)
                for (name in c("components", "paired", "scenarios")) {
                  path <- file.path(artifacts, paste0(toupper(name), ".csv"))
                  gs100_assert(!file.exists(path), "PRESERVE_EXISTING_AGGREGATE_RESULT_TABLE")
                  utils::write.csv(integral[[name]], path, row.names = FALSE)
                }
                integral_file <- gs100_new_rds(integral, file.path(private, "INTEGRATION.private.rds"))
                progress("REPLICATE_COMPLETE", list(fresh_A = TRUE, fresh_Y = TRUE, passed = TRUE))
                list(schema = "GLYCAEMIA_SENSITIVITY_BASE11_B100_REPLICATE_V1", status = "SUCCESS", utc = gs100_utc(), 
                  b = b, model_id = "S44_BASE11", seed = seed, binding = binding, identity_hash = identity_hash, 
                  N = prepared$manifest[c("N_target", "N_frozen_target", "N_fit", "N_standardization", "endpoints", 
                    "IPW")], components = integral$components, paired = integral$paired, scenarios = integral$scenarios, 
                  stage = "REPLICATE_COMPLETE", reason = "COMPLETE_FORWARD_K50", safe_code = "SUCCESS", order_success = order_success, 
                  audit = list(kernel = gs100_worker_kernel_audit(integral$audit), fresh_all_A = TRUE, fresh_four_Y = TRUE, 
                    original_seed = seed, rng_spec = rng_spec, first_real_risk_block = validation, no_B0_K3_distribution_diagnostic = TRUE, 
                    no_slow_quadrature_probe = TRUE, old_A_Y_or_prediction_cache_reuse = FALSE, own_same_attempt_A_resume = reused_A, 
                    own_same_attempt_Y_resume = reused_Y, native_prediction = predictor$stats(), selection_refits = 1L, 
                    bootstrap_sampling_performed = TRUE, full_followup_preserved = TRUE, report_horizon_years = 11L), 
                  artifacts = list(A = A_record$files$A, Y = Y_record$files$Y, integration = integral_file), participant_rows_in_reports = 0L)
            }, warning = function(warning) {
                warnings[[length(warnings) + 1L]] <<- list(stage = stage, message = conditionMessage(warning))
                invokeRestart("muffleWarning")
            }), error = function(error) {
                gs100_new_rds(list(stage = stage, class = class(error), message = conditionMessage(error)), file.path(private, 
                  "FAILURE_DETAIL.private.rds"))
                status <- gs100_worker_failure_class(error)
                list(schema = "GLYCAEMIA_SENSITIVITY_BASE11_B100_REPLICATE_V1", status = status, utc = gs100_utc(), 
                  b = b, model_id = "S44_BASE11", seed = seed, stage = stage, reason = if (identical(status, "STATISTICAL_FAILURE")) "BOOTSTRAP_STATISTICAL_STAGE_FAILURE" else "WORKER_TECHNICAL_STAGE_FAILURE", 
                  safe_code = paste0(status, "_REVIEW_PRIVATE_DETAIL"), automatic_retry = FALSE, resume_possible = identical(status, 
                    "TECHNICAL_FAILURE"), participant_rows_in_reports = 0L)
            })
            gs100_new_rds(warnings, file.path(private, "WARNINGS.private.rds"))
            result$warning_count <- length(warnings)
            result$elapsed_seconds <- proc.time()[["elapsed"]] - started
            result$launch <- launch_id
            result$resume_requested <- resume
            gs100_new_json(result, file.path(invocation, "RESULT.json"))
            if (result$status %in% c("SUCCESS", "STATISTICAL_FAILURE")) gs100_new_json(result, file.path(attempt, 
                "RESULT.json"))
            if (identical(result$status, "SUCCESS")) 0L else if (identical(result$status, "STATISTICAL_FAILURE")) 2L else 1L
        }
        if (sys.nframe() == 0L) {
            status <- tryCatch(gs100_worker_main(), error = function(error) {
                cat("{\"status\":\"TECHNICAL_FAILURE\",\"stage\":\"SETUP\",\"reason\":\"WORKER_SETUP_FAILURE\",\"safe_code\":\"WORKER_SETUP_FAILURE_REVIEW_SCOPED_PRIVATE_STATE\"}\n")
                1L
            })
            quit(status = status)
        }
    }), envir = context)
}

