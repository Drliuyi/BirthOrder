# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

runtime_files <- c("phase2_worker.R", "phase2_policy.R", "phase1_helpers.R")

first_failure_models_saved <- FALSE

emit <- function(stage, value = list()) {
    cat("PHASE2_JSON:", jsonlite::toJSON(c(list(stage = stage, worker_id = worker_id, host = host), value), auto_unbox = TRUE, 
        digits = 17, na = "null", null = "null"), "\n", sep = "")
    flush.console()
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        args <- commandArgs(TRUE)
        if (length(args) != 4L) stop("ROOT NEW_G_REPAIR_RUN GROUP_MANIFEST.json worker_id required", call. = FALSE)
        script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        if (length(script_arg) != 1L) stop("Rscript --file execution required", call. = FALSE)
        script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
        HERE <- dirname(script_path)
        source(file.path(HERE, "phase1_helpers.R"))
        source(file.path(HERE, "phase2_policy.R"))
        for (key in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")) do.call(Sys.setenv, 
            setNames(list("1"), key))
        for (package in c("jsonlite", "digest", "survival", "nnet")) {
            rw_assert(requireNamespace(package, quietly = TRUE), paste("DEPENDENCY_MISSING", package))
            suppressPackageStartupMessages(library(package, character.only = TRUE))
        }
        ROOT <- normalizePath(args[1L], mustWork = TRUE)
        manifest_path <- normalizePath(p2_rebase(args[3L]), mustWork = TRUE)
        manifest_sha <- rw_sha(manifest_path)
        manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
        worker_id <- as.integer(args[4L])
        host <- Sys.getenv("FF_HOST")
        rw_assert(identical(manifest$schema, "FOUR_PROTEIN_WARNING_PHASE2_V1") && as.integer(manifest$worker_count) == 
            20L && worker_id %in% 1:20, "P2_FIXED20_GROUP_MANIFEST_REQUIRED")
        rw_assert(identical(host, if (worker_id <= 16L) "winpc" else "macmini"), "P2_STATIC_HOST_WORKER_PARTITION_CHANGED")
        rw_assert(identical(Sys.getenv("FF_RECOVERY_EXECUTION_AUTHORIZED"), "YES_PHASE2_SELECTIVE_WARNING_RECOVERY"), 
            "P2_EXPLICIT_SELECTIVE_RECOVERY_AUTHORIZATION_REQUIRED")
        policy_sha <- rw_sha(file.path(HERE, "phase2_policy.R"))
        rw_assert(identical(policy_sha, manifest$recovery_policy_sha256), "P2_REVIEWED_RECOVERY_POLICY_CHANGED")
        policy_spec_sha <- digest::digest(p2_policy_spec(), algo = "sha256", serialize = FALSE)
        if (!is.null(manifest$method_review_sha256)) rw_assert(identical(manifest$method_review_sha256, p2_method_review_sha256()), 
            "P2_METHOD_REVIEW_DECISION_CHANGED")
        runtime_files <- c("phase2_worker.R", "phase2_policy.R", "phase1_helpers.R")
        runtime_pins <- stats::setNames(vapply(file.path(HERE, runtime_files), rw_sha, ""), runtime_files)
        canonical_out <- args[2L]
        rw_assert(grepl("^/mnt/g/UKB_results/birth/border_fourprotein_warning_recovery_20261031_v1/run_[A-Za-z0-9_.-]+$", 
            canonical_out), "P2_CANONICAL_NEW_G_REPAIR_RUN_REQUIRED")
        OUT <- p2_rebase(canonical_out)
        if (!dir.exists(OUT)) rw_assert(dir.create(OUT, recursive = TRUE, showWarnings = FALSE), "P2_REPAIR_RUN_CREATE_FAILED")
        OUT <- normalizePath(OUT, mustWork = TRUE)
        allowed_parent <- normalizePath(p2_rebase(birth_private_setting("PRIVATE_BORDER_FOURPROTEIN_WARNING_RECOVERY_20261031_V1")), 
            mustWork = TRUE)
        rw_assert(rw_inside(OUT, allowed_parent), "P2_G_REPAIR_SYMLINK_BOUNDARY_ESCAPE")
        groups <- manifest$groups
        group_ids <- vapply(groups, function(g) g$group_id, "")
        rw_assert(length(groups) > 0L && !anyDuplicated(group_ids) && all(grepl("^[A-Za-z0-9_.-]+$", group_ids)), 
            "P2_UNIQUE_GROUP_IDS_REQUIRED")
        all_workers <- as.integer(vapply(groups, function(g) g$assigned_worker, 0))
        rw_assert(all(all_workers %in% 1:20) && all(vapply(groups, function(g) identical(g$assigned_host, if (as.integer(g$assigned_worker) <= 
            16L) "winpc" else "macmini"), logical(1L))), "P2_GROUP_STATIC_ASSIGNMENT_INVALID")
        mine <- Filter(function(g) as.integer(g$assigned_worker) == worker_id, groups)
        groups_total_assigned <- length(mine)
        canary_value <- Sys.getenv("FF_PHASE2_MAX_GROUPS")
        rw_assert(!nzchar(canary_value) || identical(canary_value, "1"), "P2_ONLY_FIXED_ONE_GROUP_CANARY_LIMIT_ALLOWED")
        canary_mode <- identical(canary_value, "1")
        if (canary_mode) mine <- head(mine, 1L)
        singular <- manifest$known_singular_cases
        rw_assert(length(singular) == 125L, "P2_EXACT_KNOWN125_SINGULAR_LEDGER_REQUIRED")
        singular_keys <- vapply(singular, function(z) p2_key(z$indicator_id, z$bootstrap_id, z$endpoint_id), "")
        rw_assert(!anyDuplicated(singular_keys), "P2_KNOWN_SINGULAR_LEDGER_DUPLICATED")
        source_names <- c("10.prepare.R", "20.a_reuse.R", "30.outcome.R", "31.aj_backend.R", "61.center_strata.R", 
            "40.worker.R", "50.panel.R", "52.a_task_cache.R", "53.a_derived_reader.R", "54.incremental_a.R")
        source_pins <- stats::setNames(vapply(file.path(ROOT, "f", source_names), rw_sha, ""), source_names)
        if (!is.null(manifest$source_files_sha256)) for (k in names(manifest$source_files_sha256)) rw_assert(k %in% 
            names(source_pins) && identical(source_pins[[k]], manifest$source_files_sha256[[k]]), paste("P2_FROZEN_SOURCE_CHANGED", 
            k))
        for (path in file.path(ROOT, "f", source_names)) source(path)
        b16_a_derived_install(globalenv(), file.path(ROOT, "f/53.a_derived_reader.R"))
        backend_work <- Sys.getenv("FF_WORK")
        if (!is.null(manifest$backend_work_roots[[host]])) backend_work <- manifest$backend_work_roots[[host]]
        rw_assert(nzchar(backend_work) && dir.exists(backend_work), "P2_EXISTING_AJ_BACKEND_WORK_REQUIRED")
        backend_source_sha <- rw_sha(file.path(ROOT, "src/aj.c"))
        backend_runtime <- list(R = R.version$version.string, platform = R.version$platform, arch = R.version$arch, 
            pointer_size = .Machine$sizeof.pointer, extension = .Platform$dynlib.ext, isolated_build_contract = "R_CMD_SHLIB_V2_EXPLICIT_RUNTIME_COMPILER_PATH", 
            compiler_paths = Sys.which(c("cc", "clang", "gcc", "arm64-apple-darwin20.0.0-clang")))
        build <- file.path(backend_work, "software", paste0("aj_", substr(backend_source_sha, 1L, 16L), "_", substr(digest::digest(backend_runtime, 
            algo = "sha256"), 1L, 12L)))
        rw_assert(all(file.exists(file.path(build, c("BUILD_MANIFEST.rds", "aj.c", paste0("ff_aj", .Platform$dynlib.ext))))), 
            "P2_PREBUILT_AJ_REQUIRED_NO_COMPILATION")
        backend <- ff_aj_initialize(ROOT, backend_work, require_compiled = TRUE)
        namespace_rng_before <- b16_a_cache_rng()
        namespaces_before <- loadedNamespaces()
        runtime <- ff_y_runtime(require = TRUE)
        namespace_rng_after <- b16_a_cache_rng()
        namespace_audit <- list(new_namespaces = as.list(setdiff(loadedNamespaces(), namespaces_before)), RNG_existed_before = !is.null(namespace_rng_before), 
            RNG_exists_after = !is.null(namespace_rng_after), RNG_changed_on_dependency_initialization = !identical(namespace_rng_before, 
                namespace_rng_after), RNG_was_not_reset = TRUE)
        ctx <- ff_load(ROOT)
        ctx$allow_missing_a_recompute <- FALSE
        rw_assert(isTRUE(ctx$v4_decision$source_frozen_for_production) && !rw_inside(OUT, ctx$run), "P2_FROZEN_SOURCE_AND_DISJOINT_REPAIR_REQUIRED")
        counts <- list(groups_total_assigned = groups_total_assigned, groups_selected = length(mine), canary_mode = canary_mode, 
            groups_complete = 0L, cases_complete = 0L, cases_skipped_verified = 0L, old_VALID_skipped = 0L, known_singular_skipped = 0L, 
            conditional_estimates = 0L, clean_VALID_estimates = 0L, failures = 0L, Y_cause_fits = 0L, A_fits = 0L, 
            selection_fits = 0L, draws_generated = 0L)
        first_failure_models_saved <- FALSE
        emit <- function(stage, value = list()) {
            cat("PHASE2_JSON:", jsonlite::toJSON(c(list(stage = stage, worker_id = worker_id, host = host), value), 
                auto_unbox = TRUE, digits = 17, na = "null", null = "null"), "\n", sep = "")
            flush.console()
        }
        base_identity <- list(manifest_sha256 = manifest_sha, recovery_policy_sha256 = policy_sha, worker_id = worker_id, 
            recovery_runtime_sha256 = digest::digest(runtime_pins, algo = "sha256"), method_hash = ctx$ff_method_hash, 
            source_hash = ctx$source_hash, input_hash = ctx$input_lock_sha256)
        emit("READY", list(groups = length(mine), namespace_initialization = namespace_audit))
        for (group in mine) {
            id <- as.character(group$indicator_id)
            bid <- as.integer(group$bootstrap_id)
            rw_assert(id %in% c("olink__fcrl6", "olink__gzma", "olink__gzmh", "olink__klrd1") && bid %in% 0:1000, 
                "P2_EXACT_PROTEIN_DRAW_SCOPE_REQUIRED")
            group_out <- file.path(OUT, "groups", group$group_id)
            group_identity <- c(base_identity, list(group_id = group$group_id, indicator_id = id, bootstrap_id = bid))
            if (file.exists(file.path(group_out, "COMPLETE.json"))) {
                prior <- p2_verified(group_out, group_identity)
                for (entry in prior$cases) {
                  verified <- p2_verified(file.path(group_out, entry$relative_directory), c(group_identity, list(endpoint_id = entry$endpoint_id)))
                  if (!is.null(verified$representative_failure_models)) first_failure_models_saved <- TRUE
                }
                counts$groups_complete <- counts$groups_complete + 1L
                counts$cases_skipped_verified <- counts$cases_skipped_verified + length(prior$cases)
                emit("GROUP_SEALED_HASH_VERIFIED_SKIPPED", list(group_id = group$group_id))
                next
            }
            if (dir.exists(group_out)) rw_assert(identical(Sys.getenv("FF_PHASE2_RESUME_AUTHORIZED"), "YES_VERIFY_SEALED_ONLY"), 
                "P2_PARTIAL_GROUP_PRESERVED_EXPLICIT_RESUME_REQUIRED") else group_out <- rw_newdir(group_out)
            emit("GROUP_INPUT_AUTHENTICATION", list(group_id = group$group_id))
            task <- normalizePath(p2_rebase(group$task_path), mustWork = TRUE)
            rw_assert(rw_inside(task, ctx$run), "P2_TASK_OUTSIDE_CURRENT_FROZEN_PRODUCTION")
            task_seal_path <- file.path(task, "FF_COMPLETE.json")
            rw_check(task_seal_path, group$task_seal_sha256)
            seal <- jsonlite::fromJSON(task_seal_path, simplifyVector = FALSE)
            rw_assert(isTRUE(seal$complete) && identical(seal$indicator_id, id) && as.integer(seal$bootstrap_id) == 
                bid && identical(seal$source_hash, ctx$source_hash) && identical(seal$method_hash, ctx$ff_method_hash) && 
                identical(seal$input_hash, ctx$input_lock_sha256) && identical(seal$method_spec_hash, ctx$v4_decision$method_spec_sha256), 
                "P2_CURRENT_V4_TASK_BINDING_CHANGED")
            for (name in c("FF_RESULTS.csv", "DRAW.json")) rw_check(file.path(task, name), seal$artifacts_sha256[[name]])
            rows <- read.csv(file.path(task, "FF_RESULTS.csv"), colClasses = "character", check.names = FALSE)
            rw_assert(nrow(rows) == 5L && setequal(rows$endpoint_id, c("250", "250.2", "411", "411.3", "411.4")) && 
                all(rows$indicator_id == id) & all(as.integer(rows$bootstrap_id) == bid), "P2_FIVE_CURRENT_TASK_ROWS_REQUIRED")
            for (i in seq_len(nrow(rows))) ff_w_v4_row(ctx, as.list(rows[i, , drop = FALSE]))
            reference_endpoint <- group$reference_a_endpoint
            if (is.null(reference_endpoint) && is.list(group$reference_a)) reference_endpoint <- group$reference_a$endpoint_id
            rw_assert(reference_endpoint %in% c("411", "411.3", "411.4"), "P2_CORONARY_REFERENCE_A_ENDPOINT_REQUIRED")
            reference <- as.list(rows[rows$endpoint_id == reference_endpoint, , drop = FALSE])
            for (k in c("sample_hash", "target_hash", "selection_weight_hash")) rw_assert(identical(reference[[k]], 
                group[[k]]), paste("P2_SEALED_REFERENCE_A_BINDING_CHANGED", k))
            draw_receipt <- jsonlite::fromJSON(file.path(task, "DRAW.json"), simplifyVector = FALSE)
            draw_path <- p2_rebase(group$draw_rds_path)
            rw_check(draw_path, group$draw_rds_sha256)
            rw_assert(identical(draw_receipt$actual_draw_private_sha256, group$draw_rds_sha256) && identical(draw_receipt$draw_hash, 
                group$draw_hash) && identical(draw_receipt$draw_hash, seal$draw_hash) && as.integer(draw_receipt$bootstrap_id) == 
                bid && identical(draw_receipt$parent_frame_hash, ctx$parent_frame_hash) && identical(draw_receipt$parent_cluster_membership_hash, 
                ctx$parent_cluster_membership_hash), "P2_DRAW_RECEIPT_CHANGED")
            saved <- readRDS(draw_path)
            ix <- saved$sampled_parent_rows
            rw_assert(is.numeric(ix) && length(ix) > 0L && all(ix %in% seq_len(nrow(ctx$d_full))) && identical(ff_a_hash(ix), 
                draw_receipt$draw_hash) && identical(ff_a_hash(saved$cluster_multiplicity), draw_receipt$multiplicity_hash) && 
                identical(saved$parent_frame_hash, ctx$parent_frame_hash) && identical(saved$cluster_order_hash, 
                ctx$cluster_order_hash), "P2_ACTUAL_RETAINED_DRAW_CONTENT_CHANGED")
            ctx$actual_draw_hash <- draw_receipt$draw_hash
            descriptor_path <- p2_rebase(group$a_descriptors_path)
            rw_check(descriptor_path, group$a_descriptors_sha256)
            rw_assert(rw_inside(descriptor_path, file.path(ctx$run, "private/a_incremental")), "P2_A_DESCRIPTOR_OUTSIDE_CURRENT_PRIVATE")
            descriptors <- jsonlite::fromJSON(descriptor_path, simplifyVector = FALSE)$rows
            rw_assert(length(descriptors) == 5L && setequal(vapply(descriptors, function(z) z$endpoint_id, ""), 
                rows$endpoint_id) && all(vapply(descriptors, function(z) identical(z$indicator_id, id) && as.integer(z$bootstrap_id) == 
                bid && identical(z$format, "derived_group_rds") && isTRUE(z$individual_weight_vector_retained), 
                logical(1L))), "P2_EXACT_FIVE_RETAINED_A_DESCRIPTORS_REQUIRED")
            ctx$a_source_for <- local({
                ds <- descriptors
                gid <- id
                gbid <- bid
                function(id, bid, endpoint) {
                  rw_assert(identical(id, gid) && as.integer(bid) == gbid, "P2_A_GROUP_CHANGED")
                  out <- Filter(function(z) identical(z$endpoint_id, endpoint), ds)
                  rw_assert(length(out) == 1L, "P2_A_ENDPOINT_DESCRIPTOR_REQUIRED")
                  out[[1L]]
                }
            })
            masks <- lapply(rows$endpoint_id, function(e) ctx$eligibility[[id]] & ctx$endpoints[[e]]$eligibility)
            rw_assert(all(vapply(masks, function(m) identical(m, masks[[1L]]), logical(1L))), "P2_FIVE_TARGETS_DIFFER_NO_A_SHARING")
            pr <- ix[masks[[1L]][ix]]
            d <- ctx$d_full[pr, , drop = FALSE]
            ca <- ctx$features[[id]]
            column <- ca$column
            base <- unique(c(ctx$cfg$base_covariates, ca$technical_covariates))
            cache <- new.env(parent = emptyenv())
            cache$entries <- new.env(parent = emptyenv())
            cache$active <- TRUE
            cache$task <- c(id, as.character(bid))
            ctx$b16_a_task_cache <- cache
            group_cases <- list()
            endpoints <- group$endpoints
            requested <- vapply(endpoints, function(z) z$endpoint_id, "")
            rw_assert(length(requested) > 0L && !anyDuplicated(requested) && all(requested %in% rows$endpoint_id), 
                "P2_SELECTIVE_ENDPOINT_SCOPE_INVALID")
            for (item in endpoints) {
                epid <- as.character(item$endpoint_id)
                cell_dir <- file.path(group_out, paste0("endpoint_", gsub(".", "_", epid, fixed = TRUE)))
                cell_identity <- c(group_identity, list(endpoint_id = epid))
                if (dir.exists(cell_dir)) {
                  q <- p2_verified(cell_dir, cell_identity)
                  counts$cases_skipped_verified <- counts$cases_skipped_verified + 1L
                  if (!is.null(q$representative_failure_models)) first_failure_models_saved <- TRUE
                  group_cases[[epid]] <- list(endpoint_id = epid, relative_directory = basename(cell_dir), status = q$status, 
                    complete_sha256 = rw_sha(file.path(cell_dir, "COMPLETE.json")))
                  next
                }
                cell_dir <- rw_newdir(cell_dir)
                private <- rw_newdir(file.path(cell_dir, "private"))
                old <- as.list(rows[rows$endpoint_id == epid, , drop = FALSE])
                rw_assert(identical(old$status, item$original_state) && identical(old$reason, item$reason), "P2_SELECTED_ORIGINAL_STATE_OR_REASON_CHANGED")
                q <- c(ff_w_base_row(ctx, id, bid, epid), list(recovery_policy_sha256 = policy_sha, recovery_policy_spec_sha256 = policy_spec_sha, 
                  manifest_sha256 = manifest_sha, worker_id = worker_id, method_review_decision_sha256 = p2_method_review_sha256(), 
                  original_status = old$status, original_reason = old$reason, group_id = group$group_id, original_task_path = group$task_path, 
                  original_task_seal_sha256 = group$task_seal_sha256, finite_joint_MLE_certified = FALSE, conditional_recovery_not_native_warning_gate_pass = TRUE, 
                  original_VALID_outputs_modified = FALSE, production_queue_touched = FALSE, A_fits = 0L, selection_fits = 0L, 
                  draws_generated = 0L, Y_cause_fits_attempted = 0L))
                state <- new.env(parent = emptyenv())
                for (k in c("models", "warnings", "controls", "model_files", "native_errors", "warning_gates", 
                  "audit", "designs")) state[[k]] <- list()
                state$fits_attempted <- 0L
                cache_object <- NULL
                emit("CELL_START", list(group_id = group$group_id, endpoint_id = epid, counts = counts))
                tryCatch({
                  if (identical(old$status, "VALID")) {
                    q$status <- "SKIPPED_OLD_VALID_PRESERVED"
                    q$reason <- "EXISTING_CURRENT_V4_VALID_NEVER_REFIT"
                    counts$old_VALID_skipped <- counts$old_VALID_skipped + 1L
                  } else if (p2_key(id, bid, epid) %in% singular_keys) {
                    q$status <- "SKIPPED_KNOWN_SINGULAR_NO_REFIT"
                    q$reason <- old$reason
                    counts$known_singular_skipped <- counts$known_singular_skipped + 1L
                  } else {
                    rw_assert(identical(ctx$variables[[id]]$status, "AVAILABLE_FROZEN") && identical(ctx$endpoints[[epid]]$status, 
                      "AVAILABLE_FROZEN"), "P2_FROZEN_FEATURE_OR_ENDPOINT_UNAVAILABLE")
                    rng_before <- b16_a_cache_rng()
                    a <- ff_a_scenarios(ctx, id, bid, epid, pr, d, base)
                    rw_assert(identical(rng_before, b16_a_cache_rng()) && !grepl("RECOMPUTED", a$ledger$A_origin), 
                      "P2_A_RECOMPUTE_OR_RNG_CHANGE_FORBIDDEN")
                    rw_assert(!is.null(a$alpha), "P2_CONTINUOUS_PROTEIN_A_REQUIRED")
                    cc <- which(a$mi$complete)
                    rr <- ctx$endpoints[[epid]]$risk[pr, , drop = FALSE][cc, , drop = FALSE]
                    data <- a$data
                    w <- a$mi$w
                    q$target_hash <- ff_a_hash(pr)
                    q$sample_hash <- ff_a_hash(pr[cc])
                    q$selection_weight_hash <- ff_a_hash(format(w, digits = 17, scientific = TRUE))
                    for (k in c("target_hash", "sample_hash", "selection_weight_hash")) {
                      rw_assert(identical(q[[k]], group[[k]]), paste("P2_ACTUAL_A_VS_SEALED_GROUP_BINDING_CHANGED", 
                        k))
                      if (epid %in% c("411", "411.3", "411.4")) rw_assert(identical(q[[k]], old[[k]]), paste("P2_FAILED_CORONARY_BINDING_CHANGED", 
                        k))
                    }
                    q$actual_draw_hash <- draw_receipt$draw_hash
                    q$draw_private_sha256 <- group$draw_rds_sha256
                    q$A_reuse <- as.list(a$ledger[1L, , drop = FALSE])
                    q$A_descriptor_sha256 <- group$a_descriptors_sha256
                    q$actual_basis_reference_hash <- v5_hash(as.character(a$reference_basis))
                    q$N_target <- nrow(d)
                    q$N_complete <- nrow(data)
                    if (isTRUE(item$phase1_reuse)) {
                      p1 <- p2_import_phase1(item, q, group, policy_sha)
                      q <- p1$result
                      cache_object <- p1$cache
                      emit("CELL_REUSED_AUTHENTICATED_PHASE1_NO_Y_FIT", list(group_id = group$group_id, endpoint_id = epid))
                    } else {
                      X <- ctx$engine$b3_design(data, column, base, a$reference_basis)$X
                      design <- ff_y_matrix(X)
                      endpoint <- ff_y_risk(rr, ctx$tau_report)
                      centres <- ff_ce_centers(design$X, centre_actual = data$center)
                      changed <- ff_y_matrix(ff_ce_columns(design$X, centres$reference$columns))
                      state$designs <- list(disease = if (identical(epid, "411.3")) changed else design, death = changed)
                      audit <- function(rec) state$audit[[rec$cause]] <- rec[setdiff(names(rec), c("original_coefficients", 
                        "original_vcov"))]
                      fit_fn <- p2_clone_fit(state, column, private)
                      emit("CELL_NATIVE_Y_FIT", list(group_id = group$group_id, endpoint_id = epid))
                      warm <- fit_fn(X, rr, w, tau_report = ctx$tau_report, endpoint_id = epid, centre_actual = data$center, 
                        diagnostics = FALSE, audit_callback = audit)
                      originals <- stats::setNames(lapply(c("disease", "death"), function(cause) state$models[[paste0("original_", 
                        cause)]]), c("disease", "death"))
                      original <- rw_rebuild_family(warm, originals, X, rr, w, centres, ctx$tau_report)
                      diagnostics <- p2_diagnostics(state, column, endpoint, w)
                      rw_assert(all(vapply(diagnostics, function(z) isTRUE(z$core_stability_passed) && isTRUE(z$core_positive_information_both) && 
                        isTRUE(z$original$finite_coefficients_and_covariance) && isTRUE(z$warm$finite_coefficients_and_covariance) && 
                        isTRUE(z$original$warning_map$qualified) && isTRUE(z$warm$warning_map$qualified) && isTRUE(z$original$warning_support_qualified) && 
                        isTRUE(z$warm$warning_support_qualified), logical(1L))), "P2_CORE_INFORMATION_STABILITY_OR_RECOGNIZED_WARNING_SUPPORT_RULE_FAILED")
                      rw_assert(isTRUE(original$native_validation$passed) && isTRUE(warm$native_validation$passed), 
                        "P2_NATIVE_VALIDATION_FAILED")
                      compact_original <- p2_compact_family(original)
                      compact_warm <- p2_compact_family(warm)
                      compact_checks <- list()
                      probe <- unique(as.integer(round(seq(1, nrow(data), length.out = min(16L, nrow(data))))))
                      compare_compact <- function(SX, label) {
                        p0 <- ff_ce_predict_y_matrix(original, SX[probe, , drop = FALSE], ctx$tau_report)
                        c0 <- ff_ce_predict_y_matrix(compact_original, SX[probe, , drop = FALSE], ctx$tau_report)
                        p1 <- ff_ce_predict_y_matrix(warm, SX[probe, , drop = FALSE], ctx$tau_report)
                        c1 <- ff_ce_predict_y_matrix(compact_warm, SX[probe, , drop = FALSE], ctx$tau_report)
                        difference <- max(abs(p0 - c0), abs(p1 - c1))
                        rw_assert(identical(unname(p0), unname(c0)) && identical(unname(p1), unname(c1)), "P2_COMPACT_NATIVE_PREDICTION_NOT_BITWISE_IDENTICAL")
                        compact_checks[[label]] <<- list(rows = length(probe), max_abs_difference = difference, 
                          bitwise_identical = TRUE)
                      }
                      emit("CELL_NATIVE_CIF_PREDICTION", list(group_id = group$group_id, endpoint_id = epid))
                      predictions <- list(observed = rw_probability_pair(original, warm, X, w, ctx$tau_report))
                      compare_compact(X, "observed")
                      risks <- warm_risks <- c(R00 = 0, R10 = 0, R11 = 0)
                      for (scenario in names(risks)) {
                        nd <- data
                        nd[[column]] <- if (scenario == "R11") a$M1 else a$M0
                        nd$B <- if (scenario == "R00") 0 else 1
                        SX <- ctx$engine$b3_design(nd, column, base, a$reference_basis)$X
                        pred <- rw_probability_pair(original, warm, SX, w, ctx$tau_report)
                        predictions[[scenario]] <- pred
                        risks[[scenario]] <- pred$original_weighted_disease_risk
                        warm_risks[[scenario]] <- pred$warm_weighted_disease_risk
                        compare_compact(SX, scenario)
                      }
                      rw_assert(all(vapply(predictions, function(p) isTRUE(p$original$legal) && isTRUE(p$warm$legal), 
                        logical(1L))), "P2_OBSERVED_OR_SCENARIO_PROBABILITY_ILLEGAL")
                      rw_assert(isTRUE(ff_y_effect_closure(risks)$passed) && isTRUE(ff_y_effect_closure(warm_risks)$passed), 
                        "P2_EFFECT_CLOSURE_FAILED")
                      q$original_effects <- rw_effects(risks)
                      q$warm_effects <- rw_effects(warm_risks)
                      for (k in c("R00", "R10", "R11", "TE", "DE", "IE", "PM")) q[[k]] <- q$original_effects[[k]]
                      q$coefficient_diagnostics <- diagnostics
                      q$predictions <- predictions
                      q$native_validation <- list(original = original$native_validation, warm = warm$native_validation)
                      q$event_weight_support <- list(disease_events = sum(rr$event), death_events = sum(rr$death_event), 
                        center = warm$center_support, weights = list(min = min(w), max = max(w), sum = sum(w), 
                          ESS = sum(w)^2/sum(w^2)))
                      q$individual_probability_1e9_diagnostic_passed <- all(vapply(predictions, function(p) isTRUE(p$individual_probability_stability_passed), 
                        logical(1L)))
                      q$individual_probability_1e9_is_acceptance_gate <- FALSE
                      q$individual_probability_1e9_failed_flag_was_not_renamed_or_relaxed <- TRUE
                      warned <- any(vapply(state$warnings, length, 0L) > 0L)
                      q$status <- if (warned) "WARNING_QUALIFIED_CONDITIONAL_ESTIMATE" else "VALID"
                      q$reason <- if (warned) "AUTHORIZED_RECOGNIZED_NUISANCE_WARNING_CORE_AND_NATIVE_PROBABILITY_GUARDS_PASS" else "CLEAN_ORIGINAL_NATIVE_V4_ALL_GUARDS_PASS"
                      q$estimate_family <- "NATIVE_ORIGINAL_FAMILY_NO_IE_SIGN_P_SELECTION"
                      q$compact_prediction_validation <- compact_checks
                      cache_object <- list(schema = "P2_PRIVATE_COMPACT_PREDICTION_CACHE_V1", identity = cell_identity, 
                        original = compact_original, warm = compact_warm, original_effects = q$original_effects, 
                        warm_effects = q$warm_effects, predictions = predictions, diagnostics = diagnostics, compact_validation = compact_checks, 
                        frozen_sources = list(task = group$task_path, task_seal_sha256 = group$task_seal_sha256, 
                          draw_rds_path = group$draw_rds_path, draw_rds_sha256 = group$draw_rds_sha256, A_descriptors_path = group$a_descriptors_path, 
                          A_descriptors_sha256 = group$a_descriptors_sha256), binding = list(actual_draw_hash = q$actual_draw_hash, 
                          target_hash = q$target_hash, sample_hash = q$sample_hash, selection_weight_hash = q$selection_weight_hash, 
                          basis_reference_hash = q$actual_basis_reference_hash), Y_refits_required_for_prediction = FALSE, 
                        contains_subject_model_x_y_or_fitting_environment = FALSE)
                    }
                    rw_assert(identical(rng_before, b16_a_cache_rng()), "P2_WARM_RUNTIME_A_Y_OR_PREDICTION_CHANGED_RNG")
                  }
                }, error = function(e) {
                  q$status <<- ff_w_classify(e, bid)
                  q$reason <<- conditionMessage(e)
                  q$error_class <<- class(e)
                  q$native_singular_error <<- grepl("singular", conditionMessage(e), ignore.case = TRUE)
                  q$no_fallback_after_failure <<- TRUE
                  cache_object <<- NULL
                })
                q$Y_cause_fits_attempted <- state$fits_attempted
                q$native_errors <- state$native_errors
                q$native_warning_gate_exceptions <- state$warning_gates
                q$native_attempt_audit <- state$audit
                q$clone_body_sha256 <- state$clone_body_sha256
                if (length(state$designs) && exists("endpoint", inherits = FALSE) && exists("w", inherits = FALSE)) q$coefficient_diagnostics <- tryCatch(p2_diagnostics(state, 
                  column, endpoint, w), error = function(e) list(error = conditionMessage(e)))
                artifacts <- list()
                if (is.null(cache_object) && length(state$models)) {
                  failed_aggregate <- lapply(names(state$models), function(key) {
                    m <- state$models[[key]]
                    cause <- sub("^(original|warm)_", "", key)
                    list(attempt = key, cause = cause, actual_terms = state$designs[[cause]]$reference$columns, 
                      internal_terms = names(stats::coef(m)), coefficients = stats::coef(m), covariance = stats::vcov(m), 
                      warnings = state$warnings[[key]], controls = state$controls[[key]], iterations = m$iter, 
                      formula = paste(deparse(stats::formula(m)), collapse = " "))
                  })
                  fa <- file.path(private, "FAILED_MODEL_AGGREGATES.private.rds")
                  saveRDS(list(schema = "P2_FAILED_MODEL_AGGREGATES_V1", identity = cell_identity, attempts = failed_aggregate, 
                    native_errors = state$native_errors, reason = q$reason), fa, compress = FALSE, version = 2)
                  q$failed_model_aggregate_cache <- list(relative_path = "private/FAILED_MODEL_AGGREGATES.private.rds", 
                    sha256 = rw_sha(fa), bytes = file.info(fa)$size)
                  artifacts[["private/FAILED_MODEL_AGGREGATES.private.rds"]] <- rw_sha(fa)
                }
                if (!is.null(cache_object)) {
                  raw_cache <- serialize(cache_object, NULL, version = 2)
                  max_cache <- if (is.null(manifest$max_compact_cache_bytes_per_cell)) 32 * 1024^2 else as.numeric(manifest$max_compact_cache_bytes_per_cell)
                  rw_assert(length(raw_cache) <= max_cache, "P2_COMPACT_CACHE_SIZE_CAP_EXCEEDED_PRESERVE_FULL_MODELS_IN_MEMORY_NO_WRITE")
                  cp <- file.path(private, "POSTPROCESS_CACHE.private.rds")
                  saveRDS(cache_object, cp, compress = FALSE, version = 2)
                  q$compact_private_cache <- list(relative_path = "private/POSTPROCESS_CACHE.private.rds", sha256 = rw_sha(cp), 
                    bytes = file.info(cp)$size)
                  artifacts[["private/POSTPROCESS_CACHE.private.rds"]] <- q$compact_private_cache$sha256
                } else if (length(state$models) && !first_failure_models_saved) {
                  max_failure <- if (is.null(manifest$max_representative_failure_model_bytes_per_worker)) 128 * 
                    1024^2 else as.numeric(manifest$max_representative_failure_model_bytes_per_worker)
                  raw_failure <- serialize(state$models, NULL, version = 2)
                  if (length(raw_failure) <= max_failure) {
                    fp <- file.path(private, "REPRESENTATIVE_NEW_FAILURE_MODELS.private.rds")
                    saveRDS(state$models, fp, compress = FALSE, version = 2)
                    first_failure_models_saved <- TRUE
                    q$representative_failure_models <- list(relative_path = "private/REPRESENTATIVE_NEW_FAILURE_MODELS.private.rds", 
                      bytes = file.info(fp)$size, sha256 = rw_sha(fp))
                    artifacts[["private/REPRESENTATIVE_NEW_FAILURE_MODELS.private.rds"]] <- rw_sha(fp)
                  } else q$full_failure_model_retention_skipped_byte_cap <- TRUE
                }
                q$full_successful_subject_models_retained <- FALSE
                rw_json(q, file.path(cell_dir, "RESULT.json"))
                artifacts[["RESULT.json"]] <- rw_sha(file.path(cell_dir, "RESULT.json"))
                rw_json(c(cell_identity, list(complete = TRUE, status = q$status, artifacts_sha256 = artifacts)), 
                  file.path(cell_dir, "COMPLETE.json"))
                group_cases[[epid]] <- list(endpoint_id = epid, relative_directory = basename(cell_dir), status = q$status, 
                  complete_sha256 = rw_sha(file.path(cell_dir, "COMPLETE.json")))
                counts$cases_complete <- counts$cases_complete + 1L
                counts$Y_cause_fits <- counts$Y_cause_fits + state$fits_attempted
                if (identical(q$status, "WARNING_QUALIFIED_CONDITIONAL_ESTIMATE")) counts$conditional_estimates <- counts$conditional_estimates + 
                  1L else if (identical(q$status, "VALID")) counts$clean_VALID_estimates <- counts$clean_VALID_estimates + 
                  1L else if (!startsWith(q$status, "SKIPPED_")) counts$failures <- counts$failures + 1L
                emit("CELL_SEALED", list(group_id = group$group_id, endpoint_id = epid, status = q$status, reason = q$reason, 
                  counts = counts))
                release <- intersect(c("state", "a", "X", "SX", "data", "w", "rr", "endpoint", "design", "changed", 
                  "centres", "warm", "original", "originals", "diagnostics", "compact_original", "compact_warm", 
                  "predictions", "nd", "pred", "cache_object", "raw_cache", "raw_failure", "fit_fn", "audit", "compare_compact", 
                  "p1", "failed_aggregate"), ls(envir = environment()))
                rm(list = release, envir = environment())
                invisible(gc(FALSE))
            }
            rw_check(task_seal_path, group$task_seal_sha256)
            rw_check(draw_path, group$draw_rds_sha256)
            rw_check(descriptor_path, group$a_descriptors_sha256)
            group_result <- c(group_identity, list(schema = "P2_GROUP_RESULT_V1", cases = group_cases, target_masks_identical_all_five = TRUE, 
                production_modified = FALSE, A_fits = 0L, selection_fits = 0L, draws_generated = 0L))
            rw_json(group_result, file.path(group_out, "RESULT.json"))
            group_artifacts <- list(RESULT.json = rw_sha(file.path(group_out, "RESULT.json")))
            for (entry in group_cases) group_artifacts[[file.path(entry$relative_directory, "COMPLETE.json")]] <- entry$complete_sha256
            rw_json(c(group_identity, list(complete = TRUE, artifacts_sha256 = group_artifacts)), file.path(group_out, 
                "COMPLETE.json"))
            counts$groups_complete <- counts$groups_complete + 1L
            ctx$b16_a_task_cache <- NULL
            rm(d, pr, ix, saved, cache, descriptors, envir = environment())
            invisible(gc(FALSE))
            emit("GROUP_SEALED", list(group_id = group$group_id, counts = counts))
        }
        rw_check(manifest_path, manifest_sha)
        for (k in names(source_pins)) rw_check(file.path(ROOT, "f", k), source_pins[[k]])
        for (k in names(runtime_pins)) rw_check(file.path(HERE, k), runtime_pins[[k]])
        summary_dir <- file.path(OUT, paste0("worker_", sprintf("%02d", worker_id), if (canary_mode) "_canary" else ""))
        if (dir.exists(summary_dir)) {
            p2_verified(summary_dir, base_identity)
        } else {
            summary_dir <- rw_newdir(summary_dir)
            summary <- c(base_identity, list(schema = "P2_WORKER_SUMMARY_V1", counts = counts, namespace_initialization = namespace_audit, 
                runtime = runtime, frozen_source_sha256 = as.list(source_pins), recovery_runtime_sha256 = as.list(runtime_pins), 
                old_VALID_overwritten = FALSE, new_draws = 0L, new_A_models = 0L, new_selection_models = 0L, conditional_estimates_are_not_finite_joint_MLE_certificates = TRUE))
            rw_json(summary, file.path(summary_dir, "RESULT.json"))
            rw_json(c(base_identity, list(complete = TRUE, artifacts_sha256 = list(RESULT.json = rw_sha(file.path(summary_dir, 
                "RESULT.json"))))), file.path(summary_dir, "COMPLETE.json"))
        }
        emit("WORKER_COMPLETE", list(counts = counts))
    }), envir = context)
}

