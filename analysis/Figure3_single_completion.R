# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

new_json <- function(z, p) {
    stopifnot(!file.exists(p))
    jsonlite::write_json(z, p, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
}

emit <- function(z) {
    cat(jsonlite::toJSON(z, auto_unbox = TRUE, digits = 17, na = "null", null = "null"), "\n")
    flush.console()
}

assemble <- function(fits, baselines, dd, risk, w) {
    names(fits) <- names(baselines) <- c("disease", "death")
    grid <- sort(unique(unlist(lapply(baselines, `[[`, "time"))))
    cumulative <- vapply(baselines, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
    hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
        , drop = FALSE])
    stopifnot(all(is.finite(hazard)), all(hazard >= 0))
    keep <- rowSums(hazard) > 0
    ans <- list(fit = fits, engine = list(time = grid[keep], hazard = hazard[keep, , drop = FALSE], coefficients = vapply(fits, 
        stats::coef, numeric(ncol(dd$X)))), matrix_reference = dd$reference, tau_report = 11, endpoint = ff_y_risk(risk, 
        11), binding = list(weight = w, cluster = NULL))
    class(ans) <- c("ff_y_matrix_fit", "ff_y_fit")
    ans
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        args <- commandArgs(TRUE)
        stopifnot(length(args) == 4L)
        ROOT <- normalizePath(args[1], mustWork = TRUE)
        OUT <- normalizePath(args[2], mustWork = TRUE)
        PLAN <- normalizePath(args[3], mustWork = TRUE)
        SLOT <- args[4]
        for (package in c("jsonlite", "digest", "survival", "nnet")) suppressPackageStartupMessages(library(package, 
            character.only = TRUE))
        for (file in c("10.prepare.R", "20.a_reuse.R", "30.outcome.R", "31.aj_backend.R", "40.worker.R")) source(file.path(ROOT, 
            "f", file))
        new_json <- function(z, p) {
            stopifnot(!file.exists(p))
            jsonlite::write_json(z, p, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
        }
        emit <- function(z) {
            cat(jsonlite::toJSON(z, auto_unbox = TRUE, digits = 17, na = "null", null = "null"), "\n")
            flush.console()
        }
        plan <- jsonlite::fromJSON(PLAN, simplifyVector = FALSE)
        stopifnot(plan$max_bootstrap_id == 100L, isTRUE(plan$no_old_output_writes))
        tasks <- Filter(function(z) identical(z$host, Sys.getenv("FF_HOST")) && identical(z$slot, SLOT), plan$tasks)
        backend <- ff_aj_initialize(ROOT, Sys.getenv("FF_WORK"), require_compiled = TRUE)
        ctx <- ff_load(ROOT)
        oldrun <- ctx$run
        ctx$run <- OUT
        ctx$ff_run <- OUT
        ctx$run_id <- plan$run_id
        ctx$allow_missing_a_recompute <- TRUE
        ctx$ff_draw_cache <- new.env(parent = emptyenv())
        ctx$ff_bridge_cache <- new.env(parent = emptyenv())
        ctx$monitor <- new.env(parent = emptyenv())
        ctx$monitor$selection_weight_hashes <- list()
        ctx$monitor$endpoint_status <- list()
        ctx$monitor$process_start <- as.numeric(Sys.time())
        state_dir <- file.path(OUT, "ops", Sys.getenv("FF_HOST"))
        dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
        Sys.setenv(V5_WORKER_STATE = file.path(state_dir, paste0(SLOT, "_STATE.json")), V5_WORKER_EVENTS = file.path(state_dir, 
            paste0(SLOT, "_EVENTS.jsonl")), V5_WORKER_SLOT = SLOT)
        for (sub in c("tasks", "private/new_A")) dir.create(file.path(OUT, sub), recursive = TRUE, showWarnings = FALSE)
        assemble <- function(fits, baselines, dd, risk, w) {
            names(fits) <- names(baselines) <- c("disease", "death")
            grid <- sort(unique(unlist(lapply(baselines, `[[`, "time"))))
            cumulative <- vapply(baselines, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
            hazard <- rbind(cumulative[1L, , drop = FALSE], cumulative[-1L, , drop = FALSE] - cumulative[-nrow(cumulative), 
                , drop = FALSE])
            stopifnot(all(is.finite(hazard)), all(hazard >= 0))
            keep <- rowSums(hazard) > 0
            ans <- list(fit = fits, engine = list(time = grid[keep], hazard = hazard[keep, , drop = FALSE], coefficients = vapply(fits, 
                stats::coef, numeric(ncol(dd$X)))), matrix_reference = dd$reference, tau_report = 11, endpoint = ff_y_risk(risk, 
                11), binding = list(weight = w, cluster = NULL))
            class(ans) <- c("ff_y_matrix_fit", "ff_y_fit")
            ans
        }
        emit(list(phase = "READY", pid = Sys.getpid(), host = Sys.getenv("FF_HOST"), slot = SLOT, units = length(tasks), 
            roots = list(code = ROOT, input = oldrun, output = OUT), method_hash = ctx$ff_method_hash, input_sha = ctx$input_lock_sha256))
        for (task in tasks) {
            id <- task$indicator_id
            bid <- as.integer(task$bootstrap_id)
            stopifnot(id %in% ctx$scope$indicator_id, bid %in% 0:100)
            unit <- file.path(OUT, "tasks", paste0(id, "_B", sprintf("%04d", bid)))
            if (file.exists(file.path(unit, "COMPLETE.json"))) {
                emit(list(phase = "REUSE_NEW_SEALED_UNIT", id = id, bid = bid))
                next
            }
            stopifnot(!dir.exists(unit), dir.create(unit))
            new_json(list(host = Sys.getenv("FF_HOST"), slot = SLOT, pid = Sys.getpid(), indicator_id = id, bootstrap_id = bid, 
                started_utc = v5_utc()), file.path(unit, "CLAIM.json"))
            ctx$monitor$indicator_id <- id
            ctx$monitor$bootstrap_id <- bid
            ctx$monitor$attempt_id <- basename(unit)
            ctx$monitor$package_started_epoch <- as.numeric(Sys.time())
            ctx$monitor$package_started_monotonic <- proc.time()[["elapsed"]]
            actual <- ff_w_draw(ctx, bid)
            ctx$actual_draw_hash <- actual$receipt$draw_hash
            rows <- list()
            for (epid in unlist(task$endpoints)) {
                started <- proc.time()[["elapsed"]]
                z <- list(indicator_id = id, endpoint_id = epid, bootstrap_id = bid, host = Sys.getenv("FF_HOST"), 
                  slot = SLOT, status = "NOT_ESTIMABLE", reason = "", TE = NULL, DE = NULL, IE = NULL, PM = NULL, 
                  R00 = NULL, R10 = NULL, R11 = NULL, draw = actual$receipt, method_hash = ctx$ff_method_hash, 
                  source_hash = ctx$source_hash, input_hash = ctx$input_lock_sha256, report_policy = "WARNING_AWARE_EXPLORATORY_FIXED_B100_NO_FINITE_JOINT_MLE_CLAIM", 
                  model_formula_changed = FALSE, estimator_changed = FALSE, seed_changed = FALSE, probability_clipping = FALSE)
                tryCatch({
                  ca <- ctx$features[[id]]
                  ep <- ctx$endpoints[[epid]]
                  pr <- actual$ix[(ctx$eligibility[[id]] & ep$eligibility)[actual$ix]]
                  d <- ctx$d_full[pr, , drop = FALSE]
                  risk <- ep$risk[pr, , drop = FALSE]
                  base <- unique(c(ctx$cfg$base_covariates, ca$technical_covariates))
                  ff_w_stage(ctx, "A_CACHE_REUSE_OR_EXPLICIT_MISSING_OBJECT_COMPLETION", epid)
                  a <- ff_a_scenarios(ctx, id, bid, epid, pr, d, base)
                  z$A_reuse <- as.list(a$ledger[1L, , drop = FALSE])
                  cc <- which(a$mi$complete)
                  rr <- risk[cc, , drop = FALSE]
                  data <- a$data
                  z$N_target <- nrow(d)
                  z$N_complete <- nrow(data)
                  z$events <- sum(rr$event)
                  z$deaths <- sum(rr$death_event)
                  z$target_hash <- ff_a_hash(pr)
                  z$sample_hash <- ff_a_hash(pr[cc])
                  z$selection_weight_hash <- ff_a_hash(format(a$mi$w, digits = 17, scientific = TRUE))
                  ap <- file.path(OUT, "private/new_A", paste0(id, "_B", sprintf("%04d", bid), "_", epid, ".rds"))
                  stopifnot(!file.exists(ap))
                  aa <- list(ledger = a$ledger, draw = actual$receipt, selection_weight = a$mi$w, selection_coefficients = if (!is.null(a$mi$fit)) coef(a$mi$fit) else a$mi$saved_coefficients, 
                    alpha = a$alpha, g0 = a$g0, g1 = a$g1, reference_basis = a$reference_basis, target_hash = z$target_hash, 
                    sample_hash = z$sample_hash, input_hash = z$input_hash)
                  saveRDS(aa, ap, compress = FALSE)
                  z$new_A_checkpoint_sha256 <- ff_a_sha(ap)
                  if (!is.null(a$alpha)) X <- ctx$engine$b3_design(data, ca$column, base, a$reference_basis)$X else {
                    recipe <- ctx$discrete$c40_recipe(data, ca$column, base)
                    X <- ctx$discrete$c40_matrix(recipe, seq_len(nrow(data)), M = data[[ca$column]])
                  }
                  dd <- ff_y_matrix(X)
                  endpoint <- ff_y_risk(rr, 11)
                  work <- dd$work
                  work$.ff_time <- endpoint$time
                  work$.ff_disease <- as.integer(endpoint$status == 1L)
                  work$.ff_death <- as.integer(endpoint$status == 2L)
                  work$.ff_weight <- a$mi$w
                  core <- which(colnames(dd$X) == "B" | grepl(ca$column, colnames(dd$X), fixed = TRUE))
                  stopifnot(length(core) >= 2L)
                  originals <- warms <- baselines <- warm_baselines <- list()
                  audit <- list()
                  ff_w_stage(ctx, "FULL_FOLLOWUP_COX_SAME_FORMULA", epid)
                  for (cause in 1:2) {
                    got <- ff_y_fit_native_one(work, dd$reference$safe_names, cause, FALSE)
                    m <- got$value
                    stopifnot(all(is.finite(coef(m))), all(is.finite(vcov(m))), all(diag(vcov(m))[core] > 0))
                    b <- ff_y_native_baseline(m, dd$reference$safe_names)
                    warm <- ff_y_capture(survival::coxph(formula(m), data = work, weights = .ff_weight, robust = FALSE, 
                      ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, model = TRUE, x = TRUE, 
                      y = TRUE, init = coef(m), control = survival::coxph.control(eps = 1e-11, iter.max = 50L)))
                    wfit <- warm$value
                    stopifnot(all(is.finite(coef(wfit))), all(is.finite(vcov(wfit))))
                    wb <- ff_y_native_baseline(wfit, dd$reference$safe_names)
                    text <- paste(got$warnings, collapse = " | ")
                    hits <- regmatches(text, gregexpr("variable[[:space:]]+[0-9]+", text))[[1L]]
                    indices <- as.integer(sub("variable[[:space:]]+", "", hits))
                    indices <- indices[!is.na(indices)]
                    core_hit <- intersect(indices, core)
                    audit[[cause]] <- list(cause = if (cause == 1L) "disease" else "death", warnings = got$warnings, 
                      warning_terms = colnames(dd$X)[indices], warning_core_terms = colnames(dd$X)[core_hit], warm_warnings = warm$warnings, 
                      core_terms = colnames(dd$X)[core], core_beta = unname(coef(m)[core]), core_SE = sqrt(diag(vcov(m)))[core], 
                      core_beta_max_change_over_original_SE = max(abs(coef(m)[core] - coef(wfit)[core])/sqrt(diag(vcov(m)))[core]), 
                      all_beta_max_abs_change = max(abs(coef(m) - coef(wfit))), core_information_min_eigenvalue = min(eigen(solve(vcov(m)[core, 
                        core, drop = FALSE]), symmetric = TRUE, only.values = TRUE)$values), native_baseline_public_validation = b$max_abs_basehaz_vs_survfit)
                    originals[[cause]] <- m
                    warms[[cause]] <- wfit
                    baselines[[cause]] <- b
                    warm_baselines[[cause]] <- wb
                  }
                  z$coefficient_diagnostics <- audit
                  stable <- all(vapply(audit, function(q) q$core_beta_max_change_over_original_SE <= 0.01 && q$core_information_min_eigenvalue > 
                    0 && !length(q$warning_core_terms), logical(1L)))
                  if (!stable) stop("CORE_COEFFICIENT_NOT_IDENTIFIABLE_OR_NOT_STABLE_OR_TARGET_WARNING")
                  object <- assemble(originals, baselines, dd, rr, a$mi$w)
                  warm_object <- assemble(warms, warm_baselines, dd, rr, a$mi$w)
                  ff_w_stage(ctx, "FULL_FOLLOWUP_CIF_INTEGRAL", epid)
                  observed <- ff_predict_y_matrix(object, dd$X, 11)
                  integral <- ff_w_integrate(ctx, a, base, object, ca$column)
                  warm_integral <- ff_w_integrate(ctx, a, base, warm_object, ca$column)
                  z$observed_probability_range <- apply(observed, 2L, range)
                  z$observed_max_closure <- max(abs(rowSums(observed) - 1))
                  z$risk_max_abs_movement <- max(abs(integral$risks - warm_integral$risks))
                  for (k in names(integral$risks)) z[[k]] <- unname(integral$risks[[k]])
                  z$TE <- z$R11 - z$R00
                  z$DE <- z$R10 - z$R00
                  z$IE <- z$R11 - z$R10
                  z$PM <- if (z$TE != 0) z$IE/z$TE else NULL
                  z$IE_abs_movement <- abs(z$IE - unname(warm_integral$risks[["R11"]] - warm_integral$risks[["R10"]]))
                  z$TE_DE_IE_closure <- z$TE - z$DE - z$IE
                  stopifnot(abs(z$TE_DE_IE_closure) <= 1e-12, all(is.finite(integral$risks)), all(integral$risks >= 
                    0 & integral$risks <= 1))
                  z$status <- "VALID_WARNED_EXPLORATORY"
                  z$reason <- "LEGAL_PREDICTIONS_STABLE_TARGET_COEFFICIENTS_NUISANCE_WARNINGS_RETAINED"
                  z$core_stable <- TRUE
                  rm(originals, warms, baselines, warm_baselines, object, warm_object, observed, integral, warm_integral, 
                    a, aa, X, dd, work)
                }, error = function(e) {
                  z$reason <<- conditionMessage(e)
                  z$status <<- "NOT_ESTIMABLE_OR_EXECUTION_FAILURE_RETAINED"
                })
                z$seconds <- proc.time()[["elapsed"]] - started
                rows[[epid]] <- z
                new_json(z, file.path(unit, paste0("RESULT_", epid, ".json")))
                emit(list(phase = "SAVED", pid = Sys.getpid(), host = Sys.getenv("FF_HOST"), slot = SLOT, id = id, 
                  bid = bid, endpoint = epid, status = z$status, reason = z$reason, seconds = z$seconds))
                invisible(gc(FALSE))
            }
            new_json(list(indicator_id = id, bootstrap_id = bid, rows = rows, draw = actual$receipt), file.path(unit, 
                "RESULTS.json"))
            new_json(list(complete = TRUE, indicator_id = id, bootstrap_id = bid, draw_hash = actual$receipt$draw_hash, 
                sha256 = ff_a_sha(file.path(unit, "RESULTS.json")), finished_utc = v5_utc()), file.path(unit, "COMPLETE.json"))
        }
        new_json(list(complete = TRUE, slot = SLOT, host = Sys.getenv("FF_HOST"), pid = Sys.getpid(), units = length(tasks), 
            finished_utc = v5_utc()), file.path(state_dir, paste0(tools::file_path_sans_ext(basename(PLAN)), "_", 
            SLOT, "_DONE.json")))
    }), envir = context)
}

