# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

stage <- "LOAD"

event <- 0L

progress <- function(stage, details = list()) {
    event <<- event + 1L
    rec <- c(list(stage = stage, model_id = model, b = b, phase = phase, variant = variant, pid = Sys.getpid(), 
        utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")), details)
    pr_json(rec, file.path(out, "PROGRESS.json"))
    cat(jsonlite::toJSON(rec, auto_unbox = TRUE, digits = 17, null = "null"), "\n")
    flush.console()
}

fractions <- c(ethnicity = 0.5, body = 0.5)

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        args <- commandArgs(TRUE)
        stopifnot(length(args) == 6L)
        root <- normalizePath(args[[1]], mustWork = TRUE)
        phase <- args[[2]]
        model <- args[[3]]
        b <- as.integer(args[[4]])
        variant <- args[[5]]
        out <- args[[6]]
        source(file.path(root, "R/runtime.R"))
        cfg <- pr_read(file.path(root, "ops/runtime.json"))
        pr_assert(phase %in% c("Y", "MC") && variant %in% c("primary", "body_ridge"), "EXACT_REPAIR_PHASE_AND_AUTHORIZED_METHOD_REQUIRED")
        pr_assert(!dir.exists(out), "NEW_EXCLUSIVE_ATTEMPT_DIRECTORY_REQUIRED")
        dir.create(out, recursive = TRUE)
        start <- proc.time()[["elapsed"]]
        stage <- "LOAD"
        event <- 0L
        progress <- function(stage, details = list()) {
            event <<- event + 1L
            rec <- c(list(stage = stage, model_id = model, b = b, phase = phase, variant = variant, pid = Sys.getpid(), 
                utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")), details)
            pr_json(rec, file.path(out, "PROGRESS.json"))
            cat(jsonlite::toJSON(rec, auto_unbox = TRUE, digits = 17, null = "null"), "\n")
            flush.console()
        }
        result <- list(status = "RUNNING", model_id = model, b = b, phase = phase, variant = variant)
        setting <- Sys.getenv("PRY_B0_SENSITIVITY", "")
        fractions <- c(ethnicity = 0.5, body = 0.5)
        if (nzchar(setting)) {
            allowed <- c(ethnicity_low = 1/3, ethnicity_high = 2/3, body_low = 1/3, body_high = 2/3)
            pr_assert(b == 0L && setting %in% names(allowed), "ONLY_PRESPECIFIED_B0_PENALTY_SENSITIVITY")
            block <- if (startsWith(setting, "ethnicity")) "ethnicity" else "body"
            pr_assert(block != "body" || variant == "body_ridge", "BODY_SENSITIVITY_ONLY_IF_AUTHORIZED_BODY_FALLBACK_ADOPTED")
            fractions[block] <- allowed[[setting]]
        }
        result$analysis_role <- if (nzchar(setting)) paste0("B0_SENSITIVITY_", setting) else "MAIN_FIXED202"
        tryCatch({
            progress(stage)
            rt <- pr_runtime(cfg)
            methods <- rt$methods
            engine <- rt$engine
            bundle <- pr_read(file.path(root, "ops/cache_bundle.json"))
            jobs <- Filter(function(j) j$ticket$identity$model_id == model && j$ticket$identity$b == b, bundle)
            pr_assert(length(jobs) == 1L, "EXACT_ONE_FROZEN_JOB_REQUIRED")
            job <- jobs[[1]]
            x <- pr_data(job, rt)
            id <- x$identity
            pr_json(list(original = job$ticket, new_method = variant, binding = x$binding, original_unchanged = TRUE), 
                file.path(out, "BINDING.json"))
            reference <- readRDS(file.path(cfg$state_root, "references", paste0(model, ".rds")))
            pieces <- list()
            validation <- list()
            if (phase == "Y") {
                for (endpoint in c("IHD", "T2D")) {
                  risk <- engine$ff_y_risk(x$risks[[endpoint]], 11)
                  parts <- list()
                  for (cause in c("disease", "death")) {
                    stage <- paste("FIT", endpoint, cause, sep = "_")
                    progress(stage)
                    piece <- methods$pry_fit_piece(x$X, risk$time, risk$status, x$d$center, x$w, reference, model_id = model, 
                      b = b, endpoint = endpoint, cause = cause, body_penalty = variant == "body_ridge", ethnic_df_fraction = fractions[["ethnicity"]], 
                      body_df_fraction = fractions[["body"]])
                    pr_rds(piece, file.path(out, paste0(endpoint, "_", cause, ".private.rds")))
                    pr_json(piece$diagnostics, file.path(out, paste0(endpoint, "_", cause, "_DIAGNOSTICS.json")))
                    parts[[cause]] <- piece
                    pr_assert(isTRUE(piece$diagnostics$passed), paste("Y_DIAGNOSTIC_FAILURE", endpoint, cause, 
                      paste(piece$diagnostics$warnings, collapse = " | ")))
                  }
                  pieces[[endpoint]] <- parts
                  stage <- paste0("VALIDATE_", endpoint)
                  progress(stage)
                  compiled <- methods$pry_compile_endpoint(parts, engine)
                  ii <- sort(unique(c(match(unique(x$d$center), x$d$center), match(unique(x$d$ethnicity), x$d$ethnicity))))
                  cc <- x$d[ii, , drop = FALSE]
                  M <- cc[, x$columns, drop = FALSE]
                  cases <- list(list(label = "actual_observed_rows", kind = "observed", X = x$X[ii, , drop = FALSE], 
                    center = cc$center, ethnicity = cc$ethnicity))
                  for (B in 0:1) cases[[length(cases) + 1L]] <- list(label = paste0("fixed_observed_M_B", B), kind = "counterfactual", 
                    X = x$design(B, M, cc), center = cc$center, ethnicity = cc$ethnicity)
                  validation[[endpoint]] <- methods$pry_validate_prediction(compiled, cases, engine, expected_ethnicity = unique(as.character(x$d$ethnicity)))
                  pr_json(validation[[endpoint]], file.path(out, paste0(endpoint, "_NATIVE_PREDICTION.json")))
                }
                payload <- list(pieces = pieces, binding = x$binding, variant = variant, model = model, b = b, 
                  reference_sha256 = reference$reference_sha256, validation = validation, sensitivity_setting = setting)
                record <- pr_rds(payload, file.path(out, "VALIDATED_Y.private.rds"))
                result$Y <- record
                result$status <- "SUCCESS"
                result$validation <- validation
            } else {
                gate <- pr_read(file.path(cfg$run_root, variant, "ALL_Y_GATE.json"))
                pr_assert(isTRUE(gate$passed) && gate$tasks == 202, "ALL_202_Y_BEFORE_EXPENSIVE_MC_REQUIRED")
                if (nzchar(setting)) {
                  path <- Sys.getenv("PRY_B0_SENSITIVITY_Y")
                  pr_assert(file.exists(path), "ACTUAL_VALIDATED_SENSITIVITY_Y_REQUIRED")
                  rec <- list(path = path, sha256 = Sys.getenv("PRY_B0_SENSITIVITY_Y_SHA256"))
                } else rec <- gate$jobs[[paste0(model, ":", b)]]
                pr_assert(!is.null(rec) && identical(pr_sha(rec$path), rec$sha256), "VALIDATED_NEW_Y_CACHE_CHANGED")
                saved <- readRDS(rec$path)
                pr_assert(identical(saved$binding, x$binding) && saved$variant == variant && saved$model == model && 
                  saved$b == b && identical(saved$sensitivity_setting, setting), "EXACT_NEW_Y_BINDING_CHANGED_NO_REFIT")
                pieces <- saved$pieces
                compiled <- lapply(pieces, methods$pry_compile_endpoint, engine = engine)
                for (e in names(compiled)) {
                  pr_assert(isTRUE(saved$validation[[e]]$passed), "SAVED_PREDICTION_NOT_VALIDATED")
                  compiled[[e]]$validation_state$passed <- TRUE
                  compiled[[e]]$validation_state$receipt <- saved$validation[[e]]
                }
                stage <- "REUSE_A"
                progress(stage)
                Apath <- file.path(cfg$state_root, "a_cache", model, sprintf("B%04d", b), "FRESH_A.private.rds")
                reused <- methods$s44pr_reuse_A(Apath, job, x$A_context, engine, methods, cfg$source_root)
                A <- reused$A
                pr_json(reused$provenance, file.path(out, "A_REUSE.json"))
                draw_joint <- function(b, C, uniforms, requested_domain = NULL) engine$c39_draw_a(A, C, b, uniforms[, 
                  seq_along(x$members), drop = FALSE])[, x$columns, drop = FALSE]
                predict_risk <- function(b, M, C) {
                  X <- x$design(b, M, C)
                  do.call(cbind, lapply(compiled, function(f) methods$pry_predict_matrix(f, X, C$center, engine)[, 
                    "disease"]))
                }
                support <- function(M, C, context) {
                  engine$c39_support(M, Filter(function(m) m$column %in% names(M), x$members))
                  TRUE
                }
                ii <- sort(unique(c(match(unique(x$d$center), x$d$center), match(unique(x$d$ethnicity), x$d$ethnicity))))
                cc <- x$d[ii, , drop = FALSE]
                u <- matrix(0.5, nrow(cc), 15)
                for (B in 0:1) {
                  M <- draw_joint(B, cc, u)
                  support(M, cc, "GENERATED_PROBE")
                  for (e in names(compiled)) {
                    cases <- list(list(label = "actual", kind = "observed", X = x$X[ii, , drop = FALSE], center = cc$center, 
                      ethnicity = cc$ethnicity))
                    for (Bcheck in 0:1) cases[[length(cases) + 1L]] <- list(label = paste0("generated_B", B, "_Y", 
                      Bcheck), kind = "counterfactual", X = x$design(Bcheck, M, cc), center = cc$center, ethnicity = cc$ethnicity)
                    methods$pry_validate_prediction(compiled[[e]], cases, engine, unique(as.character(x$d$ethnicity)))
                  }
                }
                checkpoint <- NULL
                checkpoint_root <- file.path(cfg$state_root, "mc_checkpoints", variant, if (nzchar(setting)) paste0(model, 
                  "_sensitivity_", setting) else model, sprintf("B%04d", b))
                dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
                cpfiles <- list.files(checkpoint_root, pattern = "^MC_[0-9]+[.]private[.]rds$", full.names = TRUE)
                if (length(cpfiles)) checkpoint <- readRDS(tail(sort(cpfiles), 1))
                save_checkpoint <- function(s) pr_rds(s, file.path(checkpoint_root, sprintf("MC_%02d.private.rds", 
                  s$completed_mc)))
                binding <- list(panel_id = paste0(id$model_id, "_", variant), bootstrap_id = b, sample_hash = id$sampling_hash, 
                  population_hash = id$population_hash, model_spec_hash = pr_hash(list(original = id$model_spec, 
                    Y = variant, validated_Y_sha256 = rec$sha256, reference = reference$reference_sha256)), input_hash = id$input_hash, 
                  rng_spec_hash = id$rng_spec_hash)
                RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
                stage <- "K50_INTEGRATION"
                progress(stage)
                integral <- methods$nmr_product_stream(x$d, x$w, x$domains, draw_joint, predict_risk, support, 
                  uniform_dim = 15L, n_mc = 50L, block_size = 512L, seed = as.integer(id$seed), allocation_order = names(x$domains), 
                  include_reverse = TRUE, binding = binding, checkpoint = checkpoint, save_checkpoint = save_checkpoint, 
                  progress = function(v) progress(v$stage, v[setdiff(names(v), "stage")]))
                for (name in c("components", "paired", "scenarios")) write.csv(integral[[name]], file.path(out, 
                  paste0(toupper(name), ".csv")), row.names = FALSE)
                pr_rds(integral, file.path(out, "INTEGRATION.private.rds"))
                result$status <- "SUCCESS"
                result$components <- integral$components
                result$paired <- integral$paired
                result$risks <- integral$scenarios
                result$audit <- integral$audit
                result$Y_refits <- 0L
                result$A_refits <- 0L
            }
            result$binding <- x$binding
            result$pair_binding <- x$pair_binding
            result$N <- nrow(x$d)
        }, error = function(e) {
            result$status <<- "FAILED"
            result$reason <<- conditionMessage(e)
            result$failed_stage <<- stage
        })
        result$pid <- Sys.getpid()
        result$elapsed_seconds <- proc.time()[["elapsed"]] - start
        result$utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        pr_json(result, file.path(out, "result.json"))
        progress(result$status, list(reason = result$reason, elapsed_seconds = result$elapsed_seconds))
        if (result$status != "SUCCESS") quit(status = 2)
    }), envir = context)
}

