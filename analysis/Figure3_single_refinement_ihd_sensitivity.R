# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

rows <- list()

details <- list()

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        args <- commandArgs(TRUE)
        stopifnot(length(args) == 3L, Sys.getenv("FF16_MODE") == "SMOKE", Sys.getenv("FF_REAL_EXECUTION_AUTHORIZED") == 
            "YES_FOR_THIS_RUN")
        for (p in c("jsonlite", "digest", "survival", "nnet")) suppressPackageStartupMessages(library(p, character.only = TRUE))
        oldroot <- args[1]
        out <- args[2]
        reference_path <- args[3]
        stopifnot(!dir.exists(out), dir.exists(dirname(out)), grepl("/border_fullfollowup_b1000_16_20261028_v2/run_20261028_deployment_v2/estimator/", 
            out, fixed = TRUE))
        stopifnot(dir.create(out))
        work <- file.path(out, "work")
        dir.create(work)
        for (k in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")) do.call(Sys.setenv, 
            setNames(list("1"), k))
        for (prefix in c("FF", "V5", "B3", "C40")) do.call(Sys.setenv, setNames(list(work, "M_validation", file.path(out, 
            paste0(prefix, "_state.json")), file.path(out, paste0(prefix, "_events.jsonl"))), paste0(prefix, c("_WORK", 
            "_WORKER_SLOT", "_WORKER_STATE", "_WORKER_EVENTS"))))
        for (f in c("10.prepare.R", "20.a_reuse.R", "30.outcome.R", "31.aj_backend.R", "40.worker.R")) source(file.path(oldroot, 
            "f", f))
        ff_aj_initialize(oldroot, work, require_compiled = TRUE)
        ctx <- ff_load(oldroot)
        ctx$run <- dirname(out)
        ctx$ff_run <- ctx$run
        ctx$ff_draw_cache <- new.env(parent = emptyenv())
        ctx$ff_bridge_cache <- new.env(parent = emptyenv())
        ctx$monitor$worker_state_path <- Sys.getenv("V5_WORKER_STATE")
        ctx$monitor$worker_events_path <- Sys.getenv("V5_WORKER_EVENTS")
        source(file.path(dirname(reference_path), "60.ihd_pooling.R"))
        references <- jsonlite::fromJSON(reference_path, simplifyVector = FALSE)
        rows <- list()
        counts <- list()
        details <- list()
        mapping <- NULL
        began <- proc.time()[["elapsed"]]
        plan <- switch(Sys.getenv("FF_M_TEST_PLAN"), MINI_B99_THREE = data.frame(id = c("clinical__ldl", "nmr__nmr_23404", 
            "smoking"), bid = 99L), WIN_FINAL_ONE = data.frame(id = "clinical__ldl", bid = 100L), MINI_FINAL_ONE = data.frame(id = "clinical__ldl", 
            bid = 99L), expand.grid(bid = c(0L, 99L, 100L), id = c("clinical__ldl", "nmr__nmr_23404", "smoking"), 
            stringsAsFactors = FALSE))
        jsonlite::write_json(list(pid = Sys.getpid(), planned_endpoint_branch_tasks = nrow(plan), prior_V1_endpoint_tests = 14, 
            cumulative_max = 26, formal_started = FALSE, models_not_counted_in_formal_B = TRUE), file.path(out, 
            "START.json"), auto_unbox = TRUE, pretty = TRUE)
        for (i in seq_len(nrow(plan))) {
            id <- plan$id[i]
            bid <- as.integer(plan$bid[i])
            all_warnings <- character()
            key <- paste(id, bid, sep = "|")
            prior <- references[[key]]
            started <- proc.time()[["elapsed"]]
            q <- list(indicator_id = id, endpoint_id = "411", bootstrap_id = bid, branch = "M", implementation = "survival::coxph", 
                ties = "efron", status = "NOT_RUN", reason = "", warning = "", changed_estimator = TRUE, actual_y_fits_started = FALSE)
            tryCatch(withCallingHandlers({
                actual <- ff_w_draw(ctx, bid)
                ctx$actual_draw_hash <- actual$receipt$draw_hash
                ep <- ctx$endpoints[["411"]]
                ca <- ctx$features[[id]]
                base <- unique(c(ctx$cfg$base_covariates, ca$technical_covariates))
                target <- ctx$eligibility[[id]] & ep$eligibility
                pr <- actual$ix[target[actual$ix]]
                d <- ctx$d_full[pr, , drop = FALSE]
                risk <- ep$risk[pr, , drop = FALSE]
                a <- ff_a_scenarios(ctx, id, bid, "411", pr, d, base)
                rr <- risk[which(a$mi$complete), , drop = FALSE]
                dd <- a$data
                q$draw_hash <- actual$receipt$draw_hash
                q$sample_hash <- v5_hash(pr[a$mi$complete])
                q$selection_weight_hash <- v5_hash(format(a$mi$w, digits = 17, scientific = TRUE))
                q$A_origin <- a$ledger$A_origin
                for (nm in c("draw_hash", "sample_hash", "selection_weight_hash")) {
                  expected <- prior[[nm]]
                  if (is.null(expected) && !is.null(prior$a_reuse)) expected <- prior$a_reuse[[nm]]
                  if (!is.null(expected) && nzchar(expected) && !identical(q[[nm]], expected)) stop(paste0("FROZEN_REFERENCE_BINDING_MISMATCH:", 
                    nm, ":expected=", expected, ":actual=", q[[nm]]), call. = FALSE)
                }
                lev <- levels(dd$ethnicity)
                mapped <- ifelse(lev %in% c("Asian or Asian British", "Chinese"), "Asian_or_Chinese_pooled", lev)
                stopifnot(all(c("Asian or Asian British", "Chinese") %in% lev))
                mapping <- data.frame(source_value = lev, target_value = mapped, model = "IHD_death_Cox_only", 
                  original_field = "ethnicity", new_field = "ethnicity_death_pooled", reason = "USER_AUTHORIZED_FIXED_MAP_NOT_OFFICIAL_RECLASSIFICATION", 
                  version = "IHD_DEATH_POOL_V2", official_source = "https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=1001")
                for (z in lev) counts[[length(counts) + 1L]] <- data.frame(indicator_id = id, bid = bid, source_level = z, 
                  pooled_level = mapped[match(z, lev)], N = sum(dd$ethnicity == z), disease = sum(rr$event[dd$ethnicity == 
                    z]), death = sum(rr$death_event[dd$ethnicity == z]))
                if (!is.null(a$alpha)) X <- ctx$engine$b3_design(dd, ca$column, base, a$reference_basis)$X else {
                  recipe <- ctx$discrete$c40_recipe(dd, ca$column, base)
                  X <- ctx$discrete$c40_matrix(recipe, seq_len(nrow(dd)), M = dd[[ca$column]])
                }
                q$actual_y_fits_started <- TRUE
                fit <- ff_m_fit(X, rr, a$mi$w, tau_report = ctx$tau_report)
                integral <- ff_w_integrate(ctx, a, base, fit, ca$column)
                for (nm in names(integral$risks)) q[[nm]] <- unname(integral$risks[[nm]])
                q$TE <- q$R11 - q$R00
                q$DE <- q$R10 - q$R00
                q$IE <- q$R11 - q$R10
                q$PM <- if (q$TE == 0) NA_real_ else q$IE/q$TE
                q$closure <- q$TE - q$DE - q$IE
                q$N <- nrow(dd)
                q$events <- sum(rr$event)
                q$deaths <- sum(rr$death_event)
                observed <- ff_predict_y_matrix(fit, X, ctx$tau_report)
                q$observed_death_CIF_weighted_mean <- weighted.mean(observed[, "death"], a$mi$w)
                q$observed_death_CIF_min <- min(observed[, "death"])
                q$observed_death_CIF_max <- max(observed[, "death"])
                q$noninteger_weights <- fit$weight$noninteger
                q$weight_min <- fit$weight$min
                q$weight_max <- fit$weight$max
                q$weight_ESS <- fit$weight$ESS
                q$native_LP_max_abs <- fit$native_validation$lp_max_abs
                q$basehaz_survfit_max_abs <- fit$native_validation$basehaz_survfit_max_abs
                q$probability_closure <- fit$native_validation$three_state_closure
                q$monotonicity_violation <- fit$native_validation$monotonicity_violation
                q$CIF_valid <- TRUE
                q$finite_coefficients <- all(is.finite(fit$coefficients$coefficient))
                q$warning <- paste(fit$warnings, collapse = " | ")
                q$prior_status <- prior$status
                for (nm in c("R00", "R10", "R11", "TE", "DE", "IE")) q[[paste0("difference_vs_original_", nm)]] <- if (!is.null(prior[[nm]]) && 
                  is.finite(prior[[nm]])) q[[nm]] - prior[[nm]] else NA_real_
                q$status <- "VALID"
                details[[key]] <- list(coefficients = fit$coefficients, native_checks = fit$native_validation, 
                  weights = fit$weight, A_ledger = a$ledger, actual_draw = actual$receipt, original_unchanged_design_sha = v5_hash(X))
                rm(fit, X, observed, a)
                gc(FALSE)
            }, warning = function(w) {
                all_warnings <<- c(all_warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }), error = function(e) {
                q$status <<- "BLOCKED"
                q$reason <<- conditionMessage(e)
            })
            q$execution_warning_all <- paste(unique(all_warnings), collapse = " | ")
            q$seconds <- proc.time()[["elapsed"]] - started
            q$peak_rss_gib <- ctx$engine$b3_rss()
            rows[[key]] <- q
            ff_w_json(q, file.path(out, paste0(gsub("[|]", "__", key), "_RESULT.json")))
            cat("V2_M_CASE ", jsonlite::toJSON(q, auto_unbox = TRUE, na = "null", digits = 17), "\n", sep = "")
            flush.console()
        }
        ff_w_csv(rows, file.path(out, "M_SENSITIVITY.csv"))
        ff_w_json(details, file.path(out, "M_DETAILS.json"))
        if (length(counts)) write.csv(do.call(rbind, counts), file.path(out, "ETHNICITY_EVENT_COUNTS.csv"), row.names = FALSE)
        if (!is.null(mapping)) write.csv(mapping, file.path(out, "ETHNICITY_POOLING_MAP.csv"), row.names = FALSE)
        ff_w_json(list(complete = TRUE, pid = Sys.getpid(), endpoint_branch_attempts = nrow(plan), valid = sum(vapply(rows, 
            function(z) z$status == "VALID", FALSE)), blocked = sum(vapply(rows, function(z) z$status != "VALID", 
            FALSE)), seconds = proc.time()[["elapsed"]] - began, formal_started = FALSE, new_complete_fixed_sequence_CI = "NOT_MEASURED_NO_FULL_SEQUENCE", 
            no_other_endpoint_changed = TRUE), file.path(out, "COMPLETE.json"))
    }), envir = context)
}

