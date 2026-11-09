# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ff_w_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

ff_w_json <- function(value, path) {
    ff_w_assert(!file.exists(path), paste("PRESERVE_EXISTING_ARTIFACT", path))
    jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
}

ff_w_flat <- function(row) {
    lapply(row, function(value) {
        if (is.null(value)) 
            NA
        else if (is.list(value) || length(value) != 1L) 
            as.character(jsonlite::toJSON(value, auto_unbox = TRUE, digits = 17, na = "null", null = "null"))
        else value
    })
}

ff_w_csv <- function(rows, path) {
    ff_w_assert(!file.exists(path), paste("PRESERVE_EXISTING_ARTIFACT", path))
    d <- if (length(rows)) 
        v5_bind(lapply(rows, function(r) as.data.frame(ff_w_flat(r), stringsAsFactors = FALSE, check.names = FALSE)))
    else data.frame()
    write.csv(d, path, row.names = FALSE, na = "NA")
}

ff_w_stage <- function(ctx, stage, endpoint = "SHARED", detail = NULL) {
    v5_stage(ctx, stage, endpoint, detail)
}

ff_w_draw <- function(ctx, bid) {
    ff_w_assert(bid %in% 0:100 && bid == as.integer(bid), "FIXED_IDS_ZERO_TO_100_ONLY")
    key <- paste0("draw_", bid)
    if (exists(key, ctx$ff_draw_cache, inherits = FALSE)) 
        return(get(key, ctx$ff_draw_cache))
    directory <- file.path(ctx$old_run, "private/draws", sprintf("bootstrap_%04d", bid))
    seal <- jsonlite::fromJSON(file.path(directory, "COMPLETE.json"), simplifyVector = FALSE)
    path <- file.path(directory, "DRAW.private.rds")
    v5_check(path, seal$private_sha256)
    ff_w_assert(isTRUE(seal$complete) && seal$bootstrap_id == bid && identical(seal$parent_frame_hash, ctx$parent_frame_hash) && 
        identical(seal$parent_cluster_membership_hash, ctx$parent_cluster_membership_hash), "ACTUAL_OLD_CLUSTER_DRAW_BINDING_CHANGED")
    z <- readRDS(path)
    ix <- z$sampled_parent_rows
    multiplicity <- z$cluster_multiplicity
    ff_w_assert(is.numeric(ix) && all(ix %in% seq_len(nrow(ctx$d_full))) && identical(v5_hash(ix), seal$draw_hash) && 
        identical(v5_hash(multiplicity), seal$multiplicity_hash), "ACTUAL_OLD_INDEX_OR_FREQUENCY_CHANGED")
    answer <- list(ix = ix, receipt = c(seal, list(actual_draw_private_sha256 = seal$private_sha256, cached_draw_index_sha256 = seal$private_sha256, 
        actual_index_path = path, draw_retention = "READ_ONLY_OLD_ACTUAL_INDEX_NOT_RESEEDED", cluster_order_hash = ctx$cluster_order_hash)))
    assign(key, answer, ctx$ff_draw_cache)
    old <- ls(ctx$ff_draw_cache)
    if (length(old) > 2L) 
        rm(list = setdiff(old, key)[1L], envir = ctx$ff_draw_cache)
    answer
}

ff_w_base_row <- function(ctx, id, bid, endpoint) {
    s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
    list(indicator_id = id, feature_id = id, label_zh = s$label_zh, type = s$type, domain = s$category, endpoint_id = endpoint, 
        endpoint_label = switch(endpoint, `250` = "DM", `250.2` = "T2D", `411` = "IHD", `411.3` = "AG", `411.4` = "CAS"), 
        bootstrap_id = bid, time_origin = "baseline", horizon_years = ctx$tau_report, tau_report = ctx$tau_report, 
        outcome_fit = "FULL_FOLLOWUP_NATIVE_WEIGHTED_TWO_CAUSE_COX_EFRON_CIF", endpoint_version = ctx$endpoints[[endpoint]]$endpoint_version, 
        source_hash = ctx$source_hash, method_hash = ctx$ff_method_hash, model_hash = ctx$ff_method_hash, input_hash = ctx$input_lock_sha256, 
        status = "NOT_RUN_PHASE", reason = "", TE = NA_real_, DE = NA_real_, IE = NA_real_, PM = NA_real_)
}

ff_w_classify <- function(error, bid) {
    message <- conditionMessage(error)
    technical <- grepl("cannot allocate|cannot open|permission|no space|write error|could not find|no applicable method|object .* not found|unused argument|subscript|SOURCE_SHA256|FILE_MISSING|SOURCE_NOT_AUDITED|BINDING_CHANGED|HASH_MISMATCH|DRAW.*CHANGED|FF_A_SOURCE|FF_A_SAVED_SELECTION_DESIGN", 
        message, ignore.case = TRUE)
    if (technical) 
        "FAILED_TECH"
    else if (bid == 0L) 
        "BLOCKED_POINT"
    else "FAILED_STAT"
}

ff_w_integrate <- function(ctx, a, base, fit, column) {
    dd <- a$data
    w <- a$mi$w
    n <- nrow(dd)
    if (!is.null(a$alpha)) {
        m0 <- dd
        m0[[column]] <- a$M0
        m1 <- dd
        m1[[column]] <- a$M1
        sets <- list(R00 = m0, R10 = m0, R11 = m1)
        for (name in names(sets)) sets[[name]]$B <- if (name == "R00") 
            0
        else 1
        values <- vapply(sets, function(d) {
            X <- ctx$engine$b3_design(d, column, base, a$reference_basis)$X
            pred <- ff_predict_y_matrix(fit, X, ctx$tau_report)
            weighted.mean(pred[, "disease"], w)
        }, 0)
    }
    else {
        recipe <- ctx$discrete$c40_recipe(dd, column, base)
        sums <- c(R00 = 0, R10 = 0, R11 = 0)
        for (level in a$levels) {
            X0 <- ctx$discrete$c40_matrix(recipe, seq_len(n), M = rep(level, n), B = 0)
            X1 <- ctx$discrete$c40_matrix(recipe, seq_len(n), M = rep(level, n), B = 1)
            p0 <- ff_predict_y_matrix(fit, X0, ctx$tau_report)[, "disease"]
            p1 <- ff_predict_y_matrix(fit, X1, ctx$tau_report)[, "disease"]
            sums <- sums + c(R00 = sum(w * a$g0[, level] * p0), R10 = sum(w * a$g0[, level] * p1), R11 = sum(w * 
                a$g1[, level] * p1))
        }
        values <- sums/sum(w)
    }
    closure <- ff_y_effect_closure(values)
    ff_w_assert(isTRUE(closure$passed), "FULLFOLLOWUP_EFFECT_CLOSURE_FAILED")
    list(risks = values, closure = closure, integral_method = a$ledger$simulation_rule, mc_check = list(passed = TRUE, 
        status = "DETERMINISTIC_FROZEN_EMPIRICAL_SHIFT_OR_EXACT_CATEGORY_SUM", Monte_Carlo_draws_added = 0L, mean_substitution = FALSE, 
        explanation = "The inherited V5 integral is deterministic conditional on its real outer sample/A cache; no new inner random draw or MC averaging is introduced."))
}

ff_w_bridge <- function(ctx, d, risk, base, sample_hash, endpoint) {
    key <- v5_hash(c(endpoint, sample_hash, base, ctx$ff_method_hash))
    if (exists(key, ctx$ff_bridge_cache, inherits = FALSE)) 
        return(get(key, ctx$ff_bridge_cache))
    parent <- file.path(ctx$run, "private/total_association_cache")
    if (!dir.exists(parent)) 
        dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    cache <- file.path(parent, key)
    result_path <- file.path(cache, "RESULT.json")
    seal_path <- file.path(cache, "COMPLETE.json")
    cached <- function() {
        seal <- jsonlite::fromJSON(seal_path, simplifyVector = FALSE)
        ff_w_assert(identical(seal$key, key), "TOTAL_ASSOCIATION_CACHE_KEY_CHANGED")
        v5_check(result_path, seal$sha256)
        jsonlite::fromJSON(result_path, simplifyVector = FALSE)
    }
    if (file.exists(seal_path)) 
        return(cached())
    owns <- dir.create(cache, showWarnings = FALSE)
    if (!owns) {
        deadline <- proc.time()[["elapsed"]] + 30
        while (!file.exists(seal_path) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.25)
        if (file.exists(seal_path)) 
            return(cached())
        return(list(status = "WAIT_SHARED_CACHE_OWNER_NOT_REFITTED", cache = cache, sample_hash = sample_hash, 
            formula = paste(base, collapse = " + "), not_mediation_P_or_PM_denominator = TRUE))
    }
    ff_w_json(list(key = key, pid = Sys.getpid(), host = Sys.getenv("FF_HOST"), claimed_utc = v5_utc()), file.path(cache, 
        "CLAIM.json"))
    dd <- d
    dd$.ff_time <- risk$exit
    dd$.ff_event <- risk$event
    form <- reformulate(base, response = "survival::Surv(.ff_time,.ff_event)")
    captured <- ff_y_capture(survival::coxph(form, data = dd, ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, 
        control = survival::coxph.control(iter.max = 50)))
    co <- summary(captured$value)$coefficients
    ff_w_assert("B" %in% rownames(co), "FULLFOLLOWUP_BRIDGE_B_COEFFICIENT_MISSING")
    z <- co["B", ]
    beta <- unname(z["coef"])
    se <- unname(z["se(coef)"])
    answer <- list(N = nrow(dd), events = sum(risk$event), coefficient = beta, SE = se, HR = exp(beta), CI_low = exp(beta - 
        1.96 * se), CI_high = exp(beta + 1.96 * se), P_raw = unname(z["Pr(>|z|)"]), formula = paste(deparse(form), 
        collapse = " "), ties = "efron", full_followup = TRUE, weights = "ORIGINAL_UNWEIGHTED_INCIDENT_COX_ON_IDENTICAL_COMPLETE_SAMPLE", 
        warnings = captured$warnings, sample_hash = sample_hash, not_mediation_P_or_PM_denominator = TRUE)
    ff_w_json(answer, result_path)
    ff_w_json(list(key = key, sha256 = v5_sha(result_path), complete = TRUE), seal_path)
    assign(key, answer, ctx$ff_bridge_cache)
    answer
}

ff_compute <- function(ctx, command) {
    id <- command$indicator_id
    bid <- as.integer(command$bootstrap_id)
    ff_w_assert(id %in% ctx$scope$indicator_id && bid %in% 0:100, "OUTSIDE_LOCKED_121_BY_FIXED100")
    out <- command$output_dir
    ff_w_assert(startsWith(paste0(out, "/"), paste0(ctx$run, "/")) && !grepl("(^|/)\\.\\.(/|$)", out), "TASK_OUTSIDE_NEW_RUN")
    ff_w_assert(!file.exists(out) && dir.create(out, recursive = TRUE), "PRESERVE_EXISTING_ATTEMPT")
    began <- proc.time()[["elapsed"]]
    m <- ctx$monitor
    m$indicator_id <- id
    m$bootstrap_id <- bid
    m$attempt_id <- command$attempt_id
    m$package_started_epoch <- as.numeric(Sys.time())
    m$package_started_monotonic <- began
    m$endpoint_status <- as.list(setNames(rep("NOT_RUN_PHASE", 5), ctx$cfg$outcomes))
    actual <- ff_w_draw(ctx, bid)
    ctx$actual_draw_hash <- actual$receipt$draw_hash
    rows <- list()
    ledgers <- list()
    checks <- list()
    coefficients <- list()
    warnings <- list()
    for (endpoint in ctx$cfg$outcomes) {
        prior <- Filter(function(r) identical(as.character(r$endpoint_id), endpoint), command$existing_rows)
        if (!endpoint %in% unlist(command$active_endpoints)) {
            q <- if (length(prior) == 1L) 
                prior[[1L]]
            else ff_w_base_row(ctx, id, bid, endpoint)
            if (length(prior) != 1L) 
                q$status <- "SKIPPED_SCOPE_FOR_ENGINEERING_CHECK"
            rows[[endpoint]] <- q
            next
        }
        q <- ff_w_base_row(ctx, id, bid, endpoint)
        start <- proc.time()[["elapsed"]]
        tryCatch({
            ep <- ctx$endpoints[[endpoint]]
            ca <- ctx$features[[id]]
            ff_w_assert(identical(ep$status, "AVAILABLE_FROZEN"), paste("OUTCOME_ALIGNMENT_BLOCKED", ep$reason))
            ff_w_assert(identical(ctx$variables[[id]]$status, "AVAILABLE_FROZEN"), paste("FROZEN_FEATURE_BLOCKED", 
                ctx$variables[[id]]$reason))
            target <- ctx$eligibility[[id]] & ep$eligibility
            pr <- actual$ix[target[actual$ix]]
            d <- ctx$d_full[pr, , drop = FALSE]
            risk <- ep$risk[pr, , drop = FALSE]
            base <- unique(c(ctx$cfg$base_covariates, ca$technical_covariates))
            ff_w_stage(ctx, "A_CACHE_REUSE", endpoint)
            a <- ff_a_scenarios(ctx, id, bid, endpoint, pr, d, base)
            q$a_reuse <- as.list(a$ledger[1, , drop = FALSE])
            q$a_reuse_status <- a$ledger$A_origin
            ledgers[[endpoint]] <- a$ledger
            cc <- which(a$mi$complete)
            rr <- risk[cc, , drop = FALSE]
            dd <- a$data
            q$N_target <- nrow(d)
            q$N_complete <- nrow(dd)
            q$target_hash <- v5_hash(pr)
            q$sample_hash <- v5_hash(pr[cc])
            q$standardization_sample_hash <- q$sample_hash
            q$selection_weight_hash <- v5_hash(format(a$mi$w, digits = 17, scientific = TRUE))
            q$formula_hash <- v5_hash(c(base, ca$column, if (is.null(a$alpha)) c(a$levels, "FROZEN_TREATMENT_FACTOR_REFERENCE_FIRST") else "FROZEN_V5_NS_QUANTILE_RULE_5_35_65_95", 
                "FULL_E0_EFRON_TWO_CAUSE"))
            q$actual_basis_reference_hash <- v5_hash(as.character(a$reference_basis))
            q$prepared_transform <- ca$prepared_transform
            q$prepared_unit <- ca$prepared_unit
            q$inverse_transform <- ca$inverse_transform
            q$events_full_followup <- sum(rr$event)
            q$death_full_followup <- sum(rr$death_event)
            q$censor_full_followup <- sum(rr$event == 0 & rr$death_event == 0)
            q$events_by_tau <- sum(rr$event == 1 & rr$exit <= ctx$tau_report)
            q$death_by_tau <- sum(rr$death_event == 1 & rr$exit <= ctx$tau_report)
            q$person_years_full_followup <- sum(rr$exit - rr$entry)
            ff_w_stage(ctx, "FULL_FOLLOWUP_COX", endpoint)
            if (!is.null(a$alpha)) 
                X <- ctx$engine$b3_design(dd, ca$column, base, a$reference_basis)$X
            else {
                recipe <- ctx$discrete$c40_recipe(dd, ca$column, base)
                X <- ctx$discrete$c40_matrix(recipe, seq_len(nrow(dd)), M = dd[[ca$column]])
            }
            fit <- ff_fit_y_matrix(X, rr, a$mi$w, tau_report = ctx$tau_report, ties = "efron", diagnostics = bid == 
                0L)
            representative <- id %in% ctx$ff_config$representative_indicators && endpoint %in% ctx$ff_config$representative_endpoints
            if (representative && bid %in% c(0L, 1L)) {
                ff_w_stage(ctx, "CIF_NUMERICAL_VALIDATION", endpoint)
                validation <- ff_validate_y_matrix(fit, X, rr, a$mi$w, require_rr = TRUE)
                checks[[endpoint]] <- validation
                ff_w_assert(isTRUE(validation$passed) && isTRUE(validation$certified_for_representative_gate), 
                  "REPRESENTATIVE_CIF_VALIDATION_NOT_PASSED")
            }
            else validation <- list(passed = isTRUE(fit$native_validation$passed), native_checks = fit$native_validation, 
                independent_RR_reference = "REQUIRED_BOUND_REPRESENTATIVE_B0_B1_RECEIPTS_NOT_REPEATED_EVERY_FIT")
            ff_w_stage(ctx, "FULL_FOLLOWUP_CIF_INTEGRAL", endpoint)
            integral <- ff_w_integrate(ctx, a, base, fit, ca$column)
            for (name in names(integral$risks)) q[[name]] <- unname(integral$risks[[name]])
            q$TE <- q$R11 - q$R00
            q$DE <- q$R10 - q$R00
            q$IDE <- q$DE
            q$IE <- q$R11 - q$R10
            q$IIE <- q$IE
            q$PM <- if (q$TE != 0) 
                q$IE/q$TE
            else NA_real_
            q$PM_signed <- q$PM
            q$PM_signed_percent <- 100 * q$PM
            q$weak_denominator_diagnostic <- abs(q$TE) <= 1e-08
            q$PM_status <- if (q$TE == 0) 
                "EXACT_ZERO_TE_PM_UNDEFINED"
            else "RAW_SIGNED_SAME_MODEL_RATIO_FIELLER_REQUIRED_NO_CLIPPING"
            q$closure <- q$TE - q$DE - q$IE
            q$draw_hash <- actual$receipt$draw_hash
            q$standardization_hash <- v5_hash(c(q$sample_hash, q$selection_weight_hash, base))
            q$`_draw` <- c(actual$receipt, list(standardization_hash = q$standardization_hash, standardization_sample_hash = q$sample_hash, 
                selection_weight_hash = q$selection_weight_hash))
            coefficients[[endpoint]] <- fit$coefficients
            checks[[endpoint]] <- c(validation, list(effect_closure = integral$closure, integral = integral$mc_check, 
                PH = fit$PH, weights = fit$weight))
            q$warning <- paste(fit$warnings, collapse = " | ")
            warnings[[endpoint]] <- fit$warnings
            q$engineering_checks <- list(validation_pass = isTRUE(validation$passed) && isTRUE(fit$native_validation$passed), 
                mc_error_pass = isTRUE(integral$mc_check$passed), memory_pass = is.finite(ctx$engine$b3_rss()) && 
                  ctx$engine$b3_rss() > 0, method_sample_basis_precision_unchanged = TRUE, wall_seconds = proc.time()[["elapsed"]] - 
                  start, peak_rss_gib = ctx$engine$b3_rss(), evidence = "ACTUAL_NATIVE_HAZARD_LP_VALIDATION_DETERMINISTIC_INTEGRAL_CLOSURE_AND_PROCESS_RSS", 
                memory_is_not_target_concurrency_certification = TRUE)
            if (bid == 0L) 
                q$coxtotalassociation <- tryCatch(ff_w_bridge(ctx, dd, rr, base, q$sample_hash, endpoint), error = function(e) list(status = "TOTAL_ASSOCIATION_DIAGNOSTIC_FAILED_RETAINED", 
                  reason = conditionMessage(e), not_mediation_P_or_PM_denominator = TRUE))
            q$status <- "VALID"
            q$reason <- ""
            rm(fit, X, integral, a)
            gc(FALSE)
        }, error = function(error) {
            q$status <<- ff_w_classify(error, bid)
            q$reason <<- conditionMessage(error)
            if (!is.null(error$validation)) 
                checks[[endpoint]] <<- error$validation
        })
        for (name in c("TE", "DE", "IE")) q[[paste0(name, "_per_1000")]] <- 1000 * q[[name]]
        rows[[endpoint]] <- q
        m$endpoint_status[[endpoint]] <- q$status
        ff_w_stage(ctx, "WRITE_CELL", endpoint)
    }
    ff_w_stage(ctx, "WRITE_PACKAGE")
    ff_w_csv(rows, file.path(out, "FF_RESULTS.csv"))
    ff_w_csv(lapply(ledgers, function(x) as.list(x[1, , drop = FALSE])), file.path(out, "A_REUSE.csv"))
    ff_w_json(actual$receipt, file.path(out, "DRAW.json"))
    ff_w_json(checks, file.path(out, "NUMERICAL_VALIDATION.json"))
    ff_w_json(coefficients, file.path(out, "Y_COEFFICIENTS.json"))
    ff_w_json(warnings, file.path(out, "WARNINGS.json"))
    artifacts <- list.files(out, full.names = TRUE)
    seal <- list(complete = TRUE, indicator_id = id, bootstrap_id = bid, attempt_id = command$attempt_id, source_hash = ctx$source_hash, 
        method_hash = ctx$ff_method_hash, input_hash = ctx$input_lock_sha256, active_endpoints = command$active_endpoints, 
        draw_hash = actual$receipt$draw_hash, artifacts_sha256 = as.list(setNames(vapply(artifacts, v5_sha, ""), 
            basename(artifacts))))
    ff_w_json(seal, file.path(out, "FF_COMPLETE.json"))
    ff_w_stage(ctx, "TASK_FINISHED")
    list(ok = TRUE, indicator_id = id, bootstrap_id = bid, rows = unname(rows), output_dir = out, complete_sha256 = v5_sha(file.path(out, 
        "FF_COMPLETE.json")), total_seconds = proc.time()[["elapsed"]] - began)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ff_w_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        ff_w_json <- function(value, path) {
            ff_w_assert(!file.exists(path), paste("PRESERVE_EXISTING_ARTIFACT", path))
            jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
        }
        ff_w_flat <- function(row) {
            lapply(row, function(value) {
                if (is.null(value)) NA else if (is.list(value) || length(value) != 1L) as.character(jsonlite::toJSON(value, 
                  auto_unbox = TRUE, digits = 17, na = "null", null = "null")) else value
            })
        }
        ff_w_csv <- function(rows, path) {
            ff_w_assert(!file.exists(path), paste("PRESERVE_EXISTING_ARTIFACT", path))
            d <- if (length(rows)) v5_bind(lapply(rows, function(r) as.data.frame(ff_w_flat(r), stringsAsFactors = FALSE, 
                check.names = FALSE))) else data.frame()
            write.csv(d, path, row.names = FALSE, na = "NA")
        }
        ff_w_stage <- function(ctx, stage, endpoint = "SHARED", detail = NULL) {
            v5_stage(ctx, stage, endpoint, detail)
        }
        ff_w_draw <- function(ctx, bid) {
            ff_w_assert(bid %in% 0:100 && bid == as.integer(bid), "FIXED_IDS_ZERO_TO_100_ONLY")
            key <- paste0("draw_", bid)
            if (exists(key, ctx$ff_draw_cache, inherits = FALSE)) return(get(key, ctx$ff_draw_cache))
            directory <- file.path(ctx$old_run, "private/draws", sprintf("bootstrap_%04d", bid))
            seal <- jsonlite::fromJSON(file.path(directory, "COMPLETE.json"), simplifyVector = FALSE)
            path <- file.path(directory, "DRAW.private.rds")
            v5_check(path, seal$private_sha256)
            ff_w_assert(isTRUE(seal$complete) && seal$bootstrap_id == bid && identical(seal$parent_frame_hash, 
                ctx$parent_frame_hash) && identical(seal$parent_cluster_membership_hash, ctx$parent_cluster_membership_hash), 
                "ACTUAL_OLD_CLUSTER_DRAW_BINDING_CHANGED")
            z <- readRDS(path)
            ix <- z$sampled_parent_rows
            multiplicity <- z$cluster_multiplicity
            ff_w_assert(is.numeric(ix) && all(ix %in% seq_len(nrow(ctx$d_full))) && identical(v5_hash(ix), seal$draw_hash) && 
                identical(v5_hash(multiplicity), seal$multiplicity_hash), "ACTUAL_OLD_INDEX_OR_FREQUENCY_CHANGED")
            answer <- list(ix = ix, receipt = c(seal, list(actual_draw_private_sha256 = seal$private_sha256, cached_draw_index_sha256 = seal$private_sha256, 
                actual_index_path = path, draw_retention = "READ_ONLY_OLD_ACTUAL_INDEX_NOT_RESEEDED", cluster_order_hash = ctx$cluster_order_hash)))
            assign(key, answer, ctx$ff_draw_cache)
            old <- ls(ctx$ff_draw_cache)
            if (length(old) > 2L) rm(list = setdiff(old, key)[1L], envir = ctx$ff_draw_cache)
            answer
        }
        ff_w_base_row <- function(ctx, id, bid, endpoint) {
            s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
            list(indicator_id = id, feature_id = id, label_zh = s$label_zh, type = s$type, domain = s$category, 
                endpoint_id = endpoint, endpoint_label = switch(endpoint, `250` = "DM", `250.2` = "T2D", `411` = "IHD", 
                  `411.3` = "AG", `411.4` = "CAS"), bootstrap_id = bid, time_origin = "baseline", horizon_years = ctx$tau_report, 
                tau_report = ctx$tau_report, outcome_fit = "FULL_FOLLOWUP_NATIVE_WEIGHTED_TWO_CAUSE_COX_EFRON_CIF", 
                endpoint_version = ctx$endpoints[[endpoint]]$endpoint_version, source_hash = ctx$source_hash, method_hash = ctx$ff_method_hash, 
                model_hash = ctx$ff_method_hash, input_hash = ctx$input_lock_sha256, status = "NOT_RUN_PHASE", 
                reason = "", TE = NA_real_, DE = NA_real_, IE = NA_real_, PM = NA_real_)
        }
        ff_w_classify <- function(error, bid) {
            message <- conditionMessage(error)
            technical <- grepl("cannot allocate|cannot open|permission|no space|write error|could not find|no applicable method|object .* not found|unused argument|subscript|SOURCE_SHA256|FILE_MISSING|SOURCE_NOT_AUDITED|BINDING_CHANGED|HASH_MISMATCH|DRAW.*CHANGED|FF_A_SOURCE|FF_A_SAVED_SELECTION_DESIGN", 
                message, ignore.case = TRUE)
            if (technical) "FAILED_TECH" else if (bid == 0L) "BLOCKED_POINT" else "FAILED_STAT"
        }
        ff_w_integrate <- function(ctx, a, base, fit, column) {
            dd <- a$data
            w <- a$mi$w
            n <- nrow(dd)
            if (!is.null(a$alpha)) {
                m0 <- dd
                m0[[column]] <- a$M0
                m1 <- dd
                m1[[column]] <- a$M1
                sets <- list(R00 = m0, R10 = m0, R11 = m1)
                for (name in names(sets)) sets[[name]]$B <- if (name == "R00") 0 else 1
                values <- vapply(sets, function(d) {
                  X <- ctx$engine$b3_design(d, column, base, a$reference_basis)$X
                  pred <- ff_predict_y_matrix(fit, X, ctx$tau_report)
                  weighted.mean(pred[, "disease"], w)
                }, 0)
            } else {
                recipe <- ctx$discrete$c40_recipe(dd, column, base)
                sums <- c(R00 = 0, R10 = 0, R11 = 0)
                for (level in a$levels) {
                  X0 <- ctx$discrete$c40_matrix(recipe, seq_len(n), M = rep(level, n), B = 0)
                  X1 <- ctx$discrete$c40_matrix(recipe, seq_len(n), M = rep(level, n), B = 1)
                  p0 <- ff_predict_y_matrix(fit, X0, ctx$tau_report)[, "disease"]
                  p1 <- ff_predict_y_matrix(fit, X1, ctx$tau_report)[, "disease"]
                  sums <- sums + c(R00 = sum(w * a$g0[, level] * p0), R10 = sum(w * a$g0[, level] * p1), R11 = sum(w * 
                    a$g1[, level] * p1))
                }
                values <- sums/sum(w)
            }
            closure <- ff_y_effect_closure(values)
            ff_w_assert(isTRUE(closure$passed), "FULLFOLLOWUP_EFFECT_CLOSURE_FAILED")
            list(risks = values, closure = closure, integral_method = a$ledger$simulation_rule, mc_check = list(passed = TRUE, 
                status = "DETERMINISTIC_FROZEN_EMPIRICAL_SHIFT_OR_EXACT_CATEGORY_SUM", Monte_Carlo_draws_added = 0L, 
                mean_substitution = FALSE, explanation = "The inherited V5 integral is deterministic conditional on its real outer sample/A cache; no new inner random draw or MC averaging is introduced."))
        }
        ff_w_bridge <- function(ctx, d, risk, base, sample_hash, endpoint) {
            key <- v5_hash(c(endpoint, sample_hash, base, ctx$ff_method_hash))
            if (exists(key, ctx$ff_bridge_cache, inherits = FALSE)) return(get(key, ctx$ff_bridge_cache))
            parent <- file.path(ctx$run, "private/total_association_cache")
            if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
            cache <- file.path(parent, key)
            result_path <- file.path(cache, "RESULT.json")
            seal_path <- file.path(cache, "COMPLETE.json")
            cached <- function() {
                seal <- jsonlite::fromJSON(seal_path, simplifyVector = FALSE)
                ff_w_assert(identical(seal$key, key), "TOTAL_ASSOCIATION_CACHE_KEY_CHANGED")
                v5_check(result_path, seal$sha256)
                jsonlite::fromJSON(result_path, simplifyVector = FALSE)
            }
            if (file.exists(seal_path)) return(cached())
            owns <- dir.create(cache, showWarnings = FALSE)
            if (!owns) {
                deadline <- proc.time()[["elapsed"]] + 30
                while (!file.exists(seal_path) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.25)
                if (file.exists(seal_path)) return(cached())
                return(list(status = "WAIT_SHARED_CACHE_OWNER_NOT_REFITTED", cache = cache, sample_hash = sample_hash, 
                  formula = paste(base, collapse = " + "), not_mediation_P_or_PM_denominator = TRUE))
            }
            ff_w_json(list(key = key, pid = Sys.getpid(), host = Sys.getenv("FF_HOST"), claimed_utc = v5_utc()), 
                file.path(cache, "CLAIM.json"))
            dd <- d
            dd$.ff_time <- risk$exit
            dd$.ff_event <- risk$event
            form <- reformulate(base, response = "survival::Surv(.ff_time,.ff_event)")
            captured <- ff_y_capture(survival::coxph(form, data = dd, ties = "efron", singular.ok = FALSE, na.action = stats::na.fail, 
                control = survival::coxph.control(iter.max = 50)))
            co <- summary(captured$value)$coefficients
            ff_w_assert("B" %in% rownames(co), "FULLFOLLOWUP_BRIDGE_B_COEFFICIENT_MISSING")
            z <- co["B", ]
            beta <- unname(z["coef"])
            se <- unname(z["se(coef)"])
            answer <- list(N = nrow(dd), events = sum(risk$event), coefficient = beta, SE = se, HR = exp(beta), 
                CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 1.96 * se), P_raw = unname(z["Pr(>|z|)"]), 
                formula = paste(deparse(form), collapse = " "), ties = "efron", full_followup = TRUE, weights = "ORIGINAL_UNWEIGHTED_INCIDENT_COX_ON_IDENTICAL_COMPLETE_SAMPLE", 
                warnings = captured$warnings, sample_hash = sample_hash, not_mediation_P_or_PM_denominator = TRUE)
            ff_w_json(answer, result_path)
            ff_w_json(list(key = key, sha256 = v5_sha(result_path), complete = TRUE), seal_path)
            assign(key, answer, ctx$ff_bridge_cache)
            answer
        }
        ff_compute <- function(ctx, command) {
            id <- command$indicator_id
            bid <- as.integer(command$bootstrap_id)
            ff_w_assert(id %in% ctx$scope$indicator_id && bid %in% 0:100, "OUTSIDE_LOCKED_121_BY_FIXED100")
            out <- command$output_dir
            ff_w_assert(startsWith(paste0(out, "/"), paste0(ctx$run, "/")) && !grepl("(^|/)\\.\\.(/|$)", out), 
                "TASK_OUTSIDE_NEW_RUN")
            ff_w_assert(!file.exists(out) && dir.create(out, recursive = TRUE), "PRESERVE_EXISTING_ATTEMPT")
            began <- proc.time()[["elapsed"]]
            m <- ctx$monitor
            m$indicator_id <- id
            m$bootstrap_id <- bid
            m$attempt_id <- command$attempt_id
            m$package_started_epoch <- as.numeric(Sys.time())
            m$package_started_monotonic <- began
            m$endpoint_status <- as.list(setNames(rep("NOT_RUN_PHASE", 5), ctx$cfg$outcomes))
            actual <- ff_w_draw(ctx, bid)
            ctx$actual_draw_hash <- actual$receipt$draw_hash
            rows <- list()
            ledgers <- list()
            checks <- list()
            coefficients <- list()
            warnings <- list()
            for (endpoint in ctx$cfg$outcomes) {
                prior <- Filter(function(r) identical(as.character(r$endpoint_id), endpoint), command$existing_rows)
                if (!endpoint %in% unlist(command$active_endpoints)) {
                  q <- if (length(prior) == 1L) prior[[1L]] else ff_w_base_row(ctx, id, bid, endpoint)
                  if (length(prior) != 1L) q$status <- "SKIPPED_SCOPE_FOR_ENGINEERING_CHECK"
                  rows[[endpoint]] <- q
                  next
                }
                q <- ff_w_base_row(ctx, id, bid, endpoint)
                start <- proc.time()[["elapsed"]]
                tryCatch({
                  ep <- ctx$endpoints[[endpoint]]
                  ca <- ctx$features[[id]]
                  ff_w_assert(identical(ep$status, "AVAILABLE_FROZEN"), paste("OUTCOME_ALIGNMENT_BLOCKED", ep$reason))
                  ff_w_assert(identical(ctx$variables[[id]]$status, "AVAILABLE_FROZEN"), paste("FROZEN_FEATURE_BLOCKED", 
                    ctx$variables[[id]]$reason))
                  target <- ctx$eligibility[[id]] & ep$eligibility
                  pr <- actual$ix[target[actual$ix]]
                  d <- ctx$d_full[pr, , drop = FALSE]
                  risk <- ep$risk[pr, , drop = FALSE]
                  base <- unique(c(ctx$cfg$base_covariates, ca$technical_covariates))
                  ff_w_stage(ctx, "A_CACHE_REUSE", endpoint)
                  a <- ff_a_scenarios(ctx, id, bid, endpoint, pr, d, base)
                  q$a_reuse <- as.list(a$ledger[1, , drop = FALSE])
                  q$a_reuse_status <- a$ledger$A_origin
                  ledgers[[endpoint]] <- a$ledger
                  cc <- which(a$mi$complete)
                  rr <- risk[cc, , drop = FALSE]
                  dd <- a$data
                  q$N_target <- nrow(d)
                  q$N_complete <- nrow(dd)
                  q$target_hash <- v5_hash(pr)
                  q$sample_hash <- v5_hash(pr[cc])
                  q$standardization_sample_hash <- q$sample_hash
                  q$selection_weight_hash <- v5_hash(format(a$mi$w, digits = 17, scientific = TRUE))
                  q$formula_hash <- v5_hash(c(base, ca$column, if (is.null(a$alpha)) c(a$levels, "FROZEN_TREATMENT_FACTOR_REFERENCE_FIRST") else "FROZEN_V5_NS_QUANTILE_RULE_5_35_65_95", 
                    "FULL_E0_EFRON_TWO_CAUSE"))
                  q$actual_basis_reference_hash <- v5_hash(as.character(a$reference_basis))
                  q$prepared_transform <- ca$prepared_transform
                  q$prepared_unit <- ca$prepared_unit
                  q$inverse_transform <- ca$inverse_transform
                  q$events_full_followup <- sum(rr$event)
                  q$death_full_followup <- sum(rr$death_event)
                  q$censor_full_followup <- sum(rr$event == 0 & rr$death_event == 0)
                  q$events_by_tau <- sum(rr$event == 1 & rr$exit <= ctx$tau_report)
                  q$death_by_tau <- sum(rr$death_event == 1 & rr$exit <= ctx$tau_report)
                  q$person_years_full_followup <- sum(rr$exit - rr$entry)
                  ff_w_stage(ctx, "FULL_FOLLOWUP_COX", endpoint)
                  if (!is.null(a$alpha)) X <- ctx$engine$b3_design(dd, ca$column, base, a$reference_basis)$X else {
                    recipe <- ctx$discrete$c40_recipe(dd, ca$column, base)
                    X <- ctx$discrete$c40_matrix(recipe, seq_len(nrow(dd)), M = dd[[ca$column]])
                  }
                  fit <- ff_fit_y_matrix(X, rr, a$mi$w, tau_report = ctx$tau_report, ties = "efron", diagnostics = bid == 
                    0L)
                  representative <- id %in% ctx$ff_config$representative_indicators && endpoint %in% ctx$ff_config$representative_endpoints
                  if (representative && bid %in% c(0L, 1L)) {
                    ff_w_stage(ctx, "CIF_NUMERICAL_VALIDATION", endpoint)
                    validation <- ff_validate_y_matrix(fit, X, rr, a$mi$w, require_rr = TRUE)
                    checks[[endpoint]] <- validation
                    ff_w_assert(isTRUE(validation$passed) && isTRUE(validation$certified_for_representative_gate), 
                      "REPRESENTATIVE_CIF_VALIDATION_NOT_PASSED")
                  } else validation <- list(passed = isTRUE(fit$native_validation$passed), native_checks = fit$native_validation, 
                    independent_RR_reference = "REQUIRED_BOUND_REPRESENTATIVE_B0_B1_RECEIPTS_NOT_REPEATED_EVERY_FIT")
                  ff_w_stage(ctx, "FULL_FOLLOWUP_CIF_INTEGRAL", endpoint)
                  integral <- ff_w_integrate(ctx, a, base, fit, ca$column)
                  for (name in names(integral$risks)) q[[name]] <- unname(integral$risks[[name]])
                  q$TE <- q$R11 - q$R00
                  q$DE <- q$R10 - q$R00
                  q$IDE <- q$DE
                  q$IE <- q$R11 - q$R10
                  q$IIE <- q$IE
                  q$PM <- if (q$TE != 0) q$IE/q$TE else NA_real_
                  q$PM_signed <- q$PM
                  q$PM_signed_percent <- 100 * q$PM
                  q$weak_denominator_diagnostic <- abs(q$TE) <= 1e-08
                  q$PM_status <- if (q$TE == 0) "EXACT_ZERO_TE_PM_UNDEFINED" else "RAW_SIGNED_SAME_MODEL_RATIO_FIELLER_REQUIRED_NO_CLIPPING"
                  q$closure <- q$TE - q$DE - q$IE
                  q$draw_hash <- actual$receipt$draw_hash
                  q$standardization_hash <- v5_hash(c(q$sample_hash, q$selection_weight_hash, base))
                  q$`_draw` <- c(actual$receipt, list(standardization_hash = q$standardization_hash, standardization_sample_hash = q$sample_hash, 
                    selection_weight_hash = q$selection_weight_hash))
                  coefficients[[endpoint]] <- fit$coefficients
                  checks[[endpoint]] <- c(validation, list(effect_closure = integral$closure, integral = integral$mc_check, 
                    PH = fit$PH, weights = fit$weight))
                  q$warning <- paste(fit$warnings, collapse = " | ")
                  warnings[[endpoint]] <- fit$warnings
                  q$engineering_checks <- list(validation_pass = isTRUE(validation$passed) && isTRUE(fit$native_validation$passed), 
                    mc_error_pass = isTRUE(integral$mc_check$passed), memory_pass = is.finite(ctx$engine$b3_rss()) && 
                      ctx$engine$b3_rss() > 0, method_sample_basis_precision_unchanged = TRUE, wall_seconds = proc.time()[["elapsed"]] - 
                      start, peak_rss_gib = ctx$engine$b3_rss(), evidence = "ACTUAL_NATIVE_HAZARD_LP_VALIDATION_DETERMINISTIC_INTEGRAL_CLOSURE_AND_PROCESS_RSS", 
                    memory_is_not_target_concurrency_certification = TRUE)
                  if (bid == 0L) q$coxtotalassociation <- tryCatch(ff_w_bridge(ctx, dd, rr, base, q$sample_hash, 
                    endpoint), error = function(e) list(status = "TOTAL_ASSOCIATION_DIAGNOSTIC_FAILED_RETAINED", 
                    reason = conditionMessage(e), not_mediation_P_or_PM_denominator = TRUE))
                  q$status <- "VALID"
                  q$reason <- ""
                  rm(fit, X, integral, a)
                  gc(FALSE)
                }, error = function(error) {
                  q$status <<- ff_w_classify(error, bid)
                  q$reason <<- conditionMessage(error)
                  if (!is.null(error$validation)) checks[[endpoint]] <<- error$validation
                })
                for (name in c("TE", "DE", "IE")) q[[paste0(name, "_per_1000")]] <- 1000 * q[[name]]
                rows[[endpoint]] <- q
                m$endpoint_status[[endpoint]] <- q$status
                ff_w_stage(ctx, "WRITE_CELL", endpoint)
            }
            ff_w_stage(ctx, "WRITE_PACKAGE")
            ff_w_csv(rows, file.path(out, "FF_RESULTS.csv"))
            ff_w_csv(lapply(ledgers, function(x) as.list(x[1, , drop = FALSE])), file.path(out, "A_REUSE.csv"))
            ff_w_json(actual$receipt, file.path(out, "DRAW.json"))
            ff_w_json(checks, file.path(out, "NUMERICAL_VALIDATION.json"))
            ff_w_json(coefficients, file.path(out, "Y_COEFFICIENTS.json"))
            ff_w_json(warnings, file.path(out, "WARNINGS.json"))
            artifacts <- list.files(out, full.names = TRUE)
            seal <- list(complete = TRUE, indicator_id = id, bootstrap_id = bid, attempt_id = command$attempt_id, 
                source_hash = ctx$source_hash, method_hash = ctx$ff_method_hash, input_hash = ctx$input_lock_sha256, 
                active_endpoints = command$active_endpoints, draw_hash = actual$receipt$draw_hash, artifacts_sha256 = as.list(setNames(vapply(artifacts, 
                  v5_sha, ""), basename(artifacts))))
            ff_w_json(seal, file.path(out, "FF_COMPLETE.json"))
            ff_w_stage(ctx, "TASK_FINISHED")
            list(ok = TRUE, indicator_id = id, bootstrap_id = bid, rows = unname(rows), output_dir = out, complete_sha256 = v5_sha(file.path(out, 
                "FF_COMPLETE.json")), total_seconds = proc.time()[["elapsed"]] - began)
        }
    }), envir = context)
}

