# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

filehash <- function(p) digest(file = p, algo = "sha256")

norm <- function(p) normalizePath(p, mustWork = TRUE)

diseases <- c(DM = "250", T2D = "250.2", IHD = "411", Angina = "411.3", CAS = "411.4")

near <- function(x, y, abs_tol = 1e-12, rel_tol = 1e-08) {
    length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) && all(abs(x - y) <= abs_tol + rel_tol * abs(y))
}

safe_rds <- function(object, path) {
    assert(!file.exists(path), paste("REFUSE_EXISTING_MODEL", path))
    saveRDS(object, path, compress = "gzip")
    Sys.chmod(path, mode = "0600")
}

normal_ci <- function(estimate, gradient, V) {
    variance <- drop(crossprod(gradient, V %*% gradient))
    assert(is.finite(variance) && variance >= 0, "NONFINITE_OR_NEGATIVE_DELTA_VARIANCE")
    se <- sqrt(variance)
    z <- qnorm(0.975)
    c(estimate = estimate, SE = se, CI_low = estimate - z * se, CI_high = estimate + z * se)
}

standardize <- function(fit, dd) {
    b <- coef(fit)
    V <- vcov(fit)
    assert(all(is.finite(b)) && all(is.finite(V)) && identical(names(b), colnames(V)) && identical(names(b), rownames(V)), 
        "COEFFICIENT_COVARIANCE_ALIGNMENT")
    terms0 <- delete.response(terms(fit))
    make_x <- function(B) {
        z <- dd
        z$border_relative <- B
        x <- model.matrix(terms0, z, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
        assert(identical(colnames(x), names(b)) && all(is.finite(x)), "COUNTERFACTUAL_MODEL_MATRIX_ALIGNMENT")
        x
    }
    x0 <- make_x(0)
    x1 <- make_x(1)
    p0 <- plogis(drop(x0 %*% b))
    p1 <- plogis(drop(x1 %*% b))
    assert(all(is.finite(p0)) && all(is.finite(p1)), "INVALID_STANDARDIZED_PROBABILITIES")
    g0 <- colMeans(x0 * as.numeric(p0 * (1 - p0)))
    g1 <- colMeans(x1 * as.numeric(p1 * (1 - p1)))
    gdiff <- 1000 * (g1 - g0)
    apd <- normal_ci(1000 * mean(p1 - p0), gdiff, V)
    p0_ci <- normal_ci(mean(p0), g0, V)
    p1_ci <- normal_ci(mean(p1), g1, V)
    paired_cov <- drop(crossprod(g0, V %*% g1))
    assert(near(apd["estimate"], 1000 * (p1_ci["estimate"] - p0_ci["estimate"])), "APD_MEAN_DIFFERENCE_IDENTITY_FAILED")
    assert(near(apd["SE"]^2, 1e+06 * (p0_ci["SE"]^2 + p1_ci["SE"]^2 - 2 * paired_cov)), "PAIRED_DELTA_COVARIANCE_IDENTITY_FAILED")
    dx <- x1 - x0
    other <- setdiff(colnames(dx), "border_relative")
    assert(all(dx[, "border_relative"] == 1) && all(dx[, other, drop = FALSE] == 0), "ENDPOINT_CONTRAST_NOT_EXACTLY_ORIGINAL_B_0_TO_1")
    list(apd = apd, p0 = p0_ci, p1 = p1_ci, paired_covariance = paired_cov, gradients = data.frame(term = names(b), 
        P0_gradient = g0, P1_gradient = g1, APD_per1000_gradient = gdiff), covariance = V)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))
        suppressPackageStartupMessages({
            library(data.table)
            library(jsonlite)
            library(digest)
        })
        argv <- commandArgs(TRUE)
        if (length(argv) != 3L) stop("Usage: standardize_M3.R NEW_OUT FROZEN_RUN_ANALYSIS_R ORIGINAL_RESULT_ROOT", 
            call. = FALSE)
        assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
        assert(identical(Sys.info()[["sysname"]], "Linux"), "WINPC_WSL_ONLY_NO_MAC_PARTICIPANT_DATA")
        filehash <- function(p) digest(file = p, algo = "sha256")
        norm <- function(p) normalizePath(p, mustWork = TRUE)
        SRC <- norm(argv[2])
        ORIGINAL <- norm(argv[3])
        out_parent <- norm(dirname(argv[1]))
        OUT <- file.path(out_parent, basename(argv[1]))
        assert(!file.exists(OUT) && !dir.exists(OUT), "REFUSE_EXISTING_OUTPUT_DIRECTORY")
        assert(!identical(OUT, ORIGINAL) && !startsWith(OUT, paste0(ORIGINAL, "/")), "NEW_OUTPUT_MUST_NOT_BE_INSIDE_ORIGINAL_RESULTS")
        assert(identical(filehash(SRC), birth_private_setting("PRIVATE_HASH_1")), "FROZEN_ANALYSIS_SOURCE_SHA256_MISMATCH")
        expected_files <- c(results = file.path(ORIGINAL, "MODEL_RESULTS_LONG.csv"), manifest = file.path(ORIGINAL, 
            "INPUT_MANIFEST.csv"))
        if (file.exists(file.path(ORIGINAL, "ANALYSIS_SPEC.md"))) expected_files <- c(expected_files, spec = file.path(ORIGINAL, 
            "ANALYSIS_SPEC.md"))
        assert(all(file.exists(expected_files)), "MISSING_ORIGINAL_AGGREGATE_OR_SPEC")
        assert(identical(filehash(expected_files["results"]), birth_private_setting("PRIVATE_HASH_2")), "FROZEN_MODEL_RESULTS_SHA256_MISMATCH")
        assert(identical(filehash(expected_files["manifest"]), birth_private_setting("PRIVATE_HASH_3")), "FROZEN_INPUT_MANIFEST_SHA256_MISMATCH")
        original_scopes <- file.path(ORIGINAL, c("models", "private"))
        assert(all(dir.exists(original_scopes)), "ORIGINAL_MODELS_AND_PRIVATE_SCOPES_REQUIRED")
        existing_fit_files <- unlist(lapply(original_scopes, function(p) list.files(p, pattern = "\\.(rds|rda|rdata|qs|qs2)$", 
            recursive = TRUE, full.names = TRUE, ignore.case = TRUE)), use.names = FALSE)
        assert(!length(existing_fit_files), "ORIGINAL_SAVED_FIT_FOUND_DO_NOT_REFIT")
        original_results <- read.csv(expected_files["results"], check.names = FALSE, colClasses = c(warning = "character", 
            error = "character"))
        original_manifest <- read.csv(expected_files["manifest"], check.names = FALSE)
        diseases <- c(DM = "250", T2D = "250.2", IHD = "411", Angina = "411.3", CAS = "411.4")
        targets <- original_results[original_results$model == "M3_common", , drop = FALSE]
        assert(nrow(targets) == 5L && !anyDuplicated(targets$disease) && setequal(targets$disease, names(diseases)), 
            "EXACT_FIVE_ORIGINAL_M3_ROWS_REQUIRED")
        coef_paths <- setNames(file.path(ORIGINAL, "models", paste0(names(diseases), "_M3_common_coefficients.csv")), 
            paste0("coefficients_", names(diseases)))
        assert(all(file.exists(coef_paths)), "MISSING_ORIGINAL_M3_COEFFICIENT_TABLE")
        provenance_paths <- c(frozen_source = SRC, expected_files, coef_paths)
        provenance_before <- data.frame(role = names(provenance_paths), path = unname(provenance_paths), bytes = unname(file.info(provenance_paths)$size), 
            sha256 = unname(vapply(provenance_paths, filehash, character(1))))
        source_lines <- readLines(SRC, warn = FALSE)
        boundary <- which(source_lines == "results<-flows<-states<-levels_out<-list()")
        assert(length(boundary) == 1L, "ORIGINAL_PREPARATION_BOUNDARY_NOT_UNIQUE")
        prefix <- paste(source_lines[seq_len(boundary - 1L)], collapse = "\n")
        assert(!grepl("glm\\(", prefix), "PREPARATION_MUST_NOT_CONTAIN_A_FIT")
        prefix_hash <- digest(prefix, algo = "sha256", serialize = FALSE)
        assert(dir.create(OUT, recursive = FALSE), "FAILED_TO_CREATE_NEW_OUTPUT_DIRECTORY")
        Sys.chmod(OUT, mode = "0700")
        prep <- new.env(parent = globalenv())
        prep$commandArgs <- function(trailingOnly = FALSE) OUT
        eval(parse(text = prefix, keep.source = FALSE), envir = prep)
        assert(identical(prep$OUT, OUT), "ORIGINAL_PREFIX_OUTPUT_REDIRECTION_FAILED")
        assert(identical(prep$codes, diseases), "FIVE_FROZEN_DISEASES_CHANGED")
        assert(setequal(prep$inputmeta$role, original_manifest$role), "INPUT_ROLE_SET_CHANGED")
        mi <- match(original_manifest$role, prep$inputmeta$role)
        assert(identical(as.character(prep$inputmeta$sha256[mi]), as.character(original_manifest$sha256)) && identical(as.numeric(prep$inputmeta$bytes[mi]), 
            as.numeric(original_manifest$bytes)), "FULL_INPUT_CONTENT_DIFFERS_FROM_FROZEN_MANIFEST")
        prep$csv(provenance_before, file.path(OUT, "FROZEN_SOURCE_MANIFEST.csv"))
        prep$js(list(source_sha256 = filehash(SRC), preparation_prefix_sha256 = prefix_hash, preparation_boundary_line = boundary, 
            preparation_evaluated_verbatim = TRUE, original_fit_search_scopes = original_scopes, original_saved_fits_found = length(existing_fit_files), 
            authorized_models = 5, model = "M3_common", refit_reason = "ORIGINAL_FITS_NOT_SAVED", standardization_weights = "UNIFORM_ORIGINAL_GLM_UNWEIGHTED", 
            outcome = "Frozen ever-recorded history; not certified baseline prevalence or incident risk", uncertainty = "Model-based coefficient delta method conditional on empirical covariate distribution", 
            original_outputs_untouched = TRUE, participant_export_permitted = FALSE), file.path(OUT, "STANDARDIZATION_SPEC.json"))
        Sys.chmod(file.path(OUT, "private"), mode = "0700")
        near <- function(x, y, abs_tol = 1e-12, rel_tol = 1e-08) {
            length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) && all(abs(x - y) <= abs_tol + rel_tol * 
                abs(y))
        }
        safe_rds <- function(object, path) {
            assert(!file.exists(path), paste("REFUSE_EXISTING_MODEL", path))
            saveRDS(object, path, compress = "gzip")
            Sys.chmod(path, mode = "0600")
        }
        normal_ci <- function(estimate, gradient, V) {
            variance <- drop(crossprod(gradient, V %*% gradient))
            assert(is.finite(variance) && variance >= 0, "NONFINITE_OR_NEGATIVE_DELTA_VARIANCE")
            se <- sqrt(variance)
            z <- qnorm(0.975)
            c(estimate = estimate, SE = se, CI_low = estimate - z * se, CI_high = estimate + z * se)
        }
        standardize <- function(fit, dd) {
            b <- coef(fit)
            V <- vcov(fit)
            assert(all(is.finite(b)) && all(is.finite(V)) && identical(names(b), colnames(V)) && identical(names(b), 
                rownames(V)), "COEFFICIENT_COVARIANCE_ALIGNMENT")
            terms0 <- delete.response(terms(fit))
            make_x <- function(B) {
                z <- dd
                z$border_relative <- B
                x <- model.matrix(terms0, z, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
                assert(identical(colnames(x), names(b)) && all(is.finite(x)), "COUNTERFACTUAL_MODEL_MATRIX_ALIGNMENT")
                x
            }
            x0 <- make_x(0)
            x1 <- make_x(1)
            p0 <- plogis(drop(x0 %*% b))
            p1 <- plogis(drop(x1 %*% b))
            assert(all(is.finite(p0)) && all(is.finite(p1)), "INVALID_STANDARDIZED_PROBABILITIES")
            g0 <- colMeans(x0 * as.numeric(p0 * (1 - p0)))
            g1 <- colMeans(x1 * as.numeric(p1 * (1 - p1)))
            gdiff <- 1000 * (g1 - g0)
            apd <- normal_ci(1000 * mean(p1 - p0), gdiff, V)
            p0_ci <- normal_ci(mean(p0), g0, V)
            p1_ci <- normal_ci(mean(p1), g1, V)
            paired_cov <- drop(crossprod(g0, V %*% g1))
            assert(near(apd["estimate"], 1000 * (p1_ci["estimate"] - p0_ci["estimate"])), "APD_MEAN_DIFFERENCE_IDENTITY_FAILED")
            assert(near(apd["SE"]^2, 1e+06 * (p0_ci["SE"]^2 + p1_ci["SE"]^2 - 2 * paired_cov)), "PAIRED_DELTA_COVARIANCE_IDENTITY_FAILED")
            dx <- x1 - x0
            other <- setdiff(colnames(dx), "border_relative")
            assert(all(dx[, "border_relative"] == 1) && all(dx[, other, drop = FALSE] == 0), "ENDPOINT_CONTRAST_NOT_EXACTLY_ORIGINAL_B_0_TO_1")
            list(apd = apd, p0 = p0_ci, p1 = p1_ci, paired_covariance = paired_cov, gradients = data.frame(term = names(b), 
                P0_gradient = g0, P1_gradient = g1, APD_per1000_gradient = gdiff), covariance = V)
        }
        all_results <- verification <- list()
        for (disease in names(diseases)) {
            prep$log_event("RESTORE_M3_STARTED", disease = disease)
            code <- diseases[[disease]]
            dd <- prep$d
            dd$y <- dd[[paste0("Y_", code)]]
            original <- prep$complete(dd, c("y", prep$baseline))
            common <- original & prep$complete(dd, prep$extras)
            dd <- droplevels(dd[common, c("eid", "y", prep$baseline, prep$extras), drop = FALSE])
            old <- targets[targets$disease == disease, , drop = FALSE]
            sample_hash <- prep$idhash(dd$eid)
            assert(nrow(dd) == old$N && sum(dd$y == 1) == old$cases && sum(dd$y == 0) == old$controls && identical(sample_hash, 
                as.character(old$sample_hash)), paste(disease, "ORIGINAL_SAMPLE_IDENTITY_FAILED"))
            f <- prep$forms$M3
            assert(identical(paste(deparse(f), collapse = " "), as.character(old$formula)), paste(disease, "ORIGINAL_FORMULA_CHANGED"))
            fac <- names(dd)[vapply(dd, is.factor, logical(1))]
            refs <- paste(c(vapply(fac, function(k) paste0(k, "=", levels(dd[[k]])[1]), character(1)), paste0("S=", 
                min(dd$S))), collapse = ";")
            assert(identical(refs, as.character(old$reference_levels)), paste(disease, "REFERENCE_LEVELS_CHANGED"))
            x <- model.matrix(f, dd)
            assert(all(is.finite(x)) && ncol(x) == old$columns && qr(x)$rank == old$rank && qr(x)$rank == ncol(x), 
                paste(disease, "DESIGN_OR_RANK_CHANGED"))
            rm(x)
            warnings <- character()
            fit <- withCallingHandlers(glm(f, data = dd, family = binomial(link = "logit"), na.action = na.fail, 
                control = glm.control(maxit = 25), model = TRUE, y = FALSE), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
            assert(isTRUE(fit$converged) && nobs(fit) == nrow(dd) && fit$iter == old$iterations, paste(disease, 
                "FIT_CONVERGENCE_N_OR_ITERATIONS_CHANGED"))
            assert(is.null(fit$call$weights) && all(fit$prior.weights == 1), paste(disease, "ORIGINAL_UNWEIGHTED_FIT_CHANGED"))
            cf <- coef(summary(fit))
            saved <- read.csv(coef_paths[[paste0("coefficients_", disease)]], check.names = FALSE)
            assert(all(saved$disease == disease) & all(saved$model == "M3_common") && !anyDuplicated(saved$term) && 
                identical(rownames(cf), as.character(saved$term)), paste(disease, "COEFFICIENT_TERM_SET_CHANGED"))
            assert(near(cf[, 1], saved$beta) && near(cf[, 2], saved$SE) && near(cf[, 4], saved$P), paste(disease, 
                "ALL_COEFFICIENTS_NOT_REPRODUCED"))
            b <- unname(cf["border_relative", 1])
            se <- unname(cf["border_relative", 2])
            now <- c(beta = b, SE = se, OR = exp(b), CI_low = exp(b - qnorm(0.975) * se), CI_high = exp(b + qnorm(0.975) * 
                se), P = unname(cf["border_relative", 4]))
            assert(near(now, as.numeric(old[1, names(now)])), paste(disease, "TARGET_AGGREGATE_NOT_REPRODUCED"))
            restored_status <- if (length(warnings)) "OK_WARNING_DESCRIPTIVE_PENDING" else "OK_DESCRIPTIVE_PENDING"
            assert(identical(restored_status, as.character(old$status)) && identical(paste(unique(warnings), collapse = "; "), 
                as.character(old$warning)), paste(disease, "ORIGINAL_STATUS_OR_WARNINGS_CHANGED"))
            ans <- standardize(fit, dd)
            row <- data.frame(disease = disease, phecode = code, model = "M3_common", N = nrow(dd), cases = sum(dd$y == 
                1), controls = sum(dd$y == 0), sample_hash = sample_hash, observed_recorded_history_probability = mean(dd$y), 
                standardized_recorded_history_P0 = unname(ans$p0["estimate"]), P0_SE = unname(ans$p0["SE"]), P0_CI_low = unname(ans$p0["CI_low"]), 
                P0_CI_high = unname(ans$p0["CI_high"]), standardized_recorded_history_P1 = unname(ans$p1["estimate"]), 
                P1_SE = unname(ans$p1["SE"]), P1_CI_low = unname(ans$p1["CI_low"]), P1_CI_high = unname(ans$p1["CI_high"]), 
                APD_per1000 = unname(ans$apd["estimate"]), APD_SE = unname(ans$apd["SE"]), APD_CI_low = unname(ans$apd["CI_low"]), 
                APD_CI_high = unname(ans$apd["CI_high"]), P0_P1_covariance = ans$paired_covariance, CI_method = "Paired_coefficient_delta_method_normal_95CI_fixed_empirical_distribution", 
                scale = "1000*(mean_P_at_B1-mean_P_at_B0)", weights = "uniform_1_over_N", original_status = as.character(old$status), 
                original_warning = as.character(old$warning), status = "PASS_ORIGINAL_M3_RESTORED_STANDARDIZATION_COMPLETE")
            clean_env <- new.env(parent = baseenv())
            environment(fit$formula) <- clean_env
            attr(fit$terms, ".Environment") <- clean_env
            attr(fit$model, "terms") <- fit$terms
            fit_path <- file.path(OUT, "private", paste0(disease, "_M3_common_restored_fit.rds"))
            safe_rds(list(fit = fit, sample_ids = dd$eid, original_aggregate = old, sample_hash = sample_hash, 
                source_sha256 = filehash(SRC), prepared_model_frame_in_fit = TRUE, participant_data_remote_only = TRUE), 
                fit_path)
            row$P0 <- row$standardized_recorded_history_P0
            row$P1 <- row$standardized_recorded_history_P1
            row$APD_per_1000 <- row$APD_per1000
            row$SE_APD_per_1000 <- row$APD_SE
            row$CI_low_per_1000 <- row$APD_CI_low
            row$CI_high_per_1000 <- row$APD_CI_high
            row$model_source <- "Original_M3_common_restored_and_all_coefficients_verified"
            row$model_RDS_sha256 <- filehash(fit_path)
            row$frozen_analysis_source_sha256 <- filehash(SRC)
            prep$csv(ans$gradients, file.path(OUT, "models", paste0(disease, "_M3_standardization_gradients.csv")))
            prep$csv(data.frame(term = rownames(ans$covariance), ans$covariance, check.names = FALSE), file.path(OUT, 
                "models", paste0(disease, "_M3_covariance.csv")))
            vr <- data.frame(disease = disease, N = nrow(dd), cases = sum(dd$y == 1), controls = sum(dd$y == 0), 
                sample_hash = sample_hash, coefficient_terms = nrow(cf), all_coefficients_match = TRUE, max_abs_beta_difference = max(abs(cf[, 
                  1] - saved$beta)), max_abs_SE_difference = max(abs(cf[, 2] - saved$SE)), max_abs_P_difference = max(abs(cf[, 
                  4] - saved$P)), target_metrics_match = TRUE, reference_levels_match = TRUE, original_status_match = TRUE, 
                beta = b, SE = se, OR = now["OR"], CI_low = now["CI_low"], CI_high = now["CI_high"], P = now["P"], 
                iterations = fit$iter, rank = fit$rank, converged = fit$converged)
            prep$csv(row, file.path(OUT, "models", paste0(disease, "_M3_standardization.csv")))
            prep$csv(vr, file.path(OUT, "models", paste0(disease, "_M3_RESTORE_VERIFICATION.csv")))
            all_results[[disease]] <- row
            verification[[disease]] <- vr
            prep$log_event("M3_STANDARDIZATION_COMPLETE", disease = disease, N = nrow(dd), cases = sum(dd$y == 
                1), APD_per1000 = row$APD_per1000, CI_low = row$APD_CI_low, CI_high = row$APD_CI_high)
            rm(dd, fit, ans)
            gc(FALSE)
        }
        assert(length(all_results) == 5L, "EXACTLY_FIVE_M3_RESULTS_REQUIRED")
        st1 <- file.info(prep$paths)[, c("size", "mtime")]
        assert(identical(prep$st0, st1), "FROZEN_INPUT_METADATA_CHANGED_DURING_STANDARDIZATION")
        assert(identical(unname(vapply(prep$paths, filehash, character(1))), as.character(prep$inputmeta$sha256)), 
            "FROZEN_INPUT_HASH_CHANGED_DURING_STANDARDIZATION")
        assert(identical(unname(vapply(provenance_paths, filehash, character(1))), provenance_before$sha256), "FROZEN_SOURCE_OR_RESULTS_CHANGED_DURING_STANDARDIZATION")
        prep$csv(as.data.frame(rbindlist(all_results)), file.path(OUT, "PANEL_C_M3_STANDARDIZATION.csv"))
        prep$csv(as.data.frame(rbindlist(verification)), file.path(OUT, "M3_RESTORE_VERIFICATION.csv"))
        session_path <- file.path(OUT, "SESSION_INFO.txt")
        assert(!file.exists(session_path), "REFUSE_EXISTING_SESSION_INFO")
        writeLines(capture.output(sessionInfo()), session_path)
        prep$js(list(status = "PASS", models_refitted = 5, models_standardized = 5, other_models_refitted = 0, 
            bootstrap_replicates = 0, all_coefficient_verifications = "PASS", all_sample_hashes_match = TRUE, full_input_hashes_match_before_and_after = TRUE, 
            source_and_original_result_hashes_unchanged = TRUE, source_sha256 = filehash(SRC), preparation_prefix_sha256 = prefix_hash, 
            R_version = R.version.string, original_outputs_overwritten = FALSE, participant_data_exported = FALSE, 
            model_RDS_location = "PRIVATE_WINPC_ONLY_DO_NOT_COPY_TO_MAC", finished = as.character(Sys.time()), 
            elapsed_seconds = proc.time()[["elapsed"]] - prep$t0), file.path(OUT, "RUN_COMPLETE.json"))
        prep$log_event("RUN_COMPLETE", models = 5, standardized = 5)
    }), envir = context)
}

