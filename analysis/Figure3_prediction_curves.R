# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

filehash <- function(p) digest(file = p, algo = "sha256")

idhash <- function(x, ordered = FALSE) {
    v <- as.character(x)
    if (!ordered) 
        v <- sort(v, method = "radix")
    digest(enc2utf8(paste0(paste(v, collapse = "\n"), "\n")), algo = "sha256", serialize = FALSE)
}

csv <- function(z, p) {
    assert(!file.exists(p), paste("REFUSE_OVERWRITE", p))
    write.csv(z, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

js <- function(z, p) {
    assert(!file.exists(p), paste("REFUSE_OVERWRITE", p))
    write_json(z, p, pretty = TRUE, auto_unbox = TRUE, digits = NA, na = "null", null = "null")
}

near <- function(x, y, tol = 1e-08) length(x) == length(y) && all(is.finite(c(x, y))) && all(abs(x - y) <= tol)

log_event <- function(event, ...) {
    cat(toJSON(c(list(time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), event = event), list(...)), auto_unbox = TRUE, 
        digits = NA), "\n")
    flush.console()
}

codes <- c(DM = "250", T2D = "250.2", IHD = "411", Angina = "411.3", CAS = "411.4")

fullnames <- c(DM = "Diabetes mellitus", T2D = "Type 2 diabetes", IHD = "Ischemic heart disease", Angina = "Angina pectoris", 
    CAS = "Coronary atherosclerosis")

ci_method <- "Paired_coefficient_delta_method_normal_95CI_fixed_empirical_distribution_pointwise"

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))
        suppressPackageStartupMessages({
            library(jsonlite)
            library(digest)
        })
        args <- commandArgs(TRUE)
        if (length(args) != 2L) stop("Usage: predict_curve.R SAVED_STANDARDIZATION_ROOT NEW_OUTPUT_ROOT")
        assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
        assert(identical(Sys.info()[["sysname"]], "Linux"), "SAVED_MODEL_PREDICTION_RUNS_ON_WINPC_WSL_ONLY")
        SRC <- normalizePath(args[1], mustWork = TRUE)
        OUT <- normalizePath(args[2], mustWork = TRUE)
        assert(!identical(SRC, OUT) && !startsWith(OUT, paste0(SRC, "/")), "OUTPUT_MUST_BE_NEW_AND_ISOLATED")
        assert(!file.exists(file.path(OUT, "RUN_STARTED.json")) && !file.exists(file.path(OUT, "RUN_COMPLETE.json")), 
            "DO_NOT_REPEAT_RUN")
        filehash <- function(p) digest(file = p, algo = "sha256")
        idhash <- function(x, ordered = FALSE) {
            v <- as.character(x)
            if (!ordered) v <- sort(v, method = "radix")
            digest(enc2utf8(paste0(paste(v, collapse = "\n"), "\n")), algo = "sha256", serialize = FALSE)
        }
        csv <- function(z, p) {
            assert(!file.exists(p), paste("REFUSE_OVERWRITE", p))
            write.csv(z, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
        }
        js <- function(z, p) {
            assert(!file.exists(p), paste("REFUSE_OVERWRITE", p))
            write_json(z, p, pretty = TRUE, auto_unbox = TRUE, digits = NA, na = "null", null = "null")
        }
        near <- function(x, y, tol = 1e-08) length(x) == length(y) && all(is.finite(c(x, y))) && all(abs(x - y) <= 
            tol)
        log_event <- function(event, ...) {
            cat(toJSON(c(list(time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), event = event), list(...)), auto_unbox = TRUE, 
                digits = NA), "\n")
            flush.console()
        }
        started <- Sys.time()
        source_files <- c(endpoint = file.path(SRC, "PANEL_C_M3_STANDARDIZATION.csv"), restoration = file.path(SRC, 
            "M3_RESTORE_VERIFICATION.csv"), completion = file.path(SRC, "RUN_COMPLETE.json"), specification = file.path(SRC, 
            "STANDARDIZATION_SPEC.json"))
        assert(all(file.exists(source_files)), "MISSING_SAVED_STANDARDIZATION_PROVENANCE")
        assert(identical(filehash(source_files["endpoint"]), birth_private_setting("PRIVATE_HASH_1")), "OLD_ENDPOINT_TABLE_CHANGED")
        assert(identical(filehash(source_files["restoration"]), birth_private_setting("PRIVATE_HASH_2")), "OLD_RESTORATION_VERIFICATION_CHANGED")
        old <- read.csv(source_files["endpoint"], check.names = FALSE)
        restored <- read.csv(source_files["restoration"], check.names = FALSE)
        done <- read_json(source_files["completion"], simplifyVector = TRUE)
        spec <- read_json(source_files["specification"], simplifyVector = TRUE)
        codes <- c(DM = "250", T2D = "250.2", IHD = "411", Angina = "411.3", CAS = "411.4")
        fullnames <- c(DM = "Diabetes mellitus", T2D = "Type 2 diabetes", IHD = "Ischemic heart disease", Angina = "Angina pectoris", 
            CAS = "Coronary atherosclerosis")
        assert(nrow(old) == 5L && !anyDuplicated(old$disease) && setequal(old$disease, names(codes)), "EXACT_FIVE_OLD_ENDPOINTS_REQUIRED")
        assert(all(old$model == "M3_common") && all(old$weights == "uniform_1_over_N"), "MODEL_OR_STANDARDIZATION_WEIGHT_CHANGED")
        assert(done$status == "PASS" && isTRUE(done$all_sample_hashes_match) && isTRUE(done$all_coefficient_verifications == 
            "PASS"), "ORIGINAL_STANDARDIZATION_NOT_VERIFIED")
        assert(spec$model == "M3_common" && spec$standardization_weights == "UNIFORM_ORIGINAL_GLM_UNWEIGHTED", 
            "SAVED_MODEL_SPECIFICATION_CHANGED")
        assert(all(old$frozen_analysis_source_sha256 == birth_private_setting("PRIVATE_HASH_3")), "ORIGINAL_ANALYSIS_SOURCE_CHANGED")
        model_paths <- setNames(file.path(SRC, "private", paste0(names(codes), "_M3_common_restored_fit.rds")), 
            names(codes))
        cov_paths <- setNames(file.path(SRC, "models", paste0(names(codes), "_M3_covariance.csv")), names(codes))
        assert(all(file.exists(model_paths)) && all(file.exists(cov_paths)), "MISSING_SAVED_MODEL_OR_COVARIANCE_NO_AUTOMATIC_REFIT")
        all_inputs <- c(source_files, setNames(model_paths, paste0("model_", names(codes))), setNames(cov_paths, 
            paste0("covariance_", names(codes))))
        before <- data.frame(role = names(all_inputs), path = unname(all_inputs), bytes = unname(file.info(all_inputs)$size), 
            sha256 = unname(vapply(all_inputs, filehash, character(1))))
        csv(before, file.path(OUT, "INPUT_BINDINGS.csv"))
        js(list(started = format(started, "%Y-%m-%dT%H:%M:%S%z"), model = "M3_common", model_reuse = TRUE, refits_permitted = 0, 
            bootstrap_permitted = FALSE, grid = seq(0, 1, by = 0.05), participant_export_permitted = FALSE), file.path(OUT, 
            "RUN_STARTED.json"))
        grid <- seq(0, 1, by = 0.05)
        ci_method <- "Paired_coefficient_delta_method_normal_95CI_fixed_empirical_distribution_pointwise"
        curves <- bindings <- endpoint_checks <- gradients <- list()
        for (d in names(codes)) {
            log_event("SAVED_MODEL_REUSE_START", disease = d)
            target <- old[old$disease == d, , drop = FALSE]
            assert(identical(filehash(model_paths[d]), as.character(target$model_RDS_sha256)), paste(d, "MODEL_FILE_SHA_MISMATCH"))
            z <- readRDS(model_paths[d])
            fit <- z$fit
            assert(is.list(z) && inherits(fit, "glm") && all(c("sample_ids", "sample_hash", "original_aggregate", 
                "source_sha256") %in% names(z)), "INVALID_SAVED_MODEL_BUNDLE")
            assert(identical(z$source_sha256, as.character(target$frozen_analysis_source_sha256)), paste(d, "SOURCE_BINDING_MISMATCH"))
            assert(identical(z$sample_hash, as.character(target$sample_hash)) && identical(idhash(z$sample_ids), 
                z$sample_hash), paste(d, "SAVED_PARTICIPANT_ID_HASH_MISMATCH"))
            assert(!anyNA(z$sample_ids) && !anyDuplicated(z$sample_ids), "INVALID_SAVED_IDS")
            original <- z$original_aggregate
            assert(nrow(original) == 1L && original$model == "M3_common" && original$disease == d && original$phecode == 
                codes[d] && original$sample_hash == target$sample_hash, "ORIGINAL_AGGREGATE_BINDING_MISMATCH")
            assert(isTRUE(fit$converged) && fit$family$family == "binomial" && fit$family$link == "logit", "SAVED_GLM_FAMILY_OR_CONVERGENCE_CHANGED")
            assert(is.null(fit$call$weights) && all(fit$prior.weights == 1), "ORIGINAL_UNWEIGHTED_GLM_REQUIRED")
            assert(is.null(fit$offset) && is.null(fit$call$offset), "UNEXPECTED_MODEL_OFFSET")
            assert(identical(paste(deparse(formula(fit)), collapse = " "), as.character(original$formula)), "SAVED_FORMULA_MISMATCH")
            mf <- model.frame(fit)
            N <- nrow(mf)
            assert(N == target$N && N == length(z$sample_ids) && N == nobs(fit) && N == original$N, "SAVED_SAMPLE_SIZE_BINDING_MISMATCH")
            assert(sum(mf$y == 1) == target$cases && sum(mf$y == 0) == target$controls, "SAVED_OUTCOME_COUNTS_MISMATCH")
            assert(identical(levels(mf[["factor(S)"]]), as.character(2:8)), "SAVED_SIBSHIP_SUPPORT_NOT_2_TO_8")
            labels <- attr(terms(fit), "term.labels")
            assert(sum(labels == "border_relative") == 1L && sum(grepl("border_relative", labels, fixed = TRUE)) == 
                1L, "B_MUST_REMAIN_ORIGINAL_LINEAR_ADDITIVE_TERM")
            beta <- coef(fit)
            V <- vcov(fit)
            assert(all(is.finite(beta)) && all(is.finite(V)) && identical(names(beta), rownames(V)) && identical(names(beta), 
                colnames(V)), "FULL_COEFFICIENT_COVARIANCE_ALIGNMENT")
            original_v <- read.csv(cov_paths[d], check.names = FALSE)
            assert(identical(original_v$term, names(beta)) && identical(names(original_v)[-1L], names(beta)) && 
                near(as.numeric(as.matrix(original_v[, -1L])), as.numeric(V), 1e-12), "ORIGINAL_MODEL_BASED_COVARIANCE_NOT_REUSED")
            verify <- restored[restored$disease == d, , drop = FALSE]
            assert(nrow(verify) == 1L && verify$sample_hash == target$sample_hash && isTRUE(verify$all_coefficients_match) && 
                near(unname(beta["border_relative"]), verify$beta, 1e-12), "ORIGINAL_RESTORATION_BINDING_FAILED")
            X0 <- model.matrix(fit)
            assert(nrow(X0) == N && identical(colnames(X0), names(beta)) && all(is.finite(X0)) && fit$rank == length(beta), 
                "SAVED_DESIGN_ALIGNMENT")
            ib <- which(colnames(X0) == "border_relative")
            assert(length(ib) == 1L && near(X0[, ib], mf$border_relative, 1e-12) && all(X0[, ib] >= 0 & X0[, ib] <= 
                1), "SAVED_EXPOSURE_ENCODING_MISMATCH")
            X0[, ib] <- 0
            prob0 <- plogis(drop(X0 %*% beta))
            assert(all(is.finite(prob0)) && all(prob0 >= 0 & prob0 <= 1), "ILLEGAL_REFERENCE_PREDICTIONS")
            pbar0 <- mean(prob0)
            g0 <- colMeans(X0 * as.numeric(prob0 * (1 - prob0)))
            out_rows <- list()
            grad_rows <- list()
            all_prob_min <- 1
            all_prob_max <- 0
            for (k in seq_along(grid)) {
                b <- grid[k]
                Xb <- X0
                Xb[, ib] <- b
                prob <- plogis(drop(Xb %*% beta))
                assert(all(is.finite(prob)) && all(prob >= 0 & prob <= 1), "ILLEGAL_GRID_PREDICTIONS")
                all_prob_min <- min(all_prob_min, min(prob))
                all_prob_max <- max(all_prob_max, max(prob))
                pbar <- mean(prob)
                gb <- colMeans(Xb * as.numeric(prob * (1 - prob)))
                gd <- 1000 * (gb - g0)
                delta <- 1000 * (pbar - pbar0)
                variance <- drop(crossprod(gd, V %*% gd))
                assert(is.finite(variance) && variance >= 0, "INVALID_PAIRED_DELTA_VARIANCE")
                se <- sqrt(variance)
                lo <- delta - qnorm(0.975) * se
                hi <- delta + qnorm(0.975) * se
                assert(all(is.finite(c(pbar0, pbar, delta, se, lo, hi))) && lo <= delta && delta <= hi, "INVALID_GRID_DIFFERENCE_OR_CI")
                if (b == 0) assert(identical(delta, 0) && identical(lo, 0) && identical(hi, 0) && all(gd == 0), 
                  "REFERENCE_DIFFERENCE_MUST_BE_EXACTLY_ZERO")
                out_rows[[k]] <- data.frame(disease = unname(fullnames[d]), disease_code = d, phecode = unname(codes[d]), 
                  B = b, N = N, p_B0 = pbar0, p_B = pbar, difference_per1000 = delta, CI_lower = lo, CI_upper = hi, 
                  CI_method = ci_method, SE_difference_per1000 = se, model = "M3_common", sample_hash = z$sample_hash, 
                  ordered_sample_hash = idhash(z$sample_ids, TRUE), weights = "uniform_1_over_N", covariance = "Original stats::vcov.glm model-based full coefficient covariance", 
                  model_RDS_sha256 = as.character(target$model_RDS_sha256), model_RDS_source = unname(model_paths[d]), 
                  frozen_analysis_source_sha256 = z$source_sha256, outcome_scope = "Frozen OLD ever-recorded mapped disease-history state, not certified baseline-only prevalence or incident risk")
                grad_rows[[k]] <- data.frame(disease = d, B = b, term = names(beta), P0_gradient = g0, PB_gradient = gb, 
                  difference_per1000_gradient = gd)
                rm(Xb, prob, gb, gd)
            }
            result <- do.call(rbind, out_rows)
            ep <- result[result$B == 1, , drop = FALSE]
            assert(near(c(ep$p_B0, ep$p_B, ep$difference_per1000, ep$SE_difference_per1000, ep$CI_lower, ep$CI_upper), 
                c(target$P0, target$P1, target$APD_per_1000, target$SE_APD_per_1000, target$CI_low_per_1000, target$CI_high_per_1000), 
                1e-08), paste(d, "ORIGINAL_ENDPOINT_NOT_REPRODUCED_DO_NOT_ADJUST_MODEL"))
            endpoint_checks[[d]] <- data.frame(disease = unname(fullnames[d]), disease_code = d, N = N, sample_hash = z$sample_hash, 
                old_difference = target$APD_per_1000, new_difference = ep$difference_per1000, absolute_difference_error = abs(ep$difference_per1000 - 
                  target$APD_per_1000), old_CI_lower = target$CI_low_per_1000, new_CI_lower = ep$CI_lower, old_CI_upper = target$CI_high_per_1000, 
                new_CI_upper = ep$CI_upper, max_absolute_endpoint_error = max(abs(c(ep$p_B0 - target$P0, ep$p_B - 
                  target$P1, ep$difference_per1000 - target$APD_per_1000, ep$SE_difference_per1000 - target$SE_APD_per_1000, 
                  ep$CI_lower - target$CI_low_per_1000, ep$CI_upper - target$CI_high_per_1000))), status = "PASS")
            bindings[[d]] <- data.frame(disease = unname(fullnames[d]), disease_code = d, model = "M3_common", 
                N = N, cases = sum(mf$y == 1), controls = sum(mf$y == 0), sample_hash = z$sample_hash, ordered_sample_hash = idhash(z$sample_ids, 
                  TRUE), sample_ID_binding = "Exact saved RDS-bound IDs, matched original sorted-ID SHA; not N-only", 
                formula = paste(deparse(formula(fit)), collapse = " "), rank = fit$rank, coefficient_terms = length(beta), 
                weights = "uniform_1_over_N", covariance = "Original stats::vcov.glm model-based full coefficient covariance", 
                B_term = "Original linear additive border_relative, no interaction/nonlinear term", S_levels = paste(levels(mf[["factor(S)"]]), 
                  collapse = ";"), min_individual_predicted_probability = all_prob_min, max_individual_predicted_probability = all_prob_max, 
                P0_SE = sqrt(drop(crossprod(g0, V %*% g0))), model_RDS_sha256 = as.character(target$model_RDS_sha256), 
                model_RDS_source = unname(model_paths[d]), frozen_analysis_source_sha256 = z$source_sha256, original_inference_status = as.character(target$original_status), 
                model_reused = TRUE, model_refitted = FALSE)
            curves[[d]] <- result
            gradients[[d]] <- do.call(rbind, grad_rows)
            log_event("GRID_PREDICTION_COMPLETE", disease = d, N = N, grid_points = nrow(result), endpoint = ep$difference_per1000, 
                max_endpoint_error = endpoint_checks[[d]]$max_absolute_endpoint_error)
            rm(z, fit, mf, X0, V, beta, prob0, g0, out_rows, grad_rows, result)
            gc(FALSE)
        }
        curve_table <- do.call(rbind, curves)
        bind_table <- do.call(rbind, bindings)
        assert(nrow(curve_table) == 105L && !anyDuplicated(paste(curve_table$disease, curve_table$B)), "CURVE_GRID_INCOMPLETE")
        assert(length(unique(bind_table$sample_hash)) == 1L && length(unique(bind_table$ordered_sample_hash)) == 
            1L && length(unique(bind_table$N)) == 1L, "COMMON_ANALYSIS_SAMPLE_NOT_IDENTICAL_ACROSS_DISEASES")
        csv(curve_table, file.path(OUT, "PREDICTION_CURVES.csv"))
        csv(bind_table, file.path(OUT, "MODEL_BINDINGS.csv"))
        csv(do.call(rbind, endpoint_checks), file.path(OUT, "ENDPOINT_CHECK.csv"))
        csv(do.call(rbind, gradients), file.path(OUT, "CURVE_GRADIENTS.csv"))
        after <- data.frame(role = names(all_inputs), path = unname(all_inputs), bytes = unname(file.info(all_inputs)$size), 
            sha256 = unname(vapply(all_inputs, filehash, character(1))))
        assert(identical(before, after), "SAVED_INPUT_CHANGED_DURING_PREDICTION")
        capture.output(sessionInfo(), file = file.path(OUT, "SESSION_INFO.txt"))
        js(list(status = "PASS", models_reused = 5, models_refitted = 0, other_models_refitted = 0, new_bootstrap_jobs = 0, 
            grid_points_per_disease = length(grid), aggregate_rows = nrow(curve_table), common_N = unique(bind_table$N), 
            sample_hash = unique(bind_table$sample_hash), ordered_sample_hash = unique(bind_table$ordered_sample_hash), 
            exact_saved_RDS_sample_model_binding_verified = TRUE, weights = "uniform_1_over_N", CI_method = ci_method, 
            paired_full_covariance = TRUE, pointwise_not_simultaneous_CI = TRUE, B0_difference_and_CI_exactly_zero = TRUE, 
            old_endpoint_reproduction = "PASS", source_files_unchanged = TRUE, participant_data_exported = FALSE, 
            scientific_scope = "Associational prediction of frozen ever-recorded history; not incident risk or causal effect", 
            elapsed_seconds = unname(proc.time()[["elapsed"]]), finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), 
            file.path(OUT, "RUN_COMPLETE.json"))
        log_event("RUN_COMPLETE", models_reused = 5, models_refitted = 0, rows = nrow(curve_table))
    }), envir = context)
}

