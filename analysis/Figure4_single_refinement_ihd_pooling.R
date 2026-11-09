# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ff_m_original_predict_at <- ff_y_predict_at

ff_m_columns <- function(X) {
    a <- "ethnicityAsian or Asian British"
    c <- "ethnicityChinese"
    ff_y_assert(all(c(a, c) %in% colnames(X)), "IHD_POOL_EXPECTED_FROZEN_ETHNICITY_COLUMNS_ABSENT")
    ff_y_assert(all(X[, a] %in% 0:1) && all(X[, c] %in% 0:1) && all((X[, a] + X[, c]) %in% 0:1), "IHD_POOL_NOT_DISJOINT_FROZEN_TREATMENT_DUMMIES")
    z <- X[, colnames(X) != c, drop = FALSE]
    z[, a] <- X[, a] + X[, c]
    colnames(z)[colnames(z) == a] <- "ethnicityAsian_or_Chinese_pooled"
    z
}

ff_m_fit <- function(X, risk, w, tau_report, cluster = NULL, endpoint_id = "411") {
    ff_y_assert(identical(endpoint_id, "411"), "IHD_POOL_FORBIDDEN_OTHER_ENDPOINT")
    runtime <- ff_y_runtime(require = TRUE)
    original <- ff_y_matrix(X)
    designs <- list(disease = original, death = ff_y_matrix(ff_m_columns(original$X)))
    ep <- ff_y_risk(risk, tau_report)
    ff_y_assert(nrow(original$X) == ep$n && length(w) == ep$n && all(is.finite(w) & w > 0), "IHD_POOL_WEIGHT_OR_SAMPLE_INVALID")
    fits <- baselines <- vector("list", 2)
    warnings <- character()
    lp_difference <- numeric(2)
    for (j in 1:2) {
        dd <- designs[[j]]
        work <- dd$work
        work$.ff_time <- ep$time
        work$.ff_disease <- as.integer(ep$status == 1)
        work$.ff_death <- as.integer(ep$status == 2)
        work$.ff_weight <- w
        if (!is.null(cluster)) 
            work$.ff_cluster <- as.character(cluster)
        got <- ff_y_fit_native_one(work, dd$reference$safe_names, j, !is.null(cluster))
        m <- got$value
        warnings <- c(warnings, got$warnings)
        ff_y_assert(!any(grepl("infinite|converg|singular|iteration", got$warnings, ignore.case = TRUE)) && inherits(m, 
            "coxph") && identical(m$method, "efron") && all(is.finite(coef(m))) && all(is.finite(vcov(m))), paste0("IHD_POOL_CAUSE_", 
            j, "_NONCONVERGED"), paste(got$warnings, collapse = " | "))
        mf <- model.frame(m)
        ff_y_assert(m$n == ep$n && identical(as.numeric(m$y[, "time"]), as.numeric(ep$time)) && identical(as.numeric(m$y[, 
            "status"]), as.numeric(ep$status == j)) && identical(as.numeric(model.weights(mf)), as.numeric(w)) && 
            isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), "IHD_POOL_FROZEN_TIME_EVENT_WEIGHTS_OR_CAUSE_DESIGN_CHANGED")
        fits[[j]] <- m
        baselines[[j]] <- ff_y_native_baseline(m, dd$reference$safe_names)
        probe <- unique(as.integer(round(seq(1, ep$n, length.out = min(16, ep$n)))))
        lp <- as.numeric(predict(m, newdata = dd$work[probe, , drop = FALSE], type = "lp", reference = "zero"))
        lp_difference[j] <- max(abs(lp - as.numeric(dd$X[probe, , drop = FALSE] %*% coef(m))))
        ff_y_assert(lp_difference[j] <= 1e-10, "IHD_POOL_NATIVE_LP_CHECK_FAILED")
    }
    names(fits) <- names(baselines) <- c("disease", "death")
    grid <- sort(unique(unlist(lapply(baselines, `[[`, "time"))))
    cumulative <- vapply(baselines, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
    hazard <- rbind(cumulative[1, , drop = FALSE], cumulative[-1, , drop = FALSE] - cumulative[-nrow(cumulative), 
        , drop = FALSE])
    keep <- rowSums(hazard) > 0
    ff_y_assert(all(is.finite(hazard)) && all(hazard >= 0), "IHD_POOL_INVALID_NATIVE_HAZARD")
    coefficients <- do.call(rbind, lapply(1:2, function(j) data.frame(cause = c("disease", "death")[j], term = designs[[j]]$reference$columns, 
        coefficient = unname(coef(fits[[j]])), SE = sqrt(diag(vcov(fits[[j]]))), row.names = NULL)))
    out <- list(fit = fits, engine = list(time = grid[keep], hazard = hazard[keep, , drop = FALSE], baseline_public_validation = baselines), 
        matrix_reference = original$reference, death_reference = designs$death$reference, coefficients = coefficients, 
        endpoint = ep, tau_report = tau_report, ties = "efron", runtime = runtime, ihd_pooling = TRUE, binding = list(weight = w, 
            cluster = cluster), warnings = unique(c(warnings, unlist(lapply(baselines, `[[`, "warnings")))), method = "IHD_DEATH_FIXED_ASIAN_CHINESE_POOL_V2", 
        weight = list(min = min(w), max = max(w), sum = sum(w), ESS = sum(w)^2/sum(w^2), noninteger = sum(abs(w - 
            round(w)) > 1e-10)), native_validation = list(passed = TRUE, lp_max_abs = max(lp_difference), basehaz_survfit_max_abs = max(vapply(baselines, 
            `[[`, 0, "max_abs_basehaz_vs_survfit"))))
    class(out) <- c("ff_y_matrix_fit", "ff_y_fit")
    curve <- lapply(c(0, 0.25, 0.5, 0.75, 1) * tau_report, function(t) ff_y_predict_at(out, original$X[probe, , 
        drop = FALSE], t))
    violation <- 0
    for (k in 2:length(curve)) {
        delta <- curve[[k]] - curve[[k - 1]]
        violation <- max(violation, -delta[, "disease"], -delta[, "death"], delta[, "survival"])
    }
    ff_y_assert(violation <= 1e-12, "IHD_POOL_CIF_TIME_MONOTONICITY_FAILURE")
    out$native_validation$monotonicity_violation <- violation
    out$native_validation$three_state_closure <- max(vapply(curve, function(z) max(abs(rowSums(z) - 1)), 0))
    out
}

ff_y_predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
    64)) {
    if (!isTRUE(object$ihd_pooling)) 
        return(ff_m_original_predict_at(object, X, tau_report, chunk_size, max_prediction_mb))
    ff_y_assert(tau_report >= 0 && tau_report <= object$endpoint$max_time, "IHD_POOL_REPORT_OUTSIDE_FOLLOWUP_SUPPORT")
    dd <- ff_y_matrix(X, object$matrix_reference)
    dm <- ff_y_matrix(ff_m_columns(dd$X), object$death_reference)
    step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
    p <- matrix(NA_real_, nrow(dd$X), 3, dimnames = list(NULL, c("disease", "death", "survival")))
    h <- object$engine$hazard[object$engine$time <= tau_report, , drop = FALSE]
    for (start in seq.int(1L, nrow(dd$X), by = step)) {
        ii <- seq.int(start, min(nrow(dd$X), start + step - 1L))
        lp <- cbind(as.numeric(predict(object$fit$disease, newdata = dd$work[ii, , drop = FALSE], type = "lp", 
            reference = "zero")), as.numeric(predict(object$fit$death, newdata = dm$work[ii, , drop = FALSE], type = "lp", 
            reference = "zero")))
        rr <- exp(lp)
        ff_y_assert(all(is.finite(rr) & rr > 0), "IHD_POOL_RELATIVE_RISK_INVALID")
        p[ii, ] <- ff_aj_backend(h, rr)
        ff_y_probability_check(p[ii, , drop = FALSE])
    }
    ff_y_probability_check(p)
    p
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ff_m_original_predict_at <- ff_y_predict_at
        ff_m_columns <- function(X) {
            a <- "ethnicityAsian or Asian British"
            c <- "ethnicityChinese"
            ff_y_assert(all(c(a, c) %in% colnames(X)), "IHD_POOL_EXPECTED_FROZEN_ETHNICITY_COLUMNS_ABSENT")
            ff_y_assert(all(X[, a] %in% 0:1) && all(X[, c] %in% 0:1) && all((X[, a] + X[, c]) %in% 0:1), "IHD_POOL_NOT_DISJOINT_FROZEN_TREATMENT_DUMMIES")
            z <- X[, colnames(X) != c, drop = FALSE]
            z[, a] <- X[, a] + X[, c]
            colnames(z)[colnames(z) == a] <- "ethnicityAsian_or_Chinese_pooled"
            z
        }
        ff_m_fit <- function(X, risk, w, tau_report, cluster = NULL, endpoint_id = "411") {
            ff_y_assert(identical(endpoint_id, "411"), "IHD_POOL_FORBIDDEN_OTHER_ENDPOINT")
            runtime <- ff_y_runtime(require = TRUE)
            original <- ff_y_matrix(X)
            designs <- list(disease = original, death = ff_y_matrix(ff_m_columns(original$X)))
            ep <- ff_y_risk(risk, tau_report)
            ff_y_assert(nrow(original$X) == ep$n && length(w) == ep$n && all(is.finite(w) & w > 0), "IHD_POOL_WEIGHT_OR_SAMPLE_INVALID")
            fits <- baselines <- vector("list", 2)
            warnings <- character()
            lp_difference <- numeric(2)
            for (j in 1:2) {
                dd <- designs[[j]]
                work <- dd$work
                work$.ff_time <- ep$time
                work$.ff_disease <- as.integer(ep$status == 1)
                work$.ff_death <- as.integer(ep$status == 2)
                work$.ff_weight <- w
                if (!is.null(cluster)) work$.ff_cluster <- as.character(cluster)
                got <- ff_y_fit_native_one(work, dd$reference$safe_names, j, !is.null(cluster))
                m <- got$value
                warnings <- c(warnings, got$warnings)
                ff_y_assert(!any(grepl("infinite|converg|singular|iteration", got$warnings, ignore.case = TRUE)) && 
                  inherits(m, "coxph") && identical(m$method, "efron") && all(is.finite(coef(m))) && all(is.finite(vcov(m))), 
                  paste0("IHD_POOL_CAUSE_", j, "_NONCONVERGED"), paste(got$warnings, collapse = " | "))
                mf <- model.frame(m)
                ff_y_assert(m$n == ep$n && identical(as.numeric(m$y[, "time"]), as.numeric(ep$time)) && identical(as.numeric(m$y[, 
                  "status"]), as.numeric(ep$status == j)) && identical(as.numeric(model.weights(mf)), as.numeric(w)) && 
                  isTRUE(all.equal(unname(m$x), unname(dd$X), tolerance = 1e-14, check.attributes = FALSE)), "IHD_POOL_FROZEN_TIME_EVENT_WEIGHTS_OR_CAUSE_DESIGN_CHANGED")
                fits[[j]] <- m
                baselines[[j]] <- ff_y_native_baseline(m, dd$reference$safe_names)
                probe <- unique(as.integer(round(seq(1, ep$n, length.out = min(16, ep$n)))))
                lp <- as.numeric(predict(m, newdata = dd$work[probe, , drop = FALSE], type = "lp", reference = "zero"))
                lp_difference[j] <- max(abs(lp - as.numeric(dd$X[probe, , drop = FALSE] %*% coef(m))))
                ff_y_assert(lp_difference[j] <= 1e-10, "IHD_POOL_NATIVE_LP_CHECK_FAILED")
            }
            names(fits) <- names(baselines) <- c("disease", "death")
            grid <- sort(unique(unlist(lapply(baselines, `[[`, "time"))))
            cumulative <- vapply(baselines, function(z) c(0, z$cumulative)[findInterval(grid, z$time) + 1L], numeric(length(grid)))
            hazard <- rbind(cumulative[1, , drop = FALSE], cumulative[-1, , drop = FALSE] - cumulative[-nrow(cumulative), 
                , drop = FALSE])
            keep <- rowSums(hazard) > 0
            ff_y_assert(all(is.finite(hazard)) && all(hazard >= 0), "IHD_POOL_INVALID_NATIVE_HAZARD")
            coefficients <- do.call(rbind, lapply(1:2, function(j) data.frame(cause = c("disease", "death")[j], 
                term = designs[[j]]$reference$columns, coefficient = unname(coef(fits[[j]])), SE = sqrt(diag(vcov(fits[[j]]))), 
                row.names = NULL)))
            out <- list(fit = fits, engine = list(time = grid[keep], hazard = hazard[keep, , drop = FALSE], baseline_public_validation = baselines), 
                matrix_reference = original$reference, death_reference = designs$death$reference, coefficients = coefficients, 
                endpoint = ep, tau_report = tau_report, ties = "efron", runtime = runtime, ihd_pooling = TRUE, 
                binding = list(weight = w, cluster = cluster), warnings = unique(c(warnings, unlist(lapply(baselines, 
                  `[[`, "warnings")))), method = "IHD_DEATH_FIXED_ASIAN_CHINESE_POOL_V2", weight = list(min = min(w), 
                  max = max(w), sum = sum(w), ESS = sum(w)^2/sum(w^2), noninteger = sum(abs(w - round(w)) > 1e-10)), 
                native_validation = list(passed = TRUE, lp_max_abs = max(lp_difference), basehaz_survfit_max_abs = max(vapply(baselines, 
                  `[[`, 0, "max_abs_basehaz_vs_survfit"))))
            class(out) <- c("ff_y_matrix_fit", "ff_y_fit")
            curve <- lapply(c(0, 0.25, 0.5, 0.75, 1) * tau_report, function(t) ff_y_predict_at(out, original$X[probe, 
                , drop = FALSE], t))
            violation <- 0
            for (k in 2:length(curve)) {
                delta <- curve[[k]] - curve[[k - 1]]
                violation <- max(violation, -delta[, "disease"], -delta[, "death"], delta[, "survival"])
            }
            ff_y_assert(violation <= 1e-12, "IHD_POOL_CIF_TIME_MONOTONICITY_FAILURE")
            out$native_validation$monotonicity_violation <- violation
            out$native_validation$three_state_closure <- max(vapply(curve, function(z) max(abs(rowSums(z) - 1)), 
                0))
            out
        }
        ff_y_predict_at <- function(object, X, tau_report, chunk_size = getOption("ff.y.chunk_size", 4096L), max_prediction_mb = getOption("ff.y.max_prediction_mb", 
            64)) {
            if (!isTRUE(object$ihd_pooling)) return(ff_m_original_predict_at(object, X, tau_report, chunk_size, 
                max_prediction_mb))
            ff_y_assert(tau_report >= 0 && tau_report <= object$endpoint$max_time, "IHD_POOL_REPORT_OUTSIDE_FOLLOWUP_SUPPORT")
            dd <- ff_y_matrix(X, object$matrix_reference)
            dm <- ff_y_matrix(ff_m_columns(dd$X), object$death_reference)
            step <- as.integer(min(chunk_size, max(1L, floor(max_prediction_mb * 1024^2/(8 * 64)))))
            p <- matrix(NA_real_, nrow(dd$X), 3, dimnames = list(NULL, c("disease", "death", "survival")))
            h <- object$engine$hazard[object$engine$time <= tau_report, , drop = FALSE]
            for (start in seq.int(1L, nrow(dd$X), by = step)) {
                ii <- seq.int(start, min(nrow(dd$X), start + step - 1L))
                lp <- cbind(as.numeric(predict(object$fit$disease, newdata = dd$work[ii, , drop = FALSE], type = "lp", 
                  reference = "zero")), as.numeric(predict(object$fit$death, newdata = dm$work[ii, , drop = FALSE], 
                  type = "lp", reference = "zero")))
                rr <- exp(lp)
                ff_y_assert(all(is.finite(rr) & rr > 0), "IHD_POOL_RELATIVE_RISK_INVALID")
                p[ii, ] <- ff_aj_backend(h, rr)
                ff_y_probability_check(p[ii, , drop = FALSE])
            }
            ff_y_probability_check(p)
            p
        }
    }), envir = context)
}

