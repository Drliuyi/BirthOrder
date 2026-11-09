# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

bl2_capture <- function(expr) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = identity)
    list(value = value, warnings = unique(warnings))
}

bl2_assert <- function(ok, message) {
    if (!isTRUE(ok)) 
        stop(message, call. = FALSE)
}

bl2_status_factor <- function(x) {
    factor(as.character(x), levels = c("0", "1", "2"))
}

bl2_formula <- function() {
    y ~ border_relative + age + sex + ethnicity + center + assessment_year + factor(S)
}

bl2_contrasts <- function(formula, data) {
    mf <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
    is_factor <- vapply(mf, is.factor, logical(1))
    is_factor[1L] <- FALSE
    lapply(mf[is_factor], function(x) stats::contr.treatment(levels(x), base = 1L))
}

bl2_sample <- function(data, outcome) {
    bl2_assert(is.data.frame(data), "MODEL_READY_DATA_FRAME_REQUIRED")
    required <- c("eid", outcome, "border_relative", "age", "sex", "ethnicity", "center", "assessment_year", "S")
    missing <- setdiff(required, names(data))
    bl2_assert(!length(missing), paste("MISSING_MODEL_FIELDS", paste(missing, collapse = ",")))
    bl2_assert(!anyNA(data$eid) && !anyDuplicated(as.character(data$eid)), "DUPLICATE_OR_MISSING_ID")
    bl2_assert(is.numeric(data$border_relative), "NUMERIC_BORDER_RELATIVE_REQUIRED")
    bl2_assert(is.numeric(data$age), "LINEAR_NUMERIC_AGE_REQUIRED")
    vars <- setdiff(required, "eid")
    keep <- stats::complete.cases(data[, vars, drop = FALSE])
    for (v in vars) if (is.numeric(data[[v]])) 
        keep <- keep & is.finite(data[[v]])
    rows <- which(keep)
    d <- droplevels(as.data.frame(data[rows, , drop = FALSE]))
    d$y <- data[[outcome]][rows]
    if (is.factor(d$y)) 
        d$y <- factor(as.character(d$y), levels = levels(data[[outcome]]))
    rownames(d) <- as.character(rows)
    bl2_assert(nrow(d) > 0L, "NO_COMPLETE_CASES")
    bl2_assert(all(d$border_relative >= 0 & d$border_relative <= 1), "BORDER_RELATIVE_OUTSIDE_0_1")
    s <- suppressWarnings(as.numeric(as.character(d$S)))
    bl2_assert(all(is.finite(s) & s %in% 2:8), "S2_TO_8_REQUIRED")
    bl2_assert(length(unique(d$border_relative)) > 1L, "CONSTANT_EXPOSURE")
    formula <- bl2_formula()
    labels <- attr(stats::terms(formula), "term.labels")
    removed <- character()
    for (term in setdiff(labels, "border_relative")) {
        v <- all.vars(stats::as.formula(paste("~", term)))
        if (length(v) == 1L && length(unique(d[[v]])) <= 1L) 
            removed <- c(removed, term)
    }
    if (length(removed)) 
        formula <- stats::reformulate(setdiff(labels, removed), response = "y")
    null_formula <- stats::reformulate(setdiff(attr(stats::terms(formula), "term.labels"), "border_relative"), 
        response = "y")
    contrasts <- bl2_contrasts(formula, d)
    x <- stats::model.matrix(formula, d, contrasts.arg = contrasts)
    x0 <- stats::model.matrix(null_formula, d, contrasts.arg = contrasts)
    list(data = d, formula = formula, null_formula = null_formula, rows = rows, eids = as.character(d$eid), removed_terms = removed, 
        contrasts = contrasts, x = x, x0 = x0, rank = qr(x)$rank, columns = ncol(x), rank_null = qr(x0)$rank, columns_null = ncol(x0))
}

bl2_global_template <- function(data, outcome, model) {
    data.frame(outcome = outcome, model = model, exposure = "border_relative", N_input = nrow(data), N = NA_integer_, 
        N_excluded = NA_integer_, status = "NOT_ESTIMABLE", error = "", warning = "", formula = paste(deparse(bl2_formula()), 
            collapse = " "), null_formula = "", removed_terms = "", rank = NA_integer_, columns = NA_integer_, 
        rank_null = NA_integer_, columns_null = NA_integer_, beta = NA_real_, SE = NA_real_, CI_low = NA_real_, 
        CI_high = NA_real_, beta_SD = NA_real_, SE_SD = NA_real_, CI_low_SD = NA_real_, CI_high_SD = NA_real_, 
        outcome_SD = NA_real_, df_residual = NA_real_, p = NA_real_, p_raw = NA_real_, raw_beta = NA_real_, raw_SE = NA_real_, 
        raw_CI_low = NA_real_, raw_CI_high = NA_real_, raw_beta_SD = NA_real_, raw_SE_SD = NA_real_, raw_CI_low_SD = NA_real_, 
        raw_CI_high_SD = NA_real_, LR = NA_real_, LR_raw = NA_real_, df = if (model == "lm") 
            1L
        else 2L, df_raw = NA_real_, logLik_full = NA_real_, logLik_null = NA_real_, convergence_full = NA_integer_, 
        convergence_null = NA_integer_, numeric_reliable = FALSE, inference = "ASSOCIATION_NONCAUSAL", stringsAsFactors = FALSE)
}

bl2_add_sample_metadata <- function(global, sample) {
    global$N <- nrow(sample$data)
    global$N_excluded <- global$N_input - global$N
    global$formula <- paste(deparse(sample$formula), collapse = " ")
    global$null_formula <- paste(deparse(sample$null_formula), collapse = " ")
    global$removed_terms <- paste(sample$removed_terms, collapse = ";")
    for (nm in c("rank", "columns", "rank_null", "columns_null")) global[[nm]] <- sample[[nm]]
    global
}

bl2_sample_diagnostics <- function(sample) {
    list(sample_rows = sample$rows, sample_eids = sample$eids, removed_constant_terms = sample$removed_terms, full_design_columns = colnames(sample$x), 
        null_design_columns = colnames(sample$x0), rank_full = sample$rank, columns_full = sample$columns, rank_null = sample$rank_null, 
        columns_null = sample$columns_null)
}

bl2_assert_design <- function(sample) {
    bl2_assert(all(is.finite(sample$x)) && all(is.finite(sample$x0)), "NONFINITE_DESIGN")
    bl2_assert(sample$rank == sample$columns && sample$rank_null == sample$columns_null, "RANK_DEFICIENT_DESIGN")
    bl2_assert("border_relative" %in% colnames(sample$x), "TARGET_COLUMN_MISSING")
}

bl2_fit_bmi <- function(data, outcome = "bmi") {
    global <- bl2_global_template(data, outcome, "lm")
    diagnostics <- list(raw_coefficients = data.frame(), warnings = character())
    tryCatch({
        bl2_assert(outcome %in% names(data) && is.numeric(data[[outcome]]), "NUMERIC_BMI_OUTCOME_REQUIRED")
        sample <- bl2_sample(data, outcome)
        global <- bl2_add_sample_metadata(global, sample)
        diagnostics <- c(diagnostics, bl2_sample_diagnostics(sample))
        bl2_assert_design(sample)
        bl2_assert(nrow(sample$data) > sample$columns, "NO_RESIDUAL_DEGREES_OF_FREEDOM")
        global$outcome_SD <- stats::sd(sample$data$y)
        bl2_assert(is.finite(global$outcome_SD) && global$outcome_SD > 0, "CONSTANT_BMI_OUTCOME")
        captured <- bl2_capture(stats::lm(sample$formula, data = sample$data, contrasts = sample$contrasts, na.action = stats::na.fail))
        diagnostics$warnings <- captured$warnings
        global$warning <- paste(captured$warnings, collapse = "; ")
        if (inherits(captured$value, "error")) 
            stop(conditionMessage(captured$value), call. = FALSE)
        fit <- captured$value
        captured_summary <- bl2_capture(summary(fit))
        diagnostics$warnings <- unique(c(diagnostics$warnings, captured_summary$warnings))
        global$warning <- paste(diagnostics$warnings, collapse = "; ")
        if (inherits(captured_summary$value, "error")) 
            stop(conditionMessage(captured_summary$value), call. = FALSE)
        cf <- captured_summary$value$coefficients
        diagnostics$raw_coefficients <- data.frame(term = rownames(cf), raw_beta = cf[, 1L], raw_SE = cf[, 2L], 
            raw_t = cf[, 3L], raw_p = cf[, 4L], row.names = NULL)
        global$df_residual <- fit$df.residual
        global$raw_beta <- cf["border_relative", "Estimate"]
        global$raw_SE <- cf["border_relative", "Std. Error"]
        global$p_raw <- cf["border_relative", "Pr(>|t|)"]
        crit <- stats::qt(0.975, df = fit$df.residual)
        global$raw_CI_low <- global$raw_beta - crit * global$raw_SE
        global$raw_CI_high <- global$raw_beta + crit * global$raw_SE
        global$raw_beta_SD <- global$raw_beta/global$outcome_SD
        global$raw_SE_SD <- global$raw_SE/global$outcome_SD
        global$raw_CI_low_SD <- global$raw_CI_low/global$outcome_SD
        global$raw_CI_high_SD <- global$raw_CI_high/global$outcome_SD
        diagnostics$target_numeric_reliable <- all(is.finite(c(global$raw_beta, global$raw_SE, global$p_raw, global$raw_CI_low, 
            global$raw_CI_high))) && global$raw_SE > 0
        bl2_assert(all(is.finite(cf)) && diagnostics$target_numeric_reliable, "NONFINITE_OR_INVALID_LM_INFERENCE")
        global$beta <- global$raw_beta
        global$SE <- global$raw_SE
        global$CI_low <- global$raw_CI_low
        global$CI_high <- global$raw_CI_high
        global$beta_SD <- global$raw_beta_SD
        global$SE_SD <- global$raw_SE_SD
        global$CI_low_SD <- global$raw_CI_low_SD
        global$CI_high_SD <- global$raw_CI_high_SD
        global$p <- global$p_raw
        global$numeric_reliable <- TRUE
        global$status <- if (length(diagnostics$warnings)) 
            "FIT_WARNING_REVIEW"
        else "OK"
    }, error = function(e) {
        global$error <<- conditionMessage(e)
        global$status <<- "NOT_ESTIMABLE"
    })
    list(global = global, contrasts = data.frame(), diagnostics = diagnostics)
}

bl2_numerical_scale <- function(data, formula) {
    d <- data
    transform <- data.frame(variable = character(), center = numeric(), scale = numeric(), original_unit = character(), 
        stringsAsFactors = FALSE)
    used <- all.vars(formula)
    for (v in intersect(c("age", "assessment_year"), used)) {
        if (!is.numeric(d[[v]])) 
            next
        center <- min(d[[v]])
        scale <- diff(range(d[[v]]))
        bl2_assert(is.finite(scale) && scale > 0, "INVALID_NUMERICAL_SCALE")
        d[[v]] <- (d[[v]] - center)/scale
        transform <- rbind(transform, data.frame(variable = v, center = center, scale = scale, original_unit = "year"))
    }
    list(data = d, transform = transform)
}

bl2_backmap <- function(coef_names, transform, classes = 2L) {
    tmat <- diag(length(coef_names))
    dimnames(tmat) <- list(coef_names, coef_names)
    for (i in seq_len(nrow(transform))) {
        v <- transform$variable[i]
        bl2_assert(v %in% coef_names && "(Intercept)" %in% coef_names, "AFFINE_BACKMAP_COLUMN_MISSING")
        tmat[v, v] <- 1/transform$scale[i]
        tmat["(Intercept)", v] <- -transform$center[i]/transform$scale[i]
    }
    list(one_class = tmat, all_classes = kronecker(diag(classes), tmat))
}

bl2_hessian_diagnostics <- function(hessian) {
    result <- list(finite = FALSE, dimension = NROW(hessian), rank = NA_integer_, minimum_eigenvalue = NA_real_, 
        maximum_eigenvalue = NA_real_, condition_number = Inf, tolerance = sqrt(.Machine$double.eps), reliable = FALSE, 
        error = "")
    if (!is.matrix(hessian) || !nrow(hessian) || nrow(hessian) != ncol(hessian) || any(!is.finite(hessian))) {
        result$error <- "NONFINITE_OR_INVALID_HESSIAN"
        return(result)
    }
    result$finite <- TRUE
    captured <- bl2_capture(eigen((hessian + t(hessian))/2, symmetric = TRUE, only.values = TRUE)$values)
    if (inherits(captured$value, "error")) {
        result$error <- conditionMessage(captured$value)
        return(result)
    }
    ev <- captured$value
    result$minimum_eigenvalue <- min(ev)
    result$maximum_eigenvalue <- max(ev)
    result$rank <- sum(ev > max(ev) * result$tolerance)
    if (min(ev) > 0) 
        result$condition_number <- max(ev)/min(ev)
    result$reliable <- max(ev) > 0 && min(ev) > 0 && result$rank == nrow(hessian)
    if (!result$reliable) 
        result$error <- "SINGULAR_OR_SEVERELY_ILL_CONDITIONED_HESSIAN"
    result
}

bl2_multinom_diagnostics <- function(fit, formula, data, contrasts, transform) {
    cf <- stats::coef(fit)
    bl2_assert(is.matrix(cf) && identical(rownames(cf), c("1", "2")), "THREE_CLASS_COEFFICIENT_LAYOUT_MISMATCH")
    backmap <- bl2_backmap(colnames(cf), transform, nrow(cf))
    raw_cf <- cf %*% t(backmap$one_class)
    vc <- bl2_capture(stats::vcov(fit))
    raw_vcov <- if (!inherits(vc$value, "error")) 
        backmap$all_classes %*% vc$value %*% t(backmap$all_classes)
    else matrix(NA_real_, nrow(cf) * ncol(cf), nrow(cf) * ncol(cf))
    variances <- diag(raw_vcov)
    raw_se <- matrix(sqrt(ifelse(variances >= 0, variances, NA_real_)), nrow = nrow(cf), byrow = TRUE, dimnames = dimnames(cf))
    raw_table <- data.frame(outcome_level = rep(rownames(raw_cf), each = ncol(raw_cf)), term = rep(colnames(raw_cf), 
        times = nrow(raw_cf)), raw_beta = as.vector(t(raw_cf)), raw_SE = as.vector(t(raw_se)), stringsAsFactors = FALSE)
    raw_table$raw_z <- raw_table$raw_beta/raw_table$raw_SE
    raw_table$raw_p <- 2 * stats::pnorm(abs(raw_table$raw_z), lower.tail = FALSE)
    target <- raw_table[raw_table$term == "border_relative", , drop = FALSE]
    hessian <- bl2_hessian_diagnostics(fit$Hessian)
    prob <- fit$fitted.values
    x <- stats::model.matrix(formula, data, contrasts.arg = contrasts)
    response <- sapply(c("1", "2"), function(k) as.numeric(data$y == k))
    score <- crossprod(x, response - prob[, c("1", "2"), drop = FALSE])
    target_reliable <- nrow(target) == 2L && all(is.finite(target$raw_beta)) && all(is.finite(target$raw_SE) & 
        target$raw_SE > 0) && all(is.finite(target$raw_p))
    issues <- character()
    if (!identical(fit$convergence, 0L)) 
        issues <- c(issues, "MULTINOM_MAXIT_OR_NONCONVERGENCE")
    if (!all(is.finite(cf)) || !all(is.finite(raw_vcov)) || !all(is.finite(prob)) || !is.finite(as.numeric(stats::logLik(fit)))) 
        issues <- c(issues, "NONFINITE_MULTINOM_FIT")
    if (!hessian$reliable) 
        issues <- c(issues, hessian$error)
    if (!target_reliable) 
        issues <- c(issues, "INVALID_TARGET_ESTIMATE")
    if (inherits(vc$value, "error")) 
        issues <- c(issues, paste("VCOV_ERROR", conditionMessage(vc$value)))
    list(convergence = fit$convergence, hessian = hessian, raw_coefficients = raw_table, raw_vcov = raw_vcov, coefficient_backmap = backmap$one_class, 
        target_numeric_reliable = target_reliable, max_abs_score_internal = max(abs(score)), target_score_internal = if ("border_relative" %in% 
            rownames(score)) score["border_relative", ] else rep(NA_real_, 2L), minimum_fitted_probability = min(prob), 
        maximum_fitted_probability = max(prob), extreme_probability_flag = any(prob < 1e-12 | prob > 1 - 1e-12), 
        vcov_warnings = vc$warnings, reliable = !length(issues), issues = unique(issues), nobs = nrow(prob), sample_rows = rownames(prob), 
        logLik = as.numeric(stats::logLik(fit)), edf = fit$edf)
}

bl2_contrast_template <- function(outcome) {
    data.frame(outcome = outcome, exposure = "border_relative", reference_level = "0", outcome_level = c("1", "2"), 
        contrast = c("Previous_vs_Never", "Current_vs_Never"), N = NA_integer_, beta = NA_real_, SE = NA_real_, 
        relative_odds_ratio = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, raw_beta = NA_real_, 
        raw_SE = NA_real_, raw_relative_odds_ratio = NA_real_, raw_CI_low = NA_real_, raw_CI_high = NA_real_, p_raw = NA_real_, 
        status = "NOT_ESTIMABLE", error = "", warning = "", stringsAsFactors = FALSE)
}

bl2_fit_multinomial <- function(data, outcome, maxit = 1000L) {
    global <- bl2_global_template(data, outcome, "nnet::multinom")
    contrasts <- bl2_contrast_template(outcome)
    diagnostics <- list(raw_coefficients = data.frame(), warnings = character(), class_levels = c("0", "1", "2"), 
        reference_level = "0", optimizer = list(maxit = maxit, decay = 0, reltol = 1e-10, Hess = TRUE, summ = 0L), 
        numeric_reparameterization = "affine nuisance-only age/assessment_year; backmapped")
    tryCatch({
        bl2_assert(requireNamespace("nnet", quietly = TRUE), "NNET_PACKAGE_REQUIRED")
        bl2_assert(outcome %in% names(data) && is.factor(data[[outcome]]) && identical(levels(data[[outcome]]), 
            c("0", "1", "2")), "VERIFIED_THREE_LEVEL_FACTOR_0_REFERENCE_REQUIRED")
        bl2_assert(length(maxit) == 1L && is.finite(maxit) && maxit >= 1 && maxit == as.integer(maxit), "INVALID_MAXIT")
        sample <- bl2_sample(data, outcome)
        global <- bl2_add_sample_metadata(global, sample)
        contrasts$N <- global$N
        diagnostics <- c(diagnostics, bl2_sample_diagnostics(sample))
        diagnostics$class_counts <- table(sample$data$y)
        bl2_assert(all(diagnostics$class_counts > 0L), "MISSING_OUTCOME_CLASS_NO_BINARY_FALLBACK")
        bl2_assert_design(sample)
        scaled <- bl2_numerical_scale(sample$data, sample$formula)
        diagnostics$numerical_transform <- scaled$transform
        d <- scaled$data
        fit_one <- function(formula) bl2_capture(nnet::multinom(formula, data = d, contrasts = sample$contrasts, 
            na.action = stats::na.fail, Hess = TRUE, model = FALSE, summ = 0L, trace = FALSE, maxit = as.integer(maxit), 
            decay = 0, reltol = 1e-10, MaxNWts = max(1000L, as.integer(4L * (sample$columns + 1L)))))
        full <- fit_one(sample$formula)
        null <- fit_one(sample$null_formula)
        diagnostics$warnings_full <- full$warnings
        diagnostics$warnings_null <- null$warnings
        diagnostics$warnings <- unique(c(full$warnings, null$warnings))
        global$warning <- paste(diagnostics$warnings, collapse = "; ")
        full_error <- inherits(full$value, "error")
        null_error <- inherits(null$value, "error")
        if (!full_error) {
            diagnostics$full <- bl2_multinom_diagnostics(full$value, sample$formula, d, sample$contrasts, scaled$transform)
            diagnostics$raw_coefficients <- diagnostics$full$raw_coefficients
            global$convergence_full <- diagnostics$full$convergence
            global$logLik_full <- diagnostics$full$logLik
            target <- diagnostics$full$raw_coefficients
            target <- target[target$term == "border_relative", , drop = FALSE]
            ix <- match(contrasts$outcome_level, target$outcome_level)
            contrasts$raw_beta <- target$raw_beta[ix]
            contrasts$raw_SE <- target$raw_SE[ix]
            contrasts$p_raw <- target$raw_p[ix]
            contrasts$raw_relative_odds_ratio <- exp(contrasts$raw_beta)
            contrasts$raw_CI_low <- exp(contrasts$raw_beta - stats::qnorm(0.975) * contrasts$raw_SE)
            contrasts$raw_CI_high <- exp(contrasts$raw_beta + stats::qnorm(0.975) * contrasts$raw_SE)
        }
        else diagnostics$error_full <- conditionMessage(full$value)
        if (!null_error) {
            diagnostics$null <- bl2_multinom_diagnostics(null$value, sample$null_formula, d, sample$contrasts, 
                scaled$transform)
            diagnostics$null$issues <- setdiff(diagnostics$null$issues, "INVALID_TARGET_ESTIMATE")
            diagnostics$null$target_numeric_reliable <- NA
            diagnostics$null$reliable <- !length(diagnostics$null$issues)
            global$convergence_null <- diagnostics$null$convergence
            global$logLik_null <- diagnostics$null$logLik
        }
        else diagnostics$error_null <- conditionMessage(null$value)
        if (full_error || null_error) 
            stop(paste(c(if (full_error) paste("FULL_FIT_ERROR", diagnostics$error_full), if (null_error) paste("NULL_FIT_ERROR", 
                diagnostics$error_null)), collapse = "; "), call. = FALSE)
        diagnostics$warnings <- unique(c(diagnostics$warnings, diagnostics$full$vcov_warnings, diagnostics$null$vcov_warnings))
        global$warning <- paste(diagnostics$warnings, collapse = "; ")
        diagnostics$same_complete_sample <- identical(diagnostics$full$sample_rows, diagnostics$null$sample_rows) && 
            identical(diagnostics$full$sample_rows, as.character(sample$rows)) && diagnostics$full$nobs == diagnostics$null$nobs && 
            diagnostics$full$nobs == nrow(d)
        global$df_raw <- diagnostics$full$edf - diagnostics$null$edf
        global$LR_raw <- 2 * (global$logLik_full - global$logLik_null)
        diagnostics$LR_roundoff_tolerance <- max(1e-06, 100 * .Machine$double.eps * (1 + abs(global$logLik_full) + 
            abs(global$logLik_null)))
        if (is.finite(global$LR_raw) && global$LR_raw >= -diagnostics$LR_roundoff_tolerance && is.finite(global$df_raw) && 
            global$df_raw > 0) 
            global$p_raw <- stats::pchisq(max(0, global$LR_raw), df = global$df_raw, lower.tail = FALSE)
        issues <- c(diagnostics$full$issues, diagnostics$null$issues)
        if (!diagnostics$same_complete_sample) 
            issues <- c(issues, "FULL_NULL_SAMPLE_MISMATCH")
        if (!isTRUE(global$df_raw == 2L)) 
            issues <- c(issues, "GLOBAL_LRT_NOT_TWO_DF")
        if (!is.finite(global$p_raw) || !is.finite(global$LR_raw) || global$LR_raw < -diagnostics$LR_roundoff_tolerance) 
            issues <- c(issues, "INVALID_LIKELIHOOD_RATIO")
        if (any(!is.finite(as.matrix(contrasts[, c("raw_relative_odds_ratio", "raw_CI_low", "raw_CI_high")])) | 
            contrasts$raw_CI_low <= 0)) 
            issues <- c(issues, "NONFINITE_OR_UNDERFLOWED_TARGET_OR_CI")
        diagnostics$reliability_issues <- unique(issues)
        bl2_assert(!length(issues), paste(unique(issues), collapse = "; "))
        global$LR <- max(0, global$LR_raw)
        global$p <- global$p_raw
        global$numeric_reliable <- TRUE
        global$status <- if (length(diagnostics$warnings)) 
            "FIT_WARNING_REVIEW"
        else "OK"
        contrasts$beta <- contrasts$raw_beta
        contrasts$SE <- contrasts$raw_SE
        contrasts$relative_odds_ratio <- contrasts$raw_relative_odds_ratio
        contrasts$CI_low <- contrasts$raw_CI_low
        contrasts$CI_high <- contrasts$raw_CI_high
        contrasts$p <- contrasts$p_raw
    }, error = function(e) {
        global$error <<- conditionMessage(e)
        global$status <<- "NOT_ESTIMABLE"
    })
    contrasts$status <- global$status
    contrasts$error <- global$error
    contrasts$warning <- global$warning
    list(global = global, contrasts = contrasts, diagnostics = diagnostics)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        bl2_capture <- function(expr) {
            warnings <- character()
            value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }), error = identity)
            list(value = value, warnings = unique(warnings))
        }
        bl2_assert <- function(ok, message) {
            if (!isTRUE(ok)) stop(message, call. = FALSE)
        }
        bl2_status_factor <- function(x) {
            factor(as.character(x), levels = c("0", "1", "2"))
        }
        bl2_formula <- function() {
            y ~ border_relative + age + sex + ethnicity + center + assessment_year + factor(S)
        }
        bl2_contrasts <- function(formula, data) {
            mf <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
            is_factor <- vapply(mf, is.factor, logical(1))
            is_factor[1L] <- FALSE
            lapply(mf[is_factor], function(x) stats::contr.treatment(levels(x), base = 1L))
        }
        bl2_sample <- function(data, outcome) {
            bl2_assert(is.data.frame(data), "MODEL_READY_DATA_FRAME_REQUIRED")
            required <- c("eid", outcome, "border_relative", "age", "sex", "ethnicity", "center", "assessment_year", 
                "S")
            missing <- setdiff(required, names(data))
            bl2_assert(!length(missing), paste("MISSING_MODEL_FIELDS", paste(missing, collapse = ",")))
            bl2_assert(!anyNA(data$eid) && !anyDuplicated(as.character(data$eid)), "DUPLICATE_OR_MISSING_ID")
            bl2_assert(is.numeric(data$border_relative), "NUMERIC_BORDER_RELATIVE_REQUIRED")
            bl2_assert(is.numeric(data$age), "LINEAR_NUMERIC_AGE_REQUIRED")
            vars <- setdiff(required, "eid")
            keep <- stats::complete.cases(data[, vars, drop = FALSE])
            for (v in vars) if (is.numeric(data[[v]])) keep <- keep & is.finite(data[[v]])
            rows <- which(keep)
            d <- droplevels(as.data.frame(data[rows, , drop = FALSE]))
            d$y <- data[[outcome]][rows]
            if (is.factor(d$y)) d$y <- factor(as.character(d$y), levels = levels(data[[outcome]]))
            rownames(d) <- as.character(rows)
            bl2_assert(nrow(d) > 0L, "NO_COMPLETE_CASES")
            bl2_assert(all(d$border_relative >= 0 & d$border_relative <= 1), "BORDER_RELATIVE_OUTSIDE_0_1")
            s <- suppressWarnings(as.numeric(as.character(d$S)))
            bl2_assert(all(is.finite(s) & s %in% 2:8), "S2_TO_8_REQUIRED")
            bl2_assert(length(unique(d$border_relative)) > 1L, "CONSTANT_EXPOSURE")
            formula <- bl2_formula()
            labels <- attr(stats::terms(formula), "term.labels")
            removed <- character()
            for (term in setdiff(labels, "border_relative")) {
                v <- all.vars(stats::as.formula(paste("~", term)))
                if (length(v) == 1L && length(unique(d[[v]])) <= 1L) removed <- c(removed, term)
            }
            if (length(removed)) formula <- stats::reformulate(setdiff(labels, removed), response = "y")
            null_formula <- stats::reformulate(setdiff(attr(stats::terms(formula), "term.labels"), "border_relative"), 
                response = "y")
            contrasts <- bl2_contrasts(formula, d)
            x <- stats::model.matrix(formula, d, contrasts.arg = contrasts)
            x0 <- stats::model.matrix(null_formula, d, contrasts.arg = contrasts)
            list(data = d, formula = formula, null_formula = null_formula, rows = rows, eids = as.character(d$eid), 
                removed_terms = removed, contrasts = contrasts, x = x, x0 = x0, rank = qr(x)$rank, columns = ncol(x), 
                rank_null = qr(x0)$rank, columns_null = ncol(x0))
        }
        bl2_global_template <- function(data, outcome, model) {
            data.frame(outcome = outcome, model = model, exposure = "border_relative", N_input = nrow(data), N = NA_integer_, 
                N_excluded = NA_integer_, status = "NOT_ESTIMABLE", error = "", warning = "", formula = paste(deparse(bl2_formula()), 
                  collapse = " "), null_formula = "", removed_terms = "", rank = NA_integer_, columns = NA_integer_, 
                rank_null = NA_integer_, columns_null = NA_integer_, beta = NA_real_, SE = NA_real_, CI_low = NA_real_, 
                CI_high = NA_real_, beta_SD = NA_real_, SE_SD = NA_real_, CI_low_SD = NA_real_, CI_high_SD = NA_real_, 
                outcome_SD = NA_real_, df_residual = NA_real_, p = NA_real_, p_raw = NA_real_, raw_beta = NA_real_, 
                raw_SE = NA_real_, raw_CI_low = NA_real_, raw_CI_high = NA_real_, raw_beta_SD = NA_real_, raw_SE_SD = NA_real_, 
                raw_CI_low_SD = NA_real_, raw_CI_high_SD = NA_real_, LR = NA_real_, LR_raw = NA_real_, df = if (model == 
                  "lm") 1L else 2L, df_raw = NA_real_, logLik_full = NA_real_, logLik_null = NA_real_, convergence_full = NA_integer_, 
                convergence_null = NA_integer_, numeric_reliable = FALSE, inference = "ASSOCIATION_NONCAUSAL", 
                stringsAsFactors = FALSE)
        }
        bl2_add_sample_metadata <- function(global, sample) {
            global$N <- nrow(sample$data)
            global$N_excluded <- global$N_input - global$N
            global$formula <- paste(deparse(sample$formula), collapse = " ")
            global$null_formula <- paste(deparse(sample$null_formula), collapse = " ")
            global$removed_terms <- paste(sample$removed_terms, collapse = ";")
            for (nm in c("rank", "columns", "rank_null", "columns_null")) global[[nm]] <- sample[[nm]]
            global
        }
        bl2_sample_diagnostics <- function(sample) {
            list(sample_rows = sample$rows, sample_eids = sample$eids, removed_constant_terms = sample$removed_terms, 
                full_design_columns = colnames(sample$x), null_design_columns = colnames(sample$x0), rank_full = sample$rank, 
                columns_full = sample$columns, rank_null = sample$rank_null, columns_null = sample$columns_null)
        }
        bl2_assert_design <- function(sample) {
            bl2_assert(all(is.finite(sample$x)) && all(is.finite(sample$x0)), "NONFINITE_DESIGN")
            bl2_assert(sample$rank == sample$columns && sample$rank_null == sample$columns_null, "RANK_DEFICIENT_DESIGN")
            bl2_assert("border_relative" %in% colnames(sample$x), "TARGET_COLUMN_MISSING")
        }
        bl2_fit_bmi <- function(data, outcome = "bmi") {
            global <- bl2_global_template(data, outcome, "lm")
            diagnostics <- list(raw_coefficients = data.frame(), warnings = character())
            tryCatch({
                bl2_assert(outcome %in% names(data) && is.numeric(data[[outcome]]), "NUMERIC_BMI_OUTCOME_REQUIRED")
                sample <- bl2_sample(data, outcome)
                global <- bl2_add_sample_metadata(global, sample)
                diagnostics <- c(diagnostics, bl2_sample_diagnostics(sample))
                bl2_assert_design(sample)
                bl2_assert(nrow(sample$data) > sample$columns, "NO_RESIDUAL_DEGREES_OF_FREEDOM")
                global$outcome_SD <- stats::sd(sample$data$y)
                bl2_assert(is.finite(global$outcome_SD) && global$outcome_SD > 0, "CONSTANT_BMI_OUTCOME")
                captured <- bl2_capture(stats::lm(sample$formula, data = sample$data, contrasts = sample$contrasts, 
                  na.action = stats::na.fail))
                diagnostics$warnings <- captured$warnings
                global$warning <- paste(captured$warnings, collapse = "; ")
                if (inherits(captured$value, "error")) stop(conditionMessage(captured$value), call. = FALSE)
                fit <- captured$value
                captured_summary <- bl2_capture(summary(fit))
                diagnostics$warnings <- unique(c(diagnostics$warnings, captured_summary$warnings))
                global$warning <- paste(diagnostics$warnings, collapse = "; ")
                if (inherits(captured_summary$value, "error")) stop(conditionMessage(captured_summary$value), call. = FALSE)
                cf <- captured_summary$value$coefficients
                diagnostics$raw_coefficients <- data.frame(term = rownames(cf), raw_beta = cf[, 1L], raw_SE = cf[, 
                  2L], raw_t = cf[, 3L], raw_p = cf[, 4L], row.names = NULL)
                global$df_residual <- fit$df.residual
                global$raw_beta <- cf["border_relative", "Estimate"]
                global$raw_SE <- cf["border_relative", "Std. Error"]
                global$p_raw <- cf["border_relative", "Pr(>|t|)"]
                crit <- stats::qt(0.975, df = fit$df.residual)
                global$raw_CI_low <- global$raw_beta - crit * global$raw_SE
                global$raw_CI_high <- global$raw_beta + crit * global$raw_SE
                global$raw_beta_SD <- global$raw_beta/global$outcome_SD
                global$raw_SE_SD <- global$raw_SE/global$outcome_SD
                global$raw_CI_low_SD <- global$raw_CI_low/global$outcome_SD
                global$raw_CI_high_SD <- global$raw_CI_high/global$outcome_SD
                diagnostics$target_numeric_reliable <- all(is.finite(c(global$raw_beta, global$raw_SE, global$p_raw, 
                  global$raw_CI_low, global$raw_CI_high))) && global$raw_SE > 0
                bl2_assert(all(is.finite(cf)) && diagnostics$target_numeric_reliable, "NONFINITE_OR_INVALID_LM_INFERENCE")
                global$beta <- global$raw_beta
                global$SE <- global$raw_SE
                global$CI_low <- global$raw_CI_low
                global$CI_high <- global$raw_CI_high
                global$beta_SD <- global$raw_beta_SD
                global$SE_SD <- global$raw_SE_SD
                global$CI_low_SD <- global$raw_CI_low_SD
                global$CI_high_SD <- global$raw_CI_high_SD
                global$p <- global$p_raw
                global$numeric_reliable <- TRUE
                global$status <- if (length(diagnostics$warnings)) "FIT_WARNING_REVIEW" else "OK"
            }, error = function(e) {
                global$error <<- conditionMessage(e)
                global$status <<- "NOT_ESTIMABLE"
            })
            list(global = global, contrasts = data.frame(), diagnostics = diagnostics)
        }
        bl2_numerical_scale <- function(data, formula) {
            d <- data
            transform <- data.frame(variable = character(), center = numeric(), scale = numeric(), original_unit = character(), 
                stringsAsFactors = FALSE)
            used <- all.vars(formula)
            for (v in intersect(c("age", "assessment_year"), used)) {
                if (!is.numeric(d[[v]])) next
                center <- min(d[[v]])
                scale <- diff(range(d[[v]]))
                bl2_assert(is.finite(scale) && scale > 0, "INVALID_NUMERICAL_SCALE")
                d[[v]] <- (d[[v]] - center)/scale
                transform <- rbind(transform, data.frame(variable = v, center = center, scale = scale, original_unit = "year"))
            }
            list(data = d, transform = transform)
        }
        bl2_backmap <- function(coef_names, transform, classes = 2L) {
            tmat <- diag(length(coef_names))
            dimnames(tmat) <- list(coef_names, coef_names)
            for (i in seq_len(nrow(transform))) {
                v <- transform$variable[i]
                bl2_assert(v %in% coef_names && "(Intercept)" %in% coef_names, "AFFINE_BACKMAP_COLUMN_MISSING")
                tmat[v, v] <- 1/transform$scale[i]
                tmat["(Intercept)", v] <- -transform$center[i]/transform$scale[i]
            }
            list(one_class = tmat, all_classes = kronecker(diag(classes), tmat))
        }
        bl2_hessian_diagnostics <- function(hessian) {
            result <- list(finite = FALSE, dimension = NROW(hessian), rank = NA_integer_, minimum_eigenvalue = NA_real_, 
                maximum_eigenvalue = NA_real_, condition_number = Inf, tolerance = sqrt(.Machine$double.eps), reliable = FALSE, 
                error = "")
            if (!is.matrix(hessian) || !nrow(hessian) || nrow(hessian) != ncol(hessian) || any(!is.finite(hessian))) {
                result$error <- "NONFINITE_OR_INVALID_HESSIAN"
                return(result)
            }
            result$finite <- TRUE
            captured <- bl2_capture(eigen((hessian + t(hessian))/2, symmetric = TRUE, only.values = TRUE)$values)
            if (inherits(captured$value, "error")) {
                result$error <- conditionMessage(captured$value)
                return(result)
            }
            ev <- captured$value
            result$minimum_eigenvalue <- min(ev)
            result$maximum_eigenvalue <- max(ev)
            result$rank <- sum(ev > max(ev) * result$tolerance)
            if (min(ev) > 0) result$condition_number <- max(ev)/min(ev)
            result$reliable <- max(ev) > 0 && min(ev) > 0 && result$rank == nrow(hessian)
            if (!result$reliable) result$error <- "SINGULAR_OR_SEVERELY_ILL_CONDITIONED_HESSIAN"
            result
        }
        bl2_multinom_diagnostics <- function(fit, formula, data, contrasts, transform) {
            cf <- stats::coef(fit)
            bl2_assert(is.matrix(cf) && identical(rownames(cf), c("1", "2")), "THREE_CLASS_COEFFICIENT_LAYOUT_MISMATCH")
            backmap <- bl2_backmap(colnames(cf), transform, nrow(cf))
            raw_cf <- cf %*% t(backmap$one_class)
            vc <- bl2_capture(stats::vcov(fit))
            raw_vcov <- if (!inherits(vc$value, "error")) backmap$all_classes %*% vc$value %*% t(backmap$all_classes) else matrix(NA_real_, 
                nrow(cf) * ncol(cf), nrow(cf) * ncol(cf))
            variances <- diag(raw_vcov)
            raw_se <- matrix(sqrt(ifelse(variances >= 0, variances, NA_real_)), nrow = nrow(cf), byrow = TRUE, 
                dimnames = dimnames(cf))
            raw_table <- data.frame(outcome_level = rep(rownames(raw_cf), each = ncol(raw_cf)), term = rep(colnames(raw_cf), 
                times = nrow(raw_cf)), raw_beta = as.vector(t(raw_cf)), raw_SE = as.vector(t(raw_se)), stringsAsFactors = FALSE)
            raw_table$raw_z <- raw_table$raw_beta/raw_table$raw_SE
            raw_table$raw_p <- 2 * stats::pnorm(abs(raw_table$raw_z), lower.tail = FALSE)
            target <- raw_table[raw_table$term == "border_relative", , drop = FALSE]
            hessian <- bl2_hessian_diagnostics(fit$Hessian)
            prob <- fit$fitted.values
            x <- stats::model.matrix(formula, data, contrasts.arg = contrasts)
            response <- sapply(c("1", "2"), function(k) as.numeric(data$y == k))
            score <- crossprod(x, response - prob[, c("1", "2"), drop = FALSE])
            target_reliable <- nrow(target) == 2L && all(is.finite(target$raw_beta)) && all(is.finite(target$raw_SE) & 
                target$raw_SE > 0) && all(is.finite(target$raw_p))
            issues <- character()
            if (!identical(fit$convergence, 0L)) issues <- c(issues, "MULTINOM_MAXIT_OR_NONCONVERGENCE")
            if (!all(is.finite(cf)) || !all(is.finite(raw_vcov)) || !all(is.finite(prob)) || !is.finite(as.numeric(stats::logLik(fit)))) issues <- c(issues, 
                "NONFINITE_MULTINOM_FIT")
            if (!hessian$reliable) issues <- c(issues, hessian$error)
            if (!target_reliable) issues <- c(issues, "INVALID_TARGET_ESTIMATE")
            if (inherits(vc$value, "error")) issues <- c(issues, paste("VCOV_ERROR", conditionMessage(vc$value)))
            list(convergence = fit$convergence, hessian = hessian, raw_coefficients = raw_table, raw_vcov = raw_vcov, 
                coefficient_backmap = backmap$one_class, target_numeric_reliable = target_reliable, max_abs_score_internal = max(abs(score)), 
                target_score_internal = if ("border_relative" %in% rownames(score)) score["border_relative", ] else rep(NA_real_, 
                  2L), minimum_fitted_probability = min(prob), maximum_fitted_probability = max(prob), extreme_probability_flag = any(prob < 
                  1e-12 | prob > 1 - 1e-12), vcov_warnings = vc$warnings, reliable = !length(issues), issues = unique(issues), 
                nobs = nrow(prob), sample_rows = rownames(prob), logLik = as.numeric(stats::logLik(fit)), edf = fit$edf)
        }
        bl2_contrast_template <- function(outcome) {
            data.frame(outcome = outcome, exposure = "border_relative", reference_level = "0", outcome_level = c("1", 
                "2"), contrast = c("Previous_vs_Never", "Current_vs_Never"), N = NA_integer_, beta = NA_real_, 
                SE = NA_real_, relative_odds_ratio = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, 
                raw_beta = NA_real_, raw_SE = NA_real_, raw_relative_odds_ratio = NA_real_, raw_CI_low = NA_real_, 
                raw_CI_high = NA_real_, p_raw = NA_real_, status = "NOT_ESTIMABLE", error = "", warning = "", stringsAsFactors = FALSE)
        }
        bl2_fit_multinomial <- function(data, outcome, maxit = 1000L) {
            global <- bl2_global_template(data, outcome, "nnet::multinom")
            contrasts <- bl2_contrast_template(outcome)
            diagnostics <- list(raw_coefficients = data.frame(), warnings = character(), class_levels = c("0", 
                "1", "2"), reference_level = "0", optimizer = list(maxit = maxit, decay = 0, reltol = 1e-10, Hess = TRUE, 
                summ = 0L), numeric_reparameterization = "affine nuisance-only age/assessment_year; backmapped")
            tryCatch({
                bl2_assert(requireNamespace("nnet", quietly = TRUE), "NNET_PACKAGE_REQUIRED")
                bl2_assert(outcome %in% names(data) && is.factor(data[[outcome]]) && identical(levels(data[[outcome]]), 
                  c("0", "1", "2")), "VERIFIED_THREE_LEVEL_FACTOR_0_REFERENCE_REQUIRED")
                bl2_assert(length(maxit) == 1L && is.finite(maxit) && maxit >= 1 && maxit == as.integer(maxit), 
                  "INVALID_MAXIT")
                sample <- bl2_sample(data, outcome)
                global <- bl2_add_sample_metadata(global, sample)
                contrasts$N <- global$N
                diagnostics <- c(diagnostics, bl2_sample_diagnostics(sample))
                diagnostics$class_counts <- table(sample$data$y)
                bl2_assert(all(diagnostics$class_counts > 0L), "MISSING_OUTCOME_CLASS_NO_BINARY_FALLBACK")
                bl2_assert_design(sample)
                scaled <- bl2_numerical_scale(sample$data, sample$formula)
                diagnostics$numerical_transform <- scaled$transform
                d <- scaled$data
                fit_one <- function(formula) bl2_capture(nnet::multinom(formula, data = d, contrasts = sample$contrasts, 
                  na.action = stats::na.fail, Hess = TRUE, model = FALSE, summ = 0L, trace = FALSE, maxit = as.integer(maxit), 
                  decay = 0, reltol = 1e-10, MaxNWts = max(1000L, as.integer(4L * (sample$columns + 1L)))))
                full <- fit_one(sample$formula)
                null <- fit_one(sample$null_formula)
                diagnostics$warnings_full <- full$warnings
                diagnostics$warnings_null <- null$warnings
                diagnostics$warnings <- unique(c(full$warnings, null$warnings))
                global$warning <- paste(diagnostics$warnings, collapse = "; ")
                full_error <- inherits(full$value, "error")
                null_error <- inherits(null$value, "error")
                if (!full_error) {
                  diagnostics$full <- bl2_multinom_diagnostics(full$value, sample$formula, d, sample$contrasts, 
                    scaled$transform)
                  diagnostics$raw_coefficients <- diagnostics$full$raw_coefficients
                  global$convergence_full <- diagnostics$full$convergence
                  global$logLik_full <- diagnostics$full$logLik
                  target <- diagnostics$full$raw_coefficients
                  target <- target[target$term == "border_relative", , drop = FALSE]
                  ix <- match(contrasts$outcome_level, target$outcome_level)
                  contrasts$raw_beta <- target$raw_beta[ix]
                  contrasts$raw_SE <- target$raw_SE[ix]
                  contrasts$p_raw <- target$raw_p[ix]
                  contrasts$raw_relative_odds_ratio <- exp(contrasts$raw_beta)
                  contrasts$raw_CI_low <- exp(contrasts$raw_beta - stats::qnorm(0.975) * contrasts$raw_SE)
                  contrasts$raw_CI_high <- exp(contrasts$raw_beta + stats::qnorm(0.975) * contrasts$raw_SE)
                } else diagnostics$error_full <- conditionMessage(full$value)
                if (!null_error) {
                  diagnostics$null <- bl2_multinom_diagnostics(null$value, sample$null_formula, d, sample$contrasts, 
                    scaled$transform)
                  diagnostics$null$issues <- setdiff(diagnostics$null$issues, "INVALID_TARGET_ESTIMATE")
                  diagnostics$null$target_numeric_reliable <- NA
                  diagnostics$null$reliable <- !length(diagnostics$null$issues)
                  global$convergence_null <- diagnostics$null$convergence
                  global$logLik_null <- diagnostics$null$logLik
                } else diagnostics$error_null <- conditionMessage(null$value)
                if (full_error || null_error) stop(paste(c(if (full_error) paste("FULL_FIT_ERROR", diagnostics$error_full), 
                  if (null_error) paste("NULL_FIT_ERROR", diagnostics$error_null)), collapse = "; "), call. = FALSE)
                diagnostics$warnings <- unique(c(diagnostics$warnings, diagnostics$full$vcov_warnings, diagnostics$null$vcov_warnings))
                global$warning <- paste(diagnostics$warnings, collapse = "; ")
                diagnostics$same_complete_sample <- identical(diagnostics$full$sample_rows, diagnostics$null$sample_rows) && 
                  identical(diagnostics$full$sample_rows, as.character(sample$rows)) && diagnostics$full$nobs == 
                  diagnostics$null$nobs && diagnostics$full$nobs == nrow(d)
                global$df_raw <- diagnostics$full$edf - diagnostics$null$edf
                global$LR_raw <- 2 * (global$logLik_full - global$logLik_null)
                diagnostics$LR_roundoff_tolerance <- max(1e-06, 100 * .Machine$double.eps * (1 + abs(global$logLik_full) + 
                  abs(global$logLik_null)))
                if (is.finite(global$LR_raw) && global$LR_raw >= -diagnostics$LR_roundoff_tolerance && is.finite(global$df_raw) && 
                  global$df_raw > 0) global$p_raw <- stats::pchisq(max(0, global$LR_raw), df = global$df_raw, lower.tail = FALSE)
                issues <- c(diagnostics$full$issues, diagnostics$null$issues)
                if (!diagnostics$same_complete_sample) issues <- c(issues, "FULL_NULL_SAMPLE_MISMATCH")
                if (!isTRUE(global$df_raw == 2L)) issues <- c(issues, "GLOBAL_LRT_NOT_TWO_DF")
                if (!is.finite(global$p_raw) || !is.finite(global$LR_raw) || global$LR_raw < -diagnostics$LR_roundoff_tolerance) issues <- c(issues, 
                  "INVALID_LIKELIHOOD_RATIO")
                if (any(!is.finite(as.matrix(contrasts[, c("raw_relative_odds_ratio", "raw_CI_low", "raw_CI_high")])) | 
                  contrasts$raw_CI_low <= 0)) issues <- c(issues, "NONFINITE_OR_UNDERFLOWED_TARGET_OR_CI")
                diagnostics$reliability_issues <- unique(issues)
                bl2_assert(!length(issues), paste(unique(issues), collapse = "; "))
                global$LR <- max(0, global$LR_raw)
                global$p <- global$p_raw
                global$numeric_reliable <- TRUE
                global$status <- if (length(diagnostics$warnings)) "FIT_WARNING_REVIEW" else "OK"
                contrasts$beta <- contrasts$raw_beta
                contrasts$SE <- contrasts$raw_SE
                contrasts$relative_odds_ratio <- contrasts$raw_relative_odds_ratio
                contrasts$CI_low <- contrasts$raw_CI_low
                contrasts$CI_high <- contrasts$raw_CI_high
                contrasts$p <- contrasts$p_raw
            }, error = function(e) {
                global$error <<- conditionMessage(e)
                global$status <<- "NOT_ESTIMABLE"
            })
            contrasts$status <- global$status
            contrasts$error <- global$error
            contrasts$warning <- global$warning
            list(global = global, contrasts = contrasts, diagnostics = diagnostics)
        }
    }), envir = context)
}

