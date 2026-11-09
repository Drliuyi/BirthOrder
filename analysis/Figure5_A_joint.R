# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

fd_validate_members <- function(ctx, model) {
    ids <- as.character(unlist(model$generation_order))
    ans <- lapply(ids, function(id) {
        ca <- ctx$features[[id]]
        s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
        fd_assert(is.list(ca) && nrow(s) == 1L && identical(ctx$variables[[id]]$status, "AVAILABLE_FROZEN"), paste("UNVERIFIED_MEMBER", 
            id))
        type <- as.character(s$type)
        fd_assert(type %in% c("continuous", "binary", "categorical", "ordinal", "nominal"), paste("UNSUPPORTED_TYPE", 
            id, type))
        fd_assert(ca$column %in% names(ctx$d_full), paste("MEMBER_COLUMN_MISSING", id))
        if (type == "continuous") 
            fd_assert(identical(ca$prepared_transform, "identity"), paste("UNREVIEWED_CONTINUOUS_TRANSFORM", id))
        x <- ctx$d_full[[ca$column]]
        list(indicator_id = id, column = ca$column, type = type, levels = if (type != "continuous") levels(x) else NULL, 
            lower = if (is.null(ca$physical_lower_bound)) -Inf else ca$physical_lower_bound, upper = if (is.null(ca$physical_upper_bound)) Inf else ca$physical_upper_bound, 
            prepared_transform = ca$prepared_transform, technical_covariates = ca$technical_covariates)
    })
    names(ans) <- ids
    fd_assert(!anyDuplicated(vapply(ans, `[[`, "", "column")), "DUPLICATE_OR_DERIVATIVE_COLUMNS_REQUIRE_REVIEW")
    ans
}

fd_support <- function(d, members) {
    for (m in members) {
        x <- d[[m$column]]
        if (m$type == "continuous") 
            fd_assert(is.numeric(x) && all(is.finite(x) & x >= m$lower & x <= m$upper), paste("PHYSICAL_SUPPORT_FAILURE_NO_CLIPPING", 
                m$indicator_id))
        else fd_assert(!anyNA(x) && all(as.character(x) %in% m$levels), paste("CATEGORY_SUPPORT_FAILURE", m$indicator_id))
    }
    TRUE
}

fd_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED", count_rule = "NONE") {
    fd_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0), "INVALID_DOMAIN_SELECTION_WEIGHT")
    fd_support(d, members)
    fitted <- list()
    previous <- character()
    warnings <- character()
    for (m in members) {
        count <- m$indicator_id %in% c("moderate_days", "walking_days") && identical(count_rule, "FIXED_SUPPORT_MULTINOMIAL_0_7")
        if (count) {
            fd_assert(all(d[[m$column]] %in% 0:7), "COUNT_OBSERVATION_OUTSIDE_FIXED_0_7")
            d[[m$column]] <- factor(as.character(d[[m$column]]), levels = as.character(0:7))
            m$levels <- as.character(0:7)
        }
        form <- stats::reformulate(unique(c(base, previous)), response = m$column)
        X <- stats::model.matrix(stats::delete.response(stats::terms(form)), d)
        fd_assert(qr(X)$rank == ncol(X), paste("A_DESIGN_RANK_DEFICIENT", m$indicator_id))
        if (m$type == "continuous" && !count) {
            fit <- stats::lm.wfit(X, d[[m$column]], w)
            fd_assert(fit$rank == ncol(X) && all(is.finite(fit$coefficients)), "A_CONTINUOUS_FIT_FAILED")
            residual <- d[[m$column]] - as.numeric(X %*% fit$coefficients)
            order <- order(residual)
            f <- list(kind = "continuous", coefficients = fit$coefficients, residuals = residual[order], cumulative = c(0, 
                cumsum(w[order])/sum(w)), cumulative_squared_mass = c(0, cumsum((w[order]/sum(w))^2)), rule = continuous_rule)
        }
        else {
            fd_assert(length(m$levels) >= 2L && all(table(d[[m$column]]) > 0), paste("FROZEN_CATEGORY_ABSENT", 
                m$indicator_id))
            got <- fd_capture(nnet::multinom(form, data = d, weights = w, trace = FALSE, maxit = 1000L, MaxNWts = 100000L, 
                model = FALSE))
            fit <- got$value
            warnings <- c(warnings, got$warnings)
            fd_assert(identical(fit$convergence, 0L) && all(is.finite(stats::coef(fit))), paste("A_MULTINOM_NONCONVERGED", 
                m$indicator_id))
            f <- list(kind = if (count) "count" else "categorical", fit = fit)
        }
        f$formula <- form
        f$columns <- colnames(X)
        f$member <- m
        fitted[[m$indicator_id]] <- f
        previous <- c(previous, m$column)
        if (count) 
            d[[m$column]] <- as.numeric(as.character(d[[m$column]]))
    }
    list(models = fitted, members = members, warnings = unique(warnings), conditional_order = names(fitted), continuous_rule = continuous_rule, 
        residual_assumption = "Independent conditional errors after conditioning on all earlier mediator values; empirical weighted residual distribution")
}

fd_probs <- function(f, d) {
    p <- predict(f$fit, newdata = d, type = "probs")
    lev <- f$member$levels
    if (is.null(dim(p))) {
        if (length(lev) == 2L) 
            p <- cbind(1 - p, p)
        else p <- matrix(p, nrow = 1L)
        colnames(p) <- lev
    }
    fd_assert(setequal(colnames(p), lev), "A_PREDICTION_CATEGORIES_CHANGED")
    p <- p[, lev, drop = FALSE]
    fd_assert(nrow(p) == nrow(d) && all(is.finite(p) & p >= 0 & p <= 1) && max(abs(rowSums(p) - 1)) < 1e-08, "A_INVALID_PROBABILITIES")
    p
}

fd_draw_a <- function(a, C, b, u, audit = NULL) {
    fd_assert(length(b) %in% c(1L, nrow(C)) && all(is.finite(b) & b >= 0 & b <= 1) && is.matrix(u) && nrow(u) == 
        nrow(C) && ncol(u) == length(a$models) && all(u >= 0 & u < 1), "INVALID_JOINT_DRAW_RANDOM_INPUT")
    d <- C
    d$B <- b
    for (j in seq_along(a$models)) {
        f <- a$models[[j]]
        m <- f$member
        X <- model.matrix(delete.response(terms(f$formula)), d)
        fd_assert(identical(colnames(X), f$columns), "A_CONDITIONAL_DESIGN_DRIFT")
        if (f$kind == "continuous") {
            mu <- as.numeric(X %*% f$coefficients)
            q <- u[, j]
            if (f$rule == "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL") {
                low <- findInterval(m$lower - mu, f$residuals, left.open = TRUE)
                high <- findInterval(m$upper - mu, f$residuals)
                lo <- f$cumulative[low + 1L]
                hi <- f$cumulative[high + 1L]
                square_mass <- f$cumulative_squared_mass[high + 1L] - f$cumulative_squared_mass[low + 1L]
                pool_ESS <- ifelse(square_mass > 0, (hi - lo)^2/square_mass, 0)
                if (is.function(audit)) 
                  audit(m$indicator_id, high - low, hi - lo, pool_ESS)
                fd_assert(all(hi > lo), paste("CONDITIONAL_SUPPORT_EMPTY", m$indicator_id))
                q <- lo + q * (hi - lo)
            }
            else fd_assert(f$rule == "EMPIRICAL_RESIDUAL_UNTRUNCATED", "UNREVIEWED_CONTINUOUS_CONDITIONAL_LAW")
            ix <- findInterval(q, f$cumulative, rightmost.closed = TRUE)
            fd_assert(all(ix >= 1 & ix <= length(f$residuals)), "RESIDUAL_INVERSE_CDF_INVALID")
            d[[m$column]] <- mu + f$residuals[ix]
        }
        else {
            p <- fd_probs(f, d)
            cum <- p
            if (ncol(p) > 1L) 
                for (k in 2:ncol(p)) cum[, k] <- cum[, k - 1L] + p[, k]
            ix <- 1L + rowSums(u[, j] >= cum[, -ncol(cum), drop = FALSE])
            d[[m$column]] <- factor(m$levels[ix], levels = m$levels)
            if (f$kind == "count") 
                d[[m$column]] <- as.numeric(as.character(d[[m$column]]))
        }
    }
    fd_support(d, a$members)
    d
}

fd_design <- function(d, members, base, reference = NULL) {
    discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
    f <- reformulate(unique(c(base, discrete)))
    X <- model.matrix(f, d)
    continuous <- Filter(function(m) m$type == "continuous", members)
    if (is.null(reference)) 
        reference <- lapply(continuous, function(m) {
            q <- unname(quantile(d[[m$column]], c(0.05, 0.35, 0.65, 0.95), type = 7))
            fd_assert(all(diff(q) > 0), paste("FROZEN_SPLINE_KNOTS_NOT_DISTINCT", m$indicator_id))
            list(knots = q[2:3], boundary = q[c(1, 4)])
        })
    for (id in names(continuous)) {
        m <- continuous[[id]]
        k <- reference[[id]]
        z <- splines::ns(d[[m$column]], knots = k$knots, Boundary.knots = k$boundary, intercept = FALSE)
        colnames(z) <- paste0(m$column, "__ns", seq_len(ncol(z)))
        X <- cbind(X, z)
    }
    list(X = X, reference = reference)
}

fd_distribution_state <- function(members) {
    p <- length(members)
    list(weight = 0, rows = 0L, sum = numeric(p), sumsq = numeric(p), cross = matrix(0, p, p), zero = numeric(p), 
        frequencies = setNames(lapply(members, function(m) numeric()), names(members)), diagnostic_category_mappings = list(), 
        zero_applicable = rep(TRUE, p), nominal_pairs = list())
}

fd_distribution_update <- function(state, d, w, members) {
    columns <- lapply(seq_along(members), function(j) {
        m <- members[[j]]
        x <- d[[m$column]]
        if (is.factor(x) || is.character(x)) {
            lev <- if (!is.null(m$levels)) 
                m$levels
            else levels(x)
            fd_assert(length(lev) > 0L && !anyNA(match(as.character(x), lev)), "DIAGNOSTIC_CATEGORY_OUTSIDE_FROZEN_LEVELS")
            values <- suppressWarnings(as.numeric(lev))
            numeric_labels <- all(is.finite(values))
            if (!numeric_labels) 
                values <- seq_along(lev)
            mapping <- list(rule = if (numeric_labels) "NUMERIC_LEVEL_LABELS_DIAGNOSTICS_ONLY" else "FROZEN_LEVEL_INDEX_DIAGNOSTICS_ONLY_NOT_ORDINAL_MODEL", 
                label_to_value = as.list(setNames(values, lev)), A_Y_encoding_unchanged = TRUE)
            key <- names(members)[j]
            if (!is.null(state$diagnostic_category_mappings[[key]])) 
                fd_assert(identical(state$diagnostic_category_mappings[[key]], mapping), "DIAGNOSTIC_CATEGORY_MAPPING_DRIFT")
            state$diagnostic_category_mappings[[key]] <<- mapping
            state$zero_applicable[j] <<- numeric_labels
            as.numeric(values[match(as.character(x), lev)])
        }
        else as.numeric(x)
    })
    num <- matrix(unlist(columns, use.names = FALSE), nrow = nrow(d), ncol = length(members))
    fd_assert(all(is.finite(num)), "NONFINITE_DIAGNOSTIC_VALUES")
    state$weight <- state$weight + sum(w)
    state$rows <- state$rows + nrow(d)
    state$sum <- state$sum + colSums(num * w)
    state$sumsq <- state$sumsq + colSums(num^2 * w)
    state$cross <- state$cross + crossprod(num * sqrt(w))
    state$zero <- state$zero + colSums((num == 0) * w)
    categorical <- names(Filter(function(m) m$type != "continuous" || m$indicator_id %in% c("moderate_days", "walking_days"), 
        members))
    merge_counts <- function(old, fresh) {
        keys <- union(names(old), names(fresh))
        ans <- setNames(numeric(length(keys)), keys)
        ans[names(old)] <- old
        ans[names(fresh)] <- ans[names(fresh)] + fresh
        ans
    }
    for (id in categorical) {
        v <- as.character(d[[members[[id]]$column]])
        state$frequencies[[id]] <- merge_counts(state$frequencies[[id]], tapply(w, v, sum))
    }
    if (length(categorical) >= 2L) 
        for (pair in combn(categorical, 2, simplify = FALSE)) {
            key <- paste(pair, collapse = "__WITH__")
            v <- paste(as.character(d[[members[[pair[1]]]$column]]), as.character(d[[members[[pair[2]]]$column]]), 
                sep = "|")
            state$nominal_pairs[[key]] <- merge_counts(state$nominal_pairs[[key]], tapply(w, v, sum))
        }
    state
}

fd_distribution_finish <- function(state, members) {
    mean <- state$sum/state$weight
    covariance <- state$cross/state$weight - tcrossprod(mean)
    variance <- state$sumsq/state$weight - mean^2
    denominator <- sqrt(outer(variance, variance))
    correlation <- covariance/denominator
    dimnames(correlation) <- list(names(members), names(members))
    list(N_simulated_rows = state$rows, total_weight = state$weight, weighted_mean = setNames(mean, names(members)), 
        weighted_variance = setNames(variance, names(members)), weighted_true_zero_fraction = setNames(ifelse(state$zero_applicable, 
            state$zero/state$weight, NA_real_), names(members)), diagnostic_category_mappings = state$diagnostic_category_mappings, 
        weighted_category_frequencies = lapply(state$frequencies, function(v) v/state$weight), numeric_code_correlation_descriptive_only_not_ordinal_coding = correlation, 
        weighted_nominal_pair_joint_frequencies = lapply(state$nominal_pairs, function(v) v/state$weight))
}

fd_integrate <- function(a, C, w, predict_fun, K, seed, chunk = 2048L, support_log = "") {
    fd_assert(K >= 1L && K == as.integer(K), "K_MUST_BE_POSITIVE_INTEGER")
    endpoints <- c("250.2", "411")
    names3 <- c("R00", "R10", "R11")
    sums <- matrix(0, 2, 3, dimnames = list(endpoints, names3))
    per_k <- array(0, c(K, 2, 3))
    if (nzchar(support_log)) 
        fd_assert(!file.exists(support_log), "PRESERVE_SUPPORT_AUDIT")
    distributions <- setNames(lapply(1:3, function(j) fd_distribution_state(a$members)), c("B0", "B1", "OBSERVED_B_RECONSTRUCTION"))
    observed <- fd_distribution_update(fd_distribution_state(a$members), C, w, a$members)
    support_hist <- list()
    audit_for <- function(scenario, k, start) function(id, counts, probability_mass, pool_ESS) {
        key <- paste(scenario, id, sep = ":")
        old <- support_hist[[key]]
        fresh <- table(counts)
        keys <- union(names(old), names(fresh))
        ans <- setNames(numeric(length(keys)), keys)
        ans[names(old)] <- old
        ans[names(fresh)] <- ans[names(fresh)] + as.numeric(fresh)
        support_hist[[key]] <<- ans
        if (nzchar(support_log)) 
            cat(jsonlite::toJSON(list(scenario = scenario, replicate = k, row_start = start, indicator_id = id, 
                n = length(counts), candidate_min = min(counts), candidate_max = max(counts), candidate_mean = mean(counts), 
                candidate_quantiles = unname(quantile(counts, c(0, 0.01, 0.5, 0.99, 1))), empty_count = sum(counts == 
                  0 | probability_mass <= 0), minimum_retained_probability_mass = min(probability_mass), retained_probability_mass_quantiles = unname(quantile(probability_mass, 
                  c(0, 0.01, 0.5, 0.99, 1))), effective_residual_pool_size_quantiles = unname(quantile(pool_ESS, 
                  c(0, 0.01, 0.5, 0.99, 1)))), auto_unbox = TRUE, digits = 17), "\n", file = support_log, append = TRUE)
    }
    fd_rng(seed, {
        for (k in seq_len(K)) for (start in seq.int(1L, nrow(C), by = chunk)) {
            ii <- seq.int(start, min(nrow(C), start + chunk - 1L))
            u <- matrix(runif(length(ii) * length(a$models)), nrow = length(ii), byrow = TRUE)
            m0 <- fd_draw_a(a, C[ii, , drop = FALSE], 0, u, audit_for("B0", k, start))
            m1 <- fd_draw_a(a, C[ii, , drop = FALSE], 1, u, audit_for("B1", k, start))
            mo <- fd_draw_a(a, C[ii, , drop = FALSE], C$B[ii], u, audit_for("OBSERVED_B_RECONSTRUCTION", k, start))
            for (name in names(distributions)) distributions[[name]] <- fd_distribution_update(distributions[[name]], 
                switch(name, B0 = m0, B1 = m1, OBSERVED_B_RECONSTRUCTION = mo), w[ii], a$members)
            scenarios <- list(R00 = m0, R10 = m0, R11 = m1)
            scenarios$R10$B <- 1
            for (e in seq_along(endpoints)) for (s in seq_along(names3)) {
                p <- predict_fun(scenarios[[s]], endpoints[e])
                fd_assert(length(p) == length(ii) && all(is.finite(p) & p >= 0 & p <= 1), "CIF_NONFINITE_OR_INVALID_NO_CLIPPING")
                per_k[k, e, s] <- per_k[k, e, s] + sum(w[ii] * p)/sum(w)
            }
        }
    })
    sums <- apply(per_k, c(2, 3), mean)
    dimnames(sums) <- list(endpoints, names3)
    IE <- sums[, "R11"] - sums[, "R10"]
    DE <- sums[, "R10"] - sums[, "R00"]
    TE <- sums[, "R11"] - sums[, "R00"]
    fd_assert(max(abs(TE - DE - IE)) <= 1e-12, "RISK_COMPONENT_CLOSURE_FAILED")
    list(risks = sums, IE = IE, DE = DE, TE = TE, delta = unname(IE["411"] - IE["250.2"]), K = K, seed = seed, 
        closure = max(abs(TE - DE - IE)), mc_precision_passed = FALSE, distribution_diagnostics = c(list(OBSERVED = fd_distribution_finish(observed, 
            a$members)), lapply(distributions, fd_distribution_finish, members = a$members), list(interpretation = "Observed-B reconstruction audits the working distribution; B0/B1 counterfactual differences are not automatic failures. Diagnostics are not proof of model fit.")), 
        support_histograms = support_hist, integration = "WEIGHTED_FULL_VECTOR_MONTE_CARLO_SHARED_RANDOM_NUMBERS_BOTH_B_AND_BOTH_DISEASES")
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        fd_validate_members <- function(ctx, model) {
            ids <- as.character(unlist(model$generation_order))
            ans <- lapply(ids, function(id) {
                ca <- ctx$features[[id]]
                s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
                fd_assert(is.list(ca) && nrow(s) == 1L && identical(ctx$variables[[id]]$status, "AVAILABLE_FROZEN"), 
                  paste("UNVERIFIED_MEMBER", id))
                type <- as.character(s$type)
                fd_assert(type %in% c("continuous", "binary", "categorical", "ordinal", "nominal"), paste("UNSUPPORTED_TYPE", 
                  id, type))
                fd_assert(ca$column %in% names(ctx$d_full), paste("MEMBER_COLUMN_MISSING", id))
                if (type == "continuous") fd_assert(identical(ca$prepared_transform, "identity"), paste("UNREVIEWED_CONTINUOUS_TRANSFORM", 
                  id))
                x <- ctx$d_full[[ca$column]]
                list(indicator_id = id, column = ca$column, type = type, levels = if (type != "continuous") levels(x) else NULL, 
                  lower = if (is.null(ca$physical_lower_bound)) -Inf else ca$physical_lower_bound, upper = if (is.null(ca$physical_upper_bound)) Inf else ca$physical_upper_bound, 
                  prepared_transform = ca$prepared_transform, technical_covariates = ca$technical_covariates)
            })
            names(ans) <- ids
            fd_assert(!anyDuplicated(vapply(ans, `[[`, "", "column")), "DUPLICATE_OR_DERIVATIVE_COLUMNS_REQUIRE_REVIEW")
            ans
        }
        fd_support <- function(d, members) {
            for (m in members) {
                x <- d[[m$column]]
                if (m$type == "continuous") fd_assert(is.numeric(x) && all(is.finite(x) & x >= m$lower & x <= m$upper), 
                  paste("PHYSICAL_SUPPORT_FAILURE_NO_CLIPPING", m$indicator_id)) else fd_assert(!anyNA(x) && all(as.character(x) %in% 
                  m$levels), paste("CATEGORY_SUPPORT_FAILURE", m$indicator_id))
            }
            TRUE
        }
        fd_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED", count_rule = "NONE") {
            fd_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0), "INVALID_DOMAIN_SELECTION_WEIGHT")
            fd_support(d, members)
            fitted <- list()
            previous <- character()
            warnings <- character()
            for (m in members) {
                count <- m$indicator_id %in% c("moderate_days", "walking_days") && identical(count_rule, "FIXED_SUPPORT_MULTINOMIAL_0_7")
                if (count) {
                  fd_assert(all(d[[m$column]] %in% 0:7), "COUNT_OBSERVATION_OUTSIDE_FIXED_0_7")
                  d[[m$column]] <- factor(as.character(d[[m$column]]), levels = as.character(0:7))
                  m$levels <- as.character(0:7)
                }
                form <- stats::reformulate(unique(c(base, previous)), response = m$column)
                X <- stats::model.matrix(stats::delete.response(stats::terms(form)), d)
                fd_assert(qr(X)$rank == ncol(X), paste("A_DESIGN_RANK_DEFICIENT", m$indicator_id))
                if (m$type == "continuous" && !count) {
                  fit <- stats::lm.wfit(X, d[[m$column]], w)
                  fd_assert(fit$rank == ncol(X) && all(is.finite(fit$coefficients)), "A_CONTINUOUS_FIT_FAILED")
                  residual <- d[[m$column]] - as.numeric(X %*% fit$coefficients)
                  order <- order(residual)
                  f <- list(kind = "continuous", coefficients = fit$coefficients, residuals = residual[order], 
                    cumulative = c(0, cumsum(w[order])/sum(w)), cumulative_squared_mass = c(0, cumsum((w[order]/sum(w))^2)), 
                    rule = continuous_rule)
                } else {
                  fd_assert(length(m$levels) >= 2L && all(table(d[[m$column]]) > 0), paste("FROZEN_CATEGORY_ABSENT", 
                    m$indicator_id))
                  got <- fd_capture(nnet::multinom(form, data = d, weights = w, trace = FALSE, maxit = 1000L, MaxNWts = 100000L, 
                    model = FALSE))
                  fit <- got$value
                  warnings <- c(warnings, got$warnings)
                  fd_assert(identical(fit$convergence, 0L) && all(is.finite(stats::coef(fit))), paste("A_MULTINOM_NONCONVERGED", 
                    m$indicator_id))
                  f <- list(kind = if (count) "count" else "categorical", fit = fit)
                }
                f$formula <- form
                f$columns <- colnames(X)
                f$member <- m
                fitted[[m$indicator_id]] <- f
                previous <- c(previous, m$column)
                if (count) d[[m$column]] <- as.numeric(as.character(d[[m$column]]))
            }
            list(models = fitted, members = members, warnings = unique(warnings), conditional_order = names(fitted), 
                continuous_rule = continuous_rule, residual_assumption = "Independent conditional errors after conditioning on all earlier mediator values; empirical weighted residual distribution")
        }
        fd_probs <- function(f, d) {
            p <- predict(f$fit, newdata = d, type = "probs")
            lev <- f$member$levels
            if (is.null(dim(p))) {
                if (length(lev) == 2L) p <- cbind(1 - p, p) else p <- matrix(p, nrow = 1L)
                colnames(p) <- lev
            }
            fd_assert(setequal(colnames(p), lev), "A_PREDICTION_CATEGORIES_CHANGED")
            p <- p[, lev, drop = FALSE]
            fd_assert(nrow(p) == nrow(d) && all(is.finite(p) & p >= 0 & p <= 1) && max(abs(rowSums(p) - 1)) < 1e-08, 
                "A_INVALID_PROBABILITIES")
            p
        }
        fd_draw_a <- function(a, C, b, u, audit = NULL) {
            fd_assert(length(b) %in% c(1L, nrow(C)) && all(is.finite(b) & b >= 0 & b <= 1) && is.matrix(u) && nrow(u) == 
                nrow(C) && ncol(u) == length(a$models) && all(u >= 0 & u < 1), "INVALID_JOINT_DRAW_RANDOM_INPUT")
            d <- C
            d$B <- b
            for (j in seq_along(a$models)) {
                f <- a$models[[j]]
                m <- f$member
                X <- model.matrix(delete.response(terms(f$formula)), d)
                fd_assert(identical(colnames(X), f$columns), "A_CONDITIONAL_DESIGN_DRIFT")
                if (f$kind == "continuous") {
                  mu <- as.numeric(X %*% f$coefficients)
                  q <- u[, j]
                  if (f$rule == "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL") {
                    low <- findInterval(m$lower - mu, f$residuals, left.open = TRUE)
                    high <- findInterval(m$upper - mu, f$residuals)
                    lo <- f$cumulative[low + 1L]
                    hi <- f$cumulative[high + 1L]
                    square_mass <- f$cumulative_squared_mass[high + 1L] - f$cumulative_squared_mass[low + 1L]
                    pool_ESS <- ifelse(square_mass > 0, (hi - lo)^2/square_mass, 0)
                    if (is.function(audit)) audit(m$indicator_id, high - low, hi - lo, pool_ESS)
                    fd_assert(all(hi > lo), paste("CONDITIONAL_SUPPORT_EMPTY", m$indicator_id))
                    q <- lo + q * (hi - lo)
                  } else fd_assert(f$rule == "EMPIRICAL_RESIDUAL_UNTRUNCATED", "UNREVIEWED_CONTINUOUS_CONDITIONAL_LAW")
                  ix <- findInterval(q, f$cumulative, rightmost.closed = TRUE)
                  fd_assert(all(ix >= 1 & ix <= length(f$residuals)), "RESIDUAL_INVERSE_CDF_INVALID")
                  d[[m$column]] <- mu + f$residuals[ix]
                } else {
                  p <- fd_probs(f, d)
                  cum <- p
                  if (ncol(p) > 1L) for (k in 2:ncol(p)) cum[, k] <- cum[, k - 1L] + p[, k]
                  ix <- 1L + rowSums(u[, j] >= cum[, -ncol(cum), drop = FALSE])
                  d[[m$column]] <- factor(m$levels[ix], levels = m$levels)
                  if (f$kind == "count") d[[m$column]] <- as.numeric(as.character(d[[m$column]]))
                }
            }
            fd_support(d, a$members)
            d
        }
        fd_design <- function(d, members, base, reference = NULL) {
            discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
            f <- reformulate(unique(c(base, discrete)))
            X <- model.matrix(f, d)
            continuous <- Filter(function(m) m$type == "continuous", members)
            if (is.null(reference)) reference <- lapply(continuous, function(m) {
                q <- unname(quantile(d[[m$column]], c(0.05, 0.35, 0.65, 0.95), type = 7))
                fd_assert(all(diff(q) > 0), paste("FROZEN_SPLINE_KNOTS_NOT_DISTINCT", m$indicator_id))
                list(knots = q[2:3], boundary = q[c(1, 4)])
            })
            for (id in names(continuous)) {
                m <- continuous[[id]]
                k <- reference[[id]]
                z <- splines::ns(d[[m$column]], knots = k$knots, Boundary.knots = k$boundary, intercept = FALSE)
                colnames(z) <- paste0(m$column, "__ns", seq_len(ncol(z)))
                X <- cbind(X, z)
            }
            list(X = X, reference = reference)
        }
        fd_distribution_state <- function(members) {
            p <- length(members)
            list(weight = 0, rows = 0L, sum = numeric(p), sumsq = numeric(p), cross = matrix(0, p, p), zero = numeric(p), 
                frequencies = setNames(lapply(members, function(m) numeric()), names(members)), diagnostic_category_mappings = list(), 
                zero_applicable = rep(TRUE, p), nominal_pairs = list())
        }
        fd_distribution_update <- function(state, d, w, members) {
            columns <- lapply(seq_along(members), function(j) {
                m <- members[[j]]
                x <- d[[m$column]]
                if (is.factor(x) || is.character(x)) {
                  lev <- if (!is.null(m$levels)) m$levels else levels(x)
                  fd_assert(length(lev) > 0L && !anyNA(match(as.character(x), lev)), "DIAGNOSTIC_CATEGORY_OUTSIDE_FROZEN_LEVELS")
                  values <- suppressWarnings(as.numeric(lev))
                  numeric_labels <- all(is.finite(values))
                  if (!numeric_labels) values <- seq_along(lev)
                  mapping <- list(rule = if (numeric_labels) "NUMERIC_LEVEL_LABELS_DIAGNOSTICS_ONLY" else "FROZEN_LEVEL_INDEX_DIAGNOSTICS_ONLY_NOT_ORDINAL_MODEL", 
                    label_to_value = as.list(setNames(values, lev)), A_Y_encoding_unchanged = TRUE)
                  key <- names(members)[j]
                  if (!is.null(state$diagnostic_category_mappings[[key]])) fd_assert(identical(state$diagnostic_category_mappings[[key]], 
                    mapping), "DIAGNOSTIC_CATEGORY_MAPPING_DRIFT")
                  state$diagnostic_category_mappings[[key]] <<- mapping
                  state$zero_applicable[j] <<- numeric_labels
                  as.numeric(values[match(as.character(x), lev)])
                } else as.numeric(x)
            })
            num <- matrix(unlist(columns, use.names = FALSE), nrow = nrow(d), ncol = length(members))
            fd_assert(all(is.finite(num)), "NONFINITE_DIAGNOSTIC_VALUES")
            state$weight <- state$weight + sum(w)
            state$rows <- state$rows + nrow(d)
            state$sum <- state$sum + colSums(num * w)
            state$sumsq <- state$sumsq + colSums(num^2 * w)
            state$cross <- state$cross + crossprod(num * sqrt(w))
            state$zero <- state$zero + colSums((num == 0) * w)
            categorical <- names(Filter(function(m) m$type != "continuous" || m$indicator_id %in% c("moderate_days", 
                "walking_days"), members))
            merge_counts <- function(old, fresh) {
                keys <- union(names(old), names(fresh))
                ans <- setNames(numeric(length(keys)), keys)
                ans[names(old)] <- old
                ans[names(fresh)] <- ans[names(fresh)] + fresh
                ans
            }
            for (id in categorical) {
                v <- as.character(d[[members[[id]]$column]])
                state$frequencies[[id]] <- merge_counts(state$frequencies[[id]], tapply(w, v, sum))
            }
            if (length(categorical) >= 2L) for (pair in combn(categorical, 2, simplify = FALSE)) {
                key <- paste(pair, collapse = "__WITH__")
                v <- paste(as.character(d[[members[[pair[1]]]$column]]), as.character(d[[members[[pair[2]]]$column]]), 
                  sep = "|")
                state$nominal_pairs[[key]] <- merge_counts(state$nominal_pairs[[key]], tapply(w, v, sum))
            }
            state
        }
        fd_distribution_finish <- function(state, members) {
            mean <- state$sum/state$weight
            covariance <- state$cross/state$weight - tcrossprod(mean)
            variance <- state$sumsq/state$weight - mean^2
            denominator <- sqrt(outer(variance, variance))
            correlation <- covariance/denominator
            dimnames(correlation) <- list(names(members), names(members))
            list(N_simulated_rows = state$rows, total_weight = state$weight, weighted_mean = setNames(mean, names(members)), 
                weighted_variance = setNames(variance, names(members)), weighted_true_zero_fraction = setNames(ifelse(state$zero_applicable, 
                  state$zero/state$weight, NA_real_), names(members)), diagnostic_category_mappings = state$diagnostic_category_mappings, 
                weighted_category_frequencies = lapply(state$frequencies, function(v) v/state$weight), numeric_code_correlation_descriptive_only_not_ordinal_coding = correlation, 
                weighted_nominal_pair_joint_frequencies = lapply(state$nominal_pairs, function(v) v/state$weight))
        }
        fd_integrate <- function(a, C, w, predict_fun, K, seed, chunk = 2048L, support_log = "") {
            fd_assert(K >= 1L && K == as.integer(K), "K_MUST_BE_POSITIVE_INTEGER")
            endpoints <- c("250.2", "411")
            names3 <- c("R00", "R10", "R11")
            sums <- matrix(0, 2, 3, dimnames = list(endpoints, names3))
            per_k <- array(0, c(K, 2, 3))
            if (nzchar(support_log)) fd_assert(!file.exists(support_log), "PRESERVE_SUPPORT_AUDIT")
            distributions <- setNames(lapply(1:3, function(j) fd_distribution_state(a$members)), c("B0", "B1", 
                "OBSERVED_B_RECONSTRUCTION"))
            observed <- fd_distribution_update(fd_distribution_state(a$members), C, w, a$members)
            support_hist <- list()
            audit_for <- function(scenario, k, start) function(id, counts, probability_mass, pool_ESS) {
                key <- paste(scenario, id, sep = ":")
                old <- support_hist[[key]]
                fresh <- table(counts)
                keys <- union(names(old), names(fresh))
                ans <- setNames(numeric(length(keys)), keys)
                ans[names(old)] <- old
                ans[names(fresh)] <- ans[names(fresh)] + as.numeric(fresh)
                support_hist[[key]] <<- ans
                if (nzchar(support_log)) cat(jsonlite::toJSON(list(scenario = scenario, replicate = k, row_start = start, 
                  indicator_id = id, n = length(counts), candidate_min = min(counts), candidate_max = max(counts), 
                  candidate_mean = mean(counts), candidate_quantiles = unname(quantile(counts, c(0, 0.01, 0.5, 
                    0.99, 1))), empty_count = sum(counts == 0 | probability_mass <= 0), minimum_retained_probability_mass = min(probability_mass), 
                  retained_probability_mass_quantiles = unname(quantile(probability_mass, c(0, 0.01, 0.5, 0.99, 
                    1))), effective_residual_pool_size_quantiles = unname(quantile(pool_ESS, c(0, 0.01, 0.5, 0.99, 
                    1)))), auto_unbox = TRUE, digits = 17), "\n", file = support_log, append = TRUE)
            }
            fd_rng(seed, {
                for (k in seq_len(K)) for (start in seq.int(1L, nrow(C), by = chunk)) {
                  ii <- seq.int(start, min(nrow(C), start + chunk - 1L))
                  u <- matrix(runif(length(ii) * length(a$models)), nrow = length(ii), byrow = TRUE)
                  m0 <- fd_draw_a(a, C[ii, , drop = FALSE], 0, u, audit_for("B0", k, start))
                  m1 <- fd_draw_a(a, C[ii, , drop = FALSE], 1, u, audit_for("B1", k, start))
                  mo <- fd_draw_a(a, C[ii, , drop = FALSE], C$B[ii], u, audit_for("OBSERVED_B_RECONSTRUCTION", 
                    k, start))
                  for (name in names(distributions)) distributions[[name]] <- fd_distribution_update(distributions[[name]], 
                    switch(name, B0 = m0, B1 = m1, OBSERVED_B_RECONSTRUCTION = mo), w[ii], a$members)
                  scenarios <- list(R00 = m0, R10 = m0, R11 = m1)
                  scenarios$R10$B <- 1
                  for (e in seq_along(endpoints)) for (s in seq_along(names3)) {
                    p <- predict_fun(scenarios[[s]], endpoints[e])
                    fd_assert(length(p) == length(ii) && all(is.finite(p) & p >= 0 & p <= 1), "CIF_NONFINITE_OR_INVALID_NO_CLIPPING")
                    per_k[k, e, s] <- per_k[k, e, s] + sum(w[ii] * p)/sum(w)
                  }
                }
            })
            sums <- apply(per_k, c(2, 3), mean)
            dimnames(sums) <- list(endpoints, names3)
            IE <- sums[, "R11"] - sums[, "R10"]
            DE <- sums[, "R10"] - sums[, "R00"]
            TE <- sums[, "R11"] - sums[, "R00"]
            fd_assert(max(abs(TE - DE - IE)) <= 1e-12, "RISK_COMPONENT_CLOSURE_FAILED")
            list(risks = sums, IE = IE, DE = DE, TE = TE, delta = unname(IE["411"] - IE["250.2"]), K = K, seed = seed, 
                closure = max(abs(TE - DE - IE)), mc_precision_passed = FALSE, distribution_diagnostics = c(list(OBSERVED = fd_distribution_finish(observed, 
                  a$members)), lapply(distributions, fd_distribution_finish, members = a$members), list(interpretation = "Observed-B reconstruction audits the working distribution; B0/B1 counterfactual differences are not automatic failures. Diagnostics are not proof of model fit.")), 
                support_histograms = support_hist, integration = "WEIGHTED_FULL_VECTOR_MONTE_CARLO_SHARED_RANDOM_NUMBERS_BOTH_B_AND_BOTH_DISEASES")
        }
    }), envir = context)
}

