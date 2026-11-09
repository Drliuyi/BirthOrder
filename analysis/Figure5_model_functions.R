# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

c39_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

c39_basis <- function(x, ref, column) {
    c39_assert(is.numeric(x) && all(is.finite(x)), paste("NONFINITE_Y_COORDINATE", column))
    if (identical(ref$mode, "constant")) {
        c39_assert(all(x == ref$constant), paste("UNESTIMABLE_CONSTANT_COORDINATE_CHANGED", column))
        z <- matrix(numeric(), nrow = length(x), ncol = 0L)
    }
    else if (identical(ref$mode, "linear")) {
        z <- matrix(x, ncol = 1L)
    }
    else {
        c39_assert(identical(ref$mode, "ns"), "UNREVIEWED_SPLINE_MODE")
        z <- splines::ns(x, knots = ref$knots, Boundary.knots = ref$boundary, intercept = FALSE)
    }
    colnames(z) <- ref$column_names
    z
}

c39_spline_reference <- function(x, column) {
    c39_assert(is.numeric(x) && length(x) && all(is.finite(x)), "INVALID_SPLINE_INPUT")
    q <- unname(stats::quantile(x, c(0.05, 0.35, 0.65, 0.95), type = 7))
    boundary <- q[c(1, 4)]
    knots <- sort(unique(q[2:3]))
    knots <- knots[knots > boundary[1] & knots < boundary[2]]
    constant <- length(unique(x)) == 1L
    mode <- if (constant) 
        "constant"
    else if (boundary[1] == boundary[2] || !length(knots)) 
        "linear"
    else "ns"
    df <- switch(mode, constant = 0L, linear = 1L, ns = length(knots) + 1L)
    list(mode = mode, knots = knots, boundary = boundary, Boundary.knots = boundary, quantile_probabilities = c(0.05, 
        0.35, 0.65, 0.95), original_quantiles = q, effective_df = df, constant = if (constant) x[1L] else NULL, 
        constant_handling = if (constant) "FIXED_COORDINATE_OMIT_EQUIVALENT_DESIGN_ASSERT_AT_PREDICTION" else "NONE", 
        column_names = if (df) paste0(column, "__ns", seq_len(df)) else character())
}

c39_design <- function(d, members, base, reference = NULL) {
    discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
    form <- stats::reformulate(unique(c(base, discrete)), env = baseenv())
    X <- stats::model.matrix(form, d)
    continuous <- Filter(function(m) m$type == "continuous", members)
    if (is.null(reference)) {
        reference <- lapply(continuous, function(m) c39_spline_reference(d[[m$column]], m$column))
        attr(reference, "base_discrete") <- list(formula = form, columns = colnames(X), factor_levels = lapply(d[unique(c(base, 
            discrete))], function(x) if (is.factor(x)) levels(x) else NULL))
    }
    else {
        c39_assert(identical(names(reference), names(continuous)), "FROZEN_SPLINE_MEMBERS_DRIFT")
        saved <- attr(reference, "base_discrete")
        c39_assert(is.list(saved) && identical(colnames(X), saved$columns), "FROZEN_BASE_DISCRETE_DESIGN_DRIFT")
        for (column in names(saved$factor_levels)) if (!is.null(saved$factor_levels[[column]])) 
            c39_assert(is.factor(d[[column]]) && identical(levels(d[[column]]), saved$factor_levels[[column]]), 
                paste("FROZEN_Y_FACTOR_LEVELS_DRIFT", column))
    }
    for (id in names(continuous)) {
        m <- continuous[[id]]
        X <- cbind(X, c39_basis(d[[m$column]], reference[[id]], m$column))
    }
    c39_assert(all(is.finite(X)), "NONFINITE_Y_DESIGN")
    if (is.null(attr(reference, "full_colnames"))) 
        attr(reference, "full_colnames") <- colnames(X)
    c39_assert(identical(colnames(X), attr(reference, "full_colnames")), "FROZEN_Y_DESIGN_COLUMNS_DRIFT")
    list(X = X, reference = reference, fullcolnames = colnames(X))
}

c39_support <- function(d, members) {
    for (m in members) {
        x <- d[[m$column]]
        if (m$type == "continuous") 
            c39_assert(is.numeric(x) && all(is.finite(x) & x >= m$lower & x <= m$upper), paste("PHYSICAL_SUPPORT_FAILURE_NO_CLIPPING", 
                m$indicator_id))
        else c39_assert(!anyNA(x) && all(as.character(x) %in% m$levels), paste("CATEGORY_SUPPORT_FAILURE", m$indicator_id))
    }
    ids <- vapply(members, `[[`, "", "indicator_id")
    if (all(c("blood_pressure__DBP", "blood_pressure__SBP") %in% ids)) {
        dbp <- members[[match("blood_pressure__DBP", ids)]]$column
        sbp <- members[[match("blood_pressure__SBP", ids)]]$column
        c39_assert(all(d[[sbp]] > d[[dbp]]), "SBP_NOT_STRICTLY_ABOVE_DBP_NO_REWRITE")
    }
    TRUE
}

c39_a_design <- function(d, base, previous, reference = NULL) {
    if (is.null(reference)) {
        constants <- previous[vapply(d[previous], function(x) is.numeric(x) && length(unique(x)) == 1L, logical(1))]
        fixed <- lapply(d[constants], function(x) x[1L])
        form <- stats::reformulate(unique(c(base, setdiff(previous, constants))), env = baseenv())
        X <- stats::model.matrix(form, d)
        reference <- list(formula = form, columns = colnames(X), fixed = fixed)
    }
    else {
        for (column in names(reference$fixed)) c39_assert(all(d[[column]] == reference$fixed[[column]]), paste("A_FIXED_PREDECESSOR_CHANGED", 
            column))
        X <- stats::model.matrix(reference$formula, d)
        c39_assert(identical(colnames(X), reference$columns), "A_CONDITIONAL_DESIGN_DRIFT")
    }
    c39_assert(all(is.finite(X)), "A_NONFINITE_DESIGN")
    list(X = X, reference = reference)
}

c39_linear_fit <- function(X, y, w, id, rule) {
    if (length(unique(y)) == 1L) 
        return(list(kind = "constant", value = y[1L], support = "FIXED_AMOUNT_SUPPORT_IN_THIS_SAMPLE", rule = rule))
    c39_assert(qr(X)$rank == ncol(X), paste("A_DESIGN_RANK_DEFICIENT", id))
    fit <- stats::lm.wfit(X, y, w)
    c39_assert(fit$rank == ncol(X) && all(is.finite(fit$coefficients)), paste("A_CONTINUOUS_FIT_FAILED", id))
    residual <- y - as.numeric(X %*% fit$coefficients)
    ix <- order(residual)
    mass <- w[ix]/sum(w)
    cumulative <- c(0, cumsum(mass))
    cumulative[length(cumulative)] <- 1
    list(kind = "continuous", coefficients = fit$coefficients, residuals = residual[ix], cumulative = cumulative, 
        cumulative_squared_mass = c(0, cumsum(mass^2)), rule = rule)
}

c39_logistic_fit <- function(X, positive, w, id) {
    if (length(unique(positive)) == 1L) 
        return(list(kind = "constant", probability = as.numeric(positive[1L]), support = "FIXED_ZERO_POSITIVE_STATUS_IN_THIS_SAMPLE", 
            warnings = character()))
    c39_assert(qr(X)$rank == ncol(X), paste("A_ZERO_STATUS_DESIGN_RANK_DEFICIENT", id))
    warnings <- character()
    fit <- withCallingHandlers(stats::glm.fit(X, as.numeric(positive), weights = w, family = stats::binomial(), 
        control = stats::glm.control(maxit = 1000L)), warning = function(e) {
        warnings <<- c(warnings, conditionMessage(e))
        invokeRestart("muffleWarning")
    })
    c39_assert(isTRUE(fit$converged) && fit$rank == ncol(X) && all(is.finite(fit$coefficients)), paste("A_ZERO_STATUS_LOGISTIC_FAILED", 
        id))
    list(kind = "logistic", coefficients = fit$coefficients, iterations = fit$iter, converged = fit$converged, 
        rank = fit$rank, design_columns = ncol(X), boundary = fit$boundary, deviance = fit$deviance, fitted_probability_range = range(fit$fitted.values), 
        warnings = unique(warnings))
}

c39_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", 
    count_rule = "FIXED_SUPPORT_MULTINOMIAL_0_7") {
    c39_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "INVALID_JOINT_SELECTION_WEIGHT")
    c39_assert(continuous_rule %in% c("EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", "EMPIRICAL_RESIDUAL_UNTRUNCATED"), 
        "UNREVIEWED_CONTINUOUS_CONDITIONAL_LAW")
    c39_support(d, members)
    ids <- vapply(members, `[[`, "", "indicator_id")
    if ("blood_pressure__SBP" %in% ids) 
        c39_assert("blood_pressure__DBP" %in% ids && match("blood_pressure__DBP", ids) < match("blood_pressure__SBP", 
            ids), "DBP_MUST_PRECEDE_SBP")
    previous <- character()
    fitted <- list()
    warnings <- character()
    for (m in members) {
        id <- m$indicator_id
        y <- d[[m$column]]
        count <- id %in% c("moderate_days", "walking_days")
        if (count) {
            c39_assert(identical(count_rule, "FIXED_SUPPORT_MULTINOMIAL_0_7"), "ACTIVITY_REQUIRES_FIXED_0_7_MULTINOMIAL")
            c39_assert(all(y %in% 0:7), "COUNT_OBSERVATION_OUTSIDE_FIXED_0_7")
            d[[m$column]] <- factor(as.character(y), levels = as.character(0:7))
            m$levels <- as.character(0:7)
        }
        design <- c39_a_design(d, base, previous)
        X <- design$X
        two_part <- id %in% c("cooked_veg", "dried_fruit", "raw_veg", "computer_hours", "tv_hours")
        if (two_part) {
            c39_assert(m$type == "continuous" && all(y >= 0), paste("INVALID_TWO_PART_OBSERVATIONS", id))
            positive <- y > 0
            status <- c39_logistic_fit(X, positive, w, id)
            amount <- if (any(positive)) 
                c39_linear_fit(X[positive, , drop = FALSE], y[positive], w[positive], id, "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL")
            else NULL
            f <- list(kind = "two_part", status = status, positive = amount, observed_zero_n = sum(!positive), 
                observed_positive_n = sum(positive), weighted_observed_zero_fraction = sum(w[!positive])/sum(w), 
                original_scale = "IDENTITY_FROZEN_UNITS", model_scale = "POSITIVE_IDENTITY_CONDITIONAL_MEAN", rule = "ZERO_ATOM_PLUS_CONDITIONALLY_LEGAL_WEIGHTED_POSITIVE_EMPIRICAL_RESIDUAL")
            warnings <- c(warnings, status$warnings)
        }
        else if (m$type == "continuous" && !count) {
            f <- c39_linear_fit(X, y, w, id, continuous_rule)
            if (id == "blood_pressure__SBP") 
                f$rule <- "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
        }
        else {
            c39_assert(length(m$levels) >= 2L && all(table(factor(as.character(d[[m$column]]), levels = m$levels)) > 
                0), paste("FROZEN_CATEGORY_ABSENT", id))
            c39_assert(qr(X)$rank == ncol(X), paste("A_DESIGN_RANK_DEFICIENT", id))
            form <- stats::reformulate(attr(stats::terms(design$reference$formula), "term.labels"), response = m$column, 
                env = baseenv())
            c39_assert(!".c39_fit_weights" %in% names(d), "RESERVED_INTERNAL_WEIGHT_COLUMN")
            fit_data <- d
            fit_data$.c39_fit_weights <- w
            fit <- withCallingHandlers(nnet::multinom(form, data = fit_data, weights = .c39_fit_weights, trace = FALSE, 
                maxit = 1000L, MaxNWts = 100000L, model = FALSE), warning = function(e) {
                warnings <<- c(warnings, conditionMessage(e))
                invokeRestart("muffleWarning")
            })
            c39_assert(identical(fit$convergence, 0L) && all(is.finite(stats::coef(fit))), paste("A_MULTINOM_NONCONVERGED", 
                id))
            f <- list(kind = if (count) "count" else "categorical", fit = fit)
        }
        f$formula <- design$reference$formula
        f$columns <- colnames(X)
        f$design_reference <- design$reference
        f$member <- m
        f$predecessor_columns <- previous
        fitted[[id]] <- f
        checkpoint <- getOption("c39.a_checkpoint_callback")
        if (is.function(checkpoint)) 
            checkpoint(id, f, length(fitted))
        previous <- c(previous, m$column)
        if (count) 
            d[[m$column]] <- y
    }
    list(models = fitted, members = members, warnings = unique(warnings), conditional_order = names(fitted), continuous_rule = continuous_rule, 
        model_revision = "CORE39_TWO_PART_SPLINE_BP_V2_1", residual_assumption = "Weighted empirical residual law conditional on legal support and every preceding coordinate; positive atoms only for two-part quantity")
}

c39_residual_draw <- function(f, X, u, lower, upper, lower_open, id, audit = NULL) {
    n <- nrow(X)
    lower <- rep_len(lower, n)
    upper <- rep_len(upper, n)
    lower_open <- rep_len(lower_open, n)
    if (identical(f$kind, "constant")) {
        value <- rep(f$value, n)
        legal <- ifelse(lower_open, value > lower, value >= lower) & value <= upper
        c39_assert(all(legal), paste("CONDITIONAL_CONSTANT_SUPPORT_EMPTY", id))
        if (is.function(audit)) 
            audit(id, rep(1L, n), rep(1, n), rep(1, n))
        return(value)
    }
    mu <- as.numeric(X %*% f$coefficients)
    q <- u
    if (identical(f$rule, "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL")) {
        low <- findInterval(lower - mu, f$residuals, left.open = TRUE)
        low[lower_open] <- findInterval((lower - mu)[lower_open], f$residuals)
        high <- findInterval(upper - mu, f$residuals)
        lo <- f$cumulative[low + 1L]
        hi <- f$cumulative[high + 1L]
        mass <- hi - lo
        square_mass <- f$cumulative_squared_mass[high + 1L] - f$cumulative_squared_mass[low + 1L]
        ess <- ifelse(square_mass > 0, mass^2/square_mass, 0)
        if (is.function(audit)) 
            audit(id, high - low, mass, ess)
        c39_assert(all(mass > 0 & high > low), paste("CONDITIONAL_SUPPORT_EMPTY", id))
        q <- lo + u * mass
        ix <- findInterval(q, f$cumulative, rightmost.closed = TRUE)
        upper_rounding <- q == hi & ix > high
        ix[upper_rounding] <- high[upper_rounding]
        c39_assert(all(ix > low & ix <= high), paste("CONDITIONAL_INVERSE_CDF_INVALID", id))
    }
    else {
        c39_assert(identical(f$rule, "EMPIRICAL_RESIDUAL_UNTRUNCATED"), "UNREVIEWED_RESIDUAL_LAW")
        ix <- findInterval(q, f$cumulative, rightmost.closed = TRUE)
    }
    c39_assert(all(ix >= 1L & ix <= length(f$residuals)), "RESIDUAL_INVERSE_CDF_INVALID")
    value <- mu + f$residuals[ix]
    legal <- ifelse(lower_open, value > lower, value >= lower) & value <= upper
    c39_assert(all(is.finite(value) & legal), paste("DRAW_ILLEGAL_NO_CLIPPING", id))
    value
}

c39_draw_a <- function(a, C, b, u, audit = NULL) {
    c39_assert(length(b) %in% c(1L, nrow(C)) && all(is.finite(b) & b >= 0 & b <= 1) && is.matrix(u) && nrow(u) == 
        nrow(C) && ncol(u) == length(a$models) && all(is.finite(u) & u >= 0 & u < 1), "INVALID_JOINT_DRAW_RANDOM_INPUT")
    d <- C
    d$B <- b
    for (j in seq_along(a$models)) {
        f <- a$models[[j]]
        m <- f$member
        id <- m$indicator_id
        X <- c39_a_design(d, NULL, NULL, f$design_reference)$X
        lower <- rep(m$lower, nrow(d))
        lower_open <- rep(FALSE, nrow(d))
        upper <- rep(m$upper, nrow(d))
        if (id == "blood_pressure__SBP") {
            dbp <- a$models[["blood_pressure__DBP"]]$member$column
            lower_open <- d[[dbp]] >= lower
            lower <- pmax(lower, d[[dbp]])
        }
        if (f$kind == "two_part") {
            p <- if (identical(f$status$kind, "constant")) 
                rep(f$status$probability, nrow(d))
            else stats::plogis(as.numeric(X %*% f$status$coefficients))
            c39_assert(all(is.finite(p) & p >= 0 & p <= 1), "INVALID_ZERO_STATUS_PROBABILITY")
            positive <- u[, j] >= 1 - p & p > 0
            value <- numeric(nrow(d))
            if (any(positive)) {
                c39_assert(is.list(f$positive), paste("POSITIVE_STATUS_WITHOUT_AMOUNT_MODEL", id))
                q <- (u[positive, j] - (1 - p[positive]))/p[positive]
                c39_assert(all(q >= 0 & q <= 1), "INVALID_POSITIVE_CONDITIONAL_UNIFORM")
                plower <- pmax(0, lower[positive])
                popen <- lower[positive] <= 0
                value[positive] <- c39_residual_draw(f$positive, X[positive, , drop = FALSE], q, plower, upper[positive], 
                  popen, id, audit)
            }
            d[[m$column]] <- value
        }
        else if (f$kind %in% c("continuous", "constant")) {
            d[[m$column]] <- c39_residual_draw(f, X, u[, j], lower, upper, lower_open, id, audit)
        }
        else {
            p <- fd_probs(f, d)
            cumulative <- p
            if (ncol(p) > 1L) 
                for (k in 2:ncol(p)) cumulative[, k] <- cumulative[, k - 1L] + p[, k]
            ix <- 1L + rowSums(u[, j] >= cumulative[, -ncol(cumulative), drop = FALSE])
            d[[m$column]] <- factor(m$levels[ix], levels = m$levels)
            if (f$kind == "count") 
                d[[m$column]] <- as.numeric(as.character(d[[m$column]]))
        }
    }
    c39_support(d, a$members)
    d
}

c39_install_models <- function(engine) {
    c39_assert(is.environment(engine) && is.function(engine$fd_probs), "ISOLATED_SOURCE_LOCKED_ENGINE_REQUIRED")
    functions <- c("c39_assert", "c39_basis", "c39_spline_reference", "c39_design", "c39_support", "c39_a_design", 
        "c39_linear_fit", "c39_logistic_fit", "c39_fit_a", "c39_residual_draw", "c39_draw_a")
    source_environment <- environment(c39_install_models)
    for (name in functions) {
        f <- get(name, envir = source_environment, inherits = FALSE)
        environment(f) <- engine
        assign(name, f, envir = engine)
    }
    engine$fd_fit_a <- engine$c39_fit_a
    engine$fd_draw_a <- engine$c39_draw_a
    engine$fd_design <- engine$c39_design
    engine$fd_support <- engine$c39_support
    invisible(engine)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        c39_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        c39_basis <- function(x, ref, column) {
            c39_assert(is.numeric(x) && all(is.finite(x)), paste("NONFINITE_Y_COORDINATE", column))
            if (identical(ref$mode, "constant")) {
                c39_assert(all(x == ref$constant), paste("UNESTIMABLE_CONSTANT_COORDINATE_CHANGED", column))
                z <- matrix(numeric(), nrow = length(x), ncol = 0L)
            } else if (identical(ref$mode, "linear")) {
                z <- matrix(x, ncol = 1L)
            } else {
                c39_assert(identical(ref$mode, "ns"), "UNREVIEWED_SPLINE_MODE")
                z <- splines::ns(x, knots = ref$knots, Boundary.knots = ref$boundary, intercept = FALSE)
            }
            colnames(z) <- ref$column_names
            z
        }
        c39_spline_reference <- function(x, column) {
            c39_assert(is.numeric(x) && length(x) && all(is.finite(x)), "INVALID_SPLINE_INPUT")
            q <- unname(stats::quantile(x, c(0.05, 0.35, 0.65, 0.95), type = 7))
            boundary <- q[c(1, 4)]
            knots <- sort(unique(q[2:3]))
            knots <- knots[knots > boundary[1] & knots < boundary[2]]
            constant <- length(unique(x)) == 1L
            mode <- if (constant) "constant" else if (boundary[1] == boundary[2] || !length(knots)) "linear" else "ns"
            df <- switch(mode, constant = 0L, linear = 1L, ns = length(knots) + 1L)
            list(mode = mode, knots = knots, boundary = boundary, Boundary.knots = boundary, quantile_probabilities = c(0.05, 
                0.35, 0.65, 0.95), original_quantiles = q, effective_df = df, constant = if (constant) x[1L] else NULL, 
                constant_handling = if (constant) "FIXED_COORDINATE_OMIT_EQUIVALENT_DESIGN_ASSERT_AT_PREDICTION" else "NONE", 
                column_names = if (df) paste0(column, "__ns", seq_len(df)) else character())
        }
        c39_design <- function(d, members, base, reference = NULL) {
            discrete <- vapply(Filter(function(m) m$type != "continuous", members), `[[`, "", "column")
            form <- stats::reformulate(unique(c(base, discrete)), env = baseenv())
            X <- stats::model.matrix(form, d)
            continuous <- Filter(function(m) m$type == "continuous", members)
            if (is.null(reference)) {
                reference <- lapply(continuous, function(m) c39_spline_reference(d[[m$column]], m$column))
                attr(reference, "base_discrete") <- list(formula = form, columns = colnames(X), factor_levels = lapply(d[unique(c(base, 
                  discrete))], function(x) if (is.factor(x)) levels(x) else NULL))
            } else {
                c39_assert(identical(names(reference), names(continuous)), "FROZEN_SPLINE_MEMBERS_DRIFT")
                saved <- attr(reference, "base_discrete")
                c39_assert(is.list(saved) && identical(colnames(X), saved$columns), "FROZEN_BASE_DISCRETE_DESIGN_DRIFT")
                for (column in names(saved$factor_levels)) if (!is.null(saved$factor_levels[[column]])) c39_assert(is.factor(d[[column]]) && 
                  identical(levels(d[[column]]), saved$factor_levels[[column]]), paste("FROZEN_Y_FACTOR_LEVELS_DRIFT", 
                  column))
            }
            for (id in names(continuous)) {
                m <- continuous[[id]]
                X <- cbind(X, c39_basis(d[[m$column]], reference[[id]], m$column))
            }
            c39_assert(all(is.finite(X)), "NONFINITE_Y_DESIGN")
            if (is.null(attr(reference, "full_colnames"))) attr(reference, "full_colnames") <- colnames(X)
            c39_assert(identical(colnames(X), attr(reference, "full_colnames")), "FROZEN_Y_DESIGN_COLUMNS_DRIFT")
            list(X = X, reference = reference, fullcolnames = colnames(X))
        }
        c39_support <- function(d, members) {
            for (m in members) {
                x <- d[[m$column]]
                if (m$type == "continuous") c39_assert(is.numeric(x) && all(is.finite(x) & x >= m$lower & x <= 
                  m$upper), paste("PHYSICAL_SUPPORT_FAILURE_NO_CLIPPING", m$indicator_id)) else c39_assert(!anyNA(x) && 
                  all(as.character(x) %in% m$levels), paste("CATEGORY_SUPPORT_FAILURE", m$indicator_id))
            }
            ids <- vapply(members, `[[`, "", "indicator_id")
            if (all(c("blood_pressure__DBP", "blood_pressure__SBP") %in% ids)) {
                dbp <- members[[match("blood_pressure__DBP", ids)]]$column
                sbp <- members[[match("blood_pressure__SBP", ids)]]$column
                c39_assert(all(d[[sbp]] > d[[dbp]]), "SBP_NOT_STRICTLY_ABOVE_DBP_NO_REWRITE")
            }
            TRUE
        }
        c39_a_design <- function(d, base, previous, reference = NULL) {
            if (is.null(reference)) {
                constants <- previous[vapply(d[previous], function(x) is.numeric(x) && length(unique(x)) == 1L, 
                  logical(1))]
                fixed <- lapply(d[constants], function(x) x[1L])
                form <- stats::reformulate(unique(c(base, setdiff(previous, constants))), env = baseenv())
                X <- stats::model.matrix(form, d)
                reference <- list(formula = form, columns = colnames(X), fixed = fixed)
            } else {
                for (column in names(reference$fixed)) c39_assert(all(d[[column]] == reference$fixed[[column]]), 
                  paste("A_FIXED_PREDECESSOR_CHANGED", column))
                X <- stats::model.matrix(reference$formula, d)
                c39_assert(identical(colnames(X), reference$columns), "A_CONDITIONAL_DESIGN_DRIFT")
            }
            c39_assert(all(is.finite(X)), "A_NONFINITE_DESIGN")
            list(X = X, reference = reference)
        }
        c39_linear_fit <- function(X, y, w, id, rule) {
            if (length(unique(y)) == 1L) return(list(kind = "constant", value = y[1L], support = "FIXED_AMOUNT_SUPPORT_IN_THIS_SAMPLE", 
                rule = rule))
            c39_assert(qr(X)$rank == ncol(X), paste("A_DESIGN_RANK_DEFICIENT", id))
            fit <- stats::lm.wfit(X, y, w)
            c39_assert(fit$rank == ncol(X) && all(is.finite(fit$coefficients)), paste("A_CONTINUOUS_FIT_FAILED", 
                id))
            residual <- y - as.numeric(X %*% fit$coefficients)
            ix <- order(residual)
            mass <- w[ix]/sum(w)
            cumulative <- c(0, cumsum(mass))
            cumulative[length(cumulative)] <- 1
            list(kind = "continuous", coefficients = fit$coefficients, residuals = residual[ix], cumulative = cumulative, 
                cumulative_squared_mass = c(0, cumsum(mass^2)), rule = rule)
        }
        c39_logistic_fit <- function(X, positive, w, id) {
            if (length(unique(positive)) == 1L) return(list(kind = "constant", probability = as.numeric(positive[1L]), 
                support = "FIXED_ZERO_POSITIVE_STATUS_IN_THIS_SAMPLE", warnings = character()))
            c39_assert(qr(X)$rank == ncol(X), paste("A_ZERO_STATUS_DESIGN_RANK_DEFICIENT", id))
            warnings <- character()
            fit <- withCallingHandlers(stats::glm.fit(X, as.numeric(positive), weights = w, family = stats::binomial(), 
                control = stats::glm.control(maxit = 1000L)), warning = function(e) {
                warnings <<- c(warnings, conditionMessage(e))
                invokeRestart("muffleWarning")
            })
            c39_assert(isTRUE(fit$converged) && fit$rank == ncol(X) && all(is.finite(fit$coefficients)), paste("A_ZERO_STATUS_LOGISTIC_FAILED", 
                id))
            list(kind = "logistic", coefficients = fit$coefficients, iterations = fit$iter, converged = fit$converged, 
                rank = fit$rank, design_columns = ncol(X), boundary = fit$boundary, deviance = fit$deviance, fitted_probability_range = range(fit$fitted.values), 
                warnings = unique(warnings))
        }
        c39_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", 
            count_rule = "FIXED_SUPPORT_MULTINOMIAL_0_7") {
            c39_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "INVALID_JOINT_SELECTION_WEIGHT")
            c39_assert(continuous_rule %in% c("EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", "EMPIRICAL_RESIDUAL_UNTRUNCATED"), 
                "UNREVIEWED_CONTINUOUS_CONDITIONAL_LAW")
            c39_support(d, members)
            ids <- vapply(members, `[[`, "", "indicator_id")
            if ("blood_pressure__SBP" %in% ids) c39_assert("blood_pressure__DBP" %in% ids && match("blood_pressure__DBP", 
                ids) < match("blood_pressure__SBP", ids), "DBP_MUST_PRECEDE_SBP")
            previous <- character()
            fitted <- list()
            warnings <- character()
            for (m in members) {
                id <- m$indicator_id
                y <- d[[m$column]]
                count <- id %in% c("moderate_days", "walking_days")
                if (count) {
                  c39_assert(identical(count_rule, "FIXED_SUPPORT_MULTINOMIAL_0_7"), "ACTIVITY_REQUIRES_FIXED_0_7_MULTINOMIAL")
                  c39_assert(all(y %in% 0:7), "COUNT_OBSERVATION_OUTSIDE_FIXED_0_7")
                  d[[m$column]] <- factor(as.character(y), levels = as.character(0:7))
                  m$levels <- as.character(0:7)
                }
                design <- c39_a_design(d, base, previous)
                X <- design$X
                two_part <- id %in% c("cooked_veg", "dried_fruit", "raw_veg", "computer_hours", "tv_hours")
                if (two_part) {
                  c39_assert(m$type == "continuous" && all(y >= 0), paste("INVALID_TWO_PART_OBSERVATIONS", id))
                  positive <- y > 0
                  status <- c39_logistic_fit(X, positive, w, id)
                  amount <- if (any(positive)) c39_linear_fit(X[positive, , drop = FALSE], y[positive], w[positive], 
                    id, "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL") else NULL
                  f <- list(kind = "two_part", status = status, positive = amount, observed_zero_n = sum(!positive), 
                    observed_positive_n = sum(positive), weighted_observed_zero_fraction = sum(w[!positive])/sum(w), 
                    original_scale = "IDENTITY_FROZEN_UNITS", model_scale = "POSITIVE_IDENTITY_CONDITIONAL_MEAN", 
                    rule = "ZERO_ATOM_PLUS_CONDITIONALLY_LEGAL_WEIGHTED_POSITIVE_EMPIRICAL_RESIDUAL")
                  warnings <- c(warnings, status$warnings)
                } else if (m$type == "continuous" && !count) {
                  f <- c39_linear_fit(X, y, w, id, continuous_rule)
                  if (id == "blood_pressure__SBP") f$rule <- "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
                } else {
                  c39_assert(length(m$levels) >= 2L && all(table(factor(as.character(d[[m$column]]), levels = m$levels)) > 
                    0), paste("FROZEN_CATEGORY_ABSENT", id))
                  c39_assert(qr(X)$rank == ncol(X), paste("A_DESIGN_RANK_DEFICIENT", id))
                  form <- stats::reformulate(attr(stats::terms(design$reference$formula), "term.labels"), response = m$column, 
                    env = baseenv())
                  c39_assert(!".c39_fit_weights" %in% names(d), "RESERVED_INTERNAL_WEIGHT_COLUMN")
                  fit_data <- d
                  fit_data$.c39_fit_weights <- w
                  fit <- withCallingHandlers(nnet::multinom(form, data = fit_data, weights = .c39_fit_weights, 
                    trace = FALSE, maxit = 1000L, MaxNWts = 100000L, model = FALSE), warning = function(e) {
                    warnings <<- c(warnings, conditionMessage(e))
                    invokeRestart("muffleWarning")
                  })
                  c39_assert(identical(fit$convergence, 0L) && all(is.finite(stats::coef(fit))), paste("A_MULTINOM_NONCONVERGED", 
                    id))
                  f <- list(kind = if (count) "count" else "categorical", fit = fit)
                }
                f$formula <- design$reference$formula
                f$columns <- colnames(X)
                f$design_reference <- design$reference
                f$member <- m
                f$predecessor_columns <- previous
                fitted[[id]] <- f
                checkpoint <- getOption("c39.a_checkpoint_callback")
                if (is.function(checkpoint)) checkpoint(id, f, length(fitted))
                previous <- c(previous, m$column)
                if (count) d[[m$column]] <- y
            }
            list(models = fitted, members = members, warnings = unique(warnings), conditional_order = names(fitted), 
                continuous_rule = continuous_rule, model_revision = "CORE39_TWO_PART_SPLINE_BP_V2_1", residual_assumption = "Weighted empirical residual law conditional on legal support and every preceding coordinate; positive atoms only for two-part quantity")
        }
        c39_residual_draw <- function(f, X, u, lower, upper, lower_open, id, audit = NULL) {
            n <- nrow(X)
            lower <- rep_len(lower, n)
            upper <- rep_len(upper, n)
            lower_open <- rep_len(lower_open, n)
            if (identical(f$kind, "constant")) {
                value <- rep(f$value, n)
                legal <- ifelse(lower_open, value > lower, value >= lower) & value <= upper
                c39_assert(all(legal), paste("CONDITIONAL_CONSTANT_SUPPORT_EMPTY", id))
                if (is.function(audit)) audit(id, rep(1L, n), rep(1, n), rep(1, n))
                return(value)
            }
            mu <- as.numeric(X %*% f$coefficients)
            q <- u
            if (identical(f$rule, "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL")) {
                low <- findInterval(lower - mu, f$residuals, left.open = TRUE)
                low[lower_open] <- findInterval((lower - mu)[lower_open], f$residuals)
                high <- findInterval(upper - mu, f$residuals)
                lo <- f$cumulative[low + 1L]
                hi <- f$cumulative[high + 1L]
                mass <- hi - lo
                square_mass <- f$cumulative_squared_mass[high + 1L] - f$cumulative_squared_mass[low + 1L]
                ess <- ifelse(square_mass > 0, mass^2/square_mass, 0)
                if (is.function(audit)) audit(id, high - low, mass, ess)
                c39_assert(all(mass > 0 & high > low), paste("CONDITIONAL_SUPPORT_EMPTY", id))
                q <- lo + u * mass
                ix <- findInterval(q, f$cumulative, rightmost.closed = TRUE)
                upper_rounding <- q == hi & ix > high
                ix[upper_rounding] <- high[upper_rounding]
                c39_assert(all(ix > low & ix <= high), paste("CONDITIONAL_INVERSE_CDF_INVALID", id))
            } else {
                c39_assert(identical(f$rule, "EMPIRICAL_RESIDUAL_UNTRUNCATED"), "UNREVIEWED_RESIDUAL_LAW")
                ix <- findInterval(q, f$cumulative, rightmost.closed = TRUE)
            }
            c39_assert(all(ix >= 1L & ix <= length(f$residuals)), "RESIDUAL_INVERSE_CDF_INVALID")
            value <- mu + f$residuals[ix]
            legal <- ifelse(lower_open, value > lower, value >= lower) & value <= upper
            c39_assert(all(is.finite(value) & legal), paste("DRAW_ILLEGAL_NO_CLIPPING", id))
            value
        }
        c39_draw_a <- function(a, C, b, u, audit = NULL) {
            c39_assert(length(b) %in% c(1L, nrow(C)) && all(is.finite(b) & b >= 0 & b <= 1) && is.matrix(u) && 
                nrow(u) == nrow(C) && ncol(u) == length(a$models) && all(is.finite(u) & u >= 0 & u < 1), "INVALID_JOINT_DRAW_RANDOM_INPUT")
            d <- C
            d$B <- b
            for (j in seq_along(a$models)) {
                f <- a$models[[j]]
                m <- f$member
                id <- m$indicator_id
                X <- c39_a_design(d, NULL, NULL, f$design_reference)$X
                lower <- rep(m$lower, nrow(d))
                lower_open <- rep(FALSE, nrow(d))
                upper <- rep(m$upper, nrow(d))
                if (id == "blood_pressure__SBP") {
                  dbp <- a$models[["blood_pressure__DBP"]]$member$column
                  lower_open <- d[[dbp]] >= lower
                  lower <- pmax(lower, d[[dbp]])
                }
                if (f$kind == "two_part") {
                  p <- if (identical(f$status$kind, "constant")) rep(f$status$probability, nrow(d)) else stats::plogis(as.numeric(X %*% 
                    f$status$coefficients))
                  c39_assert(all(is.finite(p) & p >= 0 & p <= 1), "INVALID_ZERO_STATUS_PROBABILITY")
                  positive <- u[, j] >= 1 - p & p > 0
                  value <- numeric(nrow(d))
                  if (any(positive)) {
                    c39_assert(is.list(f$positive), paste("POSITIVE_STATUS_WITHOUT_AMOUNT_MODEL", id))
                    q <- (u[positive, j] - (1 - p[positive]))/p[positive]
                    c39_assert(all(q >= 0 & q <= 1), "INVALID_POSITIVE_CONDITIONAL_UNIFORM")
                    plower <- pmax(0, lower[positive])
                    popen <- lower[positive] <= 0
                    value[positive] <- c39_residual_draw(f$positive, X[positive, , drop = FALSE], q, plower, upper[positive], 
                      popen, id, audit)
                  }
                  d[[m$column]] <- value
                } else if (f$kind %in% c("continuous", "constant")) {
                  d[[m$column]] <- c39_residual_draw(f, X, u[, j], lower, upper, lower_open, id, audit)
                } else {
                  p <- fd_probs(f, d)
                  cumulative <- p
                  if (ncol(p) > 1L) for (k in 2:ncol(p)) cumulative[, k] <- cumulative[, k - 1L] + p[, k]
                  ix <- 1L + rowSums(u[, j] >= cumulative[, -ncol(cumulative), drop = FALSE])
                  d[[m$column]] <- factor(m$levels[ix], levels = m$levels)
                  if (f$kind == "count") d[[m$column]] <- as.numeric(as.character(d[[m$column]]))
                }
            }
            c39_support(d, a$members)
            d
        }
        c39_install_models <- function(engine) {
            c39_assert(is.environment(engine) && is.function(engine$fd_probs), "ISOLATED_SOURCE_LOCKED_ENGINE_REQUIRED")
            functions <- c("c39_assert", "c39_basis", "c39_spline_reference", "c39_design", "c39_support", "c39_a_design", 
                "c39_linear_fit", "c39_logistic_fit", "c39_fit_a", "c39_residual_draw", "c39_draw_a")
            source_environment <- environment(c39_install_models)
            for (name in functions) {
                f <- get(name, envir = source_environment, inherits = FALSE)
                environment(f) <- engine
                assign(name, f, envir = engine)
            }
            engine$fd_fit_a <- engine$c39_fit_a
            engine$fd_draw_a <- engine$c39_draw_a
            engine$fd_design <- engine$c39_design
            engine$fd_support <- engine$c39_support
            invisible(engine)
        }
    }), envir = context)
}

